# coding: utf-8
"""Implementation of the 'qc' command.
"""
import logging
from atropos.commands import (
    Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin, create_reader)
from atropos.commands.stats import (
    SingleEndReadStatistics, PairedEndReadStatistics)
from atropos.util import run_interruptible

class QcPipeline(Pipeline):
    """Base Pipeline for the qc command.
    """
    def __init__(self, read_statistics_class, **kwargs):
        super().__init__()
        self.read_statistics_class = read_statistics_class
        self.stats = {}
        self.stats_kwargs = kwargs
    
    def _get_stats(self, source):
        if source not in self.stats:
            self.stats[source] = self.read_statistics_class(**self.stats_kwargs)
        return self.stats[source]
    
    def handle_reads(self, context, read1, read2=None):
        self._get_stats(context['source']).collect(read1, read2)
    
    def finish(self, summary, worker=None):
        super().finish(summary)
        summary['pre'] = dict(
            (source, stats.summarize())
            for source, stats in self.stats.items())

class SingleEndQcPipeline(SingleEndPipelineMixin, QcPipeline):
    """QcPipeline for single-end data.
    """
    def __init__(self, **kwargs):
        super().__init__(SingleEndReadStatistics, **kwargs)

class PairedEndQcPipeline(PairedEndPipelineMixin, QcPipeline):
    """QcPipeline for paired-end data.
    """
    def __init__(self, **kwargs):
        super().__init__(PairedEndReadStatistics, **kwargs)

def execute(options, summary):
    """Execute the qc command.
    
    Args:
        options: Command-line options.
        summary: The summary dict.
    """
    reader, _, qualities, _ = create_reader(options)
    if options.paired:
        pipeline_class = PairedEndQcPipeline
    else:
        pipeline_class = SingleEndQcPipeline
    pipeline_args = dict(
        qualities=qualities,
        tile_key_regexp=options.tile_key_regexp)
    
    if options.threads is None:
        summary.update(mode='serial', threads=1)
        pipeline = pipeline_class(**pipeline_args)
        retcode = run_interruptible(pipeline, reader, summary)
    else:
        summary.update(mode='parallel', threads=options.threads)
        retcode = run_parallel(
            reader, pipeline_class, pipeline_args, summary, options.threads,
            options.process_timeout, options.read_queue_size)
    
    reader.close()
    
    return retcode

def run_parallel(
        reader, pipeline_class, pipeline_args, summary, threads=2,
        timeout=30, input_queue_size=0):
    """Execute qc in parallel mode.
    
    Args:
        reader: Iterator over batches of reads (most likely a BatchIterator)
        read_stats: Template ReadStatistics object.
        threads: Number of worker threads to use; additional threads are used
            for the main proccess and the writer process (if requested).
        timeout: number of seconds after which waiting processes escalate their
            messages from DEBUG to ERROR.
        input_queue_size: max number of items that can be in the input queue,
            or 0 for no limit (be warned that this could explode memory usage)
    """
    from multiprocessing import Queue
    from atropos.commands.multicore import (
        ParallelPipelineMixin, MulticoreError, launch_workers, ensure_processes,
        enqueue_all, dequeue, wait_on, wait_on_process, RETRY_INTERVAL)
    
    logging.getLogger().debug(
        "Starting atropos qc in parallel mode with threads=%d, timeout=%d",
        threads, timeout)
    
    if threads < 2:
        raise ValueError("'threads' must be >= 2")
    
    timeout = max(timeout, RETRY_INTERVAL)
    
    # Queue by which batches of reads are sent to worker processes
    input_queue = Queue(input_queue_size)
    # Queue for processes to send summary information back to main process
    summary_queue = Queue(threads)
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    pipeline_class = type(
        'QcPipelineImpl', (ParallelPipelineMixin, pipeline_class))
    pipeline = pipeline_class(**pipeline_args)
    worker_args = (input_queue, pipeline, summary_queue, timeout)
    worker_processes = launch_workers(threads - 1, worker_args)
    
    def ensure_alive():
        """Ensure that all worker processes have not terminated.
        """
        ensure_processes(worker_processes)
    
    def _run(worker_processes):
        # Add batches of reads to the input queue. Provide a timeout callback
        # to check that subprocesses are alive.
        num_batches = enqueue_all(reader, input_queue, timeout, ensure_alive)
        logging.getLogger().debug(
            "Main loop complete; saw %d batches", num_batches)
        
        # Tell the worker processes no more input is coming
        enqueue_all((None,) * threads, input_queue, timeout, ensure_alive)
        
        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        worker_processes.extend(
            launch_workers(1, worker_args, offset=threads-1))
        
        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            """Ensure that worker processes are still alive while waiting for
            summaries.
            """
            try:
                ensure_processes(
                    worker_processes,
                    "Workers are still alive and haven't returned summaries: {}",
                    alive=False)
            except Exception as err:
                logging.getLogger().error(err)
        
        wait_on(
            lambda: summary_queue.full(),
            wait_message="Waiting on worker summaries {}",
            timeout=timeout,
            wait=True,
            timeout_callback=summary_timeout_callback)
        
        # Process summary information from worker processes
        logging.getLogger().debug(
            "Processing summary information from worker processes")
        seen_summaries = set()
        seen_batches = set()
        def summary_fail_callback():
            """Raise exception if any summaries are missing.
            """
            missing_summaries = set(range(1, threads)) - seen_summaries
            raise MulticoreError(
                "Missing summaries from processes {}".format(
                    ",".join(str(s) for s in missing_summaries)))
        
        for _ in range(1, threads+1):
            batch = dequeue(summary_queue, fail_callback=summary_fail_callback)
            worker_index, worker_batches, worker_summary = batch
            if worker_summary is None:
                raise MulticoreError(
                    "Worker process {} died unexpectedly".format(worker_index))
            else:
                logging.getLogger().debug(
                    "Processing summary for worker %d", worker_index)
            seen_summaries.add(worker_index)
            seen_batches |= worker_batches
            summary.merge(worker_summary)
        
        # Check if any batches were missed
        if num_batches > 0:
            missing_batches = set(range(1, num_batches+1)) - seen_batches
            if len(missing_batches) > 0:
                raise MulticoreError(
                    "Workers did not process batches {}".format(
                        ",".join(str(b) for b in missing_batches)))
    
    retcode = run_interruptible(_run, worker_processes)
    
    # notify all threads that they should stop
    logging.getLogger().debug("Exiting all processes")
    def kill(process):
        """Kill a process if it fails to terminate on its own.
        """
        if retcode <= 1:
            wait_on_process(process, timeout, terminate=True)
        elif process.is_alive():
            process.terminate()
    for process in worker_processes:
        kill(process)
    
    return retcode
