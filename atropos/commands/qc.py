# coding: utf-8
"""Implementation of the 'qc' command.
"""
import logging
from atropos.commands import (
    Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin, create_reader)
from atropos.reports.txt import print_read_stats
from atropos.commands.stats import Summary, ReadStatistics

class QcPipeline(Pipeline):
    def __init__(self, **kwargs):
        super().__init__()
        self.stats = {}
        self.stats_kwargs = kwargs
    
    def _get_stats(self, source):
        if source not in self.stats:
            self.stats[source] = ReadStatCollector(**self.stats_kwargs)
        return self.stats[source]

class SingleEndQcPipeline(SingleEndPipelineMixin, QcPipeline):
    def handle_reads(self, context, read):
        self._get_stats(context['batch_source']).collect(read)

class PairedEndQcPipeline(PairedEndPipelineMixin, QcPipeline):
    def handle_reads(self, context, read1, read2):
        src1, src2 = context['batch_source']
        self._get_stats(src1).collect(read1)
        self._get_stats(src2).collect(read2)

def execute(options, summary):
    reader, names, qualities, _ = create_reader(options)
    
    if options.paired:
        pipeline_class = PairedEndQcPipeline
    else:
        pipeline_class = SingleEndQcPipeline
    
    pipeline_args = dict(
        qualities=qualities,
        tile_key_regexp=options.tile_key_regexp)
    
    if options.threads is None:
        summary.update(mode='serial', threads=1)
        rc = run_interruptible(pipeline_class(**pipeline_args), reader, summary)
    else:
        summary.update(mode='parallel', threads=options.threads)
        rc = run_parallel(
            reader, pipeline_args, summary, options.threads,
            options.process_timeout, options.read_queue_size)
    
    return rc

def run_parallel(reader, pipeline_args, summary, threads=2, timeout=30, input_queue_size=0):
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
        enqueue_all, dequeue)
    
    pipeline_class = type('QcPipelineImpl', (ParallelPipelineMixin, pipeline_class))
    
    logging.getLogger().debug(
        "Starting atropos qc in parallel mode with threads={}, timeout={}".format(
        threads, timeout))
    
    if threads < 2:
        raise ValueError("'threads' must be >= 2")
    
    timeout = max(timeout, RETRY_INTERVAL)
    
    # Queue by which batches of reads are sent to worker processes
    input_queue = Queue(input_queue_size)
    # Queue for processes to send summary information back to main process
    summary_queue = Queue(threads)
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    worker_args = (input_queue, summary_queue, timeout, read_stats)
    worker_processes = launch_workers(
        threads - 1, worker_args, worker_class=QcWorker)
    
    def ensure_alive():
        ensure_processes(worker_processes)
    
    def _run():
        # Add batches of reads to the input queue. Provide a timeout callback
        # to check that subprocesses are alive.
        num_batches = enqueue_all(reader, input_queue, timeout, ensure_alive)
        logging.getLogger().debug(
            "Main loop complete; saw {} batches".format(num_batches))
        
        # Tell the worker processes no more input is coming
        enqueue_all((None,) * threads, input_queue, timeout, ensure_alive)
        
        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        worker_processes += launch_workers(
            1, worker_args, offset=threads-1, worker_class=QcWorker)
        
        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            try:
                ensure_processes(worker_processes,
                    "Workers are still alive and haven't returned summaries: {}",
                    alive=False)
            except Exception as e:
                logging.getLogger().error(e)
        
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
            missing_summaries = set(range(1, threads)) - seen_summaries
            raise MulticoreError("Missing summaries from processes {}".format(
                ",".join(str(s) for s in missing)))
        
        for i in range(1, threads+1):
            batch = dequeue(summary_queue, fail_callback=summary_fail_callback)
            worker_index, worker_batches, worker_summary = batch
            if worker_summary is None:
                raise MulticoreError(
                    "Worker process {} died unexpectedly".format(worker_index))
            else:
                logging.getLogger().debug(
                    "Processing summary for worker {}".format(worker_index))
            seen_summaries.add(worker_index)
            seen_batches |= worker_batches
            summary.merge(worker_summary)
        
        # Check if any batches were missed
        if num_batches > 0:
            missing_batches = set(range(1, num_batches+1)) - seen_batches
            if len(missing_batches) > 0:
                raise MulticoreError("Workers did not process batches {}".format(
                    ",".join(str(b) for b in missing_batches)))
    
    rc = run_interruptible(_run, worker_processes)
    
    # notify all threads that they should stop
    logging.getLogger().debug("Exiting all processes")
    def kill(process):
        if rc <= 1:
            wait_on_process(process, terminate=True)
        elif process.is_alive():
            process.terminate()
    for process in worker_processes:
        kill(process)
    
    return rc
