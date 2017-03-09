# coding: utf-8
"""Implementation of the 'qc' command.
"""
import logging
from atropos.commands.multicore import *
from atropos.commands.stats import Summary, ReadStatistics

def qc(options, parser):
    from atropos.report.text import print_read_stats
    
    reader, names, qualities, _ = create_reader(options, parser)
    stats = ReadStatistics(
        'pre', options.paired, qualities=qualities,
        tile_key_regexp=options.tile_key_regexp)
    
    if options.threads is None:
        from atropos.qc import run_serial
        rc, report, details = run_serial(reader, stats)
    else:
        from atropos.qc import run_parallel
        rc, report, details = run_parallel(
            reader, stats, options.threads, options.process_timeout,
            options.read_queue_size)
    
    print_read_stats(options, report)
    return (rc, None, details)

def run_serial(reader, read_stats):
    def _run():
        for batch_size, batch in reader:
            for record in batch:
                read_stats.pre_trim(record)
    rc = run_interruptible(_run)
    report = read_stats.finish() if rc == 0 else None
    details = dict(mode='serial', threads=1)
    return (rc, report, details)

class QcWorker(WorkerProcess):
    def __init__(self, index, input_queue, summary_queue, timeout, read_stats):
        super().__init__(index, input_queue, summary_queue, timeout)
        self.read_stats = read_stats
        
    def _handle_records(self, batch_num, batch_size, records):
        for batch_size, batch in reader:
            for record in batch:
                self.read_stats.pre_trim(record)
    
    def _get_summary(self, error=False):
        worker_stats = None if error else self.read_stats
        return (self.index, self.seen_batches, worker_stats)

def run_parallel(reader, read_stats, threads=2, timeout=30, input_queue_size=0):
    """Execute qc in parallel mode.
    
    reader 				:: iterator over batches of reads (most likely a BatchIterator)
    read_stats          :: template ReadStatistics object.
    threads				:: number of worker threads to use; additional threads are used
                        for the main proccess and the writer process (if requested).
    timeout				:: number of seconds after which waiting processes escalate their
                        messages from DEBUG to ERROR.
    input_queue_size 	:: max number of items that can be in the input queue, or 0 for
                        no limit (be warned that this could explode memory usage)
    
    TODO: There is a lot of common code between this and atropos.trim.run_parallel.
    See if I can refactor.
    """
    logging.getLogger().debug(
        "Starting atropos qc in parallel mode with threads={}, timeout={}".format(threads, timeout))
    
    if threads < 2:
        raise ValueError("'threads' must be >= 2")
    
    timeout = max(timeout, RETRY_INTERVAL)
    
    # Queue by which batches of reads are sent to worker processes
    input_queue = Queue(input_queue_size)
    # Queue for processes to send summary information back to main process
    summary_queue = Queue(threads)
    # Aggregate summary
    summary = Summary()
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    worker_args = (input_queue, summary_queue, timeout, read_stats)
    worker_processes = launch_workers(threads - 1, QcWorker, worker_args)
    
    def ensure_alive():
        ensure_processes(worker_processes)
    
    def _run():
        # Add batches of reads to the input queue. Provide a timeout callback
        # to check that subprocesses are alive.
        num_batches = enqueue_all(enumerate(reader, 1), input_queue, timeout, ensure_alive)
        logging.getLogger().debug(
            "Main loop complete; saw {} batches".format(num_batches))
        
        # Tell the worker processes no more input is coming
        enqueue_all((None,) * threads, input_queue, timeout, ensure_alive)
        
        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        worker_processes += launch_workers(1, QcWorker, worker_args, offset=threads-1)
        
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
        logging.getLogger().debug("Processing summary information from worker processes")
        seen_summaries = set()
        seen_batches = set()
        def summary_fail_callback():
            missing_summaries = set(range(1, threads)) - seen_summaries
            raise Exception("Missing summaries from processes {}".format(
                ",".join(str(s) for s in missing)))
        
        for i in range(1, threads+1):
            batch = dequeue(
                summary_queue,
                fail_callback=summary_fail_callback)
            worker_index, worker_batches, worker_stats = batch
            if worker_stats is None:
                raise Exception("Worker process {} died unexpectedly".format(worker_index))
            else:
                logging.getLogger().debug("Processing summary for worker {}".format(worker_index))
            seen_summaries.add(worker_index)
            seen_batches |= worker_batches
            summary.add_read_stats(read_stats)
        
        # Check if any batches were missed
        if num_batches > 0:
            missing_batches = set(range(1, num_batches+1)) - seen_batches
            if len(missing_batches) > 0:
                raise Exception("Workers did not process batches {}".format(
                    ",".join(str(b) for b in missing_batches)))
    
    try:
        rc = run_interruptible(_run, worker_processes)
    finally:
        # notify all threads that they should stop
        logging.getLogger().debug("Exiting all processes")
        def kill(process):
            if rc <= 1:
                wait_on_process(process, terminate=True)
            elif process.is_alive():
                process.terminate()
        for process in worker_processes:
            kill(process)
    
    report = summary.finish() if rc == 0 else None
    details = dict(mode='parallel', threads=threads)
    return (rc, report, details)
