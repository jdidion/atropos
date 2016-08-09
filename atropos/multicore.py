# Parallel implementation of run_atropos. Works as follows:
#
# 1. Main thread creates N worker processes (where N is the number of threads to be allocated) and
# (optionally) one writer process.
# 2. Main thread loads batches of reads (or read pairs) from input file(s) and adds them to a queue
# (the input queue).
# 3. Worker processes take batches from the input queue, process them as atropos normally does,
# and either add the results to the result queue (if using a writer process) or write the results
# to disk. Each result is a dict mapping output file names to strings, where each string is a
# concatenation of reads (with appropriate line endings) to be written. A parameter also controls
# whether data compression is done by the workers or the writer.
# 4. If using a writer process, it takes results from the result queue and writes each string to
# its corresponding file.
# 5. When the main process finishes loading reads from the input file(s), it sends a signal to the
# worker processes that they should complete when the input queue is empty. It also singals the
# writer process how many total batches to expect, and the writer process exits after it has
# processed that many batches.
# 6. When a worker process completes, it adds a summary of its activity to the summary queue.
# 7. The main process reads summaries from the summary queue and merges them to create the complete
# summary, which is used to generate the report.
#
# There are several possible points of failure:
#
# 1. The main process may exit due to an unexpected error, or becuase the user forces it to exit
# (Ctrl-C). In this case, an attempt is made to cancel all processes before exiting.
# 2. A worker or writer process may exit due to an unknown error. To handle this, the main process
# checks that each process is alive whenver it times out writing to the input queue, and again when
# waiting for worker summaries. If a process has died, the program exits with an error since some data
# might have gotten lost.
# 3. More commonly, process will time out blocking on reading from or writing to a queue. Size
# limits are used (optionally) for the input and result queues to prevent using lots of memory. When
# few threads are allocated, it is most likely that the main and writer processes will block, whereas
# with many threads allocated the workers are most likely to block. Also, e.g. in a cluster
# environment, I/O latency may cause a "backup" resulting in frequent blocking of the main and workers
# processes. Finally, also e.g. in a cluster environment, processes may suspended for periods of time.
# Use of a hard timeout period, after which processes are forced to exit, is thus undesirable.
# Instead, parameters are provided for the user to tune the batch size and max queue sizes to their
# particular environment. Additionally, a "soft" timeout is used, after which log messages are
# escallated from DEBUG to ERROR level. The user can then make the decision of whether or not to kill
# the program.

from collections import defaultdict
import inspect
import logging
import os
from queue import Empty, Full
import time
from multiprocessing import Process, Queue, Value, cpu_count

from .report import *
from .seqio import Writers, FormatError
from .compression import get_compressor, can_use_system_compression

__author__ = "John Didion"

# Max time to wait between retrying operations
RETRY_INTERVAL = 5

def wait_on(condition, wait_message="Waiting {}", timeout=None, fail_callback=None, wait=None, timeout_callback=None):
    if wait is True:
        wait = lambda: time.sleep(RETRY_INTERVAL)
    elif isinstance(wait, int):
        wait_time = wait
        wait = lambda: time.sleep(wait_time)
    wait_start = None
    while True:
        result = condition()
        if result is not False:
            return result
        if fail_callback:
            fail_callback()
        now = time.time()
        if not wait_start:
            wait_start = now
        else:
            waiting = now - wait_start
            msg = wait_message.format("for {} seconds".format(round(waiting, 1)))
            if timeout is not None and waiting >= timeout:
                logging.getLogger().error(msg)
                if timeout_callback:
                    if inspect.isclass(timeout_callback):
                        raise timeout_callback()
                    else:
                        timeout_callback()
            else:
                logging.getLogger().debug(msg)
            if wait:
                wait()

def enqueue(queue, item, wait_message="Waiting to enqueue item {}", block_timeout=RETRY_INTERVAL, **kwargs):
    def condition():
        try:
            queue.put(item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False
    wait_on(condition, wait_message=wait_message, **kwargs)

def dequeue(queue, wait_message="Waiting to dequeue item {}", block_timeout=RETRY_INTERVAL, **kwargs):
    def condition():
        try:
            return queue.get(block=True, timeout=block_timeout)
        except Empty:
            return False
    return wait_on(condition, wait_message=wait_message, **kwargs)

class Control(object):
    def __init__(self, initial_value):
        self.control = Value('l', initial_value)

    def check_value(self, value, lock=False):
        return self.get_value(lock=lock) == value

    def check_value_positive(self, lock=False):
        return self.get_value(lock=lock) > 0

    def get_value(self, lock=True):
        if lock:
            with self.control.get_lock():
                return self.control.value
        else:
            return self.control.value

    def set_value(self, value):
        with self.control.get_lock():
            self.control.value = value

class WorkerProcess(Process):
    def __init__(self, index, modifiers, filters, formatters, input_queue,
                 result_handler, summary_queue, timeout):
        """
        input_queue: queue with batches of records to process
        output_queue: queue where results are written
        summary_queue: queue where summary information is written
        control: shared value with value 0 unless the thread
        should die (< 0) or there are no more batches coming (> 0).
        """
        super(WorkerProcess, self).__init__(name="Worker process {}".format(index))
        self.index = index
        self.modifiers = modifiers
        self.filters = filters
        self.formatters = formatters
        self.input_queue = input_queue
        self.result_handler = result_handler
        self.summary_queue = summary_queue
        self.timeout = timeout
        self.processed_reads = 0
        self.total_bp1 = 0
        self.total_bp2 = 0
        self.seen_batches = set()
    
    def run(self):
        logging.getLogger().debug("{} running under pid {}".format(self.name, os.getpid()))
        
        def iter_batches():
            while True:
                batch = dequeue(
                    self.input_queue,
                    wait_message="{} waiting on batch {{}}".format(self.name),
                    timeout=self.timeout)
                yield batch
        
        def enqueue_summary(error=False):
            if error:
                process_stats = adapter_stats = None
            else:
                process_stats = collect_process_statistics(
                    self.processed_reads, self.total_bp1, self.total_bp2,
                    self.modifiers, self.filters, self.formatters)
                adapter_stats = summarize_adapters(self.modifiers.get_modifiers(AdapterCutter))
            enqueue(
                self.summary_queue,
                (self.index, self.seen_batches, process_stats, adapter_stats),
                wait_message="{} waiting to queue summary {{}}".format(self.name),
                timeout=self.timeout
            )
        
        try:
            self.result_handler.start(self)
            
            for batch in iter_batches():
                if batch is None:
                    break
                
                batch_num, (batch_size, records) = batch
                self.seen_batches.add(batch_num)
                self.processed_reads += batch_size
                
                logging.getLogger().debug("{} processing batch of size {}".format(self.name, batch_size))
                
                result = defaultdict(lambda: [])
                for record in records:
                    reads, bp = self.modifiers.modify(record)
                    self.total_bp1 += bp[0]
                    self.total_bp2 += bp[1]
                    dest = self.filters.filter(*reads)
                    self.formatters.format(result, dest, *reads)
                
                self.result_handler.write_result(batch_num, result)
            
            logging.getLogger().debug("{} finished; processed {} batches, {} reads ({},{} bp)".format(
                self.name, len(self.seen_batches), self.processed_reads, self.total_bp1, self.total_bp2))
            
            logging.getLogger().debug("{} sending summary".format(self.name))
            enqueue_summary()
            
            logging.getLogger().debug("{} exiting normally".format(self.name))
        except:
            logging.getLogger().error("Unexpected error in {}".format(self.name), exc_info=True)
            enqueue_summary(error=True)
        finally:
            self.result_handler.finish()

# Control values
# Writer is running normally
WRITER_ACTIVE = 0
# Writer should exit (usu. because there were no reads in the user input)
WRITER_ERROR = -1

class WriterProcess(Process):
    """
    Thread that accepts results from the worker threads
    and writes the results to output file(s). Each batch
    is expected to be (batch_num, path, records), where
    path is the destination file and records is a string.
    Not guaranteed to preserve the original order of
    sequence records.

    writers: Writers object, which manages opening and writing
    to files.
    queue: input queue.
    control: a shared value for communcation with the main
    process.
    timeout: seconds to wait for next batch before complaining
    """
    def __init__(self, result_handler, queue, control, timeout=60):
        super(WriterProcess, self).__init__(name="Writer process")
        self.result_handler = result_handler
        self.queue = queue
        self.control = control
        self.timeout = timeout
        self.seen_batches = set()
        self.num_batches = None

    def run(self):
        logging.getLogger().debug("Writer process running under pid {}".format(self.name, os.getpid()))
        
        class Done(Exception): pass
        class Killed(Exception): pass

        def fail_callback():
            if self.num_batches is None and self.control.check_value_positive():
                self.num_batches = self.control.get_value()
            if self.num_batches is not None and len(self.seen_batches) >= self.num_batches:
                raise Done()

        def timeout_callback():
            if self.num_batches is not None:
                missing = set(range(1, self.num_batches+1)) - self.seen_batches
                logging.getLogger().error("Writer thread still missing batches {} of {}".format(
                    ",".join(str(i) for i in missing), self.num_batches))

        def iter_batches():
            while True:
                batch = dequeue(
                    self.queue,
                    wait_message="Writer waiting on result {}",
                    timeout=self.timeout,
                    fail_callback=fail_callback,
                    timeout_callback=timeout_callback)
                yield batch

        try:
            self.result_handler.start(self)
            
            for batch_num, result in iter_batches():
                self.seen_batches.add(batch_num)
                self.result_handler.write_result(batch_num, result)
        except Done:
            logging.getLogger().debug("Writer process exiting normally")
        except Killed:
            logging.getLogger().debug("Writer process exited early")
        except:
            logging.getLogger().error("Unexpected error in writer process", exc_info=True)
            self.control.set_value(WRITER_ERROR)
        finally:
            num_batches = self.control.get_value(lock=True)
            self.result_handler.finish(num_batches if num_batches > 0 else None)

class WorkerResultHandler(object):
    """
    Wraps a ResultHandler and compresses results prior
    to writing.
    """
    def __init__(self, handler):
        self.handler = handler

    def start(self, worker):
        self.handler.start(worker)
    
    def write_result(self, batch_num, result):
        """
        Give a dict mapping files to lists of strings,
        join the strings and compress them (if necessary)
        and then return the property formatted result
        dict.
        """
        self.handler.write_result(batch_num, dict(self.prepare_file(*item)
            for item in result.items()))
    
    def prepare_file(self, path, strings):
        return (path, "".join(strings))
    
    def finish(self, total_batches=None):
        self.handler.finish()

class CompressingWorkerResultHandler(WorkerResultHandler):
    """
    Wraps a ResultHandler and compresses results prior
    to writing.
    """
    def start(self, worker):
        super(CompressingWorkerResultHandler, self).start(worker)
        self.file_compressors = {}
    
    def prepare_file(self, path, strings):
        compressor = self.get_compressor(path)
        if compressor:
            return ((path, 'wb'), compressor.compress(b''.join(
                s.encode() for s in strings)))
        else:
            return ((path, 'wt'), "".join(strings))
    
    def get_compressor(self, filename):
        if filename not in self.file_compressors:
            self.file_compressors[filename] = get_compressor(filename)
        return self.file_compressors[filename]

class QueueResultHandler(object):
    """
    ResultHandler that writes results to the output queue.
    """
    def __init__(self, queue):
        self.queue = queue

    def start(self, worker):
        self.message = "{} waiting to queue result {{}}".format(worker.name)
        self.timeout = worker.timeout

    def write_result(self, batch_num, result):
        enqueue(
            self.queue,
            (batch_num, result),
            wait_message=self.message,
            timeout=self.timeout
        )

    def finish(self, total_batches=None):
        pass

class WriterResultHandler(object):
    """
    ResultHandler that writes results to disk.
    """
    def __init__(self, writers, compressed=False, use_suffix=False):
        self.writers = writers
        self.compressed = compressed
        self.use_suffix = use_suffix

    def start(self, worker):
        if self.use_suffix:
            self.writers.suffix = ".{}".format(worker.index)

    def write_result(self, batch_num, result):
        self.writers.write_result(result, self.compressed)

    def finish(self, total_batches=None):
        self.writers.close()

class OrderPreservingWriterResultHandler(WriterResultHandler):
    """
    Writer thread that is less time/memory efficient, but is
    guaranteed to preserve the original order of records.
    """
    def start(self, worker):
        super(OrderPreservingWriterResultHandler, self).__init__(worker)
        self.pending = PendingQueue()
        self.cur_batch = 1

    def write_result(self, batch_num, result):
        if batch_num == self.cur_batch:
            self.writers.write_result(result, self.compressed)
            self.cur_batch += 1
            self.consume_pending()
        else:
            self.pending.push(batch_num, result)

    def finish(self, total_batches):
        if total_batches is not None:
            self.consume_pending()
            if self.cur_batch != total_batches:
                raise Exception("OrderPreservingWriterResultHandler finishing without having seen "
                                "{} batches".format(total_batches))
        self.writers.close()

    def consume_pending(self):
        while (not self.pending.empty) and (self.cur_batch == pending.min_priority):
            self.writers.write_result(pending.pop(), self.compressed)
            self.cur_batch += 1

class PendingQueue(object):
    def _init(self, max_size=None):
        self.queue = {}
        self.max_size = max_size
        self.min_priority = None

    def push(self, priority, value):
        assert not self.full
        self.queue[priority] = value
        if self.min_priority is None or priority < self.min_priority:
            self.min_priority = priority

    def pop(self):
        assert not self.empty
        value = self.queue.pop(self.min_priority)
        if self.empty:
            self.min_priority = None
        else:
            self.min_priority = min(self.queue.keys())
        return value

    @property
    def full(self):
        return self.max_size and len(self.queue) >= self.max_size

    @property
    def empty(self):
        return len(self.queue) == 0

def run_parallel(reader, modifiers, filters, formatters, writers, threads=2, timeout=30,
                 preserve_order=False, input_queue_size=0, result_queue_size=0,
                 use_writer_process=True, compression=None):
    """
    Execute atropos in parallel mode.

    reader 				:: iterator over batches of reads (most likely a BatchIterator)
    modifiers 			::
    filters 			::
    formatters 			::
    writers				::
    threads				:: number of worker threads to use; additional threads are used
                        for the main proccess and the writer process (if requested).
    timeout				:: number of seconds after which waiting processes escalate their
                        messages from DEBUG to ERROR.
    preserve_order 		:: whether to preserve the input order of reads when writing
                        (only valid when `use_writer_process=True`)
    input_queue_size 	:: max number of items that can be in the input queue, or 0 for
                        no limit (be warned that this could explode memory usage)
    result_queue_size	:: max number of items that can be in the result queue, or 0 for
                        no limit (be warned that this could explode memory usage)
    use_writer_process	:: if True, a separate thread will be used to write results to
                        disk. Otherwise, each worker thread will write its results to
                        an output file with a '.N' extension, where N is the thread index.
                        This is useful in cases where the I/O is the main bottleneck.
    compression	        If "writer", the writer process perform data compression, otherwise
                        the worker processes performs compression.
    """
    logging.getLogger().debug(
        "Starting atropos in parallel mode with threads={}, timeout={}".format(threads, timeout))
    
    assert threads >= 2
    
    # Reserve a thread for the writer process if it will be doing the compression and if one is available.
    if compression is None:
        compression = "writer" if use_writer_process and can_use_system_compression() else "worker"
    if compression == "writer" and threads > 2:
        threads -= 1
    
    timeout = max(timeout, RETRY_INTERVAL)

    # Queue by which batches of reads are sent to worker processes
    input_queue = Queue(input_queue_size)
    # Queue by which results are sent from the worker processes to the writer process
    result_queue = Queue(result_queue_size)
    # Queue for processes to send summary information back to main process
    summary_queue = Queue(threads)
    # Aggregate summary
    summary = Summary(trimmer_classes=modifiers.get_trimmer_classes())
    # Return code (0=normal, anything else is an error)
    rc = 0

    if use_writer_process:
        worker_result_handler = QueueResultHandler(result_queue)
        if compression == "writer":
            worker_result_handler = WorkerResultHandler(worker_result_handler)
        else:
            worker_result_handler = CompressingWorkerResultHandler(worker_result_handler)
        
        # Shared variable for communicating with writer thread
        writer_control = Control(WRITER_ACTIVE)
        # result handler
        if preserve_order:
            writer_result_handler = OrderPreservingWriterResultHandler(
                writers, compressed=compression == "worker")
        else:
            writer_result_handler = WriterResultHandler(
                writers, compressed=compression == "worker")
        # writer process
        writer_process = WriterProcess(writer_result_handler, result_queue, writer_control, timeout)
        writer_process.start()
    else:
        worker_result_handler = WorkerResultHandler(WriterResultHandler(writers, use_suffix=True))

    # worker processes
    def launch_workers(n, offset=0):
        logging.getLogger().info("Starting {} worker processes".format(n))
        workers = [
            WorkerProcess(
                i+offset, modifiers, filters, formatters, input_queue,
                worker_result_handler, summary_queue, timeout)
            for i in range(n)
        ]
        # start workers
        for worker in workers: worker.start()
        return workers
    
    def ensure_alive():
        alive = [worker.is_alive() for worker in worker_processes]
        if not all(alive):
            raise Exception("One or more worker process exited: {}".format(",".join(
                str(i) for i in range(len(alive)) if not alive[i])))
        if use_writer_process and not (writer_process.is_alive() and writer_control.check_value(WRITER_ACTIVE)):
            raise Exception("Writer process exited")

    def enqueue_all(iterable):
        num_items = 0
        for item in iterable:
            def condition():
                try:
                    input_queue.put(item, block=True, timeout=RETRY_INTERVAL)
                    return True
                except Full:
                    return False
            wait_on(
                condition,
                wait_message="Main process waiting to queue item {}",
                timeout=timeout,
                fail_callback=ensure_alive)
            num_items += 1
        return num_items

    def wait_on_process(process, terminate=False):
        timeout_callback = lambda: process.terminate() if terminate else None
        wait_on(
            lambda: not process.is_alive(),
            wait_message="Waiting on {} to terminate {{}}".format(process.name),
            timeout=timeout,
            wait=lambda: process.join(RETRY_INTERVAL),
            timeout_callback=timeout_callback)
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    worker_processes = launch_workers(threads - 1)
    
    try:
        # Add batches of reads to the input queue. Provide a timeout callback
        # to check that subprocesses are alive.
        num_batches = enqueue_all(enumerate(reader, 1))
        logging.getLogger().debug(
            "Main loop complete; saw {} batches".format(num_batches))
        
        # Tell the worker processes no more input is coming
        enqueue_all((None,) * threads)
        
        # Tell the writer thread the max number of batches to expect
        if use_writer_process:
            writer_control.set_value(num_batches)
        
        # Wait for workers to exit
        #for worker in worker_processes:
        #	# Wait for worker to exit
        #	wait_on_process(worker)
        
        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        worker_processes += launch_workers(1, threads-1)
        
        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            alive = [worker.is_alive() for worker in worker_processes]
            if any(alive):
                missing = ",".join(str(i) for i in range(len(alive)) if alive[i])
                logging.getLogger().error(
                    "Workers are still alive and haven't returned summaries: {}".format(missing))
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
            worker_index, worker_batches, process_stats, adapter_stats = batch
            if process_stats is None or adapter_stats is None:
                raise Exception("Worker process {} died unexpectedly".format(worker_index))
            else:
                logging.getLogger().debug("Processing summary for worker {}".format(worker_index))
            seen_summaries.add(worker_index)
            seen_batches |= worker_batches
            summary.add_process_stats(process_stats)
            summary.add_adapter_stats(adapter_stats)

        # Check if any batches were missed
        if num_batches > 0:
            missing_batches = set(range(1, num_batches+1)) - seen_batches
            if len(missing_batches) > 0:
                raise Exception("Workers did not process batches {}".format(
                    ",".join(str(b) for b in missing_batches)))

        if use_writer_process:
            # Wait for writer to complete
            wait_on_process(writer_process)

    except KeyboardInterrupt as e:
        logging.getLogger().error("Interrupted")
        rc = 130

    except IOError as e:
        if e.errno == errno.EPIPE:
            rc = 1
        else:
            raise

    except (FormatError, EOFError) as e:
        logging.getLogger().error("Cutadapt error", exc_info=True)
        rc = 1

    except Exception as e:
        logging.getLogger().error("Unknown error", exc_info=True)
        rc = 1

    finally:
        logging.getLogger().debug("Waiting for reader to close")
        reader.close()

        # notify all threads that they should stop
        logging.getLogger().debug("Exiting all processes")
        def kill(process):
            if rc <= 1:
                wait_on_process(process, terminate=True)
            elif process.is_alive():
                process.terminate()
        for process in worker_processes:
            kill(process)
        if use_writer_process:
            kill(writer_process)
    
    report = summary.finish() if rc == 0 else None
    return (rc, report)
