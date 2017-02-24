"""Classes and methods to support parallel operations.
"""
import inspect
import logging
from multiprocessing import Process, Queue, Value, cpu_count
import os
from queue import Empty, Full
import time
from atropos.util import run_interruptible

RETRY_INTERVAL = 5
"""Max time to wait between retrying operations."""

# Control values
CONTROL_ACTIVE = 0
"""Controlled process should run normally."""
CONTROL_ERROR = -1
"""Controlled process should exit."""

class Control(object):
    def __init__(self, initial_value=CONTROL_ACTIVE):
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

class PendingQueue(object):
    def _init(self, max_size=None):
        self.queue = {}
        self.max_size = max_size
        self.min_priority = None
    
    def push(self, priority, value):
        if self.full:
            raise Full()
        self.queue[priority] = value
        if self.min_priority is None or priority < self.min_priority:
            self.min_priority = priority
    
    def pop(self):
        if self.empty:
            raise Empty()
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

class WorkerProcess(Process):
    """Parent class for worker processes that operate on batches of reads.
    
    Args:
        index: A unique ID for the process
        input_queue: queue with batches of records to process
        summary_queue: queue where summary information is written
        timeout: time to wait upon queue full/empty
    """
    def __init__(self, index, input_queue, summary_queue, timeout):
        super().__init__(name="Worker process {}".format(index))
        self.index = index
        self.input_queue = input_queue
        self.summary_queue = summary_queue
        self.timeout = timeout
        self.processed_reads = 0
        self.seen_batches = set()
    
    def _on_start(self):
        pass
    
    def _on_finish(self):
        pass
    
    def _handle_records(self, batch_num, batch_size, records):
        raise NotImplementedError()
    
    def _get_summary(self, error=False):
        return (self.index, self.seen_batches)
    
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
            enqueue(
                self.summary_queue, self._get_summary(error),
                wait_message="{} waiting to queue summary {{}}".format(self.name),
                timeout=self.timeout
            )
        
        try:
            self._on_start()
            
            for batch in iter_batches():
                if batch is None:
                    break
                
                batch_num, (batch_size, records) = batch
                self.seen_batches.add(batch_num)
                self.processed_reads += batch_size
                
                logging.getLogger().debug("{} processing batch of size {}".format(self.name, batch_size))
                self._handle_records(batch_num, batch_size, records)
            
            logging.getLogger().debug("{} finished; processed {} batches, {} reads".format(
                self.name, len(self.seen_batches), self.processed_reads))
            
            logging.getLogger().debug("{} sending summary".format(self.name))
            enqueue_summary()
            
            logging.getLogger().debug("{} exiting normally".format(self.name))
        except:
            logging.getLogger().error("Unexpected error in {}".format(self.name), exc_info=True)
            enqueue_summary(error=True)
        finally:
            self._on_finish()

class ResultHandlerWorkerProcess(WorkerProcess):
    """Parent class for worker processes that operate on batches of reads and
    handles each record.
    
    Args:
        index: A unique ID for the process
        input_queue: queue with batches of records to process
        result_handler:
        summary_queue: queue where summary information is written
        timeout: time to wait upon queue full/empty
    """
    def __init__(self, index, input_queue, summary_queue, timeout, result_handler):
        super().__init__(index, input_queue, summary_queue, timeout)
        self.result_handler = result_handler
    
    def _on_start(self):
        self.result_handler.start(self)
    
    def _on_finish(self):
        self.result_handler.finish()
    
    def _handle_records(self, batch_num, batch_size, records):
        result = defaultdict(lambda: [])
        for record in records:
            self._handle_record(record, result)
        self.result_handler.write_result(batch_num, result)
    
    def _handle_record(self, record, result):
        raise NotImplementedError()

class ResultProcess(Process):
    """Thread that accepts results from the worker threads and process them
    using a ResultHandler. Each batch is expected to be
    (batch_num, path, records), where path is the destination file and records
    is a string. Not guaranteed to preserve the original order of sequence records.

    result_handler: A ResultHandler object.
    queue: input queue.
    control: a shared value for communcation with the main process.
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
                logging.getLogger().error("Result thread still missing batches {} of {}".format(
                    ",".join(str(i) for i in missing), self.num_batches))
        
        def iter_batches():
            while True:
                batch = dequeue(
                    self.queue,
                    wait_message="Result process waiting on result {}",
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
            self.control.set_value(CONTROL_ERROR)
        finally:
            num_batches = self.control.get_value(lock=True)
            self.result_handler.finish(num_batches if num_batches > 0 else None)

class ResultHandler(object):
    def start(self, worker):
        pass
    
    def finish(self, total_batches=None):
        pass
    
    def write_result(self, batch_num, result):
        raise NotImplementedError()

class QueueResultHandler(ResultHandler):
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
            timeout=self.timeout)

class ResultHandlerWrapper(ResultHandler):
    """Wraps a ResultHandler.
    """
    def __init__(self, handler):
        self.handler = handler
    
    def start(self, worker):
        self.handler.start(worker)
    
    def write_result(self, batch_num, result):
        self.handler.write_result(batch_num, result)
    
    def finish(self, total_batches=None):
        self.handler.finish(total_batches=total_batches)

def launch_workers(n, worker_class, args=(), offset=0):
    logging.getLogger().info("Starting {} worker processes".format(n))
    workers = [
        worker_class(i+offset, *args)
        for i in range(n)
    ]
    # start workers
    for worker in workers: worker.start()
    return workers

def ensure_processes(processes, message="One or more process exited: {}", alive=True):
    is_alive = [process.is_alive() for worker in processes]
    if alive != all(is_alive):
        raise Exception(message.format(",".join(
            str(i) for i in range(len(is_alive)) if not is_alive[i])))

def wait_on(condition, wait_message="Waiting {}", timeout=None, fail_callback=None,
            wait=None, timeout_callback=None):
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

def wait_on_process(process, timeout, terminate=False):
    timeout_callback = lambda: process.terminate() if terminate else None
    wait_on(
        lambda: not process.is_alive(),
        wait_message="Waiting on {} to terminate {{}}".format(process.name),
        timeout=timeout,
        wait=lambda: process.join(RETRY_INTERVAL),
        timeout_callback=timeout_callback)

def enqueue(queue, item, wait_message="Waiting to enqueue item {}", block_timeout=RETRY_INTERVAL, **kwargs):
    def condition():
        try:
            queue.put(item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False
    wait_on(condition, wait_message=wait_message, **kwargs)

def enqueue_all(iterable, queue, timeout, fail_callback):
    num_items = 0
    for item in iterable:
        def condition():
            try:
                queue.put(item, block=True, timeout=RETRY_INTERVAL)
                return True
            except Full:
                return False
        wait_on(
            condition,
            wait_message="Main process waiting to queue item {}",
            timeout=timeout,
            fail_callback=fail_callback)
        num_items += 1
    return num_items

def dequeue(queue, wait_message="Waiting to dequeue item {}", block_timeout=RETRY_INTERVAL, **kwargs):
    def condition():
        try:
            return queue.get(block=True, timeout=block_timeout)
        except Empty:
            return False
    return wait_on(condition, wait_message=wait_message, **kwargs)
