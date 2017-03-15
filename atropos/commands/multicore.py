"""Classes and methods to support parallel operations.
"""
from collections import defaultdict
import inspect
import logging
from multiprocessing import Process, Queue, Value, cpu_count
import os
from queue import Empty, Full
import time
from atropos import AtroposError
from atropos.util import run_interruptible

RETRY_INTERVAL = 5
"""Max time to wait between retrying operations."""

# Control values
CONTROL_ACTIVE = 0
"""Controlled process should run normally."""
CONTROL_ERROR = -1
"""Controlled process should exit."""

class MulticoreError(AtroposError):
    pass

class Control(object):
    """Shared (long) value for passing control information between main and
    worker threads.
    """
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
    """Queue for items with sequentially increasing priority. An item whose
    priority is below the current level is queued. Pop returns the item with
    the current priority and increments the current priority.
    
    Args:
        max_size: Maximum queue size; None == infinite.
    """
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
    """Parent class for worker processes that execute Pipelines.
    
    Args:
        index: A unique ID for the process.
        input_queue: Queue with batches of records to process.
        pipeline: The pipeline to execute.
        summary_queue: Queue where summary information is written.
        timeout: Time to wait upon queue full/empty.
    """
    def __init__(self, index, input_queue, pipeline, summary_queue, timeout):
        super().__init__(name="Worker process {}".format(index))
        self.index = index
        self.input_queue = input_queue
        self.pipeline = pipeline
        self.summary_queue = summary_queue
        self.timeout = timeout
    
    def run(self):
        logging.getLogger().debug(
            "{} running under pid {}".format(self.name, os.getpid()))
        
        def iter_batches():
            while True:
                batch = dequeue(
                    self.input_queue,
                    wait_message="{} waiting on batch {{}}".format(self.name),
                    timeout=self.timeout)
                yield batch
        
        def enqueue_summary(error=None):
            seen_batches, summary = self.pipeline.summarize(error)
            enqueue(
                self.summary_queue, (self.index, seen_batches, summary),
                wait_message="{} waiting to queue summary {{}}".format(self.name),
                timeout=self.timeout
            )
        
        try:
            self.pipeline.start(worker=self)
            
            try:
                for batch in iter_batches():
                    if batch is None:
                        break
                    logging.getLogger().debug(
                        "{} processing batch of size {}".format(
                            self.name, batch[0]['size']))
                    self.pipeline.process_batch(batch)
            finally:
                self.pipeline.finish(worker=self)
            
            logging.getLogger().debug("{} sending summary".format(self.name))
            enqueue_summary()
            
            logging.getLogger().debug("{} exiting normally".format(self.name))
        except Exception as e:
            logging.getLogger().error(
                "Unexpected error in {}".format(self.name), exc_info=True)
            enqueue_summary(error=e)

def launch_workers(n, worker_class=WorkerProcess, args=(), offset=0):
    """Launch `n` workers. Each worker is initialized with an incremental
    index starting with `offset`, followed by `args`.
    """
    logging.getLogger().info("Starting {} worker processes".format(n))
    # create workers
    workers = [worker_class(i+offset, *args) for i in range(n)]
    # start workers
    for worker in workers:
        worker.start()
    return workers

def ensure_processes(processes, message="One or more process exited: {}", alive=True):
    """Raise an exception if all processes do not have the expected status
    (alive or dead).
    """
    is_alive = [process.is_alive() for worker in processes]
    if alive != all(is_alive):
        raise MulticoreError(message.format(",".join(
            str(i) for i in range(len(is_alive)) if not is_alive[i])))

def wait_on(
        condition, wait_message="Waiting {}", timeout=None, fail_callback=None,
        wait=None, timeout_callback=None):
    """Wait on a condition to be non-False.
    
    Args:
        condition: Function that returns either False or a non-False value.
        wait_message: The message to log while `condition` is False.
        timeout: Number of seconds after which the log messages escalate from
            DEBUG to ERROR.
        fail_callback: Function that is called each time `condition` returns False.
        wait: Either a boolean or a wait function to execute. If True,
            `time.sleep(RETRY_INTERVAL)` is used as the wait function.
        timeout_callback: Eather a function that is called each time, or an
            Exception class that is raised the first time `condition`
            returns False after waiting for `timeout` seconds.
    """
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
    """Wait on a process to terminate.
    
    Args:
        process: The process on which to wait.
        timeout: Number of seconds to wait for process to terminate.
        terminate: Whether to force the process to terminate after `timeout` seconds.
    """
    timeout_callback = lambda: process.terminate() if terminate else None
    return wait_on(
        lambda: not process.is_alive(),
        wait_message="Waiting on {} to terminate {{}}".format(process.name),
        timeout=timeout,
        wait=lambda: process.join(RETRY_INTERVAL),
        timeout_callback=timeout_callback)

def enqueue(
        queue, item, wait_message="Waiting to enqueue item {}",
        block_timeout=RETRY_INTERVAL, **kwargs):
    """Enqueue an item, using `wait_on` to wait while `queue` is full.
    
    Args:
        queue: The queue to which to add the item.
        item: The item to queue.
        wait_message: The message to log while waiting.
        block_timeout: Number of seconds to wait after each `queue.put` attempt.
        kwargs: Additional arguments to `wait_on`.
    """
    def condition():
        try:
            queue.put(item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False
    wait_on(condition, wait_message=wait_message, **kwargs)

def enqueue_all(iterable, queue, timeout, fail_callback):
    """Enqueue all items in `iterable`, using `wait_on` to wait while `queue`
    is full.
    
    Args:
        iterable: Iterable of items to queue.
        queue: The queue to which to add the item.
        timeout: Number of seconds to wait after each `queue.put` attempt.
        fail_callback: Function called (or Exception raised) after timeout.
    
    Returns:
        The number of items queued.
    """
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

def dequeue(
        queue, wait_message="Waiting to dequeue item {}",
        block_timeout=RETRY_INTERVAL, **kwargs):
    """Dequeue an item, using `wait_on` to wait while `queue` is empty.
    """
    def condition():
        try:
            return queue.get(block=True, timeout=block_timeout)
        except Empty:
            return False
    return wait_on(condition, wait_message=wait_message, **kwargs)
