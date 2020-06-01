"""Classes and methods to support parallelization of operations.
"""
import inspect
import logging
from multiprocessing import Process, Value, Queue
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
    """Base error for parallel processes.
    """
    pass

class Control(object):
    """Shared (long) value for passing control information between main and
    worker threads.

    Args:
        initial_value: Initial value of the shared control variable.
    """
    def __init__(self, initial_value=CONTROL_ACTIVE):
        self.control = Value('l', initial_value)

    def check_value(self, value, lock=False):
        """Check that the current control value == `value`.

        Args:
            value: The value to check.
            lock: Whether to lock the shared variable before checking.

        Returns:
            True if the values are equal.
        """
        return self.get_value(lock=lock) == value

    def check_value_positive(self, lock=False):
        """Check that the current control value is positive.

        Args:
            lock: Whether to lock the shared variable before checking.
        """
        return self.get_value(lock=lock) > 0

    def get_value(self, lock=True):
        """Returns the current control value.

        Args:
            lock: Whether to lock the shared variable before checking.
        """
        if lock:
            with self.control.get_lock():
                return self.control.value
        else:
            return self.control.value

    def set_value(self, value):
        """Set the control value. The shared variable is always locked.

        Args:
            value: The value to set.
        """
        with self.control.get_lock():
            self.control.value = value

class PendingQueue(object):
    """Queue for items with sequentially increasing priority. An item whose
    priority is below the current level is queued. Pop returns the item with
    the current priority and increments the current priority.

    Args:
        max_size: Maximum queue size; None == infinite.
    """
    def __init__(self, max_size=None):
        self.queue = {}
        self.max_size = max_size
        self.min_priority = None

    def push(self, priority, value):
        """Add an item to the queue with priority.

        Args:
            priority: An integer that determines the placement of `value` in
                the queue. Must be unique.
            value: The value to queue.

        Raises:
            Full if queue is full.
        """
        if self.full:
            raise Full()
        if priority in self.queue:
            raise ValueError("Duplicate priority value: {}".format(priority))
        self.queue[priority] = value
        if self.min_priority is None or priority < self.min_priority:
            self.min_priority = priority

    def pop(self):
        """Remove and return the item in the queue with lowest priority.

        Raises:
            Empty if the queue is emtpy.
        """
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
        """Whether the queue is full.
        """
        return self.max_size and len(self.queue) >= self.max_size

    @property
    def empty(self):
        """Whether the queue is empty.
        """
        return len(self.queue) == 0

class ParallelPipelineMixin(object):
    """Mixin that implements the `start`, `finish`, and `process_batch` methods
    of :class:`Pipeline`.
    """
    def start(self, **kwargs):
        super().start(**kwargs)
        self.seen_batches = set()

    def process_batch(self, batch):
        self.seen_batches.add(batch[0]['index'])
        super().process_batch(batch)

    def finish(self, summary, worker=None):
        super().finish(summary, worker=worker)
        logging.getLogger().debug(
            "%s finished; processed %d batches, %d reads",
            worker.name, len(self.seen_batches),
            sum(self.record_counts.values()))

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
            "%s running under pid %d", self.name, os.getpid())

        summary = {}

        def iter_batches():
            """Deque and yield batches.
            """
            while True:
                batch = dequeue(
                    self.input_queue,
                    wait_message="{} waiting on batch {{}}".format(self.name),
                    timeout=self.timeout)
                yield batch

        def enqueue_summary():
            """Enqueue a summary dict.
            """
            enqueue(
                self.summary_queue,
                (self.index, self.pipeline.seen_batches, summary),
                wait_message="{} waiting to queue summary {{}}".format(
                    self.name),
                timeout=self.timeout
            )

        try:
            self.pipeline.start(worker=self)

            try:
                for batch in iter_batches():
                    if batch is None:
                        break
                    logging.getLogger().debug(
                        "%s processing batch %d of size %d",
                        self.name, batch[0]['index'], batch[0]['size'])
                    self.pipeline.process_batch(batch)
            finally:
                self.pipeline.finish(summary, worker=self)

            logging.getLogger().debug("%s finished normally", self.name)
        except Exception as err:
            logging.getLogger().error(
                "Unexpected error in %s", self.name, exc_info=True)
            summary['exception'] = err

        logging.getLogger().debug("%s sending summary", self.name)
        enqueue_summary()

class ParallelPipelineRunner(object):
    """Run a pipeline in parallel.

    Args:
        reader: A :class:`BatchReader`.
        pipeline: A :class:`Pipeline`.
        threads: Number of threads to use. If None, the value will be taken
            from command_runner.
    """
    def __init__(self, command_runner, pipeline, threads=None):
        self.threads = threads or command_runner.threads
        if threads < 2:
            raise ValueError("'threads' must be >= 2")
        self.command_runner = command_runner
        self.pipeline = pipeline
        self.timeout = max(command_runner.process_timeout, RETRY_INTERVAL)
        # Queue by which batches of reads are sent to worker processes
        self.input_queue = Queue(command_runner.read_queue_size)
        # Queue for processes to send summary information back to main process
        self.summary_queue = Queue(self.threads)
        self.worker_processes = None
        self.num_batches = None
        self.seen_summaries = None
        self.seen_batches = None

    def ensure_alive(self):
        """Callback when enqueue times out.
        """
        ensure_processes(self.worker_processes)

    def after_enqueue(self):
        """Called after all batches are queued.
        """
        pass

    def finish(self):
        """Called at the end of the run.
        """
        pass

    def run(self):
        """Run the pipeline.

        Returns:
            The return code.
        """
        retcode = run_interruptible(self)
        self.terminate(retcode)
        return retcode

    def terminate(self, retcode):
        """Terminates all processes.

        Args:
            retcode: The return code.
        """
        # notify all threads that they should stop
        logging.getLogger().debug("Exiting all processes")
        for process in self.worker_processes:
            kill(process, retcode, self.timeout)

    def __call__(self):
        # Start worker processes, reserve a thread for the reader process,
        # which we will get back after it completes
        worker_args = (
            self.input_queue, self.pipeline, self.summary_queue, self.timeout)
        self.worker_processes = launch_workers(self.threads - 1, worker_args)

        self.num_batches = enqueue_all(
            self.command_runner.iterator(), self.input_queue, self.timeout,
            self.ensure_alive)

        logging.getLogger().debug(
            "Main loop complete; saw %d batches", self.num_batches)

        # Tell the worker processes no more input is coming
        enqueue_all(
            (None,) * self.threads, self.input_queue, self.timeout,
            self.ensure_alive)

        self.after_enqueue()

        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        self.worker_processes.extend(
            launch_workers(1, worker_args, offset=self.threads-1))

        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            """Ensure that workers are still alive.
            """
            try:
                ensure_processes(
                    self.worker_processes,
                    "Workers are still alive and haven't returned summaries: {}",
                    alive=False)
            except Exception as err:
                logging.getLogger().error(err)

        wait_on(
            self.summary_queue.full,
            wait_message="Waiting on worker summaries {}",
            timeout=self.timeout,
            wait=True,
            timeout_callback=summary_timeout_callback)

        # Process summary information from worker processes
        logging.getLogger().debug(
            "Processing summary information from worker processes")

        self.seen_summaries = set()
        self.seen_batches = set()

        def summary_fail_callback():
            """Raises AtroposError with workers that did not report summaries.
            """
            missing_summaries = (
                set(range(1, self.threads)) - self.seen_summaries)
            raise AtroposError(
                "Missing summaries from processes %s",
                ",".join(str(summ) for summ in missing_summaries))

        for _ in range(1, self.threads+1):
            batch = dequeue(
                self.summary_queue, fail_callback=summary_fail_callback)
            worker_index, worker_batches, worker_summary = batch
            if worker_summary is None:
                raise MulticoreError(
                    "Worker process {} died unexpectedly".format(worker_index))
            elif (
                    'exception' in worker_summary and
                    worker_summary['exception'] is not None):
                raise AtroposError(
                    "Worker process {} died unexpectedly".format(worker_index),
                    worker_summary['exception'])
            else:
                logging.getLogger().debug(
                    "Processing summary for worker %d", worker_index)
            self.seen_summaries.add(worker_index)
            self.seen_batches |= worker_batches
            self.command_runner.summary.merge(worker_summary)

        # Check if any batches were missed
        if self.num_batches > 0:
            missing_batches = (
                set(range(1, self.num_batches+1)) - self.seen_batches)
            if len(missing_batches) > 0:
                raise AtroposError(
                    "Workers did not process batches {}".format(
                        ",".join(str(batch) for batch in missing_batches)))

        self.finish()

def launch_workers(num_workers, args=(), offset=0, worker_class=WorkerProcess):
    """Launch `n` workers. Each worker is initialized with an incremental
    index starting with `offset`, followed by `args`.
    """
    logging.getLogger().info("Starting %d worker processes", num_workers)
    # create workers
    workers = [worker_class(i+offset, *args) for i in range(num_workers)]
    # start workers
    for worker in workers:
        worker.start()
    return workers

def ensure_processes(
        processes, message="One or more process exited: {}", alive=True):
    """Raise an exception if all processes do not have the expected status
    (alive or dead).
    """
    is_alive = [worker.is_alive() for worker in processes]
    if alive != all(is_alive):
        raise MulticoreError(message.format(",".join(
            str(i) for i, a in enumerate(is_alive) if a != alive)))

def wait_on(
        condition, *args, wait_message="Waiting {}", timeout=None,
        fail_callback=None, wait=None, timeout_callback=None):
    """Wait on a condition to be non-False.

    Args:
        condition: Function that returns either False or a non-False value.
        args: Args with which to call `condition`.
        wait_message: The message to log while `condition` is False.
        timeout: Number of seconds after which the log messages escalate from
            DEBUG to ERROR.
        fail_callback: Function that is called each time `condition` returns
            False.
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
        result = condition(*args)
        if result is not False:
            return result
        if fail_callback:
            fail_callback()
        now = time.time()
        if not wait_start:
            wait_start = now
        else:
            waiting = now - wait_start
            msg = wait_message.format(
                "for {} seconds".format(round(waiting, 1)))
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
        terminate: Whether to force the process to terminate after `timeout`
            seconds.
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
    def condition(item):
        """Returns True if enqueing was successful.
        """
        try:
            queue.put(item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False
    wait_on(condition, item, wait_message=wait_message, **kwargs)

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
    def condition(item):
        """Returns True if enqueing was successful.
        """
        try:
            queue.put(item, block=True, timeout=RETRY_INTERVAL)
            return True
        except Full:
            return False
    for item in iterable:
        wait_on(
            condition, item,
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
        """Returns an item from the queue, or False if the queue is empty.
        """
        try:
            return queue.get(block=True, timeout=block_timeout)
        except Empty:
            return False
    return wait_on(condition, wait_message=wait_message, **kwargs)

def kill(process, retcode, timeout):
    """Kill a process if it fails to terminate on its own.
    """
    if retcode <= 1:
        wait_on_process(process, timeout, terminate=True)
    elif process.is_alive():
        process.terminate()
