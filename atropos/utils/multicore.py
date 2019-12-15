from enum import IntEnum
from functools import partial
import inspect
from multiprocessing import Process, Value
from queue import Empty, Full, Queue
from threading import Lock
import time
from typing import Any, Callable, Iterable, Optional, Union, cast

from loguru import logger

from atropos.errors import MulticoreError
from atropos.utils import ReturnCode


DEFAULT_RETRY_INTERVAL = 5
DEFAULT_BLOCK_TIMEOUT = 5


class ControlSignal(IntEnum):
    ACTIVE = 0
    """Controlled process should run normally."""
    ERROR = -1
    """Controlled process should exit."""


class Control:
    """
    Shared (long) value for passing control information between main and worker threads.

    Args:
        initial_value: Initial value of the shared control variable.
    """

    def __init__(self, initial_value=ControlSignal.ACTIVE):
        self.control = Value('l', initial_value)

    def check_value(self, value, lock=False):
        """
        Checks that the current control value == `value`.

        Args:
            value: The value to check.
            lock: Whether to lock the shared variable before checking.

        Returns:
            True if the values are equal.
        """
        return self.get_value(lock=lock) == value

    def check_value_positive(self, lock=False):
        """
        Checks that the current control value is positive.

        Args:
            lock: Whether to lock the shared variable before checking.
        """
        return self.get_value(lock=lock) > 0

    def get_value(self, lock=True):
        """
        Returns the current control value.

        Args:
            lock: Whether to lock the shared variable before checking.
        """
        if lock:
            with self.control.get_lock():
                return self.control.value

        else:
            return self.control.value

    def set_value(self, value):
        """
        Sets the control value. The shared variable is always locked.

        Args:
            value: The value to set.
        """
        with self.control.get_lock():
            self.control.value = value


class PendingQueue:
    """
    Queue for items with sequentially increasing priority. An item whose
    priority is below the current level is queued. Pop returns the item with
    the current priority and increments the current priority.
    """

    def __init__(self, max_size=None):
        """
        Args:
            max_size: Maximum queue size; None == infinite.
        """
        self.queue = {}
        self.max_size = max_size
        self.min_priority = None

    def push(self, priority, value):
        """
        Adds an item to the queue with priority.

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
            raise ValueError(f"Duplicate priority value: {priority}")

        self.queue[priority] = value
        if self.min_priority is None or priority < self.min_priority:
            self.min_priority = priority

    def pop(self):
        """
        Removes and return the item in the queue with lowest priority.

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
        """
        Whether the queue is full.
        """
        return self.max_size and len(self.queue) >= self.max_size

    @property
    def empty(self):
        """
        Whether the queue is empty.
        """
        return len(self.queue) == 0


def synchronized(f):
    """Synchronization decorator."""

    f._lock = Lock()

    def wrapper(*args, **kw):
        f._lock.acquire()

        try:
            return f(*args, **kw)
        finally:
            f._lock.release()

    return wrapper


def ensure_processes(
    processes: Iterable[Process],
    message: str = "One or more process exited: {}",
    alive: bool = True
):
    """
    Raises an exception if all processes do not have the expected status (alive or
    dead).
    """
    is_alive = [worker.is_alive() for worker in processes]
    if alive != all(is_alive):
        raise MulticoreError(
            message.format(
                ",".join(str(i) for i, a in enumerate(is_alive) if a != alive)
            )
        )


def wait_on(
    condition: Callable,
    *args,
    wait_message: str = "Waiting {}",
    timeout: Optional[int] = None,
    fail_callback: Optional[Callable] = None,
    wait: Union[bool, int, Callable, None] = None,
    timeout_callback: Optional[Callable] = None,
    retry_interval: int = DEFAULT_RETRY_INTERVAL
) -> Any:
    """
    Waits on a condition to be non-False.

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
        retry_interval:

    Returns:

    """
    if wait is True:
        wait_fn = partial(time.sleep, retry_interval)
    elif isinstance(wait, int):
        wait_time = wait
        wait_fn = partial(time.sleep, wait_time)
    elif wait:
        wait_fn = cast(Callable, wait)
    else:
        wait_fn = None

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
            msg = wait_message.format(f"for {round(waiting, 1)} seconds")

            if timeout is not None and waiting >= timeout:
                logger.error(msg)
                if timeout_callback:
                    if inspect.isclass(timeout_callback):
                        raise timeout_callback()
                    else:
                        timeout_callback()
            else:
                logger.debug(msg)

            if wait_fn:
                wait_fn()


def wait_on_process(
    process: Process,
    wait_timeout: int,
    terminate: bool = False,
    retry_interval: int = DEFAULT_RETRY_INTERVAL
) -> Any:
    """
    Waits on a process to terminate.

    Args:
        process: The process on which to wait.
        wait_timeout: Number of seconds to wait for process to terminate.
        terminate: Whether to force the process to terminate after `timeout`
            seconds.
        retry_interval:
    """
    timeout_callback = process.terminate if terminate else None
    return wait_on(
        lambda: not process.is_alive(),
        wait_message=f"Waiting on {process.name} to terminate {{}}",
        timeout=wait_timeout,
        wait=lambda: process.join(retry_interval),
        timeout_callback=timeout_callback,
    )


def enqueue(
    queue: Queue,
    item: Any,
    wait_message: str = "Waiting to enqueue item {}",
    block_timeout: int = DEFAULT_BLOCK_TIMEOUT,
    **kwargs,
) -> None:
    """
    Enqueues an item, using `wait_on` to wait while `queue` is full.

    Args:
        queue: The queue to which to add the item.
        item: The item to queue.
        wait_message: The message to log while waiting.
        block_timeout: Number of seconds to wait after each `queue.put` attempt.
        kwargs: Additional arguments to `wait_on`.
    """

    def condition(_item):
        """
        Returns True if enqueing was successful.
        """
        try:
            queue.put(_item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False

    wait_on(condition, item, wait_message=wait_message, **kwargs)


def enqueue_all(
    iterable: Iterable,
    queue: Queue,
    wait_timeout: int,
    fail_callback: Callable,
    block_timeout: int = DEFAULT_BLOCK_TIMEOUT
) -> int:
    """
    Enqueues all items in `iterable`, using `wait_on` to wait while `queue` is full.

    Args:
        iterable: Iterable of items to queue.
        queue: The queue to which to add the item.
        wait_timeout: Number of seconds to wait after each `queue.put` attempt.
        fail_callback: Function called (or Exception raised) after timeout.
        block_timeout: Number of seconds to block while attempting to enqueue an item.

    Returns:
        The number of items queued.
    """
    num_items = 0

    def condition(_item):
        """Returns True if enqueing was successful.
        """
        try:
            queue.put(_item, block=True, timeout=block_timeout)
            return True
        except Full:
            return False

    for item in iterable:
        wait_on(
            condition,
            item,
            wait_message="Main process waiting to queue item {}",
            timeout=wait_timeout,
            fail_callback=fail_callback,
        )
        num_items += 1

    return num_items


def dequeue(
    queue: Queue,
    wait_message: str = "Waiting to dequeue item {}",
    block_timeout: int = DEFAULT_BLOCK_TIMEOUT,
    **kwargs,
) -> Any:
    """
    Dequeues an item, using `wait_on` to wait while `queue` is empty.
    """

    def condition():
        """
        Returns an item from the queue, or False if the queue is empty.
        """
        try:
            return queue.get(block=True, timeout=block_timeout)
        except Empty:
            return False

    return wait_on(condition, wait_message=wait_message, **kwargs)


def kill(process: Process, retcode: ReturnCode, wait_timeout: int):
    """
    Kills a process if it fails to terminate on its own.
    """
    if retcode <= 1:
        wait_on_process(process, wait_timeout, terminate=True)
    elif process.is_alive():
        process.terminate()
