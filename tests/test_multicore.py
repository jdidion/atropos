# coding: utf-8
from pytest import raises
import io
from multiprocessing import Process, Queue
import time
from atropos.commands.multicore import *


class TimeoutException(Exception):
    pass


log_capture_string = None


def setup():
    ### Create the logger
    logger = logging.getLogger('basic_logger')
    logger.setLevel(logging.DEBUG)
    ### Setup the console handler with a StringIO object
    log_capture_string = io.StringIO()
    ch = logging.StreamHandler(log_capture_string)
    ch.setLevel(logging.DEBUG)
    ### Optionally add a formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    ch.setFormatter(formatter)
    ### Add the console handler to the logger
    logger.addHandler(ch)


def assert_log(expected):
    log_contents = log_capture_string.getvalue()
    log_capture_string.close()
    assert log_contents == expected


def test_wait_on():

    class callbacks():

        def __init__(self):
            self.i = 0
            self.j = 0

        def condition(self):
            self.i += 1
            return self.i >= 5

        def fail_callback(self):
            self.j += 1

    c = callbacks()
    wait_on(c.condition, wait_message="waiting", fail_callback=c.fail_callback)
    assert c.i == 5
    assert c.j == 4




# assert_log("DEBUG: waiting for {} seconds".format())
def test_timeout():
    with raises(TimeoutException):

        def condition():
            return False

        wait_on(condition, timeout=2, wait=1, timeout_callback=TimeoutException)


def test_enqueue_dequeue():
    q = Queue(1)
    enqueue(q, 1)
    assert dequeue(q) == 1


def test_enqueue_timeout():
    with raises(TimeoutException):
        q = Queue(1)
        q.put(1)
        enqueue(q, 2, timeout=1, block_timeout=2, timeout_callback=TimeoutException)


def test_dequeue_timeout():
    with raises(TimeoutException):
        dequeue(Queue(1), timeout=1, block_timeout=2, timeout_callback=TimeoutException)


# TODO: port tests from testparallel here
# Test worker vs writer compression
# Test without writer process
