# Tests for internal components of the atropos commands
from pytest import raises
from atropos.commands.trim.multicore import OrderPreservingWriterResultHandler
from atropos.commands.trim.writers import Writers
import tempfile
import os


def test_order_preserving_writer():
    path = tempfile.mkstemp()[1]
    try:
        writers = Writers()
        handler = OrderPreservingWriterResultHandler(writers)
        handler.start(None)
        # write three batches out of order
        result2 = "result2"
        handler.write_result(2, {path: result2})
        result3 = "result3"
        handler.write_result(3, {path: result3})
        result1 = "result1"
        handler.write_result(1, {path: result1})
        handler.finish(total_batches=3)
        # check that the results are in the right order
        with open(path, 'rt') as inp:
            assert inp.read() == (result1 + result2 + result3)
    finally:
        os.remove(path)
