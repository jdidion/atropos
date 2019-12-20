import os
import tempfile

from atropos.adapters import Adapter, AdapterType, ColorspaceAdapter
from atropos.commands.trim.modifiers import AdapterCutter
from atropos.commands.trim.multicore import OrderPreservingWriterResultHandler
from atropos.commands.trim.writers import Writers
from atropos.io.sequence import ColorspaceSequence, Sequence


def test_cs_5p():
    read = ColorspaceSequence("name", "0123", "DEFG", "T")
    adapter = ColorspaceAdapter("CG", AdapterType.PREFIX, 0.1)
    cutter = AdapterCutter([adapter])
    cutter(read)


# no assertion here, just make sure the above code runs without
# an exception
def test_metrics():
    read = Sequence("name", "AAAACCCCAAAA")
    adapters = [Adapter("CCCC", AdapterType.BACK, 0.1)]
    cutter = AdapterCutter(adapters, times=3)
    cutter(read)
    # TODO make this a lot simpler
    trimmed_bp = 0
    for adapter in adapters:
        for d in (adapter.lengths_front, adapter.lengths_back):
            trimmed_bp += sum(seqlen * count for (seqlen, count) in d.items())
    assert trimmed_bp <= len(read), trimmed_bp


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
        with open(path, "rt") as inp:
            assert inp.read() == (result1 + result2 + result3)
    finally:
        os.remove(path)
