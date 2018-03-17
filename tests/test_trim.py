# coding: utf-8
from atropos.adapters import Adapter, ColorspaceAdapter, PREFIX, BACK
from atropos.commands.trim.modifiers import AdapterCutter
from atropos.io.seqio import ColorspaceSequence, Sequence


def test_cs_5p():
    read = ColorspaceSequence("name", "0123", "DEFG", "T")
    adapter = ColorspaceAdapter("CG", PREFIX, 0.1)
    cutter = AdapterCutter([adapter])
    trimmed_read = cutter(read)




# no assertion here, just make sure the above code runs without
# an exception
def test_statistics():
    read = Sequence('name', 'AAAACCCCAAAA')
    adapters = [Adapter('CCCC', BACK, 0.1)]
    cutter = AdapterCutter(adapters, times=3)
    trimmed_read = cutter(read)
    # TODO make this a lot simpler
    trimmed_bp = 0
    for adapter in adapters:
        for d in (adapter.lengths_front, adapter.lengths_back):
            trimmed_bp += sum(seqlen * count for (seqlen, count) in d.items())
    assert trimmed_bp <= len(read), trimmed_bp
