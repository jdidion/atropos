from atropos.util import *
from unittest import TestCase

def UtilTests(TestCase):
    def test_CountingDict():
        cd = CountingDict(keys=('a', 'b'), sort_by=1)
        assert cd['a'] == cd['b'] == 1
        cd.increment('a')
        assert cd['a'] == 2
        assert cd['c'] == 0
        merged = cd.merge(CountingDict(keys=('b', 'c')))
        assert merged['a'] == 2
        assert merged['b'] == 2
        assert merged['c'] == 1
        merged.increment('b')
        assert list(merged.summarize().keys()) == ['c', 'a', 'b']

    def test_Histogram():
        h = Histogram(keys=(1, 1, 1, 2, 2, 3), sort_by=1)
        s = h.summarize()
        assert list(s['hist'].keys()) == [3, 2, 1]
        assert s['summary']['mean'] == 10 / 6
        assert s['summary']['median'] == 1.5
        assert s['summary']['modes'] == [1]

    def test_NestedDict():
        nd1 = NestedDict()
        nd1['a'].increment('x')
        nd1['a'].increment('y')
        
        with self.assertRaises(ValueError):
            nd1.merge(1)
        
        nd2 = NestedDict()
        nd2['a'].increment('x')
        nd2['a'].increment('z')
        nd2['b'].increment('y')
        
        nd1.merge(nd2)
        
        assert nd1['a']['x'] == 2
        assert nd1['a']['y'] == 1
        assert nd1['a']['z'] == 1
        assert nd1['b']['y'] == 1
        
        assert dict(nd1.summarize()) == dict(
            columns=('x', 'y', 'z'),
            rows=dict(
                a=(2,1,1),
                b=(0,1,0)
            ))
        
        nd1.shape = 'long'
        assert nd1.summarize() == (
            ('a', 'x', 2),
            ('a', 'y', 1),
            ('a', 'z', 1),
            ('b', 'y', 1)
        )

    def test_MergingDict():
        md = MergingDict()
        md['a'] = Const(1)
        md['b'] = dict(a=1, b=2)
        md['c'] = 'foo'
        md['d'] = 5
        md['e'] = (1,2,3)
        md['f'] = True
        
        md2 = MergingDict()
        md2['a'] = 1
        md2['b'] = dict(a=3, b=4)
        md2['c'] = 'foo'
        md2['d'] = 7
        md2['e'] = (4,5,6)
        md2['f'] = True
        
        md.merge(md2)
        
        assert md['a'] == 1
        assert md['b'] == dict(a=4, b=6)
        assert md['c'] == 'foo'
        assert md['d'] == 12
        assert md['e'] == [5, 7, 9]
        assert md['f'] == True
        
        md2['f'] = False
        with self.assertRaises():
            md.merge(md2)
    
    def test_