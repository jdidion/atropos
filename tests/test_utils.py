from atropos import AtroposError
from atropos.util import *
import time
from unittest import TestCase

class UtilTests(TestCase):
    def test_factorial(self):
        f = RandomMatchProbability()
        # simple test
        assert f.factorial(0) == 1
        assert f.factorial(1) == 1
        assert f.factorial(3) == 6
        assert int(f.factorial(27)) == int(math.factorial(27))
        # test big number
        assert int(f.factorial(150)) == int(math.factorial(150))
        # test really big number
        assert f.factorial(10000) is not None
        
    def test_timing(self):
        t = Timing()
        with self.assertRaises(AtroposError):
            t.summarize()
        
        with Timing() as t:
            time.sleep(1)
        summary = t.summarize()
        assert summary['start'] == t.start_time.isoformat()
        assert summary['wallclock'] >= 1
    
    def test_CountingDict(self):
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
        
        with self.assertRaises(ValueError):
            cd.merge({})

    def test_Histogram(self):
        h = Histogram(keys=(1, 1, 1, 2, 2, 3), sort_by=1)
        s = h.summarize()
        assert list(s['hist'].keys()) == [3, 2, 1]
        assert s['summary']['mean'] == 10 / 6
        assert s['summary']['median'] == 1.5
        assert s['summary']['modes'] == [1]

    def test_NestedDict(self):
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
        assert set(nd1.summarize()) == set((
            ('a', 'x', 2),
            ('a', 'y', 1),
            ('a', 'z', 1),
            ('b', 'y', 1)
        ))

    def test_MergingDict(self):
        obj = object()
        
        md = MergingDict()
        md['a'] = Const(1)
        md['b'] = dict(a=1, b=2)
        md['c'] = 'foo'
        md['d'] = 5
        md['e'] = (1,2,3)
        md['f'] = True
        md['g'] = obj
        
        md2 = MergingDict()
        md2['a'] = 1
        md2['b'] = dict(a=3, b=4)
        md2['c'] = 'foo'
        md2['d'] = 7
        md2['e'] = (4,5,6)
        md2['f'] = True
        md2['g'] = obj
        
        md.merge(md2)
        
        assert md['a'] == 1
        assert md['b'] == dict(a=4, b=6)
        assert md['c'] == 'foo'
        assert md['d'] == 12
        assert md['e'] == [5, 7, 9]
        assert md['f'] == True
        
        failures = [
            ('a', 2, ValueError),
            ('b', 'foo', TypeError),
            ('c', 'bar', ValueError),
            ('f', False, ValueError),
            ('g', object(), ValueError)
        ]
        
        for k, v, err_type in failures:
            md3 = md2.copy()
            md3[k] = v
            with self.assertRaises(err_type):
                md.merge(md3)
    
    def test_complement(self):
        assert complement('ACCTGCCA') == 'TGGACGGT'
    
    def test_reverse_complement(self):
        assert reverse_complement('ACCTGCCA') == 'TGGCAGGT'
    
    def test_sequence_complexity(self):
        assert sequence_complexity('AAAA') == 0
        assert sequence_complexity('ACGT') == 2
    
    def test_qual2int(self):
        assert qual2int('H') == 39
    
    def test_quals2ints(self):
        assert tuple(quals2ints('HIJ')) == (39,40,41)
    
    def test_qual2prob(self):
        assert qual2prob('+') == 0.1
    
    def test_enumerate_range(self):
        assert list(enumerate_range((1,2,3), 5, 7)) == [(5, 1), (6, 2)]
    
    def test_mean(self):
        with self.assertRaises(ValueError):
            mean([])
        assert mean(range(10)) == 4.5
    
    def test_weighted_mean(self):
        with self.assertRaises(ValueError):
            weighted_mean([], [])
        with self.assertRaises(ValueError):
            weighted_mean([1,2], [1,2,3])
        
    
    def test_stdev(self):
        with self.assertRaises(ValueError):
            stdev([])
        assert stdev((1,)) == 0
        assert stdev((1,2,3)) == 1
        assert stdev((1,3,5)) == 2
    
    def test_weighted_stdev(self):
        with self.assertRaises(ValueError):
            weighted_stdev([], [])
        with self.assertRaises(ValueError):
            weighted_stdev([1,2], [1,2,3])
        assert weighted_stdev([1], [10]) == 0
        assert round(weighted_stdev([1,2,3], [1,2,3]), 4) == 0.8165
    
    def test_median(self):
        with self.assertRaises(ValueError):
            median([])
        assert median((1,2,3)) == 2
        assert median((1,2,3,4)) == 2.5
    
    def test_weighted_median(self):
        with self.assertRaises(ValueError):
            weighted_median([], [])
        with self.assertRaises(ValueError):
            weighted_median([1,2], [1,2,3])
        assert weighted_median([1], [0]) == None
        assert weighted_median([1,2,3],[2,2,2]) == 2
        assert weighted_median([1,2,3,4],[2,2,2,2]) == 2.5
        assert weighted_median([1,2,3],[2,1,2]) == 2
        assert weighted_median([1,2,3],[1,2,3]) == 2.5
    
    def test_modes(self):
        with self.assertRaises(ValueError):
            modes([])
        assert modes([1]) == [1]
        assert modes([1,2]) == [1,2]
        assert modes([1,1,2]) == [1]
    
    def test_weighted_modes(self):
        with self.assertRaises(ValueError):
            weighted_modes([], [])
        with self.assertRaises(ValueError):
            weighted_modes([1,2], [1,2,3])
        assert weighted_modes([1], [1]) == [1]
        assert weighted_modes([1,2,3], [1,2,3]) == [3]
        assert weighted_modes([1,2,3], [1,1,1]) == [1,2,3]
    
    def test_truncate_string(self):
        assert truncate_string(None) is None
        assert truncate_string('a' * 100) == ('a' * 100)
        assert truncate_string('a' * 101) == ('a' * 97) + '...'
    
    def test_run_interruptible(self):
        def raise_interrupted(x=False):
            if x is True:
                raise KeyboardInterrupt()
        assert run_interruptible(raise_interrupted, False) == 0
        assert run_interruptible(raise_interrupted, True) == 130
        
        import errno
        def raise_io_error(x=False, y=errno.EPIPE):
            if x is True:
                err = IOError()
                err.errno = y
                raise err
        assert run_interruptible(raise_io_error, False) == 0
        assert run_interruptible(raise_io_error, True) == 1
        with self.assertRaises(IOError):
            run_interruptible(raise_io_error, True, 5)
        
        def raises_exception(x):
            if x is True:
                raise Exception()
        assert run_interruptible(raises_exception, False) == 0
        assert run_interruptible(raises_exception, True) == 1
