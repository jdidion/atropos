# cython: profile=False, emit_code_comments=False
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc

from cpython.array cimport array, clone
cdef array ld_array = array('d', [])

cdef class FactorialCache:
    cdef array factorials
    cdef int cur_array_size
    cdef int max_n
    cdef bint nan_limit
    
    def __cinit__(self, bint init_size=150):
        self.factorials = clone(ld_array, init_size, zero=True)
        self.factorials[0] = 1
        self.factorials[1] = 1
        self.max_n = 1
        self.cur_array_size = init_size
        self.nan_limit = False
        
    def factorial(self, n):
        if n > self.max_n:
            if self.nan_limit:
                return -1
            if n >= self.cur_array_size:
                self._extend(n)
            self._fill_upto(n)
            if self.nan_limit:
                return -1
        return self.factorials[n]
    
    def can_compute(self, n):
        if n <= self.max_n:
            return True
        if not self.nan_limit:
            self._fill_upto(n)
        return n <= self.max_n
    
    def _extend(self, n):
        cdef int extension_size = n - self.cur_array_size + 1
        cdef array extension = clone(ld_array, extension_size, zero=False)
        array.extend(self.factorials, extension)
        self.cur_array_size += extension_size
    
    def _fill_upto(self, n):
        cdef int i = self.max_n
        cdef int next_i = i + 1
        cdef double next_val
        while i < n:
            next_val = next_i * self.factorials[i]
            if next_val == float('inf'):
                self.nan_limit = True
                break
            else:
                self.factorials[next_i] = next_val
                i = next_i
                next_i += 1
        self.max_n = i

def match_probability(int matches, int mismatches, object factorial_cache):
    cdef int count = matches + mismatches
    cdef double p = 0.0
    cdef int i
    
    while not factorial_cache.can_compute(count):
        count = matches // 2
        mismatches = mismatches // 2
        count = matches + mismatches
    
    i = matches
    nfac = factorial_cache.factorial(count)
    while i <= count:
        p += (
            (0.75 ** (count - i)) *
            (0.25 ** i) *
            nfac /
            factorial_cache.factorial(i) /
            factorial_cache.factorial(count - i)
        )
        i += 1
    
    return p
