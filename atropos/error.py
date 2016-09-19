# Estimate the empircal error rate

from .seqio import FastqReader
from .util import enumerate_range

def qual2prob(qchar):
    q = ord(qchar) - 33
    return 10 ** (-q / 10)

class ErrorEstimator(object):
    def __init__(self):
        self.total_qual = 0.0
        self.total_len = 0
            
    def consume(self, read):
        self.total_qual += sum(qual2prob(qchar) for qchar in read.qualities)
        self.toal_len += len(read.qualities)
    
    def estimate(self):
        return self.total_qual / self.total_len
