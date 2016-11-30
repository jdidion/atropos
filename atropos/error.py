# Estimate the empircal error rate

from atropos.seqio import FastqReader
from atropos.util import enumerate_range

def qual2prob(qchar):
    q = ord(qchar) - 33
    return 10 ** (-q / 10)

class ErrorEstimator(object):
    def __init__(self, max_read_len=None):
        self.total_qual = 0.0
        self.total_len = 0
        self.max_read_len = max_read_len
    
    def consume_all_batches(self, batch_iterator):
        for batch_num, batch in batch_iterator:
            for read in batch:
                self.consume(read)
    
    def consume(self, read):
        quals = read.qualities
        readlen = len(read.qualities)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            quals = quals[:readlen]
        self.total_qual += sum(qual2prob(qchar) for qchar in quals)
        self.total_len += readlen
    
    def estimate(self):
        return self.total_qual / self.total_len
    
    def summarize(self, outstream, name=None):
        print("", file=outstream)
        if name is not None:
            header = "File: {}".format(name)
            print(header)
            print('-' * len(header), file=outstream)
        print("Error rate: {:.2%}".format(self.estimate()))

class PairedErrorEstimator(object):
    def __init__(self, max_read_len=None):
        self.e1 = ErrorEstimator(max_read_len)
        self.e2 = ErrorEstimator(max_read_len)
    
    def consume_all_batches(self, batch_iterator):
        for batch_num, batch in batch_iterator:
            for read1, read2 in batch:
                self.e1.consume(read1)
                self.e2.consume(read2)
    
    def estimate(self):
        return (self.e1.estimate(), self.e2.estimate())
    
    def summarize(self, outstream, names=None):
        print("", file=outstream)
        if names:
            header = "File 1: {}".format(names[0])
            print(header)
            print('-' * len(header), file=outstream)
        print("Error rate: {:.2%}\n".format(self.e1.estimate()))
        
        if names:
            header = "File 2: {}".format(names[1])
            print(header)
            print('-' * len(header), file=outstream)
        print("Error rate: {:.2%}\n".format(self.e2.estimate()))
        
        print("Overall error rate: {:.2%}".format(
            (self.e1.total_qual + self.e2.total_qual) / (self.e1.total_len + self.e2.total_len)))
