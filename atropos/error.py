# Estimate the empircal error rate

from .seqio import FastqReader
from .util import enumerate_range

def main(options):
    n_reads = options.max_reads or 10000
    
    infiles = [f for f in (
        options.single_input, options.input1, options.input2, options.interleaved_input
    ) if f is not None]
    
    def qual2prob(qchar):
        q = ord(qchar) - 33
        return 10**(-q/10)
    
    total_qual = 0.0
    total_len = 0
    
    for i, f in enumerate(infiles, 1):
        fq = FastqReader(f)
        fq_iter = iter(fq)
        try:
            if options.progress:
                from .progress import create_progress_reader
                fq_iter = create_progress_reader(fq_iter, max_items=n_reads,
                    counter_magnitude="K", values_have_size=False)
            
            print("File {}: {}".format(i, f))
            
            cum_qual = 0.0
            cum_len = 0
            for i, read in enumerate_range(fq_iter, 0, n_reads):
                cum_qual += sum(qual2prob(qchar) for qchar in read.qualities)
                cum_len += len(read.qualities)
            
            print("  Error rate: {:.2%}".format(cum_qual / cum_len))
            total_qual += cum_qual
            total_len += cum_len
        finally:
            fq.close()
        
    print("Overall error rate: {:.2%}".format(total_qual / total_len))
