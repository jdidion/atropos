# Estimate the empircal error rate

from collections import defaultdict
import csv
import re
from atropos.seqio import FastqReader
from atropos.util import enumerate_range

class ErrorEstimator(object):
    def consume_all_batches(self, batch_iterator):
        for batch_num, batch in batch_iterator:
            for read in batch:
                self.consume(read)
    
    def summarize(self, outstream, name=None, show_details=True, indent="  "):
        print("", file=outstream)
        if name is not None:
            header = "File: {}".format(name)
            print(header, file=outstream)
            print('-' * len(header), file=outstream)
        overall, details = self.estimate()
        print("Error rate: {:.2%}".format(overall), file=outstream)
        if details and show_details:
            print("Details:\n", file=outstream)
            self.print_details(details, outstream, indent)
    
    def print_details(self, details, outstream, indent):
        pass

# Error estimation using base qualities

def qual2prob(qchar):
    q = ord(qchar) - 33
    return 10 ** (-q / 10)

class BaseQualityErrorEstimator(ErrorEstimator):
    """Simple error estimation using base qualities. It is well-known that base
    qualities significantly overestimate true error rates, so take the estimate
    with a grain of salt.
    """
    def __init__(self, max_read_len=None):
        self.total_qual = 0.0
        self.total_len = 0
        self.max_read_len = max_read_len
    
    def consume(self, read):
        quals = read.qualities
        readlen = len(quals)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            quals = quals[:readlen]
        self.total_qual += sum(qual2prob(qchar) for qchar in quals)
        self.total_len += readlen
    
    def estimate(self):
        return (self.total_qual / self.total_len, None)

# Error estimation using shadow counts

FILTER_RE = re.compile("A+|C+|G+|T+|.*N.*")

class ShadowRegressionErrorEstimator(ErrorEstimator):
    """Re-implementation of the shadow regression method described in:
    Wang et al., "Estimation of sequencing error rates in short reads",
        BMC Bioinformatics 2012 13:185, DOI: 10.1186/1471-2105-13-185
    
    Args:
        method: The differences that are considered in the error rate
            calculation; sub = substitutions, indel = insertions and
            deletions; all = both substitutions and indels.
        max_read_len: The maximum number of bases (starting from the 5' end)
            to consider from each read.
    """
    def __init__(self, method='sub', max_read_len=None, rscript_exe="Rscript"):
        self.seqs = defaultdict(lambda: 0)
        self.total_len = 0
        self.method = method
        self.max_read_len = max_read_len
        self.rscript_exe = rscript_exe
    
    def consume(self, read):
        seq = read.sequence
        readlen = len(seq)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            seq = seq[:readlen]
        if FILTER_RE.fullmatch(seq):
            return
        self.seqs[seq] += 1
        self.total_len += readlen
    
    def write_read_counts(self, fileobj):
        writer = csv.writer(fileobj, delimiter=" ")
        writer.writerows(sorted(
            self.seqs.items(), reverse=True, key=lambda i: i[1]))
    
    def estimate(self):
        # This is a temporary solution that requires R and the
        # ShadowRegression package. Eventually, this will be
        # replaced with a pure-python implementation.
        import os
        import subprocess
        import tempfile
        
        tempfiles = (tempfile.mkstemp()[1] for i in range(4))
        read_counts, per_read, per_cycle, script_file = tempfiles
        try:
            # write counts to a file
            with open(read_counts, 'wt') as o:
                self.write_read_counts(o)
            
            # execute R script
            script = """
            library(ShadowRegression)
            errorRates = getErrorRates("{reads}", type="{method}")
            write.table(errorRates$perReadER, "{per_read}", sep="\t", quote=F, col.names=F, row.names=T)
            write.table(errorRates$cycleER, "{per_cycle}", sep="\t", quote=F, col.names=F, row.names=T)
            """.format(
                reads=read_counts,
                method=self.method,
                per_read=per_read,
                per_cycle=per_cycle)
            with open(script_file, 'wt') as o:
                o.write(script)
            with subprocess.Popen(
                    [self.rscript_exe, "--vanilla", script_file],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
                stdout, stderr = p.communicate()
                if p.returncode != 0:
                    raise Exception(
                        "R script failed: rc={}; stdout={}; stderr={}".format(
                        p.returncode, stdout, stderr))
            
            # read the results
            with open(per_read, 'rt') as i:
                reader = csv.reader(i, delimiter="\t")
                per_read_error = dict(reader)
                if len(per_read_error) != 4:
                    raise Exception("Invalid output from R script")
            with open(per_cycle, 'rt') as i:
                reader = csv.reader(i, delimiter="\t")
                per_cycle_error = list(row[0:3] for row in reader)
                if not per_cycle_error:
                    raise Exception("Invalid output from R script")
            
            return (per_read_error["error rate"], dict(
                per_read=per_read_error,
                per_cycle=per_cycle_error))
        finally:
            for f in tempfiles:
                os.remove(f)
    
    def print_details(self, outstream, details, indent):
        per_read = details['per_read']
        per_cycle = details['per_cycle']
        
        print("{}StdErr: {:.2%}".format(indent, per_read['standard error']), file=outstream)
        print("{}Per-cycle rates:".format(indent), file=outstream)
        for cycle in per_cycle:
            print("{}Cycle: {}, Error: {:.2%}, StdErr: {:.2%}".format(indent*2, *cycle))

# Estimator for a pair of input files

class PairedErrorEstimator(object):
    def __init__(self, max_read_len=None,
                 estimator_class=BaseQualityErrorEstimator):
        self.e1 = estimator_class(max_read_len=max_read_len)
        self.e2 = estimator_class(max_read_len=max_read_len)
    
    def estimate(self):
        return (self.e1.estimate(), self.e2.estimate())
    
    def consume_all_batches(self, batch_iterator):
        for batch_num, batch in batch_iterator:
            for read1, read2 in batch:
                self.e1.consume(read1)
                self.e2.consume(read2)
    
    def summarize(self, outstream, names=None, show_details=True, indent="  "):
        print("", file=outstream)
        if names:
            header = "File 1: {}".format(names[0])
            print(header, file=outstream)
            print('-' * len(header), file=outstream)
        err1, details1 = self.e1.estimate()
        print("Error rate: {:.2%}".format(err1), file=outstream)
        if show_details and details1:
            print("Details:\n", file=outstream)
            self.e1.print_details(details1, outstream, indent)
        print("", file=outstream)
        
        if names:
            header = "File 2: {}".format(names[1])
            print(header, file=outstream)
            print('-' * len(header), file=outstream)
        err2, details2 = self.e2.estimate()
        print("Error rate: {:.2%}".format(err2), file=outstream)
        if show_details and details2:
            print("Details:\n", file=outstream)
            self.e2.print_details(details2, outstream, indent)
        print("", file=outstream)
        
        l1 = self.e1.total_len
        l2 = self.e2.total_len
        overall_err = ((err1 * l1) + (err2 * l2)) / (l1+l2)
        print("Overall error rate: {:.2%}".format(overall_err), file=outstream)
