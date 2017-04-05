"""Estimate the empircal error rate.
"""
from collections import defaultdict
import csv
import re
from atropos import AtroposError
from atropos.commands.base import BaseCommandRunner
from atropos.io import open_output
from atropos.util import qual2prob

class CommandRunner(BaseCommandRunner):
    name = 'error'
    
    def __call__(self):
        if not self.delivers_qualities:
            raise ValueError(
                "Cannot estimate error rate without base qualities")
        
        if self.algorithm == 'quality':
            estimator_class = BaseQualityErrorEstimator
        elif self.algorithm == 'shadow':
            estimator_class = ShadowRegressionErrorEstimator
        
        if self.paired:
            estimator = PairedErrorEstimator(
                max_read_len=self.max_bases,
                estimator_class=estimator_class)
        else:
            estimator = estimator_class(max_read_len=self.max_bases)
        
        self.summary['error_estimator'] = estimator
        
        estimator.consume_all_batches(self)
        
        return 0

class ErrorEstimator(object):
    """Base class for error estimators.
    """
    def consume_all_batches(self, batch_iterator):
        """Consume all batches in `batch_iterator`.
        """
        for batch_meta, batch in batch_iterator:
            if batch_meta['size'] == 0:
                continue
            for read in batch:
                self.consume(read)
    
    def consume(self, read):
        """Consume a sequence record.
        """
        raise NotImplementedError()
    
    def estimate(self):
        """Returns an estimate of the error rate.
        """
        raise NotImplementedError()
    
    # TODO: move report to the report package, have summarize update the
    # summary dict.
    
    def summarize(self, outstream, name=None, show_details=True, indent="  "):
        """Print a summary report.
        """
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
        """Print additional details.
        """
        pass

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

SHADOW_REGRESSION_R_SCRIPT = """
library(ShadowRegression)
errorRates = getErrorRates("{reads}", type="{method}")
write.table(errorRates$perReadER, "{per_read}", sep="\t", quote=F, col.names=F,
            row.names=T)
write.table(errorRates$cycleER, "{per_cycle}", sep="\t", quote=F, col.names=F,
            row.names=T)
"""

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
    
    def _write_read_counts(self, fileobj):
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
            with open(read_counts, 'wt') as out:
                self._write_read_counts(out)
            
            # execute R script
            script = SHADOW_REGRESSION_R_SCRIPT.format(
                reads=read_counts, method=self.method, per_read=per_read,
                per_cycle=per_cycle)
            with open(script_file, 'wt') as out:
                out.write(script)
            proc = subprocess.Popen(
                [self.rscript_exe, "--vanilla", script_file],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with proc:
                stdout, stderr = proc.communicate()
                if proc.returncode != 0:
                    raise AtroposError(
                        "R script failed: rc={}; stdout={}; stderr={}".format(
                            proc.returncode, stdout, stderr))
            
            # read the results
            with open(per_read, 'rt') as i:
                reader = csv.reader(i, delimiter="\t")
                per_read_error = dict(reader)
                if len(per_read_error) != 4:
                    raise AtroposError("Invalid output from R script")
            with open(per_cycle, 'rt') as i:
                reader = csv.reader(i, delimiter="\t")
                per_cycle_error = list(row[0:3] for row in reader)
                if not per_cycle_error:
                    raise AtroposError("Invalid output from R script")
            
            return (
                per_read_error["error rate"],
                dict(per_read=per_read_error, per_cycle=per_cycle_error))
        finally:
            for path in tempfiles:
                os.remove(path)
    
    def print_details(self, outstream, details, indent):
        per_read = details['per_read']
        per_cycle = details['per_cycle']
        
        print(
            "{}StdErr: {:.2%}".format(indent, per_read['standard error']),
            file=outstream)
        print("{}Per-cycle rates:".format(indent), file=outstream)
        for cycle in per_cycle:
            print("{}Cycle: {}, Error: {:.2%}, StdErr: {:.2%}".format(
                indent*2, *cycle))

class PairedErrorEstimator(object):
    """Estimator for a pair of input files.
    """
    def __init__(self, max_read_len=None,
                 estimator_class=BaseQualityErrorEstimator):
        self.estimator1 = estimator_class(max_read_len=max_read_len)
        self.estimator2 = estimator_class(max_read_len=max_read_len)
    
    def consume_all_batches(self, batch_iterator):
        """Consume all batches in iterator.
        """
        for batch_meta, batch in batch_iterator:
            if batch_meta['size'] == 0:
                continue
            for read1, read2 in batch:
                self.estimator1.consume(read1)
                self.estimator2.consume(read2)
    
    def estimate(self):
        """Estimate error rates.
        
        Returns:
            Tuple (read1_estimate, read2_estimate).
        """
        return (self.estimator1.estimate(), self.estimator2.estimate())
    
    def summarize(self, outstream, names=None, show_details=True, indent="  "):
        """Print a summary of the error detection.
        """
        print("", file=outstream)
        if names:
            header = "File 1: {}".format(names[0])
            print(header, file=outstream)
            print('-' * len(header), file=outstream)
        err1, details1 = self.estimator1.estimate()
        print("Error rate: {:.2%}".format(err1), file=outstream)
        if show_details and details1:
            print("Details:\n", file=outstream)
            self.estimator1.print_details(details1, outstream, indent)
        print("", file=outstream)
        
        if names:
            header = "File 2: {}".format(names[1])
            print(header, file=outstream)
            print('-' * len(header), file=outstream)
        err2, details2 = self.estimator2.estimate()
        print("Error rate: {:.2%}".format(err2), file=outstream)
        if show_details and details2:
            print("Details:\n", file=outstream)
            self.estimator2.print_details(details2, outstream, indent)
        print("", file=outstream)
        
        len1 = self.estimator1.total_len
        len2 = self.estimator2.total_len
        overall_err = ((err1 * len1) + (err2 * len2)) / (len1+len2)
        print("Overall error rate: {:.2%}".format(overall_err), file=outstream)
