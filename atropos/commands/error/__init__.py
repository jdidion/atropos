"""Estimate the empircal error rate.
"""
from collections import defaultdict
import csv
import re
from atropos import AtroposError
from atropos.commands.base import (
    BaseCommandRunner, Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin)
from atropos.io import open_output
from atropos.util import qual2prob, run_interruptible

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
        
        estimator_args = dict(max_read_len=self.max_bases)
        if self.paired:
            estimator = PairedErrorEstimator(
                estimator_class=estimator_class, **estimator_args)
        else:
            estimator = estimator_class(**estimator_args)
        
        self.summary['errorrate'] = estimator_args
        
        # currently only single-threaded operation is supproted
        self.summary.update(mode='serial', threads=1)
        return run_interruptible(estimator, self, raise_on_error=True)

class ErrorEstimator(SingleEndPipelineMixin, Pipeline):
    """Base class for error estimators.
    """
    def __init__(self, max_read_len):
        super().__init__()
        self.total_len = 0
        self.max_read_len = max_read_len
    
    def handle_reads(self, context, read1, read2=None):
        raise NotImplementedError()
    
    def estimate(self):
        """Returns an estimate of the error rate.
        """
        raise NotImplementedError()
    
    def finish(self, summary, **kwargs):
        super().finish(summary)
        estimate, details = self.estimate()
        summary['errorrate'].update(
            estimate=(estimate,), 
            total_len=(self.total_len,),
            details=(details,))

class BaseQualityErrorEstimator(ErrorEstimator):
    """Simple error estimation using base qualities. It is well-known that base
    qualities significantly overestimate true error rates, so take the estimate
    with a grain of salt.
    """
    def __init__(self, max_read_len=None):
        super().__init__(max_read_len)
        self.total_qual = 0.0
    
    def handle_reads(self, context, read1, read2=None):
        quals = read1.qualities
        readlen = len(quals)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            quals = quals[:readlen]
        self.total_qual += sum(qual2prob(qchar) for qchar in quals)
        self.total_len += readlen
    
    def estimate(self):
        return (self.total_qual / self.total_len, None)

# Error estimation using shadow counts
# DEPRECATED
# TODO: Replace with Zhu et al. cubic splines method or something faster
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1052-3

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
        super().__init__(max_read_len)
        self.seqs = defaultdict(lambda: 0)
        self.method = method
        self.rscript_exe = rscript_exe
    
    def handle_reads(self, context, read1, read2=None):
        seq = read1.sequence
        readlen = len(seq)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            seq = seq[:readlen]
        if FILTER_RE.fullmatch(seq):
            return
        self.seqs[seq] += 1
        self.total_len += readlen
    
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
    
    def _write_read_counts(self, fileobj):
        writer = csv.writer(fileobj, delimiter=" ")
        writer.writerows(sorted(
            self.seqs.items(), reverse=True, key=lambda i: i[1]))

class PairedErrorEstimator(PairedEndPipelineMixin, Pipeline):
    """Estimator for a pair of input files.
    """
    def __init__(
            self, estimator_class=BaseQualityErrorEstimator, **kwargs):
        super().__init__()
        self.estimator1 = estimator_class(**kwargs)
        self.estimator2 = estimator_class(**kwargs)
    
    def handle_reads(self, context, read1, read2):
        self.estimator1.handle_reads(context, read1)
        self.estimator2.handle_reads(context, read2)
    
    def finish(self, summary, **kwargs):
        """Estimate error rates.
        
        Returns:
            Tuple (read1_estimate, read2_estimate).
        """
        super().finish(summary)
        estimate1, details1 = self.estimator1.estimate()
        estimate2, details2 = self.estimator2.estimate()
        summary['errorrate'].update(
            estimate=(estimate1, estimate2),
            total_len=(self.estimator1.total_len, self.estimator2.total_len),
            details=(details1, details2))
