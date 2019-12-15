from abc import ABCMeta, abstractmethod
from collections import defaultdict
import csv
import re
from pathlib import Path
from typing import Optional, Tuple, Union

from xphyle import open_

from atropos.errors import AtroposError
from atropos.commands import (
    BaseCommand,
    Pipeline,
    SingleEndPipelineMixin,
    PairedEndPipelineMixin,
)
from atropos.io.sequence import Sequence
from atropos.utils import ReturnCode, classproperty, run_interruptible
from atropos.utils.ngs import qual2prob


class ErrorCommand(BaseCommand):
    @classproperty
    def name(cls) -> str:
        return "error"

    def __call__(self) -> ReturnCode:
        if not self.get_option("delivers_qualities"):
            raise ValueError("Cannot estimate error rate without base qualities")

        algorithm = self.get_option("algorithm")
        if algorithm == "quality":
            estimator_class = BaseQualityErrorEstimator
        elif algorithm == "shadow":
            estimator_class = ShadowRegressionErrorEstimator
        else:
            raise ValueError(f"Invalid algorithm: {algorithm}")

        estimator_args = dict(max_read_len=self.get_option("max_bases"))

        if self.get_option("paired"):
            estimator = PairedErrorEstimator(
                estimator_class=estimator_class, **estimator_args
            )
        else:
            estimator = estimator_class(**estimator_args)

        self.summary["errorrate"] = estimator_args
        # currently only single-threaded operation is supproted
        self.summary.update(mode="serial", threads=1)

        return run_interruptible(estimator, self, raise_on_error=True)


class ErrorEstimator(SingleEndPipelineMixin, Pipeline, metaclass=ABCMeta):
    """
    Base class for error estimators.
    """

    def __init__(self, max_read_len: int):
        super().__init__()
        self.total_len = 0
        self.max_read_len = max_read_len

    @abstractmethod
    def estimate(self) -> Tuple[float, Optional[dict]]:
        """
        Estimate the error rate.
        """

    def finish(self, summary: dict, **kwargs):
        super().finish(summary)
        estimate, details = self.estimate()
        summary["errorrate"].update(
            estimate=(estimate,), total_len=(self.total_len,), details=(details,)
        )


class BaseQualityErrorEstimator(ErrorEstimator):
    """
    Simple error estimation using base qualities. It is well-known that base
    qualities significantly overestimate true error rates, so take the estimate with
    a grain of salt.
    """

    def __init__(self, max_read_len=None):
        super().__init__(max_read_len)
        self.total_qual = 0.0

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        quals = read1.qualities
        readlen = len(quals)
        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            quals = quals[:readlen]
        self.total_qual += sum(qual2prob(qchar) for qchar in quals)
        self.total_len += readlen

    def estimate(self) -> Tuple[float, Optional[dict]]:
        return self.total_qual / self.total_len, None


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
    """
    Re-implementation of the shadow regression method described in:
    Wang et al., "Estimation of sequencing error rates in short reads",
        BMC Bioinformatics 2012 13:185, DOI: 10.1186/1471-2105-13-185
    """

    def __init__(
        self,
        method: str = "sub",
        max_read_len: Optional[int] = None,
        rscript_exe: Union[str, Path] = "Rscript"
    ):
        """
        Args:
            method: The differences that are considered in the error rate
                calculation; sub = substitutions, indel = insertions and
                deletions; all = both substitutions and indels.
            max_read_len: The maximum number of bases (starting from the 5' end)
                to consider from each read.
            rscript_exe:
        """
        super().__init__(max_read_len)
        self.seqs = defaultdict(lambda: 0)
        self.method = method
        self.rscript_exe = rscript_exe

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        seq = read1.sequence
        readlen = len(seq)

        if self.max_read_len and self.max_read_len < readlen:
            readlen = self.max_read_len
            seq = seq[:readlen]

        if FILTER_RE.fullmatch(seq):
            return

        self.seqs[seq] += 1
        self.total_len += readlen

    def estimate(self) -> Tuple[float, Optional[dict]]:
        # This is a temporary solution that requires R and the
        # ShadowRegression package. Eventually, this will be
        # replaced with a pure-python implementation.
        import subprocess
        import tempfile

        tempfiles = (Path(tempfile.mkstemp()[1]) for _ in range(4))
        read_counts, per_read, per_cycle, script_file = tempfiles

        try:
            # write counts to a file
            with open_(read_counts, "wt") as out:
                self._write_read_counts(out)

            # execute R script
            script = SHADOW_REGRESSION_R_SCRIPT.format(
                reads=read_counts,
                method=self.method,
                per_read=per_read,
                per_cycle=per_cycle,
            )

            with open_(script_file, "wt") as out:
                out.write(script)

            proc = subprocess.Popen(
                [str(self.rscript_exe), "--vanilla", script_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            with proc:
                stdout, stderr = proc.communicate()
                if proc.returncode != 0:
                    raise AtroposError(
                        f"R script failed: rc={proc.returncode}; stdout={stdout}; "
                        f"stderr={stderr}"
                    )

            # read the results
            with open_(per_read, "rt") as i:
                reader = csv.reader(i, delimiter="\t")
                per_read_error = dict(reader)
                if len(per_read_error) != 4:
                    raise AtroposError("Invalid output from R script")

            with open_(per_cycle, "rt") as i:
                reader = csv.reader(i, delimiter="\t")
                per_cycle_error = list(row[0:3] for row in reader)
                if not per_cycle_error:
                    raise AtroposError("Invalid output from R script")

            return (
                per_read_error["error rate"],
                dict(per_read=per_read_error, per_cycle=per_cycle_error),
            )

        finally:
            for path in tempfiles:
                path.unlink()

    def _write_read_counts(self, fileobj):
        writer = csv.writer(fileobj, delimiter=" ")
        writer.writerows(sorted(self.seqs.items(), reverse=True, key=lambda i: i[1]))


class PairedErrorEstimator(PairedEndPipelineMixin, Pipeline):
    """
    Estimator for a pair of input files.
    """

    def __init__(self, estimator_class=BaseQualityErrorEstimator, **kwargs):
        super().__init__()
        self.estimator1 = estimator_class(**kwargs)
        self.estimator2 = estimator_class(**kwargs)

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        self.estimator1.handle_reads(context, read1)
        self.estimator2.handle_reads(context, read2)

    def finish(self, summary, **kwargs):
        """
        Estimates error rates.

        Returns:
            Tuple (read1_estimate, read2_estimate).
        """
        super().finish(summary)
        estimate1, details1 = self.estimator1.estimate()
        estimate2, details2 = self.estimator2.estimate()
        summary["errorrate"].update(
            estimate=(estimate1, estimate2),
            total_len=(self.estimator1.total_len, self.estimator2.total_len),
            details=(details1, details2),
        )
