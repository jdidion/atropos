from abc import ABCMeta

from loguru import logger

from atropos.commands import (
    BaseCommand,
    Pipeline,
    SingleEndPipelineMixin,
    PairedEndPipelineMixin,
)
from atropos.commands.metrics import SingleEndReadMetrics, PairedEndReadMetrics
from atropos.utils import ReturnCode, classproperty, run_interruptible


class QcPipeline(Pipeline, metaclass=ABCMeta):
    """
    Base Pipeline for the qc command.
    """

    def __init__(self, read_statistics_class, **kwargs):
        super().__init__()
        self.read_statistics_class = read_statistics_class
        self._metrics = {}
        self._metrics_kwargs = kwargs

    def _get_metrics(self, source):
        if source not in self._metrics:
            self._metrics[source] = self.read_statistics_class(**self._metrics_kwargs)
        return self._metrics[source]

    def handle_reads(self, context, read1, read2=None):
        self._get_metrics(context["source"]).collect(read1, read2)

    def finish(self, summary, **kwargs):
        super().finish(summary)
        summary["pre"] = dict(
            (source, metrics.summarize()) for source, metrics in self._metrics.items()
        )


class SingleEndQcPipeline(SingleEndPipelineMixin, QcPipeline):
    """
    QcPipeline for single-end data.
    """

    def __init__(self, **kwargs):
        super().__init__(SingleEndReadMetrics, **kwargs)


class PairedEndQcPipeline(PairedEndPipelineMixin, QcPipeline):
    """
    QcPipeline for paired-end data.
    """

    def __init__(self, **kwargs):
        super().__init__(PairedEndReadMetrics, **kwargs)


class QcCommand(BaseCommand):
    @classproperty
    def name(cls) -> str:
        return "qc"

    def __call__(self) -> ReturnCode:
        if self.get_option("paired"):
            pipeline_class = PairedEndQcPipeline
        else:
            pipeline_class = SingleEndQcPipeline

        pipeline_args = dict(
            qualities=self.get_option("delivers_qualities"),
            quality_base=self.get_option("quality_base"),
        )
        if self.get_option("metrics"):
            pipeline_args.update(self.get_option("metrics"))

        threads = self.get_option("threads", 1)

        if threads <= 1:
            self.summary.update(mode="serial", threads=1)
            pipeline = pipeline_class(**pipeline_args)
            return run_interruptible(pipeline, self)
        else:
            self.summary.update(mode="parallel", threads=threads)
            return self.run_parallel(pipeline_class, pipeline_args)

    def run_parallel(self, pipeline_class, pipeline_args):
        """
        Executes qc in parallel mode.

        Args:
            pipeline_class: Pipeline class to instantiate.
            pipeline_args: Arguments to pass to Pipeline constructor.

        Returns:
            The return code.
        """
        from atropos.commands.multicore import (
            ParallelPipelineMixin,
            ParallelPipelineRunner,
        )

        logger.debug(
            f"Starting atropos qc in parallel mode with "
            f"threads={self.get_option('threads')}, "
            f"timeout={self.get_option('timeout')}"
        )

        if self.get_option("threads") < 2:
            raise ValueError("'threads' must be >= 2")

        # Start worker processes, reserve a thread for the reader process,
        # which we will get back after it completes
        pipeline_class = type("QcPipelineImpl", (ParallelPipelineMixin, pipeline_class))
        pipeline = pipeline_class(**pipeline_args)

        runner = ParallelPipelineRunner(self, pipeline)

        return runner.run()
