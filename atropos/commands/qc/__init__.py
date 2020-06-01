# coding: utf-8
"""Implementation of the 'qc' command.
"""
import logging
from atropos.commands.base import (
    BaseCommandRunner, Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin)
from atropos.commands.stats import (
    SingleEndReadStatistics, PairedEndReadStatistics)
from atropos.util import run_interruptible

class QcPipeline(Pipeline):
    """Base Pipeline for the qc command.
    """
    def __init__(self, read_statistics_class, **kwargs):
        super().__init__()
        self.read_statistics_class = read_statistics_class
        self.stats = {}
        self.stats_kwargs = kwargs

    def _get_stats(self, source):
        if source not in self.stats:
            self.stats[source] = self.read_statistics_class(**self.stats_kwargs)
        return self.stats[source]

    def handle_reads(self, context, read1, read2=None):
        self._get_stats(context['source']).collect(read1, read2)

    def finish(self, summary, **kwargs):
        super().finish(summary)
        summary['pre'] = dict(
            (source, stats.summarize())
            for source, stats in self.stats.items())

class SingleEndQcPipeline(SingleEndPipelineMixin, QcPipeline):
    """QcPipeline for single-end data.
    """
    def __init__(self, **kwargs):
        super().__init__(SingleEndReadStatistics, **kwargs)

class PairedEndQcPipeline(PairedEndPipelineMixin, QcPipeline):
    """QcPipeline for paired-end data.
    """
    def __init__(self, **kwargs):
        super().__init__(PairedEndReadStatistics, **kwargs)

class CommandRunner(BaseCommandRunner):
    name = 'qc'

    def __call__(self):
        if self.paired:
            pipeline_class = PairedEndQcPipeline
        else:
            pipeline_class = SingleEndQcPipeline
        pipeline_args = dict(
            qualities=self.delivers_qualities,
            quality_base=self.quality_base)
        if self.stats:
            pipeline_args.update(self.stats)

        if self.threads is None:
            self.summary.update(mode='serial', threads=1)
            pipeline = pipeline_class(**pipeline_args)
            return run_interruptible(pipeline, self)
        else:
            self.summary.update(mode='parallel', threads=self.threads)
            return self.run_parallel(pipeline_class, pipeline_args)

    def run_parallel(self, pipeline_class, pipeline_args):
        """Execute qc in parallel mode.

        Args:
            pipeline_class: Pipeline class to instantiate.
            pipeline_args: Arguments to pass to Pipeline constructor.

        Returns:
            The return code.
        """
        from atropos.commands.multicore import (
            ParallelPipelineMixin, ParallelPipelineRunner, RETRY_INTERVAL)

        pipeline_class = type(
            'QcPipelineImpl', (ParallelPipelineMixin, pipeline_class))
        pipeline = pipeline_class(**pipeline_args)
        runner = ParallelPipelineRunner(self, pipeline)

        logging.getLogger().debug(
            "Starting atropos qc in parallel mode with threads=%d, timeout=%d",
            runner.threads, runner.timeout)

        return runner.run()
