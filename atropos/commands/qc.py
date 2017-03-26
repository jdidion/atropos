# coding: utf-8
"""Implementation of the 'qc' command.
"""
import logging
from atropos.commands import (
    Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin, create_reader)
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
    
    def finish(self, summary, worker=None):
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

def execute(options, summary):
    """Execute the qc command.
    
    Args:
        options: Command-line options.
        summary: The summary dict.
    """
    reader, _, qualities, _ = create_reader(options)
    if options.paired:
        pipeline_class = PairedEndQcPipeline
    else:
        pipeline_class = SingleEndQcPipeline
    pipeline_args = dict(
        qualities=qualities,
        tile_key_regexp=options.tile_key_regexp)
    
    if options.threads is None:
        summary.update(mode='serial', threads=1)
        pipeline = pipeline_class(**pipeline_args)
        retcode = run_interruptible(pipeline, reader, summary)
    else:
        summary.update(mode='parallel', threads=options.threads)
        retcode = run_parallel(
            reader, pipeline_class, pipeline_args, summary, options.threads,
            options.process_timeout, options.read_queue_size)
    
    reader.close()
    
    return retcode

def run_parallel(
        reader, pipeline_class, pipeline_args, summary, threads=2,
        timeout=30, input_queue_size=0):
    """Execute qc in parallel mode.
    
    Args:
        reader: Iterator over batches of reads (most likely a BatchIterator)
        read_stats: Template ReadStatistics object.
        threads: Number of worker threads to use; additional threads are used
            for the main proccess and the writer process (if requested).
        timeout: number of seconds after which waiting processes escalate their
            messages from DEBUG to ERROR.
        input_queue_size: max number of items that can be in the input queue,
            or 0 for no limit (be warned that this could explode memory usage)
    """
    from atropos.commands.multicore import (
        ParallelPipelineMixin, ParallelPipelineRunner)
    
    logging.getLogger().debug(
        "Starting atropos qc in parallel mode with threads=%d, timeout=%d",
        threads, timeout)
    
    if threads < 2:
        raise ValueError("'threads' must be >= 2")
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    pipeline_class = type(
        'QcPipelineImpl', (ParallelPipelineMixin, pipeline_class))
    pipeline = pipeline_class(**pipeline_args)
    runner = ParallelPipelineRunner(
        reader, pipeline, threads, input_queue_size, timeout)
    return runner.run(summary)
