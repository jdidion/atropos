"""Implementation of the 'trim' command.
"""
from collections import Sequence, defaultdict
import logging
import os
import sys
import textwrap
from atropos.commands.base import (
    BaseCommandRunner, Summary, Pipeline, SingleEndPipelineMixin,
    PairedEndPipelineMixin)
from atropos.commands.stats import (
    SingleEndReadStatistics, PairedEndReadStatistics)
from atropos.adapters import AdapterParser, BACK
from atropos.io import STDOUT
from atropos.util import RandomMatchProbability, Const, run_interruptible
from .modifiers import (
    AdapterCutter, DoubleEncoder, InsertAdapterCutter, LengthTagModifier,
    MergeOverlapping, MinCutter, NEndTrimmer, NextseqQualityTrimmer,
    NonDirectionalBisulfiteTrimmer, OverwriteRead, PairedEndModifiers,
    PrefixSuffixAdder, PrimerTrimmer, QualityTrimmer, RRBSTrimmer,
    SingleEndModifiers, SuffixRemover, SwiftBisulfiteTrimmer,
    UnconditionalCutter, ZeroCapper)
from .filters import (
    FilterFactory, Filters, MergedReadFilter, NContentFilter, NoFilter,
    TooLongReadFilter, TooShortReadFilter, TrimmedFilter, UntrimmedFilter)
from .writers import (
    Formatters, InfoFormatter, RestFormatter, WildcardFormatter, Writers)

class TrimPipeline(Pipeline):
    """Base trimming pipeline.
    
    Args:
        record_handler:
        result_handler:
    """
    def __init__(self, record_handler, result_handler):
        super().__init__()
        self.record_handler = record_handler
        self.result_handler = result_handler
    
    def start(self, worker=None):
        self.result_handler.start(worker)
    
    def add_to_context(self, context):
        context['results'] = defaultdict(lambda: [])
    
    def handle_records(self, context, records):
        super().handle_records(context, records)
        self.result_handler.write_result(context['index'], context['results'])
    
    def handle_reads(self, context, read1, read2=None):
        return self.record_handler.handle_record(context, read1, read2)
    
    def finish(self, summary, **kwargs):
        self.result_handler.finish()
        super().finish(summary)
        summary.update(self.record_handler.summarize())

class RecordHandler(object):
    """Base class for record handlers.
    """
    def __init__(self, modifiers, filters, formatters):
        self.modifiers = modifiers
        self.filters = filters
        self.formatters = formatters
    
    def handle_record(self, context, read1, read2=None):
        """Handle a pair of reads.
        """
        reads = self.modifiers.modify(read1, read2)
        dest = self.filters.filter(*reads)
        self.formatters.format(context['results'], dest, *reads)
        return (dest, reads)
    
    def summarize(self):
        """Returns a summary dict.
        """
        return dict(trim=dict(
            modifiers=self.modifiers.summarize(),
            filters=self.filters.summarize(),
            formatters=self.formatters.summarize()))

class StatsRecordHandlerWrapper(object):
    """Wrapper around a record handler that collects read statistics
    before and/or after trimming.
    
    Args:
        record_handler:
        paired: Whether reads are paired-end.
        mode: Collection mode; pre=only collect pre-trim stats; post=only
            collect post-trim stats; both=collect both pre- and post-trim stats.
    """
    def __init__(self, record_handler, paired, stats_args, **kwargs):
        self.record_handler = record_handler
        self.read_statistics_class = (
            PairedEndReadStatistics if paired else SingleEndReadStatistics)
        self.pre = self.post = None
        if 'pre' in stats_args:
            self.pre = {}
            self.pre_kwargs = kwargs.copy()
            self.pre_kwargs.update(stats_args['pre'])
        if 'post' in stats_args:
            self.post = {}
            self.post_kwargs = kwargs.copy()
            self.post_kwargs.update(stats_args['post'])
    
    def handle_record(self, context, read1, read2=None):
        """Handle a pair of reads.
        """
        if self.pre is not None:
            self.collect(
                self.pre, context['source'], read1, read2, **self.pre_kwargs)
        dest, reads = self.record_handler.handle_record(context, read1, read2)
        if self.post is not None:
            if dest not in self.post:
                self.post[dest] = {}
            self.collect(
                self.post[dest], context['source'], *reads, **self.post_kwargs)
        return (dest, reads)
    
    def collect(self, stats, source, read1, read2=None, **kwargs):
        """Collect stats on a pair of reads.
        
        Args:
            stats: The :class:`ReadStatistics` object.
            source: The source file(s).
            read1, read2: The reads.
        """
        if source not in stats:
            stats[source] = self.read_statistics_class(**kwargs)
        stats[source].collect(read1, read2)
    
    def summarize(self):
        """Returns a summary dict.
        """
        summary = self.record_handler.summarize()
        if self.pre is not None:
            summary['pre'] = dict(
                (source, stats.summarize())
                for source, stats in self.pre.items())
        if self.post is not None:
            summary['post'] = {}
            for dest, stats_dict in self.post.items():
                summary['post'][dest.name] = dict(
                    (source, stats.summarize())
                    for source, stats in stats_dict.items())
        return summary

class ResultHandler(object):
    """Base class for result handlers.
    """
    def start(self, worker=None):
        """Start the result handler.
        """
        pass
    
    def finish(self, total_batches=None):
        """Finish the result handler.
        
        Args:
            total_batches: Total number of batches processed.
        """
        pass
    
    def write_result(self, batch_num, result):
        """Write a batch of results to output.
        
        Args:
            batch_num: The batch number.
            result: The result to write.
        """
        raise NotImplementedError()

class ResultHandlerWrapper(ResultHandler):
    """Wraps a ResultHandler.
    """
    def __init__(self, handler):
        self.handler = handler
    
    def start(self, worker):
        self.handler.start(worker)
    
    def write_result(self, batch_num, result):
        self.handler.write_result(batch_num, result)
    
    def finish(self, total_batches=None):
        self.handler.finish(total_batches=total_batches)

class WorkerResultHandler(ResultHandlerWrapper):
    """Wraps a ResultHandler and compresses results prior to writing.
    """
    def write_result(self, batch_num, result):
        """Given a dict mapping files to lists of strings, join the strings and
        compress them (if necessary) and then return the property formatted
        result dict.
        """
        self.handler.write_result(
            batch_num, dict(
                self.prepare_file(*item)
                for item in result.items()))
    
    def prepare_file(self, path, strings):
        """Prepare data for writing.
        
        Returns:
            Tuple (path, data).
        """
        return (path, "".join(strings))

class WriterResultHandler(ResultHandler):
    """ResultHandler that writes results to disk.
    
    Args:
        writers: :class:`Writers` object.
        compressed: Whether the data is compressed.
        use_suffix: Whether to add the worker index as a file suffix. Used for
            parallel-write mode.
    """
    def __init__(self, writers, compressed=False, use_suffix=False):
        self.writers = writers
        self.compressed = compressed
        self.use_suffix = use_suffix
    
    def start(self, worker=None):
        if self.use_suffix:
            if worker is None:
                raise ValueError("")
            self.writers.suffix = ".{}".format(worker.index)
    
    def write_result(self, batch_num, result):
        self.writers.write_result(result, self.compressed)
    
    def finish(self, total_batches=None):
        self.writers.close()

class TrimSummary(Summary):
    """Summary that adds aggregate values for record and bp stats.
    """
    def _post_process_other(self, dict_val, key, value):
        """For trim stats, any value with a name that starts with 'records_'
        will have 'fraction_<var>' computed as value / total_records, and any
        value with a name that starts with 'bp_' will have 'fraction_<var>' and
        'total_<var>' computed. We also replace any `Const`s with their values.
        """
        if self.has_exception:
            return
        
        def frac(val, total):
            """Compute fraction of total.
            """
            return (val / total) if val and total != 0 else 0
        
        if isinstance(key, str):
            if key.startswith('records_'):
                frac_key = "fraction_{}".format(key)
                total_records = self['total_record_count']
                if isinstance(value, Sequence):
                    dict_val[frac_key] = [
                        frac(val, total_records) for val in value]
                    total = sum(val for val in value if val)
                    dict_val["total_{}".format(key)] = total
                else:
                    dict_val[frac_key] = frac(value, total_records)
            elif key.startswith('bp_'):
                frac_key = "fraction_{}".format(key)
                sum_total_bp = self['sum_total_bp_count']
                if isinstance(value, Sequence):
                    dict_val[frac_key] = [
                        frac(val, bps)
                        for val, bps in zip(value, self['total_bp_counts'])]
                    total = sum(val for val in value if val)
                    dict_val["total_{}".format(key)] = total
                    dict_val["fraction_total_{}".format(key)] = \
                        frac(total, sum_total_bp)
                else:
                    dict_val[frac_key] = frac(value, sum_total_bp)

class CommandRunner(BaseCommandRunner):
    name = 'trim'
    
    def __init__(self, options):
        super().__init__(options, TrimSummary)
    
    def __call__(self):
        options = self.options
        match_probability = RandomMatchProbability()
        
        # Create Adapters
        
        has_adapters1 = options.adapters or options.anywhere or options.front
        has_adapters2 = options.adapters2 or options.anywhere2 or options.front2
        
        adapters1 = adapters2 = []
        if has_adapters1 or has_adapters2:
            adapter_cache = super().load_known_adapters()
            parser_args = dict(
                colorspace=options.colorspace,
                max_error_rate=options.error_rate,
                min_overlap=options.overlap,
                read_wildcards=options.match_read_wildcards,
                adapter_wildcards=options.match_adapter_wildcards,
                indels=options.indels, indel_cost=options.indel_cost,
                cache=adapter_cache, gc_content=options.gc_content,
                match_probability=match_probability, alphabet=options.alphabet)
            if options.adapter_max_rmp:
                parser_args['max_rmp'] = options.adapter_max_rmp
            adapter_parser = AdapterParser(**parser_args)
            
            if has_adapters1:
                adapters1 = adapter_parser.parse_multi(
                    options.adapters, options.anywhere, options.front)
            if has_adapters2:
                adapters2 = adapter_parser.parse_multi(
                    options.adapters2, options.anywhere2, options.front2)
            
            if options.cache_adapters:
                adapter_cache.save()
        
        # Create Modifiers
        
        # TODO: can this be replaced with an argparse required group?
        if (
                not adapters1 and not adapters2 and
                not options.quality_cutoff and
                options.nextseq_trim is None and
                options.cut == [] and options.cut2 == [] and
                options.cut_min == [] and options.cut_min2 == [] and
                (
                    options.minimum_length is None or
                    options.minimum_length <= 0) and
                options.maximum_length == sys.maxsize and not options.trim_n and
                not self.has_qualfile and options.max_n is None and
                (not options.paired or options.overwrite_low_quality is None)):
            raise ValueError(
                "You need to provide at least one adapter sequence.")
        
        if (
                options.aligner == 'insert' and any(
                    not a or len(a) != 1 or a[0].where != BACK
                    for a in (adapters1, adapters2))):
            raise ValueError(
                "Insert aligner requires a single 3' adapter for each read")
        
        if options.debug:
            for adapter in adapters1 + adapters2:
                adapter.enable_debug()
        
        if options.paired:
            modifiers = PairedEndModifiers(options.paired)
        else:
            modifiers = SingleEndModifiers()
        
        for oper in options.op_order:
            if oper == 'W' and options.overwrite_low_quality:
                lowq, highq, window = options.overwrite_low_quality
                modifiers.add_modifier(
                    OverwriteRead,
                    worse_read_min_quality=lowq, better_read_min_quality=highq,
                    window_size=window, base=options.quality_base)
                
            elif oper == 'A' and (adapters1 or adapters2):
                # TODO: generalize this using some kind of factory class
                if options.aligner == 'insert':
                    # Use different base probabilities if we're trimming
                    # bisulfite data.
                    # TODO: this doesn't seem to help things, so commenting it
                    # out for now
                    #if options.bisulfite:
                    #   base_probs = dict(match_prob=0.33, mismatch_prob=0.67)
                    # else:
                    #   base_probs = dict(match_prob=0.25, mismatch_prob=0.75)
                    modifiers.add_modifier(
                        InsertAdapterCutter,
                        adapter1=adapters1[0], adapter2=adapters2[0],
                        action=options.action,
                        mismatch_action=options.correct_mismatches,
                        max_insert_mismatch_frac=\
                            options.insert_match_error_rate,
                        max_adapter_mismatch_frac=\
                            options.insert_match_adapter_error_rate,
                        match_probability=match_probability,
                        insert_max_rmp=options.insert_max_rmp,
                        read_wildcards=options.match_read_wildcards,
                        adapter_wildcards=options.match_adapter_wildcards)
                else:
                    a1_args = dict(
                        adapters=adapters1,
                        times=options.times,
                        action=options.action) if adapters1 else None
                    a2_args = dict(
                        adapters=adapters2,
                        times=options.times,
                        action=options.action) if adapters2 else None
                    modifiers.add_modifier_pair(AdapterCutter, a1_args, a2_args)
            elif oper == 'C' and (options.cut or options.cut2):
                modifiers.add_modifier_pair(
                    UnconditionalCutter,
                    dict(lengths=options.cut),
                    dict(lengths=options.cut2))
            elif oper == 'G' and (options.nextseq_trim is not None):
                modifiers.add_modifier(
                    NextseqQualityTrimmer,
                    cutoff=options.nextseq_trim,
                    base=options.quality_base)
            elif oper == 'Q' and options.quality_cutoff:
                modifiers.add_modifier(
                    QualityTrimmer,
                    cutoff_front=options.quality_cutoff[0],
                    cutoff_back=options.quality_cutoff[1],
                    base=options.quality_base)
        
        if options.bisulfite:
            if isinstance(options.bisulfite, str):
                if "non-directional" in options.bisulfite:
                    modifiers.add_modifier(
                        NonDirectionalBisulfiteTrimmer,
                        rrbs=options.bisulfite=="non-directional-rrbs")
                elif options.bisulfite == "rrbs":
                    modifiers.add_modifier(RRBSTrimmer)
                elif options.bisulfite in ("epignome", "truseq"):
                    # Trimming leads to worse results
                    #modifiers.add_modifier(TruSeqBisulfiteTrimmer)
                    pass
                elif options.bisulfite == "swift":
                    modifiers.add_modifier(SwiftBisulfiteTrimmer)
            else:
                if options.bisulfite[0]:
                    modifiers.add_modifier(
                        MinCutter, read=1, **(options.bisulfite[0]))
                if len(options.bisulfite) > 1 and options.bisulfite[1]:
                    modifiers.add_modifier(
                        MinCutter, read=2, **(options.bisulfite[1]))
        
        if options.trim_n:
            modifiers.add_modifier(NEndTrimmer)
        
        if options.cut_min or options.cut_min2:
            modifiers.add_modifier_pair(
                MinCutter,
                dict(lengths=options.cut_min),
                dict(lengths=options.cut_min2))
        
        if options.length_tag:
            modifiers.add_modifier(
                LengthTagModifier, length_tag=options.length_tag)
        
        if options.strip_suffix:
            modifiers.add_modifier(SuffixRemover, suffixes=options.strip_suffix)
        
        if options.prefix or options.suffix:
            modifiers.add_modifier(
                PrefixSuffixAdder, prefix=options.prefix, suffix=options.suffix)
        
        if options.double_encode:
            modifiers.add_modifier(DoubleEncoder)
        
        if options.zero_cap and self.delivers_qualities:
            modifiers.add_modifier(
                ZeroCapper, quality_base=options.quality_base)
        
        if options.trim_primer:
            modifiers.add_modifier(PrimerTrimmer)
        
        if options.merge_overlapping:
            modifiers.add_modifier(
                MergeOverlapping,
                min_overlap=options.merge_min_overlap,
                error_rate=options.merge_error_rate,
                mismatch_action=options.correct_mismatches)
        
        # Create Filters and Formatters
        
        min_affected = 2 if options.pair_filter == 'both' else 1
        filters = Filters(FilterFactory(options.paired, min_affected))
        
        output1 = output2 = None
        interleaved = False
        if options.interleaved_output:
            output1 = options.interleaved_output
            interleaved = True
        else:
            output1 = options.output
            output2 = options.paired_output
        
        seq_formatter_args = dict(
            qualities=self.delivers_qualities,
            colorspace=options.colorspace,
            interleaved=interleaved
        )
        formatters = Formatters(output1, seq_formatter_args)
        force_create = []
            
        if options.merge_overlapping:
            filters.add_filter(MergedReadFilter)
            if options.merged_output:
                formatters.add_seq_formatter(
                    MergedReadFilter, options.merged_output)
            
        if options.minimum_length is not None and options.minimum_length > 0:
            filters.add_filter(TooShortReadFilter, options.minimum_length)
            if options.too_short_output:
                formatters.add_seq_formatter(
                    TooShortReadFilter,
                    options.too_short_output, options.too_short_paired_output)
        
        if options.maximum_length < sys.maxsize:
            filters.add_filter(TooLongReadFilter, options.maximum_length)
            if options.too_long_output is not None:
                formatters.add_seq_formatter(
                    TooLongReadFilter,
                    options.too_long_output, options.too_long_paired_output)
        
        if options.max_n is not None:
            filters.add_filter(NContentFilter, options.max_n)
        
        if options.discard_trimmed:
            filters.add_filter(TrimmedFilter)
        
        if not formatters.multiplexed:
            if output1 is not None:
                formatters.add_seq_formatter(NoFilter, output1, output2)
                if output1 != STDOUT and options.writer_process:
                    force_create.append(output1)
                    if output2 is not None:
                        force_create.append(output2)
            elif not (options.discard_trimmed and options.untrimmed_output):
                formatters.add_seq_formatter(NoFilter, options.default_outfile)
                if options.default_outfile != STDOUT and options.writer_process:
                    force_create.append(options.default_outfile)
        
        if options.discard_untrimmed or options.untrimmed_output:
            filters.add_filter(UntrimmedFilter)
        
        if not options.discard_untrimmed:
            if formatters.multiplexed:
                untrimmed = options.untrimmed_output or output1.format(
                    name='unknown')
                formatters.add_seq_formatter(UntrimmedFilter, untrimmed)
                formatters.add_seq_formatter(NoFilter, untrimmed)
            elif options.untrimmed_output:
                formatters.add_seq_formatter(
                    UntrimmedFilter,
                    options.untrimmed_output, options.untrimmed_paired_output)
        
        if options.rest_file:
            formatters.add_info_formatter(RestFormatter(options.rest_file))
        if options.info_file:
            formatters.add_info_formatter(InfoFormatter(options.info_file))
        if options.wildcard_file:
            formatters.add_info_formatter(
                WildcardFormatter(options.wildcard_file))
        
        if options.paired:
            mixin_class = PairedEndPipelineMixin
        else:
            mixin_class = SingleEndPipelineMixin
        writers = Writers(force_create)
        record_handler = RecordHandler(modifiers, filters, formatters)
        if options.stats:
            record_handler = StatsRecordHandlerWrapper(
                record_handler, options.paired, options.stats,
                qualities=self.delivers_qualities,
                quality_base=self.quality_base)
        
        logger = logging.getLogger()
        num_adapters = sum(len(a) for a in modifiers.get_adapters())
        logger.info(
            "Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
            num_adapters,
            's' if num_adapters > 1 else '', options.error_rate * 100,
            {
                False: 'single-end',
                'first': 'paired-end legacy',
                'both': 'paired-end'
            }[options.paired])
        if (
                options.paired == 'first' and (
                    len(record_handler.modifiers.get_modifiers(read=2)) > 0 or
                    options.quality_cutoff)):
            logger.warning('\n'.join(textwrap.wrap(
                'Requested read modifications are applied only to the '
                'first read since backwards compatibility mode is enabled. '
                'To modify both reads, also use any of the -A/-B/-G/-U '
                'options. Use a dummy adapter sequence when necessary: '
                '-A XXX')))
        
        if options.threads is None:
            # Run single-threaded version
            result_handler = WorkerResultHandler(WriterResultHandler(writers))
            pipeline_class = type(
                'TrimPipelineImpl', (mixin_class, TrimPipeline), {})
            pipeline = pipeline_class(record_handler, result_handler)
            self.summary.update(mode='serial', threads=1)
            return run_interruptible(pipeline, self, raise_on_error=True)
        else:
            # Run multiprocessing version
            self.summary.update(mode='parallel', threads=options.threads)
            return self.run_parallel(record_handler, writers, mixin_class)
    
    def run_parallel(self, record_handler, writers, mixin_class):
        """Parallel implementation of run_atropos. Works as follows:
        
        1. Main thread creates N worker processes (where N is the number of
        threads to be allocated) and (optionally) one writer process.
        2. Main thread loads batches of reads (or read pairs) from input file(s)
        and adds them to a queue (the input queue).
        3. Worker processes take batches from the input queue, process them as
        Atropos normally does, and either add the results to the result queue
        (if using a writer process) or write the results to disk. Each result is
        a dict mapping output file names to strings, where each string is a
        concatenation of reads (with appropriate line endings) to be written.
        A parameter also controls whether data compression is done by the
        workers or the writer.
        4. If using a writer process, it takes results from the result queue and
        writes each string to its corresponding file.
        5. When the main process finishes loading reads from the input file(s),
        it sends a signal to the worker processes that they should complete when
        the input queue is empty. It also singals the writer process how many
        total batches to expect, and the writer process exits after it has
        processed that many batches.
        6. When a worker process completes, it adds a summary of its activity to
        the summary queue.
        7. The main process reads summaries from the summary queue and merges
        them to create the complete summary, which is used to generate the
        report.
        
        There are several possible points of failure:
        
        1. The main process may exit due to an unexpected error, or becuase the
        user forces it to exit (Ctrl-C). In this case, an attempt is made to
        cancel all processes before exiting.
        2. A worker or writer process may exit due to an unknown error. To
        handle this, the main process checks that each process is alive whenver
        it times out writing to the input queue, and again when waiting for
        worker summaries. If a process has died, the program exits with an error
        since some data might have gotten lost.
        3. More commonly, process will time out blocking on reading from or
        writing to a queue. Size limits are used (optionally) for the input and
        result queues to prevent using lots of memory. When few threads are
        allocated, it is most likely that the main and writer processes will
        block, whereas with many threads allocated the workers are most likely
        to block. Also, e.g. in a cluster environment, I/O latency may cause a
        "backup" resulting in frequent blocking of the main and workers
        processes. Finally, also e.g. in a cluster environment, processes may
        suspended for periods of time. Use of a hard timeout period, after which
        processes are forced to exit, is thus undesirable. Instead, parameters
        are provided for the user to tune the batch size and max queue sizes to
        their particular environment. Additionally, a "soft" timeout is used,
        after which log messages are escallated from DEBUG to ERROR level. The
        user can then make the decision of whether or not to kill the program.
        
        Args:
            record_handler: RecordHandler object.
            writers: Writers object.
            mixin_class: Mixin to use for creating pipeline class.
        
        Returns:
            The return code.
        """
        # We do all the multicore imports and class definitions within the
        # run_parallel method to avoid extra work if only running in serial
        # mode.
        from multiprocessing import Queue
        from atropos.commands.multicore import (
            ParallelPipelineMixin, RETRY_INTERVAL)
        from atropos.commands.trim.multicore import (
            Done, Killed, ParallelTrimPipelineRunner, QueueResultHandler,
            CompressingWorkerResultHandler, OrderPreservingWriterResultHandler,
            ResultProcess, WriterManager)
        from atropos.io.compression import can_use_system_compression
        
        # Main process
        
        timeout = max(self.process_timeout, RETRY_INTERVAL)
        threads = self.threads
        
        logging.getLogger().debug(
            "Starting atropos in parallel mode with threads=%d, timeout=%d",
            threads, timeout)
        
        if threads < 2:
            raise ValueError("'threads' must be >= 2")
        
        # Reserve a thread for the writer process if it will be doing the
        # compression and if one is available.
        compression = self.compression
        if compression is None:
            compression = "worker"
            if self.writer_process and can_use_system_compression():
                compression = "writer"
        if compression == "writer" and threads > 2:
            threads -= 1
        
        # Queue by which results are sent from the worker processes to the
        # writer process
        result_queue = Queue(self.result_queue_size)
        writer_manager = None
        
        if self.writer_process:
            if compression == "writer":
                worker_result_handler = WorkerResultHandler(
                    QueueResultHandler(result_queue))
            else:
                worker_result_handler = CompressingWorkerResultHandler(
                    QueueResultHandler(result_queue))
            writer_manager = WriterManager(
                writers, compression, self.preserve_order, result_queue,
                timeout)
        else:
            worker_result_handler = WorkerResultHandler(
                WriterResultHandler(writers, use_suffix=True))
        
        pipeline_class = type(
            'TrimPipelineImpl',
            (ParallelPipelineMixin, mixin_class, TrimPipeline), {})
        pipeline = pipeline_class(record_handler, worker_result_handler)
        runner = ParallelTrimPipelineRunner(
            self, pipeline, threads, writer_manager)
        return runner.run()
