"""Implementation of the 'trim' command.
"""
from collections import defaultdict
import logging
import sys
import textwrap

from atropos.commands import (
    Pipeline, SingleEndPipelineMixin, PairedEndPipelineMixin, create_reader)
from atropos.commands.stats import *
from atropos.adapters import AdapterParser, BACK, load_known_adapters
from atropos.trim.modifiers import *
from atropos.trim.filters import *
from atropos.trim.writers import *
from atropos.io import STDOUT
from atropos.reports import generate_reports
from atropos.util import RandomMatchProbability, Timing, run_interruptible_with_result

class TrimPipeline(Pipeline):
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
    
    def finish(self, worker=None):
        self.result_handler.finish()
    
    def summarize(self, error=None):
        summary = super().summarize()
        if not error:
            summary['trim'] = self.record_handler.summarize()
        return summary

class SerialTrimPipeline(TrimPipeline):
    def __init__(self, record_handler, writers):
        result_handler = WorkerResultHandler(WriterResultHandler(writers))
        super().__init__(record_handler, result_handler)

class ParallelTrimPipeline(TrimPipeline):
    def __init__(self, record_handler, result_handler):
        super().__init__(record_handler, result_handler)
        self.seen_batches = set()
    
    def process_batch(self, batch):
        batch_num, batch_data = batch
        self.seen_batches.add(batch_num)
        super().process_batch(batch_data)
    
    def finish(self, worker=None):
        super().finish(worker)
        logging.getLogger().debug("{} finished; processed {} batches, {} reads".format(
            worker.name, len(self.seen_batches), sum(self.record_counts.values())))
    
    def summarize(self, error=None):
        return (self.seen_batches, super().summarize())

class RecordHandler(object):
    def __init__(self, modifiers, filters, formatters):
        self.modifiers = modifiers
        self.filters = filters
        self.formatters = formatters
    
    def handle_record(self, context, read1, read2=None):
        reads = self.modifiers.modify(read1, read2)
        dest = self.filters.filter(*reads)
        self.formatters.format(context['results'], dest, *reads)
        return (dest, reads)
    
    def summarize(self):
        return dict(
            modifiers=self.modifiers.summarize(),
            filters=self.filters.summarize(),
            formatters=self.formatters.summarize())

class StatsRecordHandlerWrapper(object):
    def __init__(self, record_handler, paired, mode='both', **kwargs):
        self.record_handler = record_handler
        self.paired = paired
        self.pre = self.post = None
        if mode in ('pre', 'both'):
            self.pre = {}
        if mode in ('post', 'both'):
            self.post = {}
        self.stats_kwargs = kwargs
    
    def handle_record(self, context, read1, read2=None):
        if self.pre:
            self.collect(self.pre, context['source'], read1, read2)
        dest, reads = self.record_handler.handle_record(context, read1, read2)
        if self.post:
            if dest not in self.post:
                self.post[dest] = {}
            self.collect(self.post[dest], context['source'], *reads)
        return (dest, reads)
    
    def collect(self, stats, source, read1, read2=None):
        if paired:
            self._collect(stats, source[0], read1)
            self._collect(stats, source[1], read2)
        else:
            self._collect(stats, source, read1)
    
    def _collect(self, stats, source, read):
        if source not in self.stats:
            stats[source] = ReadStatCollector(**self.stats_kwargs)
        self.stats[source].collect(read)
    
    def summarize(self):
        summary = self.record_handler.summarize()
        if self.pre:
            summary['pre'] = self.pre
        if self.post:
            summary['post'] = self.post
        return summary

class ResultHandler(object):
    def start(self, worker):
        pass
    
    def finish(self, total_batches=None):
        pass
    
    def write_result(self, batch_num, result):
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
        """
        Given a dict mapping files to lists of strings,
        join the strings and compress them (if necessary)
        and then return the property formatted result
        dict.
        """
        self.handler.write_result(
            batch_num, dict(
                self.prepare_file(*item)
                for item in result.items()))
    
    def prepare_file(self, path, strings):
        return (path, "".join(strings))

class WriterResultHandler(ResultHandler):
    """
    ResultHandler that writes results to disk.
    """
    def __init__(self, writers, compressed=False, use_suffix=False):
        self.writers = writers
        self.compressed = compressed
        self.use_suffix = use_suffix
    
    def start(self, worker):
        if self.use_suffix:
            self.writers.suffix = ".{}".format(worker.index)
    
    def write_result(self, batch_num, result):
        self.writers.write_result(result, self.compressed)
    
    def finish(self, total_batches=None):
        self.writers.close()

def execute(options):
    reader, input_names, qualities, has_qual_file = create_reader(options)
    
    if options.adapter_max_rmp or options.aligner == 'insert':
        match_probability = RandomMatchProbability()
    
    # Create Adapters
    
    has_adapters1 = options.adapters or options.anywhere or options.front
    has_adapters2 = options.adapters2 or options.anywhere2 or options.front2
    
    adapters1 = adapters2 = []
    if has_adapters1 or has_adapters2:
        adapter_cache = load_known_adapters(options)
        parser_args = dict(
            colorspace=options.colorspace,
            max_error_rate=options.error_rate,
            min_overlap=options.overlap,
            read_wildcards=options.match_read_wildcards,
            adapter_wildcards=options.match_adapter_wildcards,
            indels=options.indels, indel_cost=options.indel_cost,
            cache=adapter_cache
        )
        if options.adapter_max_rmp:
            parser_args['match_probability'] = match_probability
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
    if not adapters1 and not adapters2 and not options.quality_cutoff and \
            options.nextseq_trim is None and \
            options.cut == [] and options.cut2 == [] and \
            options.cut_min == [] and options.cut_min2 == [] and \
            (options.minimum_length is None or options.minimum_length <= 0) and \
            options.maximum_length == sys.maxsize and \
            not has_qual_file and options.max_n is None and not options.trim_n \
            and (not options.paired or options.overwrite_low_quality is None):
        raise ValueError("You need to provide at least one adapter sequence.")
    
    if options.aligner == 'insert' and any(
            not a or len(a) != 1 or a[0].where != BACK
            for a in (adapters1, adapters2)):
        raise ValueError("Insert aligner requires a single 3' adapter for each read")
    
    if options.debug:
        for adapter in adapters1 + adapters2:
            adapter.enable_debug()
    
    if options.paired:
        modifiers = PairedEndModifiers(options.paired)
    else:
        modifiers = SingleEndModifiers()
            
    for op in options.op_order:
        if op == 'W' and options.overwrite_low_quality:
            lowq, highq, window = options.overwrite_low_quality
            modifiers.add_modifier(OverwriteRead,
                worse_read_min_quality=lowq, better_read_min_quality=highq,
                window_size=window, base=options.quality_base)
            
        elif op == 'A' and (adapters1 or adapters2):
            # TODO: generalize this using some kind of factory class
            if options.aligner == 'insert':
                # Use different base probabilities if we're trimming bisulfite data.
                # TODO: this doesn't seem to help things, so commenting it out for now
                #base_probs = dict(p1=0.33, p2=0.67) if options.bisulfite else dict(p1=0.25, p2=0.75)
                modifiers.add_modifier(InsertAdapterCutter,
                    adapter1=adapters1[0], adapter2=adapters2[0], action=options.action,
                    mismatch_action=options.correct_mismatches,
                    max_insert_mismatch_frac=options.insert_match_error_rate,
                    max_adapter_mismatch_frac=options.insert_match_adapter_error_rate,
                    match_probability=match_probability,
                    insert_max_rmp=options.insert_max_rmp)
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
        elif op == 'C' and (options.cut or options.cut2):
            modifiers.add_modifier_pair(UnconditionalCutter,
                dict(lengths=options.cut),
                dict(lengths=options.cut2))
        elif op == 'G' and (options.nextseq_trim is not None):
            modifiers.add_modifier(NextseqQualityTrimmer,
                read=1, cutoff=options.nextseq_trim, base=options.quality_base)
        elif op == 'Q' and options.quality_cutoff:
            modifiers.add_modifier(QualityTrimmer,
                cutoff_front=options.quality_cutoff[0],
                cutoff_back=options.quality_cutoff[1],
                base=options.quality_base)
    
    if options.bisulfite:
        if isinstance(options.bisulfite, str):
            if "non-directional" in options.bisulfite:
                modifiers.add_modifier(NonDirectionalBisulfiteTrimmer,
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
                modifiers.add_modifier(MinCutter, read=1, **(options.bisulfite[0]))
            if len(options.bisulfite) > 1 and options.bisulfite[1]:
                modifiers.add_modifier(MinCutter, read=2, **(options.bisulfite[1]))
    
    if options.trim_n:
        modifiers.add_modifier(NEndTrimmer)
    
    if options.cut_min or options.cut_min2:
        modifiers.add_modifier_pair(MinCutter,
            dict(lengths=options.cut_min),
            dict(lengths=options.cut_min2)
        )
    
    if options.length_tag:
        modifiers.add_modifier(LengthTagModifier, length_tag=options.length_tag)
    
    if options.strip_suffix:
        modifiers.add_modifier(SuffixRemover, suffixes=options.strip_suffix)
    
    if options.prefix or options.suffix:
        modifiers.add_modifier(PrefixSuffixAdder, prefix=options.prefix, suffix=options.suffix)
    
    if options.double_encode:
        modifiers.add_modifier(DoubleEncoder)
    
    if options.zero_cap and qualities:
        modifiers.add_modifier(ZeroCapper, quality_base=options.quality_base)
    
    if options.trim_primer:
        modifiers.add_modifier(PrimerTrimmer)
    
    if options.merge_overlapping:
        modifiers.add_modifier(MergeOverlapping,
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
        qualities=qualities,
        colorspace=options.colorspace,
        interleaved=interleaved
    )
    formatters = Formatters(output1, seq_formatter_args)
    force_create = []
        
    if options.merge_overlapping:
        filters.add_filter(MergedReadFilter)
        if options.merged_output:
            formatters.add_seq_formatter(MergedReadFilter, options.merged_output)
        
    if options.minimum_length is not None and options.minimum_length > 0:
        filters.add_filter(TooShortReadFilter, options.minimum_length)
        if options.too_short_output:
            formatters.add_seq_formatter(TooShortReadFilter,
                options.too_short_output, options.too_short_paired_output)

    if options.maximum_length < sys.maxsize:
        filters.add_filter(TooLongReadFilter, options.maximum_length)
        if options.too_long_output is not None:
            formatters.add_seq_formatter(TooLongReadFilter,
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
            formatters.add_seq_formatter(NoFilter, default_outfile)
            if default_outfile != STDOUT and options.writer_process:
                force_create.append(default_outfile)
    
    if options.discard_untrimmed or options.untrimmed_output:
        filters.add_filter(UntrimmedFilter)

    if not options.discard_untrimmed:
        if formatters.multiplexed:
            untrimmed = options.untrimmed_output or output1.format(name='unknown')
            formatters.add_seq_formatter(UntrimmedFilter, untrimmed)
            formatters.add_seq_formatter(NoFilter, untrimmed)
        elif options.untrimmed_output:
            formatters.add_seq_formatter(UntrimmedFilter,
                options.untrimmed_output, options.untrimmed_paired_output)

    if options.rest_file:
        formatters.add_info_formatter(RestFormatter(options.rest_file))
    if options.info_file:
        formatters.add_info_formatter(InfoFormatter(options.info_file))
    if options.wildcard_file:
        formatters.add_info_formatter(WildcardFormatter(options.wildcard_file))
    
    mixin_class = PairedEndPipelineMixin if options.paired else SingleEndPipelineMixin
    writers = Writers(force_create)
    record_handler = RecordHandler(modifiers, filters, formatters)
    if options.stats:
        record_handler = StatsRecordHandlerWrapper(
            record_handler, options.paired, mode=options.stats,
            qualities=qualities, tile_key_regexp=options.tile_key_regexp)
    
    logger = logging.getLogger()
    num_adapters = sum(len(a) for a in modifiers.get_adapters())
    logger.info("Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
        num_adapters, 's' if num_adapters > 1 else '', options.error_rate * 100,
        { False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end' }[options.paired])
    if options.paired == 'first' and (
            len(record_handler.modifiers.get_modifiers(read=2)) > 0 or
            options.quality_cutoff):
        logger.warning('\n'.join(textwrap.wrap('WARNING: Requested read '
            'modifications are applied only to the first '
            'read since backwards compatibility mode is enabled. '
            'To modify both reads, also use any of the -A/-B/-G/-U options. '
            'Use a dummy adapter sequence when necessary: -A XXX')))
    
    with Timing() as timing:
        if options.threads is None:
            # Run single-threaded version
            pipeline_class = type('TrimPipelineImpl', (SerialTrimPipeline, mixin_class), {})
            pipeline = pipeline_class(record_handler, writers)
            rc, summary = run_interruptible_with_result(pipeline, reader)
        else:
            # Run multiprocessing version
            pipeline_class = type('TrimPipelineImpl', (ParallelTrimPipeline, mixin_class), {})
            rc, summary = run_parallel(
                reader, record_handler, writers, pipeline_class,
                options.threads, options.process_timeout, options.preserve_order,
                options.read_queue_size, options.result_queue_size,
                options.writer_process, options.compression)
    
    reader.close()
    
    # For trim stats, any value with a name that starts with 'records_' will
    # have 'fraction_<var>' computed as value / total_reads, and any value with a
    # name that starts with 'bp_' will have 'fraction_<var>' and 'total_<var>'
    # computed.
    
    if rc == 0:
        summary['timing'] = timing.summarize()
        
        record_counts = summary['record_counts']
        total_records = sum(record_counts.values())
        
        bp_counts = summary['bp_counts']
        total_bp = [sum(counts) for counts in zip(*bp_counts.values())]
        
        def _recurse(d):
            for key, value in tuple(d.items()):
                if isinstance(value, dict):
                    _recurse(value)
                elif isinstance(key, str):
                    if key.startswith('records_'):
                        d["fraction_{}".format(key)] = (
                            (value / total_records) if total_records != 0 else 0)
                    elif key.startswith('bp_'):
                        print("{}=<{}>".format(key,value))
                        d["fraction_{}".format(key)] = [
                            (v / b) if b != 0 else 0
                            for v, b in zip(value, total_bp)]
                        d["total_{}".format(key)] = sum(value)
        
        _recurse(summary['trim'])
        
        generate_reports(options, summary)
    
    return (rc, None, summary)

# stats["written_fraction"] = 0
# stats["too_short_fraction"] = 0
# stats["too_long_fraction"] = 0
# stats["too_many_n_fraction"] = 0
# stats["with_adapters_fraction"] = [0, 0]
# stats["total_written_bp_fraction"] = 0.0
# for modifier_class in self.trimmer_classes:
#     name = modifier_class.__name__
#     stats["{}_fraction".format(name)] = 0.0
#
# N = stats["N"]
# if N > 0:
#     stats["written_fraction"] = stats["written"] / N if stats["written"] else 0
#     stats["too_short_fraction"] = stats["too_short"] / N if stats["too_short"] else 0
#     stats["too_long_fraction"] = stats["too_long"] / N if stats["too_long"] else 0
#     stats["too_many_n_fraction"] = stats["too_many_n"] / N if stats["too_many_n"] else 0
#     stats["with_adapters_fraction"] = [ (v / N) for v in stats["with_adapters"] ]
#     if "corrected" in stats:
#         stats["corrected_fraction"] = stats["corrected"] / N
#
# if stats["total_bp"] > 0:
#     N = stats["total_bp"]
#     stats["total_written_bp_fraction"] = (stats["total_written_bp"] / N) if stats["total_written_bp"] else 0
#     if "corrected" in stats:
#         stats["corrected_bp_fraction"] = [ (c / N) for c in stats["corrected_bp"] ]
#         stats["total_corrected_bp_fraction"] = stats["total_corrected_bp"] / N
#     for modifier_class in self.trimmer_classes:
#         name = modifier_class.__name__
#         if stats[name]:
#             stats["{}_fraction".format(name)] = (stats[name] / N)

def run_parallel(
        reader, record_handler, writers, pipeline_class, threads=2, timeout=30,
        preserve_order=False, input_queue_size=0, result_queue_size=0,
        use_writer_process=True, compression=None):
    """Parallel implementation of run_atropos.
    
    Works as follows:
    
    1. Main thread creates N worker processes (where N is the number of threads
    to be allocated) and (optionally) one writer process.
    2. Main thread loads batches of reads (or read pairs) from input file(s) and
    adds them to a queue (the input queue).
    3. Worker processes take batches from the input queue, process them as
    Atropos normally does, and either add the results to the result queue (if
    using a writer process) or write the results to disk. Each result is a dict
    mapping output file names to strings, where each string is a concatenation
    of reads (with appropriate line endings) to be written. A parameter also
    controls whether data compression is done by the workers or the writer.
    4. If using a writer process, it takes results from the result queue and
    writes each string to its corresponding file.
    5. When the main process finishes loading reads from the input file(s), it
    sends a signal to the worker processes that they should complete when the
    input queue is empty. It also singals the writer process how many total
    batches to expect, and the writer process exits after it has processed that
    many batches.
    6. When a worker process completes, it adds a summary of its activity to the
    summary queue.
    7. The main process reads summaries from the summary queue and merges them
    to create the complete summary, which is used to generate the report.
    
    There are several possible points of failure:
    
    1. The main process may exit due to an unexpected error, or becuase the user
    forces it to exit (Ctrl-C). In this case, an attempt is made to cancel all
    processes before exiting.
    2. A worker or writer process may exit due to an unknown error. To handle
    this, the main process checks that each process is alive whenver it times
    out writing to the input queue, and again when waiting for worker summaries.
    If a process has died, the program exits with an error since some data might
    have gotten lost.
    3. More commonly, process will time out blocking on reading from or writing
    to a queue. Size limits are used (optionally) for the input and result
    queues to prevent using lots of memory. When few threads are allocated, it
    is most likely that the main and writer processes will block, whereas with
    many threads allocated the workers are most likely to block. Also, e.g. in a
    cluster environment, I/O latency may cause a "backup" resulting in frequent
    blocking of the main and workers processes. Finally, also e.g. in a cluster
    environment, processes may suspended for periods of time. Use of a hard
    timeout period, after which processes are forced to exit, is thus
    undesirable. Instead, parameters are provided for the user to tune the batch
    size and max queue sizes to their particular environment. Additionally, a
    "soft" timeout is used, after which log messages are escallated from DEBUG
    to ERROR level. The user can then make the decision of whether or not to
    kill the program.
    
    Args:
        reader: Iterator over batches of reads (most likely a BatchIterator).
        record_handler: RecordHandler object.
        writers: Writers object.
        pipeline_class: Class of pipeline to instantiate.
        threads: Number of worker threads to use; additional threads are used
            for the main proccess and the writer process (if requested).
        timeout: Number of seconds after which waiting processes escalate their
            messages from DEBUG to ERROR.
        preserve_order: Whether to preserve the input order of reads when
            writing (only valid when `use_writer_process=True`)
        input_queue_size: Max number of items that can be in the input queue,
            or 0 for no limit (be warned that this could explode memory usage)
        result_queue_size: Max number of items that can be in the result queue,
            or 0 for no limit (be warned that this could explode memory usage)
        use_writer_process: If True, a separate thread will be used to write
            results to disk. Otherwise, each worker thread will write its
            results to an output file with a '.N' extension, where N is the
            thread index. This is useful in cases where the I/O is the main
            bottleneck.
        compression: If "writer", the writer process perform data compression,
            otherwise the worker processes performs compression.
    
    Returns:
        Tuple of (return_code, {summary}).
    """
    # We do all the multicore imports and class definitions within the
    # run_parallel method to avoid extra work if only running in serial mode.
    from multiprocessing import (
        Process, Queue, MulticoreError, launch_workers, ensure_processes,
        wait_on, enqueue, enqueue_all, dequeue)
    from atropos.commands.multicore import Control, PendingQueue
    from atropos.io.compression import get_compressor, can_use_system_compression
    
    class Done(MulticoreError):
        pass
    
    class Killed(MulticoreError):
        pass
    
    class QueueResultHandler(ResultHandler):
        """ResultHandler that writes results to the output queue.
        """
        def __init__(self, queue):
            self.queue = queue
        
        def start(self, worker):
            self.message = "{} waiting to queue result {{}}".format(worker.name)
            self.timeout = worker.timeout
        
        def write_result(self, batch_num, result):
            enqueue(
                self.queue,
                (batch_num, result),
                wait_message=self.message,
                timeout=self.timeout)

    class CompressingWorkerResultHandler(WorkerResultHandler):
        """Wraps a ResultHandler and compresses results prior to writing.
        """
        def start(self, worker):
            super().start(worker)
            self.file_compressors = {}
        
        def prepare_file(self, path, strings):
            compressor = self.get_compressor(path)
            if compressor:
                return ((path, 'wb'), compressor.compress(b''.join(
                    s.encode() for s in strings)))
            else:
                return ((path, 'wt'), "".join(strings))
        
        def get_compressor(self, filename):
            if filename not in self.file_compressors:
                self.file_compressors[filename] = get_compressor(filename)
            return self.file_compressors[filename]

    class OrderPreservingWriterResultHandler(WriterResultHandler):
        """
        Writer thread that is less time/memory efficient, but is
        guaranteed to preserve the original order of records.
        """
        def start(self, worker):
            super().__init__(worker)
            self.pending = PendingQueue()
            self.cur_batch = 1
        
        def write_result(self, batch_num, result):
            if batch_num == self.cur_batch:
                self.writers.write_result(result, self.compressed)
                self.cur_batch += 1
                self.consume_pending()
            else:
                self.pending.push(batch_num, result)
        
        def finish(self, total_batches):
            if total_batches is not None:
                self.consume_pending()
                if self.cur_batch != total_batches:
                    raise AtroposError(
                        "OrderPreservingWriterResultHandler finishing without "
                        "having seen {} batches".format(total_batches))
            self.writers.close()
        
        def consume_pending(self):
            while (not self.pending.empty) and (self.cur_batch == pending.min_priority):
                self.writers.write_result(pending.pop(), self.compressed)
                self.cur_batch += 1
    
    class ResultProcess(Process):
        """Thread that accepts results from the worker threads and process them
        using a ResultHandler. Each batch is expected to be
        (batch_num, path, records), where path is the destination file and records
        is a string. Not guaranteed to preserve the original order of sequence records.

        result_handler: A ResultHandler object.
        queue: Input queue.
        control: A shared value for communcation with the main process.
        timeout: Seconds to wait for next batch before complaining.
        """
        def __init__(self, result_handler, queue, control, timeout=60):
            super().__init__(name="Result process")
            self.result_handler = result_handler
            self.queue = queue
            self.control = control
            self.timeout = timeout
            self.seen_batches = set()
            self.num_batches = None
        
        def run(self):
            logging.getLogger().debug(
                "Writer process running under pid {}".format(
                    self.name, os.getpid()))
            
            def fail_callback():
                if self.num_batches is None and self.control.check_value_positive():
                    self.num_batches = self.control.get_value()
                if self.num_batches is not None and len(self.seen_batches) >= self.num_batches:
                    raise Done()
            
            def timeout_callback():
                if self.num_batches is not None:
                    missing = set(range(1, self.num_batches+1)) - self.seen_batches
                    logging.getLogger().error(
                        "Result thread still missing batches {} of {}".format(
                            ",".join(str(i) for i in missing), self.num_batches))
            
            def iter_batches():
                while True:
                    batch = dequeue(
                        self.queue,
                        wait_message="Result process waiting on result {}",
                        timeout=self.timeout,
                        fail_callback=fail_callback,
                        timeout_callback=timeout_callback)
                    yield batch

            try:
                self.result_handler.start(self)
                
                for batch_num, result in iter_batches():
                    self.seen_batches.add(batch_num)
                    self.result_handler.write_result(batch_num, result)
            except Done:
                logging.getLogger().debug("Writer process exiting normally")
            except Killed:
                logging.getLogger().debug("Writer process exited early")
            except:
                logging.getLogger().error(
                    "Unexpected error in writer process", exc_info=True)
                self.control.set_value(CONTROL_ERROR)
            finally:
                num_batches = self.control.get_value(lock=True)
                self.result_handler.finish(num_batches if num_batches > 0 else None)
    
    # Main process
    
    logging.getLogger().debug(
        "Starting atropos in parallel mode with threads={}, timeout={}".format(
            threads, timeout))
    
    if threads < 2:
        raise ValueError("'threads' must be >= 2")
    
    # Reserve a thread for the writer process if it will be doing the compression and if one is available.
    if compression is None:
        compression = "writer" if use_writer_process and can_use_system_compression() else "worker"
    if compression == "writer" and threads > 2:
        threads -= 1
    
    timeout = max(timeout, RETRY_INTERVAL)
    
    # Queue by which batches of reads are sent to worker processes
    input_queue = Queue(input_queue_size)
    # Queue by which results are sent from the worker processes to the writer process
    result_queue = Queue(result_queue_size)
    # Queue for processes to send summary information back to main process
    summary_queue = Queue(threads)
    # Aggregate summary
    summary = MergingDict()
    
    if use_writer_process:
        worker_result_handler = QueueResultHandler(result_queue)
        if compression == "writer":
            worker_result_handler = WorkerResultHandler(worker_result_handler)
        else:
            worker_result_handler = CompressingWorkerResultHandler(worker_result_handler)
        
        # Shared variable for communicating with writer thread
        writer_control = Control(CONTROL_ACTIVE)
        
        # result handler
        if preserve_order:
            writer_result_handler = OrderPreservingWriterResultHandler(
                writers, compressed=compression == "worker")
        else:
            writer_result_handler = WriterResultHandler(
                writers, compressed=compression == "worker")
        
        # writer process
        writer_process = ResultProcess(
            writer_result_handler, result_queue, writer_control, timeout)
        writer_process.start()
    else:
        worker_result_handler = WorkerResultHandler(
            WriterResultHandler(writers, use_suffix=True))
    
    # Start worker processes, reserve a thread for the reader process,
    # which we will get back after it completes
    pipeline = pipeline_class(record_handler, worker_result_handler)
    worker_args = (input_queue, pipeline, summary_queue, timeout)
    worker_processes = launch_workers(threads - 1, worker_args)
    
    def ensure_alive():
        ensure_processes(worker_processes)
        if (use_writer_process and not (
                writer_process.is_alive() and
                writer_control.check_value(CONTROL_ACTIVE))):
            raise AtroposError("Writer process exited")
    
    def _run(worker_processes):
        # Add batches of reads to the input queue. Provide a timeout callback
        # to check that subprocesses are alive.
        num_batches = enqueue_all(reader, input_queue, timeout, ensure_alive)
        logging.getLogger().debug(
            "Main loop complete; saw {} batches".format(num_batches))
        
        # Tell the worker processes no more input is coming
        enqueue_all((None,) * threads, input_queue, timeout, ensure_alive)
        
        # Tell the writer thread the max number of batches to expect
        if use_writer_process:
            writer_control.set_value(num_batches)
        
        # Now that the reader process is done, it essentially
        # frees up another thread to use for a worker
        worker_processes.extend(launch_workers(1, worker_args, offset=threads-1))
        
        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            try:
                ensure_processes(worker_processes,
                    "Workers are still alive and haven't returned summaries: {}",
                    alive=False)
            except Exception as e:
                logging.getLogger().error(e)
            
        wait_on(
            lambda: summary_queue.full(),
            wait_message="Waiting on worker summaries {}",
            timeout=timeout,
            wait=True,
            timeout_callback=summary_timeout_callback)
        
        # Process summary information from worker processes
        logging.getLogger().debug(
            "Processing summary information from worker processes")
        seen_summaries = set()
        seen_batches = set()
        
        def summary_fail_callback():
            missing_summaries = set(range(1, threads)) - seen_summaries
            raise AtroposError("Missing summaries from processes {}".format(
                ",".join(str(s) for s in missing)))
        
        for i in range(1, threads+1):
            batch = dequeue(summary_queue, fail_callback=summary_fail_callback)
            worker_index, worker_batches, stats = batch
            if stats is None:
                raise AtroposError(
                    "Worker process {} died unexpectedly".format(worker_index))
            else:
                logging.getLogger().debug(
                    "Processing summary for worker {}".format(worker_index))
            seen_summaries.add(worker_index)
            seen_batches |= worker_batches
            summary.merge(stats)
        
        # Check if any batches were missed
        if num_batches > 0:
            missing_batches = set(range(1, num_batches+1)) - seen_batches
            if len(missing_batches) > 0:
                raise AtroposError("Workers did not process batches {}".format(
                    ",".join(str(b) for b in missing_batches)))
        
        if use_writer_process:
            # Wait for writer to complete
            wait_on_process(writer_process, timeout)
    
    try:
        rc = run_interruptible_with_result(_run, worker_processes)
    finally:
        # notify all threads that they should stop
        logging.getLogger().debug("Exiting all processes")
        def kill(process):
            if rc <= 1:
                wait_on_process(process, timeout, terminate=True)
            elif process.is_alive():
                process.terminate()
        for process in worker_processes:
            kill(process)
        if use_writer_process:
            kill(writer_process)
    
    return (rc, summary)
