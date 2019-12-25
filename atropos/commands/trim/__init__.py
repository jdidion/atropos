from abc import ABCMeta, abstractmethod
from collections import Sequence as SequenceCollection
from multiprocessing import Process
from pathlib import Path
import sys
from typing import (
    Dict,
    Iterable,
    Optional,
    Sequence as SequenceType,
    Tuple,
    Type,
    Union,
)

from loguru import logger
from xphyle import STDOUT, STDERR

from atropos.adapters import AdapterParser, AdapterType
from atropos.commands import (
    BaseCommand,
    Summary,
    Pipeline,
    SingleEndPipelineMixin,
    PairedEndPipelineMixin,
)
from atropos.commands.multicore import WorkerProcess
from atropos.commands.metrics import (
    MetricsMode,
    SingleEndReadMetrics,
    PairedEndReadMetrics,
)
from atropos.commands.trim.modifiers import (
    AdapterCutter,
    DoubleEncoder,
    InsertAdapterCutter,
    LengthTagModifier,
    MergeOverlapping,
    MinCutter,
    NEndTrimmer,
    TwoColorQualityTrimmer,
    UmiTrimmer,
    SyncUmi,
    AddUmi,
    NonDirectionalBisulfiteTrimmer,
    OverwriteRead,
    PairedEndModifiers,
    PrefixSuffixAdder,
    PrimerTrimmer,
    QualityTrimmer,
    RRBSTrimmer,
    SingleEndModifiers,
    SuffixRemover,
    SwiftBisulfiteTrimmer,
    UnconditionalCutter,
    ZeroCapper,
    Modifiers,
)
from atropos.commands.trim.filters import (
    Filter,
    FilterFactory,
    Filters,
    MergedReadFilter,
    NContentFilter,
    NoFilter,
    TooLongReadFilter,
    TooShortReadFilter,
    TrimmedFilter,
    UntrimmedFilter,
)
from atropos.commands.trim.writers import (
    Formatters,
    InfoFormatter,
    RestFormatter,
    WildcardFormatter,
    Writers,
)
from atropos.io.sequence import Sequence
from atropos.utils import ReturnCode, classproperty, run_interruptible
from atropos.utils.argparse import Namespace
from atropos.utils.collections import Summarizable
from atropos.utils.statistics import RandomMatchProbability


class RecordHandler(Summarizable):
    """
    Base class for record handlers.
    """

    def __init__(self, modifiers: Modifiers, filters: Filters, formatters: Formatters):
        self.modifiers = modifiers
        self.filters = filters
        self.formatters = formatters

    def handle_record(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ) -> Tuple[Type[Filter], Tuple[Sequence, Optional[Sequence]]]:
        """
        Handles a pair of reads.
        """
        reads = self.modifiers.modify(read1, read2)
        dest = self.filters.filter(*reads)
        self.formatters.format(context["results"], dest, *reads)
        return dest, reads

    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """
        return dict(
            trim=dict(
                modifiers=self.modifiers.summarize(),
                filters=self.filters.summarize(),
                formatters=self.formatters.summarize(),
            )
        )


class MetricsRecordHandlerWrapper(Summarizable):
    """
    Wrapper around a record handler that collects read metrics before and/or after
    trimming.
    """

    def __init__(
        self,
        record_handler: RecordHandler,
        paired: bool,
        metrics_args: Dict[MetricsMode, dict],
        **kwargs,
    ):
        """
        Args:
            record_handler:
            paired: Whether reads are paired-end.
            metrics_args: Sequence; when to collect metrics ('pre', 'post')
            kwargs:
                mode: Collection mode; pre=only collect pre-trim metrics; post=only
                collect post-trim metrics; both=collect both pre- and post-trim metrics.
        """
        self.record_handler = record_handler
        self.read_metrics_class = (
            PairedEndReadMetrics if paired else SingleEndReadMetrics
        )

        if MetricsMode.PRE in metrics_args:
            self.pre = {}
            self.pre_kwargs = kwargs.copy()
            self.pre_kwargs.update(metrics_args[MetricsMode.PRE])
        else:
            self.pre = None

        if MetricsMode.POST in metrics_args:
            self.post = {}
            self.post_kwargs = kwargs.copy()
            self.post_kwargs.update(metrics_args[MetricsMode.POST])
        else:
            self.post = None

    def handle_record(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ) -> Tuple[Type[Filter], Tuple[Sequence, Optional[Sequence]]]:
        """
        Handles a pair of reads.

        Args:
            context:
            read1:
            read2:
        """
        if self.pre is not None:
            self.collect(self.pre, context["source"], read1, read2, **self.pre_kwargs)
        dest, reads = self.record_handler.handle_record(context, read1, read2)
        if self.post is not None:
            if dest not in self.post:
                self.post[dest] = {}
            self.collect(self.post[dest], context["source"], *reads, **self.post_kwargs)
        return dest, reads

    def collect(
        self,
        metrics: dict,
        source: Union[str, Tuple[str, str]],
        read1: Sequence,
        read2: Optional[Sequence] = None,
        **kwargs,
    ):
        """
        Collects metrics on a pair of reads.

        Args:
            metrics: The :class:`ReadStatistics` object.
            source: The source file(s).
            read1
            read2: The reads.
            kwargs:
        """
        if source not in metrics:
            metrics[source] = self.read_metrics_class(**kwargs)

        metrics[source].collect(read1, read2)

    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """
        summary = self.record_handler.summarize()

        if self.pre is not None:
            summary["pre"] = dict(
                (source, metrics.summarize()) for source, metrics in self.pre.items()
            )

        if self.post is not None:
            summary["post"] = {}
            for dest, metrics_dict in self.post.items():
                summary["post"][dest.name] = dict(
                    (source, metrics.summarize())
                    for source, metrics in metrics_dict.items()
                )

        return summary


class ResultHandler(metaclass=ABCMeta):
    """
    Base class for result handlers.
    """

    def start(self, worker: Optional[Process] = None):
        """
        Starts the result handler.
        """

    def finish(self, total_batches: Optional[int] = None):
        """
        Finishes the result handler.

        Args:
            total_batches: Total number of batches processed.
        """

    @abstractmethod
    def write_result(self, batch_num: int, result: dict):
        """
        Writes a batch of results to output.

        Args:
            batch_num: The batch number.
            result: The result to write.
        """


class ResultHandlerWrapper(ResultHandler):
    """
    Wraps a ResultHandler.
    """

    def __init__(self, handler: ResultHandler):
        self.handler = handler

    def start(self, worker: Optional[Process] = None):
        self.handler.start(worker)

    def write_result(self, batch_num: int, result: dict):
        self.handler.write_result(batch_num, result)

    def finish(self, total_batches: Optional[int] = None):
        self.handler.finish(total_batches=total_batches)


class WorkerResultHandler(ResultHandlerWrapper):
    """
    Wraps a ResultHandler and compresses results prior to writing.
    """

    def write_result(self, batch_num: int, result: dict):
        """
        Given a dict mapping files to lists of strings, joins the strings and
        compresses them (if necessary) and then returns the properly formatted
        result dict.
        """
        self.handler.write_result(
            batch_num, dict(self.prepare_file(*item) for item in result.items())
        )

    def prepare_file(self, path: str, strings: Iterable[str]):
        """
        Prepares data for writing.

        Returns:
            Tuple (path, data).
        """
        return path, "".join(strings)


class WriterResultHandler(ResultHandler):
    """
    ResultHandler that writes results to disk.
    """

    def __init__(
        self, writers: Writers, compressed: bool = False, use_suffix: bool = False
    ):
        """
        Args:
            writers: :class:`Writers` object.
            compressed: Whether the data is compressed.
            use_suffix: Whether to add the worker index as a file suffix. Used for
                parallel-write mode.
        """
        self.writers = writers
        self.compressed = compressed
        self.use_suffix = use_suffix

    def start(self, worker: Optional[WorkerProcess] = None):
        if self.use_suffix:
            if worker is None:
                raise ValueError("worker must not be None")

            self.writers.suffix = f".{worker.index}"

    def write_result(self, batch_num: int, result: dict):
        self.writers.write_result(result, self.compressed)

    def finish(self, total_batches: Optional[int] = None):
        self.writers.close()


class TrimPipeline(Pipeline, metaclass=ABCMeta):
    """
    Base trimming pipeline.

    Args:
        record_handler:
        result_handler:
    """

    def __init__(self, record_handler: RecordHandler, result_handler: ResultHandler):
        super().__init__()
        self.record_handler = record_handler
        self.result_handler = result_handler

    def start(self, worker: Optional[Process] = None):
        self.result_handler.start(worker)

    def add_to_context(self, context: dict):
        context["results"] = {}

    def handle_records(self, context: dict, records: SequenceType[Sequence]):
        super().handle_records(context, records)
        self.result_handler.write_result(context["index"], context["results"])

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        return self.record_handler.handle_record(context, read1, read2)

    def finish(self, summary: Summary, **kwargs):
        self.result_handler.finish()
        super().finish(summary)
        summary.update(self.record_handler.summarize())


class TrimSummary(Summary):
    """
    Summary that adds aggregate values for record and bp metrics.
    """

    def _post_process_other(self, dict_val: dict, key, value) -> None:
        """
        For trim metrics, any value with a name that starts with 'records_' will have
        'fraction_<var>' computed as value / total_records, and any value with a name
        that starts with 'bp_' will have 'fraction_<var>' and 'total_<var>' computed.
        We also replace any `Const`s with their values.
        """
        if self.has_exception:
            return

        def frac(val: int, _total: int):
            """
            Computes fraction of total.
            """
            return (val / _total) if val and _total != 0 else 0

        if isinstance(key, str):
            if key.startswith("records_"):
                frac_key = f"fraction_{key}"
                total_records = self["total_record_count"]

                if isinstance(value, SequenceCollection):
                    dict_val[frac_key] = [frac(val, total_records) for val in value]
                    total = sum(val for val in value if val)
                    dict_val[f"total_{key}"] = total
                else:
                    dict_val[frac_key] = frac(value, total_records)
            elif key.startswith("bp_"):
                frac_key = f"fraction_{key}"
                sum_total_bp = self["sum_total_bp_count"]

                if isinstance(value, SequenceCollection):
                    dict_val[frac_key] = [
                        frac(val, bps)
                        for val, bps in zip(value, self["total_bp_counts"])
                    ]
                    total = sum(val for val in value if val)
                    dict_val[f"total_{key}"] = total
                    dict_val[f"fraction_total_{key}"] = frac(
                        total, sum_total_bp
                    )
                else:
                    dict_val[frac_key] = frac(value, sum_total_bp)


class TrimCommand(BaseCommand):
    @classproperty
    def name(cls) -> str:
        return "trim"

    @classmethod
    def _create_summary(cls) -> dict:
        return TrimSummary()

    def __call__(self) -> ReturnCode:
        options = self.options
        match_probability = RandomMatchProbability()

        # Create Adapters
        adapters1 = []
        adapters2 = []
        has_adapters1 = bool(options.adapters or options.anywhere or options.front)
        has_adapters2 = bool(options.adapters2 or options.anywhere2 or options.front2)

        if has_adapters1 or has_adapters2:
            adapter_cache = self._load_known_adapters()

            parser_args = dict(
                colorspace=options.colorspace,
                max_error_rate=options.error_rate,
                min_overlap=options.overlap,
                read_wildcards=options.match_read_wildcards,
                adapter_wildcards=options.match_adapter_wildcards,
                indels=options.indels,
                indel_cost=options.indel_cost,
                cache=adapter_cache,
                gc_content=options.gc_content,
                match_probability=match_probability,
                alphabet=options.alphabet,
            )
            if options.adapter_max_rmp:
                parser_args["max_rmp"] = options.adapter_max_rmp

            adapter_parser = AdapterParser(**parser_args)

            if has_adapters1:
                adapters1 = adapter_parser.parse_multi(
                    options.adapters, options.anywhere, options.front
                )

            if has_adapters2:
                adapters2 = adapter_parser.parse_multi(
                    options.adapters2, options.anywhere2, options.front2
                )

            if options.cache_adapters:
                adapter_cache.save()

        # Create Modifiers
        # TODO: can this be replaced with an argparse required group?
        if (
            not adapters1
            and not adapters2
            and not options.quality_cutoff
            and options.twocolor_trim is None
            and options.cut == []
            and options.cut2 == []
            and options.cut_min == []
            and options.cut_min2 == []
            and (options.minimum_length is None or options.minimum_length <= 0)
            and options.maximum_length == sys.maxsize
            and not options.trim_n
            and not self.get_option("has_qualfile")
            and options.max_n is None
            and (not options.paired or options.overwrite_low_quality is None)
        ):
            raise ValueError("You need to provide at least one adapter sequence.")

        if options.aligner == "insert" and any(
            not a or len(a) != 1 or a[0].where != AdapterType.BACK
            for a in (adapters1, adapters2)
        ):
            raise ValueError(
                "Insert aligner requires a single 3' adapter for each read"
            )

        if options.debug:
            for adapter in adapters1 + adapters2:
                adapter.enable_debug()

        if options.paired:
            modifiers = PairedEndModifiers(options.paired)
        else:
            modifiers = SingleEndModifiers()

        if options.read1_umi or options.read2_umi:
            modifiers.add_modifier_pair(
                UmiTrimmer,
                dict(number_of_bases=options.read1_umi),
                dict(number_of_bases=options.read2_umi),
            )
            if options.paired:
                modifiers.add_modifier(SyncUmi, delim=options.umi_delim)
            else:
                modifiers.add_modifier(AddUmi, delim=options.umi_delim)

        for oper in options.op_order:
            if oper == "W" and options.overwrite_low_quality:
                lowq, highq, window = options.overwrite_low_quality
                modifiers.add_modifier(
                    OverwriteRead,
                    worse_read_min_quality=lowq,
                    better_read_min_quality=highq,
                    window_size=window,
                    base=options.quality_base,
                )
            elif oper == "A" and (adapters1 or adapters2):
                # TODO: generalize this using some kind of factory class
                if options.aligner == "insert":
                    # Use different base probabilities if we're trimming
                    # bisulfite data.
                    # TODO: this doesn't seem to help things, so commenting it
                    #  out for now
                    # if options.bisulfite:
                    #   base_probs = dict(match_prob=0.33, mismatch_prob=0.67)
                    # else:
                    #   base_probs = dict(match_prob=0.25, mismatch_prob=0.75)
                    modifiers.add_modifier(
                        InsertAdapterCutter,
                        adapter1=adapters1[0],
                        adapter2=adapters2[0],
                        action=options.action,
                        mismatch_action=options.correct_mismatches,
                        max_insert_mismatch_frac=options.insert_match_error_rate,
                        max_adapter_mismatch_frac=
                        options.insert_match_adapter_error_rate,
                        match_probability=match_probability,
                        insert_max_rmp=options.insert_max_rmp,
                        read_wildcards=options.match_read_wildcards,
                        adapter_wildcards=options.match_adapter_wildcards,
                    )
                else:
                    a1_args = (
                        dict(
                            adapters=adapters1,
                            times=options.times,
                            action=options.action,
                        )
                        if adapters1
                        else None
                    )
                    a2_args = (
                        dict(
                            adapters=adapters2,
                            times=options.times,
                            action=options.action,
                        )
                        if adapters2
                        else None
                    )
                    modifiers.add_modifier_pair(AdapterCutter, a1_args, a2_args)
            elif oper == "C" and (options.cut or options.cut2):
                modifiers.add_modifier_pair(
                    UnconditionalCutter,
                    dict(lengths=options.cut),
                    dict(lengths=options.cut2),
                )
            elif oper == "G" and (options.twocolor_trim is not None):
                modifiers.add_modifier(
                    TwoColorQualityTrimmer,
                    cutoff=options.twocolor_trim,
                    base=options.quality_base,
                )
            elif oper == "Q" and options.quality_cutoff:
                modifiers.add_modifier(
                    QualityTrimmer,
                    cutoff_front=options.quality_cutoff[0],
                    cutoff_back=options.quality_cutoff[1],
                    base=options.quality_base,
                )

        if options.bisulfite:
            if isinstance(options.bisulfite, str):
                if "non-directional" in options.bisulfite:
                    modifiers.add_modifier(
                        NonDirectionalBisulfiteTrimmer,
                        rrbs=options.bisulfite == "non-directional-rrbs",
                    )
                elif options.bisulfite == "rrbs":
                    modifiers.add_modifier(RRBSTrimmer)
                elif options.bisulfite in ("epignome", "truseq"):
                    # Trimming leads to worse results
                    # modifiers.add_modifier(TruSeqBisulfiteTrimmer)
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
            modifiers.add_modifier_pair(
                MinCutter, dict(lengths=options.cut_min), dict(lengths=options.cut_min2)
            )

        if options.length_tag:
            modifiers.add_modifier(LengthTagModifier, length_tag=options.length_tag)

        if options.strip_suffix:
            modifiers.add_modifier(SuffixRemover, suffixes=options.strip_suffix)

        if options.prefix or options.suffix:
            modifiers.add_modifier(
                PrefixSuffixAdder, prefix=options.prefix, suffix=options.suffix
            )

        if options.double_encode:
            modifiers.add_modifier(DoubleEncoder)

        delivers_qualities = self.get_option("delivers_qualities")

        if options.zero_cap and delivers_qualities:
            modifiers.add_modifier(ZeroCapper, quality_base=options.quality_base)

        if options.trim_primer:
            modifiers.add_modifier(PrimerTrimmer)

        if options.merge_overlapping:
            modifiers.add_modifier(
                MergeOverlapping,
                min_overlap=options.merge_min_overlap,
                error_rate=options.merge_error_rate,
                mismatch_action=options.correct_mismatches,
            )

        # Create Formatters

        interleaved = options.use_interleaved_output

        if interleaved:
            output1 = options.interleaved_output
            output2 = None
        else:
            output1 = options.output
            output2 = options.paired_output

        if options.output_format is None:
            if delivers_qualities:
                options.output_format = "fastq"
            else:
                options.output_format = "fasta"

        formatters = Formatters(
            output1,
            dict(
                file_format=options.output_format,
                qualities=delivers_qualities,
                colorspace=options.colorspace,
                interleaved=interleaved,
            ),
        )

        # Create filters
        min_affected = 2 if options.pair_filter == "both" else 1
        filters = Filters(FilterFactory(options.paired, min_affected))
        force_create = []

        if options.merge_overlapping:
            filters.add_filter(MergedReadFilter)
            if options.merged_output:
                formatters.add_seq_formatter(MergedReadFilter, options.merged_output)

        if options.minimum_length is not None and options.minimum_length > 0:
            filters.add_filter(TooShortReadFilter, options.minimum_length)
            if options.too_short_output:
                formatters.add_seq_formatter(
                    TooShortReadFilter,
                    options.too_short_output,
                    options.too_short_paired_output,
                )

        if options.maximum_length < sys.maxsize:
            filters.add_filter(TooLongReadFilter, options.maximum_length)
            if options.too_long_output is not None:
                formatters.add_seq_formatter(
                    TooLongReadFilter,
                    options.too_long_output,
                    options.too_long_paired_output,
                )

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
                if (
                    options.default_outfile not in (STDOUT, STDERR) and
                    options.writer_process
                ):
                    force_create.append(options.default_outfile)

        if options.discard_untrimmed or options.untrimmed_output:
            filters.add_filter(UntrimmedFilter)

        if not options.discard_untrimmed:
            if formatters.multiplexed:
                untrimmed = (
                    options.untrimmed_output or
                    Path(str(output1).format(name="unknown"))
                )
                formatters.add_seq_formatter(UntrimmedFilter, untrimmed)
                formatters.add_seq_formatter(NoFilter, untrimmed)
            elif options.untrimmed_output:
                formatters.add_seq_formatter(
                    UntrimmedFilter,
                    options.untrimmed_output,
                    options.untrimmed_paired_output,
                )

        if options.rest_file:
            formatters.add_detail_formatter(RestFormatter(options.rest_file))

        if options.info_file:
            formatters.add_detail_formatter(InfoFormatter(options.info_file))

        if options.wildcard_file:
            formatters.add_detail_formatter(WildcardFormatter(options.wildcard_file))

        # Create writers
        writers = Writers(force_create, options.compression_format)

        # Create record handler
        record_handler = RecordHandler(modifiers, filters, formatters)

        if options.metrics:
            record_handler = MetricsRecordHandlerWrapper(
                record_handler,
                options.paired,
                options.metrics,
                qualities=delivers_qualities,
                quality_base=options.quality_base,
            )

        num_adapters = sum(len(a) for a in modifiers.get_adapters())
        trim_mode = {
            False: "single-end",
            "first": "paired-end legacy",
            "both": "paired-end",
        }[options.paired]
        logger.info(
            f"Trimming {num_adapters} adapter{'s' if num_adapters > 1 else ''} "
            f"with at most {options.error_rate * 100:.1f} errors in {trim_mode} mode..."
        )

        if options.paired == "first" and (
            len(record_handler.modifiers.get_modifiers(read=2)) > 0
            or options.quality_cutoff
        ):
            logger.warning(
                """Requested read modifications are applied only to the first read 
                since backwards compatibility mode is enabled. To modify both reads, 
                also use any of the -A/-B/-G/-U options. Use a dummy adapter sequence 
                when necessary: -A XXX"""
            )

        if options.paired:
            mixin_class = PairedEndPipelineMixin
        else:
            mixin_class = SingleEndPipelineMixin

        if options.threads is None:
            # Run single-threaded version
            result_handler = WorkerResultHandler(WriterResultHandler(writers))
            pipeline_class = type("TrimPipelineImpl", (mixin_class, TrimPipeline), {})
            pipeline = pipeline_class(record_handler, result_handler)
            self.summary.update(mode="serial", threads=1)
            return run_interruptible(pipeline, self, raise_on_error=True)
        else:
            # Run multiprocessing version
            self.summary.update(mode="parallel", threads=options.threads)
            return self.run_parallel(record_handler, writers, mixin_class)

    def run_parallel(
        self,
        record_handler: RecordHandler,
        writers: Writers,
        mixin_class: Type[Union[SingleEndPipelineMixin, PairedEndPipelineMixin]],
    ):
        """
        Parallel implementation of run_atropos. Works as follows:

        1. Main thread creates N worker processes (where N is the number of
        threads to be allocated) and (optionally) one writer process.
        2. Main thread loads batches of reads (or read pairs) from input file(s)
        and adds them to a queue (the input queue).
        3. Worker processes take batches from the input queue, process them as
        Atropos normally does, and either add the results to the result queue
        (if using a writer process) or write the results to disk. Each result is
        a dict mapping output file names to strings, where each string is a
        concatenation of reads (with appropriate line endings) to be written.
        A parameter also controls whether data compression is done by the workers or
        the writer.
        4. If using a writer process, it takes results from the result queue and
        writes each string to its corresponding file.
        5. When the main process finishes loading reads from the input file(s),
        it sends a signal to the worker processes that they should complete when
        the input queue is empty. It also singals the writer process how many
        total batches to expect, and the writer process exits after it has processed
        that many batches.
        6. When a worker process completes, it adds a summary of its activity to
        the summary queue.
        7. The main process reads summaries from the summary queue and merges
        them to create the complete summary, which is used to generate the report.

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
        # run_parallel method to avoid extra overhead if only running in serial mode.
        from multiprocessing import Queue
        from atropos.commands.multicore import ParallelPipelineMixin
        from atropos.commands.trim.multicore import (
            CompressingWorkerResultHandler,
            Done,
            Killed,
            OrderPreservingWriterResultHandler,
            ParallelTrimPipelineRunner,
            QueueResultHandler,
            ResultProcess,
            WriterManager,
        )
        from atropos.utils.multicore import DEFAULT_RETRY_INTERVAL

        # Main process
        timeout = max(self.get_option("process_timeout"), DEFAULT_RETRY_INTERVAL)
        threads = self.get_option("threads")
        logger.debug(
            f"Starting atropos in parallel mode with threads={threads}, "
            f"timeout={timeout}"
        )

        if threads < 2:
            raise ValueError("'threads' must be >= 2")

        compression_mode = self.get_option("compression_mode")

        if compression_mode is None:
            compression_mode = choose_compression_mode(self.options)

        # Reserve a thread for the writer process if it will be doing the compression
        # and if one is available.
        if compression_mode == "writer" and threads > 2:
            threads -= 1

        # Queue by which results are sent from the worker processes to the writer
        # process
        result_queue = Queue(self.get_option("result_queue_size"))
        writer_manager = None

        if self.get_option("writer_process"):
            if compression_mode == "writer":
                worker_result_handler = WorkerResultHandler(
                    QueueResultHandler(result_queue)
                )
            else:
                worker_result_handler = CompressingWorkerResultHandler(
                    QueueResultHandler(result_queue),
                    compression_format=self.get_option("compression_format")
                )

            writer_manager = WriterManager(
                writers,
                compression_mode,
                self.get_option("preserve_order"),
                result_queue,
                timeout,
            )
        else:
            worker_result_handler = WorkerResultHandler(
                WriterResultHandler(writers, use_suffix=True)
            )

        pipeline_class = type(
            "TrimPipelineImpl", (ParallelPipelineMixin, mixin_class, TrimPipeline), {}
        )
        pipeline = pipeline_class(record_handler, worker_result_handler)

        runner = ParallelTrimPipelineRunner(self, pipeline, threads, writer_manager)

        return runner.run()


def choose_compression_mode(options: Namespace) -> str:
    if (
        options.output is None or
        options.output in (STDOUT, STDERR) or (
            options.writer_process and
            2 < options.threads < 8 and
            options.can_use_system_compression
        )
    ):
        # We must use a writer process to write to stdout/stderr. Otherwise, our
        # tests show that with 8 or more threads, worker compression is more efficient.
        return "writer"
    else:
        return "worker"
