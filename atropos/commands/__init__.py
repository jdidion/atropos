from abc import ABCMeta, abstractmethod
from collections.abc import Sequence as SequenceCollection
import copy
from functools import lru_cache
import platform
import sys
from typing import Any, Dict, Iterator, Optional, Sequence as SequenceType, Tuple, cast

from atropos import __version__
from atropos.adapters import AdapterCache
from atropos.errors import AtroposError
from atropos.io.readers import Sequence, open_reader
from atropos.utils import ReturnCode, classproperty, create_progress_reader
from atropos.utils.argparse import Namespace
from atropos.utils.collections import Const, MergingDict, Summarizable, Timing


EXCEPTION_KEY = "exception"

Batch = Tuple[dict, SequenceType[Sequence]]


class Summary(MergingDict):
    """
    Container for summary information.
    """

    @property
    def has_exception(self) -> bool:
        return EXCEPTION_KEY in self

    def finish(self) -> None:
        self._post_process_dict(self)

    def _post_process_dict(self, dict_val: Optional[Dict[str, Any]]) -> None:
        """
        Replaces `Summarizable` members with their summaries, and computes some
        aggregate values.

        Args:
            dict_val:
        """
        if dict_val is None:
            return

        for key, value in tuple(dict_val.items()):
            if value is None:
                continue

            if isinstance(value, Summarizable):
                dict_val[key] = value = value.summarize()

            if isinstance(value, dict):
                self._post_process_dict(value)
            elif (
                isinstance(value, SequenceCollection)
                and len(value) > 0
                and all(val is None or isinstance(val, dict) for val in value)
            ):
                for val in value:
                    self._post_process_dict(val)
            else:
                if isinstance(value, Const):
                    dict_val[key] = value = value.value
                self._post_process_other(dict_val, key, value)

    def _post_process_other(self, parent: dict, key: str, value: Any) -> None:
        """
        Override this method to perform additional post-processing on non-dict values.

        Args:
            parent: The parent dict
            key: The key
            value: The current value
        """


class Command(Iterator[Batch], metaclass=ABCMeta):
    """
    An Atropos command.
    """

    @abstractmethod
    def get_option(self, name: str, default: Optional[Any] = None) -> Any:
        """
        Gets the value of the specified option.

        Args:
            name: The option name.
            default: The default value to return if option `name` is not set.

        Returns:
            The option value.
        """

    @abstractmethod
    def iterator(self) -> Iterator[Batch]:
        """
        Returns an iterator over input batches.
        """

    @abstractmethod
    def merge_summary(self, summary: dict) -> None:
        """
        Merge a summary dict into the command's current summary.

        Args:
            summary: The summary to merge
        """


class BaseCommand(Command, metaclass=ABCMeta):
    """
    Base class for Atropos commands.
    """

    @classproperty
    @abstractmethod
    def name(cls) -> str:
        """
        The command name.
        """

    @classmethod
    def _create_summary(cls) -> dict:
        return Summary()

    def __init__(self, options: Namespace):
        """
        Args:
            options: Command-line options.
        """
        self.options = options
        self.summary: dict = self._create_summary()
        self.timing = Timing()
        self.return_code: Optional[ReturnCode] = None
        self.size = options.batch_size or 1000
        self.batches = 0
        self.done = False
        self._empty_batch = [None] * self.size
        self._progress_options = None

        open_args = dict(
            file_format=options.input_format,
            quality_base=options.quality_base,
            colorspace=options.colorspace,
            input_read=options.input_read,
            alphabet=options.alphabet,
        )

        if options.ngstream_reader:
            open_args["ngstream_reader"] = options.ngstream_reader
            options.ngstream_reader = None
        else:
            interleaved = bool(options.interleaved_input)
            input1 = options.interleaved_input if interleaved else options.input1
            input2 = qualfile = None

            if options.paired and not interleaved:
                input2 = options.input2
            else:
                qualfile = options.input2

            open_args.update(
                file1=input1,
                file2=input2,
                qualfile=qualfile,
                interleaved=interleaved,
            )

        self.reader = reader = open_reader(**open_args)

        # Wrap reader in subsampler
        if options.subsample:
            import random

            if options.subsample_seed:
                random.seed(options.subsample_seed)

            def subsample(_reader, frac):
                """Generator that yields a random subsample of records.

                Args:
                    _reader: The reader from which to sample.
                    frac: The fraction of records to yield.
                """
                for reads in _reader:
                    if random.random() < frac:
                        yield reads

            reader = subsample(reader, options.subsample)

        self.iterable = enumerate(reader, 1)

        if options.progress:
            max_reads = self.get_option("max_reads")

            if not max_reads:
                max_reads = self.reader.estimate_num_records()
                if max_reads and options.subsample:
                    max_reads = int(max_reads * options.subsample)

            self._progress_options = (
                options.progress,
                self.size,
                max_reads,
                options.counter_magnitude,
            )

        self._init_summary()

    @lru_cache(maxsize=None)
    def get_option(self, name: str, default: Optional[Any] = None) -> Any:
        if hasattr(self.reader, name):
            return getattr(self.reader, name)
        elif hasattr(self.options, name):
            return getattr(self.options, name)
        else:
            return default

    def iterator(self) -> Iterator[Batch]:
        """
        Returns an iterator (an object with the __iter__ method) over input batches.
        BaseCommand is itself an iterable, and will be wrapped with a progress bar if
        `_progress_options` is True.
        """
        if not self._progress_options:
            return self

        # Wrap iterator in progress bar
        return create_progress_reader(self, *self._progress_options)

    def merge_summary(self, summary: dict) -> None:
        if isinstance(self.summary, MergingDict):
            cast(MergingDict, self.summary).merge(summary)
        else:
            self.summary.update(summary)

    def __iter__(self):
        return self

    def __next__(self) -> Batch:
        if self.done:
            raise StopIteration()

        try:
            read_index, record = next(self.iterable)
        except StopIteration:
            self.finish()
            raise

        batch = copy.copy(self._empty_batch)
        batch[0] = record
        batch_index = 1
        max_reads = self.get_option("max_reads")
        if max_reads:
            max_size = min(self.size, max_reads - read_index + 1)
        else:
            max_size = self.size

        while batch_index < max_size:
            try:
                read_index, record = next(self.iterable)
                batch[batch_index] = record
                batch_index += 1
            except StopIteration:
                self.finish()
                break
            except:
                self.finish()
                raise

        if max_reads and read_index >= max_reads:
            self.finish()

        self.batches += 1

        batch_meta = dict(
            index=self.batches,
            # TODO: When multi-file input is supported, 'source' will need to
            # be the index of the current file/pair from which records are
            # being read.
            source=0,
            size=batch_index,
        )

        if batch_index == self.size:
            return batch_meta, batch
        else:
            return batch_meta, batch[0:batch_index]

    def _init_summary(self) -> None:
        """
        Initializes the summary dict with general information.
        """
        self.summary["program"] = "Atropos"
        self.summary["version"] = __version__
        self.summary["python"] = platform.python_version()
        self.summary["command"] = self.name
        self.summary["options"] = self.options.__dict__.copy()
        self.summary["timing"] = self.timing
        self.summary["sample_id"] = self.options.sample_id
        self.summary["input"] = self.reader.summarize()
        self.summary["input"].update(
            batch_size=self.size,
            max_reads=self.get_option("max_reads"),
            batches=self.batches
        )

    def run(self) -> Tuple[ReturnCode, dict]:
        """
        Runs the command, wrapping it in a Timing, catching any exceptions,
        and finally closing the reader.

        Returns:
            The tuple (retcode, summary).
        """
        with self.timing:
            try:
                self.return_code = self()
            except Exception as err:  # pylint: disable=broad-except
                self.summary["exception"] = dict(
                    message=str(err), details=sys.exc_info()
                )
                self.return_code = 1
            finally:
                self.finish()
        return self.return_code, self.summary

    @abstractmethod
    def __call__(self) -> ReturnCode:
        """
        Executes the command.

        Returns:
            The return code.
        """

    def finish(self) -> None:
        """
        Finishes the command.
        """
        # Close the underlying reader.
        if not self.done:
            self.done = True
            self.reader.close()
        if isinstance(self.summary, Summary):
            cast(Summary, self.summary).finish()

    def _load_known_adapters(self) -> AdapterCache:
        """
        Load known adapters based on setting in command-line options.
        """
        cache_file = None
        if self.options.cache_adapters:
            cache_file = self.options.adapter_cache_file
        adapter_cache = AdapterCache(cache_file)
        if adapter_cache.empty and self.options.default_adapters:
            adapter_cache.load_default()
        if self.options.known_adapter:
            for known in self.options.known_adapter:
                name, seq = known.split("=")
                adapter_cache.add(name, seq)
        if self.options.known_adapters_file:
            for known_file in self.options.known_adapters_file:
                adapter_cache.load_from_url(known_file)
        if self.options.cache_adapters:
            adapter_cache.save()
        return adapter_cache


class Pipeline(metaclass=ABCMeta):
    """
    Base class for analysis pipelines.
    """

    def __init__(self):
        self.record_counts = {}
        self.bp_counts = {}

    def __call__(
        self, command_runner: BaseCommand, raise_on_error: bool = False, **kwargs
    ) -> None:
        """
        Executes the pipeline.

        Args:
            command_runner:
            raise_on_error:
            **kwargs:

        Raises:

        """
        self.start(**kwargs)
        try:
            for batch in command_runner.iterator():
                self.process_batch(batch)
        except Exception as err:
            if raise_on_error:
                raise
            else:
                command_runner.summary["exception"] = dict(
                    message=str(err), details=sys.exc_info()
                )
        finally:
            self.finish(command_runner.summary, **kwargs)

    def start(self, **kwargs) -> None:
        """
        Starts the pipeline.
        """

    def process_batch(self, batch: Batch) -> None:
        """
        Runs the pipeline on a batch of records.

        Args:
            batch: A batch of reads. A batch has the format
            ({batch_metadata}, [records]).
        """
        batch_meta, records = batch
        context = batch_meta.copy()
        if context["source"] not in self.record_counts:
            self.record_counts[context["source"]] = 0
        self.record_counts[context["source"]] += context["size"]
        if not context["source"] in self.bp_counts:
            self.bp_counts[context["source"]] = [0, 0]
        context["bp"] = self.bp_counts[context["source"]]
        self.add_to_context(context)
        self.handle_records(context, records)

    def add_to_context(self, context: dict) -> None:
        """
        Adds items to the batch context.
        """

    def handle_records(self, context: dict, records: SequenceType[Sequence]) -> None:
        """
        Handles a sequence of records.

        Args:
            context: The pipeline context (dict).
            records: The sequence of records.

        Raises:
            AtroposError
        """
        for idx, record in enumerate(records):
            try:
                self.handle_record(context, record)
            except Exception as err:
                raise AtroposError(
                    f"An error occurred at record {idx} of batch {context['index']}"
                ) from err

    @abstractmethod
    def handle_record(self, context: dict, record: Sequence) -> None:
        """
        Handles a single record.

        Args:
            context: The pipeline context (dict).
            record: The record.
        """

    @abstractmethod
    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ) -> None:
        """
        Handles a read or read-pair.

        Args:
            context: The pipeline context (dict).
            read1: The first read in the pair.
            read2: The second read in the pair; will be `None` for single-end data.
        """

    def finish(self, summary: dict, **kwargs) -> None:
        """
        Finishes the pipeline, including adding information to the summary.

        Args:
            summary: Summary dict to update.
        """
        total_bp_counts = tuple(sum(b) for b in zip(*self.bp_counts.values()))
        summary.update(
            record_counts=self.record_counts,
            total_record_count=sum(self.record_counts.values()),
            bp_counts=self.bp_counts,
            total_bp_counts=total_bp_counts,
            sum_total_bp_count=sum(total_bp_counts),
        )


class SingleEndPipelineMixin:
    """
    Mixin for pipelines that implements `handle_record` for single-end data.
    """

    def handle_record(self, context: dict, record: Sequence) -> None:
        context["bp"][0] += len(record)
        cast(Pipeline, self).handle_reads(context, record)


class PairedEndPipelineMixin:
    """
    Mixin for pipelines that implements `handle_record` for paired-end data.
    """

    def handle_record(self, context: dict, record: Sequence) -> None:
        read1, read2 = record
        bps = context["bp"]
        bps[0] += len(read1.sequence)
        bps[1] += len(read2.sequence)
        cast(Pipeline, self).handle_reads(context, read1, read2)
