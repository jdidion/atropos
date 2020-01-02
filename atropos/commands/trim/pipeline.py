from abc import ABCMeta, abstractmethod
from enum import Enum
from multiprocessing import Process
from typing import (
    Iterable,
    Optional,
    Sequence as SequenceType,
    Tuple,
    Type,
)

from atropos.commands import Pipeline, Summary
from atropos.commands.multicore import WorkerProcess
from atropos.commands.trim.modifiers import Modifiers
from atropos.commands.trim.filters import Filter, Filters
from atropos.commands.trim.writers import Formatters, Writers
from atropos.io.sequence import Sequence
from atropos.utils.collections import Summarizable


class CompressionMode(Enum):
    WORKER = "worker"
    WRITER = "writer"


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
