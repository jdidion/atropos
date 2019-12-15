from multiprocessing import Process, Queue
import os
from typing import Dict, Iterable, Optional

from loguru import logger
import xphyle
from xphyle.formats import CompressionFormat

from atropos.commands.multicore import ParallelPipelineRunner, WorkerProcess
from atropos.commands.trim import (
    ResultHandler,
    WorkerResultHandler,
    WriterResultHandler,
    TrimCommand,
    TrimPipeline,
)
from atropos.commands.trim.writers import Writers
from atropos.utils import ReturnCode
from atropos.utils.multicore import (
    Control,
    ControlSignal,
    MulticoreError,
    PendingQueue,
    enqueue,
    dequeue,
    kill,
    wait_on_process,
)


class Done(MulticoreError):
    """
    Raised when process exits normally.
    """

    pass


class Killed(MulticoreError):
    """
    Raised when process exits is killed.
    """

    pass


class WriterManager:
    """
    Manager for a writer process and control variable.
    """

    def __init__(
        self,
        writers: Writers,
        compression_mode: str,
        preserve_order: bool,
        result_queue: Queue,
        timeout: int,
    ):
        # result handler
        if preserve_order:
            writer_result_handler = OrderPreservingWriterResultHandler(
                writers, compressed=compression_mode == "worker"
            )
        else:
            writer_result_handler = WriterResultHandler(
                writers, compressed=compression_mode == "worker"
            )

        self.timeout = timeout
        # Shared variable for communicating with writer thread
        self.writer_control = Control(ControlSignal.ACTIVE)
        # writer process
        self.writer_process = ResultProcess(
            writer_result_handler, result_queue, self.writer_control, timeout
        )
        self.writer_process.start()

    @property
    def is_active(self) -> bool:
        """
        Returns True if the writer process is alive and the control value is
        ControlSignal.ACTIVE.
        """
        return self.writer_process.is_alive() and self.writer_control.check_value(
            ControlSignal.ACTIVE
        )

    def set_num_batches(self, num_batches: int):
        """
        Sets the number of batches to the control variable.
        """
        self.writer_control.set_value(num_batches)

    def wait(self):
        """
        Waits for the writer process to terminate.
        """
        wait_on_process(self.writer_process, self.timeout)

    def terminate(self, retcode: ReturnCode):
        """
        Forces the writer process to terminate.
        """
        kill(self.writer_process, retcode, self.timeout)


class ParallelTrimPipelineRunner(ParallelPipelineRunner):
    """
    ParallelPipelineRunner for a TrimPipeline.
    """

    def __init__(
        self,
        command: TrimCommand,
        pipeline: TrimPipeline,
        threads: int,
        writer_manager: Optional[WriterManager] = None,
    ):
        super().__init__(command, pipeline, threads)
        self.writer_manager = writer_manager

    def ensure_alive(self):
        super().ensure_alive()
        if self.writer_manager and not self.writer_manager.is_active:
            raise MulticoreError("Writer process exited")

    def after_enqueue(self):
        # Tell the writer thread the max number of batches to expect
        if self.writer_manager:
            self.writer_manager.set_num_batches(self.num_batches)

    def finish(self):
        if self.writer_manager:
            # Wait for writer to complete
            self.writer_manager.wait()

    def terminate(self, retcode: ReturnCode):
        super().terminate(retcode)
        if self.writer_manager:
            self.writer_manager.terminate(retcode)


class QueueResultHandler(ResultHandler):
    """
    ResultHandler that writes results to the output queue.
    """

    def __init__(self, queue: Queue):
        self.queue = queue
        self.message = None
        self.timeout = None

    def start(self, worker: Optional[WorkerProcess] = None):
        if worker is None:
            raise RuntimeError("worker must not be None")

        self.message = f"{worker.name} waiting to queue result {{}}"
        self.timeout = worker.timeout

    def write_result(self, batch_num: int, result: dict):
        enqueue(
            self.queue,
            (batch_num, result),
            wait_message=self.message,
            timeout=self.timeout,
        )


class CompressingWorkerResultHandler(WorkerResultHandler):
    """
    Wraps a ResultHandler and compresses results prior to writing.
    """

    def __init__(self, *args, compression_format: Optional[str] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.compression_format = compression_format
        self.file_compressors: Optional[Dict[str, CompressionFormat]] = None

    def start(self, worker: Optional[WorkerProcess] = None):
        super().start(worker)
        self.file_compressors = {}

    def prepare_file(self, path: str, strings: Iterable[str]):
        compressor = self.get_compressor(path)

        if compressor:
            return (
                (path, "wb"),
                compressor.compress(b"".join(s.encode() for s in strings)),
            )
        else:
            return (path, "wt"), "".join(strings)

    def get_compressor(self, filename: str) -> CompressionFormat:
        """
        Returns the file compressor based on the file extension.
        """
        if filename not in self.file_compressors:
            self.file_compressors[filename] = xphyle.get_compressor(
                self.compression_format or filename
            )

        return self.file_compressors[filename]


class OrderPreservingWriterResultHandler(WriterResultHandler):
    """
    Writer thread that is less time/memory efficient, but is guaranteed to preserve
    the original order of records.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pending = None
        self.cur_batch = None

    def start(self, worker: Optional[WorkerProcess] = None):
        super().start(worker)
        self.pending = PendingQueue()
        self.cur_batch = 1

    def write_result(self, batch_num: int, result: dict):
        if batch_num == self.cur_batch:
            self.writers.write_result(result, self.compressed)
            self.cur_batch += 1
            self.consume_pending()
        else:
            self.pending.push(batch_num, result)

    def finish(self, total_batches: Optional[int] = None):
        if total_batches is not None:
            self.consume_pending()

            if self.cur_batch != total_batches + 1:
                raise MulticoreError(
                    f"OrderPreservingWriterResultHandler finishing "
                    f"without having seen {total_batches + 1 - self.cur_batch} of "
                    f"{total_batches} batches"
                )

        super().finish(total_batches=total_batches)

    def consume_pending(self):
        """
        Consumes any remaining items in the queue.
        """
        while (not self.pending.empty) and (
            self.cur_batch == self.pending.min_priority
        ):
            self.writers.write_result(self.pending.pop(), self.compressed)
            self.cur_batch += 1


class ResultProcess(Process):
    """
    Thread that accepts results from the worker threads and process them using a
    ResultHandler. Each batch is expected to be (batch_num, path, records),
    where path is the destination file and records is a string. Not guaranteed to
    preserve the original order of sequence records.
    """

    def __init__(
        self,
        result_handler: ResultHandler,
        queue: Queue,
        control: Control,
        timeout: int = 60,
    ):
        """
        Args:
            result_handler: A ResultHandler object.
            queue: Input queue.
            control: A shared value for communcation with the main process.
            timeout: Seconds to wait for next batch before complaining.
        """
        super().__init__(name="Result process")
        self.result_handler = result_handler
        self.queue = queue
        self.control = control
        self.timeout = timeout
        self.seen_batches = set()
        self.num_batches = None

    def run(self):
        logger.debug(f"Writer process {self.name} running under pid {os.getpid()}")

        def fail_callback():
            """
            Raises Done if the expected number of batches has been seen.
            """
            if self.num_batches is None and self.control.check_value_positive():
                self.num_batches = self.control.get_value()

            if (
                self.num_batches is not None
                and len(self.seen_batches) >= self.num_batches
            ):
                raise Done()

        def timeout_callback():
            """
            Logs an error with the missing batches.
            """
            if self.num_batches is not None:
                missing = set(range(1, self.num_batches + 1)) - self.seen_batches
                logger.error(
                    f"Result thread still missing batches"
                    f"{','.join(str(i) for i in missing)} of {self.num_batches}"
                )

        def iter_batches():
            """
            Dequeues and yields batches.
            """
            while True:
                yield dequeue(
                    self.queue,
                    wait_message="Result process waiting on result {}",
                    timeout=self.timeout,
                    fail_callback=fail_callback,
                    timeout_callback=timeout_callback,
                )

        try:
            self.result_handler.start(self)
            for batch_num, result in iter_batches():
                self.seen_batches.add(batch_num)
                self.result_handler.write_result(batch_num, result)
        except Done:
            logger.debug("Writer process exiting normally")
        except Killed:
            logger.debug("Writer process exited early")
        except:
            logger.exception("Unexpected error in writer process")
            self.control.set_value(ControlSignal.ERROR)
        finally:
            num_batches = self.control.get_value(lock=True)
            self.result_handler.finish(num_batches if num_batches > 0 else None)
