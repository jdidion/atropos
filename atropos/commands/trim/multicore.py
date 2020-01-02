from multiprocessing import Process, Queue
import os
from typing import Dict, Iterable, Optional

from loguru import logger
import xphyle
from xphyle.formats import CompressionFormat

from atropos.commands import Command, PairedEndPipelineMixin, SingleEndPipelineMixin
from atropos.commands.multicore import (
    MulticorePipelineMixin, ParallelPipelineRunner, WorkerProcess
)
from atropos.commands.trim.pipeline import (
    CompressionMode,
    RecordHandler,
    ResultHandler,
    WorkerResultHandler,
    WriterResultHandler,
    TrimPipeline,
)
from atropos.commands.trim.writers import Writers
from atropos.utils import ReturnCode
from atropos.utils.multicore import (
    DEFAULT_RETRY_INTERVAL,
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
        compression_mode: CompressionMode,
        preserve_order: bool,
        result_queue: Queue,
        timeout: int,
    ):
        compressed = (compression_mode == CompressionMode.WORKER)

        if preserve_order:
            writer_result_handler = OrderPreservingWriterResultHandler(
                writers, compressed=compressed
            )
        else:
            writer_result_handler = WriterResultHandler(
                writers, compressed=compressed
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
        command: Command,
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


class SingleEndMulticoreTrimPipeline(
    MulticorePipelineMixin, SingleEndPipelineMixin, TrimPipeline
):
    pass


class PairedEndMulticoreTrimPipeline(
    MulticorePipelineMixin, PairedEndPipelineMixin, TrimPipeline
):
    pass


def run_multicore(
    trim_command: Command,
    paired: bool,
    record_handler: RecordHandler,
    writers: Writers,
    threads: int,
    result_queue_size: int,
    writer_process: bool,
    preserve_order: bool,
    compression_mode: CompressionMode,
    compression_format: str,
    process_timeout: int,
) -> ReturnCode:
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
        threads:
        process_timeout:
        compression_mode:
        result_queue_size:
        writer_process:
        compression_format:
        writers: Writers object
        preserve_order:
        paired: Whether the input is paired-end.
        record_handler: RecordHandler object.
        trim_command:

    Returns:
        The return code.
    """
    timeout = max(process_timeout, DEFAULT_RETRY_INTERVAL)
    logger.debug(
        f"Starting atropos in multicore mode with threads={threads}, "
        f"timeout={timeout}"
    )

    if threads < 2:
        raise ValueError("'threads' must be >= 2")

    # Reserve a thread for the writer process if it will be doing the compression
    # and if one is available.
    if compression_mode == CompressionMode.WRITER and threads > 2:
        threads -= 1

    # Queue by which results are sent from the worker processes to the writer
    # process
    result_queue = Queue(result_queue_size)
    writer_manager = None

    if writer_process:
        if compression_mode == CompressionMode.WRITER:
            worker_result_handler = WorkerResultHandler(
                QueueResultHandler(result_queue)
            )
        else:
            worker_result_handler = CompressingWorkerResultHandler(
                QueueResultHandler(result_queue),
                compression_format=compression_format
            )

        writer_manager = WriterManager(
            writers,
            compression_mode,
            preserve_order,
            result_queue,
            timeout,
        )
    else:
        worker_result_handler = WorkerResultHandler(
            WriterResultHandler(writers, use_suffix=True)
        )

    if paired:
        pipeline_class = PairedEndMulticoreTrimPipeline
    else:
        pipeline_class = SingleEndMulticoreTrimPipeline

    pipeline = pipeline_class(record_handler, worker_result_handler)

    runner = ParallelTrimPipelineRunner(
        trim_command, pipeline, threads, writer_manager
    )

    return runner.run()
