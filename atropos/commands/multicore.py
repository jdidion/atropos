from multiprocessing import Process, Queue
import os
from typing import Iterator, Optional, Type, cast

from loguru import logger

from atropos.commands import Batch, Command, Pipeline
from atropos.errors import AtroposError, MulticoreError
from atropos.utils import ReturnCode, run_interruptible
from atropos.utils.multicore import (
    DEFAULT_RETRY_INTERVAL,
    dequeue,
    enqueue,
    enqueue_all,
    ensure_processes,
    kill,
    wait_on,
)


class ParallelPipelineMixin:
    """
    Mixin that implements the `start`, `finish`, and `process_batch` methods of
    :class:`Pipeline`.
    """

    seen_batches = None

    def start(self, **kwargs):
        cast(Pipeline, super()).start(**kwargs)
        self.seen_batches = set()

    def process_batch(self, batch: Batch):
        self.seen_batches.add(batch[0]["index"])
        cast(Pipeline, super()).process_batch(batch)

    def finish(self, summary: dict, worker: "WorkerProcess" = None):
        cast(Pipeline, super()).finish(summary, worker=worker)
        logger.debug(
            f"{worker.name} finished; processed {len(self.seen_batches)} batches, "
            f"{sum(cast(Pipeline, self).record_counts.values())} reads"
        )


class WorkerProcess(Process):
    """
    Parent class for worker processes that execute Pipelines.
    """

    def __init__(
        self,
        index: int,
        input_queue: Queue,
        pipeline: Pipeline,
        summary_queue: Queue,
        timeout: int,
    ):
        """
        Args:
            index: A unique ID for the process.
            input_queue: Queue with batches of records to process.
            pipeline: The pipeline to execute.
            summary_queue: Queue where summary information is written.
            timeout: Time to wait upon queue full/empty.
        """
        super().__init__(name=f"Worker process {index}")
        self.index = index
        self.input_queue = input_queue
        self.pipeline = pipeline
        self.summary_queue = summary_queue
        self.timeout = timeout

    def run(self) -> None:
        logger.debug(f"{self.name} running under pid {os.getpid()}")
        summary = {}

        def iter_batches() -> Iterator:
            """
            Deques and yields batches.
            """
            while True:
                yield dequeue(
                    self.input_queue,
                    wait_message=f"{self.name} waiting on batch {{}}",
                    block_timeout=self.timeout,
                )

        def enqueue_summary():
            """
            Enqueues a summary dict.
            """
            enqueue(
                self.summary_queue,
                (
                    self.index,
                    cast(ParallelPipelineMixin, self.pipeline).seen_batches,
                    summary,
                ),
                wait_message=f"{self.name} waiting to queue summary {{}}",
                block_timeout=self.timeout,
            )

        try:
            self.pipeline.start(worker=self)

            try:
                for batch in iter_batches():
                    if batch is None:
                        break

                    logger.debug(
                        f"{self.name} processing batch {batch[0]['index']} of size "
                        f"{batch[0]['size']}"
                    )
                    self.pipeline.process_batch(batch)
            finally:
                self.pipeline.finish(summary, worker=self)

            logger.debug(f"{self.name} finished normally")
        except Exception as err:
            logger.exception(f"Unexpected error in {self.name}")
            summary["exception"] = err

        logger.debug(f"{self.name} sending summary")
        enqueue_summary()


class ParallelPipelineRunner:
    """
    Runs a pipeline in parallel.
    """

    def __init__(
        self, command: Command, pipeline: Pipeline, threads: Optional[int] = None
    ):
        """
        Args:
            command:
            pipeline: A :class:`Pipeline`.
            threads: Number of threads to use. If None, the value will be taken
                from command_runner.
        """
        self.command = command
        self.pipeline = pipeline
        self.threads = threads or command.get_option("threads")
        self.timeout = max(
            command.get_option("process_timeout", 0), DEFAULT_RETRY_INTERVAL
        )
        # Queue by which batches of reads are sent to worker processes
        self.input_queue = Queue(command.get_option("read_queue_size"))
        # Queue for processes to send summary information back to main process
        self.summary_queue = Queue(self.threads)
        self.worker_processes = None
        self.num_batches = None
        self.seen_summaries = None
        self.seen_batches = None

    def ensure_alive(self) -> None:
        """
        Called when enqueue times out.
        """
        ensure_processes(self.worker_processes)

    def after_enqueue(self) -> None:
        """
        Called after all batches are queued.
        """

    def finish(self) -> None:
        """
        Called at the end of the run.
        """

    def run(self) -> ReturnCode:
        """
        Runs the pipeline.

        Returns:
            The return code.
        """
        retcode = run_interruptible(self)
        self.terminate(retcode)
        return retcode

    def terminate(self, retcode: ReturnCode) -> None:
        """
        Terminates all processes.

        Args:
            retcode: The return code.
        """
        # notify all threads that they should stop
        logger.debug("Exiting all processes")
        for process in self.worker_processes:
            kill(process, retcode, self.timeout)

    def __call__(self) -> None:
        # Start worker processes, reserve a thread for the reader process,
        # which we will get back after it completes
        worker_args = (
            self.input_queue,
            self.pipeline,
            self.summary_queue,
            self.timeout,
        )
        self.worker_processes = launch_workers(self.threads - 1, worker_args)
        self.num_batches = enqueue_all(
            self.command.iterator(), self.input_queue, self.timeout, self.ensure_alive,
        )

        logger.debug(f"Main loop complete; saw {self.num_batches} batches")
        # Tell the worker processes no more input is coming
        enqueue_all(
            (None,) * self.threads, self.input_queue, self.timeout, self.ensure_alive
        )
        self.after_enqueue()
        # Now that the reader process is done, it essentially frees up another thread
        # to use for a worker
        self.worker_processes.extend(
            launch_workers(1, worker_args, offset=self.threads - 1)
        )

        # Wait for all summaries to be available on queue
        def summary_timeout_callback():
            """
            Ensures that workers are still alive.
            """
            try:
                ensure_processes(
                    self.worker_processes,
                    "Workers are still alive and haven't returned summaries: {}",
                    alive=False,
                )
            except:
                logger.exception("Error calling ensure_processes")

        wait_on(
            self.summary_queue.full,
            wait_message="Waiting on worker summaries {}",
            timeout=self.timeout,
            wait=True,
            timeout_callback=summary_timeout_callback,
        )

        logger.debug("Processing summary information from worker processes")
        self.seen_summaries = set()
        self.seen_batches = set()

        def summary_fail_callback():
            """
            Raises AtroposError with workers that did not report summaries.
            """
            missing_summaries = set(range(1, self.threads)) - self.seen_summaries
            raise AtroposError(
                f"Missing summaries from processes "
                f"{','.join(str(summ) for summ in missing_summaries)}"
            )

        for _ in range(1, self.threads + 1):
            batch = dequeue(self.summary_queue, fail_callback=summary_fail_callback)
            worker_index, worker_batches, worker_summary = batch
            if worker_summary is None:
                raise MulticoreError(f"Worker process {worker_index} died unexpectedly")
            elif (
                "exception" in worker_summary
                and worker_summary["exception"] is not None
            ):
                raise AtroposError(
                    f"Worker process {worker_index} died unexpectedly",
                    worker_summary["exception"],
                )
            else:
                logger.debug(f"Processing summary for worker {worker_index}")

            self.seen_summaries.add(worker_index)
            self.seen_batches |= worker_batches
            self.command.merge_summary(worker_summary)

        # Check if any batches were missed
        if self.num_batches > 0:
            missing_batches = set(range(1, self.num_batches + 1)) - self.seen_batches
            if len(missing_batches) > 0:
                raise AtroposError(
                    f"Workers did not process batches "
                    f"{','.join(str(batch) for batch in missing_batches)}"
                )

        self.finish()


def launch_workers(
    num_workers: int,
    args=(),
    offset: int = 0,
    worker_class: Type[WorkerProcess] = WorkerProcess
):
    """
    Launches `num_workers` workers. Each worker is initialized with an incremental
    index starting with `offset`, followed by `args`.
    """
    logger.info(f"Starting {num_workers} worker processes")
    # create workers
    workers = [worker_class(i + offset, *args) for i in range(num_workers)]
    # start workers
    for worker in workers:
        worker.start()
    return workers
