import logging
from multiprocessing import Process
import os
from atropos.commands.trim import (
    ResultHandler, WorkerResultHandler, WriterResultHandler)
from atropos.commands.multicore import (
    Control, PendingQueue, ParallelPipelineRunner, MulticoreError, 
    wait_on_process, enqueue, dequeue, kill, CONTROL_ACTIVE, CONTROL_ERROR)
from atropos.io.compression import get_compressor

class Done(MulticoreError):
    """Raised when process exits normally.
    """
    pass

class Killed(MulticoreError):
    """Raised when process exits is killed.
    """
    pass

class ParallelTrimPipelineRunner(ParallelPipelineRunner):
    """ParallelPipelineRunner for a TrimPipeline.
    """
    def __init__(
            self, command_runner, pipeline, threads, writer_manager=None):
        super().__init__(command_runner, pipeline, threads)
        self.writer_manager = writer_manager
    
    def ensure_alive(self):
        super().ensure_alive()
        if self.writer_manager and not self.writer_manager.is_active():
            raise MulticoreError("Writer process exited")
    
    def after_enqueue(self):
        # Tell the writer thread the max number of batches to expect
        if self.writer_manager:
            self.writer_manager.set_num_batches(self.num_batches)
    
    def finish(self):
        if self.writer_manager:
            # Wait for writer to complete
            self.writer_manager.wait()
    
    def terminate(self, retcode):
        super().terminate(retcode)
        if self.writer_manager:
            self.writer_manager.terminate(retcode)

class QueueResultHandler(ResultHandler):
    """ResultHandler that writes results to the output queue.
    """
    def __init__(self, queue):
        self.queue = queue
        self.message = None
        self.timeout = None
    
    def start(self, worker):
        self.message = "{} waiting to queue result {{}}".format(
            worker.name)
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_compressors = None
    
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
        """Returns the file compressor based on the file extension.
        """
        if filename not in self.file_compressors:
            self.file_compressors[filename] = get_compressor(filename)
        return self.file_compressors[filename]

class OrderPreservingWriterResultHandler(WriterResultHandler):
    """Writer thread that is less time/memory efficient, but is
    guaranteed to preserve the original order of records.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pending = None
        self.cur_batch = None
    
    def start(self, worker=None):
        super().start(worker)
        self.pending = PendingQueue()
        self.cur_batch = 1
    
    def write_result(self, batch_num, result):
        if batch_num == self.cur_batch:
            self.writers.write_result(result, self.compressed)
            self.cur_batch += 1
            self.consume_pending()
        else:
            self.pending.push(batch_num, result)
    
    def finish(self, total_batches=None):
        if total_batches is not None:
            self.consume_pending()
            if self.cur_batch != total_batches + 1:
                raise MulticoreError(
                    "OrderPreservingWriterResultHandler finishing "
                    "without having seen {} of {} batches".format(
                        total_batches + 1 - self.cur_batch,
                        total_batches))
        super().finish(total_batches=total_batches)
    
    def consume_pending(self):
        """Consume any remaining items in the queue.
        """
        while (
                (not self.pending.empty) and
                (self.cur_batch == self.pending.min_priority)):
            self.writers.write_result(
                self.pending.pop(), self.compressed)
            self.cur_batch += 1

class ResultProcess(Process):
    """Thread that accepts results from the worker threads and process
    them using a ResultHandler. Each batch is expected to be
    (batch_num, path, records), where path is the destination file and
    records is a string. Not guaranteed to preserve the original order
    of sequence records.
    
    Args:
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
            "Writer process %s running under pid %d",
            self.name, os.getpid())
        
        def fail_callback():
            """Raises Done if the expected number of batches has been
            seen.
            """
            if (
                    self.num_batches is None and
                    self.control.check_value_positive()):
                self.num_batches = self.control.get_value()
            if (
                    self.num_batches is not None and
                    len(self.seen_batches) >= self.num_batches):
                raise Done()
        
        def timeout_callback():
            """Logs an error with the missing batches.
            """
            if self.num_batches is not None:
                missing = (
                    set(range(1, self.num_batches+1)) -
                    self.seen_batches)
                logging.getLogger().error(
                    "Result thread still missing batches %s of %d",
                    ",".join(str(i) for i in missing),
                    self.num_batches)
        
        def iter_batches():
            """Dequeue and yield batches.
            """
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
            self.result_handler.finish(
                num_batches if num_batches > 0 else None)

class WriterManager(object):
    """Manager for a writer process and control variable.
    """
    def __init__(
            self, writers, compression, preserve_order, result_queue,
            timeout):
        # result handler
        if preserve_order:
            writer_result_handler = OrderPreservingWriterResultHandler(
                writers, compressed=compression == "worker")
        else:
            writer_result_handler = WriterResultHandler(
                writers, compressed=compression == "worker")
        
        self.timeout = timeout
        # Shared variable for communicating with writer thread
        self.writer_control = Control(CONTROL_ACTIVE)
        # writer process
        self.writer_process = ResultProcess(
            writer_result_handler, result_queue, self.writer_control,
            timeout)
        self.writer_process.start()
    
    def is_active(self):
        """Returns True if the writer process is alive and the control
        value is CONTROL_ACTIVE.
        """
        return (
            self.writer_process.is_alive() and
            self.writer_control.check_value(CONTROL_ACTIVE))
    
    def set_num_batches(self, num_batches):
        """Set the number of batches to the control variable.
        """
        self.writer_control.set_value(num_batches)
    
    def wait(self):
        """Wait for the writer process to terminate.
        """
        wait_on_process(self.writer_process, self.timeout)
    
    def terminate(self, retcode):
        """Force the writer process to terminate.
        """
        kill(self.writer_process, retcode, self.timeout)