from atropos.commands.trim.pipeline import (
    RecordHandler, WorkerResultHandler, WriterResultHandler,
)
from atropos.commands import Command, PairedEndPipelineMixin, SingleEndPipelineMixin
from atropos.commands.trim.pipeline import TrimPipeline
from atropos.commands.trim.writers import Writers
from atropos.utils import ReturnCode, run_interruptible


class SingleEndSerialTrimPipeline(SingleEndPipelineMixin, TrimPipeline):
    pass


class PairedEndSerialTrimPipeline(PairedEndPipelineMixin, TrimPipeline):
    pass


def run_serial(
    command: Command,
    paired: bool,
    record_handler: RecordHandler,
    writers: Writers,
) -> ReturnCode:
    """

    Args:
        command:
        paired:
        record_handler:
        writers:

    Returns:
        ReturnCode
    """
    if paired:
        pipeline_class = PairedEndSerialTrimPipeline
    else:
        pipeline_class = SingleEndSerialTrimPipeline

    result_handler = WorkerResultHandler(WriterResultHandler(writers))

    return run_interruptible(
        pipeline_class(record_handler, result_handler),
        command,
        raise_on_error=True
    )
