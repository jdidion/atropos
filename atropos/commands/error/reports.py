from itertools import repeat
from typing import IO

from atropos.commands.reports import BaseReportGenerator, BaseReportWriter, ReportWriter
from atropos.commands.legacy_reports import Printer, TitlePrinter
from atropos.utils.argparse import Namespace


class ErrorTextReportWriter(BaseReportWriter):
    def __init__(self, options: Namespace):
        super().__init__("txt", options)

    def serialize(self, summary: dict, stream: IO) -> None:
        names = summary['input']['input_names'] or repeat(None)
        estimates = summary['errorrate']['estimate']
        _print = Printer(stream)
        _print_title = TitlePrinter(stream)

        input_idx = 0

        for input_idx, (estimate, details, name) in enumerate(
            zip(estimates, summary['errorrate']['details'], names), 1
        ):
            self._generate_estimator_report(
                stream, input_idx, estimate, details, _print, _print_title, name
            )

        if input_idx > 1:
            _print.newline()
            _print_title("Overall", level=0)
            total_lens = summary['errorrate']['total_len']
            overall_err = (
                sum(err * total_len for err, total_len in zip(estimates, total_lens)) /
                sum(total_lens)
            )
            print(f"Error rate: {overall_err:.2%}", file=stream)

    @staticmethod
    def _generate_estimator_report(
        stream: IO,
        input_idx: int,
        estimate: float,
        details: dict,
        _print: Printer,
        _print_title: Printer,
        input_name: str = None
    ):
        _print_indent = Printer(stream, indent='  ')
        _print.newline()
        _print_title(f"Input {input_idx}", level=0)

        if input_name:
            _print(f"File: {input_name}")

        _print(f"Error rate: {estimate:.2%}")

        if details:
            _print("Details:\n")
            per_read = details['per_read']
            per_cycle = details['per_cycle']
            _print_indent(f"StdErr: {per_read['standard error']:.2%}")
            _print_indent("Per-cycle rates:")

            for cycle in per_cycle:
                _print_indent(
                    "Cycle: {}, Error: {:.2%}, StdErr: {:.2%}".format(*cycle), indent=2
                )


class ErrorReportGenerator(BaseReportGenerator):
    @classmethod
    def _create_report_writer(cls, fmt: str, options: Namespace) -> ReportWriter:
        if fmt == "txt":
            return ErrorTextReportWriter(options)
        else:
            return super()._create_report_writer(fmt, options)
