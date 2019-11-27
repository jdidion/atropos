"""Interface to report generation.
"""
from abc import ABCMeta, abstractmethod
from importlib import import_module
from pathlib import Path
from typing import Dict, IO, Optional, Sequence, Type

from xphyle import STDOUT, STDERR, open_

from atropos.io.seqio import InputRead


class ReportWriter(metaclass=ABCMeta):
    def __init__(self, output_file: Path, options=None, **kwargs):
        self.output_file = output_file

    def __call__(self, summary: dict):
        with open_(self.output_file, "wt") as stream:
            self.serialize(summary, stream)

    @abstractmethod
    def serialize(self, summary: dict, stream: IO):
        pass


class ReportWriterFactory:
    def __init__(
        self,
        name: str,
        cls: Type[ReportWriter],
        ext: Optional[str] = None,
        default: bool = False,
        **kwargs
    ):
        self.name = name
        self.cls = cls
        self.ext = ext if ext is not None else name
        self.default = default
        self.kwargs = kwargs

    def __call__(self, output_file: Path, options):
        return self.cls(output_file, name=self.name, options=options, **self.kwargs)


class TemplateReportWriter(ReportWriter):
    """
    Args:
        output_file: Where to write the report.
        name: Report format. Must match the file extension on a discoverable
            template.
        template_name: A template name to use, rather than auto-discover.
        template_paths: Sequence of paths to search for templates.
        template_globals: Dict of additional globals to add to the template
            environment.
    """
    def __init__(
        self,
        output_file: Path,
        name: str = "txt",
        template_name: Optional[str] = None,
        template_paths: Sequence[Path] = None,
        template_globals: Optional[dict] = None,
        **_
    ):
        super().__init__(output_file)
        self.name = name
        self.template_name = template_name
        self.template_paths = template_paths
        self.template_globals = template_globals

    def serialize(self, summary: dict, stream: IO):
        """
        Generate a report using Jinja2.
        """
        import jinja2

        template_name = self.template_name or f"template.{self.name}"
        template_paths = self.template_paths or []
        if hasattr(self, "template_path"):
            template_paths.append(self.template_path)

        try:
            env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_paths))
            if self.template_globals:
                env.globals.update(self.template_globals)
            template = env.get_template(template_name)
        except:
            raise IOError(f"Could not load template file '{template_name}'")

        stream.write(template.render(summary=summary))


class DumpReportWriter(ReportWriter):
    def __init__(self, output_file: Path, name: str, options, **kwargs):
        super().__init__(output_file)
        self.mod_name = name
        self.kwargs = kwargs

    def serialize(self, summary: dict, stream: IO):
        mod = import_module(self.mod_name)
        mod.dump(summary, stream, **self.kwargs)


class BaseReportGenerator:
    """
    Base class for command report generators.

    Args:
        options: Command-line options.
    """

    def __init__(self, options, available_formats: Dict[str, ReportWriterFactory]):
        report_file = options.report_file
        report_formats = options.report_formats
        default_formats = [
            fmt for fmt, factory in available_formats.items() if factory.default
        ]

        if report_file in (STDOUT, STDERR):
            if not report_formats:
                report_formats = default_formats
            self.report_formats = [
                available_formats[fmt](report_file, options)
                for fmt in report_formats
            ]
        else:
            if not report_formats:
                ext = report_file.suffix
                if ext:
                    report_formats = (ext[1:],)
                else:
                    report_formats = default_formats

            if len(report_formats) == 1:
                self.report_formats = [
                    available_formats[report_formats[0]](report_file, options)
                ]
            else:
                self.report_formats = [
                    factory(Path(f"{report_file}.{factory.ext}"), options)
                    for factory in (available_formats[fmt] for fmt in report_formats)
                ]

    def generate_reports(self, summary: dict):
        """
        Generate report(s) from a summary.

        Args:
            summary: The summary dict.
        """
        self.add_derived_data(summary)
        for fmt in self.report_formats:
            fmt(summary)

    @staticmethod
    def add_derived_data(summary: dict) -> None:
        """
        Get some additional fields that are useful in the report.
        """
        derived = {
            "mean_sequence_lengths": tuple(
                None if bp is None else bp / summary["total_record_count"]
                for bp in summary["total_bp_counts"]
            )
        }
        inp = summary["input"]
        fmt = inp["file_format"]
        if inp["input_read"] == InputRead.PAIRED:
            fmt += ", Paired"
        else:
            fmt += ", Read {}".format(inp["input_read"])
        if inp["colorspace"]:
            fmt += ", Colorspace"
        if inp["interleaved"]:
            fmt += ", Interleaved"
        if inp["delivers_qualities"]:
            fmt += ", w/ Qualities"
        else:
            fmt += ", w/o Qualities"
        derived["input_format"] = fmt
        summary["derived"] = derived


def prettyprint_summary(summary: dict, outfile: Path = Path("summary.dump.txt")):
    """Pretty-print the summary to a file. Mostly used for debugging.
    """
    from pprint import pprint

    with open_(outfile, "w") as out:
        pprint(summary, out)
