from abc import ABCMeta, abstractmethod
from enum import Enum
from importlib import import_module
from pathlib import Path
from typing import Any, Callable, Dict, IO, Optional, Sequence, cast

from xphyle import STDOUT, STDERR, open_

from atropos.io import InputRead
from atropos.utils import classproperty
from atropos.utils.argparse import Namespace


class ReportWriter(metaclass=ABCMeta):
    def __init__(self, options: Namespace):
        self.options = options

    @property
    @abstractmethod
    def extensions(self) -> Sequence[str]:
        pass

    @abstractmethod
    def write_report(self, summary: dict, output_file: Path) -> None:
        pass

    def _get_output_file(
        self, output_file: Path, extensions: Optional[Sequence[str]] = None
    ) -> Path:
        if not extensions:
            extensions = self.extensions
        if output_file.suffix and output_file.suffix[1:] in extensions:
            return output_file
        else:
            return Path(str(output_file) + f".{extensions[0]}")


class BaseReportWriter(ReportWriter, metaclass=ABCMeta):
    def write_report(self, summary: dict, output_file: Path):
        with open_(self._get_output_file(output_file), "wt") as stream:
            self.serialize(summary, stream)

    @abstractmethod
    def serialize(self, summary: dict, stream: IO) -> None:
        pass


class TemplateReportWriter(BaseReportWriter):
    def __init__(
        self,
        name: str,
        options: Namespace,
        template_name: Optional[str] = None,
        template_paths: Sequence[Path] = None,
        template_globals: Optional[dict] = None
    ):
        """
        Args:
            name: Report format. Must match the file extension on a discoverable
                template.
            template_name: A template name to use, rather than auto-discover.
            template_paths: Sequence of paths to search for templates.
            template_globals: Dict of additional globals to add to the template
                environment.
        """
        super().__init__(options)
        self.name = name
        self.template_name = template_name
        self.template_paths = template_paths
        self.template_globals = template_globals

    @property
    def extensions(self) -> Sequence[str]:
        return self.name,

    def serialize(self, summary: dict, stream: IO):
        """
        Generate a report using Jinja2.
        """
        import jinja2

        template_name = self.template_name or f"template.{self.name}"
        template_paths = self.template_paths or []

        try:
            env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_paths))
            if self.template_globals:
                env.globals.update(self.template_globals)
            template = env.get_template(template_name)
        except:
            raise IOError(f"Could not load template file '{template_name}'")

        stream.write(template.render(summary=summary))


class DumpReportWriter(BaseReportWriter, metaclass=ABCMeta):
    def __init__(self, options: Namespace):
        super().__init__(options)

    @classproperty
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    def extensions(self) -> Sequence[str]:
        return self.name,

    def serialize(self, summary: dict, stream: IO):
        mod = import_module(self.name)
        mod.dump(summary, stream)


class PickleReportWriter(DumpReportWriter):
    @classproperty
    def name(self) -> str:
        return "pickle"


class DumpEncodedReportWriter(DumpReportWriter, metaclass=ABCMeta):
    def serialize(self, summary: dict, stream: IO):
        super().serialize(self._make_serializable(summary), stream)

    def _make_serializable(self, obj: Any) -> Any:
        if obj is None:
            return obj
        if isinstance(obj, (str, int, float, bool)):
            return obj
        if isinstance(obj, dict):
            return dict(
                (key, self._make_serializable(value))
                for key, value in cast(dict, obj).items()
            )
        if isinstance(obj, Sequence):
            return list(self._make_serializable(item) for item in cast(Sequence, obj))
        if isinstance(obj, Enum):
            return self._make_serializable(cast(Enum, obj).value)
        return str(obj)


class JsonReportWriter(DumpEncodedReportWriter):
    @classproperty
    def name(self) -> str:
        return "json"


class YamlReportWriter(DumpEncodedReportWriter):
    @classproperty
    def name(self) -> str:
        return "yaml"


class DefaultReportGenerator(metaclass=ABCMeta):
    """
    Base class for command report generators.
    """

    @classproperty
    def report_formats(cls) -> Dict[str, Callable[[Namespace], ReportWriter]]:
        return dict((cls.name, cls) for cls in (
            JsonReportWriter, YamlReportWriter, PickleReportWriter
        ))

    @classmethod
    def list_report_formats(cls) -> Sequence[str]:
        return tuple(sorted(cls.report_formats.keys()))

    @classproperty
    def default_report_formats(cls) -> Sequence[str]:
        return ()

    @classmethod
    def generate_reports(cls, summary: dict, options: Namespace):
        """
        Generates report(s) from a summary.

        Args:
            summary: The summary dict.
            options:
        """
        cls._add_derived_data(summary)

        report_file = options.report_file
        report_formats = options.report_formats

        if not report_formats:
            if report_file not in (STDOUT, STDERR) and report_file.suffix:
                report_formats = (report_file.suffix[1:],)
            else:
                report_formats = cls.default_report_formats

        for fmt in report_formats:
            writer = cls._create_report_writer(fmt, options)
            writer.write_report(summary, report_file)

    @classmethod
    def _add_derived_data(cls, summary: dict) -> None:
        """
        Adds to the summary some additional fields that are useful in the report.
        """
        derived = {
            "mean_sequence_lengths": tuple(
                None if bp is None else bp / summary["total_record_count"]
                for bp in summary["total_bp_counts"]
            )
        }
        inp = summary["input"]
        fmt = [inp["file_format"].name]
        if inp["input_read"] == InputRead.PAIRED:
            fmt.append("Paired")
        else:
            fmt.append(f"Read {inp['input_read']}")
        if inp["colorspace"]:
            fmt.append("Colorspace")
        if inp["interleaved"]:
            fmt.append("Interleaved")
        if inp["delivers_qualities"]:
            fmt.append("w/ Qualities")
        else:
            fmt.append("w/o Qualities")
        derived["input_format"] = ", ".join(fmt)
        summary["derived"] = derived

    @classmethod
    def _create_report_writer(cls, fmt: str, options: Namespace) -> ReportWriter:
        report_formats = cls.report_formats

        if fmt in report_formats:
            report_writer_class = report_formats[fmt]
            return report_writer_class(options)
        else:
            try:
                return TemplateReportWriter(fmt, options)
            except IOError:
                raise ValueError(f"Unsupported format {fmt}")


def prettyprint_summary(summary: dict, outfile: Path = Path("summary.dump.txt")):
    """
    Pretty-prints the summary to a file. Mostly used for debugging.
    """
    from pprint import pprint

    with open_(outfile, "w") as out:
        pprint(summary, out)
