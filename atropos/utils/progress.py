from enum import Enum
from typing import Iterable, Iterator, Optional, TypeVar

import pokrok


T = TypeVar("T")


class ProgressBarType(Enum):
    BAR = "bar"
    MSG = "msg"
    ANY = "any"
    NONE = "none"


def create_progress_reader(
    itr: Iterable[T],
    progress_type: ProgressBarType = ProgressBarType.BAR,
    batch_size: int = 1,
    max_items: Optional[int] = None,
    **kwargs
) -> Iterator[T]:
    """
    Wraps an iterable in a progress bar of the specified type.

    Args:
        itr: The iterable to wrap.
        progress_type: msg = a custom progress bar that reports via log
            messages; bar = use a ProgressBar (from the progressbar library)
            or tqdm.
        max_items: Max number of items, if known in advance.
        batch_size: The number of records in each iterable item (iterable is
            typically a BatchReader).
        kwargs: Additional arguments to pass to the progress bar constructor.

    Returns:
        A wrapped iterable. If `progress_type == 'bar'` and neither of the
        supported libraries are available, a warning is logged and the unwrapped
        reader is returned.
    """
    if progress_type == ProgressBarType.NONE:
        return iter(itr)

    if progress_type == ProgressBarType.MSG:
        pokrok.set_plugins(["logging"])
    elif progress_type == ProgressBarType.BAR:
        pokrok.set_plugins(["tqdm", "progressbar2"])

    if max_items:
        widgets = [
            pokrok.Widget.COUNTER,
            pokrok.Widget.PERCENT,
            pokrok.Widget.ELAPSED,
            pokrok.Widget.BAR,
            pokrok.Widget.ETA,
        ]
    else:
        widgets = [pokrok.Widget.COUNTER, pokrok.Widget.ELAPSED, pokrok.Widget.SPINNER]

    return pokrok.progress_iter(
        itr,
        desc="Processed",
        size=max_items,
        unit="records",
        multiplier=batch_size,
        style=pokrok.Style(widgets),
        **kwargs,
    )
