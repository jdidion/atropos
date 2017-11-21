import progressbar
import progressbar.widgets
from atropos.util import magnitude_formatter

# TODO: put these configurations in the .pokrok file
# if max_reads:
#     reader = ProgressBarReader(reader, [
#         MagCounter(mag_format), " Reads (", progressbar.Percentage(), ") ",
#         progressbar.Timer(), " ", progressbar.Bar(),
#         progressbar.AdaptiveETA()
#     ], batch_size, max_reads)
# else:
#     reader = ProgressBarReader(reader, batch_size, [
#         MagCounter(mag_format), " Reads", progressbar.Timer(),
#         progressbar.AnimatedMarker()
#     ], batch_size)

class ProgressBarReader(progressbar.ProgressBar):
    """Extension of ProgressBar that supports starting and stopping the
    BatchReader.
    """
    def __init__(self, iterable, widgets, size=None, batch_size=1):
        super(ProgressBarReader, self).__init__(
            widgets=widgets,
            max_value=size or progressbar.UnknownLength)
        self.batch_size = batch_size
        self._iterable = iterable
        self.done = False

    def __next__(self):
        try:
            value = next(self._iterable)
            if self.start_time is None:
                self.start()
            self.update(self.value + (value[0] * self.batch_size))
            return value
        except StopIteration:
            self.close()
            raise

    def close(self):
        """Finish the progress bar and close the underlying iterator.
        """
        if not self.done:
            self.finish()
            self.done = True
        try:
            self._iterable.close()
        except:
            pass

class MagCounter(progressbar.widgets.WidgetBase):
    """Custom widget that formats the value using a specified magnitude.
    """
    def __init__(self, mag_format):
        super().__init__()
        self._format = mag_format

    def __call__(self, progress, data):
        return self._format(data["value"])