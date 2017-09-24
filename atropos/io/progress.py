"""Implementation of python and system-level progress bar wrappers for
iterators.
"""
import logging
from atropos.util import MAGNITUDE, Timestamp

class ProgressMessageReader(object):
    """Wrapper around an iterable that logs an INFO message at a specified
    interval.
    
    Args:
        iterable: The iterable to wrap.
        batch_size: The number of records in each iterable item (iterable is
            typically a BatchReader).
        interval: The reporting interval.
        max_items: Max number of items, if known in advance.
        mag_format: Function that formats an integer as a string with magnitude
            (e.g. 1000000 => 1M).
    """
    def __init__(
            self, iterable, batch_size, interval=1000000, max_items=None,
            mag_format=None):
        self.iterable = iterable
        self.batch_size = batch_size
        self.interval = interval
        self.ctr = 0
        self.mag_format = mag_format
        self.start = None
        if max_items:
            if mag_format:
                max_items = mag_format(max_items)
            else:
                max_items = str(max_items)
            self.msg = "Read {0}/" + max_items + " records in {1:.1f} seconds"
        else:
            self.msg = "Read {0} records in {1:.1f} seconds"
        
    def __next__(self):
        value = next(self.iterable)
        if value:
            self.ctr += value[0]
            if self.ctr % self.interval < self.batch_size:
                duration = Timestamp() - self.start
                ctr = self.ctr
                if self.mag_format:
                    ctr = self.mag_format(ctr)
                logging.getLogger().info(
                    self.msg.format(ctr, duration['wallclock']))
        return value
    
    next = __next__
    
    def __iter__(self):
        self.start = Timestamp()
        return self
    
    def close(self):
        """Log the total number of records iterated and close the underlying
        iterable.
        """
        logging.getLogger().info("Read a total of %s records", self.ctr)
        self.iterable.close()

def create_progress_reader(
        reader, progress_type="bar", batch_size=1, max_items=None,
        counter_magnitude="M", **kwargs):
    """Wrap an iterable in a progress bar of the specified type.
    
    Args:
        reader: The iterable to wrap.
        progress_type: msg = a custom progress bar that reports via log
            messages; bar = use a ProgressBar (from the progressbar library)
            or tqdm.
        max_items: Max number of items, if known in advance.
        mag_format: Function that formats an integer as a string with magnitude
            (e.g. 1000000 => 1M).
        batch_size: The number of records in each iterable item (iterable is
            typically a BatchReader).
        kwargs: Additional arguments to pass to the progress bar constructor.
    
    Returns:
        A wrapped iterable. If `progress_type == 'bar'` and neither of the
        supported libraries are available, a warning is logged and the unwrapped
        reader is returned.
    """
    mag_format = magnitude_formatter(counter_magnitude)
    
    if progress_type == "msg":
        return ProgressMessageReader(
            reader, batch_size, max_items=max_items, mag_format=mag_format,
            **kwargs)
            
    try:
        return create_progressbar_reader(
            reader, max_items, mag_format, **kwargs)
    except:
        pass
    
    try:
        return create_tqdm_reader(reader, max_items, **kwargs)
    except:
        pass
    
    logging.getLogger().warning("No progress bar library available")
    return reader

def magnitude_formatter(magnitude):
    """Returns a function that formats integers as magnitude strings.
    """
    suffix = ""
    if magnitude is None:
        div = 1.0
    else:
        div = float(MAGNITUDE[magnitude.upper()])
        suffix = magnitude
    return lambda val: "{:.1f} {}".format(val / div, suffix)

def create_progressbar_reader(reader, max_reads=None, mag_format=None):
    """Wrap an iterable in a ProgressBar.
    
    Args:
        max_reads: Max number of items, if known in advance.
        mag_format: Function that formats an integer as a string with magnitude
            (e.g. 1000000 => 1M).
    """
    import progressbar
    import progressbar.widgets
    
    class ProgressBarReader(progressbar.ProgressBar):
        """Extension of ProgressBar that supports starting and stopping the
        BatchReader.
        """
        def __init__(self, iterable, widgets, max_value=None):
            super(ProgressBarReader, self).__init__(
                widgets=widgets,
                max_value=max_value or progressbar.UnknownLength)
            self._iterable = iterable
            self.done = False
        
        def __next__(self):
            try:
                value = next(self._iterable)
                if self.start_time is None:
                    self.start()
                self.update(self.value + value[0]["size"])
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
        
    if max_reads:
        reader = ProgressBarReader(reader, [
            MagCounter(mag_format), " Reads (", progressbar.Percentage(), ") ",
            progressbar.Timer(), " ", progressbar.Bar(),
            progressbar.AdaptiveETA()
        ], max_reads)
    else:
        reader = ProgressBarReader(reader, [
            MagCounter(mag_format), " Reads", progressbar.Timer(),
            progressbar.AnimatedMarker()
        ])
    
    return reader

def create_tqdm_reader(reader, max_reads=None):
    """Wrap an iterable in a tqdm progress bar.
    
    Args:
        reader: The iterable to wrap.
        max_reads: Max number of items, if known in advance.
    
    Returns:
        The wrapped iterable.
    """
    import tqdm
    return tqdm.tqdm(reader, total=max_reads)
