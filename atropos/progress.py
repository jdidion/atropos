import logging
from atropos.util import magnitude_formatter

def create_progress_reader(reader, progress_type="bar", batch_size=1, max_items=None,
                           counter_magnitude="M", **kwargs):
    mag_format = magnitude_formatter(counter_magnitude)
    
    if progress_type == "msg":
        return ProgressMessageReader(reader, batch_size, max_items, mag_format, **kwargs)
    elif progress_type != "bar":
        raise Exception("Unsupported progress type {}".format(progress_type))
            
    try:
        return create_progressbar_reader(reader, max_items, mag_format, **kwargs)
    except:
        pass
    
    try:
        return create_tqdm_reader(reader, max_items, **kwargs)
    except:
        pass
    
    logging.getLogger().warn("No progress bar library available")
    return reader

class ProgressMessageReader(object):
    def __init__(self, iterable, batch_size, interval=1000000, max_items=None, mag_format=None):
        self.iterable = iterable
        self.batch_size = batch_size
        self.interval = interval
        self.ctr = 0
        self.mat_format = mag_format
        if max_items:
            if mag_format:
                max_items = mag_format(max_items)
            else:
                max_items = str(max_items)
            self.msg = "Read {0}/" + max_items + " records in {1:.1f} seconds"
        else:
            self.msg = "Read {0} records in {1:.1f} seconds"
        
    def __next__(self):
        value = self.iterable.next()
        if value:
            self.ctr += value[0]
            if self.ctr % self.interval < self.batch_size:
                now = time.time()
                ctr = self.mag_format(self.ctr) if self.mag_format else self.ctr
                logging.getLogger().info(self.msg.format(ctr, now - self.start))
        return value
    
    next = __next__
    def __iter__(self):
        self.start = time.time()
        return self
    
    def close(self):
        logging.getLogger().info("Read a total of {} records".format(self.ctr))
        self.iterable.close()

def create_progressbar_reader(reader, max_reads=None, mag_format=None, values_have_size=True):
    import progressbar
    import progressbar.widgets
    import math

    class ProgressBarReader(progressbar.ProgressBar):
        def __init__(self, iterable, widgets, max_value=None, values_have_size=True):
            super(ProgressBarReader, self).__init__(
                widgets=widgets, max_value=max_value or progressbar.UnknownLength)
            self._iterable = iterable
            self.values_have_size = values_have_size
            self.done = False
        
        def __next__(self):
            try:
                value = next(self._iterable)
                if self.start_time is None:
                    self.start()
                if self.values_have_size:
                    self.update(self.value + value[0])
                else:
                    self.update(self.value + 1)
                return value
            except StopIteration:
                self.close()
                raise
        
        def close(self):
            if not self.done:
                self.finish()
                self.done = True
            try:
                self._iterable.close()
            except:
                pass
    
    class MagCounter(progressbar.widgets.WidgetBase):
        def __init__(self, mag_format):
            self._format = mag_format
    
        def __call__(self, progress, data):
            return self._format(data["value"])
        
    if max_reads:
        reader = ProgressBarReader(reader, [
            MagCounter(mag_format), " Reads (", progressbar.Percentage(), ") ",
            progressbar.Timer(), " ", progressbar.Bar(), progressbar.AdaptiveETA()
        ], max_reads, values_have_size=values_have_size)
    else:
        reader = ProgressBarReader(reader, [
            MagCounter(mag_format), " Reads", progressbar.Timer(),
            progressbar.AnimatedMarker()
        ], values_have_size=values_have_size)
    
    return reader

def create_tqdm_reader(reader, max_reads=None):
    import tqdm
    return tqdm.tqdm(reader)
