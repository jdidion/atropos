"""Implementation of python and system-level progress bar wrappers for
iterators.
"""
import logging
from pokrok.plugins import DefaultProgressMeterFactory, BaseProgressMeter
from atropos.util import Timestamp, magnitude_formatter

class AtroposMsgProgressMeterFactory(DefaultProgressMeterFactory):
    def __init__(self):
        super().__init__('atropos_msg', AtroposMsgProgressMeter)

    def iterate(
            self, iterable, size=None, widgets=None, desc=None, start=None,
            **kwargs):
        try:
            with self.create(size, widgets, desc, start, **kwargs) as pbar:
                for value in iterable:
                    yield value
                    if value:
                        pbar.increment(value[0]['size'])
        finally:
            iterable.close()

class AtroposMsgProgressMeter(BaseProgressMeter):
    """Wrapper around an iterable that logs an INFO message at a specified
    interval.

    Args:
        iterable: The iterable to wrap.
        size: Max number of items, if known in advance.
        batch_size: The number of records in each iterable item (iterable is
            typically a BatchReader).
        interval: The reporting interval.
        counter_magnitude: The magnitude to use when formatting numbers.
    """
    def __init__(
            self, size=None, widgets=None, desc=None, start=None,
            batch_size=None, interval=1000000, counter_magnitude=None,
            **kwargs):
        super().__init(size)
        self.batch_size = batch_size
        self.interval = interval
        self.ctr = 0
        self.mag_format = magnitude_formatter(counter_magnitude)
        self.start_time = None
        if self.size:
            if self.mag_format:
                max_items = self.mag_format(self.size)
            else:
                max_items = str(self.size)
            self.msg = "Read {0}/" + self.size + " records in {1:.1f} seconds"
        else:
            self.msg = "Read {0} records in {1:.1f} seconds"

    def increment(self, n=1):
        self.ctr += n
        if self.ctr % self.interval < self.batch_size:
            duration = Timestamp() - self.start
            ctr = self.ctr
            if self.mag_format:
                ctr = self.mag_format(ctr)
            logging.getLogger().info(
                self.msg.format(ctr, duration['wallclock']))

    def start(self):
        super().start()
        self.start_time = Timestamp()

    def finish(self):
        """Log the total number of records iterated and close the underlying
        iterable.
        """
        super().finish()
        logging.getLogger().info("Read a total of %s records", self.ctr)
