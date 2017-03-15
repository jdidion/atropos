import copy
import importlib
import logging
from atropos.io.seqio import open_reader

class Pipeline(object):
    """Base class for analysis pipelines.
    """
    def __init__(self):
        self.record_counts = {}
        self.bp_counts = {}
    
    def __call__(self, reader, **kwargs):
        self.start(**kwargs)
        error = None
        try:
            for batch in reader:
                self.process_batch(batch)
        except Exception as error:
            pass
        finally:
            self.finish(**kwargs)
        return self.summarize(error=error)
    
    def start(self, **kwargs):
        pass
    
    def process_batch(self, batch):
        """Run the pipeline on a batch of records.
        
        Args:
            batch: A batch of reads. A batch has the format
            ({batch_metadata}, [records]).
        """
        batch_meta, records = batch
        context = batch_meta.copy()
        
        if not context['source'] in self.record_counts:
            self.record_counts[context['source']] = 0
        self.record_counts[context['source']] += context['size']
        
        if not context['source'] in self.bp_counts:
            self.bp_counts[context['source']] = [0, 0]
        context['bp'] = self.bp_counts[context['source']]
        
        self.add_to_context(context)
        self.handle_records(context, records)
    
    def add_to_context(self, context):
        """Add items to the batch context.
        """
        pass
    
    def handle_records(self, context, records):
        for record in records:
            self.handle_record(context, record)
    
    def handle_reads(self, context, read1, read2=None):
        raise NotImplementedError()
    
    def finish(self, **kwargs):
        pass
    
    def summarize(self, error=None):
        return dict(
            record_counts=self.record_counts,
            bp_counts=self.bp_counts)

class SingleEndPipelineMixin(object):
    def handle_record(self, context, record):
        context['bp'][0] += len(record)
        return self.handle_reads(context, record)

class PairedEndPipelineMixin(object):
    def handle_record(self, context, record):
        read1, read2 = record
        bp = context['bp']
        bp[0] += len(read1.sequence)
        bp[1] += len(read2.sequence)
        return self.handle_reads(context, read1, read2)

class BatchIterator(object):
    def __init__(self, reader, size, max_reads=None):
        self.reader = reader
        self.iterable = enumerate(reader, 1)
        self.size = size
        self.max_reads = max_reads
        self.batches = 0
        self.done = False
        self._empty_batch = [None] * size
        self._source = None
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.done:
            raise StopIteration()
        
        try:
            read_index, record = next(self.iterable)
        except:
            self.close()
            raise
        
        batch = copy.copy(self._empty_batch)
        batch[0] = record
        batch_index = 1
        max_size = self.size
        if self.max_reads:
            max_size = min(max_size, self.max_reads - read_index + 1)
        
        while batch_index < max_size:
            try:
                read_index, record = next(self.iterable)
                batch[batch_index] = record
                batch_index += 1
            except StopIteration:
                self.close()
                break
            except:
                self.close()
                raise
        
        if self.max_reads and read_index >= self.max_reads:
            self.close()
        
        self.batches += 1
        
        batch_meta = dict(
            index=self.batches,
            source=self._source,
            size=batch_index)
        
        if batch_index == self.size:
            return (batch_meta, batch)
        else:
            return (batch_meta, batch[0:batch_index])
    
    def close(self):
        self.done = True
        self.reader.close()

def execute_command(name, options):
    mod = importlib.import_module("atropos.commands.{}".format(name))
    return mod.execute(options)

def create_reader(options, counter_magnitude="M"):
    """Create sequence reader based on configured options.
    
    Args:
        options: Namespace-like object with configuration options.
        counter_magnitude: Magnitutde to use for progress bar.
    
    Returns:
        BatchIterator, possibly wrapped in progress bar.
    """
    interleaved = bool(options.interleaved_input)
    input1 = options.interleaved_input if interleaved else options.input1
    input2 = qualfile = None
    if options.paired and not interleaved:
        input2 = options.input2
    else:
        qualfile = options.input2
    
    reader = open_reader(input1, file2=input2, qualfile=qualfile,
        colorspace=options.colorspace, fileformat=options.format,
        interleaved=interleaved, single_input_read=options.single_input_read)
    
    qualities = reader.delivers_qualities
    
    # Wrap reader in subsampler
    if options.subsample:
        reader = subsample(reader, options.subsample)
    
    # Wrap reader in batch iterator
    batch_size = options.batch_size or 1000
    reader = BatchIterator(reader, batch_size, options.max_reads)
    
    # HACK: This is temporary until multi-file input is supported, at which
    # point the reader will keep track of the current source files.
    if input2:
        reader._source = (input1, input2)
    else:
        reader._source = input1
    
    # Wrap iterator in progress bar
    if options.progress:
        from atropos.io.progress import create_progress_reader
        reader = create_progress_reader(
            reader, options.progress, batch_size, options.max_reads,
            counter_magnitude)
    
    return (reader, (input1, input2), qualities, qualfile is not None)

def subsample(reader, frac):
    from random import random
    for reads in reader:
        if random() < frac:
            yield reads
