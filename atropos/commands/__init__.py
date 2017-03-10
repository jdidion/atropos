import importlib
from atropos.io.seqio import UnknownFileType, BatchIterator, open_reader

class Pipeline(object):
    """
    Args:
        paired: Whether each record will contain a pair of reads.
    """
    def __init__(self):
        self.record_counts = {}
        self.bp_counts = {}
    
    def __call__(self, reader, **kwargs):
        self.start(**kwargs)
        for batch in reader:
            self.process_batch(batch)
        return self.finish()
    
    def start(self, **kwargs):
        pass
    
    def process_batch(self, batch):
        """Run the pipeline on a batch of records.
        
        Args:
            batch: A batch of reads. A batch has the format
            (batch_source, batch_size, records).
        
        Returns:
            
        """
        batch_source, batch_size, records = batch
        if not batch_source in record_count:
            self.record_counts[batch_source] = 0
            self.bp_counts[batch_source] = [0, 0]
        self.record_counts[batch_source] += batch_size
        context = self.get_context(batch_num, batch_source)
        self.handle_records(context, records)
    
    def get_context(self, batch_num, batch_source):
        return dict(
            batch_num=batch_num,
            batch_source=batch_source,
            bp=self.bp_counts[batch_source])
    
    def handle_records(self, context, records):
        for record in records:
            self.handle_record(context, record)
    
    def handle_reads(self, context, read1, read2=None):
        raise NotImplementedError()
    
    def finish(self):
        pass
    
    def summarize(self):
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

def execute_command(name, options, parser):
    mod = importlib.import_module("atropos.commands.{}".format(name))
    return mod.execute(options, parser)

def create_reader(options, parser, counter_magnitude="M"):
    interleaved = bool(options.interleaved_input)
    input1 = options.interleaved_input if interleaved else options.input1
    input2 = qualfile = None
    if options.paired and not interleaved:
        input2 = options.input2
    else:
        qualfile = options.input2
    
    try:
        reader = open_reader(input1, file2=input2, qualfile=qualfile,
            colorspace=options.colorspace, fileformat=options.format,
            interleaved=interleaved, single_input_read=options.single_input_read)
    except (UnknownFileType, IOError) as e:
        parser.error(e)
    
    qualities = reader.delivers_qualities
    
    # Wrap reader in subsampler
    if options.subsample:
        reader = subsample(reader, options.subsample)
    
    # Wrap reader in batch iterator
    batch_size = options.batch_size or 1000
    reader = BatchIterator(reader, batch_size, options.max_reads)
    
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
