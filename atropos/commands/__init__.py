"""Common classes/functions used in commands.
"""
import copy
import importlib
import logging
import platform
import random
import sys
from atropos import get_package_path, __version__
from atropos.adapters import AdapterCache
from atropos.io.seqio import open_reader
from atropos.util import MergingDict, Timing

class Pipeline(object):
    """Base class for analysis pipelines.
    """
    def __init__(self):
        self.record_counts = {}
        self.bp_counts = {}
    
    def __call__(self, reader, summary, raise_on_error=False, **kwargs):
        self.start(**kwargs)
        try:
            for batch in reader:
                self.process_batch(batch)
        except Exception as err:
            if raise_on_error:
                raise
            else:
                summary['error'] = dict(
                    message=str(err),
                    details=sys.exc_info())
        finally:
            self.finish(summary, **kwargs)
    
    def start(self, **kwargs):
        """Start the pipeline.
        """
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
        """Handle a sequence of records.
        
        Args:
            context: The pipeline context (dict).
            records: The sequence of records.
        """
        for record in records:
            self.handle_record(context, record)
    
    def handle_record(self, context, record):
        """Handle a single record.
        
        Args:
            context: The pipeline context (dict).
            record: The record.
        """
        raise NotImplementedError()
    
    def handle_reads(self, context, read1, read2=None):
        """Handle a read or read-pair.
        
        Args:
            context: The pipeline context (dict).
            read1, read2: The read pair; read2 will be None for single-end data.
        """
        raise NotImplementedError()
    
    def finish(self, summary, **kwargs):
        """Finish the pipeline, including adding information to the summary.
        
        Args:
            summary: Summary dict to update.
        """
        summary.update(
            record_counts=self.record_counts,
            total_record_count=sum(self.record_counts.values()),
            bp_counts=self.bp_counts,
            total_bp_counts=tuple(
                sum(b) for b in zip(*self.bp_counts.values())))

class SingleEndPipelineMixin(object):
    """Mixin for pipelines that implements `handle_record` for single-end data.
    """
    def handle_record(self, context, record):
        context['bp'][0] += len(record)
        return self.handle_reads(context, record)

class PairedEndPipelineMixin(object):
    """Mixin for pipelines that implements `handle_record` for paired-end data.
    """
    def handle_record(self, context, record):
        read1, read2 = record
        bps = context['bp']
        bps[0] += len(read1.sequence)
        bps[1] += len(read2.sequence)
        return self.handle_reads(context, read1, read2)

class BatchReader(object):
    """An iterator that yields batches of sequence records.
    
    Args:
        options: Namespace-like object with configuration options.
    """
    def __init__(self, options):
        interleaved = bool(options.interleaved_input)
        input1 = options.interleaved_input if interleaved else options.input1
        input2 = None
        qualfile = None
        if options.paired and not interleaved:
            input2 = options.input2
        else:
            qualfile = options.input2
        
        self.reader = reader = open_reader(
            input1, file2=input2, qualfile=qualfile,
            colorspace=options.colorspace, file_format=options.format,
            interleaved=interleaved,
            input_read=options.input_read)
        
        # Wrap reader in subsampler
        if options.subsample:
            reader = subsample(reader, options.subsample)
        
        self.iterable = enumerate(reader, 1)
        self.size = options.batch_size or 1000
        self.max_reads = options.max_reads
        self.batches = 0
        self.done = False
        self._empty_batch = [None] * self.size
        
        if options.progress:
            # Wrap iterator in progress bar
            from atropos.io.progress import create_progress_reader
            self.iterator = create_progress_reader(
                self, options.progress, self.size, self.max_reads,
                options.counter_magnitude)
        else:
            self.iterator = self
    
    def __getattr__(self, name):
        return getattr(self.reader, name)
    
    def __iter__(self):
        return self.iterator
    
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
            # HACK: Using input_names is temporary until multi-file input is
            # supported, at which point the reader will keep track of the
            # current source files.
            source=self.input_names,
            size=batch_index)
        
        if batch_index == self.size:
            return (batch_meta, batch)
        else:
            return (batch_meta, batch[0:batch_index])
    
    def close(self):
        """Close the underlying reader.
        """
        self.done = True
        self.reader.close()
    
    def summarize(self):
        summary = self.reader.summarize()
        summary.update(
            batch_size=self.size,
            max_reads=self.max_reads,
            batches=self.batches)
        return summary

def execute_command(name, options, orig_args):
    """Execute a subcommand. Loads module `name` within the atropos.commands
    package and calls that module's `execute` method.
    
    Args:
        options: Command-line options.
        summary: The summary dict.
    """
    retcode = 0
    summary = MergingDict()
    summary['program'] = 'Atropos'
    summary['version'] = __version__
    summary['python'] = platform.python_version()
    summary['command_line'] = [name] + list(orig_args)
    summary['options'] = options.__dict__.copy()
    summary['sample_id'] = options.sample_id
    
    reader = BatchReader(options)
    
    with Timing() as timing:
        try:
            mod = importlib.import_module("atropos.commands.{}".format(name))
            retcode = mod.execute(reader, options, summary)
        except Exception as err: # pylint: disable=broad-except
            summary['error'] = dict(
                message=str(err),
                details=sys.exc_info())
            retcode = 1
        finally:
            reader.close()
            summary['input'] = reader.summarize()
            summary['timing'] = timing.summarize()
    
    return (retcode, summary)

def subsample(reader, frac):
    """Generator that yields a random subsample of records.
    
    Args:
        reader: The reader from which to sample.
        frac: The fraction of records to yield.
    """
    for reads in reader:
        if random.random() < frac:
            yield reads

def load_known_adapters(options):
    """Load known adapters based on setting in command-line options.
    
    Args:
        options: Command-line options.
    """
    cache_file = options.adapter_cache_file if options.cache_adapters else None
    adapter_cache = AdapterCache(cache_file)
    if adapter_cache.empty and options.default_adapters:
        adapter_cache.load_default()
    if options.known_adapter:
        for known in options.known_adapter:
            name, seq = known.split('=')
            adapter_cache.add(name, seq)
    if options.known_adapters_file:
        for known_file in options.known_adapters_file:
            adapter_cache.load_from_url(known_file)
    if options.cache_adapters:
        adapter_cache.save()
    return adapter_cache
