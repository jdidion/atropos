"""Common classes/functions used in commands.
"""
from collections import Sequence
import copy
import platform
import sys
from atropos import __version__, AtroposError
from atropos.adapters import AdapterCache
from atropos.io.seqio import open_reader, sra_reader
from atropos.util import MergingDict, Const, Summarizable, Timing

class Pipeline(object):
    """Base class for analysis pipelines.
    """
    def __init__(self):
        self.record_counts = {}
        self.bp_counts = {}
    
    def __call__(self, command_runner, raise_on_error=False, **kwargs):
        self.start(**kwargs)
        try:
            for batch in command_runner.iterator():
                self.process_batch(batch)
        except Exception as err:
            if raise_on_error:
                raise
            else:
                command_runner.summary['exception'] = dict(
                    message=str(err),
                    details=sys.exc_info())
        finally:
            self.finish(command_runner.summary, **kwargs)
    
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
        for idx, record in enumerate(records):
            try:
                self.handle_record(context, record)
            except Exception as err:
                raise AtroposError(
                    "An error occurred at record {} of batch {}".format(
                    idx, context['index'])) from err
    
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
        total_bp_counts = tuple(sum(b) for b in zip(*self.bp_counts.values()))
        summary.update(
            record_counts=self.record_counts,
            total_record_count=sum(self.record_counts.values()),
            bp_counts=self.bp_counts,
            total_bp_counts=total_bp_counts,
            sum_total_bp_count=sum(total_bp_counts))

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

class Summary(MergingDict):
    """Contains summary information.
    """
    @property
    def has_exception(self):
        return 'exception' in self
    
    def finish(self):
        """Replaces Summarizable members with their summaries, and computes
        some aggregate values.
        """
        self._post_process_dict(self)
    
    def _post_process_dict(self, dict_val):
        if dict_val is None:
            return
        for key, value in tuple(dict_val.items()):
            if value is None:
                continue
            if isinstance(value, Summarizable):
                dict_val[key] = value = value.summarize()
            if isinstance(value, dict):
                self._post_process_dict(value)
            elif (
                    isinstance(value, Sequence) and
                    len(value) > 0 and
                    all(
                        val is None or isinstance(val, dict)
                        for val in value)):
                for val in value:
                    self._post_process_dict(val)
            else:
                if isinstance(value, Const):
                    dict_val[key] = value = value.value
                self._post_process_other(dict_val, key, value)
    
    def _post_process_other(self, parent, key, value):
        pass

class BaseCommandRunner(object):
    """Base class for command executors.
    
    Args:
        options: Command-line options.
    """
    def __init__(self, options, summary_class=Summary):
        self.options = options
        self.summary = summary_class()
        self.timing = Timing()
        self.return_code = None
        self.size = options.batch_size or 1000
        self.batches = 0
        self.done = False
        self._empty_batch = [None] * self.size
        self._progress_options = None
        
        if options.sra_reader:
            self.reader = reader = sra_reader(
                reader=options.sra_reader, 
                quality_base=options.quality_base, 
                colorspace=options.colorspace, 
                input_read=options.input_read,
                alphabet=options.alphabet)
            options.sra_reader = None
        else:
            interleaved = bool(options.interleaved_input)
            input1 = options.interleaved_input if interleaved else options.input1
            input2 = qualfile = None
            if options.paired and not interleaved:
                input2 = options.input2
            else:
                qualfile = options.input2
            self.reader = reader = open_reader(
                file1=input1, file2=input2, file_format=options.format, 
                qualfile=qualfile, quality_base=options.quality_base, 
                colorspace=options.colorspace, interleaved=interleaved, 
                input_read=options.input_read, alphabet=options.alphabet)
        
        # Wrap reader in subsampler
        if options.subsample:
            import random
            if options.subsample_seed:
                random.seed(options.subsample_seed)
            
            def subsample(reader, frac):
                """Generator that yields a random subsample of records.
                
                Args:
                    reader: The reader from which to sample.
                    frac: The fraction of records to yield.
                """
                for reads in reader:
                    if random.random() < frac:
                        yield reads
            
            reader = subsample(reader, options.subsample)
        
        self.iterable = enumerate(reader, 1)
        
        if options.progress:
            self._progress_options = (
                options.progress, self.size, self.max_reads,
                options.counter_magnitude)
        
        self.init_summary()
    
    def __getattr__(self, name):
        if hasattr(self.reader, name):
            return getattr(self.reader, name)
        elif hasattr(self.options, name):
            return getattr(self.options, name)
        else:
            raise ValueError("Unknown attribute: {}".format(name))
    
    def iterator(self):
        """Returns an iterator (an object with the __iter__ method) over
        input batches. BaseCommandRunner is itself an iterable, and will be
        wrapped with a progress bar if _progress_options is set.
        """
        if self._progress_options:
            # Wrap iterator in progress bar
            from atropos.io.progress import create_progress_reader
            itr = create_progress_reader(self, *self._progress_options)
            # itr may be none if there are no progress bar libraries available
            if itr is not None:
                return itr
        return self

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.done:
            raise StopIteration()
        
        try:
            read_index, record = next(self.iterable)
        except:
            self.finish()
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
                self.finish()
                break
            except:
                self.finish()
                raise
        
        if self.max_reads and read_index >= self.max_reads:
            self.finish()
        
        self.batches += 1
        
        batch_meta = dict(
            index=self.batches,
            # TODO: When multi-file input is supported, 'source' will need to
            # be the index of the current file/pair from which records are
            # being read.
            source=0,
            size=batch_index)
        
        if batch_index == self.size:
            return (batch_meta, batch)
        else:
            return (batch_meta, batch[0:batch_index])
    
    def init_summary(self):
        """Initialize the summary dict with general information.
        """
        self.summary['program'] = 'Atropos'
        self.summary['version'] = __version__
        self.summary['python'] = platform.python_version()
        self.summary['command'] = self.name
        self.summary['options'] = self.options.__dict__.copy()
        self.summary['timing'] = self.timing
        self.summary['sample_id'] = self.options.sample_id
        self.summary['input'] = self.reader.summarize()
        self.summary['input'].update(
            batch_size=self.size,
            max_reads=self.max_reads,
            batches=self.batches)
    
    def run(self):
        """Run the command, wrapping it in a Timing, catching any exceptions,
        and finally closing the reader.
        
        Returns:
            The tuple (retcode, summary).
        """
        with self.timing:
            try:
                self.return_code = self()
            except Exception as err: # pylint: disable=broad-except
                self.summary['exception'] = dict(
                    message=str(err),
                    details=sys.exc_info())
                self.return_code = 1
            finally:
                self.finish()
        
        return (self.return_code, self.summary)
    
    def __call__(self):
        """Execute the command. Must be implemented within the command
        module.
        
        Returns:
            The return code.
        """
        raise NotImplementedError()
    
    def finish(self):
        """Finish the command.
        """
        # Close the underlying reader.
        if not self.done:
            self.done = True
            self.reader.close()
        self.summary.finish()
    
    def load_known_adapters(self):
        """Load known adapters based on setting in command-line options.
        
        Args:
            options: Command-line options.
        """
        cache_file = None
        if self.options.cache_adapters:
            cache_file = self.options.adapter_cache_file
        adapter_cache = AdapterCache(cache_file)
        if adapter_cache.empty and self.options.default_adapters:
            adapter_cache.load_default()
        if self.options.known_adapter:
            for known in self.options.known_adapter:
                name, seq = known.split('=')
                adapter_cache.add(name, seq)
        if self.options.known_adapters_file:
            for known_file in self.options.known_adapters_file:
                adapter_cache.load_from_url(known_file)
        if self.options.cache_adapters:
            adapter_cache.save()
        return adapter_cache
