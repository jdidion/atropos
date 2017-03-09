
class Pipeline(object):
    """
    Args:
        paired: Whether each record will contain a pair of reads.
    """
    def __init__(self):
        self.seen_batches = set()
        self.record_counts = {}
        self.bp_counts = {}
    
    def process_batch(self, batch):
        """Run the pipeline on a batch of records.
        
        Args:
            batch: A batch of reads. A batch has the format
            (batch_index, (batch_source, batch_size, records)).
        
        Returns:
            
        """
        batch_num, (batch_source, batch_size, records) = batch
        self.seen_batches.add(batch_num)
        if not batch_source in record_count:
            self.record_counts[batch_source] = 0
            self.bp_counts[batch_source] = [0, 0]
        self.record_counts[batch_source] += batch_size
        
        context = self.get_context(batch_num, batch_source)
        for record in records:
            self.handle_record(context, record)
        self.handle_results(context)
    
    def get_context(self, batch_num, batch_source):
        return dict(
            batch_num=batch_num,
            batch_source=batch_source,
            bp=self.bp_counts[batch_source])
    
    def handle_reads(self, context, read1, read2=None):
        raise NotImplementedError()
    
    def handle_results(self, context):
        pass

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

class QcPipeline(Pipeline):
    def __init__(self, **kwargs):
        super().__init__()
        self.stats = {}
        self.stats_kwargs = kwargs
    
    def get_stats(self, source):
        if source not in self.stats:
            self.stats[source] = ReadStatCollector(**self.stats_kwargs)
        return self.stats[source]

class SingleEndQcPipeline(QcPipeline, SingleEndPipelineMixin):
    def handle_reads(self, context, read):
        self.get_stats(context['batch_source']).collect(read)

class PairedEndQcPipeline(QcPipeline, PairedEndPipelineMixin):
    def handle_reads(self, context, read1, read2):
        src1, src2 = context['batch_source']
        self.get_stats(src1).collect(read1)
        self.get_stats(src2).collect(read2)

class TrimPipeline(Pipeline):
    def __init__(self, record_handler, result_handler):
        super().__init__()
        self.record_handler = record_handler
        self.result_handler = result_handler
    
    def get_context(self, batch_num, batch_source):
        context = super().get_context(batch_num, batch_source)
        context['results'] = defaultdict(lambda: [])
        return context
    
    def handle_reads(self, context, read1, read2):
        self.record_handler.handle_record(context, read1, read2)
    
    def handle_results(self, context):
        self.result_handler.write_result(context['batch_num'], context['results'])

class SerialTrimPipeline(TrimPipeline):
    def __init__(self, record_handler, writers):
        result_handler = WorkerResultHandler(WriterResultHandler(writers))
        super().__init__(record_handler, result_handler)

class RecordHandler(object):
    def __init__(self, modifiers, filters, formatters):
        self.modifiers = modifiers
        self.filters = filters
        self.formatters = formatters
    
    def handle_record(self, context, read1, read2=None):
        read1, read2 = self.modifiers.modify(read1, read2)
        dest = self.filters.filter(read1, read2)
        self.formatters.format(context['results'], dest, read1, read2)
        return (dest, read1, read2)

class StatsRecordHandlerWrapper(object):
    def __init__(self, record_handler, paired, mode='both', **kwargs):
        self.record_handler = record_handler
        self.paired = paired
        self.pre = self.post = None
        if mode in ('pre', 'both'):
            self.pre = {}
        if mode in ('post', 'both'):
            self.post = {}
        self.stats_kwargs = kwargs
    
    def handle_record(self, context, read1, read2=None):
        if self.pre:
            self.collect(self.pre, context['batch_source'], read1, read2)
        dest, read1, read2 = self.record_handler.handle_record(context, read1, read2)
        if self.post:
            if dest not in self.post:
                self.post[dest] = {}
            self.collect(self.post[dest], context['batch_source'], read1, read2)
        return (dest, read1, read2)
    
    def collect(self, stats, source, read1, read2=None):
        if paired:
            self._collect(stats, source[0], read1)
            self._collect(stats, source[1], read2)
        else:
            self._collect(stats, source, read1)
    
    def _collect(self, stats, source, read):
        if source not in self.stats:
            stats[source] = ReadStatCollector(**self.stats_kwargs)
        self.stats[source].collect(read)
