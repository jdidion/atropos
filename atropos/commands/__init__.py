import importlib
from atropos.io.seqio import UnknownFileType, BatchIterator, open_reader

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
