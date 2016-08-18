from argparse import ArgumentParser
from atropos.modifiers import AdapterCutter, SeqPurgeAdapterCutter
from atropos.adapters import AdapterParser
from atropos.seqio import open_reader

adapter_parser = AdapterParser()
A1 = adapter_parser.parse('AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG')
A2 = adapter_parser.parse('AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')

def create_cutadapt():
    ac1 = AdapterCutter([A1])
    ac2 = AdapterCutter([A2])
    return lambda read1,read2: (ac1(read1), ac2(read2))

def create_seqpurge():
    return SeqPurgeAdapterCutter(A1, A2)

def time_process(reader, aligner):
    from datetime import datetime
    start = datetime.now()
    i = run_process(reader, aligner)
    stop = datetime.now()
    print("Time to process {} read pairs: {}".format(i, stop-start))

def run_process(reader, aligner):
    for i, (read1, read2) in enumerate(reader, 1):
        read1_new, read2_new = aligner(read1, read2)
    return i

def main():
    parser = ArgumentParser()
    parser.add_argument("-a", "--aligner")
    parser.add_argument("-m", "--mode", choices=('time','run'), default='run')
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-1", "--fastq1", default='test1.fq')
    parser.add_argument("-2", "--fastq2", default='test2.fq')
    args = parser.parse_args()
    
    reader = open_reader(args.fastq1, args.fastq2)
    aligner = create_cutadapt() if args.aligner == 'cutadapt' else create_seqpurge()
    
    if args.mode == 'time':
        time_process(reader, aligner)
    else:
        run_process(reader, aligner)

if __name__ == "__main__":
    main()
