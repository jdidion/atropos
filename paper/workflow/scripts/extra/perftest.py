from argparse import ArgumentParser
from atropos.modifiers import AdapterCutter, InsertAdapterCutter
from atropos.adapters import AdapterParser
from atropos.seqio import open_reader

adapter_parser = AdapterParser()
A1 = adapter_parser.parse('AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG')
A2 = adapter_parser.parse('AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')

def create_adapter():
    ac1 = AdapterCutter([A1])
    ac2 = AdapterCutter([A2])
    return lambda read1,read2: (ac1(read1), ac2(read2))

def create_insert():
    return InsertAdapterCutter(A1, A2)

def time_process(reader, aligner):
    from datetime import datetime
    start = datetime.now()
    i = run_process(reader, aligner)
    stop = datetime.now()
    t = stop-start
    print("Time to process {} read pairs: {}".format(i, t))
    return t.total_seconds()

def run_process(reader, aligner):
    for i, (read1, read2) in enumerate(reader, 1):
        read1_new, read2_new = aligner(read1, read2)
    return i

def main():
    parser = ArgumentParser()
    parser.add_argument("-a", "--aligner", choices=('adapter', 'insert'))
    parser.add_argument("-m", "--mode", choices=('time','run'), default='run')
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-1", "--fastq1", default='test1.fq')
    parser.add_argument("-2", "--fastq2", default='test2.fq')
    parser.add_argument("-i", "--iters", default=10)
    args = parser.parse_args()
    
    aligner = create_adapter() if args.aligner == 'adapter' else create_insert()
    
    if args.mode == 'time':
        t = 0
        for i in range(args.iters):
            reader = open_reader(args.fastq1, args.fastq2)
            t += time_process(reader, aligner)
            reader.close()
        print("Average time: {}".format(t / args.iters))
    else:
        run_process(reader, aligner)

if __name__ == "__main__":
    main()
