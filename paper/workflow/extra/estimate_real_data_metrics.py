import argparse
from common import *
import re
from xphyle import xopen
import tqdm

ADAPTER1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG" # TruSeq index 7
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" # TruSeq universal

def estimate_metrics(fq1, fq2, adapter1, adapter2):
    adapter1_re = re.compile(adapter1 + '(A{8}|A*$)')
    adapter2_re = re.compile(adapter2 + '(A{8}|A*$)')
    adapter_len1 = len(adapter1)
    adapter_len2 = len(adapter2)
    
    cum_qual = 0.0
    cum_len = 0
    num_adapters1 = 0
    adapter_bases1 = 0
    num_adapters2 = 0
    adapter_bases2 = 0
    
    def qual2prob(qchar):
        q = ord(qchar) - 33
        return 10**(-q/10)
        
    for read1, read2 in tqdm.tqdm(zip(fq_iterator(fq1, 1), fq_iterator(fq2, 2))):
        cum_qual += sum(qual2prob(qchar) for qchar in read1[3])
        cum_qual += sum(qual2prob(qchar) for qchar in read2[3])
        cum_len += len(read1[2]) + len(read2[3])
        if adapter1_re.search(read1[2]):
            num_adapters1 += 1
            adapter_bases1 += adapter_len1
        if adapter2_re.search(read2[2]):
            num_adapters2 += 1
            adapter_bases2 += adapter_len2
    
    return (cum_qual / cum_len, num_adapters1, adapter_bases1, num_adapters2, adapter_bases2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--adapter1", default=ADAPTER1)
    parser.add_argument("-A", "--adapter2", default=ADAPTER2)
    parser.add_argument("-1", "--fastq1")
    parser.add_argument("-2", "--fastq2")
    parser.add_argument("-o", "--output", default="-")
    args = parser.parse_args()
    
    with xopen(args.fastq1) as fq1, xopen(args.fastq2) as fq2:
        metrics = estimate_metrics(fq1, fq2, args.adapter1, args.adapter2)
    
    with xopen(args.output, 'w') as o:
        print("Avg error prob: {}".format(metrics[0]), file=o)
        print("Read 1 with full-length adapters: {}".format(metrics[1]), file=o)
        print("Read 1 full-length adapter bases: {}".format(metrics[2]), file=o)
        print("Read 2 with full-length adapters: {}".format(metrics[3]), file=o)
        print("Read 2 full-length adapter bases: {}".format(metrics[4]), file=o)
            
if __name__ == "__main__":
    main()
