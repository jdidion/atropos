#!/usr/bin/env python
# Summarizes trimming accuracy for simulated reads.

from argparse import ArgumentParser
from atropos import align, xopen
import csv
from common import *

def aln_iterator(i):
    for line in i:
        if line[0] in ('@','#'):
            continue
        assert line[0] == '>'
        chrm, name, pos, strand = line[1:].rstrip().split('\t')
        if name.endswith('-1'):
            name = name[:-2]
        ref = next(i).rstrip()
        actual = next(i).rstrip()
        yield (name, chrm, pos, strand, ref, actual)

def fq_iterator(i, mate):
    for read in zip(*[i] * 4):
        name = read[0].rstrip()[1:]
        mate = name[-1]
        name = name[:-2]
        yield (name, mate, read[1].rstrip())

REF_ALIGN_FLAGS = align.START_WITHIN_SEQ2 | align.STOP_WITHIN_SEQ2 | align.START_WITHIN_SEQ1

def summarize_accuracy(a, r, w, adapters):
    adapter_aligners = [align.Aligner(adapters[i], 0.1) for i in (0,1)]
    adapter_lengths = [len(adapters[i]) for i in (0,1)]
    
    cache = {}
    for reads in r:
        read_id = reads[0][0]
        assert read_id == reads[1][0], "Read IDs differ - {} != {}".format(read_id, reads[1][0])
        assert int(reads[0][1]) == 1 and int(reads[1][1]) == 2, "Mate identifiers are incorrect for {}".format(read_id)
        if read_id not in cache:
            for aln in a:
                if read_id == aln[0][0]:
                    break
                else:
                    cache[aln[0][0]] = aln
        else:
            aln = cache.pop(read_id)
        
        try:
            for i in (0,1):
                aa = aln[i]
                rr = reads[i]
                
                # identify each part of the read
                ref_seq = aa[4]
                read_seq = aa[5]
                read_len = len(read_seq)
                ref_aligner = align.Aligner(ref_seq, 0.1, REF_ALIGN_FLAGS)
                ref_match = ref_aligner.locate(read_seq)
                if ref_match is None:
                    raise Exception("reference sequence doesn't align to read")
                if ref_match[1] < len(ref_seq):
                    raise Exception("reference sequence doesn't fully align to read")
                ref_end = ref_match[3]
                has_adapter = ref_end < read_len
                has_polyA = False
                if has_adapter:
                    adapter_match = adapter_aligners[i].locate(read_seq)
                    if adapter_match is None:
                        raise Exception("adapter sequence doesn't align")
                    if adapter_match[0] > 0 or adapter_match[1] < min(read_len - ref_end, adapter_lengths[i]):
                        raise Exception("adapter sequence doesn't fully align")
                    if adapter_match[2] != ref_end:
                        raise Exception("ambiguous adapter start position")
                    if has_polyA and not all(read_seq[b] == 'A' for b in range(adapter_match[3], read_len)):
                        raise Exception("base other than A found in tail sequence")
                
                trimmed_len = len(rr[2])
                if rr[2] != read_seq[:trimmed_len]:
                    raise Exception("mismatch between raw and trimmed read sequences")
                
                w.writerow((read_id, i+1, int(has_adapter), int(has_polyA), ref_end, trimmed_len, 'OK'))
                
            
        except Exception as e:
            print("Cannot evaluate {}: {}".format(read_id, e.args[0]))
            w.writerow((read_id, i+1, '', '', '', '', e.args[0]))
    
    # all remaining alignments represent reads that were discarded
    
    def handle_discarded(aln):
        read_id = aln[0][0]
        w.writerow((read_id, 1, '', '', '', '', 'discarded'))
        w.writerow((read_id, 2, '', '', '', '', 'discarded'))
        
    for aln in cache.values():
        handle_discarded(aln)
    for aln in a:
        handle_discarded(aln)

def main():
    parser = ArgumentParser()
    parser.add_argument('-a1', '--aln1', help=".aln file associated with simulated read1")
    parser.add_argument('-a2', '--aln2', help=".aln file associated with simulated read2")
    parser.add_argument('-r1', '--reads1', help="trimmed fastq file read1")
    parser.add_argument('-r2', '--reads2', help="trimmed fastq file read1")
    parser.add_argument('-o', '--output', default='-')
    parser.add_argument("--adapters", nargs=2, default=DEFAULT_ADAPTERS)
    args = parser.parse_args()
    
    with open(args.aln1, 'rt') as a1, open(args.aln2, 'rt') as a2:
        aln_pair_iterator = zip(aln_iterator(a1), aln_iterator(a2))
        
        with xopen.xopen(args.reads1, 'rt') as r1, xopen.xopen(args.reads2, 'rt') as r2:
            read_pair_iterator = zip(fq_iterator(r1, 1), fq_iterator(r2, 2))
            
            with fileoutput(args.output) as o:
                w = csv.writer(o, delimiter="\t")
                summarize_accuracy(aln_pair_iterator, read_pair_iterator, w, args.adapters)

if __name__ == "__main__":
    main()
