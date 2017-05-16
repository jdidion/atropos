#!/usr/bin/env python
"""Summarizes trimming accuracy for simulated reads."""

from argparse import ArgumentParser
import csv
import os
from common import *

# Either of these will work
from editdistance import eval as editdistance
#from Levenshtein import distance as editdistance

summary_fields = (
    "program", "program2", "threads", "dataset", "qcut", "retained reads", 
    "mismatch reads", "discarded reads", "total reads", "reads with adapters", 
    "retained reads with adapters", "non-adapter reads trimmed", 
    "adapter reads untrimmed", "adapter reads undertrimmed", 
    "adapter reads overtrimmed", "total ref bases", "total ref edit distance", 
    "total adapter bases", "total retained adapter bases", 
    "total adapter edit dist", "overtrimmed bases", "undertrimmed bases"
)

def summarize_accuracy(aln_iter, read_iter, w, read_length, adapters):
    adapter_lengths = [len(adapters[i]) for i in (0,1)]
    
    debug = False
    
    def summarize_alignment(a):
        ref_seq = a[5]
        if debug: print(ref_seq)
        ref_len = len(ref_seq)
        if debug: print(ref_len)
        read_seq = a[6]
        if debug: print(read_seq)
        read_len = len(read_seq)
        if debug: print(read_len)
        
        read_ref = read_seq[:ref_len]
        
        ref_ins = [i for i in range(len(ref_seq)) if ref_seq[i] == '-']
        expected_read = "".join(b for b in read_ref if b != '-')
        expected_read_len = len(expected_read)
        ref_del = ref_len - expected_read_len
        
        has_adapter = ref_len < read_len
        adapter_seq = []
        adapter_len = adapter_ins = adapter_del = polyA = 0
        if debug: print(has_adapter)
        
        if has_adapter:
            for b in read_seq[ref_len:]:
                if adapter_len >= adapter_lengths[i] and b == 'A':
                    polyA += 1
                else:
                    if b == '-':
                        adapter_del += 1
                    else:
                        adapter_seq.append(b)
                        adapter_len += 1
            
            adapter_ref_len = adapter_len + adapter_del
            if adapter_ref_len > adapter_lengths[i]:
                adapter_ins = adapter_ref_len - adapter_lengths[i]
        
        ref_edit_dist = editdistance(ref_seq, read_ref)
        adapter_edit_dist = editdistance(adapters[i][:adapter_len], "".join(adapter_seq))
        return [expected_read, expected_read_len,
            (ref_len, ref_edit_dist),
            (int(has_adapter), adapter_len, adapter_edit_dist, adapter_ins, adapter_del, polyA)
        ]
    
    cache = {}
    
    total_ref_bp = 0
    total_ref_edit_dist = 0
    total_adapter_bp = 0
    total_adapter_edit_dist = 0
    overtrimmed_bp = 0
    undertrimmed_bp = 0
    raw_trimmed_mismatch = 0
    num_adapter_reads = 0
    adapter_reads_untrimmed = 0
    adapter_reads_undertrimmed = 0
    adapter_reads_overtrimmed = 0
    non_adapter_reads_trimmed = 0
    num_aln = 0
    num_discarded_adapter_reads = 0
    num_discarded_adapter_bp = 0
    
    for num_reads, reads in enumerate(read_iter, 1):
        read_id = reads[0][0]
        aln = None
        
        assert read_id == reads[1][0], "Read IDs differ - {} != {}".format(read_id, reads[1][0])
        assert int(reads[0][1]) == 1 and int(reads[1][1]) == 2, "Mate identifiers are incorrect for {}".format(read_id)
        if read_id not in cache:
            for aln in aln_iter:
                num_aln += 1
                if read_id == aln[0][0]:
                    break
                else:
                    cache[aln[0][0]] = aln
        else:
            aln = cache.pop(read_id)
        
        if debug: print(reads)
        if debug: print(aln)
        
        if aln is None:
            raise Exception("No alignment for read {}".format(read_id))
        
        for i in (0, 1):
            expected_read, expected_read_len, ref_info, adapter_info = summarize_alignment(aln[i])
            total_ref_bp += ref_info[0]
            total_ref_edit_dist += ref_info[1]
            has_adapter = adapter_info[0] == 1
            total_adapter_bp += adapter_info[1]
            total_adapter_edit_dist += adapter_info[2]
            
            if has_adapter:
                num_adapter_reads += 1
            
            r = reads[i]
            trimmed_len = len(r[2])
            
            if debug: print(trimmed_len)
            if debug: print(r[2])
            if debug: print(expected_read)
            
            status = 'OK'
            common_len = min(trimmed_len, expected_read_len)
            if expected_read_len > trimmed_len:
                overtrimmed_bp += expected_read_len - trimmed_len
                status = 'OVERTRIMMED'
                if has_adapter:
                    adapter_reads_overtrimmed += 1
                else:
                    non_adapter_reads_trimmed += 1
            elif expected_read_len < trimmed_len:
                undertrimmed_bp += trimmed_len - expected_read_len
                status = 'UNDERTRIMMED'
                if not has_adapter:
                    raise Exception("Expected read length is less than trimmed "
                                    "length for read without an adapter: {}".format(read_id))
                elif trimmed_len == read_length:
                    adapter_reads_untrimmed += 1
                else:
                    adapter_reads_undertrimmed += 1
            if r[2][:common_len] != expected_read[:common_len]:
                # This happens with Skewer results due to automatic error correction
                raw_trimmed_mismatch += 1
                status = 'MISMATCH'
            
            w.writerow((read_id, i+1, expected_read_len, trimmed_len, status) + ref_info + adapter_info)

    # all remaining alignments represent reads that were discarded
    
    def handle_discarded(aln):
        read_id = aln[0][0]
        for i in (0, 1):
            expected_read, expected_read_len, ref_info, adapter_info = summarize_alignment(aln[i])
            w.writerow((read_id, i+1, expected_read_len, '', 'DISCARDED') + ref_info + adapter_info)
            if adapter_info[0] == 1:
                nonlocal num_discarded_adapter_reads, num_discarded_adapter_bp
                num_discarded_adapter_reads += 1
                num_discarded_adapter_bp += adapter_info[1]
    
    num_discarded = len(cache)
    for aln in cache.values():
        handle_discarded(aln)
    for aln in aln_iter:
        handle_discarded(aln)
        num_discarded += 1
    
    return (
        num_reads,
        raw_trimmed_mismatch,
        num_discarded,
        num_reads + num_discarded,
        num_adapter_reads + num_discarded_adapter_reads,
        num_adapter_reads,
        non_adapter_reads_trimmed,
        adapter_reads_untrimmed,
        adapter_reads_undertrimmed,
        adapter_reads_overtrimmed,
        total_ref_bp,
        total_ref_edit_dist,
        total_adapter_bp + num_discarded_adapter_bp,
        total_adapter_bp,
        total_adapter_edit_dist,
        overtrimmed_bp,
        undertrimmed_bp
    )

def main():
    parser = ArgumentParser()
    parser.add_argument('-a1', '--aln1', help=".aln file associated with simulated read1")
    parser.add_argument('-a2', '--aln2', help=".aln file associated with simulated read2")
    parser.add_argument('-r1', '--reads1', help="trimmed fastq file read1")
    parser.add_argument('-r2', '--reads2', help="trimmed fastq file read2")
    parser.add_argument('-l', '--read-length', type=int, default=125)
    parser.add_argument('-o', '--output', default='-')
    parser.add_argument('-s', '--summary', default='-')
    parser.add_argument('-t', '--table', default=None)
    parser.add_argument("-p", "--profile", default=None)
    parser.add_argument("--adapters", nargs=2, default=DEFAULT_ADAPTERS)
    parser.add_argument("--no-progress", action="store_true", default=False)
    args = parser.parse_args()
    
    profile = parse_profile(args.profile)
    
    with fileopen(args.aln1, 'rt') as a1, fileopen(args.aln2, 'rt') as a2:
        aln_pair_iterator = zip(aln_iterator(a1), aln_iterator(a2))
        
        with fileopen(args.reads1, 'rt') as r1, fileopen(args.reads2, 'rt') as r2:
            read_pair_iterator = zip(fq_iterator(r1, 1, True), fq_iterator(r2, 2, True))
            
            if not args.no_progress:
                try:
                    import tqdm
                    aln_pair_iterator = iter(tqdm.tqdm(aln_pair_iterator))
                except:
                    print("tqdm library is required for a progress bar")
            
            with fileopen(args.output, 'wt') as o:
                w = csv.writer(o, delimiter="\t")
                w.writerow((
                    'read_id','mate','expected_len','actual_len','status','has_adapter',
                    'adapter_len','adapter_edit_dist','adapter_ins','adapter_del','polyA'))
                summary = summarize_accuracy(
                    aln_pair_iterator, read_pair_iterator, w, args.read_length, 
                    args.adapters)
            
            with fileopen(args.summary, 'wt') as s:
                for field, value in zip(("{} " + field for field in summary_fields), summary):
                    print(field.format(value), file=s)
            
            if args.table:
                with fileopen(args.table, "wt") as t:
                    w = csv.writer(t, delimiter="\t")
                    w.writerow(profile + summary)

if __name__ == "__main__":
    main()
