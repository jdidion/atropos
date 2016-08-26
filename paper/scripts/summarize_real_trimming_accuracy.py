from common import *
import csv

nuc = ('A','C','G','T','N')

def summarize(bam, w, adapter1, adapter2, mapqs, n_bases=20):
    # position base histograms
    read1_hist = [
        dict((base, [0] * n_bases) for base in ('A','C','G','T','N'))
        for i in range(2)
    ]
    read2_hist = [
        dict((base, [0] * n_bases) for base in ('A','C','G','T','N'))
        for i in range(2)
    ]
    summaries = dict((q, [0, 0, 0, 0]))
    mapqs = sorted(mapqs, reverse=True)
    
    def summarize_read(r, hist):
        for i, b in enuermate_range(r.query_sequence.upper(), 0, n_bases):
            hist[b][i] += 1
        for i, b in enuermate_range(reversed(r.query_sequence.upper()), 0, n_bases):
            hist[b][i] += 1
        
        return (
            r.is_proper_pair,
            not r.is_unmapped,
            r.mapping_quality,
            len(r.query_sequence),
            len(r.query_alignment_sequence)
        )
    
    for r1, r2 in bam:
        r1_summary = summarize_read(r1, read1_hist)
        r2_summary = summarize_read(r2, read2_hist)
        proper = r1_summary[0]
        assert proper == r2_summary[0]
        w.writerow((r1.query_name, proper) + r1_summary[1:] + r2_summary[1:])
        
        minq = min(r1_summary[2], r2_summary[2])
        for q in mapqs:
            if minq >= q:
                summaries[q][0] += 1
                summaries[q][1] += int(proper)
                summaries[q][2] += r1_summary[3] - r1_summary[4]
                summaries[q][3] += r2_summary[3] - r2_summary[4]
    
    return (summaries, (read1_hist, read2_hist))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--adapter1", default=ADAPTER1)
    parser.add_argument("-A", "--adapter2", default=ADAPTER2)
    parser.add_argument("-b", "--bam")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-s", "--summary", default="-")
    paraser.add_argument("-h", "--hist", default="trimmed_hists.txt")
    args = parser.parse_args()
    
    mapqs = (0,10,20,30,40)
    
    with BAMReader(args.bam) as bam, open_output(args.output) as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(('read_id', 'proper_pair') + (
            'read{}_'.format(i) + field
            for field in ('mapped', 'mapq', 'trimmed_len', 'clipped_len')
            for i in (1,2)
        ))
        summaries, hists = summarize(bam, w, args.adapter1, args.adapter2, mapqs)
    
    with open_output(args.summary) as s:
        w = csv.writer(s, delimiter="\t")
        for q in mapqs:
            w.writerow([q] + summaries[q])
    
    with open_output(args.hist) as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(('read', 'side', 'pos', 'base', 'count'))
        for i, h in enumerate(hists, 1):
            for j in range(2):
                for b in nuc:
                    for k, count in enumerate(h[j][b], 1):
                        w.writerow((i, j, k, b, count))
            
if __name__ == "__main__":
    main()
