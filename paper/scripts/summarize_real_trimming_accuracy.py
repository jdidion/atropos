from common import *

nuc = ('A','C','G','T','N')

def summarize(bam, adapter1, adapter2, n_bases=20):
    # position base histograms
    read1_hist = [
        dict((base, [0] * n_bases) for base in ('A','C','G','T','N'))
        for i in range(2)
    ]
    read2_hist = [
        dict((base, [0] * n_bases) for base in ('A','C','G','T','N'))
        for i in range(2)
    ]
    
    def summarize_read(r, hist):
        mapped = not r.is_unmapped
        paired = r.is_paired
        proper = r.is_proper_pair
        mapq = r.mapping_quality
        trimmed_len = len(r.query_sequence)
        clipped_len = len(r.query_alignment_sequence)
        
        for i, b in enuermate_range(r.query_sequence.upper(), 0, n_bases):
            hist[b][i] += 1
        for i, b in enuermate_range(reversed(r.query_sequence.upper()), 0, n_bases):
            hist[b][i] += 1
    
    for r1, r2 in bam:
        summarize_read(r1, read1_hist)
        summarize_read(r2, read2_hist)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--adapter1", default=ADAPTER1)
    parser.add_argument("-A", "--adapter2", default=ADAPTER2)
    parser.add_argument("-b", "--bam")
    parser.add_argument("-o", "--output", default="-")
    paraser.add_argument("-h", "--hist", default="trimmed_hists.txt")
    args = parser.parse_args()
    
    with BAMReader(args.bam) as bam:
        summarize(bam, args.adapter1, args.adapter2)
    
    with open_output(args.output) as o:
        pass
    
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
