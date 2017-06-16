#!/usr/bin/env python
"""Summarize the mapping results of trimmed versus untrimmed reads.
This requires that both input bams be name-sorted.
"""
import argparse
from collections import defaultdict
from common import *
import csv
from glob import glob
import os

# There is the option to compute edit distance between the untrimmed and 
# trimmed reads. We did not use those metrics in the paper. If you want to 
# enable the edit distance calculation, you need to install the 'editdistance'
# pyton library.
#
# pip install editdistance
#
# # Note: if this doesn't work for you, you'll need to checkout the 
# editdistance repository and edit setup.py to enable the 'cythonize' command, 
# which will recompile the cython code for your local environment. So you would 
# run:
# 
# python setup.py build_ext -i && python setup.py install

nuc = ('A','C','G','T','N')

prog_header = ['prog', 'prog2', 'threads', 'dataset', 'qcut']

pair_fields = [
    ('prog', None),
    ('read_name', None),
    ('read_idx', None),
    ('failed', False),
    ('discarded', False),
    ('split', False),
    ('skipped', False),
    ('proper', False),
    ('valid', False)
]
read_fields = [
    ('mapped', False),
    ('quality', -1),
    ('score', -1),
    ('chrm', -1),
    ('start', -1),
    ('end', -1),
    ('clipped_front', -1),
    ('clipped_back', -1),
    ('trimmed_front', 0),
    ('trimmed_back', 0),
    ('in_region', -1)
]

def get_header():
    rf = [key for key, val in read_fields]
    r1 = ['read1_{}'.format(key) for key in rf]
    r2 = ['read2_{}'.format(key) for key in rf]
    return prog_header + [key for key, val in pair_fields[1:]] + r1 + r2

HEADER = get_header()

class TableRead(object):
    def __init__(self):
        for key, val in read_fields:
            setattr(self, key, val)
    
    def set_from_read(self, r):
        self.mapped = not r.is_unmapped
        self.quality = r.mapping_quality
        if r.has_tag('AS'):
            self.score = r.get_tag('AS')
        
        if self.mapped:
            self.chrm = r.reference_name
            self.start = r.reference_start
            self.end = r.reference_end
        
            def count_soft_clipped(seq):
                clipped = 0
                for pos in seq:
                    if pos:
                        break
                    else:
                        clipped += 1
                return clipped

            ref_pos = r.get_reference_positions(full_length=True)
            self.clipped_front = count_soft_clipped(ref_pos)
            self.clipped_back = count_soft_clipped(reversed(ref_pos))
    
    def set_trimming(self, u, t, use_edit_distance=True):
        untrimmed = u.query_sequence.upper()
        untrimmed_len = len(untrimmed)
        trimmed = t.query_sequence.upper()
        trimmed_len = len(trimmed)
        
        trimmed_front = 0 if use_edit_distance else -1
        if use_edit_distance and (untrimmed_len > trimmed_len):
            for i in range(untrimmed_len - trimmed_len + 1):
                if untrimmed[i:(i+trimmed_len)] == trimmed:
                    trimmed_front = i
                    break
            else:
                # Since Skewer performs automatic error correction, the trimmed and
                # untrimmed reads may not match, so in that case we find the closest
                # match by Levenshtein distance.
                dist = None
                for i in range(untrimmed_len - trimmed_len + 1):
                    d = editdistance.eval(untrimmed[i:(i+trimmed_len)], trimmed)
                    if not dist:
                        dist = d
                    elif d < dist:
                        trimmed_front = i
                        dist = d
        
        self.trimmed_front = trimmed_front
        self.trimmed_back = untrimmed_len - (trimmed_len + trimmed_front)
    
class TableRow(object):
    def __init__(self, prog, read_name, read_idx, regions=None):
        for key, val in pair_fields:
            setattr(self, key, val)
        self.prog = prog
        self.prog_fields = list(parse_profile(prog))
        self.read_name = read_name
        self.read_idx = read_idx
        self.read1 = TableRead()
        self.read2 = TableRead()
        self.regions = regions
    
    def skip(self, split=False, failed=False, discarded=False):
        self.skipped = True
        if split:
            self.split = True
        if failed:
            self.failed = True
        if discarded:
            self.discarded = True
        if self.regions and not discarded:
            self.regions.skip(self)
    
    def set_from_pair(self, r1, r2):
        assert r1.is_proper_pair == r2.is_proper_pair
        self.proper = r1.is_proper_pair
        # A valid pair is both proper (i.e. both reads are mapped)
        # and read1 and read2 are in opposite orientation
        self.valid = self.proper and (r1.is_reverse != r2.is_reverse)
        self.read1.set_from_read(r1)
        self.read2.set_from_read(r2)
        if self.regions:
            self.regions.set_overlaps(self)
    
    def set_trimming(self, u1, r1, u2, r2, use_edit_distance):
        self.read1.set_trimming(u1, r1, use_edit_distance)
        self.read2.set_trimming(u2, r2, use_edit_distance)
    
    def write(self, w):
        w.writerow(self.prog_fields + [
            getattr(self, key) for key, val in pair_fields[1:]
        ] + [
            getattr(self.read1, key) for key, val in read_fields
        ] + [
            getattr(self.read2, key) for key, val in read_fields
        ])

class Bed(object):
    def __init__(self, path, slop=200):
        from intervaltree import IntervalTree
        from csv import reader
        self.trees = {}
        self.slop = slop
        with open(path, "rt") as bed:
            chrm = None
            ivls = []
            for row in reader(bed, delimiter="\t"):
                if row[0] != chrm:
                    if len(ivls) > 0:
                        self.trees[chrm] = IntervalTree.from_tuples(ivls)
                    chrm = row[0]
                    ivls = []
                ivls.append((
                    max(1, int(row[1]) - slop + 1),
                    int(row[2]) + slop + 1
                ))
            if len(ivls) > 0:
                self.trees[chrm] = IntervalTree.from_tuples(ivls)
    
    def set_overlaps(self, read_pair, min_overlap=1.0):
        for read in (read_pair.read1, read_pair.read2):
            if read.mapped:
                read.in_region = self.overlaps(
                    read.chrm, read.start, read.end, min_overlap)

    def overlaps(self, chrm, start, end, min_overlap=1.0):
        if chrm not in self.trees:
            return 0
        
        tree = self.trees[chrm]
        hits = tree.search(start, end+1)
        
        if len(hits) == 0:
            return 0
        
        length = end - start + 1
        
        for h in hits:
            o = min(end, h.end) - max(start, h.begin) + 1
            l = min(length, h.end - h.begin + 1)
            if o / l >= min_overlap:
                return 1
        
        return 0
    
    def skip(self, read_pair):
        pass
    
    def close(self):
        pass

class Annotations(object):
    def __init__(self, bed_files_or_dir, pattern, untrimmed_name, trimmed_names):
        self.pattern = pattern
        self.annot_beds = {}
        if isinstance(bed_files_or_dir, str):
            self.bed_dir = bed_dir
            self.add(untrimmed_name)
            for name in trimmed_names:
                self.add(name)
        else:
            for path in bed_files_or_dir:
                basename = os.path.basename(path)
                if basename.startswith(untrimmed_name):
                    self.add(untrimmed_name, path)
                    continue
                for name in trimmed_names:
                    if basename.startswith(name):
                        self.add(name, path)
                        break
                else:
                    raise ValueError(
                        "No matching BAM file for BED path {}".format(path))
            if untrimmed_name not in self.annot_beds:
                raise ValueError("No untrimmed BED file")
    
    def add(self, name, path=None):
        if path:
            self.annot_beds[name] = AnnotBed(path)
        else:
            self.annot_beds[name] = AnnotBed(os.path.join(
                self.bed_dir, self.pattern.format(name=name)))
    
    def set_overlaps(self, read_pair, min_overlap=1.0):
        annot = self.annot_beds[read_pair.prog]
        annot.set_overlaps(read_pair, min_overlap)
    
    def skip(self, read_pair):
        annot = self.annot_beds[read_pair.prog]
        annot.skip(read_pair)
    
    def close(self):
        for annots in self.annot_beds.values():
            annots.close()

class AnnotBed(object):
    """
    The result of intersecting a name-sorted BAM file of read mappings with an
    annotation database:
    bam2bed --all-reads --do-not-sort < <bam> | cut -f 1-6 | bedmap \
      --delim '\t' --echo --echo-map-id - <annotations> > out.
    """
    def __init__(self, path):
        import csv
        self.fh = open(path, 'rt')
        self.reader = csv.reader(self.fh, delimiter="\t")
    
    def set_overlaps(self, read_pair, min_overlap=1.0):
        a1, a2 = self.next_annotations(read_pair)
        
        if read_pair.proper:
            # gene names should be semicolon-delimited in the 7th column
            g1 = set(a1[6].split(';'))
            g2 = set(a2[6].split(';'))
            
            # reads are considered "in region" if they are both mapped to
            # the same gene
            intersection = g1 & g2
            if len(intersection) > 0:
                read_pair.read1.in_region = read_pair.read2.in_region = 1
                return
        
        read_pair.read1.in_region = read_pair.read2.in_region = 0
    
    def skip(self, read_pair):
        self.next_annotations(read_pair)
    
    def next_annotations(self, read_pair):
        a1 = next(self.reader)
        a2 = next(self.reader)
        if a1[3] != read_pair.read_name or a2[3] != read_pair.read_name:
            raise Exception("Read {} does not match next annotation {} for sample {}".format(
                read_pair.read_name, a1[3], self.fh.name))
        
        return (a1, a2)
    
    def close(self):
        self.fh.close()

class Hist(object):
    def __init__(self, prog, n_bases):
        # position base histograms
        self.prog = prog
        self.prog_fields = list(parse_profile(prog))
        self.n_bases = n_bases
        self.hists = [
            [
                dict((base, [0] * n_bases) for base in ('A','C','G','T','N'))
                for side in range(2)
            ]
            for read in range(2)
        ]
    
    def add_reads(self, r1, r2):
        # add start and end bases to histogram
        for seq, hist in zip((r1.query_sequence, r2.query_sequence), (self.hists)):
            for i, b in enumerate_range(seq, 0, self.n_bases):
                hist[0][b][i] += 1
            for i, b in enumerate_range(reversed(seq), 0, self.n_bases):
                hist[1][b][i] += 1
    
    def write(self, w):
        for read, read_hists in enumerate(self.hists, 1):
            for side in range(2):
                for base in nuc:
                    for pos, count in enumerate(read_hists[side][base], 1):
                        w.writerow(self.prog_fields + [read, side, pos, base, count])

def summarize(untrimmed, trimmed, ow, hw, mode=None, regions=None,
              n_hist_bases=20, max_reads=None, use_edit_distance=True,
              progress=True):
    progs = ('untrimmed',) + tuple(trimmed.keys())
    hists = dict((prog, Hist(prog, n_hist_bases)) for prog in progs)
    itr = enumerate(untrimmed, 1)
    if progress:
        try:
            import tqdm
            itr = tqdm.tqdm(itr)
        except:
            print("tqdm library not available; not showing progress bar")
    
    for i, (name, u1, u2) in itr:
        assert len(u1) > 0 and len(u2) > 0
        
        rows = dict(untrimmed=TableRow('untrimmed', name, i, regions))
        valid = {}
        fail = False
        skip = False
        
        if len(u1) > 1 or len(u2) > 1:
            rows['untrimmed'].skip(split=True)
            skip = True
        else:
            u1 = u1[0]
            u2 = u2[0]
            if u1.is_qcfail or u2.is_qcfail:
                rows['untrimmed'].skip(failed=True)
                skip = fail = True
            else:
                hists['untrimmed'].add_reads(u1, u2)
                rows['untrimmed'].set_from_pair(u1, u2)
        
        for prog, bam in trimmed.items():
            rows[prog] = TableRow(prog, name, i, regions)
            
            if not bam.finished and bam.peek().query_name == name:
                _, t1, t2 = next(bam)
                assert len(t1) > 0 and len(t2) > 0
                if len(t1) > 1 or len(t2) > 1:
                    rows[prog].skip(split=True, failed=fail)
                elif skip:
                    rows[prog].skip(failed=fail)
                else:
                    t1 = t1[0]
                    t2 = t2[0]
                    hists[prog].add_reads(t1, t2)
                    valid[prog] = (t1, t2)
            else:
                rows[prog].skip(discarded=True, failed=fail)
        
        if len(valid) > 0:
            for prog, (r1, r2) in valid.items():
                row = rows[prog]
                row.set_from_pair(r1, r2)
                row.set_trimming(u1, r1, u2, r2, use_edit_distance)
        
        for prog in progs:
            rows[prog].write(ow)
        
        if max_reads and i >= max_reads:
            break
    
    for h in hists.values():
        h.write(hw)

def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(command=None)
    parser.add_argument("-i", "--bam-files", nargs="+", default=None)
    parser.add_argument("-d", "--bam-dir")
    parser.add_argument("-x", "--bam-extension", default=".name_sorted.bam")
    parser.add_argument("-p", "--bam-pattern", default=None)
    parser.add_argument("-u", "--untrimmed-name", default="untrimmed")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-H", "--hist", default="trimmed_hists.txt")
    parser.add_argument("-m", "--max-reads", type=int, default=None)
    parser.add_argument("--no-edit-distance", action="store_false",
                        default=True, dest="edit_distance",
                        help="Don't try to match by editdistance.")
    parser.add_argument("--no-progress", action="store_false", dest="progress",
                        default=True)
    sp = parser.add_subparsers()
    
    amplicon = sp.add_parser('amplicon')
    amplicon.set_defaults(command='amplicon')
    amplicon.add_argument(
        "-b", "--bed", default=None,
        help="Sorted .bed file of regions where reads should map.")
    amplicon.add_argument(
        "--min-overlap", type=float, default=1.0,
        help="When a .bed file is specified, this is the minimum "
            "fraction of a mapped read that must overlap a selected "
            "interval for that mapping to be considered valid. (1.0)")
    amplicon.add_argument(
        "--slop", type=int, default=None,
        help="When a .bed file is specified, this is the number of bp each "
            "region is extended. This is often necessary with amplicon data "
            "because enrichment can capture sequences that only paritally "
            "overlap the probes.")
    
    mrna = sp.add_parser('mrna')
    mrna.set_defaults(command='mrna')
    mrna.add_argument(
        "-b", "--bed-files", nargs="+", default=None,
        help="BED file names corresponding to BAM files.")
    mrna.add_argument(
        "-D", "--bed-dir", default=None,
        help="Directory where to find annotation bed files. Defaults to "
            "--bam-dir.")
    mrna.add_argument(
        "-B", "--bed-pattern", default='{name}.bed',
        help="String template for bed file names. \{name\} is replaced with the "
            "BAM file name with extension (--bam-extension) removed.")
    
    args = parser.parse_args()
    
    if args.edit_distance:
        import editdistance
    
    trimmed = {}
    untrimmed = None
    if args.bam_files:
        bam_files = args.bam_files
    else:
        pattern = (args.bam_pattern or "*{}").format(args.bam_extension)
        bam_files = list(glob(os.path.join(args.bam_dir, pattern)))
    for path in bam_files:
        name = os.path.basename(path)[:-len(args.bam_extension)]
        if name == args.untrimmed_name:
            untrimmed = BAMReader(path)
        else:
            trimmed[name] = BAMReader(path)
    if untrimmed is None:
        untrimmed = BAMReader(os.path.join(
            args.bam_dir,
            "{}{}".format(args.untrimmed_name, args.bam_extension)))
    
    regions = None
    if args.command == 'amplicon':
        regions = Bed(args.bed, args.slop or 200)
    elif args.command == 'mrna':
        regions = Annotations(
            args.bed_files or args.bed_dir or args.bam_dir,
            args.bed_pattern,
            args.untrimmed_name,
            trimmed.keys())
    
    try:
        with fileopen(args.output, 'wt') as o, fileopen(args.hist, 'wt') as h:
            ow = csv.writer(o, delimiter="\t")
            ow.writerow(HEADER)
            hw = csv.writer(h, delimiter="\t")
            hw.writerow(prog_header + ['read', 'side', 'pos', 'base', 'count'])
            summarize(untrimmed, trimmed, ow, hw,
                mode=args.command, regions=regions, max_reads=args.max_reads,
                use_edit_distance=args.edit_distance, progress=args.progress)
    finally:
        if untrimmed:
            untrimmed.close()
        for t in trimmed.values():
            t.close()
        if regions:
            regions.close()

if __name__ == "__main__":
    main()

# This is code to compare the mapped reads against the reference.
# from pysam import FastaFile
#parser.add_argument("-r", "--reference",
#    help="Reference sequence. If reads are bisulfite converted, then this must be a "
#         "bisulfite-converted reference.")
#parser.add_argument("-a", "--adapter1", default=ADAPTER1)
#parser.add_argument("-A", "--adapter2", default=ADAPTER2)
#parser.add_argument("-q", "--mapqs", nargs='*', default=(0,10,20,30,40))
# with FastaFile(args.reference) as ref,
# chrm = t.reference_name
# if bisulfite:
#     chrm = ('f','r')[int(t.is_rever)] + chrm
# chrm_len = t.get_reference_length(chrm)
# ref_start = t.reference_start - offset_front
# ref_end = t.reference_end + offset_back
# if ref_start < 0:
#     untrimmed = untrimmed[abs(ref_start):]
#     offset_front += ref_start
#     ref_start = 0
# if ref_end > chrm_len:
#     overhang = chrm_len - ref_end
#     untrimmed = untrimmed[:overhang]
#     offset_back += overhang
#     ref_end = chrm_len
#
# ref_seq = ref.fetch(chrm, ref_start, ref_end)
# assert len(ref_seq) == len(untrimmed)
#
# if bisulfite:
#     untrimmed = bisulfite_convert(untrimmed, r)
#     trimmed = bisulfite_convert(trimmed, r)
#
# def compare_seqs(s1, s2):
#     n = 0
#     for i, (b1, b2) in enumerate(zip(s1, s2)):
#         if b1 == 'N' or b2 == 'N':
#             n += 1
#         elif b1 == b2:
#             n = 0
#         else:
#             break
#     return i - n
#
# if undertrimmed_front == 0 and offset_start > 0:
#     overtrimmed_front = compare_seqs(untrimmed[offset_front:0:-1], ref_seq[offset_front:0:-1])
#
# if undertrimmed_back == 0 and offset_back > 0:
#     overtrimmed_back = compare_seqs(untrimmed[offset_back:], ref_seq[offset_back:]))
#
# def bisulfite_convert(seq, read):
#     if read == 0:
#         return seq.replace('C', 'T')
#     else:
#         return seq.replace('G', 'A')
#
# ADAPTER1="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG" # TruSeq index 7
# ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" # TruSeq universal
