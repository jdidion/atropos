import sys

DEFAULT_ADAPTERS = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
]

class fileoutput(object):
    def __init__(self, path, mode='wt'):
        self.close = False
        if path == '-':
            self.fh = sys.stdout
        else:
            self.fh = open(path, mode)
            self.close = True
    
    def __enter__(self):
        return self.fh
    
    def __exit__(self, type, value, traceback):
        if self.close:
            self.fh.close()

def find_best_alignment(ref, query, side, min_match=1, cache=None, start=0, end=0.4, inc=0.01):
    from atropos import align
    best_match = None
    best_alternate = None
    
    for err in seq(start, end, inc):
        if cache is not None and err in cache:
            aligner = cache[err]
        else:
            aligner = make_aligner(ref, err, side)
            if cache is not None:
                cache[err] = aligner
        match = aligner.locate(query)
        if match is not None:
            ref_match = match[1] - match[0]
            if ref_match >= min_match:
                if match[5] == 0:
                    return match
                elif best_match is None or match[5] < best_match[5]:
                    best_match = match
            elif (best_alternate is None or ref_match > (best_alternate[1] - best_alternate[0]) or
                        (ref_match == (best_alternate[1] - best_alternate[0]) and match[5] < best_alternate[5])):
                best_alternate = match
    
    if best_match is not None:
        return best_match
    else:
        return best_alternate

def seq(start, end, inc):
    i = start
    while i <= end:
        yield i
        i += inc
