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
