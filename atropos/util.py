import math

complement = {
    'A' : 'T',
    'C' : 'G',
    'R' : 'Y',
    'S' : 'S',
    'W' : 'W',
    'K' : 'M',
    'B' : 'V',
    'D' : 'H',
    'N' : 'N'
}
for k,v in list(complement.items()):
    complement[v] = k
    complement[k.lower()] = v.lower()
    complement[v.lower()] = k.lower()

def reverse_complement(seq):
    return "".join(complement[base] for base in reversed(seq))

log2 = math.log(2)
def sequence_complexity(seq):
    """A simple measure of sequence complexity"""
    seq = seq.upper()
    seqlen = float(len(seq))
    term = 0
    for base in ('A','C','G','T'):
        if base in seq:
            d = seq.count(base) / seqlen
            term += d * math.log(d) / log2
    return -term
