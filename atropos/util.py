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

MAGNITUDE = dict(
    G=(1E9, "000000000"),
    M=(1E6, "000000"),
    K=(1E3, "000")
)

def magnitude_formatter(magnitude):
    suffix = ""
    if magnitude is None:
        div = 1.0
    else:
        div = float(MAGNITUDE[magnitude.upper()][0])
        suffix = magnitude
    return lambda val: "{:.1f} {}".format(val / div, suffix)

def int_or_str(x):
    if x is None or isinstance(x, int):
        return x
    elif isinstance(x, str):
        x = x.upper()
        for a, mag in MAGNITUDE.items():
            x = x.replace(a, mag[1])
        return int(x)
    else:
        raise Exception("Unsupported type {}".format(x))

def enumerate_range(collection, start, end):
        'Generates an indexed series:  (0,coll[0]), (1,coll[1]) ...'
        i = start
        it = iter(collection)
        while i < end:
            yield (i, it.next())
            i += 1
