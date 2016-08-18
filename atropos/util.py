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
