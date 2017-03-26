# coding: utf-8
"""
Colorspace conversion routines.

Inspired by agapython/util/Dibase.py from Corona lite,
but reimplemented to avoid licensing issues.

Encoding Table

  A C G T
A 0 1 2 3
C 1 0 3 2
G 2 3 0 1
T 3 2 1 0
"""
def _initialize_dicts():
    """
    Create the colorspace encoding and decoding dictionaries.
    """
    enc = {}
    for i, char1 in enumerate("ACGT"):
        enc['N' + char1] = '4'
        enc[char1 + 'N'] = '4'
        enc['.' + char1] = '4'
        enc[char1 + '.'] = '4'
        for j, char2 in enumerate("ACGT"):
            # XOR of nucleotides gives color
            enc[char1 + char2] = chr(ord('0') + (i ^ j))
    enc.update({ 'NN': '4', 'N.': '4', '.N': '4', '..': '4'})

    dec = {}
    for i, char1 in enumerate("ACGT"):
        dec['.' + str(i)] = 'N'
        dec['N' + str(i)] = 'N'
        dec[char1 + '4'] = 'N'
        dec[char1 + '.'] = 'N'
        for j, char2 in enumerate("ACGT"):
            # XOR of nucleotides gives color
            dec[char1 + chr(ord('0') + (i ^ j))] = char2
    dec['N4'] = 'N'

    return (enc, dec)

ENCODE, DECODE = _initialize_dicts()

def encode(nucs):
    """Given a sequence of nucleotides, convert them to colorspace. Only
    uppercase characters are allowed.
    
    Examples:
        >>> encode("ACGGTC")
        "A13012"
    """
    if not nucs:
        return nucs
    encoded = nucs[0:1]
    for idx in range(len(nucs) - 1):
        encoded += ENCODE[nucs[idx:idx+2]]
    return encoded

def decode(colors):
    """Decode a sequence of colors to nucleotide space. The first character in
    `colors` must be a nucleotide. Only uppercase characters are allowed.
    
    Examples:
        >>> decode("A13012")
        "ACGGTC"
    """
    if len(colors) < 2:
        return colors
    result = base = colors[0]
    for col in colors[1:]:
        base = DECODE[base + col]
        result += base
    return result
