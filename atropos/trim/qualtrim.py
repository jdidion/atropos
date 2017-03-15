# coding: utf-8
"""Quality trimming.
"""
import logging

# Import cythonized functions, defaulting to pure python implementations.
try:
    from atropos.trim._qualtrim import quality_trim_index, nextseq_trim_index
except:
    logging.getLogger().debug("Import failed for cythonized qualtrim functions")
    
    import sys
    from atropos.util import qual2int
    
    def quality_trim_index(qualities, cutoff_front, cutoff_back, base=33):
        """
        Find the position at which to trim a low-quality end from a nucleotide
        sequence.

        Qualities are assumed to be ASCII-encoded as chr(qual + base).

        The algorithm is the same as the one used by BWA within the function
        'bwa_trim_read':
        - Subtract the cutoff value from all qualities.
        - Compute partial sums from all indices to the end of the sequence.
        - Trim sequence at the index at which the sum is minimal.
        """
        start = 0
        stop = max_i = len(qualities)
        
        # find trim position for 5' end
        s = 0
        max_qual = 0
        for i in range(max_i):
            q = qual2int(qualities[i], base)
            s += cutoff_front - (q - base)
            if s < 0:
                break
            if s > max_qual:
                max_qual = s
                start = i + 1

        # same for 3' end
        max_qual = 0
        s = 0
        for i in reversed(range(max_i)):
            q = qual2int(qualities[i], base)
            s += cutoff_back - (q - base)
            if s < 0:
                break
            if s > max_qual:
                max_qual = s
                stop = i
        
        if start >= stop:
            start, stop = 0, 0
        
        return (start, stop)
    
    def nextseq_trim_index(sequence, cutoff, base=33):
        """
        Variant of the above quality trimming routine that works on NextSeq data.
        With Illumina NextSeq, bases are encoded with two colors. 'No color' (a
        dark cycle) usually means that a 'G' was sequenced, but that also occurs
        when sequencing falls off the end of the fragment. The read then contains
        a run of high-quality G bases in the end.

        This routine works as the one above, but counts qualities belonging to 'G'
        bases as being equal to cutoff - 1.
        """
        bases = sequence.sequence
        qualities = sequence.qualities
        s = 0
        max_qual = 0
        max_i = len(qualities)
        for i in reversed(range(max_i)):
            q = qual2int(qualities[i], base)
            if bases[i] == 'G':
                q = cutoff - 1
            s += cutoff - q
            if s < 0:
                break
            if s > max_qual:
                max_qual = s
                max_i = i
        return max_i
