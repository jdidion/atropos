# coding: utf-8
"""Quality trimming.
"""
# Import cythonized functions, defaulting to pure python implementations.
try:
    from ._qualtrim import quality_trim_index, nextseq_trim_index

except:
    import logging
    from atropos.util import qual2int
    
    logging.getLogger().debug("Import failed for cythonized qualtrim functions")
    
    def quality_trim_index(qualities, cutoff_front, cutoff_back, base=33):
        """Find the position at which to trim a low-quality end from a
        nucleotide sequence.
        
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
        score = 0
        max_qual = 0
        for idx in range(max_i):
            qual = qual2int(qualities[idx], base)
            score += cutoff_front - (qual - base)
            if score < 0:
                break
            if score > max_qual:
                max_qual = score
                start = idx + 1

        # same for 3' end
        max_qual = 0
        score = 0
        for idx in reversed(range(max_i)):
            qual = qual2int(qualities[idx], base)
            score += cutoff_back - (qual - base)
            if score < 0:
                break
            if score > max_qual:
                max_qual = score
                stop = idx
        
        if start >= stop:
            start, stop = 0, 0
        
        return (start, stop)
    
    def nextseq_trim_index(sequence, cutoff, base=33):
        """Variant of the above quality trimming routine that works on NextSeq
        data. With Illumina NextSeq, bases are encoded with two colors.
        'No color' (a dark cycle) usually means that a 'G' was sequenced, but
        that also occurs when sequencing falls off the end of the fragment. The
        read then contains a run of high-quality G bases in the end.

        This routine works as the one above, but counts qualities belonging to
        'G' bases as being equal to cutoff - 1.
        """
        bases = sequence.sequence
        qualities = sequence.qualities
        score = 0
        max_qual = 0
        max_i = len(qualities)
        for idx in reversed(range(max_i)):
            qual = qual2int(qualities[idx], base)
            if bases[idx] == 'G':
                qual = cutoff - 1
            score += cutoff - qual
            if score < 0:
                break
            if score > max_qual:
                max_qual = score
                max_i = idx
        return max_i
