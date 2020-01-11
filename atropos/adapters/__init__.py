# coding: utf-8
"""
Adapters
"""
import itertools
import logging
import os
import pickle
import re
import sys
from urllib.error import URLError
from urllib.request import urlopen
from atropos import align
from atropos.align import Match
from atropos.io.seqio import ColorspaceSequence, FastaReader
from atropos.util import (
    IUPAC_BASES, GC_BASES, MergingDict, NestedDict, CountingDict, Const,
    reverse_complement, ALPHABETS)
from atropos.util import colorspace as cs

class AdapterType(object):
    """Struct for adapter type information.

    Args:
        name: Adapter type name.
        desc: Adapter type description.
        flags: Alignment flags.
    """
    def __init__(self, name, desc, *flags):
        self.name = name
        self.desc = desc
        self.flags = flags[0]
        for i in range(1, len(flags)):
            self.flags |= flags[i]

    def asdict(self):
        """Returns AdapterType fields in a dict.
        """
        return dict(name=self.name, desc=self.desc, flags=Const(self.flags))

ADAPTER_TYPES = dict(
    back=AdapterType(
        'back', "regular 3'",
        align.START_WITHIN_SEQ2,
        align.STOP_WITHIN_SEQ2,
        align.STOP_WITHIN_SEQ1),
    front=AdapterType(
        'front', "regular 5'",
        align.START_WITHIN_SEQ2,
        align.STOP_WITHIN_SEQ2,
        align.START_WITHIN_SEQ1),
    prefix=AdapterType('prefix', "anchored 5'", align.STOP_WITHIN_SEQ2),
    suffix=AdapterType('suffix', "anchored 3'", align.START_WITHIN_SEQ2),
    anywhere=AdapterType('anywhere', "variable 5'/3'", align.SEMIGLOBAL),
    linked=AdapterType('linked', 'linked', 'linked')
)

def where_int_to_dict(where):
    """Convert a "where" flags integer to a dictionary of values from the
    corresponding AdapterType.
    """
    for adapter_type in ADAPTER_TYPES.values():
        if where == adapter_type.flags:
            return adapter_type.asdict()
    raise ValueError("Invalid WHERE value: {}".format(where))

# TODO get rid of these constants
BACK = ADAPTER_TYPES['back'].flags
FRONT = ADAPTER_TYPES['front'].flags
PREFIX = ADAPTER_TYPES['prefix'].flags
SUFFIX = ADAPTER_TYPES['suffix'].flags
ANYWHERE = ADAPTER_TYPES['anywhere'].flags
LINKED = ADAPTER_TYPES['linked'].flags

# TODO: specify this externally rather than hard-coding
DEFAULT_ADAPTERS_URL = "https://raw.githubusercontent.com/jdidion/atropos/master/atropos/adapters/sequencing_adapters.fa"
DEFAULT_ADAPTERS_PATH = os.path.join(
    os.path.dirname(__file__), 'sequencing_adapters.fa')

class AdapterParser(object):
    """Factory for Adapter classes that all use the same parameters (error rate,
    indels etc.). The given **kwargs will be passed to the Adapter constructors.

    Args:
        colorspace: Whether the adapter sequences are in colorspace.
        cache: The AdapterCache to use for fetching/storing named adapters.
    """
    def __init__(self, colorspace=False, cache=None, **kwargs):
        self.colorspace = colorspace
        self.cache = cache
        self.constructor_args = kwargs
        self.adapter_class = ColorspaceAdapter if colorspace else Adapter

    def parse(self, spec, cmdline_type='back'):
        """Parse an adapter specification not using ``file:`` notation and return
        an object of an appropriate Adapter class. The notation for anchored
        5' and 3' adapters is supported. If the name parameter is None, then
        an attempt is made to extract the name from the specification
        (If spec is 'name=ADAPTER', name will be 'name'.)

        Args:
            spec: The adapter spec.
            name: The adapter name. If not provided, one is automatically
                generated.
            cmdline_type: describes which commandline parameter was used (``-a``
                is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').

        TODO: describe the adapter spec format

        Returns:
            An :class:`Adapter` instance.
        """
        if spec.startswith('file:'):
            # read adapter sequences from a file
            with FastaReader(spec[5:]) as fasta:
                for record in fasta:
                    name = record.name.split(None, 1)[0]
                    yield self.parse_from_spec(
                        record.sequence, cmdline_type, name)
        else:
            yield self.parse_from_spec(spec, cmdline_type)

    def parse_from_spec(self, spec, cmdline_type='back', name=None):
        if cmdline_type not in ADAPTER_TYPES:
            raise ValueError(
                'cmdline_type cannot be {0!r}'.format(cmdline_type))
        orig_spec = spec
        where = ADAPTER_TYPES[cmdline_type].flags

        if name is None and spec is None:
            raise ValueError('Either name or spec must be given')
        elif name is None:
            if self.cache and self.cache.has_name(spec):
                name = spec
                spec = self.cache.get_for_name(name)
        elif spec is None:
            if self.cache and self.cache.has_name(name):
                spec = self.cache.get_for_name(name)

        if spec is None:
            raise ValueError('Name not found: {}'.format(name))
        elif name is None:
            name, spec = _extract_name_from_spec(spec)

        if self.cache and name is not None:
            self.cache.add(name, spec)

        front_anchored, back_anchored = False, False
        if spec.startswith('^'):
            spec = spec[1:]
            front_anchored = True
        if spec.endswith('$'):
            spec = spec[:-1]
            back_anchored = True

        sequence1, middle, sequence2 = spec.partition('...')

        if where == ANYWHERE:
            if front_anchored or back_anchored:
                raise ValueError("'anywhere' (-b) adapters may not be anchored")
            if middle == '...':
                raise ValueError("'anywhere' (-b) adapters may not be linked")
            return self.adapter_class(
                sequence=spec, where=where, name=name, **self.constructor_args)

        assert where == FRONT or where == BACK
        if middle == '...':
            if not sequence1:
                if where == BACK:  # -a ...ADAPTER
                    spec = sequence2
                else:  # -g ...ADAPTER
                    raise ValueError('Invalid adapter specification')
            elif not sequence2:
                if where == BACK:  # -a ADAPTER...
                    spec = sequence1
                    where = FRONT
                    front_anchored = True
                else:  # -g ADAPTER...
                    spec = sequence1
            else:
                # linked adapter
                if self.colorspace:
                    raise NotImplementedError(
                        'Using linked adapters in colorspace is not supported')
                # automatically anchor 5' adapter if -a is used
                if where == BACK:
                    front_anchored = True

                return LinkedAdapter(sequence1, sequence2, name=name,
                    front_anchored=front_anchored, back_anchored=back_anchored,
                    **self.constructor_args)

        if front_anchored and back_anchored:
            raise ValueError(
                'Trying to use both "^" and "$" in adapter specification '
                '{!r}'.format(orig_spec))
        if front_anchored:
            if where == BACK:
                raise ValueError("Cannot anchor the 3' adapter at its 5' end")
            where = PREFIX
        elif back_anchored:
            if where == FRONT:
                raise ValueError("Cannot anchor 5' adapter at 3' end")
            where = SUFFIX

        return self.adapter_class(
            sequence=spec, where=where, name=name, **self.constructor_args)

    def parse_multi(self, back=None, anywhere=None, front=None):
        """Parse all three types of commandline options that can be used to
        specify adapters. back, anywhere and front are lists of strings,
        corresponding to the respective commandline types (-a, -b, -g).

        Args:
            back: Back-adapter specs.
            anywhere: Anywhere-adapter specs.
            front: Front-adapter specs.

        Return:
            A list of appropriate Adapter classes.
        """
        adapters = []
        for specs, cmdline_type in (
                (back, 'back'), (anywhere, 'anywhere'), (front, 'front')):
            if not specs:
                continue
            for spec in specs:
                adapters.extend(self.parse(spec, cmdline_type))
        return adapters

class Adapter(object):
    """An adapter knows how to match itself to a read. In particular, it knows
    where it should be within the read and how to interpret wildcard characters.

    Args:
        sequence: The adapter sequence as string. Will be converted to
            uppercase. Also, Us will be converted to Ts.
        where: One of the BACK, FRONT, PREFIX, SUFFIX or ANYWHERE constants.
            This influences where the adapter is allowed to appear within in the
            read and also which part of the read is removed.
        max_error_rate: Maximum allowed error rate. The error rate is the number
            of errors in the alignment divided by the length of the part of the
            alignment that matches the adapter.
        min_overlap: Minimum length of the part of the alignment that
            matches the adapter.
        read_wildcards: Whether IUPAC wildcards in the read are allowed.
        adapter_wildcards: Whether IUPAC wildcards in the adapter are allowed.
        name: optional name of the adapter. If not provided, the name is set to
            a unique number.
        indels: Whether indels are allowed in adapter matches.
        indel_cost: Cost of indel bases in alignment scoring.
        match_probability: Callabale that computes the probability that a match
            could happen by random chance.
        max_rmp: Maximum random-match probability below which a match is
            considered genuine.
        gc_content: Expected GC content of sequences.
        alphabet: The Alphabet to use to validate the adapter sequence.
    """
    def __init__(
            self, sequence, where, max_error_rate=0.1, min_overlap=3,
            read_wildcards=False, adapter_wildcards=True, name=None,
            indels=True, indel_cost=1, match_probability=None, max_rmp=None,
            gc_content=0.5, alphabet=None):
        if len(sequence) == 0:
            raise ValueError("Empty adapter sequence")
        # TODO: all of this validation code should be wrapped up in Alphabet
        sequence = parse_braces(sequence.upper().replace('U', 'T'))
        seq_set = set(sequence)
        if seq_set <= set('ACGT'):
            adapter_wildcards = False
        if adapter_wildcards and not seq_set <= IUPAC_BASES:
            raise ValueError(
                "Invalid character(s) in adapter sequence: {}".format(
                    ','.join(seq_set - IUPAC_BASES)))
        if alphabet:
            if isinstance(alphabet, str):
                alphabet = ALPHABETS[alphabet]
            alphabet.validate_string(sequence)

        self.debug = False
        self.name = _generate_adapter_name() if name is None else name
        self.sequence = sequence
        self.where = where
        self.max_error_rate = max_error_rate
        self.min_overlap = min(min_overlap, len(self.sequence))
        self.match_probability = match_probability
        self.max_rmp = max_rmp
        self.gc_content = gc_content
        self.indels = indels
        self.adapter_wildcards = adapter_wildcards
        self.read_wildcards = read_wildcards
        # redirect trimmed() to appropriate function depending on adapter type
        trimmers = {
            FRONT: self._trimmed_front,
            PREFIX: self._trimmed_front,
            BACK: self._trimmed_back,
            SUFFIX: self._trimmed_back,
            ANYWHERE: self._trimmed_anywhere
        }
        self.trimmed = trimmers[where]
        if where == ANYWHERE:
            self._front_flag = None  # means: guess
        else:
            self._front_flag = where not in (BACK, SUFFIX)
        # statistics about length of removed sequences
        self.lengths_front = CountingDict()
        self.lengths_back = CountingDict()
        self.errors_front = NestedDict()
        self.errors_back = NestedDict()
        self.adjacent_bases = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, '': 0 }
        self.aligner = align.Aligner(
            self.sequence, self.max_error_rate, flags=self.where,
            wildcard_ref=self.adapter_wildcards,
            wildcard_query=self.read_wildcards)
        self.aligner.min_overlap = self.min_overlap
        if self.indels:
            self.aligner.indel_cost = indel_cost
        else:
            # TODO
            # When indels are disallowed, an entirely different algorithm
            # should be used.
            self.aligner.indel_cost = 100000

    def __repr__(self):
        return '<Adapter(name="{name}", sequence="{sequence}", where={where}, '\
            'max_error_rate={max_error_rate}, min_overlap={min_overlap}, '\
            'read_wildcards={read_wildcards}, '\
            'adapter_wildcards={adapter_wildcards}, '\
            'indels={indels})>'.format(**vars(self))

    def enable_debug(self):
        """Print out the dynamic programming matrix after matching a read to an
        adapter.
        """
        self.debug = True
        self.aligner.enable_debug()

    def match_to(self, read):
        """Attempt to match this adapter to the given read.

        Args:
            read: A :class:`Sequence` instance.

        Returns:
            A :class:`Match` instance if a match was found; return None if no
            match was found given the matching criteria (minimum overlap length,
            maximum error rate).
        """
        read_seq = read.sequence.upper()

        # try to find an exact match first unless wildcards are allowed
        pos = -1
        if not self.adapter_wildcards:
            if self.where == PREFIX:
                if read_seq.startswith(self.sequence):
                    pos = 0
            elif self.where == SUFFIX:
                if read_seq.endswith(self.sequence):
                    pos = (len(read_seq) - len(self.sequence))
            else:
                pos = read_seq.find(self.sequence)

        if pos >= 0:
            seqlen = len(self.sequence)
            return Match(
                0, seqlen, pos, pos + seqlen, seqlen, 0, self._front_flag,
                self, read)

        # try approximate matching
        if not self.indels and self.where in (PREFIX, SUFFIX):
            if self.where == PREFIX:
                alignment = align.compare_prefixes(
                    self.sequence, read_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards)
            else:
                alignment = align.compare_suffixes(
                    self.sequence, read_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards)
        else:
            alignment = self.aligner.locate(read_seq)
            if self.debug:
                print(self.aligner.dpmatrix)  # pragma: no cover

        if alignment:
            astart, astop, rstart, rstop, matches, errors = alignment
            size = astop - astart
            if ((
                    size >=
                    self.min_overlap and errors / size <=
                    self.max_error_rate
                ) and (
                    self.max_rmp is None or
                    self.match_probability(matches, size) <= self.max_rmp)):
                return Match(
                    astart, astop, rstart, rstop, matches, errors,
                    self._front_flag, self, read)

        return None

    def _trimmed_anywhere(self, match):
        """Trims an adapter from either the front or back of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        if match.front:
            return self._trimmed_front(match)
        else:
            return self._trimmed_back(match)

    def _trimmed_front(self, match):
        """Trims an adapter from the front of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        # TODO move away
        self.lengths_front[match.rstop] += 1
        self.errors_front[match.rstop][match.errors] += 1
        return match.read[match.rstop:]

    def _trimmed_back(self, match):
        """Trims an adapter from the back of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        self.lengths_back[len(match.read) - match.rstart] += 1
        self.errors_back[len(match.read) - match.rstart][match.errors] += 1
        adjacent_base = match.read.sequence[match.rstart-1:match.rstart]
        if adjacent_base not in 'ACGT':
            adjacent_base = ''
        self.adjacent_bases[adjacent_base] += 1
        return match.read[:match.rstart]

    def __len__(self):
        return len(self.sequence)

    def random_match_probabilities(self):
        """Estimate probabilities that this adapter matches a random sequence.
        Indels are not taken into account.

        Returns:
            A list of probabilities the same length as this adapter's sequence,
            where the value at position 'i' is the probability that i bases of
            this adapter match a random sequence.
        """
        if self._front_flag:
            seq = self.sequence[::-1]
        else:
            seq = self.sequence

        #matches = 0
        base_probs = (self.gc_content / 2.0, (1 - self.gc_content) / 2.0)
        probabilities = [1.0] + ([0] * len(seq))
        c_bases = frozenset(GC_BASES if self.adapter_wildcards else 'GC')

        # TODO: this doesn't work; need to figure out if RandomMatchProbability
        # can be used for this.
        #for idx, base in enumerate(seq, 1):
        #    if base in c_bases:
        #        matches += 1
        #    probabilities[idx] = self.match_probability(
        #        matches, idx, *base_probs)

        cur_p = 1.0
        for idx, base in enumerate(seq, 1):
            cur_p *= base_probs[0 if base in c_bases else 1]
            probabilities[idx] = cur_p
        return probabilities

    def summarize(self):
        """Summarize the activities of this :class:`Adapter`.
        """
        total_front = sum(self.lengths_front.values())
        total_back = sum(self.lengths_back.values())

        stats = MergingDict(
            adapter_class=self.__class__.__name__,
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back,
            match_probabilities=Const(self.random_match_probabilities()))

        where = self.where
        assert (
            where in (ANYWHERE, LINKED) or
            (where in (BACK, SUFFIX) and total_front == 0) or
            (where in (FRONT, PREFIX) and total_back == 0))

        stats["where"] = where_int_to_dict(where)
        stats["sequence"] = Const(self.sequence)
        stats["max_error_rate"] = Const(self.max_error_rate)
        if where in (ANYWHERE, FRONT, PREFIX):
            stats["lengths_front"] = self.lengths_front
            stats["errors_front"] = self.errors_front
        if where in (ANYWHERE, BACK, SUFFIX):
            stats["lengths_back"] = self.lengths_back
            stats["errors_back"] = self.errors_back
        if where in (BACK, SUFFIX):
            stats["adjacent_bases"] = self.adjacent_bases

        return stats

class ColorspaceAdapter(Adapter):
    """An adapter for a colorspace sequence.

    Args:
        args, kwargs: Arguments to pass to :class:`Adapter` constructor.
    """
    def __init__(self, *args, **kwargs):
        if kwargs.get('adapter_wildcards', False):
            raise ValueError("Wildcards not supported for colorspace adapters")
        kwargs['adapter_wildcards'] = False
        super().__init__(*args, **kwargs)
        has_nucleotide_seq = False
        if set(self.sequence) <= set('ACGT'):
            # adapter was given in basespace
            self.nucleotide_sequence = self.sequence
            has_nucleotide_seq = True
            self.sequence = cs.encode(self.sequence)[1:]
        if self.where in (PREFIX, FRONT) and not has_nucleotide_seq:
            raise ValueError(
                "A 5' colorspace adapter needs to be given in nucleotide space")
        self.aligner.reference = self.sequence

    def match_to(self, read):
        """Attempt to match this adapter to the given read.

        Args:
            read: A :class:`Sequence` instance.

        Returns:
            A :class:`Match` instance if a match was found; return None if no
            match was found given the matching criteria (minimum overlap length,
            maximum error rate).
        """
        if self.where != PREFIX:
            return super().match_to(read)
        # create artificial adapter that includes a first color that encodes the
        # transition from primer base into adapter
        asequence = (
            cs.ENCODE[read.primer + self.nucleotide_sequence[0:1]] +
            self.sequence)

        pos = 0 if read.sequence.startswith(asequence) else -1
        if pos >= 0:
            match = Match(
                0, len(asequence), pos, pos + len(asequence),
                len(asequence), 0, self._front_flag, self, read)
        else:
            # try approximate matching
            self.aligner.reference = asequence
            alignment = self.aligner.locate(read.sequence)
            if self.debug:
                print(self.aligner.dpmatrix)  # pragma: no cover
            if alignment is not None:
                match = Match(*(alignment + (self._front_flag, self, read)))
            else:
                match = None

        if match is None:
            return None
        assert (
            match.length > 0 and
            match.errors / match.length <= self.max_error_rate)
        assert match.length >= self.min_overlap
        return match

    def _trimmed_front(self, match):
        """Trims an adapter from the front of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        read = match.read
        self.lengths_front[match.rstop] += 1
        self.errors_front[match.rstop][match.errors] += 1
        # to remove a front adapter, we need to re-encode the first color
        # following the adapter match
        color_after_adapter = read.sequence[match.rstop:match.rstop + 1]
        if not color_after_adapter:
            # the read is empty
            return read[match.rstop:]
        base_after_adapter = cs.DECODE[
            self.nucleotide_sequence[-1:] + color_after_adapter]
        new_first_color = cs.ENCODE[read.primer + base_after_adapter]
        new_read = read[:]
        new_read.sequence = new_first_color + read.sequence[(match.rstop + 1):]
        new_read.qualities = None
        if read.qualities:
            new_read.qualities = read.qualities[match.rstop:]
        return new_read

    def _trimmed_back(self, match):
        """Trims an adapter from the back of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        # trim one more color if long enough
        adjusted_rstart = max(match.rstart - 1, 0)
        self.lengths_back[len(match.read) - adjusted_rstart] += 1
        self.errors_back[len(match.read) - adjusted_rstart][match.errors] += 1
        return match.read[:adjusted_rstart]

    def __repr__(self):
        return '<ColorspaceAdapter(sequence={0!r}, where={1})>'.format(
            self.sequence, self.where)

# TODO Consolidate Match and a LinkedMatch class.

class LinkedMatch(object):
    """Represent a match of a LinkedAdapter.

    Args:
        front_match: The match to the front of the sequence.
        back_match: The match to the back of the sequence.
        adapter: The matched adapter.
    """
    def __init__(self, front_match, back_match, adapter):
        self.front_match = front_match
        self.back_match = back_match
        self.adapter = adapter
        assert front_match is not None

    def get_info_record(self):
        """Returns the info record for the either the back or forward match.
        """
        if self.back_match:
            return self.back_match.get_info_record()
        else:
            return self.front_match.get_info_record()

class LinkedAdapter(object):
    """An adapter with linked front and back sequences.

    Args:
        front_sequence: Front adapter sequence.
        back_sequence: Back adapter sequence.
        front_anchored: Whether the front adapter is anchored.
        back_anchored: Whether the back adapter is anchored.
        name: Adapter name.
        kwargs: passed on to individual :class:`Adapter` constructors.
    """
    def __init__(
            self, front_sequence, back_sequence, front_anchored=True,
            back_anchored=False, name=None, **kwargs):
        assert front_anchored and not back_anchored
        where1 = PREFIX if front_anchored else FRONT
        where2 = SUFFIX if back_anchored else BACK
        self.front_anchored = front_anchored
        self.back_anchored = back_anchored

        # The following attributes are needed for the report
        self.where = LINKED
        self.name = _generate_adapter_name() if name is None else name
        self.front_adapter = Adapter(
            front_sequence, where=where1, name=None, **kwargs)
        self.back_adapter = Adapter(
            back_sequence, where=where2, name=None, **kwargs)

    def enable_debug(self):
        """Enable debug on adapters.
        """
        self.front_adapter.enable_debug()
        self.back_adapter.enable_debug()

    def match_to(self, read):
        """Match the linked adapters against the given read. If the 'front'
        adapter is not found, the 'back' adapter is not searched for.

        Args:
            read: A :class:`Sequence` instance.

        Returns:
            A :class:`Match` instance if a match was found; return None if no
            match was found given the matching criteria (minimum overlap length,
            maximum error rate).
        """
        front_match = self.front_adapter.match_to(read)
        if front_match is None:
            return None
        # TODO use match.trimmed() instead as soon as that does not update
        # statistics anymore
        read = read[front_match.rstop:]
        back_match = self.back_adapter.match_to(read)
        return LinkedMatch(front_match, back_match, self)

    def trimmed(self, match):
        """Returns the read trimmed with the front and/or back adapter
        trimmer(s).

        Args:
            match: The match to trim.

        Returns:
            The trimmed read.
        """
        front_trimmed = self.front_adapter.trimmed(match.front_match)
        if match.back_match:
            return self.back_adapter.trimmed(match.back_match)
        else:
            return front_trimmed

    def summarize(self):
        """Returns the summary dict for this adapter.
        """
        total_front = sum(self.front_adapter.lengths_front.values())
        total_back = sum(self.back_adapter.lengths_back.values())

        stats = MergingDict(
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back)

        where = self.where
        assert (
            where in (ANYWHERE, LINKED) or
            (where in (BACK, SUFFIX) and total_front == 0) or
            (where in (FRONT, PREFIX) and total_back == 0))

        stats["where"] = where_int_to_dict(where)
        stats["front_sequence"] = Const(self.front_adapter.sequence)
        stats["front_match_probabilities"] = Const(
            self.front_adapter.random_match_probabilities())
        stats["back_sequence"] = Const(self.back_adapter.sequence)
        stats["back_match_probabilities"] = Const(
            self.back_adapter.random_match_probabilities())
        stats["front_max_error_rate"] = Const(self.front_adapter.max_error_rate)
        stats["back_max_error_rate"] = Const(self.back_adapter.max_error_rate)
        stats["front_lengths_front"] = self.front_adapter.lengths_front
        stats["front_lengths_back"] = self.front_adapter.lengths_back
        stats["back_lengths_front"] = self.back_adapter.lengths_front
        stats["back_lengths_back"] = self.back_adapter.lengths_back
        # have to clone these nested dicts and set them
        # up with a custom merge function
        stats["front_errors_front"] = self.front_adapter.errors_front
        stats["front_errors_back"] = self.front_adapter.errors_back
        stats["back_errors_front"] = self.back_adapter.errors_front
        stats["back_errors_back"] = self.back_adapter.errors_back

        return stats

class AdapterCache(object):
    """Cache for known adapters.

    Args:
        path: Path to file where cache is written.
        auto_reverse_complement: Whether adapter reverse-complements should
            automatically be added to the cache.
    """
    def __init__(self, path=".adapters", auto_reverse_complement=False):
        self.path = path
        self.auto_reverse_complement = auto_reverse_complement
        if path and os.path.exists(path):
            with open(path, "rb") as cache:
                try:
                    self.seq_to_name, self.name_to_seq = pickle.load(cache)
                    return
                except:
                    # It is possible for the cache file to get corrupted - see #71
                    pass

        self.seq_to_name = {}
        self.name_to_seq = {}

    @property
    def empty(self):
        """Whether the cache is empty.
        """
        return len(self.seq_to_name) == 0

    def save(self):
        """Save the cache to file.
        """
        if self.path is not None:
            with open(self.path, "wb") as cache:
                pickle.dump((self.seq_to_name, self.name_to_seq), cache)

    def add(self, name, seq):
        """Add a sequence to the cache.

        Args:
            name: Adapter name.
            seq: Adapter sequence.
        """
        self._add(name, seq)
        if self.auto_reverse_complement:
            self._add("{}_rc".format(name), reverse_complement(seq))

    def _add(self, name, seq):
        if seq not in self.seq_to_name:
            self.seq_to_name[seq] = set()
        self.seq_to_name[seq].add(name)
        self.name_to_seq[name] = seq

    def load_from_file(self, path=DEFAULT_ADAPTERS_PATH):
        """Load cached data from a file.

        Args:
            path: Path from which to load.
        """
        with open(path, "rt") as i:
            return self.load_from_fasta(i)

    def load_from_url(self, url=DEFAULT_ADAPTERS_URL):
        """Load adapter data from a URL.

        Args:
            url: URL from which to load.
        """
        logging.getLogger().info(
            "Loading list of known contaminants from %s", url)
        try:
            fasta = urlopen(url).read().decode().split("\n")
            return self.load_from_fasta(fasta)
        except URLError:
            if url.startswith("file:"):
                url = url[5:]
            return self.load_from_file(url)

    def load_from_fasta(self, fasta):
        """Load adapter data from a FASTA file.

        Args:
            fasta: FASTA file.
        """
        close = False
        if isinstance(fasta, str):
            fasta = open(fasta, 'rt')
            close = True
        num_records = None
        with FastaReader(fasta) as fasta:
            for num_records, record in enumerate(fasta, 1):
                name = record.name.split(None, 1)[0]
                seq = record.sequence
                self.add(name, seq)
        if close:
            fasta.close()
        return num_records

    def load_default(self):
        """Tries to load from default URL first, then from default path.
        """
        try:
            return self.load_from_url()
        except (OSError, IOError):
            logging.getLogger().warning(
                "Error loading adapters from URL %s; loading from file",
                DEFAULT_ADAPTERS_URL)
        try:
            return self.load_from_file()
        except IOError:
            logging.getLogger().warning(
                "Error loading adapters from file %s; loading from file",
                DEFAULT_ADAPTERS_PATH)

    @property
    def names(self):
        """Sequence of adapter names.
        """
        return list(self.name_to_seq.keys())

    @property
    def sequences(self):
        """Sequence of adapter sequences.
        """
        return list(self.seq_to_name.keys())

    def iter_names(self):
        """Returns an iterator over adapter names.
        """
        return self.name_to_seq.items()

    def iter_sequences(self):
        """Returns an iterator over adapter sequences.
        """
        return self.seq_to_name.items()

    def has_name(self, name):
        """Returns whether this cache contains the specified name.

        Args:
            name: The adapter name.
        """
        return name in self.name_to_seq

    def get_for_name(self, name):
        """Returns the sequence associated with a name.

        Args:
            name: The name to fetch.

        Returns:
            The sequence.
        """
        return self.name_to_seq[name]

    def has_seq(self, seq):
        """Tests whether a sequence is in the cache.

        Args:
            seq: The sequence to check.

        Returns:
            True if the sequence is in the cache.
        """
        return seq in self.seq_to_name

    def get_for_seq(self, seq):
        """Returns the name associated with a given sequence.

        Args:
            seq: The sequence to fetch.

        Returns:
            The name associated with the sequence.
        """
        return list(self.seq_to_name[seq])

    def summarize(self):
        """Returns a summary dict. Does *not* add sequence info.
        """
        return dict(
            path=self.path,
            auto_reverse_complement=self.auto_reverse_complement,
            num_adapter_names=len(self.name_to_seq),
            num_adapter_seqs=len(self.seq_to_name))

def parse_braces(sequence):
    """Replace all occurrences of ``x{n}`` (where x is any character) with n
    occurrences of x. Raise ValueError if the expression cannot be parsed.

    Examples:
        >>> parse_braces('TGA{5}CT')
        TGAAAAACT
    """
    # Simple DFA with four states, encoded in prev
    result = ''
    prev = None
    for char in re.split(r'(\{|\})', sequence):
        if char == '':
            continue
        if prev is None:
            if char == '{':
                raise ValueError('"{" must be used after a character')
            if char == '}':
                raise ValueError('"}" cannot be used here')
            prev = char
            result += char
        elif prev == '{':
            prev = int(char)
            if not 0 <= prev <= 10000:
                raise ValueError('Value {} invalid'.format(prev))
        elif isinstance(prev, int):
            if char != '}':
                raise ValueError('"}" expected')
            result = result[:-1] + result[-1] * prev
            prev = None
        else:
            if char != '{':
                raise ValueError('Expected "{"')
            prev = '{'
    # Check if we are in a non-terminating state
    if isinstance(prev, int) or prev == '{':
        raise ValueError("Unterminated expression")
    return result

def _extract_name_from_spec(spec):
    """Parse an adapter specification given as 'name=adapt' into 'name'
    and 'adapt'.

    Args:
        spec: Adapter spec.

    Returns:
        (name, spec)
    """
    fields = spec.split('=', 1)
    if len(fields) > 1:
        name, spec = fields
        name = name.strip()
    else:
        name = None
    spec = spec.strip()
    return name, spec

ADAPTER_ID_GENERATOR = itertools.count(1)

def _generate_adapter_name():
    return str(next(ADAPTER_ID_GENERATOR))
