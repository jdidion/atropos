from collections import namedtuple
from enum import Enum
import itertools
from pathlib import Path
import pickle
import re
from typing import (
    Callable,
    Dict,
    Generic,
    IO,
    Iterable,
    Optional,
    Sequence as SequenceType,
    Set,
    Tuple,
    TypeVar,
    Union,
)
from urllib.error import URLError
from urllib.request import urlopen

from loguru import logger
from xphyle import open_
from xphyle.types import PathLike

from atropos.aligner import (
    Aligner, GapRule, MatchTuple, MultiAligner, compare_prefixes, compare_suffixes
)
from atropos.io.sequence import ColorspaceSequence, Sequence
from atropos.io.readers import FastaReader
from atropos.utils import colorspace
from atropos.utils.collections import Const, CountingDict, MergingDict, NestedDict
from atropos.utils.ngs import (
    ALPHABETS,
    GC_BASES,
    IUPAC_BASES,
    Alphabet,
    reverse_complement,
)
from atropos.utils.paths import as_readable_path
from atropos.utils.statistics import RandomMatchProbability


# TODO: specify this externally rather than hard-coding
DEFAULT_ADAPTERS_URL = "https://raw.githubusercontent.com/jdidion/atropos/master/atropos/adapters/sequencing_adapters.fa"
DEFAULT_ADAPTERS_PATH = Path(__file__).parent / "sequencing_adapters.fa"
DEFAULT_ADAPTER_CACHE_FILE = ".adapters"

DEFAULT_INSERT_MAX_RMP = 1e-6
"""Default value for InsertAligner max insert random match probability"""
DEFAULT_ADAPTER_MAX_RMP = 0.001
"""Default value for InsertAligner max adapter random match probability"""
DEFAULT_MIN_INSERT_OVERLAP = 1
"""Default value for InsertAligner min insert overlap"""
DEFAULT_MAX_INSERT_MISMATCH_FRAC = 0.2
"""Default value for InsertAligner max insert mismatch fraction"""
DEFAULT_MIN_ADAPTER_OVERLAP = 1
"""Default value for InsertAligner min adapter overlap"""
DEFAULT_MAX_ADAPTER_MISMATCH_FRAC = 0.2
"""Default value for InsertAligner max adapter mismatch fraction"""
DEFAULT_ADAPTER_CHECK_CUTOFF = 9
"""Default value for InsertAligner adapter check cutoff"""

M = TypeVar("M", bound="Match")
S = TypeVar("S", bound=Sequence)


class AdapterType(int, Enum):
    """
    Enuemerator of adapter types.
    """

    def __new__(cls, desc: str, flags: int):
        """
        Args:
            desc: Adapter type description.
            flags: Alignment flags.
        """
        obj = int.__new__(cls, flags)
        obj._value_ = flags
        obj.desc = desc
        return obj

    BACK = (
        "regular 3'",
        (
            GapRule.START_WITHIN_SEQ2
            | GapRule.STOP_WITHIN_SEQ2
            | GapRule.STOP_WITHIN_SEQ1
        ),
    )
    FRONT = (
        "regular 5'",
        (
            GapRule.START_WITHIN_SEQ2
            | GapRule.STOP_WITHIN_SEQ2
            | GapRule.START_WITHIN_SEQ1
        ),
    )
    PREFIX = ("anchored 5'", GapRule.STOP_WITHIN_SEQ2)
    SUFFIX = ("anchored 3'", GapRule.START_WITHIN_SEQ2)
    ANYWHERE = ("variable 5'/3'", GapRule.SEMIGLOBAL)
    LINKED = ("linked", 0)

    def as_dict(self) -> dict:
        """
        Returns AdapterType fields in a dict.
        """
        return dict(name=self.name, desc=self.desc, flags=Const(self.value))


ADAPTER_TYPE_NAMES = set(AdapterType.__members__.keys())


MatchInfo = namedtuple(
    "MatchInfo",
    (
        "read_name",
        "errors",
        "rstart",
        "rstop",
        "seq_before",
        "seq_adapter",
        "seq_after",
        "adapter_name",
        "qual_before",
        "qual_adapter",
        "qual_after",
        "is_front",
        "asize",
        "rsize_adapter",
        "rsize_total",
    ),
)


class Match:
    """
    Common match-result object returned by aligners.
    """

    __slots__ = [
        "astart",
        "astop",
        "rstart",
        "rstop",
        "matches",
        "errors",
        "front",
        "length",
    ]

    def __init__(
        self,
        astart: int,
        astop: int,
        rstart: int,
        rstop: int,
        matches: int,
        errors: int,
        front: Optional[bool] = None,
    ):
        """
        Args:
            astart: Starting position of the match within the adapter.
            astop: Ending position of the match within the adapter.
            rstart: Starting position of the match within the read.
            rstop: Ending position of the match within the read.
            matches: Number of matching bases.
            errors: Number of mismatching bases.
            front: Whether the match is to the front of the read.
        """
        self.astart = astart
        self.astop = astop
        self.rstart = rstart
        self.rstop = rstop
        self.matches = matches
        self.errors = errors
        self.front = self._guess_is_front() if front is None else front
        # Number of aligned characters in the adapter. If there are indels,
        # this may be different from the number of characters in the read.
        self.length = self.astop - self.astart
        if self.length <= 0:
            raise ValueError("Match length must be >= 0")
        if self.length - self.errors <= 0:
            raise ValueError("A Match requires at least one matching position.")
        # TODO: this assertion may not always hold now that we use both error
        #  rates and probabilities
        # if self.adapter:
        #    assert self.errors / self.length <= self.adapter.max_error_rate

    def __repr__(self) -> str:
        return (
            f"Match(astart={self.astart}, astop={self.astop}, "
            f"rstart={self.rstart}, rstop={self.rstop}, matches={self.matches}, "
            f"errors={self.errors})"
        )

    def copy(self: M) -> M:
        """
        Create a copy of this Match.
        """
        return Match(
            self.astart,
            self.astop,
            self.rstart,
            self.rstop,
            self.matches,
            self.errors,
            self.front,
        )

    def _guess_is_front(self) -> bool:
        """
        Guesses whether this is match is for a front adapter.

        The match is assumed to be a front adapter when the first base of
        the read is involved in the alignment to the adapter.
        """
        return self.rstart == 0


class AdapterMatch(Match, Generic[S]):
    __slots__ = ["_adapter", "_read"]

    def __init__(
        self,
        *args,
        adapter: Optional["Adapter"] = None,
        read: Optional[S] = None,
        **kwargs,
    ):
        """
        Args:
            adapter: The :class:`Adapter`.
            read: The :class:`Sequence`.
            *args: Positional arguments to pass to `Match` initializer.
            **kwargs: Keyword arguments to pass to `Match` initializer.
        """
        super().__init__(*args, **kwargs)
        self._adapter = adapter
        self._read = read

    def update(self, adapter: "Adapter", read: S, front: bool):
        self._adapter = adapter
        self._read = read
        self.front = front

    @property
    def adapter(self):
        return self._adapter

    @property
    def read(self):
        return self._read

    def copy(self: M) -> M:
        """
        Create a copy of this Match.
        """
        return AdapterMatch(
            self.astart,
            self.astop,
            self.rstart,
            self.rstop,
            self.matches,
            self.errors,
            self.front,
            adapter=self._adapter,
            read=self._read,
        )

    def wildcards(self, wildcard_char: str = "N") -> str:
        """
        Returns a string that contains, for each wildcard character,
        the character that it matches. For example, if the adapter
        ATNGNA matches ATCGTA, then the string 'CT' is returned.

        If there are indels, this is not reliable as the full alignment
        is not available.
        """
        return "".join(
            (
                self.read.sequence[self.rstart + i]
                for i in range(self.length)
                if (
                    self._adapter.sequence[self.astart + i] == wildcard_char
                    and self.rstart + i < len(self._read.sequence)
                )
            )
        )

    def rest(self) -> str:
        """
        Returns the part of the read before this match if this is a 'front' (5')
        adapter, or the part after the match if this is not a 'front' adapter (3').
        This can be an empty string.
        """
        if self.front:
            return self._read.sequence[: self.rstart]
        else:
            return self._read.sequence[self.rstop :]

    def get_info_record(self) -> MatchInfo:  # TODO: write test
        """
        Returns a :class:`MatchInfo`, which contains information about the match to
        write into the info file.
        """
        seq = self._read.sequence
        qualities = self._read.qualities

        if qualities is None:
            qualities = ""

        rsize = rsize_total = self.rstop - self.rstart

        if self.front and self.rstart > 0:
            rsize_total = self.rstop
        elif not self.front and self.rstop < len(seq):
            rsize_total = len(seq) - self.rstart

        return MatchInfo(
            self._read.name,
            self.errors,
            self.rstart,
            self.rstop,
            seq[0 : self.rstart],
            seq[self.rstart : self.rstop],
            seq[self.rstop :],
            self._adapter.name,
            qualities[0 : self.rstart],
            qualities[self.rstart : self.rstop],
            qualities[self.rstop :],
            self.front,
            self.astop - self.astart,
            rsize,
            rsize_total,
        )


InsertMatchResult = Tuple[MatchTuple, Optional[AdapterMatch], Optional[AdapterMatch]]


class InsertAligner:
    """
    Implementation of an insert matching algorithm.

    If the inserts align, the overhangs are searched for the adapter sequences.
    Otherwise, each read is search for its adapter separately.

    This only works with paired-end reads with 3' adapters.
    """

    def __init__(
        self,
        adapter1: str,
        adapter2: str,
        match_probability: Optional[RandomMatchProbability] = None,
        insert_max_rmp: float = DEFAULT_INSERT_MAX_RMP,
        adapter_max_rmp: float = DEFAULT_ADAPTER_MAX_RMP,
        min_insert_overlap: int = DEFAULT_MIN_INSERT_OVERLAP,
        max_insert_mismatch_frac: Union[str, float] = DEFAULT_MAX_INSERT_MISMATCH_FRAC,
        min_adapter_overlap: int = DEFAULT_MIN_ADAPTER_OVERLAP,
        max_adapter_mismatch_frac: Union[str, float] =
        DEFAULT_MAX_ADAPTER_MISMATCH_FRAC,
        adapter_check_cutoff: int = DEFAULT_ADAPTER_CHECK_CUTOFF,
        base_probs: Optional[Dict[str, float]] = None,
        adapter_wildcards: bool = True,
        read_wildcards: bool = False,
    ):
        """
        Args:
            adapter1: read1 adapter.
            adapter2: read2 adapter.
            match_probability: Callable that calculates random match probability
                given arguments (num_matches, match_len).
            insert_max_rmp: Max random match probability for the insert match.
            adapter_max_rmp: Max random match probability for the adapter match.
            min_insert_overlap: Minimum number of bases the inserts must overlap
                to be considered matching.
            max_insert_mismatch_frac: Maximum fraction of mismatching bases between
                the inserts to be considered matching.
            min_adapter_overlap: Minimum number of bases the adapter must overlap
                the read to be considered matching.
            max_adapter_mismatch_frac: Maximum fraction of mismatching bases between
                the adapter and the read to be considered matching.
            adapter_check_cutoff: Threshold number of matching bases required before
                adapter matches are checked against random match probability.
            base_probs: Dict of (match_prob, mismatch_prob), which are the
                probabilities passed to
                :method:`atropos.util.RandomMatchProbability.__call__()`.
        """

        self.adapter1 = adapter1
        self.adapter1_len = len(adapter1)
        self.adapter2 = adapter2
        self.adapter2_len = len(adapter2)
        self.match_probability = match_probability or RandomMatchProbability()
        self.insert_max_rmp = insert_max_rmp
        self.adapter_max_rmp = adapter_max_rmp
        self.min_insert_overlap = min_insert_overlap
        self.max_insert_mismatch_frac = float(max_insert_mismatch_frac)
        self.min_adapter_overlap = min_adapter_overlap
        self.max_adapter_mismatch_frac = float(max_adapter_mismatch_frac)
        self.adapter_check_cutoff = adapter_check_cutoff
        self.base_probs = base_probs or dict(match_prob=0.25, mismatch_prob=0.75)
        self.adapter_wildcards = adapter_wildcards
        self.read_wildcards = read_wildcards
        self.aligner = MultiAligner(
            max_insert_mismatch_frac,
            GapRule.START_WITHIN_SEQ1 | GapRule.STOP_WITHIN_SEQ2,
            min_insert_overlap,
        )

    def __repr__(self) -> str:
        return f"InsertAligner<adapter1={self.adapter1}, adapter2={self.adapter2}, " \
               f"match_probability={self.match_probability}, " \
               f"insert_max_rmp={self.insert_max_rmp}, " \
               f"adapter_max_rmp={self.adapter_max_rmp}, " \
               f"min_insert_overlap={self.min_insert_overlap}, " \
               f"max_insert_mismatch_frac={self.max_insert_mismatch_frac}, " \
               f"min_adapter_overlap={self.min_adapter_overlap}, " \
               f"max_adapter_mismatch_frac={self.max_adapter_mismatch_frac}, " \
               f"adapter_check_cutoff={self.adapter_check_cutoff}, " \
               f"base_probs={self.base_probs}, " \
               f"adapter_wildcards={self.adapter_wildcards}, " \
               f"read_wildcards={self.read_wildcards}, " \
               f"aligner={self.aligner}>"

    def match_insert(self, seq1: str, seq2: str) -> Optional[InsertMatchResult]:
        """
        Uses cutadapt aligner for insert and adapter matching.

        Args:
            seq1, seq2: Sequences to match.

        Returns:
            A :class:`Match` object, or None if there is no match.
        """
        seq_len1 = len(seq1)
        seq_len2 = len(seq2)

        if seq_len1 > seq_len2:
            seq1 = seq1[:seq_len2]
            seq_len = seq_len2
        else:
            seq_len = seq_len1
            if seq_len2 > seq_len1:
                seq2 = seq2[:seq_len1]

        seq2_rc = reverse_complement(seq2)

        def _match(
            _insert_match: MatchTuple, _offset: int, _insert_match_size: int
        ) -> Optional[InsertMatchResult]:
            if _offset < self.min_adapter_overlap:
                # The reads are mostly overlapping, to the point where
                # there's not enough overhang to do a confident adapter
                # match. We return just the insert match to signal that
                # error correction can be done even though no adapter
                # trimming is required.
                return _insert_match, None, None

            # TODO: this is very sensitive to the exact correct choice of adapter.
            #  For example, if you specifiy GATCGGAA... and the correct adapter is
            #  AGATCGGAA..., the prefixes will not match exactly and the alignment
            #  will fail. We need to use a comparison that is a bit more forgiving.
            def _match_adapter(
                insert_seq: str, adapter_seq: str, adapter_len: int
            ) -> Tuple[MatchTuple, int, float]:
                amatch = compare_prefixes(
                    insert_seq[_insert_match_size:],
                    adapter_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards
                )
                alen = min(_offset, adapter_len)
                return amatch, alen, round(alen * self.max_adapter_mismatch_frac)

            a1_match, a1_length, a1_max_mismatches = _match_adapter(
                seq1, self.adapter1, self.adapter1_len
            )
            a2_match, a2_length, a2_max_mismatches = _match_adapter(
                seq2, self.adapter2, self.adapter2_len
            )

            if (
                a1_match[5] > a1_max_mismatches and
                a2_match[5] > a2_max_mismatches
            ):
                return None

            if min(a1_length, a2_length) > self.adapter_check_cutoff:
                a1_prob = self.match_probability(a1_match[4], a1_length)
                a2_prob = self.match_probability(a2_match[4], a2_length)
                if (a1_prob * a2_prob) > self.adapter_max_rmp:
                    return None

            mismatches = min(a1_match[5], a2_match[5])

            def _create_adapter_match(alen: int, slen: int) -> AdapterMatch:
                alen = min(alen, slen - _insert_match_size)
                _mismatches = min(alen, mismatches)
                _matches = alen - _mismatches
                return AdapterMatch(
                    0, alen, _insert_match_size, slen, _matches, _mismatches
                )

            return (
                _insert_match,
                _create_adapter_match(a1_length, seq_len1),
                _create_adapter_match(a2_length, seq_len2)
            )

        insert_matches = self.aligner.locate(seq2_rc, seq1)

        if insert_matches:
            # Filter by random-match probability
            filtered_matches = []

            for insert_match in insert_matches:
                offset = min(insert_match[0], seq_len - insert_match[3])
                insert_match_size = seq_len - offset
                prob = self.match_probability(
                    insert_match[4], insert_match_size, **self.base_probs
                )
                if prob <= self.insert_max_rmp:
                    filtered_matches.append(
                        (insert_match, offset, insert_match_size, prob)
                    )

            if filtered_matches:
                if len(filtered_matches) == 1:
                    return _match(*filtered_matches[0][0:3])
                else:
                    # Test matches in order of random-match probability.
                    # TODO: compare against sorting by length (which is how
                    #  SeqPurge essentially does it).
                    #  filtered_matches.sort(key=lambda x: x[2], reverse=True)
                    filtered_matches.sort(key=lambda x: x[3])
                    for match_args in filtered_matches:
                        match = _match(*match_args[0:3])
                        if match:
                            return match


class AdapterCache:
    """
    Cache for known adapters.
    """

    def __init__(
        self,
        path: Optional[Union[str, PathLike]] = Path(DEFAULT_ADAPTER_CACHE_FILE),
        auto_reverse_complement: bool = False,
    ):
        """
        Args:
            path: Path to file where cache is written.
            auto_reverse_complement: Whether adapter reverse-complements should
                automatically be added to the cache.
        """
        self.path = as_readable_path(path, False)
        self.auto_reverse_complement = auto_reverse_complement

        if self.path and path.exists():
            with open_(self.path, "rb") as cache:
                self.seq_to_name, self.name_to_seq = pickle.load(cache)
        else:
            self.seq_to_name = {}
            self.name_to_seq = {}

    @property
    def empty(self) -> bool:
        """
        Whether the cache is empty.
        """
        return len(self.seq_to_name) == 0

    def save(self):
        """
        Save the cache to disk.
        """
        if self.path:
            with open_(self.path, "wb") as cache:
                pickle.dump((self.seq_to_name, self.name_to_seq), cache)

    def add(self, name: str, seq: str):
        """
        Adds a sequence to the cache.

        Args:
            name: Adapter name.
            seq: Adapter sequence.
        """
        self._add(name, seq)
        if self.auto_reverse_complement:
            self._add(f"{name}_rc", reverse_complement(seq))

    def _add(self, name: str, seq: str):
        if seq not in self.seq_to_name:
            self.seq_to_name[seq] = set()
        self.seq_to_name[seq].add(name)
        self.name_to_seq[name] = seq

    def load_from_file(self, path: Union[str, Path] = DEFAULT_ADAPTERS_PATH) -> int:
        """
        Loads cached data from a file.

        Args:
            path: Path from which to load.
        """
        p = as_readable_path(path, True)
        with open_(p, "rt") as inp:
            return self.load_from_fasta(inp)

    def load_from_url(self, url: str = DEFAULT_ADAPTERS_URL) -> int:
        """
        Loads adapter data from a URL.

        Args:
            url: URL from which to load.
        """
        logger.info(f"Loading list of known contaminants from {url}")

        try:
            fasta = urlopen(url).read().decode().split("\n")
            return self.load_from_fasta(fasta)
        except URLError:
            if url.startswith("file:"):
                url = url[5:]
            return self.load_from_file(url)

    def load_from_fasta(self, fasta: Union[str, IO]) -> int:
        """
        Loads adapter data from a FASTA file.

        Args:
            fasta: FASTA file.
        """
        with FastaReader(fasta) as fasta:
            for num_records, record in enumerate(fasta, 1):
                name = record.name.split(None, 1)[0]
                seq = record.sequence
                self.add(name, seq)

        return num_records

    def load_default(self) -> int:
        """
        Tries to load from default URL first, then from default path.
        """
        try:
            return self.load_from_url()
        except (OSError, IOError):
            logger.opt(exception=True).warning(
                f"Error loading adapters from URL {DEFAULT_ADAPTERS_URL}; loading from"
                f"file instead"
            )

        try:
            return self.load_from_file()
        except IOError:
            logger.opt(exception=True).warning(
                f"Error loading adapters from file {DEFAULT_ADAPTERS_PATH}"
            )

    @property
    def names(self) -> SequenceType[str]:
        """
        Sequence of adapter names.
        """
        return tuple(self.name_to_seq.keys())

    @property
    def sequences(self) -> SequenceType[str]:
        """
        Sequence of adapter sequences.
        """
        return tuple(self.seq_to_name.keys())

    @property
    def sequence_name_pairs(self) -> SequenceType[Tuple[str, Set[str]]]:
        """
        Sequence of (sequence, Set[name]) tuples.
        """
        return tuple(self.seq_to_name.items())

    def has_name(self, name: str) -> bool:
        """
        Returns whether this cache contains the specified name.

        Args:
            name: The adapter name.
        """
        return name in self.name_to_seq

    def get_for_name(self, name: str) -> str:
        """
        Gets the sequence associated with a name.

        Args:
            name: The name to fetch.

        Returns:
            The sequence.
        """
        return self.name_to_seq[name]

    def has_seq(self, seq: str) -> bool:
        """
        Tests whether a sequence is in the cache.

        Args:
            seq: The sequence to check.

        Returns:
            True if the sequence is in the cache.
        """
        return seq in self.seq_to_name

    def get_for_seq(self, seq: str) -> Set[str]:
        """
        Gets the name associated with a given sequence.

        Args:
            seq: The sequence to fetch.

        Returns:
            The name associated with the sequence.
        """
        return self.seq_to_name[seq]

    def summarize(self) -> dict:
        """
        Returns a summary dict. Does *not* add sequence info.
        """
        return dict(
            path=self.path,
            auto_reverse_complement=self.auto_reverse_complement,
            num_adapter_names=len(self.name_to_seq),
            num_adapter_seqs=len(self.seq_to_name),
        )


class AdapterBase:
    ADAPTER_ID_GENERATOR = itertools.count(1)

    @staticmethod
    def _generate_adapter_name():
        return str(next(Adapter.ADAPTER_ID_GENERATOR))


class Adapter(AdapterBase):
    """
    An adapter knows how to match itself to a read. In particular, it knows where it
    should be within the read and how to interpret wildcard characters.
    """

    def __init__(
        self,
        sequence: str,
        where: int,
        max_error_rate: float = 0.1,
        min_overlap: int = 3,
        read_wildcards: bool = False,
        adapter_wildcards: bool = True,
        name: Optional[str] = None,
        indels: bool = True,
        indel_cost: int = 1,
        match_probability: Optional[RandomMatchProbability] = None,
        max_rmp: Optional[float] = None,
        gc_content: float = 0.5,
        alphabet: Optional[Alphabet] = None,
    ):
        """
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
        if len(sequence) == 0:
            raise ValueError("Empty adapter sequence")

        # TODO: all of this validation code should be wrapped up in Alphabet
        sequence = parse_braces(sequence.upper().replace("U", "T"))
        seq_set = set(sequence)

        if seq_set <= set("ACGT"):
            adapter_wildcards = False

        if adapter_wildcards and not seq_set <= IUPAC_BASES:
            raise ValueError(
                f"Invalid character(s) in adapter sequence: "
                f"{','.join(seq_set - IUPAC_BASES)}"
            )

        if alphabet:
            if isinstance(alphabet, str):
                alphabet = ALPHABETS[alphabet]
            alphabet.validate_string(sequence)

        self.debug = False
        self.name = name or AdapterBase._generate_adapter_name()
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
        trimmers: Dict[int, Callable] = {
            AdapterType.FRONT: self._trimmed_front,
            AdapterType.PREFIX: self._trimmed_front,
            AdapterType.BACK: self._trimmed_back,
            AdapterType.SUFFIX: self._trimmed_back,
            AdapterType.ANYWHERE: self._trimmed_anywhere,
        }
        self.trimmed = trimmers[where]

        if where == AdapterType.ANYWHERE:
            self._front_flag = None  # means: guess
        else:
            self._front_flag = where not in (AdapterType.BACK, AdapterType.SUFFIX)

        self.aligner = Aligner(
            reference=self.sequence,
            max_error_rate=self.max_error_rate,
            flags=int(self.where),
            wildcard_ref=self.adapter_wildcards,
            wildcard_query=self.read_wildcards,
        )
        self.aligner.min_overlap = self.min_overlap
        if self.indels:
            self.aligner.indel_cost = indel_cost
        else:
            # TODO: when indels are disallowed, an entirely different algorithm
            #  should be used.
            self.aligner.indel_cost = 100000

        # statistics about length of removed sequences
        self.lengths_front = CountingDict()
        self.lengths_back = CountingDict()
        self.errors_front = NestedDict()
        self.errors_back = NestedDict()
        self.adjacent_bases = {"A": 0, "C": 0, "G": 0, "T": 0, "": 0}

    def __len__(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return (
            '<Adapter(name="{name}", sequence="{sequence}", where={where}, '
            "max_error_rate={max_error_rate}, min_overlap={min_overlap}, "
            "read_wildcards={read_wildcards}, "
            "adapter_wildcards={adapter_wildcards}, indels={indels})>".format(
                **vars(self)
            )
        )

    def enable_debug(self):
        """
        Enables logging of the dynamic programming matrix after matching a read to an
        adapter.
        """
        self.debug = True
        self.aligner.enable_debug()

    def match_to(self, read: Sequence) -> Optional[AdapterMatch]:
        """
        Attempts to match this adapter to the given read.

        Args:
            read: A :class:`Sequence` instance.

        Returns:
            A :class:`Match` instance if a match was found; return None if no match
            was found given the matching criteria (minimum overlap length, maximum
            error rate).
        """
        read_seq = read.sequence.upper()

        # try to find an exact match first unless wildcards are allowed
        pos = -1

        if not self.adapter_wildcards:
            if self.where is AdapterType.PREFIX:
                if read_seq.startswith(self.sequence):
                    pos = 0
            elif self.where is AdapterType.SUFFIX:
                if read_seq.endswith(self.sequence):
                    pos = len(read_seq) - len(self.sequence)
            else:
                pos = read_seq.find(self.sequence)

        if pos >= 0:
            seqlen = len(self.sequence)
            return AdapterMatch(
                0,
                seqlen,
                pos,
                pos + seqlen,
                seqlen,
                0,
                self._front_flag,
                adapter=self,
                read=read,
            )

        # try approximate matching
        if not self.indels and self.where in (AdapterType.PREFIX, AdapterType.SUFFIX):
            if self.where == AdapterType.PREFIX:
                alignment = compare_prefixes(
                    self.sequence,
                    read_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards,
                )
            else:
                alignment = compare_suffixes(
                    self.sequence,
                    read_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards,
                )
        else:
            alignment = self.aligner.locate(read_seq)
            if self.debug:
                print(self.aligner.dpmatrix)  # pragma: no cover

        if alignment:
            astart, astop, rstart, rstop, matches, errors = alignment
            size = astop - astart
            if (size >= self.min_overlap and errors / size <= self.max_error_rate) and (
                self.max_rmp is None
                or self.match_probability(matches, size) <= self.max_rmp
            ):
                return AdapterMatch(
                    astart,
                    astop,
                    rstart,
                    rstop,
                    matches,
                    errors,
                    self._front_flag,
                    adapter=self,
                    read=read,
                )

        return None

    def _trimmed_anywhere(self, match: AdapterMatch) -> Sequence:
        """
        Trims an adapter from either the front or back of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        if match.front:
            return self._trimmed_front(match)
        else:
            return self._trimmed_back(match)

    def _trimmed_front(self, match: AdapterMatch) -> Sequence:
        """
        Trims an adapter from the front of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        self.lengths_front[match.rstop] += 1
        self.errors_front[match.rstop][match.errors] += 1
        return match.read[match.rstop :]

    def _trimmed_back(self, match: AdapterMatch) -> Sequence:
        """
        Trims an adapter from the back of sequence.

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        self.lengths_back[len(match.read) - match.rstart] += 1
        self.errors_back[len(match.read) - match.rstart][match.errors] += 1
        adjacent_base = match.read.sequence[match.rstart - 1 : match.rstart]
        if adjacent_base not in "ACGT":
            adjacent_base = ""
        self.adjacent_bases[adjacent_base] += 1
        return match.read[: match.rstart]

    def random_match_probabilities(self) -> SequenceType[float]:
        """
        Estimates the probabilities that this adapter matches a random sequence. Indels
        are not taken into account.

        Returns:
            A list of probabilities the same length as this adapter's sequence,
            where the value at position 'i' is the probability that i bases of this
            adapter match a random sequence.
        """
        if self._front_flag:
            seq = self.sequence[::-1]
        else:
            seq = self.sequence

        base_probs = (self.gc_content / 2.0, (1 - self.gc_content) / 2.0)
        probabilities = [1.0] + ([0] * len(seq))
        c_bases = frozenset(GC_BASES if self.adapter_wildcards else "GC")

        # TODO: this doesn't work; need to figure out if RandomMatchProbability
        # can be used for this.
        # matches = 0
        # for idx, base in enumerate(seq, 1):
        #     if base in c_bases:
        #         matches += 1
        #     probabilities[idx] = self.match_probability(
        #         matches, idx, *base_probs)

        cur_p = 1.0
        for idx, base in enumerate(seq, 1):
            cur_p *= base_probs[0 if base in c_bases else 1]
            probabilities[idx] = cur_p

        return probabilities

    def summarize(self) -> dict:
        """
        Summarizes the activities of this :class:`Adapter`.
        """
        where = self.where
        total_front = sum(self.lengths_front.values())
        total_back = sum(self.lengths_back.values())

        if not (
            where in (AdapterType.ANYWHERE, AdapterType.LINKED)
            or (where in (AdapterType.BACK, AdapterType.SUFFIX) and total_front == 0)
            or (where in (AdapterType.FRONT, AdapterType.PREFIX) and total_back == 0)
        ):
            raise ValueError(
                f"Invalid combination of where={where}, total_front={total_front}, "
                f"total_back={total_back}"
            )

        metrics = MergingDict(
            adapter_class=self.__class__.__name__,
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back,
            match_probabilities=Const(self.random_match_probabilities()),
        )
        metrics["where"] = AdapterType(where).as_dict()
        metrics["sequence"] = Const(self.sequence)
        metrics["max_error_rate"] = Const(self.max_error_rate)
        if where in (AdapterType.ANYWHERE, AdapterType.FRONT, AdapterType.PREFIX):
            metrics["lengths_front"] = self.lengths_front
            metrics["errors_front"] = self.errors_front
        if where in (AdapterType.ANYWHERE, AdapterType.BACK, AdapterType.SUFFIX):
            metrics["lengths_back"] = self.lengths_back
            metrics["errors_back"] = self.errors_back
        if where in (AdapterType.BACK, AdapterType.SUFFIX):
            metrics["adjacent_bases"] = self.adjacent_bases

        return metrics


class ColorspaceAdapter(Adapter):
    """
    An adapter for a colorspace sequence.
    """

    def __init__(self, *args, **kwargs):
        """
        Args:
            args: Positional arguments passed to :class:`Adapter` constructor.
            kwargs: Keyword arguments passed to :class:`Adapter` constructor.
        """
        if kwargs.get("adapter_wildcards"):
            raise ValueError("Wildcards not supported for colorspace adapters")
        else:
            kwargs["adapter_wildcards"] = False

        super().__init__(*args, **kwargs)

        if set(self.sequence) <= set("ACGT"):
            # adapter was given in basespace
            self.nucleotide_sequence = self.sequence
            self.sequence = colorspace.encode(self.sequence)[1:]
        elif self.where in (AdapterType.PREFIX, AdapterType.FRONT):
            raise ValueError(
                "A 5' colorspace adapter needs to be given in nucleotide space"
            )

        self.aligner.reference = self.sequence

    def match_to(self, read: ColorspaceSequence) -> Optional[AdapterMatch]:
        """
        Attempts to match this adapter to the given read.

        Args:
            read: A :class:`Sequence` instance.

        Returns:
            A :class:`Match` instance if a match was found; return None if no
            match was found given the matching criteria (minimum overlap length,
            maximum error rate).
        """
        if self.where != AdapterType.PREFIX:
            return super().match_to(read)

        # create artificial adapter that includes a first color that encodes the
        # transition from primer base into adapter
        asequence = (
            colorspace.ENCODE[read.primer + self.nucleotide_sequence[0:1]]
            + self.sequence
        )

        if read.sequence.startswith(asequence):
            match = AdapterMatch(
                0,
                len(asequence),
                0,
                len(asequence),
                len(asequence),
                0,
                self._front_flag,
                adapter=self,
                read=read,
            )
        else:
            # try approximate matching
            self.aligner.reference = asequence
            alignment = self.aligner.locate(read.sequence)

            if self.debug:
                logger.debug(self.aligner.dpmatrix)

            if alignment:
                match = AdapterMatch(
                    *alignment, self._front_flag, adapter=self, read=read
                )
            else:
                return None

        # The following assertions should be guarnateed by the matching algorithm
        assert match.length > 0
        assert (match.errors / match.length) <= self.max_error_rate
        assert match.length >= self.min_overlap

        return match

    def _trimmed_front(
        self, match: AdapterMatch[ColorspaceSequence]
    ) -> ColorspaceSequence:
        """
        Trims an adapter from the front of sequence.

        Args:
            match:

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        read = match.read
        self.lengths_front[match.rstop] += 1
        self.errors_front[match.rstop][match.errors] += 1

        # to remove a front adapter, we need to re-encode the first color following
        # the adapter match
        color_after_adapter = read.sequence[match.rstop : match.rstop + 1]
        if not color_after_adapter:
            # the read is empty
            return read[match.rstop :]

        base_after_adapter = colorspace.DECODE[
            self.nucleotide_sequence[-1:] + color_after_adapter
        ]
        new_first_color = colorspace.ENCODE[read.primer + base_after_adapter]

        new_read = read[:]
        new_read.sequence = new_first_color + read.sequence[(match.rstop + 1) :]
        new_read.qualities = None
        if read.qualities:
            new_read.qualities = read.qualities[match.rstop :]

        return new_read

    def _trimmed_back(self, match: AdapterMatch) -> ColorspaceSequence:
        """
        Trims an adapter from the back of sequence.

        Args:
            match:

        Returns:
            A :class:`Sequence` instance: the trimmed read.
        """
        # trim one more color if long enough
        adjusted_rstart = max(match.rstart - 1, 0)
        self.lengths_back[len(match.read) - adjusted_rstart] += 1
        self.errors_back[len(match.read) - adjusted_rstart][match.errors] += 1
        return match.read[:adjusted_rstart]

    def __repr__(self) -> str:
        return f"<ColorspaceAdapter(sequence={self.sequence!r}, where={self.where})>"


class LinkedMatch:
    """
    Represents a match of a LinkedAdapter.

    Todo:
        Consolidate AdapterMatch and LinkedMatch classes.
    """

    def __init__(
        self,
        front_match: AdapterMatch,
        back_match: AdapterMatch,
        adapter: "LinkedAdapter",
    ):
        """
        Args:
            front_match: The match to the front of the sequence.
            back_match: The match to the back of the sequence.
            adapter: The matched adapter.
        """
        if front_match is None:
            raise ValueError("front_match cannot be None")

        self.front_match = front_match
        self.back_match = back_match
        self.adapter = adapter

    def get_info_record(self) -> MatchInfo:
        """
        Returns the info record for the either the back or forward match.
        """
        if self.back_match:
            return self.back_match.get_info_record()
        else:
            return self.front_match.get_info_record()


class LinkedAdapter(AdapterBase):
    """
    An adapter with linked front and back sequences.
    """

    def __init__(
        self,
        front_sequence: str,
        back_sequence: str,
        front_anchored: bool = True,
        back_anchored: bool = False,
        name: str = None,
        **kwargs,
    ):
        """
        Args:
            front_sequence: Front adapter sequence.
            back_sequence: Back adapter sequence.
            front_anchored: Whether the front adapter is anchored.
            back_anchored: Whether the back adapter is anchored.
            name: Adapter name.
            kwargs: passed on to individual :class:`Adapter` constructors.
        """
        if front_anchored:
            self.front_anchored = True
        else:
            raise ValueError("Front adapter must be anchored")

        if back_anchored:
            raise ValueError("Back adapater must not be anchored")
        else:
            self.back_anchored = False

        self.where = AdapterType.LINKED
        self.name = AdapterBase._generate_adapter_name() if name is None else name
        self.front_adapter = Adapter(
            front_sequence, where=AdapterType.PREFIX, name=None, **kwargs
        )
        self.back_adapter = Adapter(
            back_sequence, where=AdapterType.BACK, name=None, **kwargs
        )

    def enable_debug(self):
        """
        Enables debug on both adapters.
        """
        self.front_adapter.enable_debug()
        self.back_adapter.enable_debug()

    def match_to(self, read: Sequence) -> Optional[LinkedMatch]:
        """
        Matches the linked adapters against the given read. If the 'front' adapter is
        not found, the 'back' adapter is not searched for.

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

        # TODO use match.trimmed() instead as soon as that does not update statistics
        #  anymore
        read = read[front_match.rstop :]
        back_match = self.back_adapter.match_to(read)
        return LinkedMatch(front_match, back_match, self)

    def trimmed(self, match: LinkedMatch) -> Sequence:
        """
        Returns the read trimmed with the front and/or back adapter trimmer(s).

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

    def summarize(self) -> dict:
        """
        Returns the summary dict for this adapter.
        """
        where = self.where
        total_front = sum(self.front_adapter.lengths_front.values())
        total_back = sum(self.back_adapter.lengths_back.values())
        if not (
            where in (AdapterType.ANYWHERE, AdapterType.LINKED)
            or (where in (AdapterType.BACK, AdapterType.SUFFIX) and total_front == 0)
            or (where in (AdapterType.FRONT, AdapterType.PREFIX) and total_back == 0)
        ):
            raise ValueError(
                f"Invalid combination of where={where}, total_front={total_front}, "
                f"total_back={total_back}"
            )

        metrics = MergingDict(
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back,
        )
        metrics["where"] = AdapterType(where).as_dict()
        metrics["front_sequence"] = Const(self.front_adapter.sequence)
        metrics["front_match_probabilities"] = Const(
            self.front_adapter.random_match_probabilities()
        )
        metrics["back_sequence"] = Const(self.back_adapter.sequence)
        metrics["back_match_probabilities"] = Const(
            self.back_adapter.random_match_probabilities()
        )
        metrics["front_max_error_rate"] = Const(self.front_adapter.max_error_rate)
        metrics["back_max_error_rate"] = Const(self.back_adapter.max_error_rate)
        metrics["front_lengths_front"] = self.front_adapter.lengths_front
        metrics["front_lengths_back"] = self.front_adapter.lengths_back
        metrics["back_lengths_front"] = self.back_adapter.lengths_front
        metrics["back_lengths_back"] = self.back_adapter.lengths_back
        # have to clone these nested dicts and set them up with a custom merge function
        metrics["front_errors_front"] = self.front_adapter.errors_front
        metrics["front_errors_back"] = self.front_adapter.errors_back
        metrics["back_errors_front"] = self.back_adapter.errors_front
        metrics["back_errors_back"] = self.back_adapter.errors_back

        return metrics


class AdapterParser:
    """
    Factory for Adapter classes that all use the same parameters (error rate,
    indels etc.). The given **kwargs will be passed to the Adapter constructors.

    Todo: Use an enum for cmdline_type
    """

    def __init__(
        self, colorspace: bool = False, cache: Optional[AdapterCache] = None, **kwargs
    ):
        """
        Args:
            colorspace: Whether the adapter sequences are in colorspace.
            cache: The AdapterCache to use for fetching/storing named adapters.
        """
        self.colorspace = colorspace
        self.cache = cache
        self.constructor_args = kwargs
        self.adapter_class = ColorspaceAdapter if colorspace else Adapter

    def parse(self, spec: str, cmdline_type: str = "back") -> Iterable[Adapter]:
        """
        Parse an adapter specification not using ``file:`` notation and return an
        object of an appropriate Adapter class. The notation for anchored 5' and 3'
        adapters is supported. If the name parameter is None, then an attempt is made
        to extract the name from the specification (If spec is 'name=ADAPTER',
        name will be 'name'.)

        Args:
            spec: The adapter spec.
            cmdline_type: describes which commandline parameter was used (``-a``
                is 'back', ``-b`` is 'anywhere', and ``-g`` is 'front').

        Todo: describe the adapter spec format

        Yields:
            :class:`Adapter` instances.
        """
        if spec.startswith("file:"):
            # read adapter sequences from a file
            with FastaReader(spec[5:]) as fasta:
                for record in fasta:
                    name = record.name.split(None, 1)[0]
                    yield self.parse_from_spec(record.sequence, cmdline_type, name)
        else:
            yield self.parse_from_spec(spec, cmdline_type)

    def parse_from_spec(
        self, spec: str, cmdline_type: str = "back", name: Optional[str] = None
    ) -> Union[Adapter, LinkedAdapter]:
        if name is None and spec is None:
            raise ValueError("Either name or spec must be given")

        type_name = cmdline_type.upper()

        if type_name not in ADAPTER_TYPE_NAMES:
            raise ValueError(f"cmdline_type cannot be {cmdline_type!r}")

        orig_spec = spec
        where = AdapterType[type_name]

        if name is None:
            if spec is None:
                raise ValueError("Either name or spec must be given")

            if self.cache and self.cache.has_name(spec):
                name = spec
                spec = self.cache.get_for_name(name)

        if spec is None:
            if self.cache and self.cache.has_name(name):
                spec = self.cache.get_for_name(name)
            else:
                raise ValueError(f"Name not found: {name}")

        if name is None:
            name, spec = _extract_name_from_spec(spec)

        if self.cache and name is not None:
            self.cache.add(name, spec)

        front_anchored, back_anchored = False, False

        if spec.startswith("^"):
            spec = spec[1:]
            front_anchored = True

        if spec.endswith("$"):
            spec = spec[:-1]
            back_anchored = True

        sequence1, middle, sequence2 = spec.partition("...")

        if where is AdapterType.ANYWHERE:
            if front_anchored or back_anchored:
                raise ValueError("'anywhere' (-b) adapters may not be anchored")

            if middle == "...":
                raise ValueError("'anywhere' (-b) adapters may not be linked")

            return self.adapter_class(
                sequence=spec, where=where, name=name, **self.constructor_args
            )

        if where not in (AdapterType.FRONT, AdapterType.BACK):
            raise ValueError(f"Invalid where: {where}")

        if middle == "...":
            if not sequence1:
                if where is AdapterType.BACK:  # -a ...ADAPTER
                    spec = sequence2
                else:  # -g ...ADAPTER
                    raise ValueError(f"Invalid adapter specification: {orig_spec}")
            elif not sequence2:
                if where is AdapterType.BACK:  # -a ADAPTER...
                    spec = sequence1
                    where = AdapterType.FRONT
                    front_anchored = True
                else:  # -g ADAPTER...
                    spec = sequence1
            else:
                # linked adapter
                if self.colorspace:
                    raise NotImplementedError(
                        "Using linked adapters in colorspace is not supported"
                    )

                # automatically anchor 5' adapter if -a is used
                if where is AdapterType.BACK:
                    front_anchored = True

                return LinkedAdapter(
                    sequence1,
                    sequence2,
                    name=name,
                    front_anchored=front_anchored,
                    back_anchored=back_anchored,
                    **self.constructor_args,
                )

        if front_anchored and back_anchored:
            raise ValueError(
                f'Trying to use both "^" and "$" in adapter specification '
                f"{orig_spec!r}"
            )

        if front_anchored:
            if where is AdapterType.BACK:
                raise ValueError("Cannot anchor the 3' adapter at its 5' end")

            where = AdapterType.PREFIX
        elif back_anchored:
            if where is AdapterType.FRONT:
                raise ValueError("Cannot anchor 5' adapter at 3' end")

            where = AdapterType.SUFFIX

        return self.adapter_class(
            sequence=spec, where=where, name=name, **self.constructor_args
        )

    def parse_multi(
        self,
        back: Optional[Iterable[str]] = None,
        anywhere: Optional[Iterable[str]] = None,
        front: Optional[Iterable[str]] = None,
    ) -> SequenceType[Adapter]:
        """
        Parses all three types of commandline options that can be used to specify
        adapters. `back`, `anywhere` and `front` are lists of strings, corresponding
        to the respective commandline types (-a, -b, -g).

        Args:
            back: Back-adapter specs.
            anywhere: Anywhere-adapter specs.
            front: Front-adapter specs.

        Return:
            A list of appropriate Adapter classes.
        """
        return list(
            itertools.chain.from_iterable(
                self.parse(spec, cmdline_type)
                for specs, cmdline_type in (
                    (back, "back"),
                    (anywhere, "anywhere"),
                    (front, "front"),
                )
                for spec in specs or ()
            )
        )


BRACES = re.compile(r"([{}])")
MAX_REPEATS = 10000


def parse_braces(sequence: str) -> str:
    """
    Replaces all occurrences of ``x{n}`` (where x is any character) with n occurrences
    of x.

    Returns:
        The string with brace expressions replaced.

    Raises:
        ValueError if the expression cannot be parsed.

    Examples:
        >>> parse_braces('TGA{5}CT')
        TGAAAAACT
    """
    result = ""
    prev = None

    for char in BRACES.split(sequence):
        if char == "":
            continue

        if prev is None:
            if char == "{":
                raise ValueError(
                    f"Invalid expression: {sequence};"
                    '"{" must be used after a character'
                )

            if char == "}":
                raise ValueError(
                    f"Invalid expression: {sequence};" '"}" cannot be used here'
                )

            prev = char
            result += char
        elif prev == "{":
            prev = int(char)
            if not 0 <= prev <= MAX_REPEATS:
                raise ValueError(
                    f"Invalid expression: {sequence}; repeat value {prev} must be <= "
                    f"{MAX_REPEATS}"
                )
        elif isinstance(prev, int):
            if char != "}":
                raise ValueError(f"Invalid expression: {sequence}; " '"}" expected')

            result = result[:-1] + result[-1] * prev
            prev = None
        else:
            if char != "{":
                raise ValueError(f"Invalid expression: {sequence}; " 'Expected "{"')

            prev = "{"

    if isinstance(prev, int) or prev == "{":
        raise ValueError(f"Unterminated expression: {sequence}")

    return result


def _extract_name_from_spec(spec: str) -> Tuple[str, str]:
    """
    Parses an adapter specification given as 'name=adapt' into 'name' and 'adapt'.

    Args:
        spec: Adapter spec.

    Returns:
        The tuple (name, spec)
    """
    fields = spec.split("=", 1)
    if len(fields) > 1:
        name, spec = fields
        name = name.strip()
    else:
        name = None
    spec = spec.strip()
    return name, spec
