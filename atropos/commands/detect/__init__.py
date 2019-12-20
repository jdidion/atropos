from abc import ABCMeta, abstractmethod
from collections import defaultdict
from enum import IntFlag
import math
import re
from typing import Iterable, Optional, Sequence as SequenceType, Tuple, Union

from loguru import logger

from atropos.adapters import AdapterCache
from atropos.aligner import Aligner, GapRule
from atropos.commands import (
    BaseCommand,
    Pipeline,
    PairedEndPipelineMixin,
    SingleEndPipelineMixin,
    Summary,
)
from atropos.io.sequence import Sequence
from atropos.utils import ReturnCode, classproperty, run_interruptible
from atropos.utils.collections import Summarizable
from atropos.utils.ngs import reverse_complement, sequence_complexity


class DetectCommand(BaseCommand):
    @classproperty
    def name(cls) -> str:
        return "detect"

    def __call__(self) -> ReturnCode:
        kmer_size = self.get_option("kmer_size", 12)
        n_reads = self.get_option("max_reads")
        overrep_cutoff = 100

        if self.get_option("include_contaminants"):
            include = Include(self.get_option("include_contaminants").upper())
        else:
            include = Include.ALL

        known_contaminants = None
        if include != Include.UNKNOWN:
            known_contaminants = self._load_known_adapters()

        detector = self.get_option("detector")
        if not detector:
            if known_contaminants and include == Include.KNOWN:
                detector = "known"
            elif n_reads <= 50000:
                detector = "heuristic"
            else:
                detector = "khmer"

        detector_args = dict(known_contaminants=known_contaminants)

        if detector == "known":
            logger.debug("Detecting contaminants using the known-only algorithm")
            detector_class = KnownContaminantDetector
            detector_args["min_kmer_match_frac"] = self.get_option(
                "min_kmer_match_frac"
            )
        elif detector == "heuristic":
            logger.debug("Detecting contaminants using the heuristic algorithm")
            detector_class = HeuristicDetector
            detector_args["min_frequency"] = self.get_option("min_frequency")
            detector_args["min_contaminant_match_frac"] = self.get_option(
                "min_contaminant_match_frac"
            )
        elif detector == "khmer":
            logger.debug("Detecting contaminants using the kmer-based algorithm")
            detector_class = KhmerDetector
        else:
            raise ValueError(f"Invalid value for 'detector': {detector}")

        summary_args = dict(
            kmer_size=kmer_size,
            n_reads=n_reads,
            overrep_cutoff=overrep_cutoff,
            include=include,
            past_end_bases=self.get_option("past_end_bases"),
        )
        detector_args.update(summary_args)

        if self.get_option("paired"):
            detector = PairedDetector(detector_class, **detector_args)
        else:
            detector = detector_class(**detector_args)

        self.summary["detect"] = summary_args
        if known_contaminants:
            self.summary["detect"][
                "known_contaminants"
            ] = known_contaminants.summarize()

        # currently only single-threaded operation is supproted
        self.summary.update(mode="serial", threads=1)

        logger.info(
            "Detecting adapters and other potential contaminant "
            "sequences based on %d-mers in %d reads",
            kmer_size,
            n_reads,
        )
        return run_interruptible(detector, self, raise_on_error=True)


class ContaminantMatcher:
    """
    Matches a known contaminant against other sequences.
    """

    @staticmethod
    def create(contaminants, kmer_size) -> SequenceType["ContaminantMatcher"]:
        """
        Creates :class:`ContaminantMatcher`s from sequences.

        Args:
            contaminants: A dict of {seq:names}.
            kmer_size: The kmer size.

        Returns:
            A list of :class:`ContaminantMatcher`s.
        """
        return [
            ContaminantMatcher(seq, names, kmer_size)
            for seq, names in contaminants.sequence_name_pairs
        ]

    def __init__(self, seq: str, names, kmer_size: int):
        """
        Args:
            seq: The contaminant sequence.
            names: Sequence of names for the contaminant.
            kmer_size: kmer size.
        """
        self.seq = seq
        self.names = names
        self.kmers = set(
            seq[i:(i + kmer_size)] for i in range(len(seq) - kmer_size + 1)
        )
        self.n_kmers = len(self.kmers)
        self.kmer_size = kmer_size
        self.matches = 0

    def match(self, seq: str, seqrc: str) -> Tuple[float, float, str]:
        """
        Matches a sequence.

        Args:
            seq: The sequence to match.
            seqrc: The reverse complement of `seq`.

        Returns:
            Tuple (f1, f2, seq), where f1 is the fraction of contaminant kmers
            that match, f2 is the fraction of sequence kmers that match, and
            seq is the best matching sequence (either `seq` or `seqrc`).
        """
        fw_kmers = set(
            seq[i:(i + self.kmer_size)] for i in range(len(seq) - self.kmer_size + 1)
        )
        fw_matches = float(len(self.kmers & fw_kmers))

        rv_kmers = set(
            seqrc[i:(i + self.kmer_size)]
            for i in range(len(seqrc) - self.kmer_size + 1)
        )
        rv_matches = float(len(self.kmers & rv_kmers))

        if fw_matches >= rv_matches:
            n_matches = fw_matches
            kmers = fw_kmers
            compare_seq = seq
        else:
            n_matches = rv_matches
            kmers = rv_kmers
            compare_seq = seqrc

        self.matches += n_matches

        match_frac1 = n_matches / self.n_kmers if self.n_kmers > 0 else 0
        match_frac2 = n_matches / len(kmers) if len(kmers) > 0 else 0

        return match_frac1, match_frac2, compare_seq


class Match(Summarizable):
    """
    A contaminant match.
    """

    def __init__(
        self,
        seq_or_contam: Union[str, ContaminantMatcher],
        count: int = 0,
        names: Optional[SequenceType[str]] = None,
        match_frac: Optional[float] = None,
        match_frac2: Optional[float] = None,
        abundance: Optional[float] = None,
        reads: Optional[Iterable[str]] = None,
    ):
        """
        Args:
            seq_or_contam: The matched sequence.
            count: The number of matches.
            names: The name(s) of matching contaminant(s).
            match_frac: The fraction of contaminant kmers that match.
            match_frac2: The fraction of sequence kmers that match.
            reads: The number of reads with the contaminant.
        """
        if isinstance(seq_or_contam, ContaminantMatcher):
            self.seq = seq_or_contam.seq
            self.count = int(seq_or_contam.matches)
            self.names = tuple(seq_or_contam.names)
            self.known_seqs = [seq_or_contam.seq]
        else:
            self.seq = seq_or_contam
            self.count = count
            self.names = tuple(names) if names else None
            self.known_seqs = None

        self.match_frac = match_frac
        self.match_frac2 = match_frac2
        self.abundance = abundance
        self.longest_match = None

        if reads:
            self.set_longest_match(reads)

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:
        if self.is_known:
            return f"{self.seq} => {self.names} ({self.known_seqs}))"
        else:
            return self.seq

    @property
    def seq_complexity(self) -> float:
        """
        The complexity of the sequence (0=homopolymer, 2=random).
        """
        return sequence_complexity(self.seq)

    @property
    def count_is_frequency(self) -> bool:
        """
        Whether `self.count` represents a frequency.
        """
        return isinstance(self.count, float)

    def set_contaminant(
        self,
        contam: ContaminantMatcher,
        match_frac: float,
        match_frac2: Optional[float] = None,
    ):
        """
        Set the known contaminant from a :class:`ContaminantMatcher`.
        """
        self.set_known(contam.names, [contam.seq], match_frac, match_frac2)

    def set_known(
        self,
        names: Iterable[str],
        seqs: SequenceType[str],
        match_frac: float,
        match_frac2: Optional[float] = None,
    ):
        """
        Sets the known contaminant.
        """
        self.names = tuple(names) if names else None
        self.known_seqs = seqs
        self.match_frac = match_frac
        self.match_frac2 = match_frac2

    @property
    def is_known(self) -> bool:
        """
        Whether the matching sequence is a known contaminant.
        """
        return self.known_seqs is not None

    def set_longest_match(self, sequences: Iterable[str]):
        """
        Set the longest matching sequence from a set of sequences.
        """
        for seq in sequences:
            idx = seq.index(self.seq)
            seqlen = len(self.seq) - idx
            if self.longest_match is None or self.longest_match[1] < seqlen:
                self.longest_match = (seq[idx:], seqlen)

    def estimate_abundance(self, read_sequences: Iterable[str]):
        """
        Determines whether this match's sequence is within 'seq' by simple exact
        string comparison.
        """
        self.abundance = sum(1 for read_seq in read_sequences if self.seq in read_seq)

    def summarize(self) -> dict:
        summary = dict(
            longest_kmer=self.seq,
            kmer_freq=self.count,
            kmer_freq_type="frequency" if self.count_is_frequency else "count",
            abundance=self.abundance,
            is_known=self.is_known,
            known_to_contaminant_match_frac=None,
            contaminant_to_known_match_frac=None,
            longest_match=None,
            known_names=None,
            known_seqs=None,
        )

        if self.longest_match:
            summary.update(longest_match=self.longest_match[0])

        if self.is_known:
            summary.update(
                known_to_contaminant_match_frac=self.match_frac,
                contaminant_to_known_match_frac=self.match_frac2,
                known_names=self.names,
                known_seqs=self.known_seqs,
            )

        return summary


class Include(IntFlag):
    KNOWN = 1
    UNKNOWN = 2
    ALL = KNOWN | UNKNOWN


class Detector(SingleEndPipelineMixin, Pipeline, metaclass=ABCMeta):
    """
    Base class for contaminant detectors.
    """

    def __init__(
        self,
        kmer_size: int = 12,
        n_reads: int = 10000,
        overrep_cutoff: int = 100,
        include: Include = Include.ALL,
        known_contaminants=None,
        past_end_bases: SequenceType[str] = ("A",),
    ):
        """
        Args:
            kmer_size: Size of kmers to match.
            n_reads: Number of reads to sample.
            overrep_cutoff: Degree of overrepresentation required for a kmer to be
                considered as a contaminant.
            known_contaminants: :class:`ContaminantMatcher`s to match against.
            past_end_bases: On Illumina, long runs of A (and sometimes other bases)
                can signify that the sequencer has read past the end of a fragment
                that is shorter than the read length + adapter length. Those
                bases will be removed from any sequencers before looking for
                matching contaminants.
        """
        super().__init__()

        self.kmer_size = kmer_size
        self.n_reads = n_reads
        self.overrep_cutoff = overrep_cutoff
        self.include = include
        self.known_contaminants = known_contaminants
        self._read_length = None
        self._read_sequences = set()
        self._matches = None
        self._past_end_regexp = None

        if past_end_bases:
            if len(past_end_bases[0]) > 1:
                self._past_end_regexp = re.compile(past_end_bases[0])
            else:
                self._past_end_regexp = re.compile(
                    "|".join(
                        base + "{8,}.*|" + base + "{2,}$" for base in past_end_bases
                    )
                )

    @property
    @abstractmethod
    def min_report_freq(self):
        """
        The minimum contaminant frequency required for a contaminant to be reported.
        """

    def set_read_length(self, record: Sequence):
        if self._read_length is not None:
            raise RuntimeError("Cannot set read length once it has already been set")

        self._read_length = len(record.sequence)

    def handle_records(self, context: dict, records: SequenceType[Sequence]):
        if context["size"] == 0:
            return

        if self._read_length is None:
            self.set_read_length(records[0])

        super().handle_records(context, records)

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        seq = self._filter_seq(read1.sequence)
        if seq:
            self._read_sequences.add(seq)

    def _filter_seq(self, seq: str) -> Optional[str]:
        if sequence_complexity(seq) <= 1.0:
            return None

        if self._past_end_regexp:
            match = self._past_end_regexp.search(seq)
            if match:
                seq = seq[: match.start()]

        if len(seq) < self.kmer_size:
            return None

        return seq

    def matches(self, **kwargs) -> SequenceType[Match]:
        """
        Returns the current set of matches.
        """
        if self._matches is None or len(kwargs) > 0:
            self._filter_and_sort(**kwargs)

        return self._matches

    def _filter_and_sort(
        self,
        min_len: Optional[int] = None,
        min_complexity: float = 1.1,
        min_match_frac: float = 0.1,
        limit: int = 20,
    ) -> None:
        """
        Identify, filter, and sort contaminants.

        Args:
            min_len: Minimum contaminant length.
            min_complexity: Minimum sequence complexity.
            min_match_frac: Minimum fraction of matching kmers.
            limit: Maximum number of contaminants to return.
        """
        if min_len is None:
            min_len = self.kmer_size

        matches = self._get_contaminants()
        for match in matches:
            match.estimate_abundance(self._read_sequences)

        def _filter(_match: Match):
            if _match.count < self.min_report_freq:
                logger.debug(
                    f"Filtering {_match.seq} because frequency {_match.count} < "
                    f"{self.min_report_freq}"
                )
                return False

            if min_len and len(_match) < min_len:
                logger.debug(
                    f"Filtering {_match.seq} because it's too short (< {min_len} bp)"
                )
                return False

            if min_complexity and _match.seq_complexity < min_complexity:
                logger.debug(
                    f"Filtering {_match.seq} because its complexity "
                    f"{_match.seq_complexity} < {min_complexity}"
                )
                return False

            if self.include == "known" and not _match.is_known:
                logger.debug(f"Filtering {_match.seq} because it's not known")
                return False

            elif self.include == "unknown" and _match.is_known:
                logger.debug(f"Filtering {_match.seq} because it's known")
                return False

            if (
                min_match_frac
                and _match.is_known
                and _match.match_frac < min_match_frac
            ):
                logger.debug(
                    f"Filtering {_match.seq} because its match_frac "
                    f"{_match.match_frac} < {min_match_frac}"
                )
                return False

            return True

        matches = list(
            sorted(
                filter(_filter, matches),
                key=lambda x: len(x) * math.log(x.count),
                reverse=True,
            )
        )

        if limit is not None:
            matches = matches[:limit]

        self._matches = matches

    @abstractmethod
    def _get_contaminants(self) -> Iterable[Match]:
        """
        Implemention of contaminant matching.

        Returns:
            A list of :class:`Match`es.
        """

    def finish(self, summary: Summary, **kwargs):
        super().finish(summary)
        summary["detect"]["matches"] = (
            [match.summarize() for match in self.matches(**kwargs)],
        )


class PairedDetector(PairedEndPipelineMixin, Pipeline):
    """
    Detector for paired-end reads.
    """

    def __init__(self, detector_class, **kwargs):
        super().__init__()
        self.read1_detector = detector_class(**kwargs)
        self.read2_detector = detector_class(**kwargs)
        self._read_length_set = False

    def handle_records(self, context, records):
        if context["size"] == 0:
            return

        if not self._read_length_set:
            read1, read2 = records[0]
            self.read1_detector.set_read_length(read1)
            self.read2_detector.set_read_length(read2)
            self._read_length_set = True
        super().handle_records(context, records)

    def handle_reads(
        self, context: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        self.read1_detector.handle_reads(context, read1)
        self.read2_detector.handle_reads(context, read2)

    def finish(self, summary, **kwargs):
        super().finish(summary)
        summary["detect"]["matches"] = (
            [match.summarize() for match in self.read1_detector.matches(**kwargs)],
            [match.summarize() for match in self.read2_detector.matches(**kwargs)],
        )


class KnownContaminantDetector(Detector):
    """
    Tests known contaminants against reads. This has linear complexity and is
    more specific than the khmer matcher, but less specific than the heuristic
    matcher. It's also less sensitive since it does not try to detect unknown
    contaminants.
    """

    def __init__(
        self,
        known_contaminants: AdapterCache,
        min_kmer_match_frac: float = 0.5,
        **kwargs
    ):
        """
        Args:
            known_contaminants: List of :class:`ContaminantMatcher`s.
            min_kmer_match_frac: Minimum fraction of matching kmers required.
            kwargs: Additional arguments to pass to the :class:`Detector`
                constructor.
        """
        super().__init__(known_contaminants=known_contaminants, **kwargs)
        self.min_kmer_match_frac = min_kmer_match_frac
        self._min_k = min(len(s) for s in known_contaminants.sequences)

    @property
    def min_report_freq(self) -> float:
        return 0.1

    def _filter_seq(self, seq: str) -> Optional[str]:
        seq = super()._filter_seq(seq)
        if seq and len(seq) >= self._min_k:
            return seq

        return None

    def _get_contaminants(self) -> Iterable[Match]:
        contaminant_matchers = ContaminantMatcher.create(
            self.known_contaminants, self.kmer_size
        )
        counts = defaultdict(int)
        max_match_fracs = defaultdict(float)

        for seq in self._read_sequences:
            seqrc = reverse_complement(seq)
            for contam in contaminant_matchers:
                match = contam.match(seq, seqrc)
                if match[0] > self.min_kmer_match_frac:
                    counts[contam] += 1
                    if match[0] > max_match_fracs[contam]:
                        max_match_fracs[contam] = match[0]

        min_count = math.ceil(
            self.n_reads
            * (self._read_length - self._min_k + 1)
            * self.overrep_cutoff
            / float(4 ** self._min_k)
        )

        return [
            Match(
                c[0],
                match_frac=max_match_fracs[c[0]],
                abundance=float(c[1]) / self.n_reads,
            )
            for c in filter(lambda x: x[1] >= min_count, counts.items())
        ]


class HeuristicDetector(Detector):
    """
    Uses a heuristic iterative algorithm to arrive at likely contaminants.
    This is the most accurate algorithm overall, but it has quadratic complexity
    and becomes too slow/memory-intenstive when n_reads > 50k.
    """

    def __init__(
        self,
        min_frequency: float = 0.001,
        min_contaminant_match_frac: float = 0.9,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.min_frequency = min_frequency
        self.min_contaminant_match_frac = min_contaminant_match_frac

    @property
    def min_report_freq(self) -> float:
        return 0.1 * self.n_reads

    def _get_contaminants(self) -> Iterable[Match]:
        def _min_count(_kmer_size: int):
            return math.ceil(
                self.n_reads
                * max(
                    self.min_frequency,
                    (
                        (self._read_length - _kmer_size + 1)
                        * self.overrep_cutoff
                        / float(4 ** _kmer_size)
                    ),
                )
            )

        kmer_size = self.kmer_size
        kmers = defaultdict(lambda: [0, set()])
        for seq in self._read_sequences:
            for i in range(len(seq) - kmer_size + 1):
                kmer = seq[i:(i + kmer_size)]
                kmers[kmer][0] += 1
                kmers[kmer][1].add(seq)

        prev = None
        cur = {}
        results = {}
        result_seqs = defaultdict(set)
        min_count = _min_count(kmer_size)

        # Identify candidate kmers for increasing values of k
        while True:
            all_seqs = set()

            for kmer, (count, seqs) in kmers.items():
                if count > min_count:
                    cur[kmer] = (count, seqs)
                    all_seqs.update(seqs)

            if len(all_seqs) == 0:
                break

            if prev:
                for kmer, (count, seqs) in prev.items():
                    if (
                        not any(seq in cur for seq in seqs)
                        and sequence_complexity(kmer) > 1.0
                    ):
                        results[kmer] = count
                        result_seqs[kmer].update(seqs)

            kmer_size += 1
            kmers = defaultdict(lambda: [0, set()])

            for seq in all_seqs:
                for i in range(len(seq) - kmer_size + 1):
                    kmer = seq[i:(i + kmer_size)]
                    kmers[kmer][0] += 1
                    kmers[kmer][1].add(seq)

            min_count = _min_count(kmer_size)
            prev = cur
            cur = {}

        results = list(results.items())

        # Now merge overlapping sequences by length and frequency to eliminate
        # redundancy in the set of candidate kmers.
        results.sort(key=lambda r: len(r[0]) * math.log(r[1]), reverse=True)
        merged = []
        unmerged = []

        while len(results) > 1:
            seq1, count1 = results[0]
            for j in range(1, len(results)):
                seq2, count2 = results[j]
                if len(seq1) >= len(seq2) and seq2 in seq1:
                    count1 += count2
                elif seq1 in seq2:
                    # if they are close in count, keep the longer sequence
                    if count1 < (2 * count2):
                        seq1 = seq2
                    count1 += count2
                else:
                    unmerged.append(results[j])

            merged.append([seq1, count1])
            results = unmerged
            unmerged = []

        results = merged + results
        if len(results) == 0:
            return []

        # TODO: For each retained match, pull out the longest sequence that
        #  matches to have a better shot of identifying long adapters that
        #  appear in full very infrequently

        # Re-sort by frequency
        results.sort(key=lambda r: r[1], reverse=True)
        # Keep anything that's within 50% of the top hit
        # TODO: make this user-configurable?
        min_count = int(results[0][1] * 0.5)
        results = (x for x in results if x[1] >= min_count)
        # Convert to matches
        matches = [Match(x[0], count=x[1], reads=result_seqs[x[0]]) for x in results]

        if self.known_contaminants:
            # Match to known sequences
            contaminants = ContaminantMatcher.create(
                self.known_contaminants, self.kmer_size
            )
            known = {}
            unknown = []

            def find_best_match(_seq, _best_matches, _best_match_frac):
                """
                Finds best contaminant matches to `seq`.
                """
                seqrc = reverse_complement(_seq)

                for _contam in contaminants:
                    match_frac1, match_frac2, compare_seq = _contam.match(_seq, seqrc)
                    if match_frac1 < _best_match_frac[0]:
                        continue

                    if _contam.seq in compare_seq or align(
                        compare_seq, _contam.seq, self.min_contaminant_match_frac
                    ):
                        if match_frac1 > _best_match_frac[0] or (
                            match_frac1 == _best_match_frac[0]
                            and match_frac2 > _best_match_frac[1]
                        ):
                            _best_matches = {}
                            _best_match_frac = (match_frac1, match_frac2)
                        _best_matches[_contam] = (match, (match_frac1, match_frac2))

                return _best_matches, _best_match_frac

            for match in matches:
                best_matches, best_match_frac = find_best_match(
                    match.seq, {}, (self.min_contaminant_match_frac, 0)
                )
                if match.longest_match:
                    best_matches, best_match_frac = find_best_match(
                        match.longest_match[0], best_matches, best_match_frac
                    )
                if best_matches:
                    for contam, _match in best_matches.items():
                        if contam not in known or _match[1] > known[contam][1]:
                            known[contam] = _match
                else:
                    unknown.append(match)

            # resolve many-many relationships
            new_matches = defaultdict(list)
            for contam, (match, match_frac) in known.items():
                new_matches[match].append((contam, match_frac))

            known = []

            for match, contams in new_matches.items():
                if len(contams) == 1:
                    contam, match_frac = contams[0]
                    match.set_contaminant(contam, *match_frac)
                else:
                    contams.sort(key=lambda x: x[1], reverse=True)
                    contam, match_frac = contams[0]
                    equiv = [
                        other_contam
                        for other_contam in contams[1:]
                        if other_contam[1] == match_frac
                    ]
                    if len(equiv) == 0:
                        match.set_contaminant(contam, *match_frac)
                    else:
                        names = {contam.names}
                        seqs = {(contam.seq,)}
                        for other_contam in equiv:
                            names.update(other_contam[0].names)
                            seqs.add(other_contam[0].seq)
                        match.set_known(list(names), list(seqs), *match_frac)
                known.append(match)

            matches = known + unknown

        return matches


class KhmerDetector(Detector):
    """
    Identifies contaminants based on kmer frequency using a fast kmer counting
    approach (as implemented in the khmer library). This approach is fast but
    not as accurate as the other two.
    """

    @property
    def min_report_freq(self) -> float:
        return 0.0001

    def _get_contaminants(self) -> Iterable[Match]:
        from khmer import Countgraph, khmer_args

        # assuming all sequences are same length
        n_win = self._read_length - self.kmer_size + 1
        tablesize = self.n_reads * n_win
        n_expected = math.ceil(tablesize / float(4 ** self.kmer_size))
        min_count = n_expected * self.overrep_cutoff
        if min_count >= 2 ** 16:
            raise ValueError(
                f"The minimum count for an over-represented k-kmer {min_count} is "
                f"greater than the max khmer count (2^16)"
            )

        countgraph = Countgraph(self.kmer_size, tablesize, khmer_args.DEFAULT_N_TABLES)
        countgraph.set_use_bigcount(True)
        for seq in self._read_sequences:
            countgraph.consume_and_tag(seq)

        candidates = {}
        for tag in countgraph.get_tagset():
            count = countgraph.get(tag)
            if count >= min_count:
                candidates[tag] = count

        if self.known_contaminants:
            matches = []
            seen = set()

            def match(_kmer):
                """Returns the frequency of `kmer` in `candidates`.
                """
                freq = candidates.get(_kmer, 0)
                if freq > 0:
                    seen.add(_kmer)
                return freq

            for seq, names in self.known_contaminants.sequence_name_pairs:
                seqlen = len(seq)
                if seqlen < self.kmer_size:
                    print(
                        f"Cannot check {list(names)[0]}; sequence is shorter than "
                        f"{self.kmer_size}"
                    )
                    continue

                n_kmers = seqlen - self.kmer_size + 1
                num_matches = 0
                match_counts = []

                for idx in range(n_kmers):
                    kmer = seq[idx:(idx + self.kmer_size)]
                    kmer_count = max(match(kmer), match(reverse_complement(kmer)))
                    if kmer_count > 0:
                        num_matches += 1
                        match_counts.append(kmer_count)

                if num_matches > 0:
                    # not sure what the correct metric is to use here
                    overall_count = sum(match_counts) / float(n_kmers)
                    matches.append(
                        Match(
                            seq,
                            count=overall_count / float(tablesize),
                            names=names,
                            match_frac=float(num_matches) / n_kmers,
                        )
                    )

            # Add remaining tags
            for tag in set(candidates.keys()) - seen:
                matches.append(Match(tag, count=candidates[tag] / float(tablesize)))
        else:
            matches = [
                Match(tag, count=count / float(tablesize))
                for tag, count in candidates.items()
            ]

        return matches


def align(seq1: str, seq2: str, min_overlap_frac: float = 0.9) -> Optional[str]:
    """
    Aligns two sequences.

    Args:
        seq1: The first sequence to align.
        seq2: The second sequence to align.
        min_overlap_frac: Minimum fraction of overlapping bases required for a
            match.

    Returns:
        The matching portion of the sequence.
    """
    aligner = Aligner(seq1, 0.0, GapRule.SEMIGLOBAL, False, False)
    aligner.min_overlap = math.ceil(min(len(seq1), len(seq2)) * min_overlap_frac)
    aligner.indel_cost = 100000
    match = aligner.locate(seq2)
    if match:
        return seq1[match[0]:match[1]]
