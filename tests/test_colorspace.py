from .utils import run, datapath

from atropos.utils.colorspace import encode, decode
from atropos.console import execute_cli

# If there are any unknown characters in the test sequence,
# round tripping will only work if all characters after the
# first unknown character are also unknown:
# encode("TNGN") == "T444", but
# decode("T444") == "TNNN".
sequences = [
    "",
    "C",
    "ACGGTC",
    "TN",
    "TN.",
    "TNN.N",
    "CCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTCGGTGGTA",
    "TGCAGTTGATGATCGAAGAAAACGACATCATCAGCCAGCAAGTGC",
    "CAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGG",
]


def test_encode():
    assert encode("AA") == "A0"
    assert encode("AC") == "A1"
    assert encode("AG") == "A2"
    assert encode("AT") == "A3"
    assert encode("CA") == "C1"
    assert encode("CC") == "C0"
    assert encode("CG") == "C3"
    assert encode("CT") == "C2"
    assert encode("GA") == "G2"
    assert encode("GC") == "G3"
    assert encode("GG") == "G0"
    assert encode("GT") == "G1"
    assert encode("TA") == "T3"
    assert encode("TC") == "T2"
    assert encode("TG") == "T1"
    assert encode("TT") == "T0"
    assert encode("TN") == "T4"
    assert encode("NT") == "N4"
    assert encode("NN") == "N4"
    assert encode("ACGGTC") == "A13012"
    assert encode("TTT.N") == "T0044"
    assert encode("TTNT.N") == "T04444"


def test_decode():
    for s in sequences:
        expected = s.replace(".", "N")
        encoded = encode(s)
        assert decode(encoded) == expected
    assert decode("A.") == "AN"
    assert decode("C.") == "CN"
    assert decode("G.") == "GN"
    assert decode("T.") == "TN"


def test_qualtrim_csfastaqual(tmp_path):
    """-q with csfasta/qual files"""
    run(
        "-c -q 10", "solidqual.fastq", "solid.csfasta", qualfile="solid.qual",
        output_dir=tmp_path
    )


def test_e3m():
    """Read the E3M dataset"""
    # not really colorspace, but a fasta/qual file pair
    execute_cli(
        ["-o", "/dev/null", "-se", datapath("E3M.fasta"), "-sq", datapath("E3M.qual")]
    )


def test_bwa(tmp_path):
    """MAQ-/BWA-compatible output"""
    run(
        "-c -e 0.12 -a 330201030313112312 -x 552: --maq",
        "solidmaq.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
        output_dir=tmp_path,
    )


def test_bfast(tmp_path):
    """BFAST-compatible output"""
    run(
        "-c -e 0.12 -a 330201030313112312 -x abc: --strip-f3",
        "solidbfast.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
        output_dir=tmp_path,
    )


def test_trim_095(tmp_path):
    """some reads properly trimmed since atropos 0.9.5"""
    run(
        "-c -e 0.122 -a 330201030313112312", "solid.fasta", "solid.fasta",
        output_dir=tmp_path
    )


def test_solid(tmp_path):
    run(
        "-c -e 0.122 -a 330201030313112312", "solid.fastq", "solid.fastq",
        output_dir=tmp_path
    )


def test_solid_basespace_adapter(tmp_path):
    """colorspace adapter given in basespace"""
    run(
        "-c -e 0.122 -a CGCCTTGGCCGTACAGCAG", "solid.fastq", "solid.fastq",
        output_dir=tmp_path
    )


def test_solid5p(tmp_path):
    """test 5' colorspace adapter"""
    # this is not a real adapter, just a random string
    # in colorspace: C0302201212322332333
    run(
        "-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA",
        "solid5p.fasta",
        "solid5p.fasta",
        output_dir=tmp_path,
    )


def test_solid5p_prefix_notrim(tmp_path):
    """test anchored 5' colorspace adapter, no primer trimming"""
    run(
        "-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA",
        "solid5p-anchored.notrim.fasta",
        "solid5p.fasta",
        output_dir=tmp_path,
    )


def test_solid5p_prefix(tmp_path):
    """test anchored 5' colorspace adapter"""
    run(
        "-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA",
        "solid5p-anchored.fasta",
        "solid5p.fasta",
        output_dir=tmp_path,
    )


def test_solid5p_fastq(tmp_path):
    """test 5' colorspace adapter"""
    # this is not a real adapter, just a random string
    # in colorspace: C0302201212322332333
    run(
        "-c -e 0.1 --trim-primer -g CCGGAGGTCAGCTCGCTATA",
        "solid5p.fastq",
        "solid5p.fastq",
        output_dir=tmp_path,
    )


def test_solid5p_prefix_notrim_fastq(tmp_path):
    """test anchored 5' colorspace adapter, no primer trimming"""
    run(
        "-c -e 0.1 -g ^CCGGAGGTCAGCTCGCTATA",
        "solid5p-anchored.notrim.fastq",
        "solid5p.fastq",
        output_dir=tmp_path,
    )


def test_solid5p_prefix_fastq(tmp_path):
    """test anchored 5' colorspace adapter"""
    run(
        "-c -e 0.1 --trim-primer -g ^CCGGAGGTCAGCTCGCTATA",
        "solid5p-anchored.fastq",
        "solid5p.fastq",
        output_dir=tmp_path,
    )


def test_sra_fastq(tmp_path):
    """test SRA-formatted colorspace FASTQ"""
    run(
        "-c -e 0.1 --input-format sra-fastq -a CGCCTTGGCCGTACAGCAG",
        "sra.fastq",
        "sra.fastq",
        output_dir=tmp_path,
    )


def test_no_zero_cap(tmp_path):
    run(
        "--no-zero-cap -c -e 0.122 -a CGCCTTGGCCGTACAGCAG",
        "solid-no-zerocap.fastq",
        "solid.fastq",
        output_dir=tmp_path,
    )
