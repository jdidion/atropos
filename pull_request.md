adding feature for extracting UMIs (issue #61)

I don't know if you are working on this, but this pull request is a prototype for adding a new feature to extract UMIs from the reads (issue #61). 

Basically, I added three new options to ```trim``` command to clip off the first *N* bases from read1/read2/both reads and append to the read ID.

```
UMI options, clipping UMI from sequence and append to read name:
  --read1_umi READ1_UMI
                        How many bases on the 5' end of read 1 are UMI?
                        (default: 0)
  --read2_umi READ2_UMI
                        How many bases on the 5' end of read 2 are UMI?
                        (default: 0)
  --delim DELIM         Deliminator for separating UMI from read ID (default:
                        ':')
```

In this prototype, I am introducing a new instance variable **UMI** into the **Sequence** object from ```_seqio.pyx```, such that UMIs can be synced between read1 and read2 when the inputs are paired-end. I also implemented two classes of **Modifiers** to:

1. trim *N* bases from sequences and store as **UMI** in the **Sequence** object (```class UmiTrimmer(Trimmer)```), and
2. Appending UMI to the read ID (```class SyncUMI(ReadPairModifier)``` for paired-end, or ```class AddUMI(Modifier)``` for single-end)

I am not sure if these are the best approaches or are there any better ways I can clip and append UMI in one modifier? or syncing the UMIs of a pair of reads without an extra modifier?

Usage example of the current implementation:

1. 4-nt UMIs from both ends:
```
$ atropos trim -pe1 tests/data/paired.1.fastq -pe2 tests/data/paired.2.fastq --read1_umi 4 --read2_umi 4 -o test.1.fq -p test.2.fq  -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG --delim ':'
```

Results:
```
$ cat test.1.fq
@read1/1:TTAT:GCTG some text
TTGTCTCCAGCTTAGACATATCGCCT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/1:CAAC:TGTG
AGGCCACATTAGACATATCGGATGGT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/1:CCAA:TGTT
CTTGATATTAATAACATTAG
+
HHHHHHHHHHHHHHHHHHHH
@read4/1:GACA:CATC
GGCCGTTTGAATGTTGACGGGATGTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
```

```
$ cat test.2.fq
@read1/2:TTAT:GCTG other text
GAGACAAATAACAGTGGAGTAGTTTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/2:CAAC:TGTG
GCCTGTTGCAGTGGAGTAACTCCAGC
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/2:CCAA:TGTT
ATTAATATCAAGTTGGCAGTG
+
HHHHHHHHHHHHHHHHHHHHH
@read4/2:GACA:CATC
CCGTCAACATTCAAACGGCCTGTCCA
+
##########################
```

2. 4-nt UMIs from read 1:
```
$ atropos trim -pe1 tests/data/paired.1.fastq -pe2 tests/data/paired.2.fastq --read1_umi 4 -o test.1.fq -p test.2.fq  -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG --delim ':'
```

Results:
```
$ cat test.1.fq                                       
@read1/1:TTAT some text
TTGTCTCCAGCTTAGACATATCGCCT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/1:CAAC
AGGCCACATTAGACATATCGGATGGT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/1:CCAA
CTTGATATTAATAACATTAG
+
HHHHHHHHHHHHHHHHHHHH
@read4/1:GACA
GGCCGTTTGAATGTTGACGGGATGTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
```

```
$ cat test.2.fq                                       
@read1/2:TTAT other text
GCTGGAGACAAATAACAGTGGAGTAGTTTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/2:CAAC
TGTGGCCTGTTGCAGTGGAGTAACTCCAGC
+
###HHHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/2:CCAA
TGTTATTAATATCAAGTTGGCAGTG
+
#HHHHHHHHHHHHHHHHHHHHHHHH
@read4/2:GACA
CATCCCGTCAACATTCAAACGGCCTGTCCA
+
HH############################
```


3. 4-nt UMIs from read 2:
```
$ atropos trim -pe1 tests/data/paired.1.fastq -pe2 tests/data/paired.2.fastq --read2_umi 4 -o test.1.fq -p test.2.fq  -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG --delim ':'
```

Results:
```
$ cat test.1.fq                                       
@read1/1:GCTG some text
TTATTTGTCTCCAGCTTAGACATATCGCCT
+
##HHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/1:TGTG
CAACAGGCCACATTAGACATATCGGATGGT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/1:TGTT
CCAACTTGATATTAATAACATTAG
+
HHHHHHHHHHHHHHHHHHHHHHHH
@read4/1:CATC
GACAGGCCGTTTGAATGTTGACGGGATGTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
```

```
$ cat test.2.fq                                       
@read1/2:GCTG other text
GAGACAAATAACAGTGGAGTAGTTTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/2:TGTG
GCCTGTTGCAGTGGAGTAACTCCAGC
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/2:TGTT
ATTAATATCAAGTTGGCAGTG
+
HHHHHHHHHHHHHHHHHHHHH
@read4/2:CATC
CCGTCAACATTCAAACGGCCTGTCCA
+
##########################
```

4. 4-nt UMIs from single end:
```
$ atropos trim -se tests/data/paired.1.fastq  --read1_umi 4 -o -  -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG --delim ':' --quiet
@read1/1:TTAT some text
TTGTCTCCAGCTTAGACATATCGCCT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read2/1:CAAC
AGGCCACATTAGACATATCGGATGGT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
@read3/1:CCAA
CTTGATATTAATAACATTAG
+
HHHHHHHHHHHHHHHHHHHH
@read4/1:GACA
GGCCGTTTGAATGTTGACGGGATGTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHH
```

So I guess if this is moving forward, do you want to include these in this pull request:

1. find a way to do more flexible UMI recognition (something like the interface of [UMI-tools](https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md))?
2. maybe incorporate some quality cutoff for UMI?
3. it is currently silently ignoring ```--read2_umi``` if it is single-end input, would it be more appropriate if a warning or an error is raised?