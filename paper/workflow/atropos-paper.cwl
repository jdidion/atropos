#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: Workflow

inputs:
  - id: p1
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the first end of paired end reads in fasta or fastq format

  - id: p2
    type:
      - "null"
      - type: array
        items: File
    description: list of files containing the second end of paired end reads in fasta or fastq format

  - id: output_prefix
    type: string
    description: prefix for output files. will output prefix.aligned.bam and prefix.aligned.stats

  - id: reference
    type: File
    description: "lobSTR's bwa reference files"

  - id: rg-sample
    type: string
    description: Use this in the read group SM tag

  - id: rg-lib
    type: string
    description: Use this in the read group LB tag

  - id: strinfo
    type: File
    description: File containing statistics for each STR.

  - id: noise_model
    type: File
    description: File to read noise model parameters from (.stepmodel)
    secondaryFiles:
      - "^.stuttermodel"

outputs:
  - id: bam
    type: File
    source: "#samindex/bam_with_bai"

  - id: bam_stats
    type: File
    source: "#lobSTR/bam_stats"

  - id: vcf
    type: File
    source: "#allelotype/vcf"

  - id: vcf_stats
    type: File
    source: "#allelotype/vcf_stats"

hints:
  - class: DockerRequirement
    dockerLoad: https://workbench.qr1hi.arvadosapi.com/collections/download/qr1hi-4zz18-x2ae13tsx5jqg8d/1nduktd8dpvhdpgsva82lje0i710kgzb6rttks5jldx7s2y7k9/7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613.tar
    dockerImageId: 7e0c0ae3bf4e70442f9b8eee816ec23426d9e1169a2925316e5c932745e21613

steps:
  - id: lobSTR
    run: lobSTR-tool.cwl
    inputs:
      - { id: p1, source: "#p1" }
      - { id: p2, source: "#p2" }
      - { id: output_prefix, source: "#output_prefix" }
      - { id: reference, source: "#reference" }
      - { id: rg-sample, source: "#rg-sample" }
      - { id: rg-lib, source: "#rg-lib" }
    outputs:
      - { id: bam }
      - { id: bam_stats }

  - id: samsort
    run: samtools-sort.cwl
    inputs:
      - { id: input, source: "#lobSTR/bam" }
      - { id: output_name, default: "aligned.sorted.bam" }
    outputs:
      - { id: output_file }

  - id: samindex
    run: samtools-index.cwl
    inputs:
      - { id: input, source: "#samsort/output_file" }
    outputs:
      - { id: bam_with_bai }

  - id: allelotype
    run: allelotype.cwl
    inputs:
      - { id: bam, source: "#samindex/bam_with_bai" }
      - { id: reference, source: "#reference" }
      - { id: output_prefix, source: "#output_prefix" }
      - { id: noise_model, source: "#noise_model" }
      - { id: strinfo, source: "#strinfo" }
    outputs:
      - { id: vcf }
      - { id: vcf_stats }
