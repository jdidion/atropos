#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/cutadapt:1.16
- class: InlineJavascriptRequirement

baseCommand: [cutadapt]

inputs:

  paired_output1:
    type: File
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930", # FASTQ
    ]
    inputBinding:
      prefix: -o
  paired_output2:
    type: File?
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930"  # FASTQ
    ]
    inputBinding:
      prefix: -p
  adapter1:
    type: string
    inputBinding:
      prefix: -a
  adapter2:
    type:
    - "null"
    - string
    inputBinding:
      prefix: -A
  threads:
    type: int?
    inputBinding:
      prefix: -j
    default: 1
  quality:
    type: int?
    inputBinding:
      prefix: -q
    default: 0
  trim_n:
    type: boolean?
    inputBinding:
      prefix: --trim-n
  min_len:
    type: int?
    inputBinding:
      prefix: -m
    default: 0
  error_rate:
    type: double
    inputBinding:
      prefix: -e
    default: 0.1
  paired_input1:
    type: File
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930", # FASTQ
    ]
    streamable: true
    inputBinding:
      position: 99
  paired_input2:
    type: File?
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930"  # FASTQ
    ]
    streamable: true
    inputBinding:
      position: 100

outputs:
  trimmed1:
    type: File
    streamable: true
    outputBinding:
      glob: $(inputs.paired_output1)
  trimmed2:
    type: File?
    streamable: true
    outputBinding:
      glob: $(inputs.paired_output2)
