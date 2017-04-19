#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/seqpurge
- class: InlineJavascriptRequirement

baseCommand: [SeqPurge]

inputs:
  paired_input1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: -in1
  paired_output1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -out1
  paired_input2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: -in2
  paired_output2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -out2
  adapter1:
    type: string
    inputBinding:
      prefix: -a1
  adapter2:
    type: string?
    inputBinding:
      prefix: -a2
  qcut:
    type: int?
    inputBinding:
      prefix: -qcut
  min_len:
    type: int?
    inputBinding:
      prefix: -min_len
  threads:
    type: int?
    inputBinding:
      prefix: -threads
    default: 1
  prefetch:
    type: int?
    inputBinding:
      prefix: -prefetch
  match_perc:
    type: double?
    inputBinding:
      prefix: -match_perc
  summary_file:
    type: File?
    inputBinding:
      prefix: -summary

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
  summary:
    type: File?
    outputBinding:
      glob: $(inputs.summary_file)
