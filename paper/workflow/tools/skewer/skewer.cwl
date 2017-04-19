#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/skewer
- class: InlineJavascriptRequirement

baseCommand: [skewer]

inputs:
  mode:
    type: string
    inputBinding:
      prefix: -m
      position: 0
    default: pe
  adapter1:
    type: string
    inputBinding:
      prefix: -x
      position: 0
  adapter2:
    type: string?
    inputBinding:
      prefix: -y
      position: 0
  min_len:
    type: int?
    inputBinding:
      prefix: -l
      position: 0
  qcut:
    type: int?
    inputBinding:
      prefix: -q
      position: 0
  filter_ns:
    type: boolean?
    inputBinding:
      prefix: -n
      position: 0
  error_rate:
    type: double?
    inputBinding:
      prefix: -r
      position: 0
  threads:
    type: int?
    inputBinding:
      prefix: -t
      position: 0
    default: 1
  zipped:
    type: boolean?
    inputBinding:
      prefix: -z
      position: 0
  quiet:
    type: boolean?
    inputBinding:
      prefix: --quiet
      position: 0
  output_prefix:
    type: File?
    inputBinding:
      prefix: -o
      position: 0
  paired_input1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      position: 1
  paired_input2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      position: 2

outputs:
  trimmed1:
    type: File
    streamable: true
    outputBinding:
      glob: $(inputs.output_prefix || inputs.paired_input1.replace(/([_\.]1)(\.fq|\.fastq)(\.gz)?/g,'')+"*-trimmed")).fastq1.gz
  trimmed2:
    type: File?
    streamable: true
    outputBinding:
      glob: $(inputs.output_prefix || inputs.paired_input2.replace(/([_\.]2)(\.fq|\.fastq)(\.gz)?/g,'')+"*-trimmed")).fastq2.gz
