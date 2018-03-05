#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/fastp:0.12.3
- class: InlineJavascriptRequirement

baseCommand: [fastp]

inputs:
  paired_input1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: -i
  paired_output1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -o
  paired_input2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: -I
  paired_output2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -O
  quality:
    type: int?
    inputBinding:
      prefix: -q
    default: 0
  cut_mean_quality:
    type: int?
    inputBinding:
      prefix: --cut_mean_quality
    default: 20
  cut_quality5:
    type: int?
    inputBinding:
      prefix: --cut_by_quality5
  cut_quality3:
    type: int?
    inputBinding:
      prefix: --cut_by_quality3
  disable_length_filtering:
    type: boolean
    inputBinding:
      prefix: --disable_length_filtering
  adapter1:
    type: string
    inputBinding:
      prefix: --adapter_sequence
  adapter2:
    type: string?
    inputBinding:
      prefix: --adapter_sequence_r2
  correction:
    type: bool
    inputBinding:
      prefix: --correction
  threads:
    type: int?
    inputBinding:
      prefix: --thread
    default: 1

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
