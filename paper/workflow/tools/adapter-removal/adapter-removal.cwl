#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/seqpurge
- class: InlineJavascriptRequirement

baseCommand: [AdapterRemoval]

inputs:
  paired_input1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: --file1
  paired_output1:
    type: File
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -output1
  paired_input2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    streamable: true
    inputBinding:
      prefix: --file2
  paired_output2:
    type: File?
    format: "http://edamontology.org/format_1930" # FASTQ
    inputBinding:
      prefix: -output2
  gzip:
    type: boolean
    inputBinding:
      prefix: --gzip
  adapter1:
    type: string
    inputBinding:
      prefix: -adapter1
  adapter2:
    type: string?
    inputBinding:
      prefix: -adapter2
  mismatch_rate:
    type: double?
    inputBinding:
      prefix: --mm
  trimns:
    type: boolean
    inputBinding:
      prefix: --trimns
  trimqualities:
    type: boolean
    inputBinding:
      prefix: --trimqualities
  minquality:
    type: int?
    inputBinding:
      prefix: --minquality
  minlength:
    type: int?
    inputBinding:
      prefix: --minlength
  threads:
    type: int?
    inputBinding:
      prefix: --threads
    default: 1
  settings_file:
    type: File?
    inputBinding:
      prefix: --settings

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
  settings:
    type: File?
    outputBinding:
      glob: $(inputs.settings_file)
