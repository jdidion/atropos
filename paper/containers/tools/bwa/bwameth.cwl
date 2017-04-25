#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/bwaindex:latest
- class: InlineJavascriptRequirement

baseCommand: [bwameth.py]

inputs:
  reference:
    type: File?
    inputBinding:
      prefix: --reference
      position: 0
    default: /genomes/hg38.fa
  read_group:
    type: string?
    inputBinding:
      position: 0
      prefix: --read-group
  output_filename:
    type: File
    inputBinding:
      position: 0
      prefix: -o
  threads:
    type: int?
    inputBinding:
      position: 0
      prefix: -t
  nosort:
    type: boolean
    inputBinding:
      position: 0
      prefix: --nosort
  reads:
    type: File[]
    inputBinding:
      position: 1

stdout: $(inputs.output_filename)

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
