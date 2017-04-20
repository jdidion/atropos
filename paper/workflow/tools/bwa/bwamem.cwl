#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/bwaindex:latest
- class: InlineJavascriptRequirement

baseCommand: [bwa, mem]

inputs:
  read_group:
    type: string?
    inputBinding:
        position: 0
        prefix: -R
  mark_secondary:
    type: boolean
    inputBinding:
      position: 0
      prefix: -M
  threads:
    type: int?
    inputBinding:
      position: 0
      prefix: -t
  output_filename:
    type: File
  reference:
    type: File?
    inputBinding:
      position: 1
    default: /genomes/hg38.fa
  reads:
    type: File[]
    inputBinding:
      position: 2

stdout: $(inputs.output_filename)

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
