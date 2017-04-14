#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/atropos-1.1.1

baseCommand: [atropos, trim]

inputs:
  paired-input1:
    type: File
    format: [
      http://edamontology.org/format_1929, # FASTA
      http://edamontology.org/format_1930, # FASTQ
      http://edamontology.org/format_2572, # BAM
      http://edamontology.org/format_2573  # SAM
    ]
    streamable: true
    inputBinding:
      prefix: -pe1
  paired-input2:
    type: File
    format: [
      http://edamontology.org/format_1929, # FASTA
      http://edamontology.org/format_1930  # FASTQ
    ]
    streamable: true
    inputBinding:
      prefix: -pe2
  adapter1:
    type: [File, string]
    format: http://edamontology.org/format_1929 # FASTA
    streamable: true
  adapter2:
    type: [File, string]
    format: http://edamontology.org/format_1929 # FASTA
    streamable: true
  threads:
    type: int?
    inputBinding:
      prefix: -T
    default: 1
  

outputs:
  trimmed1:
    type: File
    format: [
      http://edamontology.org/format_1929, # FASTA
      http://edamontology.org/format_1930, # FASTQ
      http://edamontology.org/format_2572, # BAM
      http://edamontology.org/format_2573  # SAM
    ]
    streamable: true
    outputBinding:
      glob: $(inputs.output_name1)
  trimmed2:
    type: File
    format: [
      http://edamontology.org/format_1929, # FASTA
      http://edamontology.org/format_1930  # FASTQ
    ]
    streamable: true
    outputBinding:
      glob: $(inputs.output_name2)
  report:
    type: [File]
    format: [
      http://edamontology.org/format_2330, # txt
      http://edamontology.org/format_3464  # JSON
    ]
    outputBinding:
      glob: $(inputs.report_name)
