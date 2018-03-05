#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/atropos:1.1.17
- class: InlineJavascriptRequirement

baseCommand: [atropos, trim]

inputs:
  paired_input1:
    type: File
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930", # FASTQ
      "http://edamontology.org/format_2572", # BAM
      "http://edamontology.org/format_2573"  # SAM
    ]
    streamable: true
    inputBinding:
      prefix: -pe1
  paired_output1:
    type: File
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930", # FASTQ
      "http://edamontology.org/format_2572", # BAM
      "http://edamontology.org/format_2573"  # SAM
    ]
    inputBinding:
      prefix: -o
  paired_input2:
    type: File?
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930"  # FASTQ
    ]
    streamable: true
    inputBinding:
      prefix: -pe2
  paired_output2:
    type: File?
    format: [
      "http://edamontology.org/format_1929", # FASTA
      "http://edamontology.org/format_1930"  # FASTQ
    ]
    inputBinding:
      prefix: -p
  adapter1:
    type: [File, string]
    format: "http://edamontology.org/format_1929" # FASTA
    streamable: true
    inputBinding:
      prefix: -a
  adapter2:
    type:
    - "null"
    - File
    - string
    format: "http://edamontology.org/format_1929" # FASTA
    streamable: true
    inputBinding:
      prefix: -A
  threads:
    type: int?
    inputBinding:
      prefix: -T
    default: 1
  aligner:
    type: string?
    inputBinding:
      prefix: --aligner
    default: insert
  op_order:
    type: string?
    inputBinding:
      prefix: --op-order
    default: GACQW
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
  batch_size:
    type: int?
    inputBinding:
      prefix: --batch-size
    default: 1000
  report_file:
    type: File?
    format: [
      "http://edamontology.org/format_2330", # txt
      "http://edamontology.org/format_3464"  # JSON
    ]
    inputBinding:
      prefix: --report-file
  no_default_adatpers:
    type: boolean?
    inputBinding:
      prefix: --no-default-adapters
  no_cache_adapters:
    type: boolean?
    inputBinding:
      prefix: --no-cache-adapters
  log_level:
    type: string
    inputBinding:
      prefix: --log-level
    default: ERROR
  quiet:
    type: boolean?
    inputBinding:
      prefix: --quiet
  error_rate:
    type: double
    inputBinding:
      prefix: -e
    default: 0.1
  insert_match_error_rate:
    type: double?
    inputBinding:
      prefix: --insert-match-error-rate
    default: 0.2
  correct_mismatches:
    type: string?
    inputBinding:
      prefix: --correct-mismatches
  overwrite_low_quality:
    type: string?
    inputBinding:
      prefix: -w

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
  report:
    type: File?
    outputBinding:
      glob: $(inputs.report_file)
