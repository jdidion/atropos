#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: ngs-bits/seqpurge:1.1.1

baseCommand: [SeqPurge]