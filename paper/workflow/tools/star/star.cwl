#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: jdidion/starindex:latest
- class: InlineJavascriptRequirement

baseCommand: [STAR]

inputs:
  outMultimapperOrder:
    type: string?
    inputBinding:
      prefix: --outMultimapperOrder
    default: Random
  outSAMattrRGline:
    type: string?
    inputBinding:
      prefix: --outSAMattrRGline
  genomeDir:
    type: Directory
    inputBinding:
      prefix: --genomeDir
    default: /genomes/hg38
  outSAMunmapped:
    type: string[]?
    inputBinding:
      itemSeparator: ' '
      prefix: --outSAMunmapped
    default: [Within, KeepPairs]
  outFilterMultimapNmax:
    type: int?
    inputBinding:
      prefix: --outFilterMultimapNmax
    default: 100000
  outFileNamePrefix:
    type: string?
    inputBinding:
      prefix: --outFileNamePrefix
  readNameSeparator:
    type: string?
    inputBinding:
      prefix: --readNameSeparator
    default: "/"
  genomeLoad:
    type: string?
    inputBinding:
      prefix: --genomeLoad
  readFilesIn:
    type: File[]?
    inputBinding:
      itemSeparator: ' '
      prefix: --readFilesIn
  outSAMmultNmax:
    type: int?
    inputBinding:
      prefix: --outSAMmultNmax
    default: 1
  runThreadN:
    type: int?
    inputBinding:
      prefix: --runThreadN
    default: 1
  outSAMprimaryFlag:
    type: string?
    inputBinding:
      prefix: --outSAMprimaryFlag
  bamRemoveDuplicatesType:
    type: string?
    inputBinding:
      prefix: --bamRemoveDuplicatesType
  outSAMtype:
    type: string[]?
    inputBinding:
      itemSeparator: ' '
      prefix: --outSAMtype
    default: [BAM, Unsorted]
  readFilesCommand:
    type: string?
    inputBinding:
      prefix: --readFilesCommand
  outStd:
    type: string
    default: Log
    inputBinding:
      prefix: --outStd

outputs:
  aligned:
    type: File?
    outputBinding:
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          if (inputs.outSAMtype.indexOf("SAM") > -1) {
              return p+"Aligned.out.sam";
          } else {
           if ( inputs.outSAMtype.indexOf("SortedByCoordinate") > -1 )
              return p+"Aligned.sortedByCoord.out.bam";
            else
              return p+"Aligned.out.bam";
          }
        }
    secondaryFiles: |
      ${
         var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
         return [
           {"path": p+"Log.final.out", "class":"File"},
           {"path": p+"SJ.out.tab", "class":"File"},
           {"path": p+"Log.out", "class":"File"}
         ];
      }

  mappingstats:
    type: string?
    outputBinding:
      loadContents: true
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.final.out";
        }
      outputEval: |
        ${
          if (inputs.runMode == "genomeGenerate")
            return "";

          var s = self[0].contents.replace(/[ ]+.*?:\n|[ ]{2,}|\n$/g,"").
              split(/\n{1,2}/g).map(function(v){var s=v.split(/\|\t/g); var o={}; o[s[0]]=s[1]; return o;})
          return JSON.stringify(s);
        }
