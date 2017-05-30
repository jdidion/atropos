#!/usr/bin/env nextflow

/* Atropos paper workflow for simulated DNA-Seq reads
 * --------------------------------------------------
 * Configuration is externalized into separate files for desktop and cluster.
 * The main difference is that Docker is used to manage images on desktop, 
 * verusus Singularity on the cluster.
 *
 * The read data is contained within a Docker container. During local execution, 
 * we attach that container directly to the tool container (using --from-volumes 
 * <data container>). Singularity does not allow for this, so when running on 
 * the cluster we use a separate process to unpack the data from the container 
 * into $params.storeDir. Note that storeDir persists beyond the end of the 
 * script, so we delete it in the shell script that runs the workflow.
 *
 * Parameters
 * ----------
 * The following are expected to be defined in params:
 * - threadCounts: The numbers of threads (cpus) to test
 * - errorRates: The error rates of the simulated data sets
 * - minLength: minimum read length
 * - batchSize: read batch size
 * - quals: quality thresholds for trimming
 * - aligners: Atropos aligner algorithms to test
 * - compressionSchemes: Atropos compression schemes to test
 * - adapter1, adapter2: Adapter sequence
 */

// An absolute path to the container image is required for Singularity but
// not Docker
params.containerPrefix = ""
params.containerSuffix = ""

// variables for all tools
params.errorRates = [ '001', '005', '01' ]
params.quals = [ 0 ]
params.adapter1 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATATCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000

// atropos-specific variables
params.aligners = [ 'insert', 'adapter' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

process ExtractReads {
  container {
    "${params.containerPrefix}jdidion/atropos_rnaseq${params.containerSuffix}"
  }
  
  output:
  set val("untrimmed"), file("rna.{1,2}.fq.gz") into rnaseqReads
  
  script:
  """
  gunzip -c /data/rna/rna.1.fq.gz | head -40 | gzip > rna.1.fq.gz
  gunzip -c /data/rna/rna.2.fq.gz | head -40 | gzip > rna.2.fq.gz
  """
}

rnaseqReads.into {
  untrimmedRnaseqReads
  atroposRnaseqReads
}

process ExtractAnnotations {
  container {
    "${params.containerPrefix}jdidion/hg38_reference${params.containerSuffix}"
  }
  
  output:
  file "gencode.v26.annotation.gtf" into annotations
  
  script:
  """
  head /data/annotations/hg38/gencode.v26.annotation.gtf > ./gencode.v26.annotation.gtf
  """
}

process Atropos {
  //tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  tag { "atropos_${task.cpus}_rnaseq_q0_insert_writer" }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/atropos_paper${params.containerSuffix}"
  }
  
  input:
  set val(_ignore_), file(reads) from atroposRnaseqReads
  each threads from params.threadCounts
  //each qcut from params.quals
  //each aligner from params.aligners
  //each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  set val("${task.tag}"), val("trim"), file("${task.tag}.machine_info.txt") into machineAtropos
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt \
  && /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    --op-order GACQW -T 4 \
    -e 0.20 --aligner insert --insert-match-error-rate 0.30 -q 0 --trim-n \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -m 25 --no-default-adapters --no-cache-adapters --log-level ERROR \
    --correct-mismatches liberal \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz \
    --report-file ${task.tag}.report.txt --quiet \
    -pe1 ${reads[0]} -pe2 ${reads[1]}
  """
}

// concatenate all of the timing results into a single channel
Channel
  .empty()
  .concat(
    untrimmedRnaseqReads,
    trimmedAtropos
  )
  .set { trimmedMerged }

process StarAlign {
  tag { "${name}.star" }
  cpus { params.alignThreads }
  
  input:
  set val(name), file(fastq) from trimmedMerged
  
  output:
  file("${name}_rnaseq_Aligned.{bam,bam.bai}")
  set val(name), file("${name}.name_sorted.bam") into sorted
  set val(name), file("${name}.star.timing.txt") into timingStar
  set val("${name}"), val("star"), file("${task.tag}.machine_info.txt") into machineBwameth
  
  script:
  """
  touch ${task.tag}.machine_info.txt
  touch "${name}_rnaseq_Aligned.bam"
  touch "${name}_rnaseq_Aligned.bam.bai"
  touch "${name}.name_sorted.bam"
  touch "${name}.star.timing.txt"
  touch "${task.tag}.machine_info.txt"
  """
}

// clone sorted bams
sorted.into {
  sortedBam2Bed
  sortedEffectiveness
}

/* Process: Convert aligned reads into table format
 * ------------------------------------------------
 */
process Bam2Bed {
  container {
    "${params.containerPrefix}jdidion/atropos_paper_analysis${params.containerSuffix}"
  }
  
  input:
  set val(name), file(sortedBam) from sortedBam2Bed
  file annoFile from annotations
  
  output:
  file "${name}.overlap.txt" into overlap
  
  script:
  """
  bam2bed --all-reads --do-not-sort < $sortedBam \
    | cut -f 1-6 | bedmap --delim '\t' --echo --echo-map-id - $annoFile \
    > ${name}.overlap.txt
  """
}

Channel
  .empty()
  .concat(timingAtropos)
  .set { timingMerged }

process ParseTiming {
    container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
    
    input:
    set val(item), file(timing) from timingMerged
    
    output:
    stdout timingParsed
    
    script:
    """
    parse_gtime.py -i $timing -p $item
    """
}

/* Channel: display names for tools
 * --------------------------------
 */
toolNames = Channel.fromPath(
  "${workflow.projectDir}/../containers/tools/tool-names.txt")

process ShowPerformance {
    container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
    
    input:
    val timingRows from timingParsed.toList()
    
    output:
    file "timing.tex"
    file "timing.svg"
    
    script:
    data = timingRows.join("")
    """
    echo '$data' | show_performance.py -o timing -f tex svg pickle
    """
}

/* Process: Summarize trimming effectiveness
 * -----------------------------------------
 */
process ComputeEffectiveness {
  input:
  val bamFileList from sortedEffectiveness.toList()
  val bedFileList from overlap.toList()
  
  output:
  file "effectiveness.txt" into effectiveness
  
  script:
  bamMap = bamFileList.collectEntries()
  bedMap = bedFileList.collectEntries()
  names = bamMap.keySet().collect()
  bamFiles = names.collect { name -> bamMap[name] }
  bedFiles = names.collect { name -> bedMap[name] }
  """
  compute_real_effectiveness.py \
    -i $bamFiles -o effectiveness.txt \
    --no-edit-distance --no-progress \
    mrna -b $bedFiles
  """
}

/* Process: Generate plot
 * ----------------------
 */
process ShowEffectiveness {
  input:
  file effData from effectiveness
  file toolNamesFile from toolNames
  
  output:
  file "rnaseq_effectiveness.svg"
  
  script:
  """
  show_rnaseq_effectiveness.py \
    -i $effData -o rnaseq_effectiveness \
    -t $toolNamesFile --exclude-discarded -f svg pickle
  """
}

/* Channel: merged machine info
 * ----------------------------
 * If you add a tool process, make sure to add the machine channel
 * into the concat list here.
 */
Channel
  .empty()
  .concat(
    machineAtropos,
    machineBwameth
  )
  .set { machineMerged }

/* Process: summarize machine info
 * -------------------------------
 */
process SummarizeMachine {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
    
  input:
  set val(name), val(analysis), file(machine) from machineMerged

  output:
  stdout machineParsed

  script:
  """
  parse_machine.py -i $machine -p $name $analysis
  """
}

process CreateMachineTable {
  input:
  val parsedRows from machineParsed.toList()
    
  output:
  file "machine_info.txt"
  
  script:
  data = parsedRows.join("")
  """
  echo -e "prog\tprog2\tthreads\tdataset\tqcut\tanalysis\tcpus\tmemory\tcpu_details" > machine_info.txt
  echo '$data' >> machine_info.txt
  """
}
