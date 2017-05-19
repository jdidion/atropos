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

// variables for all tools
params.errorRates = [ '001', '005', '01' ]
params.quals = [ 0 ]
params.adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000

// atropos-specific variables
params.aligners = [ 'insert', 'adapter' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

process ExtractReads {
  container "jdidion/atropos_wgbs"
  
  output:
  set val("untrimmed"), file("wgbs.{1,2}.fq") into wgbsReads
  
  script:
  """
  gunzip -c /data/wgbs/wgbs.1.fq.gz | head -40 > ./wgbs.1.fq
  gunzip -c /data/wgbs/wgbs.2.fq.gz | head -40 > ./wgbs.2.fq
  """
}

wgbsReads.into {
  untrimmedWgbsReads
  atroposWgbsReads
}

process Atropos {
  //tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  tag { "atropos_${task.cpus}_wgbs_q0_insert_writer" }
  cpus { threads }
  container "jdidion/atropos_paper"
  
  input:
  set val(_ignore_), file(reads) from atroposWgbsReads
  each threads from params.threadCounts
  //each qcut from params.quals
  //each aligner from params.aligners
  //each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineAtropos
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
    untrimmedWgbsReads,
    trimmedAtropos
  )
  .set { trimmedMerged }

process BwamethAlign {
  tag { "${name}.bwameth" }
  container "jdidion/bwameth_hg38index"
  cpus { params.alignThreads }
  
  input:
  set val(name), file(fastq) from trimmedMerged
  
  output:
  file("${name}.sam")
  file("${name}.name_sorted.bam") into sortedBams
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineBwameth
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingBwameth
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt \
  && /usr/bin/time -v -o ${name}.bwameth.timing.txt bwameth.py \
    -t ${params.alignThreads} --read-group '${task.ext.readGroup}' \
    --reference /data/index/bwameth/hg38/hg38.fa \
    ${fastq[0]} ${fastq[1]} > ${name}.sam \
  && samtools view -Shb ${name}.sam | \
     samtools sort -n -O bam -@ ${params.alignThreads} \
       -o ${name}.name_sorted.bam -
  """
}

Channel
  .empty()
  .concat(timingAtropos)
  .set { timingMerged }

process ParseTiming {
    container "jdidion/python_bash"
    
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
    container "jdidion/python_bash"
    
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
  val bamFileList from sortedBams.toList()
  
  output:
  file "effectiveness.txt" into effectiveness
  
  script:
  bamFiles = bamFileList.join(" ")
  """
  compute_real_effectiveness.py \
    -i $bamFiles -o effectiveness.txt \
    --no-edit-distance --no-progress
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
  file "wgbs_effectiveness.svg"
  
  script:
  """
  show_wgbs_effectiveness.py \
    -i $effData -o wgbs_effectiveness \
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
  container "jdidion/python_bash"
    
  input:
  set val(name), file(machine) from machineMerged

  output:
  stdout machineParsed

  script:
  """
  parse_machine.py -i $machine -p $name
  """
}

process CreateMachineTable {
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  val parsedRows from machineParsed.toList()
    
  output:
  file "machine_info.txt"
  
  script:
  data = parsedRows.join("")
  """
  echo -e "prog\tprog2\tthreads\tdataset\tqcut\tcpus\tmemory\tcpu_details" > machine_info.txt
  echo '$data' >> machine_info.txt
  """
}