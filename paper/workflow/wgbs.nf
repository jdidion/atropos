#!/usr/bin/env nextflow

/* Atropos paper workflow for real WGBS reads
 * ------------------------------------------
 * Configuration is externalized into separate profiles for local and cluster.
 * The main difference is that Docker is used to manage images on desktop, 
 * verusus Singularity on the cluster.
 *
 * Parameters
 * ----------
 * The following are expected to be defined in params:
 * - threadCounts: The numbers of threads (cpus) to test
 * - minLength: minimum read length
 * - batchSize: read batch size
 * - quals: quality thresholds for trimming
 * - aligners: Atropos aligner algorithms to test
 * - compressionSchemes: Atropos compression schemes to test
 * - adapter1, adapter2: Adapter sequence
 * 
 * Note: The use of /proc/cpuinfo and /proc/meminfo is linux-specific, which
 * is fine when running in Docker/Singularity. If you disable containers,
 * you'll need to change the command.
 * 
 * Note: If Nextflow exits with error 137, this is due to insufficient
 * memory. There are a few things you can try to work around this:
 * - Set executor.cpus to the max value in params.threadCounts, so only one
 *   process will run at a time
 * - Reduce params.batchSize
 * - If you have all of the software installed locally, you can disable
 *   Docker/Singularity by commenting out the 'container' directives.e.
 */

// variables for all tools
params.publishDir = "${params.resultsDir}/${workflow.profile}/wgbs"
params.quals = [ 0, 20 ]
params.adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000

// atropos-specific variables
params.aligners = [ 'insert' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

/* Process: Extract data from container
 * ------------------------------------
 * The read data is contained within a Docker container. Since Singularity does
 * not support attaching volumes, we use a process to unpack the data from the 
 * container.
 */
process ExtractReads {
  container "jdidion/atropos_wgbs"
  
  output:
  set val("untrimmed"), file("wgbs.{1,2}.fq.gz") into wgbsReads
  
  script:
  """
  cp \
    /data/wgbs/wgbs.1.fq.gz \
    /data/wgbs/wgbs.2.fq.gz \
    .
  """
}

// Split the input reads channel into one per tool process.
wgbsReads.into {
  untrimmedWgbsReads
  atroposWgbsReads
  skewerWgbsReads
  seqPurgeWgbsReads
  adapterRemovalWgbsReads
}

/* Process: Atropos adapter trimming
 * ---------------------------------
 */
process Atropos {
  tag { "atropos_${task.cpus}_wgbs_q${qcut}_${aligner}_${compression}" }
  cpus { threads }
  container "jdidion/atropos_paper"

  input:
  set val(_ignore_), file(reads) from atroposWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  set val("${task.tag}"), val("trim"), file("${task.tag}.machine_info.txt") into machineAtropos
  file "${task.tag}.report.txt"
  
  script:
  mergeCmd = ''
  if (compression == 'nowriter') {
    mergeCmd = """
    zcat ${task.tag}.1.*.fq.gz | gzip > ${task.tag}.1.fq.gz
    zcat ${task.tag}.2.*.fq.gz | gzip > ${task.tag}.2.fq.gz
    """
  }
  // Tune different parameters based on aligner. For adapter-match, minimum
  // overlap is important to prevent over-trimming. For insert-match, the
  // error rate is important to allow for mismatches between overlaps.
  alignerArgs = null
  if (aligner == 'insert') {
    alignerArgs = "--insert-match-error-rate 0.30"
  }
  else {
    alignerArgs = "-O 7"
  }
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    --op-order GACQW -T $task.cpus --batch-size $params.batchSize \
    -e 0.20 --aligner $aligner $alignerArgs -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2  -m $params.minLength \
    --no-default-adapters --no-cache-adapters --log-level ERROR \
    --correct-mismatches liberal \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz \
    --report-file ${task.tag}.report.txt --quiet \
    $task.ext.compressionArg -pe1 ${reads[0]} -pe2 ${reads[1]}
  $mergeCmd
  """
}

/* Process: Skewer adapter trimming
 * ---------------------------------
 * Skewer forces a specific naming convention, so we rename them to be 
 * consistent with the other tools.
 */
process Skewer {
  tag { "skewer_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/skewer"

  input:
  set val(_ignore_), file(reads) from skewerWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedSkewer
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSkewer
  set val("${task.tag}"), val("trim"), file("${task.tag}.machine_info.txt") into machineSkewer
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt skewer \
    -m pe -l $params.minLength -r 0.3 \
    -o $task.tag -z --quiet \
    -x $params.adapter1 -y $params.adapter2 -t $task.cpus \
    -q $qcut -n ${reads[0]} ${reads[1]} \
    > ${task.tag}.report.txt && \
  mv ${task.tag}-trimmed-pair1.fastq.gz ${task.tag}.1.fq.gz && \
  mv ${task.tag}-trimmed-pair2.fastq.gz ${task.tag}.2.fq.gz
  """
}

/* Process: SeqPurge adapter trimming
 * ----------------------------------
 */
process SeqPurge {
  tag { "seqpurge_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/seqpurge"

  input:
  set val(_ignore_), file(reads) from seqPurgeWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSeqPurge
  set val("${task.tag}"), val("trim"), file("${task.tag}.machine_info.txt") into machineSeqPurge
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 ${task.tag}.1.fq.gz -out2 ${task.tag}.2.fq.gz \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -threads $task.cpus -prefetch $params.batchSize \
    -ec -match_perc 70 -summary ${task.tag}.report.txt
  """
}

/* Process: AdapterRemoval2 adapter trimming
 * -----------------------------------------
 */
process AdapterRemoval {
  tag { "adapterremoval_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/adapterremoval"

  input:
  set val(_ignore_), file(reads) from adapterRemovalWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAdapterRemoval
  set val("${task.tag}"), val("trim"), file("${task.tag}.machine_info.txt") into machineAdapterRemoval
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${task.tag}.1.fq.gz --output2 ${task.tag}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut --mm 0.3 \
    --minlength $params.minLength --threads $task.cpus
    > "${task.tag}.report.txt"
  """
}

/* Channel: merged tool outputs
 * ----------------------------
 * We also add in the untrimmed output, which we need to align for evaluation of
 * effectiveness. If you add a tool process, make sure to add the trimmedXXX 
 * channel into the concat list here.
 */
Channel
  .empty()
  .concat(
    untrimmedWgbsReads,
    trimmedAtropos,
    trimmedSkewer,
    trimmedSeqPurge,
    trimmedAdapterRemoval
  )
  .set { trimmedMerged }

/* Process: Align reads using bwa-meth
 * -----------------------------------
 * Alignments are written to an unsorted BAM file and then name sorted by
 * samtools.
 */
process BwamethAlign {
  tag { "${name}.bwameth" }
  cpus { params.alignThreads }
  container "jdidion/bwameth_hg38index"
  
  input:
  set val(name), file(fastq) from trimmedMerged
  
  output:
  file("${name}.sam")
  file("${name}.name_sorted.bam") into sortedBams
  set val(name), file("${task.tag}.timing.txt") into timingBwameth
  set val("${name}"), val("bwameth"), file("${task.tag}.machine_info.txt") into machineBwameth
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt bwameth.py \
    -t ${params.alignThreads} --read-group '${task.ext.readGroup}' \
    --reference /data/index/bwameth/hg38/hg38.fa \
    ${fastq[0]} ${fastq[1]} > ${name}.sam \
  && samtools view -Shb ${name}.sam | \
     samtools sort -n -O bam -@ ${params.alignThreads} \
       -o ${name}.name_sorted.bam -
  """
}

/* Channel: display names for tools
 * --------------------------------
 */
toolNames = Channel.fromPath(
  "${workflow.projectDir}/../containers/tools/tool-names.txt")

/* Process: Summarize trimming effectiveness
 * -----------------------------------------
 */
process ComputeEffectiveness {
  container "jdidion/python_bash"
  
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
  container "jdidion/python_bash"
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  file effData from effectiveness
  file toolNamesFile from toolNames
  
  output:
  file "wgbs_effectiveness.svg"
  file "wgbs_effectiveness.pickle"
  
  script:
  """
  show_wgbs_effectiveness.py \
    -i $effData -o wgbs_effectiveness \
    -t $toolNamesFile --exclude-discarded -f svg pickle
  """
}

/* Channel: merged performance outputs
 * -----------------------------------
 * If you add a tool process, make sure to add the stdout channel
 * into the concat list here.
 */
Channel
  .empty()
  .concat(
    timingAtropos,
    timingSkewer,
    timingSeqPurge,
    timingAdapterRemoval
  )
  .set { timingMerged }

/* Process: parse performance data for each tool
 * ---------------------------------------------
 */
process ParseTrimmingTiming {
    container "jdidion/python_bash"
    
    input:
    set val(name), file(timing) from timingMerged
    
    output:
    stdout timingMergedParsed
    
    script:
    """
    parse_gtime.py -i $timing -p $name
    """
}

/* Process: generate performance figure/table
 * ------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowTrimmingPerformance {
    container "jdidion/python_bash"
    publishDir "$params.publishDir", mode: 'copy', overwrite: true
    
    input:
    val parsedRows from timingMergedParsed.toList()
    
    output:
    file "trim_performance.tex"
    file "trim_performance.svg"
    file "trim_performance.pickle"
    
    script:
    data = parsedRows.join("")
    """
    echo '$data' | show_performance.py -o trim_performance -f tex svg pickle
    """
}

/* Process: parse performance for bwa-meth alignment
 * -------------------------------------------------
 */
process ParseBwamethTiming {
    container "jdidion/python_bash"
    
    input:
    set val(name), file(timing) from timingBwameth
    
    output:
    stdout timingBwamethParsed
    
    script:
    """
    parse_gtime.py -i $timing -p $name
    """
}

/* Process: generate bwa-meth alignment performance figure/table
 * -------------------------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowBwamethPerformance {
    container "jdidion/python_bash"
    publishDir "$params.publishDir", mode: 'copy', overwrite: true
    
    input:
    val parsedRows from timingBwamethParsed.toList()
    
    output:
    file "bwameth_performance.tex"
    file "bwameth_performance.svg"
    file "bwameth_performance.pickle"
    
    script:
    data = parsedRows.join("")
    """
    echo '$data' | show_performance.py -o bwameth_performance -f tex svg pickle
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
    machineSkewer,
    machineSeqPurge,
    machineAdapterRemoval,
    machineBwameth
  )
  .set { machineMerged }

/* Process: summarize machine info
 * -------------------------------
 */
process SummarizeMachine {
  container "jdidion/python_bash"
    
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
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
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