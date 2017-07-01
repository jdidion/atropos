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

// An absolute path to the container image is required for Singularity but
// not Docker
params.containerPrefix = ""
params.containerSuffix = ""

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
  container {
    "${params.containerPrefix}jdidion/atropos_wgbs${params.containerSuffix}"
  }
  
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
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/atropos_paper${params.containerSuffix}"
  }

  input:
  set val(_ignore_), file(reads) from atroposWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  output:
  val(taskId) into atroposTaskId
  set val("${taskId}"), file("${taskId}.{1,2}.fq.gz") into trimmedAtropos
  set val("${taskId}"), file("${taskId}.timing.txt") into timingAtropos
  set val("${taskId}"), val("trim"), file("${taskId}.machine_info.txt") into machineAtropos
  set val("${taskId}"), file(".jobid") into jobAtropos
  file "${taskId}.report.txt"
  
  script:
  taskId = "atropos_${task.cpus}_wgbs_q${qcut}_${aligner}_${compression}"
  mergeCmd = ''
  if (compression == 'nowriter') {
    mergeCmd = """
    zcat ${taskId}.1.*.fq.gz | gzip > ${taskId}.1.fq.gz
    zcat ${taskId}.2.*.fq.gz | gzip > ${taskId}.2.fq.gz
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
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt atropos \
    --op-order GACQW -T $task.cpus --batch-size $params.batchSize \
    -e 0.20 --aligner $aligner $alignerArgs -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2  -m $params.minLength \
    --no-default-adapters --no-cache-adapters --log-level ERROR \
    --correct-mismatches liberal \
    -o ${taskId}.1.fq.gz -p ${taskId}.2.fq.gz \
    --report-file ${taskId}.report.txt --quiet \
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
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/skewer${params.containerSuffix}"
  }

  input:
  set val(_ignore_), file(reads) from skewerWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into skewerTaskId
  set val("${taskId}"), file("${taskId}.{1,2}.fq.gz") into trimmedSkewer
  set val("${taskId}"), file("${taskId}.timing.txt") into timingSkewer
  set val("${taskId}"), val("trim"), file("${taskId}.machine_info.txt") into machineSkewer
  set val("${taskId}"), file(".jobid") into jobASkewer
  file "${taskId}.report.txt"
  
  script:
  taskId = "skewer_${task.cpus}_wgbs_q${qcut}"
  """
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt skewer \
    -m pe -l $params.minLength -r 0.3 \
    -o $task.tag -z --quiet \
    -x $params.adapter1 -y $params.adapter2 -t $task.cpus \
    -q $qcut -n ${reads[0]} ${reads[1]} \
    > ${taskId}.report.txt && \
  mv ${taskId}-trimmed-pair1.fastq.gz ${taskId}.1.fq.gz && \
  mv ${taskId}-trimmed-pair2.fastq.gz ${taskId}.2.fq.gz
  """
}

/* Process: SeqPurge adapter trimming
 * ----------------------------------
 */
process SeqPurge {
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/seqpurge${params.containerSuffix}"
  }

  input:
  set val(_ignore_), file(reads) from seqPurgeWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into seqPurgeTaskId
  set val("${taskId}"), file("${taskId}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${taskId}"), file("${taskId}.timing.txt") into timingSeqPurge
  set val("${taskId}"), val("trim"), file("${taskId}.machine_info.txt") into machineSeqPurge
  set val("${taskId}"), file(".jobid") into jobSeqPurge
  file "${taskId}.report.txt"
  
  script:
  taskId = "seqpurge_${task.cpus}_wgbs_q${qcut}"
  """
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 ${taskId}.1.fq.gz -out2 ${taskId}.2.fq.gz \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -threads $task.cpus -prefetch $params.batchSize \
    -ec -match_perc 70 -summary ${taskId}.report.txt
  """
}

/* Process: AdapterRemoval2 adapter trimming
 * -----------------------------------------
 */
process AdapterRemoval {
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/adapterremoval${params.containerSuffix}"
  }

  input:
  set val(_ignore_), file(reads) from adapterRemovalWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into adapterRemovalTaskId
  set val("${taskId}"), file("${taskId}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${taskId}"), file("${taskId}.timing.txt") into timingAdapterRemoval
  set val("${taskId}"), val("trim"), file("${taskId}.machine_info.txt") into machineAdapterRemoval
  set val("${taskId}"), file(".jobid") into jobAdapterRemoval
  file "${taskId}.report.txt"
  
  script:
  taskId = "adapterremoval_${task.cpus}_wgbs_q${qcut}"
  """
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${taskId}.1.fq.gz --output2 ${taskId}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut --mm 0.3 \
    --minlength $params.minLength --threads $task.cpus
    > "${taskId}.report.txt"
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
  tag { taskId }
  cpus { params.alignThreads }
  container {
    "${params.containerPrefix}jdidion/bwameth_hg38index${params.containerSuffix}"
  }
  
  input:
  set val(name), file(fastq) from trimmedMerged
  
  output:
  val(taskId) into bwamethTaskId
  file("${name}.sam")
  file("${name}.name_sorted.bam") into sortedBams
  set val(name), file("${taskId}.timing.txt") into timingBwameth
  set val("${name}"), val("bwameth"), file("${taskId}.machine_info.txt") into machineBwameth
  
  script:
  taskId = "${name}.bwameth"
  """
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt bwameth.py \
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
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  
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
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  memory "32 GB"
  
  input:
  file effData from effectiveness
  file toolNamesFile from toolNames
  
  output:
  file "wgbs_effectiveness.svg"
  
  script:
  """
  show_wgbs_effectiveness.py \
    -i $effData -o wgbs_effectiveness \
    -t $toolNamesFile --exclude-discarded -f svg
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
    container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
    
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
    container {
      "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
    }
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
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  
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
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
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
  echo '$data' | show_performance.py -t $params.alignThreads \
    -o bwameth_performance -f tex svg pickle
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

/* Channel: merged job outputs
 * ---------------------------
 * If you add a tool process, make sure to add the trimmedXXX channel
 * into the concat list here.
 */
Channel
  .empty()
  .concat(
    jobAtropos,
    jobSkewer,
    jobSeqPurge,
    jobAdapterRemoval
  )
  .set { jobMerged }

/* Process: summarize job info
 * ---------------------------
 * Generates a summary of each job using a job scheduler-specific
 * script. We have to run this outside the container to access the job 
 * scheduler.
 */
process SummarizeJob {
  input:
  set val(name), file(jobid) from jobMerged
  
  when:
  workflow.profile == 'cluster'
  
  output:
  set val(name), file("${name}.jobSummary.txt") into jobSummary

  script:
  """
  $params.summarizeJobScript $jobid > ${name}.jobSummary.txt
  """
}

/* Process: parse job info
 * -----------------------
 * Parses each job summary using a job scheduler-specific script.
 */
process ParseJob {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  
  input:
  set val(name), file(jobSummaryFile) from jobSummary
  
  output:
  stdout jobParsed
  
  script:
  """
  $params.parseJobScript -i $jobSummaryFile -p $name
  """
}

/* Process: show job memory usage
 * ------------------------------
 * Generates memory usage table/figure from job summaries.
 */
process ShowJobMemoryUsage {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  val parsedJobs from jobParsed.toList()
  
  output:
  file "job.mem.tex"
  file "job.mem.pickle"
  file "job.mem.svg"
  
  script:
  data = parsedJobs.join("")
  """
  echo '$data' | show_job_info.py -o job \
    -m mem=${params.jobMemoryMetric} -f tex pickle svg
  """
}
