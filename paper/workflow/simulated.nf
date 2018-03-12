#!/usr/bin/env nextflow

/* Atropos paper workflow for simulated DNA-Seq reads
 * --------------------------------------------------
 * Configuration is externalized into separate profiles for local and cluster.
 * The main difference is that Docker is used to manage images on desktop, 
 * verusus Singularity on the cluster.
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
 * 
 * Note: If Nextflow exits with error 137, this is due to insufficient
 * memory. There are a few things you can try to work around this:
 * - Set executor.cpus to the max value in params.threadCounts, so only one
 *   process will run at a time
 * - Reduce params.batchSize
 * - If you have all of the software installed locally, you can disable
 *   Docker/Singularity by commenting out the 'container' directives.
 *
 * Note: When running the 'local' profile, a fake .jobid file is generated in
 * order to satisfy the output file requirement, even though the tasks that
 * process .jobid files are not run. Hopefully a future version of nextflow
 * will support conditional output files.
 */

// An absolute path to the container image is required for Singularity but
// not Docker
params.containerPrefix = ""
params.containerSuffix = ""

// variables for all tools
params.publishDir = "${params.resultsDir}/${workflow.profile}/simulated"
params.errorRates = [ '001', '005', '01' ]
params.quals = [ 0 ]
params.adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000

// atropos-specific variables
params.aligners = [ 'insert', 'adapter' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

/* Process: Extract data from container
 * ------------------------------------
 * The read data is contained within a Docker container. Since Singularity does
 * not support attaching volumes, we use a process to unpack the data from the 
 * container.
 */
process ExtractReads {
  container { 
    "${params.containerPrefix}jdidion/atropos_simulated${params.containerSuffix}"
  }
  
  input:
  each err from params.errorRates
  
  output:
  set err, file("sim_${err}.{1,2}.fq"), file("sim_${err}.{1,2}.aln") into simReads
  
  script:
  """
  cp \
    /data/simulated/sim_${err}.1.fq  \
    /data/simulated/sim_${err}.2.fq  \
    /data/simulated/sim_${err}.1.aln \
    /data/simulated/sim_${err}.2.aln \
    .
  """
}

// Split the input reads channel into one per tool process.
simReads.into {
  atroposSimReads
  skewerSimReads
  seqPurgeSimReads
  adapterRemovalSimReads
  fastpSimReads
  cutadaptSimReads
}

/* Process: Atropos adapter trimming
 * ---------------------------------
 * Post-processing is done when aligner == 'nowriter' to concatenate all the
 * output files for each pair because the process is expecting a single pair of 
 * output files.
 */
process Atropos {
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/atropos_paper${params.containerSuffix}"
  }

  input:
  set err, file(reads), file(alns) from atroposSimReads
  each threads from params.threadCounts
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  output:
  val(taskId) into atroposTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedAtropos
  set val("${taskId}"), file("${taskId}.timing.txt") into timingAtropos
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineAtropos
  set val("${taskId}"), file(".jobid") into jobAtropos
  file "${taskId}.report.txt"
  
  script:
  taskId = "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}"
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
    alignerArgs = "--insert-match-error-rate 0.20"
  }
  else {
    alignerArgs = "-O 7"
  }
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt atropos trim \
    --op-order GACQW -T $task.cpus --batch-size $params.batchSize \
    -e 0.10 --aligner $aligner $alignerArgs -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2 -m $params.minLength \
    --no-default-adapters --no-cache-adapters --log-level ERROR \
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
  set err, file(reads), file(alns) from skewerSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into skewerTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedSkewer
  set val("${taskId}"), file("${taskId}.timing.txt") into timingSkewer
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineSkewer
  set val("${taskId}"), file(".jobid") into jobSkewer
  file "${taskId}.report.txt"
  
  script:
  taskId = "skewer_${task.cpus}_${err}_q${qcut}"
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt skewer \
    -m pe -l $params.minLength -r 0.2 \
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
  set err, file(reads), file(alns) from seqPurgeSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into seqPurgeTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${taskId}"), file("${taskId}.timing.txt") into timingSeqPurge
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineSeqPurge
  set val("${taskId}"), file(".jobid") into jobSeqPurge
  file "${taskId}.report.txt"
  
  script:
  taskId = "seqpurge_${task.cpus}_${err}_q${qcut}"
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 "${taskId}.1.fq.gz" -out2 "${taskId}.2.fq.gz" \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -threads $task.cpus -prefetch $params.batchSize \
    -match_perc 80 -summary ${taskId}.report.txt
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
  set err, file(reads), file(alns) from adapterRemovalSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into adapterRemovalTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${taskId}"), file("${taskId}.timing.txt") into timingAdapterRemoval
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineAdapterRemoval
  set val("${taskId}"), file(".jobid") into jobAdapterRemoval
  file "${taskId}.report.txt"
  
  script:
  taskId = "adapterremoval_${task.cpus}_${err}_q${qcut}"
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${taskId}.1.fq.gz --output2 ${taskId}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut \
    --minlength $params.minLength --threads $task.cpus \
    > "${taskId}.report.txt"
  """
}

/* Process: fastp adapter trimming
 * -------------------------------
 */
process fastp {
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/fastp${params.containerSuffix}"
  }

  input:
  set err, file(reads), file(alns) from fastpSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into fastpTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedFastp
  set val("${taskId}"), file("${taskId}.timing.txt") into timingFastp
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineFastp
  set val("${taskId}"), file(".jobid") into jobFastp
  file "${taskId}.report.txt"

  script:
  taskId = "fastp_${task.cpus}_${err}_q${qcut}"
  cutArgs = ""
  if (qcut > 0) {
    cutArgs = "--cut_by_quality3 --cut_mean_quality $qcut"
  }
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt fastp \
    -i ${reads[0]} -I ${reads[1]} -o ${taskId}.1.fq.gz -O ${taskId}.2.fq.gz \
    --adapter_sequence $params.adapter1 --adapter_sequence_r2 $params.adapter2 \
    --thread $task.cpus $cutArgs \
    --length_required $params.minLength --disable_quality_filtering \
    > "${taskId}.report.txt"
  """
}

/* Process: Cutadapt adapter trimming
 * ----------------------------------
 */
process Cutadapt {
  tag { taskId }
  cpus { threads }
  container {
    "${params.containerPrefix}jdidion/cutadapt${params.containerSuffix}"
  }

  input:
  set err, file(reads), file(alns) from cutadaptSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  val(taskId) into cutadaptTaskId
  set val("${taskId}"), val(err), file(alns), file("${taskId}.{1,2}.fq.gz") into trimmedCutadapt
  set val("${taskId}"), file("${taskId}.timing.txt") into timingCutadapt
  set val("${taskId}"), file("${taskId}.machine_info.txt") into machineCutadapt
  set val("${taskId}"), file(".jobid") into jobCutadapt
  file "${taskId}.report.txt"

  script:
  taskId = "cutadapt_${task.cpus}_${err}_q${qcut}"
  """
  if [ "${workflow.profile}" == "local" ]; then
    touch .jobid
  fi
  cat /proc/cpuinfo /proc/meminfo > ${taskId}.machine_info.txt
  /usr/bin/time -v -o ${taskId}.timing.txt cutadapt \
    -j $task.cpus -O 7 -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2 -m $params.minLength \
    -o ${taskId}.1.fq.gz -p ${taskId}.2.fq.gz ${reads[0]} ${reads[1]} \
    > "${taskId}.report.txt"
  """
}

/* Channel: merged tool outputs
 * ----------------------------
 * If you add a tool process, make sure to add the trimmedXXX channel
 * into the concat list here.
 */
Channel
  .empty()
  .concat(
    trimmedAtropos,
    trimmedSkewer,
    trimmedSeqPurge,
    trimmedAdapterRemoval,
    trimmedFastp,
    trimmedCutadapt
  )
  .set { trimmedMerged }

/* Process: compute accuracy for each tool
 * ---------------------------------------
 */
process ComputeSimulatedAccuracy {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  
  input:
  set val(name), val(err), file(alns), file(trimmed) from trimmedMerged
  
  output:
  file "${name}.txt"
  file "${name}.summary.txt"
  stdout tableFile

  script:
  """
  compute_simulated_accuracy.py \
    -a1 ${alns[0]} -a2 ${alns[1]} \
    -r1 ${trimmed[0]} -r2 ${trimmed[1]} \
    -p ${name} --no-progress \
    -o ${name}.txt -s ${name}.summary.txt -t -
  """
}

/* Channel: display names for error rates
 * --------------------------------------
 */
errorRates = Channel.fromPath(
  "${workflow.projectDir}/../containers/data/simulated/art_profiles.txt")

/* Channel: display names for tools
 * --------------------------------
 */
toolNames = Channel.fromPath(
  "${workflow.projectDir}/../containers/tools/tool-names.txt")

/* Process: generate accuracy table
 * --------------------------------
 * Aggregate all the accuracy rows from each tool and pass it to stdin of
 * bin/show_simulated_accuracy.py.
 */
process ShowSimulatedAccuracy {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  val accuracyRows from tableFile.toList()
  file errorRatesFile from errorRates
  file toolNamesFile from toolNames
  
  output:
  file "accuracy.tex"
  file "accuracy.pickle"
  
  script:
  data = accuracyRows.join("")
  """
  echo '$data' | show_simulated_accuracy.py -o accuracy -f txt tex pickle \
    -e $errorRatesFile -t $toolNamesFile
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
    timingAdapterRemoval,
    timingFastp,
    timingCutadapt
  )
  .set { timingMerged }

/* Process: parse performance data for each tool
 * ---------------------------------------------
 */
process ParseSimualtedTiming {
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

/* Process: generate performance figure/table
 * ------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowSimulatedPerformance {
  container {
    "${params.containerPrefix}jdidion/python_bash${params.containerSuffix}"
  }
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  val parsedRows from timingParsed.toList()
  
  output:
  file "performance.tex"
  file "performance.svg"
  file "performance.pickle"
  
  script:
  data = parsedRows.join("")
  """
  echo '$data' | show_performance.py -o performance -f tex svg pickle
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
    machineFastp,
    machineCutadapt
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
  set val(name), file(machine) from machineMerged

  output:
  stdout machineParsed

  script:
  """
  parse_machine.py -i $machine -p $name
  """
}

/* Process: print machine info table
 * ---------------------------------
 */
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
    jobAdapterRemoval,
    jobFastp,
    jobCutadapt
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

  when:
  workflow.profile == 'cluster'

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

  when:
  workflow.profile == 'cluster'

  output:
  file "job.mem.pickle"
  file "job.mem.tex"
  file "job.mem.svg"
  
  script:
  data = parsedJobs.join("")
  """
  echo '$data' | show_job_info.py -o job \
    -m mem=${params.jobMemoryMetric} -f pickle tex svg
  """
}
