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
 */

// variables for all tools
params.publishDir = "results/${workflow.profile}/simulated"
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
  container "jdidion/atropos_simulated"
  
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
}

/* Process: Atropos adapter trimming
 * ---------------------------------
 * Post-processing is done when aligner == 'nowriter' to concatenate all the
 * output files for each pair because the process is expecting a single pair of 
 * output files.
 */
process Atropos {
  tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  cpus { threads }
  container "jdidion/atropos_paper"

  input:
  set err, file(reads), file(alns) from atroposSimReads
  each threads from params.threadCounts
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), val(err), file(alns), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineAtropos
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
    alignerArgs = "--insert-match-error-rate 0.20"
  }
  else {
    alignerArgs = "-O 7"
  }
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos trim \
    --op-order GACQW -T $task.cpus --batch-size $params.batchSize \
    -e 0.10 --aligner $aligner $alignerArgs -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2 -m $params.minLength \
    --no-default-adapters --no-cache-adapters --log-level ERROR \
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
  tag { "skewer_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }
  container "jdidion/skewer"

  input:
  set err, file(reads), file(alns) from skewerSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), val(err), file(alns), file("${task.tag}.{1,2}.fq.gz") into trimmedSkewer
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSkewer
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineSkewer
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt skewer \
    -m pe -l $params.minLength -r 0.2 \
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
  tag { "seqpurge_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }
  container "jdidion/seqpurge"

  input:
  set err, file(reads), file(alns) from seqPurgeSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), val(err), file(alns), file("${task.tag}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSeqPurge
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineSeqPurge
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 "${task.tag}.1.fq.gz" -out2 "${task.tag}.2.fq.gz" \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -threads $task.cpus -prefetch $params.batchSize \
    -match_perc 80 -summary ${task.tag}.report.txt
  """
}

/* Process: AdapterRemoval2 adapter trimming
 * -----------------------------------------
 */
process AdapterRemoval {
  tag { "adapterremoval_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }
  container "jdidion/adapterremoval"

  input:
  set err, file(reads), file(alns) from adapterRemovalSimReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), val(err), file(alns), file("${task.tag}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAdapterRemoval
  set val("${task.tag}"), file("${task.tag}.machine_info.txt") into machineAdapterRemoval
  file "${task.tag}.report.txt"
  
  script:
  """
  cat /proc/cpuinfo /proc/meminfo > ${task.tag}.machine_info.txt
  /usr/bin/time -v -o ${task.tag}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${task.tag}.1.fq.gz --output2 ${task.tag}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut \
    --minlength $params.minLength --threads $task.cpus \
    > "${task.tag}.report.txt"
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
    trimmedAdapterRemoval
  )
  .set { trimmedMerged }

/* Process: compute accuracy for each tool
 * ---------------------------------------
 */
process ComputeSimulatedAccuracy {
  container "jdidion/python_bash"
  
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
  container "jdidion/python_bash"
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
    timingAdapterRemoval
  )
  .set { timingMerged }

/* Process: parse performance data for each tool
 * ---------------------------------------------
 */
process ParseSimualtedTiming {
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

/* Process: generate performance figure/table
 * ------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowSimulatedPerformance {
    container "jdidion/python_bash"
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
    machineAdapterRemoval
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
  parse_machine.py -i $timing -p $name
  """
}

/* Process: print machine info table
 * ---------------------------------
 */
process CreateMachineTable {
  publishDir "$params.publishDir", mode: 'copy', overwrite: true
  
  input:
  val parsedRows from machineMerged.toList()
    
  output:
  file "machine_info.txt"
  
  script:
  data = parsedRows.join("")
  """
  echo -e "prog\tprog2\tthreads\tdataset\tqcut\tcpus\tmemory\tcpu_details" > machine_info.txt
  echo '$data' >> machine_info.txt
  """
}