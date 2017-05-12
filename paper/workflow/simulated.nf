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
 * The figure names and captions are also defined in the config file.
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
process Extract {
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
  set val("${task.tag}"), err, alns, file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    -T $task.cpus --aligner $aligner --op-order GACQW \
    -a $params.adapter1 -A $params.adapter2 -q $qcut --trim-n \
    -m $params.minLength --batch-size $params.batchSize \
    --no-default-adapters --no-cache-adapters --log-level ERROR --quiet \
    --insert-match-error-rate 0.20 -e 0.10 \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz  \
    --report-file ${task.tag}.report.txt --quiet \
    $task.ext.compressionArg -pe1 ${reads[0]} -pe2 ${reads[1]}
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
  each err from params.errorRates
  each qcut from params.quals

  output:
  set val("${task.tag}"), err, alns, file("${task.tag}.{1,2}.fq.gz") into trimmedSkewer
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSkewer
  file "${task.tag}.report.txt"
  
  script:
  """
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
  each err from params.errorRates
  each qcut from params.quals

  output:
  set val("${task.tag}"), err, alns, file("${task.tag}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSeqPurge
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 ${task.tag}-trimmed-pair1.pair1.fastq.gz \
    -out2 ${task.tag}-trimmed-pair1.pair2.fastq.gz \
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
  each err from params.errorRates
  each qcut from params.quals

  output:
  set val("${task.tag}"), err, alns, file("${task.tag}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAdapterRemoval
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${task.tag}.1.fq.gz --output2 ${task.tag}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut \
    --minlength $params.minLength --threads $task.cpus
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
    --name ${name} --no-progress \
    -o ${name}.txt -s ${name}.summary.txt -t -
  """
}

/* Process: generate accuracy table
 * --------------------------------
 * Aggregate all the accuracy rows from each tool and pass it to stdin of
 * bin/show_simulated_accuracy.py.
 */
process ShowSimulatedAccuracy {
  container "jdidion/python_bash"
  publishDir "$publishDir", mode: 'copy', overwrite: true
  
  input:
  val accuracyRows from tableFile.toList()
  
  output:
  file "accuracy.tex"
  
  script:
  data = accuracyRows.join("")
  """
  echo '$data' | show_simulated_accuracy.py -n $task.ext.name -c $task.ext.caption -o accuracy
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
    publishDir "$publishDir", mode: 'copy', overwrite: true
    
    input:
    val parsedRows from timingParsed.toList()
    
    output:
    file "performance.tex"
    file "performance.svg"
    
    script:
    data = parsedRows.join("")
    if (workflow.profile == "local") {
      name = task.ext.local_name
      caption = task.ext.local_caption
    } else {
      name = task.ext.cluster_name
      caption = task.ext.cluster_caption
    }
    
    """
    echo '$data' | show_performance.py -n $name -c $caption -o performance
    """
}
