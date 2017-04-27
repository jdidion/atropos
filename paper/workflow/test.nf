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

process Extract {
  container true
  
  input:
  each err from params.errorRates
  
  output:
  set err, "sim_${err}.{1,2}.fq", "sim_${err}.{1,2}.aln" into simReads
  
  script:
  """
  jdidion/atropos_simulated cp \
    /data/simulated/sim_${err}.1.fq  \
    /data/simulated/sim_${err}.2.fq  \
    /data/simulated/sim_${err}.1.aln \
    /data/simulated/sim_${err}.2.aln \
    .
  """
}

process Atropos {
  //tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  tag { "atropos_${task.cpus}_${err}" }
  cpus { threads }
  container "jdidion/atropos_paper"
  
  input:
  set err, file(reads), file(alns) from simReads
  each threads from params.threadCounts
  //each qcut from params.quals
  //each aligner from params.aligners
  //each compression from params.compressionSchemes
  
  output:
  set err, "${task.tag}.{1,2}.fq.gz" into trimmedReads
  file "${task.tag}.report.txt" into reportFile
  file "${task.tag}.timing.txt" into timingAtropos

  script:
  """
  >&2 echo ${task.tag} && \
  /usr/bin/time -v atropos \
    -pe1 ${reads[0]} -pe2 ${reads[1]} \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz \
    --report-file ${task.tag}.report.txt --quiet \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    2> ${task.tag}.timing.txt
  """
}

// concatenate all of the timing results into a single channel
Channel
  .empty()
  .concat(timingAtropos)
  .set { timing }

process Timing {
  input:
  val timingFiles from timing.toList()

  output:
  file 'timing.txt'

  script:
  """
  cat $timingFiles > timing.txt
  """
}
