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
 * - dataDir: The local directory where data from the container will be copied
 *   (for Singularity execution only)
 * - dataContainer: The name of the Docker Hub repository from which the data
 *   data container should be pulled.
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
  storeDir { dataDir }
  
  input:
  each err from params.errorRates
  each readPair from {[1, 2]}
  each fileExt from {['fq', 'aln']}
  val dataDir from params.dataDir
  
  output:
  file "sim_${err}.${readPair}.${fileExt}" into destFile
  
  script:
  """
  jdidion/atropos_simulated cp /data/simulated/sim_${err}.${readPair}.${fileExt} .
  """
}

process Atropos {
  echo true
  tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  cpus { threads }
  container "jdidion/atropos_paper"
  
  input:
  each threads from params.threadCounts
  each err from params.errorRates
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  val dataDir from params.dataDir
  file input1 from "${dataDir}/sim_${err}.1.fq"
  file input2 from "${dataDir}/sim_${err}.2.fq"
  
  script:
  """
  echo $task.tag $threads $compression $task.ext.compressionArg $input1 $input2
  """
}
