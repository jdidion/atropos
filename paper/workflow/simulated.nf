#!/usr/bin/env nextflow

/* Atropos paper workflow for simulated DNA-Seq reads
 * --------------------------------------------------
 * Configuration is externalized into separate files for desktop and cluster.
 * The main difference is that Docker images are used on desktop, verusus
 * Singularity on the cluster.
 */

/* Parameters
 * ----------
 * The following are expected to be defined in params:
 * - minLength: minimum read length
 * - batchSize: read batch size
 * - quals: quality thresholds for trimming
 * - aligners: Atropos aligner algorithms to use
 * - adapter1, adapter2: Adapter sequence
 * - dataDir: The local directory where data from the container will be copied
 */

// variables for all tools
params.errorRates = [ '001', '005', '01' ]
params.quals = [ 0 ]
params.adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000
params.dataDir = '.'

// atropos-specific variables
params.aligners = [ 'insert', 'adapter' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

err = Channel.from(params.errorRates)

/* The read data is contained within a Docker container 
 * (jdidion/atropos_simulated). During local execution, we attach that container
 * directly to the tool container (using --from-volumes <data container>).
 * Singularity does not allow for this, so when running on the cluster we use a 
 * separate process to unpack the data from the container into $params.storeDir.
 * Note that storeDir persists beyond the end of the script, so we delete it in 
 * the shell script that runs the workflow.
 */

process Extract {
  storeDir "$params.dataDir"
  
  when:
  singularity.enabled is true

  input:
  each err from params.errorRates
  val fqPrefix from "sim_${err}"

  output:
  file "${fqPrefix}.1.fq"
  file "${fqPrefix}.1.aln"
  file "${fqPrefix}.2.fq"
  file "${fqPrefix}.2.aln"

  script:
  """
  cp /data/simulated/${fqPrefix}* .
  """
}

process Atropos {
  tag { name }

  input:
  each err from params.errorRates
  each aligner from params.aligners
  each qcut from params.quals
  val compressionArg { compression == "nowriter" ? "--no-writer-process" : "--compression $compression" }
  val name from "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}"
  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile

  script:
  """
  /usr/bin/time -v atropos \
    -T $task.cpus --aligner $aligner --op-order GACQW \
    "-a $params.adapter1 -A $params.adapter2 -q $qcut --trim-n \
    "-m $params.minLength --batch-size $params.batchSize \
    "--report-file $reportFile \
    "--no-default-adapters --no-cache-adapters \
    "-o $output1 -p $output2 \
    "--log-level ERROR --quiet \
    --insert-match-error-rate 0.20 -e 0.10 \
    "$compressionArg -pe1 $input1 -pe2 $input2
  """
}

process Skewer {
  tag { name }

  input:
  each err from params.errorRates
  each qcut from params.quals
  val name from "skewer_${task.cpus}_${err}_q${qcut}"
  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}-trimmed-pair1.fastq.gz" into output1
  file "${name}-trimmed-pair2.fastq.gz" into output2
  
  script:
  """
  /usr/bin/time -v skewer \
    -m pe -l $params.minLength -match_perc 80 \
    -o $name -z --quiet \
    -x $params.adapter1 -y $params.adapter2 -t $task.cpus \
    -q $qcut -n $input1 $input2
  """

  output1.onComplete: { it.renameTo("${name}.1.fq.gz") }
  output2.onComplete: { it.renameTo("${name}.2.fq.gz") }
}

process SeqPurge {
  tag { name }

  input:
  each err from params.errorRates
  each qcut from params.quals
  val name from "seqpurge_${task.cpus}_${err}_q${qcut}"
  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile

  script:
  """
  /usr/bin/time -v SeqPurge \
    -in1 $input1 -in2 $input2 \
    -out1 $output1 -out2 $output2 \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -task.cpus $task.cpus -prefetch $params.batchSize \
    -r 0.20 -summary $reportFile
  """
}

process AdapterRemoval {
  tag { name }

  input:
  each err from params.errorRates
  each qcut from params.quals
  val name from "adapterremoval_${task.cpus}_${err}_q${qcut}"
  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile

  script:
  """
  /usr/bin/time -v AdapterRemoval \
    --file1 $input1 --file2 $input2 \
    --output1 $output1 --output2 $output2 --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut \
    --minLengthgth $params.minLength --task.cpus $task.cpus
  """
}
