#!/usr/bin/env nextflow

/* Atropos paper workflow
 * ----------------------
 * Configuration is externalized into separate files for desktop and cluster.
 * The main difference is that Docker images are used on desktop, verusus
 * Singularity on the cluster.
 */

/* Parameters
 * ----------
 * The following are expected to be defined externally:
 * - minLen: minimum read length
 * - batchSize: read batch size
 * - quals: quality thresholds for trimming
 * - aligners: Atropos aligner algorithms to use
 * - adapter1, adapter2: Adapter sequence
 * - extra: map of extra arguments for each program
 */

process atropos {
  tag { name }

  input:
  val compressionArg { compression == "nowriter" ? "--no-writer-process" : "--compression $compression" }
  val name from "atropos_${threads}_${err}_q${qcut}_${aligner}_${compression}"
  file input1 from fq1
  file input2 from fq2

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile

  script:
  """
  /usr/bin/time -v atropos \
    -T $threads --aligner $aligner --op-order GACQW
    "-a $adapter1 -A $adapter2 -q $qcut --trim-n
    "-m $minLen --batch-size $batchSize
    "--report-file $reportFile
    "--no-default-adapters --no-cache-adapters
    "-o $output1 -p $output2
    "--log-level ERROR --quiet $extra.atropos
    "$compressionArg -pe1 $input1 -pe2 $input2
  """
}

process skewer {
  tag { name }

  input:
  val name from "skewer_${threads}_${err}_q${qcut}"
  file input1 from fq1
  file input2 from fq2

  output:

  script:
  """
  /usr/bin/time -v skewer \
    -m pe -l $MIN_LEN $extra.skewer \
    -o $name -z --quiet \
    -x $adapter1 -y $adapter2 -t $threads \
    -q $qcut -n $input1 $input2
  """
}