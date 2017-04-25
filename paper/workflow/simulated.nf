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
params.dataContainer = "jdidion/atropos_simulated"

// atropos-specific variables
params.aligners = [ 'insert', 'adapter' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

process Extract {
  storeDir { params.dataDir }
  
  when:
  workflow.profile == 'cluster'

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

// channel for output summaries
results = Channel.create()
// channel for outputs from /usr/bin/time
perfs = Channel.create()

process Atropos {
  tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  cpus { threads }

  input:
  each threads from params.threadCounts
  each err from params.errorRates
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  file input1 from "${params.dataDir}/${task.ext.fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${task.ext.fqPrefix}.2.fq"
  
  output:
  file "${task.tag}.1.fq.gz" into output1
  file "${task.tag}.2.fq.gz" into output2
  file "${task.tag}.report.txt" into reportFile
  stdout timing

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
    "$task.ext.compressionArg -pe1 $input1 -pe2 $input2
  """
  
  timing.onComplete {
    perfs << [ 
      "name" : name,
      "value" : it
    ]
  }

  afterScript {
    results << [ 
      "name" : name,
      "trimmed" : [ output1, output2 ]
    ]
  }
}

process Skewer {
  tag { "skewer_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }

  input:
  each threads from params.threadCounts
  each err from params.errorRates
  each qcut from params.quals

  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"
  stdout timing

  output:
  file "${name}-trimmed-pair1.fastq.gz" into output1
  file "${name}-trimmed-pair2.fastq.gz" into output2
  
  script:
  """
  /usr/bin/time -v skewer \
    -m pe -l $params.minLength -match_perc 80 \
    -o $task.tag -z --quiet \
    -x $params.adapter1 -y $params.adapter2 -t $task.cpus \
    -q $qcut -n $input1 $input2
  """

  timing.onComplete {
    perfs << [ 
      "name" : name,
      "value" : it
    ]
  }

  afterScript {
    output1.renameTo("${name}.1.fq.gz")
    output2.renameTo("${name}.2.fq.gz")
    results << [ 
      "name" : name,
      "trimmed" : [ "${name}.1.fq.gz", "${name}.2.fq.gz" ]
    ]
  }
}

process SeqPurge {
  tag { "seqpurge_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }

  input:
  each threads from params.threadCounts
  each err from params.errorRates
  each qcut from params.quals

  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile
  stdout timing

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

  timing.onComplete {
    perfs << [ 
      "name" : name,
      "value" : it
    ]
  }

  afterScript {
    results << [ 
      "name" : name,
      "trimmed" : [ output1, output2 ]
    ]
  }
}

process AdapterRemoval {
  tag { "adapterremoval_${task.cpus}_${err}_q${qcut}" }
  cpus { threads }

  input:
  each threads from params.threadCounts
  each err from params.errorRates
  each qcut from params.quals

  val fqPrefix from "sim_${err}"
  file input1 from "${params.dataDir}/${fqPrefix}.1.fq"
  file input2 from "${params.dataDir}/${fqPrefix}.2.fq"

  output:
  file "${name}.1.fq.gz" into output1
  file "${name}.2.fq.gz" into output2
  file "${name}.report.txt" into reportFile
  stdout timing

  script:
  """
  usr/bin/time -v AdapterRemoval \
    --file1 $input1 --file2 $input2 \
    --output1 $output1 --output2 $output2 --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut \
    --minLengthgth $params.minLength --task.cpus $task.cpus
  """

  timing.onComplete {
    perfs << [ 
      "name" : name,
      "value" : it
    ]
  }

  afterScript {
    results << [ 
      "name" : name,
      "trimmed" : [ output1, output2 ]
    ]
  }
}

process Accuracy {
  input:
  each result from results

  output:
  file "${result.name}.txt" into resultFile
  file "${result.name}.summary.txt" into summaryFile
  file "${result.name}.table.txt" into tableFile

  script:
  """
  python scripts/summarize_simulated_trimming_accuracy.py \
    -a1 ${params.dataDir}/${fqPrefix}.1.aln \
    -a2 ${params.dataDir}/${fqPrefix}.2.aln \
    -r1 ${result.trimmed[0]} -r2 ${result.trimmed[1]} \
    --name ${result.name} \
    -o $resultFile -s $summaryFile -t $tableFile
  """
}

process Timing {
  input:
  stdin timingInfo from perfs

  output:
  file "${process.executor}.timing.tex" into table

  script:
  """
  python scripts/summarize_timing_info.py -o $table -f latex
  """
}
