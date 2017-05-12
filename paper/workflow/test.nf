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
  container "jdidion/atropos_simulated"
  
  input:
  each err from params.errorRates
  
  output:
  set err, "sim_${err}.{1,2}.fq", "sim_${err}.{1,2}.aln" into simReads
  
  script:
  """
  head -40 /data/simulated/sim_${err}.1.fq > ./sim_${err}.1.fq && \
  head -40 /data/simulated/sim_${err}.2.fq > ./sim_${err}.2.fq && \
  head -400 /data/simulated/sim_${err}.1.aln > ./sim_${err}.1.aln && \
  head -400 /data/simulated/sim_${err}.2.aln > ./sim_${err}.2.aln
  """
}

process Atropos {
  //tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  tag { "atropos_${task.cpus}_${err}_q0_insert_writer" }
  cpus { threads }
  container "jdidion/atropos_paper"
  
  input:
  set err, file(reads), file(alns) from simReads
  each threads from params.threadCounts
  //each qcut from params.quals
  //each aligner from params.aligners
  //each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), val(err), file(alns), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    -pe1 ${reads[0]} -pe2 ${reads[1]} \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz \
    --report-file ${task.tag}.report.txt --quiet \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
  """
}

// concatenate all of the timing results into a single channel
Channel
  .empty()
  .concat(trimmedAtropos)
  .set { trimmedMerged }
Channel
  .empty()
  .concat(timingAtropos)
  .set { timingMerged }


process ParseTiming {
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

process ShowPerformance {
    container "jdidion/python_bash"
    
    input:
    val timingRows from timingParsed.toList()
    
    output:
    file "timing.tex"
    file "timing.svg"
    
    script:
    data = timingRows.join("")
    if (workflow.profile == "local") {
      name = task.ext.local_name
      caption = task.ext.local_caption
    } else {
      name = task.ext.cluster_name
      caption = task.ext.cluster_caption
    }
    """
    echo '$data' | show_performance.py -n $name -c $caption -o timing -f tex svg pickle
    """
}

process ComputeAccuracy {
  container "jdidion/python_bash"
  
  input:
  set val(name), val(err), file(alns), file(trimmed) from trimmedMerged
  
  output:
  file "${name}.txt" into resultFile
  file "${name}.summary.txt" into summaryFile
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

process ShowAccuracy {
  echo true
  container "jdidion/python_bash"
  
  input:
  val accuracyRows from tableFile.toList()
  
  output:
  file "accuracy.tex"
  
  script:
  data = accuracyRows.join("")
  """
  echo '$data' | show_simulated_accuracy.py -n foo -c bar -o accuracy
  """
}
