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

process ExtractReads {
  container "jdidion/atropos_wgbs"
  
  output:
  file("wgbs.{1,2}.fq") into wgbsReads
  
  script:
  """
  cp \
    zcat /data/wgbs/wgbs.1.fq.gz | head -40 > ./wgbs.1.fq \
    zcat /data/wgbs/wgbs.2.fq.gz | head -40 > ./wgbs.2.fq \
    .
  """
}

process Atropos {
  //tag { "atropos_${task.cpus}_${err}_q${qcut}_${aligner}_${compression}" }
  tag { "atropos_${task.cpus}_${err}_q0_insert_writer" }
  cpus { threads }
  container "jdidion/atropos_paper"
  
  input:
  file(reads) from wgbsReads
  each threads from params.threadCounts
  //each qcut from params.quals
  //each aligner from params.aligners
  //each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    --op-order GACQW -T 4 \
    -e 0.20 --aligner insert --insert-match-error-rate 0.30 -q 0 --trim-n \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -m 25 --no-default-adapters --no-cache-adapters --log-level ERROR \
    --correct-mismatches liberal \
    -o ${task.tag}.1.fq.gz -p ${task.tag}.2.fq.gz \
    --report-file ${task.tag}.report.txt --quiet \
    -pe1 ${reads[0]} -pe2 ${reads[1]}
  """
}

// concatenate all of the timing results into a single channel
Channel
  .empty()
  .concat(trimmedAtropos)
  .set { trimmedMerged }

process BwamethAlign {
  container "jdidion/bwa_hg38index"
  cpus { params.alignThreads }
  
  input:
  set val(name), file(fastq) from trimmedMerged
  
  output:
  file("${name}_wgbs.{bam,bam.bai}")
  set val(name), file("${name}_wgbs.name_sorted.bam") into sorted
  set val(name), file("${name}.bwameth.timing.txt") into timingBwameth
  
  script:
  """
  /usr/bin/time -v -o ${name}.star.timing.txt bwameth.py \
    -z -t ${params.alignThreads} --read-group '${task.ext.readGroup}' \
    --reference /data/index/bwameth/hg38/hg38
    -o ${name}_wgbs.bam ${fastq[0]} ${fastq[1]} \
  && samtools sort -n -O bam -@ ${params.alignThreads} \
    -o ${name}_wgbs.name_sorted.bam ${name}_wgbs.bam
  """
}

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
    """
    echo '$data' | show_performance.py -o timing -f tex svg pickle
    """
}
