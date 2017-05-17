#!/usr/bin/env nextflow

/* Atropos paper workflow for real WGBS reads
 * ------------------------------------------
 * Configuration is externalized into separate profiles for local and cluster.
 * The main difference is that Docker is used to manage images on desktop, 
 * verusus Singularity on the cluster.
 *
 * Parameters
 * ----------
 * The following are expected to be defined in params:
 * - threadCounts: The numbers of threads (cpus) to test
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
 *   Docker/Singularity by commenting out the 'container' directives.e.
 */

// variables for all tools
params.publishDir = "results/${workflow.profile}/wgbs"
params.quals = [ 0, 20 ]
params.adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
params.adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
params.minLength = 25
params.batchSize = 5000

// atropos-specific variables
params.aligners = [ 'insert' ]
params.compressionSchemes = [ 'worker', 'writer', 'nowriter' ]

/* Process: Extract data from container
 * ------------------------------------
 * The read data is contained within a Docker container. Since Singularity does
 * not support attaching volumes, we use a process to unpack the data from the 
 * container.
 */
process ExtractReads {
  container "jdidion/atropos_wgbs"
  
  output:
  file("wgbs.{1,2}.fq.gz") into wgbsReads
  
  script:
  """
  cp \
    /data/wgbs/wgbs.1.fq.gz \
    /data/wgbs/wgbs.2.fq.gz \
    .
  """
}

// Split the input reads channel into one per tool process.
wgbsReads.into {
  untrimmedWgbsReads
  atroposWgbsReads
  skewerWgbsReads
  seqPurgeWgbsReads
  adapterRemovalWgbsReads
}

/* Process: Atropos adapter trimming
 * ---------------------------------
 */
process Atropos {
  tag { "atropos_${task.cpus}_wgbs_q${qcut}_${aligner}_${compression}" }
  cpus { threads }
  container "jdidion/atropos_paper"

  input:
  file(reads) from atroposWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals
  each aligner from params.aligners
  each compression from params.compressionSchemes
  
  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAtropos
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAtropos
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
    alignerArgs = "--insert-match-error-rate 0.30"
  }
  else {
    alignerArgs = "-O 7"
  }
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt atropos \
    --op-order GACQW -T $task.cpus --batch-size $params.batchSize \
    -e 0.20 --aligner $aligner $alignerArgs -q $qcut --trim-n \
    -a $params.adapter1 -A $params.adapter2  -m $params.minLength \
    --no-default-adapters --no-cache-adapters --log-level ERROR \
    --correct-mismatches liberal \
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
  tag { "skewer_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/skewer"

  input:
  file(reads) from skewerWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedSkewer
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSkewer
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt skewer \
    -m pe -l $params.minLength -r 0.3 \
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
  tag { "seqpurge_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/seqpurge"

  input:
  file(reads) from seqPurgeWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedSeqPurge
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingSeqPurge
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt SeqPurge \
    -in1 ${reads[0]} -in2 ${reads[1]} \
    -out1 ${task.tag}.1.fastq.gz -out2 ${task.tag}.2.fastq.gz \
    -a1 $params.adapter1 -a2 $params.adapter2 \
    -qcut $qcut -min_len $params.minLength \
    -threads $task.cpus -prefetch $params.batchSize \
    -ec -match_perc 70 -summary ${task.tag}.report.txt
  """
}

/* Process: AdapterRemoval2 adapter trimming
 * -----------------------------------------
 */
process AdapterRemoval {
  tag { "adapterremoval_${task.cpus}_wgbs_q${qcut}" }
  cpus { threads }
  container "jdidion/adapterremoval"

  input:
  file(reads) from adapterRemovalWgbsReads
  each threads from params.threadCounts
  each qcut from params.quals

  output:
  set val("${task.tag}"), file("${task.tag}.{1,2}.fq.gz") into trimmedAdapterRemoval
  set val("${task.tag}"), file("${task.tag}.timing.txt") into timingAdapterRemoval
  file "${task.tag}.report.txt"
  
  script:
  """
  /usr/bin/time -v -o ${task.tag}.timing.txt AdapterRemoval \
    --file1 ${reads[0]} --file2 ${reads[1]} \
    --output1 ${task.tag}.1.fq.gz --output2 ${task.tag}.2.fq.gz --gzip \
    --adapter1 $params.adapter1 --adapter2 $params.adapter2 \
    --trimns --trimqualities --minquality $qcut --mm 0.3 \
    --minlength $params.minLength --threads $task.cpus
    > "${task.tag}.report.txt"
  """
}

/* Channel: Untrimmed
 * ------------------
 * Need to create a set from the untrimmed reads so it can be merged with
 * the trimmed reads.
 */
untrimmed = Channel.value(["untrimmed", untrimmedWgbsReads])

/* Channel: merged tool outputs
 * ----------------------------
 * We also add in the untrimmed output, which we need to align for evaluation of
 * effectiveness. If you add a tool process, make sure to add the trimmedXXX 
 * channel into the concat list here.
 */
Channel
  .empty()
  .concat(
    untrimmed,
    trimmedAtropos,
    trimmedSkewer,
    trimmedSeqPurge,
    trimmedAdapterRemoval
  )
  .set { trimmedMerged }

/* Process: Align reads using STAR
 * -------------------------------
 * Alignments are written to an unsorted BAM file and then name sorted by
 * samtools.
 */
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

/* Process: Summarize trimming effectiveness
 * -----------------------------------------
 */
process ComputeEffectiveness {
  input:
  val bamFileList from sortedEffectiveness.toList()
  val bedFileList from overlap.toList()
  
  output:
  file "effectiveness.txt" into effectiveness
  
  script:
  bamFiles = bamFileList.join(" ")
  bedFiles = bedFileList.join(" ")
  """
  compute_real_effectiveness.py \
    -i $bamFileList -o effectiveness.txt \
    --no-edit-distance --no-progress
  """
}

/* Process: Generate plot
 * ----------------------
 */
process ShowEffectiveness {
  input:
  file effData from effectiveness
  
  output:
  file "wgbs_effectiveness.svg"
  
  script:
  """
  show_wgbs_effectiveness.py \
    -i $effData -o wgbs_effectiveness.svg \
    --exclude-discarded
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
process ParseTrimmingTiming {
    container "jdidion/python_bash"
    
    input:
    set val(name), file(timing) from timingMerged
    
    output:
    stdout timingMergedParsed
    
    script:
    """
    parse_gtime.py -i $timing -p $name
    """
}

/* Process: generate performance figure/table
 * ------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowTrimmingPerformance {
    container "jdidion/python_bash"
    publishDir "$publishDir", mode: 'copy', overwrite: true
    
    input:
    val parsedRows from timingMergedParsed.toList()
    
    output:
    file "trim_performance.tex"
    file "trim_performance.svg"
    file "trim_performance.pickle"
    
    script:
    data = parsedRows.join("")
    """
    echo '$data' | show_performance.py -o trim_performance -f tex svg pickle
    """
}

/* Process: parse performance for STAR alignment
 * ---------------------------------------------
 */
process ParseStarTiming {
    container "jdidion/python_bash"
    
    input:
    set val(name), file(timing) from timingBwameth
    
    output:
    stdout timingBwamethParsed
    
    script:
    """
    parse_gtime.py -i $timing -p $name
    """
}

/* Process: generate STAR alignment performance figure/table
 * ---------------------------------------------------------
 * Aggregate all the parsed performance data and pass it to stdin of the
 * bin/show_performance.py script.
 */
process ShowStarPerformance {
    container "jdidion/python_bash"
    publishDir "$publishDir", mode: 'copy', overwrite: true
    
    input:
    val parsedRows from timingBwameth.toList()
    
    output:
    file "bwameth_performance.tex"
    file "bwameth_performance.svg"
    
    script:
    data = parsedRows.join("")
    if (workflow.profile == "local") {
      name = task.ext.bwameth_local_name
      caption = task.ext.bwameth_local_caption
    } else {
      name = task.ext.bwameth_luster_name
      caption = task.ext.bwameth_cluster_caption
    }
    """
    echo '$data' | show_performance.py -n $name -c $caption -o bwameth_performance
    """
}
