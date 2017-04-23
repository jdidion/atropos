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
 * - min_len: minimum read length
 * - batch_size: read batch size
 * - quals: quality thresholds for trimming
 * - aligners: Atropos aligner algorithms to use
 * - adapter1, adapter2: Adapter sequence
 * - extra: map of extra arguments for each program
 */
