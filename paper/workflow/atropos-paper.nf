#!/usr/bin/env nextflow

/* Atropos paper workflow
 * ----------------------
 * Configuration is externalized into separate files for desktop and cluster.
 * The main difference is that Docker images are used on desktop, verusus
 * Singularity on the cluster.
 */

/* Parameters
 * ----------
 * The following are exptected to be defined externally:
 * - 