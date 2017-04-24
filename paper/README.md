# Overview

The scripts in this directory will enable you to re-run the analyses in the Atropos paper. The workflows defined here runs the benchmarks and generates the figures and tables shown in the paper.

We have created [Docker](https://www.docker.com/) images for all of the software tools used, as well as data volumes containing all the raw data and resources. These images can be used directly on Mac, Windows, and some linux platforms using the Docker engine. On unsupported linux platforms (namely RedHat and derivatives, such as Scientific Linux), [Singularity](http://singularity.lbl.gov/) can be used to execute the containers directly from Docker Hub. 

Our workflows are written in [Nextflow](https://www.nextflow.io/index.html), primarily because it supports both Docker and Singularity, which we need to run benchmarks on both desktop and RedHat-based HPC cluster. We also provide [CWL](http://www.commonwl.org/) tool definitions to simplify the development of alternate workflows.

Each workflow (.nf file) runs the analysis for one data type (RNA-Seq, WGBS, or simulated DNA-Seq). We provide the configuration files we used for both the local and cluster executions. Our cluster runs SGE, so you may need to alter the cluster configuration files for your environment.

# 1. Install software

* [Docker](https://www.docker.com/) and/or [Singularity](http://singularity.lbl.gov/)
* Java 7+
* [Nextflow](https://www.nextflow.io/index.html)
* If you want to generate your own simulated reads, rather than use the ones we provide, you will need:
    * [R](https://www.r-project.org/about.html)
    * A modern C++ compiler (e.g. GCC 5.x)
    * automake

# 2. Simulate reads (optional)

There are already 3 simulated data sets. If you'd like to re-create these, in the workflow/data/simulated directory run:

    ./install.sh
    ./simulate_reads.sh

# 3. Run the workflows

In the 'workflow' directory, run:

    ./run-<mode>-workflows.sh

Where <mode> is either 'local' or 'cluster'. Note that the first time you run this it will download several Docker images requiring ~XX GB of disk space.

# TODO

* These papers do a nice job of benchmarking trimmers. Consider adding some more benchmarks.
    * http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2
    * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0454-y
