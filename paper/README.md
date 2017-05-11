# Overview

The scripts in this directory will enable you to re-run the analyses in the Atropos paper. The workflows defined here run the benchmarks and generate the figures and tables shown in the paper.

We have created [Docker](https://www.docker.com/) images for all of the software tools used, as well as data volumes containing all the raw data and resources. These images can be used directly on Mac, Windows, and some linux platforms using the Docker engine. On unsupported linux platforms (namely RedHat and derivatives, such as Scientific Linux), [Singularity](http://singularity.lbl.gov/) can be used to execute the containers directly from Docker Hub. The versions of the tools used in the paper are noted in the Dockerfile headers, and also in the supplementary data.

Our workflows are written in [Nextflow](https://www.nextflow.io/index.html), primarily because it supports both Docker and Singularity, which we need to run benchmarks on both desktop and RedHat-based HPC cluster. We also provide [CWL](http://www.commonwl.org/) tool definitions to simplify the development of alternate workflows.

Each workflow (.nf file) runs the analysis for one data type (RNA-Seq, WGBS, or simulated DNA-Seq). We provide a configuration file with profiles we used for both the local and cluster executions. Our cluster runs SGE, so you may need to alter the cluster configuration files for your environment.

# 1. Install software

* You will need a [Docker](https://www.docker.com/) engine if you want to build the containers yourself. If you only want to run the containers, you can use either Docker or [Singularity](http://singularity.lbl.gov/).
* [Nextflow](https://www.nextflow.io/index.html), which requires Java 7+.

# 2. Build containers

All of the containers defined in the 'containers' subdirectory have already been built and pushed to Docker Hub, with two exceptions: the data containers for the STAR indexes (data/hg37/star-index and data/hg38/star-index) are too large to be pushed to Docker Hub or Quay.io. Thus, you will unfortunately need to build at least one of them yourself. We use GRCh38 in the paper, so to build that container, clone data/hg38/star-index and run the build.sh script in that directory.

For full reproducibility, you are free to build the containers yourself, but you'll need to create your own account on Docker Hub, and you'll need to update the scripts to push/pull containers from your own repository. Build all the tool containers, then build all the data containers.

In general, for each tool container, run the following sequence of commands:

    # Build the container
    docker build -f Dockerfile -t <your repo>/<tool name> .
    
    # Upload the container to your Docker Hub repo
    docker push <your repo>/<tool name>

For each data container, run the following sequence of commands:

    # Download the data and build the docker container
    ./build.sh
    
    # Upload the container to your Docker Hub repo
    docker push <your repo>/<data name>

On a technical note, we use Phusion (https://hub.docker.com/r/phusion/baseimage/) as the base image for the containers for the tools we benchmark. This is primarily for convenience and comparability (i.e. removing base image as a variable); you could build smaller images using Alpine with a little extra work.

# 3. Run the workflows

Clone the files in the 'workflow' directory, including the 'bin' subdirectory. In the 'workflow' directory, first edit the nextflow.config file and customize it to your own computing environment. Then run:

    ./run-workflows.sh <env>

Where <env> is either 'local' or 'cluster'. Note that the first time you run this it will download several Docker images requiring ~[XX] GB of disk space.

All results will be placed in the 'results' subdirectory. Note that for the paper we tweaked the aesthetics of some of the figures and tables, but the data was left exactly as-is. Of course there will be some variability in the performance metrics, but the relative rankings of the tools should not change significantly -- please let us know if you find otherwise!

# TODO

* Update #4 with an estimate of the total disk space requirement.
* These papers do a nice job of benchmarking trimmers. Consider adding some more benchmarks.
    * http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2
    * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0454-y
