The scripts in this directory should enable you to re-run the analyses in the Atropos paper.

# Overview

Our workflow is written in [Common Workflow Language (CWL) v1.0](http://www.commonwl.org/v1.0/). We have created Docker images for all of the software tools used. For our paper, we executed the workflow using [Toil](https://toil.readthedocs.io).

# 1. Install software

* [Docker](https://www.docker.com/)
* [Toil](https://toil.readthedocs.io), or another workflow execution system that supports CWL. Note that Toil requires python 2.7.
* If you want to generate your own simulated reads, rather than use the ones we provide, you will need:
    * [R](https://www.r-project.org/about.html)
    * A modern C++ compiler (e.g. GCC 5.x)
    * automake

# 2. Simulate reads

There are already 3 simulated data sets in the 'data' directory. If you'd like to re-create these, in the 'simulation' directory run:

    ./install.sh
    ./simulate_reads.sh

# 3. Prepare and run the workflow

In the 'workflow' directory, run:

    ./run-<mode>-workflow.sh

Where <mode> is either 'local' or 'cluster'. Note that the first time you run this it will download several Docker images requiring ~10 GB of disk space.

# TODO

* These papers do a nice job of benchmarking trimmers. Consider adding some more benchmarks.
    * http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2
    * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0454-y
