# Data containers

These are containers for the raw data and resources used in the Atropos paper. All containers are based on Busybox for the greatest space efficiency. Each folder has a build.sh script, which downloads the data and builds the image based on the Dockerfile in the same directory.

* hg37 and hg38
    * reference: Reference sequence and GENCODE v26 gene annotations.
    * bwa-index: BWA index, using jdidion/bwabase as a baseimage.
    * bwameth-index: bwa-meth index, using jdidion/bwabase as a base image.
    * star-index: STAR index, using jdidion/starbase as a base image.
* simulated, rna, wgbs: benchmark datasets used in Atropos paper
