This directory contains all the containers used in the Atropos paper. There are two types of containers:

* Data: Minimal base container (busybox) with no software installed. Each container houses a single dataset.
* Tools: Phusion base container with software required by the analysis. Containers are designed such that each Nextflow process uses a single container.
