# Pipeline for the analysis of Oxford Nanopore Technologies human plasma sequencing data

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

[![Follow on Twitter](http://img.shields.io/badge/twitter-%40JVerwilt-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/JVerwilt)

## Introduction
This pipeline is written specifically for the analysis of data originating from our Oxford Nanopore Technologies-based mRNA sequencing protocol. It is tuned towards a faithful quantification of different isoforms (including novel ones). The pipeline uses a Docker image which is pulled from Docker Hub, when run with the ```docker``` profile. A general overview of the pipeline is visualized in the following flowchart: 

![flowchart](https://github.com/OncoRNALab/plasmaONT/blob/main/flowchart.png/?)

The pipeline starts by combining all reads. Then it uses [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI) to determine the adapter sequences and trim these sequences. The primers found by Porechop_ABI are then used by [Pychopper v2](https://github.com/epi2me-labs/pychopper) to filter for the end-to-end reads (containing both adapter sequences). The resulting FASTQ files are mapped using [Minimap2](https://github.com/lh3/minimap2) in guided mode (a GFF file is provided to prioritize existing splice sites). 

