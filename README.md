# Pipeline for the analysis of Oxford Nanopore Technologies human plasma sequencing data

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

[![Follow on Twitter](http://img.shields.io/badge/twitter-%40JVerwilt-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/JVerwilt)

## Introduction
This pipeline is written specifically for the analysis of data originating from our Oxford Nanopore Technologies-based mRNA sequencing protocol. It is tuned towards a faithful quantification of different isoforms (including novel ones). The pipeline uses a Docker image which is pulled from Docker Hub, when run with the ```docker``` profile. A general overview of the pipeline is visualized in the following flowchart: 

![flowchart](https://github.com/OncoRNALab/plasmaONT/blob/main/flowchart.png/?)

The pipeline starts by combining all reads. Then it uses [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI) to determine the adapter sequences and trim these sequences. The primers found by Porechop_ABI are then used by [Pychopper v2](https://github.com/epi2me-labs/pychopper) to filter for the end-to-end reads (containing both adapter sequences). The resulting FASTQ files are mapped using [Minimap2](https://github.com/lh3/minimap2) in guided mode (a GFF file is provided to prioritize existing splice sites) and information on the mapped reads is summarized and visualized using [NanoPlot](https://github.com/wdecoster/NanoPlot). Reads generated by internal priming are then filtered out using [polyAfilter](https://github.com/MarekSvob/polyAfilter). The pipeline detects isoforms using the two most prominent and high-performing isoform detection algorithms: [IsoQuant](https://github.com/ablab/IsoQuant) and [StringTie2](https://github.com/skovaka/stringtie2). The trimmed and end-to-end reads are quantified using [HTSeq](https://htseq.readthedocs.io/en/master/) with three different GTF files: the original one containing all Ensembl exons, the one generated by StringTie2 and the one generated by IsoQuant. Last, all information from compatible software is collected and reported in a [MultiQC](https://multiqc.info/) report. 


