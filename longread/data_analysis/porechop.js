#!/bin/bash
#PBS -N porechop
#PBS -l nodes=1:ppn=12
#PBS -l walltime=24:00:00
#PBS -l mem=180gb

input_reads_sam="/scratch/gent/vo/000/gvo00027/projects/HVD/Anchoring/fastq"
out_DIR="/scratch/gent/vo/000/gvo00027/projects/HVD/Anchoring/fastq"
work_DIR=""
sample="sample"

cd $work_DIR

module purge
ml SAMtools/1.16.1-GCC-11.3.0

samtools fastq $input_reads_sam > ${out_DIR}/${sample}.fastq
gzip ${out_DIR}/${sample}.fastq

module purge
ml Porechop/0.2.4-intel-2019b-Python-3.7.4

porechop -i ${out_DIR}/${sample}.fastq.gz -b $out_DIR --untrimmed > ${out_DIR}/${sample}_barcodinglog.txt
