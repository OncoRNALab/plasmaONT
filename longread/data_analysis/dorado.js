#!/bin/sh
#PBS -N dorado_caller
#PBS -l nodes=1:ppn=10:gpus=1
#PBS -l walltime=72:00:00
#PBS -l mem=90gb
#PBS -m abe

orig_DIR="dir_where_pod5_are"
out_DIR="dir_where_sam_should_go"
work_DIR="dir_where_you_want_logfiles"
config="dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
sample="sample"

cd $work_DIR

module purge
ml dorado/0.3.1-foss-2022a-CUDA-11.7.0
mkdir $out_DIR
dorado download --model $config
dorado basecaller $config $orig_DIR > ${out_DIR}/${sample}.sam
