#!/bin/bash
#PBS -N subsamp
#PBS -l nodes=1:ppn=16
#PBS -l walltime=70:59:00

cd $PBS_O_WORKDIR

ml purge
ml seqtk/1.3-foss-2018a

for file in *_R[1-2]_001.fastq.gz
do
  for number in {1,5,10,15,20,30}
  #  for number in 34
  do
  seed=0
  seed=`expr $seed + $number`
  echo "Extracting reads: ${number}M from $file with seed $seed"
  READ=$(echo $file | cut -f 3 -d"_")
  READ=$(echo $READ | cut -f 1 -d".")
  SAMPLE=${file%%_${READ}_001.fastq.gz}
  echo "        basename = $SAMPLE and read = $READ"
  R1=${SAMPLE}_R1_001.fastq.gz
  R2=${SAMPLE}_R2_001.fastq.gz
  mkdir -p ${SAMPLE}-subsamp${number}M
  seqtk sample -s$seed $file ${number}000000 | gzip > ${SAMPLE}-subsamp${number}M/${SAMPLE}-subsamp${number}M-s${seed}_${READ}.fastq.gz
  done
done
