#!/bin/sh
var1=$(echo $1)
var2=$(echo $2)
var3=$(echo $3)
var4=$(echo $4)

paste \
    <(awk 'BEGIN{OFS="\t"}{for (i=1;i<=3;i++) {printf $0"\t";getline} print $0}' $var1) \
    <(awk 'BEGIN{OFS="\t"}{for (i=1;i<=3;i++) {printf $0"\t";getline} print $0}' $var2) | 
    awk -v var3=${var3} -v var4=${var4} 'BEGIN{FS="\t"} {b=substr($6,1,8)} b!="NNNNNNNN" && b!="GGGGGGGG" && b!="AAAAAAAA" && b!="TTTTTTTT" && b!="CCCCCCCC" {for (i=1;i<=4;i++) {print $i > var3} for (i=5;i<=8;i++) {print $i > var4}}'

lines=4
sample="$(basename -- $1)"
sample_lines=`cat ${1} | wc -l`; 
echo $((sample_lines/lines)) $(echo $sample)  >> polyNG_filter.reads.txt

sample="$(basename -- $3)"
sample_lines=`cat ${3} | wc -l`; 
echo $((sample_lines/lines)) $(echo $sample) >> polyNG_filter.reads.txt