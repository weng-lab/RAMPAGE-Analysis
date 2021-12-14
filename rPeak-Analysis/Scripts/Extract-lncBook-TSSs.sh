#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

gtf=../Input-Data/lncBook/lncRNA_LncBook_GRCh38_9.28.gtf

awk '{if ($3 == "transcript" && $1 ~ /chr/) print $1 "\t" $4 "\t" $5 "\t" $12 \
    "\t" $6 "\t" $7 "\t" $10 "\t" $18"\t"$14"\t"$20}' $gtf | \
    awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' | \
    awk '{if ($1 ~ /chr/) print $0}' > tmp.2
 
awk '{if ($6 == "+") print $1 "\t" $2-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7;\
    else print $1 "\t" $3-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' tmp.2 | \
    sort -k1,1 -k2,2n > ../Input-Data/lncBook/lncBook-lncRNA-TSSs.bed

rm tmp.2
