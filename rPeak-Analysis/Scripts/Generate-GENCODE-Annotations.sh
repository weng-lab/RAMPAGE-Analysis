#!/bin/bash

#Jill E. Moore
#Weng Lab
#UMass Medical School
#December 2021

gtf=$1
mode=$2
scriptDir=~/Projects/ENCODE/Encyclopedia/Version5/Annotation-Files

col=14

awk '{if ($3 == "transcript" && $1 ~ /chr/) print $1 "\t" $4 "\t" $5 "\t" $12 \
    "\t" $6 "\t" $7 "\t" $10 "\t" $18"\t"$14"\t"$20}' $gtf | \
    awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' | \
    awk '{if ($1 ~ /chr/) print $0}' > tmp.2

 
awk '{if ($6 == "+") print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7;\
    else print $1 "\t" $3 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' tmp.2 | \
    sort -k1,1 -k2,2n > TSS.$mode.bed

awk '{if ($2 > 2000) print $1 "\t" $2-2000 "\t" $3+2000 "\t" $4 "\t" $5 "\t" \
    $6 "\t" $7; else print $1 "\t" 0 "\t" $3+2000 "\t" $4 "\t" $5 "\t" \
    $6 "\t" $7}' TSS.$mode.bed > TSS.$mode.4K.bed

awk '{if ($9 == "protein_coding") print $0}' tmp.2 | awk '{if ($6 == "+") \
    print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7;\
    else print $1 "\t" $3 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' | \
    sort -k1,1 -k2,2n > TSS.$mode-PC.bed

awk '{if ($2 > 2000) print $1 "\t" $2-2000 "\t" $3+2000 "\t" $4 "\t" $5 "\t" \
    $6 "\t" $7; else print $1 "\t" 0 "\t" $3+2000 "\t" $4 "\t" $5 "\t" \
    $6 "\t" $7}' TSS.$mode-PC.bed > TSS.$mode-PC.4K.bed

awk '{if ($3 == "gene") print $1 "\t" $4 "\t" $5 "\t" $10 "\t" "." "\t" $7 \
    "\t" $'$col' "\t" $12 }' $gtf | awk '{gsub(/;/," ");print}' | \
    awk '{gsub(/"/,"");print}' | awk '{if ($1 ~ /chr/) print $0}' | sort -k1,1 -k2,2n > Genes.$mode.bed

awk '{if ($8 == "protein_coding") print $0}' Genes.$mode.bed | sort -k1,1 -k2,2n > Genes.$mode-PC.bed

if [[ $gtf == *"24"* ]]
then
    awk '{if ($3 == "exon") print $1 "\t" $4 "\t" $5 "\t" $28 "\t" 0 "\t" $7 "\t" $12 "\t" \
        $10 "\t" $26 "\t" $18}' $gtf | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' \
        > Exon-Annotated.$mode.bed
    awk '{if ($3 == "transcript") print $12 "\t" $10 "\t" $18 "\t" $14 "\t" $20}' $gtf \
        | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' > Transcript-Gene-Table.$mode.txt
else
    awk '{if ($3 == "exon") print $1 "\t" $4 "\t" $5 "\t" $24 "\t" 0 "\t" $7 "\t" $12 "\t" \
        $10 "\t" $22 "\t" $16}' $gtf | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' \
        > Exon-Annotated.$mode.bed
    awk '{if ($3 == "transcript") print $12 "\t" $10 "\t" $16 "\t" $18 "\t" $14}' $gtf \
        | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' > Transcript-Gene-Table.$mode.txt
fi

awk '{if ($3 == "transcript") print $1 "\t" $4 "\t" $5 "\t" $12 "\t" "." "\t" $7 "\t" $10}' $gtf \
    | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' > Transcripts.$mode.bed

awk '{if ($9 == "1") print $1 "\t" $2 "\t" $3 "\t" $7 "\t" "." "\t" $6}' Exon-Annotated.$mode.bed > \
    Exon-One.$mode.bed

awk '{if ($9 > 1) print $7}' Exon-Annotated.$mode.bed | \
    awk 'FNR==NR {x[$1];next} !($7 in x)' -  Exon-Annotated.$mode.bed | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $7 "\t" "." "\t" $6}' > Only-One-Exon.$mode.bed 

rm tmp.2
