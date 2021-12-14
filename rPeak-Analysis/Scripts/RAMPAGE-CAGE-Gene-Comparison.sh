#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

dataDir=~/Lab/ENCODE/RAMPAGE
rampage100=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks-100bp.bed
cage100=~/Lab/ENCODE/RAMPAGE/FANTOM-CAGE/hg38_fair+new_CAGE_100bp_Summits.bed
refDir=~/Lab/Reference/Human/hg38/GENCODE31
tss=$refDir/TSS.Basic.bed

mkdir -p $dataDir/Gene-Overlap
cd $dataDir/Gene-Overlap

awk '{print $1 "\t" $2-500 "\t" $2+500 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' $tss \
    > proximal500.bed

bedtools intersect -u -s -a proximal500.bed -b $rampage100 > RAMPAGE-All-TSS.txt
bedtools intersect -u -s -a proximal500.bed -b $cage100 > CAGE-All-TSS.txt

cat CAGE-All-TSS.txt RAMPAGE-All-TSS.txt | awk '{print $NF}' \
    | sort -u > All-Genes.txt
awk 'FNR==NR {x[$NF];next} !($NF in x)' CAGE-All-TSS.txt RAMPAGE-All-TSS.txt \
    | awk '{print $NF}' | sort -u > RAMPAGE-Only-Genes.txt
awk 'FNR==NR {x[$NF];next} !($NF in x)' RAMPAGE-All-TSS.txt CAGE-All-TSS.txt \
    | awk '{print $NF}' | sort -u > CAGE-Only-Genes.txt
awk 'FNR==NR {x[$NF];next} ($NF in x)' RAMPAGE-All-TSS.txt CAGE-All-TSS.txt \
    | awk '{print $NF}' | sort -u > Common-Genes.txt 

rm proximal500.bed
