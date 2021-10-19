#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#October 2021

biosample=$1
header=$2
bioLower=$(echo $biosample | awk '{print tolower($0)}')
dataDir=~/Lab/ENCODE/RAMPAGE/GROcap

cd $dataDir

~/bin/liftOver ~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/RAMPAGE-$biosample-RPM-2.peaks.bed ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain  tmp.hg19 tmp.no

totalP=$(wc -l "tss_SS_"$bioLower"_plus.bed" | awk '{print $1}')
totalN=$(wc -l "tss_SS_"$bioLower"_minus.bed" | awk '{print $1}')

bedtools intersect -u -s -a "tss_SS_"$bioLower"_plus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalP' "\t" $1/'$totalP'}'
bedtools intersect -u -s -a "tss_SS_"$bioLower"_minus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalN' "\t" $1/'$totalN'}'

totalP=$(wc -l "tss_US_"$bioLower"_plus.bed" | awk '{print $1}')
totalN=$(wc -l "tss_US_"$bioLower"_minus.bed" | awk '{print $1}')

bedtools intersect -u -s -a "tss_US_"$bioLower"_plus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalP' "\t" $1/'$totalP'}'
bedtools intersect -u -s -a "tss_US_"$bioLower"_minus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalN' "\t" $1/'$totalN'}'

totalP=$(wc -l "tss_SU_"$bioLower"_plus.bed" | awk '{print $1}')
totalN=$(wc -l "tss_SU_"$bioLower"_minus.bed" | awk '{print $1}')

bedtools intersect -u -s -a "tss_SU_"$bioLower"_plus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalP' "\t" $1/'$totalP'}'
bedtools intersect -u -s -a "tss_SU_"$bioLower"_minus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalN' "\t" $1/'$totalN'}'

totalP=$(wc -l "tss_UU_"$bioLower"_plus.bed" | awk '{print $1}')
totalN=$(wc -l "tss_UU_"$bioLower"_minus.bed" | awk '{print $1}')

bedtools intersect -u -s -a "tss_UU_"$bioLower"_plus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalP' "\t" $1/'$totalP'}'
bedtools intersect -u -s -a "tss_UU_"$bioLower"_minus.bed" -b tmp.hg19 | wc -l | awk '{print $1 "\t" '$totalN' "\t" $1/'$totalN'}'
