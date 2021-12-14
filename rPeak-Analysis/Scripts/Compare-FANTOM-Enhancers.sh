#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

dataDir=~/Lab/ENCODE/RAMPAGE/NET-CAGE/
rPeaks=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
fantomEnhancers=$dataDir/Supplementary_Data_1_Human_FANTOM-NET_enhancers.bed
hg38Enhancers=$dataDir/Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.bed

~/bin/liftOver -bedPlus=4 $fantomEnhancers ~/Lab/Reference/Human/hg19/hg19ToHg38.over.chain \
    $hg38Enhancers no-map

totalCAGE=$(awk '{if ($4 ~ /FANTOM/) print $0}' $hg38Enhancers | wc -l | awk '{print $1}')
totalNET=$(awk '{if ($4 ~ /NET/) print $0}' $hg38Enhancers | wc -l | awk '{print $1}')

bedtools intersect -u -a $dataDir/Supplementary_Data_1_Human_FANTOM-NET_enhancers.hg38.bed -b $rPeaks > tmp.out

awk '{if ($4 ~ /FANTOM/) print $0}' tmp.out | wc -l | \
    awk '{print "Overlap CAGE enhancers:\t" $1 "\t" $1/'$totalCAGE'*100}'
awk '{if ($4 ~ /NET/) print $0}' tmp.out | wc -l | \
    awk '{print "Overlap NET-CAGE enhancers:\t" $1 "\t" $1/'$totalNET'*100}'


rm tmp.out no-map
