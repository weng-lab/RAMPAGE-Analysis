#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
dataDir=~/Lab/ENCODE/RAMPAGE/Peaks
rPeakSummary=$dataDir/rPeaks-Filtered-Summary.txt
allPeaks=$dataDir/All-Peaks.bed
genomicContext=~/Lab/ENCODE/RAMPAGE/Genomic-Context/Full/All-RAMPAGE.genomic-context

cd ~/Lab/ENCODE/RAMPAGE/Boundary-Variation

bedtools intersect -wo -s -a $rPeakSummary -b $allPeaks > tmp.intersection
python ~/Projects/RAMPAGE/calculate-boundary-variation.py tmp.intersection | \
    sort -k1,1 > tmp.annotated
sort -k4,4 $genomicContext | awk '{print $4 "\t" $(NF-1)}' | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' tmp.annotated - | \
    paste tmp.annotated - > Boundary-Variation-Summary.txt

rm tmp.*
