#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

dataDir=~/Lab/ENCODE/RAMPAGE
rampage=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
cage=~/Lab/ENCODE/RAMPAGE/FANTOM-CAGE/hg38_fair+new_CAGE_peaks_phase1and2.bed

rampageTotal=$(wc -l $rampage | awk '{print $1}')
cageTotal=$(wc -l $cage | awk '{print $1}')

rOverlap=$(bedtools intersect -u -a $rampage -b $cage | wc -l | awk '{print $1}')
cOverlap=$(bedtools intersect -u -a $cage -b $rampage | wc -l | awk '{print $1}')

echo -e "Out of" $rampageTotal "RAMPAGE peaks," $rOverlap "overlap CAGE peaks"
echo -e "Out of" $cageTotal "CAGE peaks," $cOverlap "overlap RAMPAGE peaks"
