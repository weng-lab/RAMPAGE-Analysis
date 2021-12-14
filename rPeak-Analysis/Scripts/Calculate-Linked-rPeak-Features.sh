#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
peaks=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
annotations=~/Lab/ENCODE/RAMPAGE/Read-Pair-Distance/RAMPAGE-Batch-Results/Filtered-All-Summaries.txt
longReadDir=~/Lab/ENCODE/RAMPAGE/PacBio

sort -k1,1 $annotations > tmp.annotations

awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.annotations $peaks | sort -k4,4 | paste - tmp.annotations > tmp.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' tmp.bed > tmp.mini
bigWig=~/Lab/Reference/Human/hg38/Conservation/hg38.phastCons100way.bw
~/bin/bigWigAverageOverBed $bigWig tmp.mini tmp.out
sort -k1,1 tmp.out | paste tmp.bed tmp.out | awk -F "\t" '{print $4 "\t" $11 "\t" $14 "\t" $21}' > PhastCons-Summary.txt

counts=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks-Biosample-Counts.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' tmp.annotations $counts | sort -k1,1 | \
    paste tmp.annotations - | awk -F "\t" '{print $1 "\t" $5 "\t" $NF}' > Biosample-Count-Summary.txt

biosample=GM12878
awk '{if ($NF > 2) print $0}' ~/Lab/ENCODE/RAMPAGE/signal-output/ENCSR000AEI.txt | \
    awk 'FNR==NR {x[$1];next} ($4 in x)' - tmp.bed > tmp.tmp2
awk '{if ($6 == "+") print $0}' tmp.tmp2 > plus.bed
awk '{if ($6 == "-") print $0}' tmp.tmp2 > minus.bed
bedtools intersect -c -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > tmp.pacbio
bedtools intersect -c -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> tmp.pacbio
awk -F "\t" '{print $4 "\t" $14 "\t" $NF}' tmp.pacbio > PacBio-$biosample-Overlap-Summary.txt

bedtools intersect -wo -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > tmp.out1
bedtools intersect -wo -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed > tmp.out2
fullP=$longReadDir/$biosample.full.plus.bed
fullM=$longReadDir/$biosample.full.minus.bed
python $scriptDir/extract-pacbio-read-length.py tmp.out1 $fullP tmp.bed > tmp.output
python $scriptDir/extract-pacbio-read-length.py tmp.out2 $fullM tmp.bed >> tmp.output
mv tmp.output PacBio-$biosample-Length-Summary.txt

biosample=K562
awk '{if ($NF > 2) print $0}' ~/Lab/ENCODE/RAMPAGE/signal-output/ENCSR000AER.txt | \
    awk 'FNR==NR {x[$1];next} ($4 in x)' - tmp.bed > tmp.tmp2
awk '{if ($6 == "+") print $0}' tmp.tmp2 > plus.bed
awk '{if ($6 == "-") print $0}' tmp.tmp2 > minus.bed
bedtools intersect -c -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > tmp.pacbio
bedtools intersect -c -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> tmp.pacbio
awk -F "\t" '{print $4 "\t" $14 "\t" $NF}' tmp.pacbio > PacBio-$biosample-Overlap-Summary.txt

bedtools intersect -wo -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > tmp.out1
bedtools intersect -wo -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed > tmp.out2
fullP=$longReadDir/$biosample.full.plus.bed
fullM=$longReadDir/$biosample.full.minus.bed
python $scriptDir/extract-pacbio-read-length.py tmp.out1 $fullP tmp.bed > tmp.output
python $scriptDir/extract-pacbio-read-length.py tmp.out2 $fullM tmp.bed >> tmp.output
mv tmp.output PacBio-$biosample-Length-Summary.txt

rm tmp.*
