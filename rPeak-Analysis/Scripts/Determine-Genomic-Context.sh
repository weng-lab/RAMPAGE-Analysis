#!/bin/bash

#Jill E. Moore
#Weng Lab
#UMass Medical School
#2019-2020

rPeaks=$1
output=$2
refDir=~/Lab/Reference/Human/hg38/GENCODE31

tss=$refDir/TSS.Basic.bed
exon=$refDir/Exon.Basic.bed
genes=$refDir/Transcripts.Basic.bed

awk '{print $1 "\t" $2-500 "\t" $2+500 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' $tss \
    > proximal500.bed

#TSS overlap
echo "Processing TSSs..."
bedtools intersect -u -a $rPeaks -b $tss > tmp.working1
bedtools intersect -u -s -a tmp.working1 -b $tss | \
    awk '{print $0 "\t" "TSS" "\t" "same"}' > tmp.out
bedtools intersect -v -s -a tmp.working1 -b $tss > tmp.working2
bedtools intersect -u -S -a tmp.working2 -b $tss | \
    awk '{print $0 "\t" "TSS" "\t" "opposite"}' >> tmp.out
bedtools intersect -v -a $rPeaks -b $tss > tmp.1

#TSS proximal
echo "Processing proximal..."
bedtools intersect -f 0.5 -u -a tmp.1 -b proximal500.bed > tmp.working1
bedtools intersect -f 0.5 -u -s -a tmp.working1 -b proximal500.bed | \
    awk '{print $0 "\t" "Proximal" "\t" "same"}' >> tmp.out
bedtools intersect -f 0.5 -v -s -a tmp.working1 -b proximal500.bed > tmp.working2
bedtools intersect -f 0.5 -u -S -a tmp.working2 -b proximal500.bed | \
    awk '{print $0 "\t" "Proximal" "\t" "opposite"}' >> tmp.out
bedtools intersect -f 0.5 -v -a tmp.1 -b proximal500.bed > tmp.2

#Exon overlap
echo "Processing exons..."
bedtools intersect -f 0.5 -u -a tmp.2 -b $exon > tmp.working1
bedtools intersect -f 0.5 -u -s -a tmp.working1 -b $exon | \
    awk '{print $0 "\t" "Exon" "\t" "same"}' >> tmp.out
bedtools intersect -f 0.5 -v -s -a tmp.working1 -b $exon > tmp.working2
bedtools intersect -f 0.5 -u -S -a tmp.working2 -b $exon | \
    awk '{print $0 "\t" "Exon" "\t" "opposite"}' >> tmp.out
bedtools intersect -f 0.5 -v -a tmp.2 -b $exon > tmp.1

#Intron overlap
echo "Processing introns..."
bedtools intersect -f 0.5 -u -a tmp.1 -b $genes > tmp.working1
bedtools intersect -f 0.5 -u -s -a tmp.working1 -b $genes | \
    awk '{print $0 "\t" "Intron" "\t" "same"}' >> tmp.out
bedtools intersect -f 0.5 -v -s -a tmp.working1 -b $genes > tmp.working2
bedtools intersect -f 0.5 -u -S -a tmp.working2 -b $genes | \
    awk '{print $0 "\t" "Intron" "\t" "opposite"}' >> tmp.out
bedtools intersect -f 0.5 -v -a tmp.1 -b $genes > tmp.working1

#Intergenic
echo "Processing intergenic..."
bedtools closest -a tmp.working1 -b $tss > tmp.working2
awk '{if ($6 == $15) print $4}' tmp.working2 | sort -u > tmp.match
awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.match tmp.working2 | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" \
    $9 "\t" "Intergenic" "\t" "same"}' | sort -u >> tmp.out
awk 'FNR==NR {x[$1];next} !($4 in x)' tmp.match tmp.working2 | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" \
    $9 "\t" "Intergenic" "\t" "opposite"}' | sort -u >> tmp.out    

mv tmp.out $output
rm tmp.* proximal500.bed
