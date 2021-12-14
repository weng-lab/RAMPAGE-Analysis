#!/bin/bash

#Jill E. Moore
#Weng Lab
#UMass Medical School
#December 2021

rPeaks=$1
output=$2
gencode=$3
mode=$4
refDir=~/Lab/Reference/Human/hg38/$gencode

tss=$refDir/TSS.$mode.bed
exon=$refDir/Exon-Annotated.$mode.bed
genes=$refDir/Transcripts.$mode.bed

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

#Calculate enrichment
echo "Calculating enrichment..."
#echo -e "\n"
#echo -e "Group\tObserved\tExpected\tEnrichment"
#rPeakTotal=$(awk '{sum += $3-$2}END{print sum}' tmp.out)
#regions=("TSS" "Proximal" "Exon" "Intron" "Intergenic")
#for r in ${regions[@]}
#do
#    expected=$(awk '{if ($1 == "'$r'") print $3}' $refDir/Baseline-Genomic-Coverage.txt)
#    awk '{if ($10 == "'$r'") sum += $3-$2}END{print "'$r'" "\t" \
#        sum/'$rPeakTotal' "\t" '$expected' "\t" (sum/'$rPeakTotal')/'$expected'}' \
#        tmp.out
#done


mv tmp.out $output

rm tmp.* proximal500.bed
