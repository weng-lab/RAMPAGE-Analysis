#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

genome=$1
links=$2
geneID=$3
peakID=$4

rPeaks=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
genes=~/Lab/Reference/Human/hg38/GENCODE31/Genes.Basic.bed
transcripts=~/Lab/Reference/Human/hg38/GENCODE31/TSS.Basic.bed

if [[ $genome == "hg19" ]]
then
    grep $peakID $rPeaks > tmp.peaks-hg38
    awk '{if ("'$geneID'" == $7) print $4}' $genes > tmp.gene
    awk 'FNR==NR {x[$1];next} ($7 in x)' tmp.gene $transcripts | \
        awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}'> tmp.tss-hg38
    ~/bin/liftOver tmp.peaks-hg38 ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain tmp.peaks no
    ~/bin/liftOver tmp.tss-hg38 ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain tmp.tss no
else
    grep $peakID $rPeaks > tmp.peaks
    awk '{if ("'$geneID'" == $7) print $4}' $genes > tmp.gene
    awk 'FNR==NR {x[$1];next} ($7 in x)' tmp.gene $transcripts | \
        awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}'> tmp.tss
fi

col=$(awk '{if (NR == 1) print NF}' $links)
if [[ $col -eq 3 ]]
then
    awk -F "," '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "Link-"NR}' $links  > tmp.links1
    awk '{print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8}' tmp.links1 > tmp.links2
else
    awk '{print $0 "\t" "Link-"NR}' $links > tmp.links1
    awk '{print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $9}' tmp.links1 > tmp.links2
fi

bedtools intersect -u -a tmp.links1 -b tmp.peaks > tmp.out-peaks
bedtools intersect -u -a tmp.links2 -b tmp.peaks >> tmp.out-peaks

bedtools intersect -u -a tmp.links1 -b tmp.tss > tmp.out-tss
bedtools intersect -u -a tmp.links2 -b tmp.tss >> tmp.out-tss

echo "Supporting links ..."
awk 'FNR==NR {x[$9];next} ($9 in x)' tmp.out-tss tmp.out-peaks

