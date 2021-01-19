#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

biosample=$1

dataDir=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/
rampagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/RAMPAGE-$biosample-RPM-2.peaks.bed
rampage100=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks-100bp.bed
cagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/CAGE-$biosample-TPM-2.peaks.bed
cage100=~/Lab/ENCODE/RAMPAGE/FANTOM-CAGE/hg38_fair+new_CAGE_100bp_Summits.bed
pacbioDir=~/Lab/ENCODE/RAMPAGE/PacBio
refDir=~/Lab/Reference/Human/hg38/GENCODE31
tss=$refDir/TSS.Basic.bed

mkdir -p $dataDir/Gene-Overlap
cd $dataDir/Gene-Overlap

awk '{print $1 "\t" $2-500 "\t" $2+500 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' $tss \
    > proximal500.bed

awk 'FNR==NR {x[$4];next} ($4 in x)' $rampagePeaks $rampage100 > tmp.rampage100
awk 'FNR==NR {x[$4];next} ($4 in x)' $cagePeaks $cage100 > tmp.cage100
awk '{print $1 "\t" $3-50 "\t" $2+50 "\t" $4 "\t" $5 "\t" "+"}' \
    $pacbioDir/$biosample.5ends.plus.bed | \
    awk '{if ($2 > 0) print $0; else print $1 "\t" 0 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' \
    > tmp.pacbio100
awk '{print $1 "\t" $3-50 "\t" $2+50 "\t" $4 "\t" $5 "\t" "-"}' \
    $pacbioDir/$biosample.5ends.minus.bed | \
    awk '{if ($2 > 0) print $0; else print $1 "\t" 0 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' \
    >> tmp.pacbio100

bedtools intersect -u -s -a proximal500.bed -b tmp.rampage100 > RAMPAGE-All-TSS.txt
bedtools intersect -u -s -a proximal500.bed -b tmp.cage100 > CAGE-All-TSS.txt
bedtools intersect -u -s -a proximal500.bed -b tmp.pacbio100 > PacBio-All-TSS.txt

awk '{print $NF}' RAMPAGE-All-TSS.txt | sort -u > a
awk '{print $NF}' CAGE-All-TSS.txt | sort -u > b
awk '{print $NF}' PacBio-All-TSS.txt | sort -u > c

awk 'FNR==NR {x[$0];next} ($0 in x)' a b > d
all=$(awk 'FNR==NR {x[$0];next} ($0 in x)' d c | wc -l | awk '{print $1}')
ab=$(awk 'FNR==NR {x[$0];next} !($0 in x)' c d | wc -l | awk '{print $1}')

awk 'FNR==NR {x[$0];next} ($0 in x)' a c > d
ac=$(awk 'FNR==NR {x[$0];next} !($0 in x)' b d | wc -l | awk '{print $1}')

awk 'FNR==NR {x[$0];next} ($0 in x)' b c > d
bc=$(awk 'FNR==NR {x[$0];next} !($0 in x)' a d | wc -l | awk '{print $1}')

awk 'FNR==NR {x[$0];next} !($0 in x)' b a > d
a=$(awk 'FNR==NR {x[$0];next} !($0 in x)' c d | wc -l | awk '{print $1}')

awk 'FNR==NR {x[$0];next} !($0 in x)' a b > d
b=$(awk 'FNR==NR {x[$0];next} !($0 in x)' c d | wc -l | awk '{print $1}')

awk 'FNR==NR {x[$0];next} !($0 in x)' a c > d
c=$(awk 'FNR==NR {x[$0];next} !($0 in x)' b d | wc -l | awk '{print $1}')

echo -e "R" "\t" "C" "\t" "P" "\t" "RC" "\t" "RP" "\t" "CP" "\t" "RCP"
echo -e $a "\t" $b "\t" $c "\t" $ab "\t" $ac "\t" $bc "\t" $all

rm a b c d proximal500.bed tmp.*
