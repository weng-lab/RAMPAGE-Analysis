#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

biosample=$1

rampagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/RAMPAGE-$biosample-RPM-2.peaks.bed
cagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/CAGE-$biosample-TPM-2.peaks.bed
groPeaks=~/Lab/ENCODE/RAMPAGE/GRO-PRO/"TSS_All_"$biosample.hg38.bed
groCapDir=~/Lab/ENCODE/RAMPAGE/GROcap
longReadDir=~/Lab/ENCODE/RAMPAGE/PacBio

if [ $biosample == "K562" ]
then
    cutoff=0.06
elif [ $biosample == "GM12878" ]
then
    cutoff=0.08
fi

echo -e "\t RAMPAGE \t CAGE \t PacBio \t GROcap"

totalRAMPAGE=$(wc -l $rampagePeaks | awk '{print $1}')
totalCAGE=$(wc -l $cagePeaks | awk '{print $1}')
totalPacBio=$(cat $longReadDir/$biosample.5ends.minus.bed $longReadDir/$biosample.5ends.plus.bed | wc -l | awk '{print $1}')

#RAMPAGE
rampageCAGE=$(bedtools intersect -u -s -F 0.25 -a $rampagePeaks -b $cagePeaks | wc -l | \
    awk '{print $1/'$totalRAMPAGE'*100","$1}')
awk '{if ($6 == "+") print $0}' $rampagePeaks > plus.bed
awk '{if ($6 == "-") print $0}' $rampagePeaks > minus.bed
bedtools intersect -u -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> pacbio
rampagePACBIO=$(wc -l pacbio | awk '{print $1/'$totalRAMPAGE'*100","$1}')
rampageGRO=$(awk '{if ($2 > '$cutoff' || $2 < -1*'$cutoff') print $0}' $groCapDir/$biosample-GROcap-RAMPAGE.txt | \
    wc -l | awk '{print $1/'$totalRAMPAGE'*100","$1}')
echo -e "RAMPAGE \t --- \t" $rampageCAGE "\t" $rampagePACBIO "\t" $rampageGRO

#CAGE
cageRAMPAGE=$(bedtools intersect -u -s -f 0.25 -a $cagePeaks -b $rampagePeaks | wc -l | \
    awk '{print $1/'$totalCAGE'*100","$1}')
awk '{if ($6 == "+") print $0}' $cagePeaks > plus.bed
awk '{if ($6 == "-") print $0}' $cagePeaks > minus.bed
bedtools intersect -u -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> pacbio
cagePACBIO=$(wc -l pacbio | awk '{print $1/'$totalCAGE'*100","$1}')
cageGRO=$(awk '{if ($2 > '$cutoff' || $2 < -1*'$cutoff') print $0}' $groCapDir/$biosample-GROcap-CAGE.txt | \
    wc -l | awk '{print $1/'$totalCAGE'*100","$1}')
echo -e "CAGE \t" $cageRAMPAGE "\t" "---" "\t" $cagePACBIO "\t" $cageGRO

#PacBio
awk '{if ($6 == "+") print $0}' $rampagePeaks > plus.bed
awk '{if ($6 == "-") print $0}' $rampagePeaks > minus.bed
bedtools intersect -u -b plus.bed -a $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -b minus.bed -a $longReadDir/$biosample.5ends.minus.bed >> pacbio
pacbioRAMPAGE=$(wc -l pacbio | awk '{print $1/'$totalPacBio'*100","$1}')
awk '{if ($6 == "+") print $0}' $cagePeaks > plus.bed
awk '{if ($6 == "-") print $0}' $cagePeaks > minus.bed
bedtools intersect -u -b plus.bed -a $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -b minus.bed -a $longReadDir/$biosample.5ends.minus.bed >> pacbio
pacbioCAGE=$(wc -l pacbio | awk '{print $1/'$totalPacBio'*100","$1}')
pacbioGRO=$(awk '{if ($2 > '$cutoff' || $2 < -1*'$cutoff') print $0}' $groCapDir/$biosample-GROcap-PacBio.txt | \
    wc -l | awk '{print $1/'$totalPacBio'*100","$1}')
echo -e "PacBio \t" $pacbioRAMPAGE "\t" $pacbioCAGE "\t" "---" "\t" $pacbioGRO

#GroPeak
totalGRO=$(wc -l $groPeaks | awk '{print $1}')
groRAMPAGE=$(bedtools intersect -u -s -a $groPeaks -b $rampagePeaks | wc -l | \
    awk '{print $1/'$totalGRO'*100","$1}')
groCAGE=$(bedtools intersect -u -s -a $groPeaks -b $cagePeaks | wc -l | \
    awk '{print $1/'$totalGRO'*100","$1}')

awk '{if ($6 == "+") print $0}' $groPeaks > plus.bed 
awk '{if ($6 == "-") print $0}' $groPeaks > minus.bed

bedtools intersect -u -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> pacbio
groPACBIO=$(wc -l pacbio | awk '{print $1/'$totalGRO'*100","$1}')

echo -e "GROcap \t" $groRAMPAGE "\t" $groCAGE "\t" $groPACBIO  "\t" "---"

rm plus.bed minus.bed pacbio
