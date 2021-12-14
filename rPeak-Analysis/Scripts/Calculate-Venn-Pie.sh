#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

biosample=$1

rampagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/RAMPAGE-$biosample-RPM-2.peaks.bed
cagePeaks=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/CAGE-$biosample-TPM-2.peaks.bed
groCap=~/Lab/ENCODE/RAMPAGE/GROcap/$biosample-GROcap-RAMPAGE.txt
longReadDir=~/Lab/ENCODE/RAMPAGE/PacBio

if [ $biosample == "K562" ]
then
    cutoff=0.06
elif [ $biosample == "GM12878" ]
then
    cutoff=0.08
fi

accession=$(awk '{if ($2 == "'$biosample'") print $1}' ~/Lab/ENCODE/RAMPAGE/RAMPAGE-List.txt)

mkdir -p ~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/VennPie/
cd ~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/VennPie/

#CAGE
bedtools intersect -u -s -F 0.25 -a $rampagePeaks -b $cagePeaks > cage
bedtools intersect -v -s -F 0.25 -a $rampagePeaks -b $cagePeaks > nocage

#PacBio
awk '{if ($6 == "+") print $0}' $rampagePeaks > plus.bed
awk '{if ($6 == "-") print $0}' $rampagePeaks > minus.bed
bedtools intersect -u -a plus.bed -b $longReadDir/$biosample.5ends.plus.bed > pacbio
bedtools intersect -u -a minus.bed -b $longReadDir/$biosample.5ends.minus.bed >> pacbio
awk 'FNR==NR {x[$0];next} !($0 in x)' pacbio $rampagePeaks > nopacbio

#GRO-cap
awk '{if ($2 > '$cutoff' || $2 < -1*'$cutoff') print $0}' $groCap > x
awk 'FNR==NR {x[$1];next} ($4 in x)' x $rampagePeaks > grocap
awk 'FNR==NR {x[$1];next} !($4 in x)' x $rampagePeaks > nogrocap

awk 'BEGIN{print "GRO-cap\tPacBio\tCAGE\t#"}' > output
awk 'FNR==NR {x[$4];next} ($4 in x)' grocap pacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage | wc -l | \
    awk '{print "yes\tyes\tyes\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage > 3-assay-support.txt

awk 'FNR==NR {x[$4];next} ($4 in x)' grocap pacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage | wc -l | \
    awk '{print "yes\tyes\tno\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage > 2-assay-support.txt

awk 'FNR==NR {x[$4];next} ($4 in x)' grocap nopacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage | wc -l | \
    awk '{print "yes\tno\tyes\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage >> 2-assay-support.txt
    
awk 'FNR==NR {x[$4];next} ($4 in x)' grocap nopacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage | wc -l | \
    awk '{print "yes\tno\tno\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage > 1-assay-support.txt
    
awk 'FNR==NR {x[$4];next} ($4 in x)' nogrocap pacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage | wc -l | \
    awk '{print "no\tyes\tyes\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage >> 2-assay-support.txt

awk 'FNR==NR {x[$4];next} ($4 in x)' nogrocap pacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage | wc -l | \
    awk '{print "no\tyes\tno\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage >> 1-assay-support.txt

awk 'FNR==NR {x[$4];next} ($4 in x)' nogrocap nopacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage | wc -l | \
    awk '{print "no\tno\tyes\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a cage >> 1-assay-support.txt
    
awk 'FNR==NR {x[$4];next} ($4 in x)' nogrocap nopacbio > a
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage | wc -l | \
    awk '{print "no\tno\tno\t" $1}' >> output
awk 'FNR==NR {x[$4];next} ($4 in x)' a nocage > 0-assay-support.txt

totalRAMPAGE=$(wc -l $rampagePeaks | awk '{print $1}')
cageSumm=$(wc -l cage | awk '{print $1/'$totalRAMPAGE'*100}')
pacbioSumm=$(wc -l pacbio | awk '{print $1/'$totalRAMPAGE'*100}')
grocapSumm=$(wc -l grocap | awk '{print $1/'$totalRAMPAGE'*100}')

echo -e "CAGE overlap:" $cageSumm
echo -e "PacBio overlap:" $pacbioSumm
echo -e "GROcap overlap:" $grocapSumm

support=(0 1 2 3)
for group in ${support[@]}
do
    awk 'FNR==NR {x[$4];next} ($1 in x)' $group-assay-support.txt \
        ~/Lab/ENCODE/RAMPAGE/signal-output/$accession.txt \
        | awk '{print $0 "\t" "'$group'"}' >> tmp.signal
    awk 'FNR==NR {x[$4];next} ($4 in x)' $group-assay-support.txt \
        ~/Lab/ENCODE/RAMPAGE/Genomic-Context/Full/All-RAMPAGE.genomic-context \
        | awk '{print $0 "\t" "'$group'"}' >> tmp.context
done

mv tmp.signal $biosample-Assay-Support-Signal.txt
mv tmp.context $biosample-Assay-Support-Genomic-Context.txt
mv output $biosample-VennPie-Summary.txt
rm cage nocage nogrocap grocap nopacbio pacbio a x *.bed
