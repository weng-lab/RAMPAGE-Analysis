#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

biosample=$1

rampage=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/RAMPAGE-$biosample-RPM-2.peaks.bed
rampageSummits=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks-Summits.bed

cage=~/Lab/ENCODE/RAMPAGE/$biosample-Comparison/CAGE-$biosample-TPM-2.peaks.bed
pacbio=~/Lab/ENCODE/RAMPAGE/PacBio/
bw1="*"$biosample"_GROcap_wTAP_plus.bigWig"
bw2="*"$biosample"_GROcap_wTAP_minus.bigWig"

mkdir -p ~/Lab/ENCODE/RAMPAGE/GROcap
cd ~/Lab/ENCODE/RAMPAGE/GROcap

awk 'FNR==NR {x[$4];next} ($4 in x)' $rampage $rampageSummits \
    | awk '{print $1 "\t" $10-1 "\t" $10 "\t" $4 "\t" $5 "\t" $6}' > tmp.rampage
~/bin/liftOver tmp.rampage ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain hg19.rampage tmp.no
awk '{print $1 "\t" $3-25 "\t" $3+25 "\t" $4 "\t" $5 "\t" $6}' hg19.rampage > RAMPAGE-$biosample.hg19-50bp.bed

awk '{print $1 "\t" $7 "\t" $8 "\t" $4 "\t" $5 "\t" $6}' $cage > tmp.cage
~/bin/liftOver tmp.cage ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain hg19.cage tmp.no
awk '{print $1 "\t" $3-25 "\t" $3+25 "\t" $4 "\t" $5 "\t" $6}' hg19.cage > CAGE-$biosample.hg19-50bp.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "0" "\t" "+"}' $pacbio/$biosample.5ends.plus.bed > tmp.pacbio
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "0" "\t" "-"}' $pacbio/$biosample.5ends.minus.bed >> tmp.pacbio
~/bin/liftOver tmp.pacbio ~/Lab/Reference/Human/hg38/hg38ToHg19.over.chain hg19.pacbio tmp.no
awk '{print $1 "\t" $3-25 "\t" $3+25 "\t" $4 "\t" $5 "\t" $6}' hg19.pacbio | grep -v "chrM" > PacBio-$biosample.hg19-50bp.bed

rm hg19.* tmp.*

assays=("PacBio" "CAGE" "RAMPAGE")
for assay in ${assays[@]}
do
    echo $assay
    bed=$assay-$biosample.hg19-50bp.bed    
    awk '{if ($6 == "+") print $1 "\t" $2 "\t" $3 "\t" $4}' $bed > plus
    awk '{if ($6 == "-") print $1 "\t" $2 "\t" $3 "\t" $4}' $bed > minus
    ~/bin/bigWigAverageOverBed $bw1 plus plus.tab
    ~/bin/bigWigAverageOverBed $bw2 minus minus.tab
    cat plus.tab minus.tab | awk '{print $1 "\t" $5 "\t" "'$assay'"}' > $biosample-GROcap-$assay.txt
done

~/bin/bigWigAverageOverBed $bw1 Random-500k.hg19-50bp.bed plus.tab
~/bin/bigWigAverageOverBed $bw1 Random-500k.hg19-50bp.bed minus.tab
paste plus.tab minus.tab | awk '{print $1 "\t" $5 "\t" "Random"}' > $biosample-GROcap-Random.txt

rm plus.tab minus.tab plus minus
