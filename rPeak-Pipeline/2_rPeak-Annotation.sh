#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

peakDir=~/Lab/ENCODE/RAMPAGE/Peaks
scriptDir=~/Projects/RAMPAGE/

cd $peakDir

cat *RPM2-70.RNA-Filtered.bed > tmp.all

awk '{if ($6 == "+") print $0}' tmp.all > tmp.plus
awk '{if ($6 == "-") print $0}' tmp.all > tmp.minus

sort -k1,1 -k2,2n tmp.plus > tmp.splus
sort -k1,1 -k2,2n tmp.minus > tmp.sminus

rm -f tmp.rPlus
num=$(wc -l tmp.splus | awk '{print $1}')
echo -e "Merging plus peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i tmp.splus -c 4,5 -o collapse,collapse > tmp.merge-plus
    python $scriptDir/pick-best-peak.py tmp.merge-plus > tmp.peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.peak-list tmp.splus >> tmp.rPlus
    bedtools intersect -v -a tmp.splus -b tmp.rPlus > tmp.remaining
    mv tmp.remaining tmp.splus
    num=$(wc -l tmp.splus | awk '{print $1}')
done

rm -f tmp.rMinus
num=$(wc -l tmp.sminus | awk '{print $1}')
echo -e "Merging plus peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i tmp.sminus -c 4,5 -o collapse,collapse > tmp.merge-minus
    python $scriptDir/pick-best-peak.py tmp.merge-minus > tmp.peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.peak-list tmp.sminus >> tmp.rMinus
    bedtools intersect -v -a tmp.sminus -b tmp.rMinus > tmp.remaining
    mv tmp.remaining tmp.sminus
    num=$(wc -l tmp.sminus | awk '{print $1}')
done

cat tmp.rMinus tmp.rPlus | sort -k1,1 -k2,2n | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6 "\t" $7 \
    "\t" $8 "\t" $9}' > rPeaks-Raw.bed

bedtools intersect -s -c -a rPeaks-Raw.bed -b tmp.all > tmp.intersect
python $scriptDir/filter-rPeaks.py tmp.intersect > rPeaks-Filtered.bed
python $scriptDir/accession-rPeaks.py rPeaks-Filtered.bed hg38 rPeak \
    | sort -k1,1 -k2,2n | awk '{if ($1 !~ /_/ && $1 != "chrM") print $0}' \
    > rPeaks-Filtered-Summary.txt
awk '{print $1 "\t" $2 "\t" $3 "\t" $NF "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' \
    rPeaks-filtered-Summary.txt > hg38-rPeaks.bed

~/bin/bedToBigBed -type=bed9 hg38-rPeaks.bed ~/Lab/Reference/Human/hg38/chromInfo.txt hg38-rPeaks.bigBed

rm tmp*

