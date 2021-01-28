#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

dataDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Input-Data/GWAS
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
gwasBed=$dataDir/GWAS-Jan-2019.bed
gwasCatalog=$dataDir/GWAS-Catalog-Jan-2019.txt
peaks=~/Lab/ENCODE/RAMPAGE/Annotated-hg38-rPeaks.bed
geneAssignments=~/Lab/ENCODE/RAMPAGE/hg38-rPeak-Gene-Assignments.txt

echo "Summary of TSS classes ..."
bedtools intersect -wo -a $peaks -b $gwasBed > tmp.intersection
awk '{print $16 "\t" $4 "\t" $12}' tmp.intersection | sort -u > tmp.filter
python $scriptDir/count-groups.py tmp.filter 3

python $scriptDir/compare-gwas-gene-assignments.py $gwasCatalog \
    $geneAssignments tmp.intersection tmp.same tmp.other
total=$(cat tmp.same tmp.other | awk '{print $2}' | sort -u | wc -l | awk '{print $1}')
other1=$(awk '{print $2}' tmp.other | sort -u | wc -l | awk '{print $1}')
other2=$(awk 'FNR==NR {x[$2];next} !($2 in x)' tmp.same tmp.other | awk '{print $2}' \
    | sort -u | wc -l | awk '{print $1}')

echo -e "\n"
echo "Summary of gene annotations..."
awk 'BEGIN{print "% unreported original study: "'$other1'/'$total'*100 "\n" \
    "% unreported all GWAS: " '$other2'/'$total'*100}'
mv tmp.intersection rPeak-GWAS-Intersection.txt
mv tmp.same GWAS-Matching-Gene-Assignment.txt
mv tmp.other GWAS-Different-Gene-Assignment.txt
rm tmp.*

