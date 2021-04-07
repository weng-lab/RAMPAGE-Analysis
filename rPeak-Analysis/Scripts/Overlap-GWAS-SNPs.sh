#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

dataDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Input-Data/GWAS
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
gwasBed=$dataDir/GWAS-Jan-2019.bed
gwasCatalog=$dataDir/GWAS-Catalog-Jan-2019.txt
peaks=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
annotations=~/Lab/ENCODE/RAMPAGE/Read-Pair-Distance/RAMPAGE-Batch-Results/Filtered-All-Summaries.txt
geneAssignments=~/Lab/ENCODE/RAMPAGE/Read-Pair-Distance/RAMPAGE-Batch-Results/Read-Linked-Gene-Assignment.txt
genes=~/Lab/Reference/Human/hg38/GENCODE31/Genes.Basic.bed

sort -k1,1 $annotations > tmp.annotations

awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.annotations $peaks | sort -k4,4 | paste - tmp.annotations > tmp.bed
awk 'FNR==NR {x[$1];next} !($4 in x)' tmp.annotations $peaks | awk '{print $0 "\t" $4 "\t" "Exon+" "\t" \
    "Discard" "\t" "NA" "\t" "Discard" "\t" "NA" "\t" "NA"}' >> tmp.bed

echo "Summary of TSS classes ..."
bedtools intersect -wo -a tmp.bed -b $gwasBed > tmp.intersection
awk -F "\t" '{print $20 "\t" $4 "\t" $14}' tmp.intersection | sort -u > tmp.filter
python $scriptDir/count-groups.py tmp.filter 3

python $scriptDir/compare-gwas-gene-assignments.py $gwasCatalog \
    $geneAssignments tmp.intersection tmp.same tmp.other $genes
total=$(cat tmp.same tmp.other | awk -F "\t" '{print $4}' | sort -u | wc -l | awk '{print $1}')
other1=$(awk -F "\t" '{print $4}' tmp.other | sort -u | wc -l | awk '{print $1}')
other2=$(awk -F "\t" 'FNR==NR {x[$4];next} !($4 in x)' tmp.same tmp.other | awk -F "\t" '{print $4}' \
    | sort -u | wc -l | awk '{print $1}')

echo -e "\n"
echo "Summary of gene annotations..."
awk 'BEGIN{print "% unreported original study: "'$other1'/'$total'*100 "\n" \
    "% unreported all GWAS: " '$other2'/'$total'*100}'
mv tmp.intersection rPeak-GWAS-Intersection.txt
mv tmp.same GWAS-Matching-Gene-Assignment.txt
mv tmp.other GWAS-Different-Gene-Assignment.txt
rm tmp.*

