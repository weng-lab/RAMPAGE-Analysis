#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

geneAssignments=~/Lab/ENCODE/RAMPAGE/hg38-rPeak-Gene-Assignments.txt
tss=~/Lab/Reference/Human/hg38/GENCODE31/TSS.Basic.bed
rpeaks=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks.bed
rpeak100=~/Lab/ENCODE/RAMPAGE/hg38-rPeaks-100bp.bed
annotationDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Input-Data

bedtools intersect -s -u -a $tss -b $rpeaks > tmp.overlap
awk 'FNR==NR {x[$7];next} !($7 in x)' tmp.overlap $tss > tmp.no-tss
awk 'FNR==NR {x[$2];next} ($7 in x)' $geneAssignments tmp.no-tss > tmp.ver-tss
awk 'FNR==NR {x[$7];next} ($2 in x)' tmp.ver-tss $geneAssignments > tmp.ver 

awk '{print $1 "\t" $2-50 "\t" $3+50 "\t" $4 "\t" $5 "\t" $6}' tmp.ver-tss > \
    GENCODE-Gene.Comparison-Set.100bp.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.ver $rpeak100 > Verified-TSS.Comparison-Set.100bp.bed

k562Peaks=~/Lab/ENCODE/RAMPAGE/K562-Comparison/RAMPAGE-K562-RPM-2.peaks.bed
awk 'FNR==NR {x[$4];next} ($1 in x)' $k562Peaks tmp.ver > tmp.ver-k562
awk 'FNR==NR {x[$2];next} ($7 in x)' tmp.ver-k562 tmp.ver-tss | \
    awk '{print $1 "\t" $2-50 "\t" $3+50 "\t" $4 "\t" $5 "\t" $6}' > GENCODE-Gene.K562-Set.100bp.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.ver-k562 $rpeak100 > Verified-TSS.K562-Set.100bp.bed

rm tmp.*

totalGENCODE=$(wc -l GENCODE-Gene.Comparison-Set.100bp.bed | awk '{print $1}')
totalRAMPAGE=$(wc -l Verified-TSS.Comparison-Set.100bp.bed | awk '{print $1}')

totalK562GENCODE=$(wc -l GENCODE-Gene.K562-Set.100bp.bed | awk '{print $1}')
totalK562RAMPAGE=$(wc -l Verified-TSS.K562-Set.100bp.bed | awk '{print $1}')

echo "Mapping TSSs to mm10 genome ..."
#chain=~/Lab/Reference/Human/hg38/hg38ToMm10.over.chain
#~/bin/liftOver -minMatch=0.5 Verified-TSS.Comparison-Set.100bp.bed $chain \
#    Verified-TSS.Comparison-Set.100bp.mm10.bed no
#~/bin/liftOver -minMatch=0.5 GENCODE-Gene.Comparison-Set.100bp.bed $chain \
#    GENCODE-Gene.Comparison-Set.100bp.mm10.bed no

echo "Calculating phastCons conservation ..."
#bigWig=~/Lab/Reference/Human/hg38/Conservation/hg38.phastCons100way.bw
#~/bin/bigWigAverageOverBed $bigWig Verified-TSS.Comparison-Set.100bp.bed \
#    Verified-TSS.Comparison-Set.phastCons.txt
#~/bin/bigWigAverageOverBed $bigWig GENCODE-Gene.Comparison-Set.100bp.bed \
#    GENCODE-Gene.Comparison-Set.phastCons.txt

echo "Overlap with all cCREs ..."
ccres=$annotationDir/GRCh38-cCREs.bed
bedtools intersect -u -a Verified-TSS.Comparison-Set.100bp.bed -b $ccres | wc -l \
    | awk '{print "RAMPAGE verified TSSs: " $1/'$totalRAMPAGE'*100  "\t" $1 "\t" '$totalRAMPAGE'-$1}'
bedtools intersect -u -a GENCODE-Gene.Comparison-Set.100bp.bed -b $ccres | wc -l \
    | awk '{print "GENCODE TSSs: " $1/'$totalGENCODE'*100  "\t" $1 "\t" '$totalGENCODE'-$1}'

echo "Overlap with all eQTLs ..."
eqtls=$annotationDir/GTEx-eQTLs.bed
bedtools intersect -u -a Verified-TSS.Comparison-Set.100bp.bed -b $eqtls | wc -l \
    | awk '{print "RAMPAGE verified TSSs: " $1/'$totalRAMPAGE'*100  "\t" $1 "\t" '$totalRAMPAGE'-$1}'
bedtools intersect -u -a GENCODE-Gene.Comparison-Set.100bp.bed -b $eqtls | wc -l \
    | awk '{print "GENCODE TSSs: " $1/'$totalGENCODE'*100  "\t" $1 "\t" '$totalGENCODE'-$1}'


echo "Overlap with K562 cCREs ..."
ccresK562=$annotationDir/K562-cCREs.bed
bedtools intersect -u -a Verified-TSS.K562-Set.100bp.bed -b $ccresK562 | wc -l \
    | awk '{print "RAMPAGE verified TSSs: " $1/'$totalK562RAMPAGE'*100  "\t" $1 "\t" '$totalK562RAMPAGE'-$1}'
bedtools intersect -u -a GENCODE-Gene.K562-Set.100bp.bed -b $ccresK562 | wc -l \
    | awk '{print "GENCODE TSSs: " $1/'$totalK562GENCODE'*100  "\t" $1 "\t" '$totalK562GENCODE'-$1}'

echo -e "\n"

echo "Overlap with SuRE peaks ..."
surePeaks=$annotationDir/SuRE-Peaks.hg38.bed
bedtools intersect -u -a Verified-TSS.K562-Set.100bp.bed -b $surePeaks | wc -l \
    | awk '{print "RAMPAGE verified TSSs: " $1/'$totalK562RAMPAGE'*100  "\t" $1 "\t" '$totalK562RAMPAGE'-$1}'
bedtools intersect -u -a GENCODE-Gene.K562-Set.100bp.bed -b $surePeaks | wc -l \
    | awk '{print "GENCODE TSSs: " $1/'$totalK562GENCODE'*100  "\t" $1 "\t" '$totalK562GENCODE'-$1}'

echo "Overlap with 5'ends ..."
pacBioPlus=/home/moorej3/Lab/ENCODE/RAMPAGE/PacBio/K562.5ends.plus.bed
pacBioMinus=/home/moorej3/Lab/ENCODE/RAMPAGE/PacBio/K562.5ends.minus.bed

awk '{if ($6 == "+") print $0}' Verified-TSS.K562-Set.100bp.bed > plus.bed
awk '{if ($6 == "-") print $0}' Verified-TSS.K562-Set.100bp.bed > minus.bed

bedtools intersect -c -a plus.bed -b $pacBioPlus | awk '{print $0 "\t" "RAMPAGE"}' > PacBio-Summary.txt
bedtools intersect -c -a minus.bed -b $pacBioMinus | awk '{print $0 "\t" "RAMPAGE"}' >> PacBio-Summary.txt

awk '{if ($6 == "+") print $0}' GENCODE-Gene.K562-Set.100bp.bed > plus.bed
awk '{if ($6 == "-") print $0}' GENCODE-Gene.K562-Set.100bp.bed > minus.bed

bedtools intersect -c -a plus.bed -b $pacBioPlus | awk '{print $0 "\t" "GENCODE"}' >> PacBio-Summary.txt
bedtools intersect -c -a minus.bed -b $pacBioMinus | awk '{print $0 "\t" "GENCODE"}' >> PacBio-Summary.txt

