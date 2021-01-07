#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

cd ~/Lab/ENCODE/RAMPAGE/FANTOM-CAGE

tpmMatrix=hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt
peaks=hg38_fair+new_CAGE_peaks_phase1and2.bed
peaks100=hg38_fair+new_CAGE_100bp_Summits.bed

awk '{print $1 "\t" $8-50 "\t" $8+50 "\t" $4 "\t" $5 "\t" $6 "\t" $7 \
    "\t" $8 "\t" $9}' $peaks > $peaks100

#Extract K562 peaks
echo "Processing K562 ..."
grep -v "##" $tpmMatrix | awk '{print $1 "\t" $563 "\t" $564 "\t" $565}' | \
    awk '{if (NR >3) print $1 "\t" ($2+$3+$4)/3}' > K562-Averaged-TPM.txt
awk '{if ($2 > 2) print $0}' K562-Averaged-TPM.txt > K562-TPM-2.txt
awk 'FNR==NR {x[$1];next} ($4 in x)' K562-TPM-2.txt $peaks > CAGE-K562-TPM-2.peaks.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' K562-TPM-2.txt $peaks100 > CAGE-K562-TPM-2.100bp.bed

cp CAGE-K562-TPM-2.peaks.bed CAGE-K562-TPM-2.100bp.bed ../K562-Comparison/

#Extract GM12878 peaks
echo "Processing GM12878 ..."
grep -v "##" $tpmMatrix | awk '{print $1 "\t" $171 "\t" $172 "\t" $173}' | \
    awk '{if (NR >3) print $1 "\t" ($2+$3+$4)/3}' > GM12878-Averaged-TPM.txt
awk '{if ($2 > 2) print $0}' GM12878-Averaged-TPM.txt > GM12878-TPM-2.txt
awk 'FNR==NR {x[$1];next} ($4 in x)' GM12878-TPM-2.txt $peaks > CAGE-GM12878-TPM-2.peaks.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' GM12878-TPM-2.txt $peaks100 > CAGE-GM12878-TPM-2.100bp.bed

cp CAGE-GM12878-TPM-2.peaks.bed CAGE-GM12878-TPM-2.100bp.bed ../GM12878-Comparison/
