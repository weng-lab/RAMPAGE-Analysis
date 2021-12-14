#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021


data=../RAMPAGE-List-Annotated.txt
signalDir=~/Lab/ENCODE/RAMPAGE/signal-output
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts

awk -F "\t" '{if ($3 == "tissue") print $0}' $data > tmp.list
q=$(wc -l tmp.list | awk '{print $1}')
exp=$(awk -F "\t" '{if (NR == 1) print $1}' tmp.list)
bio=$(awk -F "\t" '{if (NR == 1) print $2}' tmp.list)
awk 'BEGIN{print "rPeak" "\t" "'$bio'"}{print $1 "\t" $2}' $signalDir/$exp.txt > tmp.matrix
for j in `seq 2 1 $q`
do
    echo $j
    exp=$(awk -F "\t" '{if (NR == '$j') print $1}' tmp.list)
    bio=$(awk -F "\t" '{if (NR == '$j') print $2}' tmp.list)
    awk 'BEGIN{print "'$bio'"}{print $2}' $signalDir/$exp.txt > tmp.col
    paste tmp.matrix tmp.col > tmp.tmp
    mv tmp.tmp tmp.matrix
done
mv tmp.matrix Tissue-RPKM-Matrix.txt

cp $data tmp.list
q=$(wc -l tmp.list | awk '{print $1}')
exp=$(awk -F "\t" '{if (NR == 1) print $1}' tmp.list)
bio=$(awk -F "\t" '{if (NR == 1) print $2}' tmp.list)
awk 'BEGIN{print "rPeak" "\t" "'$bio'"}{print $1 "\t" $2}' $signalDir/$exp.txt > tmp.matrix
for j in `seq 2 1 $q`
do
    echo $j
    exp=$(awk -F "\t" '{if (NR == '$j') print $1}' tmp.list)
    bio=$(awk -F "\t" '{if (NR == '$j') print $2}' tmp.list)
    awk 'BEGIN{print "'$bio'"}{print $2}' $signalDir/$exp.txt > tmp.col
    paste tmp.matrix tmp.col > tmp.tmp
    mv tmp.tmp tmp.matrix
done
mv tmp.matrix All-Biosample-RPKM-Matrix.txt

python $scriptDir/run-umap.py Tissue-RPKM-Matrix.txt 10 > UMAP-Tissue-10.txt
python $scriptDir/run-umap.py All-Biosample-RPKM-Matrix.txt 10 > UMAP-All-10.txt
