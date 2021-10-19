#!/bin/bash

#Jill E. Moore
#Weng Lab
#UMass Medical School
#October 2021


mkdir -p /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID
cd /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID

j=$SLURM_ARRAY_TASK_ID
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
dataDir=~/Lab/ENCODE/RAMPAGE/ORF
mode=$1
inputList=$dataDir/$2
inputData=$dataDir/$3
inputSam=$dataDir/$4
bio=$5

m=$(awk '{if (NR == "'$j'") print $1}' $inputList)
grep $m $inputData > tmp.overlap

awk 'FNR==NR {x[$13];next} ($1"_'$bio'" in x)' tmp.overlap $inputSam \
    | awk '{print ">"$1 "\n" $10}' > tmp.fa

~/bin/ORFfinder -strand $mode -in tmp.fa -out tmp.orf

python $scriptDir/compile-orfs.py tmp.orf \
    | awk '{print "'$m'" "\t" $1 "\t" $2 "\t" $3 "\t" $4}'

rm -r /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID
