#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=0:5:00
#SBATCH --mem=30G
#SBATCH --array=1-22090%100
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A_%a.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A_%a.error
#SBATCH --partition=4hours

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

jid=$SLURM_ARRAY_TASK_ID

dataDir=~/Lab/ENCODE/RAMPAGE/Major-Minor-TSS
geneList=$dataDir/entropy-gene-pool
linkList=$dataDir/Read-Linked-Gene-Assignment.txt
expList=~/Lab/ENCODE/RAMPAGE/RAMPAGE-List.txt
sigDir=~/Lab/ENCODE/RAMPAGE/signal-output
outputDir=~/Lab/ENCODE/RAMPAGE/Major-Minor-TSS/Percentage-Runs
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts

mkdir -p $outputDir
mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

gene=$(awk '{if (NR == '$jid') print $1}' $geneList)
awk '{if ($2 == "'$gene'") print $0}' $linkList > tmp.linked
num=$(wc -l tmp.linked | awk '{print $1}')

rm -f tmp.total-biosample tmp.count
for i in `seq 1 1 $num`
do
    tss=$(awk '{if (NR == '$i') print $1}' tmp.linked)
    rm -f tmp.sig
    for j in `seq 1 1 115`
    do
        d=$(awk '{if (NR == '$j') print $1}' $expList)
        bio=$(awk '{if (NR == '$j') print $2}' $expList)
        grep $tss $sigDir/$d.txt | awk '{if ($NF > 2) print "'$bio'"}' >> tmp.sig
    done
    wc -l tmp.sig | awk '{print "'$tss'" "\t" $1}' >> tmp.count
    cat tmp.sig >> tmp.total-biosample
done

total=$(wc -l tmp.total-biosample | awk '{print $1}')
unique=$(sort -u tmp.total-biosample | wc -l | awk '{print $1}')
awk '{print "'$gene'" "\t" $1 "\t" $2 "\t" '$unique' "\t" $2/'$unique'*100 "\t" '$total' "\t" $2/'$total'*100}' tmp.count > tmp.output

mv tmp.output $outputDir/$gene-Summary.txt
rm -r /tmp/moorej3/$SLURM_JOBID-$jid
