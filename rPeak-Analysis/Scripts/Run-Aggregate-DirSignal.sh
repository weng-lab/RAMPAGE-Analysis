#!/bin/bash
#Jill E. Moore
#Weng Lab
#UMass Medical School
#December 2021

mode=$1
workingDir=$2
bigWigP=$3
bigWigN=$4
ccres=$workingDir/$mode.bed
echo $ccres

##Step 1 - Run first batch###
num=$(wc -l $ccres | awk '{print int($1/100)}')
sbatch --nodes 1 --array=1-$num%25 --mem=10G --time=00:30:00 \
    --output=/home/moorej3/Job-Logs/jobid_%A_%a.output \
    --error=/home/moorej3/Job-Logs/jobid_%A_%a.error \
    Aggregate-DirSignal-Batch-1.sh $mode $ccres $workingDir $bigWigP $bigWigN

##Step 2 - Run second batch###
remainder=$(wc -l $ccres | awk '{print $1-'$num'*100}')
sbatch --nodes 1 --mem=10G --time=00:30:00 \
    --output=/home/moorej3/Job-Logs/jobid_%A.output \
    --error=/home/moorej3/Job-Logs/jobid_%A.error \
    Aggregate-DirSignal-Batch-2.sh $mode $remainder \
    $ccres $workingDir $bigWigP $bigWigN

