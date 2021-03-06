#!/bin/bash
#Jill E. Moore
#Weng Lab
#UMass Medical School
#Updated March 2020

mode=$1
workingDir=$2
bigWig=$3
ccres=$workingDir/$mode.bed
echo $ccres

##Step 1 - Run first batch###
num=$(wc -l $ccres | awk '{print int($1/100)}')
sbatch --nodes 1 --array=1-$num%25 --mem=10G --time=00:30:00 \
    --output=/home/moorej3/Job-Logs/jobid_%A_%a.output \
    --error=/home/moorej3/Job-Logs/jobid_%A_%a.error \
    Aggregate-Signal-Batch-1.sh $mode $ccres $workingDir $bigWig

##Step 2 - Run second batch###
remainder=$(wc -l $ccres | awk '{print $1-'$num'*100}')
sbatch --nodes 1 --mem=10G --time=00:30:00 \
    --output=/home/moorej3/Job-Logs/jobid_%A.output \
    --error=/home/moorej3/Job-Logs/jobid_%A.error \
    Aggregate-Signal-Batch-2.sh $mode $remainder \
    $ccres $workingDir $bigWig

