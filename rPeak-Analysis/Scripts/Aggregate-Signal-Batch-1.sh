#!/bin/bash

mode=$1
ccres=$2
workingDir=$3
bigWig=$4
scriptDir=~/Projects/RAMPAGE

mkdir -p /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID
cd /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID

N=$(awk 'BEGIN {print '$SLURM_ARRAY_TASK_ID'*100}')

head -n $N $ccres | tail -n 100 > mini
awk '{if ($6 == "+") print $0}' mini > miniP
awk '{if ($6 == "-") print $0}' mini > miniN

rm -f header1 header2

for i in `seq -2000 1 2000`
do
    echo -e 0 "\t" 0 >> header1
    echo -e 0 "\t" 0 >> header2
done

f=miniP
i=$(wc -l miniP | awk '{print $1}')
for j in `seq 1 1 $i`
do
    echo $j
    chrom=$(awk '{if (NR == '$j') print $1}' $f)
    start=$(awk '{if (NR == '$j') print $2}' $f)
    stop=$(awk '{if (NR == '$j') print $3}' $f)
    python $scriptDir/per-bp.py $chrom $start $stop > mini.bed
    ~/bin/bigWigAverageOverBed $bigWig mini.bed out1.tab
    sort -k1,1g out1.tab | awk '{print $5}' > col1
    paste header1 col1 | awk '{print $1+1 "\t" $2+$3}' > tmp1
    mv tmp1 header1
done

f=miniN
i=$(wc -l miniN | awk '{print $1}')
for j in `seq 1 1 $i`
do
    echo $j
    chrom=$(awk '{if (NR == '$j') print $1}' $f)
    start=$(awk '{if (NR == '$j') print $2}' $f)
    stop=$(awk '{if (NR == '$j') print $3}' $f)
    python $scriptDir/per-bp.py $chrom $start $stop > mini.bed
    ~/bin/bigWigAverageOverBed $bigWig mini.bed out1.tab
    sort -k1,1g out1.tab | awk '{print $5}' > col1
    paste header2 col1 | awk '{print $1+1 "\t" $2+$3}' > tmp1
    mv tmp1 header2
done

sed '1!G;h;$!d' header2 > modHeader
paste header1 modHeader | awk '{print $1+$3 "\t" $2+$4}' > output

outputDir=$workingDir/Output/$mode
mkdir -p $outputDir

mv output $outputDir/agg-output.$mode.$SLURM_ARRAY_TASK_ID

rm -r /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID
