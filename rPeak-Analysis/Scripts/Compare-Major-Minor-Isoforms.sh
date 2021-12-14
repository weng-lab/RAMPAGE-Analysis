#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=0:5:00
#SBATCH --mem=30G
#SBATCH --array=1-5%100
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A_%a.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A_%a.error
#SBATCH --partition=4hours

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

jid=$SLURM_ARRAY_TASK_ID

dataDir=~/Lab/ENCODE/RAMPAGE/Major-Minor-TSS
novelList=$dataDir/working-pool-novel
verifiedList=$dataDir/working-pool-verified
expList=~/Lab/ENCODE/RAMPAGE/RAMPAGE-List.txt
sigDir=~/Lab/ENCODE/RAMPAGE/signal-output
outputDir=~/Lab/ENCODE/RAMPAGE/Major-Minor-TSS/Batch-Runs

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

novelTSS=$(awk '{if (NR == '$jid') print $1}' $novelList)
gene=$(awk '{if (NR == '$jid') print $2}' $novelList)

rm -f tmp.sig-novel
for j in `seq 1 1 115`
do
    d=$(awk '{if (NR == '$j') print $1}' $expList)
    bio=$(awk '{if (NR == '$j') print $2}' $expList)
    grep $novelTSS $sigDir/$d.txt | awk '{print "'$bio'" "\t" $1 "\t" $NF}' >> tmp.sig-novel
done


awk '{if ($2 == "'$gene'") print $0}' $verifiedList > tmp.verified
num=$(wc -l tmp.verified | awk '{print $1}')

for i in `seq 1 1 $num`
do
    verifiedTSS=$(awk '{if (NR == '$i') print $1}' tmp.verified)
    rm -f tmp.sig-verified
    for j in `seq 1 1 115`
    do
        d=$(awk '{if (NR == '$j') print $1}' $expList)
        bio=$(awk '{if (NR == '$j') print $2}' $expList)
        grep $verifiedTSS $sigDir/$d.txt | awk '{print "'$bio'" "\t" $1 "\t" $NF}' >> tmp.sig-verified
    done
    paste tmp.sig-verified tmp.sig-novel > tmp.paste
    
    verifiedNum=$(awk '{if ($3 > 2) x += 1}END{print x}' tmp.paste)
    novelNum=$(awk '{if ($6 > 2) x += 1}END{print x}' tmp.paste)

    awk '{if ($3 > 2) print $1 "\t" $2}' tmp.paste > tmp.verified-tissue
    awk '{if ($6 > 2) print $4 "\t" $5}' tmp.paste > tmp.novel-tissue

    novelHigher=$(awk '{if ($6 > $3) x += 1}END{print x}' tmp.paste)
    awk '{if ($6 > $3) print $4 "\t" $5 "\t" ($6+0.01)/($3+0.01) "\t" $3 "\t" $6}' tmp.paste > tmp.novel-higher-tissue
    
    mv tmp.verified-tissue $outputDir/$gene-$verifiedTSS-$novelTSS.Verified-Tissue.txt
    mv tmp.novel-tissue $outputDir/$gene-$verifiedTSS-$novelTSS.Novel-Tissue.txt
    mv tmp.novel-higher-tissue $outputDir/$gene-$verifiedTSS-$novelTSS.Novel-Higher-Tissue.txt
    
done

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
