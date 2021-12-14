#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --array=71
#SBATCH --output=/home/moorej3/Lab/Job-Logs/jobid_%A_%a.output
#SBATCH --error=/home/moorej3/Lab/Job-Logs/jobid_%A_%a.error
#SBATCH --partition=12hours

#Jill E Moore
#Weng Lab
#UMass Medical School
#December 2021

mode=Basic
gencode=GENCODE31

scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
j=$SLURM_ARRAY_TASK_ID

masterList=~/Lab/ENCODE/RAMPAGE/RAMPAGE-List.txt
exp=$(awk '{if (NR == '$j') print $1}' $masterList)
outputDir=~/Lab/ENCODE/RAMPAGE/Read-Pair-Distance/$gencode-$mode-Batch-Results/$exp/
bedtools=~/bin/bedtools2/bin/bedtools
rpeaks=~/Lab/ENCODE/RAMPAGE/Genomic-Context/$gencode-$mode/hg38-rPeaks-genomic-context.txt

mkdir -p /tmp/moorej3/$SLURM_JOBID-$j
cd /tmp/moorej3/$SLURM_JOBID-$j

echo "Processing Reads..."
reads=/data/zusers/zhangx/projects/rampage/0_rampage_peak/$exp/rampage_peak_uniq/rampage_link.bed
awk '{if ($6 == "+") print $1 "\t" $2+1 "\t" $2+1 "\t" $1 "_" $2 "_" $3 \
    "\t" $5 "\t" $6 "\t" $3-$2; else print $1 "\t" $3 "\t" $3 "\t" \
     $1 "_" $2 "_" $3 "\t" $5 "\t" $6 "\t" $3-$2 }' $reads > 5reads.bed
     
awk '{if ($6 == "-") print $1 "\t" $2+1 "\t" $2+1 "\t" $1 "_" $2 "_" $3 \
    "\t" $5 "\t" $6 "\t" $3-$2; else print $1 "\t" $3 "\t" $3 "\t" \
     $1 "_" $2 "_" $3 "\t" $5 "\t" $6 "\t" $3-$2 }' $reads > 3reads.bed

annotation=~/Lab/Reference/Human/hg38/$gencode/Exon-Annotated.$mode.bed
$bedtools intersect -wo -s -a $annotation -b 3reads.bed > out3
  
#TSS
echo "Processing TSS..."
annotation=~/Lab/Reference/Human/hg38/$gencode/TSS.$mode.bed
grep TSS $rpeaks | \
    awk '{if ($NF == "same") print $0}' | sort -k1,1 -k2,2n > q
$bedtools intersect -s -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "TSS"}' | sort -u > out5
    
grep TSS $rpeaks | \
    awk '{if ($NF == "opposite") print $0}' | sort -k1,1 -k2,2n > q
$bedtools intersect -S -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "TSS-AS"}' | sort -u >> out5 

#TSS-Proximal
echo "Processing Proximal..."
annotation=~/Lab/Reference/Human/hg38/$gencode/TSS.$mode.4K.bed
grep Proximal $rpeaks | \
    awk '{if ($NF == "same") print $0}' | sort -k1,1 -k2,2n > q
$bedtools intersect -s -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "Proximal"}' | sort -u >> out5
    
grep Proximal $rpeaks | \
    awk '{if ($NF == "opposite") print $0}' | sort -k1,1 -k2,2n > q
$bedtools intersect -S -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "Proximal-AS"}' | sort -u >> out5 

#Introns
echo "Processing Introns..."
annotation=~/Lab/Reference/Human/hg38/$gencode/Transcripts.$mode.bed
grep Intron $rpeaks | sort -k1,1 -k2,2n > q
$bedtools intersect -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "Intron"}' | sort -u >> out5
    
#Intergenic
echo "Processing Intergenic..."
annotation=~/Lab/Reference/Human/hg38/$gencode/TSS.$mode.bed
grep Intergenic $rpeaks | sort -k1,1 -k2,2n > q
$bedtools closest -s -t all -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $15}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $0 "\t" \
    "Intergenic"}' | sort -u >> out5
    
#Exons
echo "Processing Exons..."
annotation=~/Lab/Reference/Human/hg38/$gencode/Exon-Annotated.$mode.bed
grep Exon $rpeaks | sort -k1,1 -k2,2n > q
$bedtools intersect -wo -a q -b $annotation | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $18 "\t" "Exon-"$20}' > bed
$bedtools intersect -wo -s -a bed -b 5reads.bed | awk '{print $1 "\t" $2 "\t" \
    $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $12 "\t" \
    $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $11}' | \
    sort -u >> out5

echo "Matching reads..."
python $scriptDir/match-supporting-reads.py out5 out3 > $exp-raw.txt

echo "Producing summary..."
python $scriptDir/filter-matched-reads.py $exp-raw.txt > $exp-summary.txt

echo "Compressing raw results..."
gzip $exp-raw.txt

mkdir -p $outputDir
mv $exp-raw.txt.gz $exp-summary.txt $outputDir
#mv * $outputDir
rm -r /tmp/moorej3/$SLURM_JOBID-$j
