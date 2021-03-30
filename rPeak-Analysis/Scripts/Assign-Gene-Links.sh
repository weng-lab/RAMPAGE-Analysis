#!/bin/bash

#Jill E. Moore
#Weng Lab
#UMass Medical School
#March 2021

rpeaks=~/Lab/ENCODE/RAMPAGE/Genomic-Context/Full/All-RAMPAGE.genomic-context
workingDir=~/Lab/ENCODE/RAMPAGE/Read-Pair-Distance/RAMPAGE-Batch-Results
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts
geneDir=~/Lab/Reference/Human/hg38/GENCODE31

cd $workingDir

echo "Combining results ..."
cat ENCSR*/ENCSR*-summary.txt > Cat-All-Summaries.txt

awk '{print $1 "\t" $2-500 "\t" $2+500 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
     $geneDir/TSS.Basic.bed > proximal500.bed

grep TSS $rpeaks | bedtools intersect -wo -a - -b $geneDir/TSS.Basic.bed | \
    awk '{print $4 "\t" $(NF-1)}' | sort -u > tmp.genes
grep TSS $rpeaks | bedtools intersect -wo -a - -b $geneDir/TSS.Basic.bed | \
    awk '{print $4 "\t" $15}' | sort -u > tmp.transcripts

grep Proximal $rpeaks | bedtools intersect -wo -a - -b proximal500.bed | \
    awk '{print $4 "\t" $(NF-1)}' | sort -u >> tmp.genes
grep Proximal $rpeaks | bedtools intersect -wo -a - -b proximal500.bed | \
    awk '{print $4 "\t" $15}' | sort -u >> tmp.transcripts


grep Exon $rpeaks | bedtools intersect -wo -a - -b $geneDir/Exon-Annotated.Basic.bed | \
    awk '{print $4 "\t" $(NF-3)}' | sort -u >> tmp.genes
grep Exon $rpeaks | bedtools intersect -wo -a - -b $geneDir/Exon-Annotated.Basic.bed | \
    awk '{print $4 "\t" $(NF-4)}' | sort -u >> tmp.transcripts

grep Intron $rpeaks | bedtools intersect -wo -a - -b $geneDir/Genes.Basic.bed | \
    awk '{print $4 "\t" $(NF-5)}' | sort -u >> tmp.genes
grep Intron $rpeaks | bedtools intersect -wo -a - -b $geneDir/Transcripts.Basic.bed | \
    awk '{print $4 "\t" $15}' | sort -u >> tmp.transcripts


grep Intergenic $rpeaks | sort -k1,1 -k2,2n | bedtools closest -a - -b \
    $geneDir/TSS.Basic.bed | awk '{print $4 "\t" $(NF)}' >> tmp.genes
grep Intergenic $rpeaks | sort -k1,1 -k2,2n | bedtools closest -a - -b \
    $geneDir/TSS.Basic.bed | awk '{print $4 "\t" $15}' | sort -u >> tmp.transcripts

sort -k1,1 tmp.genes > Overlap-Gene-Assignment.txt
sort -k1,1 tmp.transcripts > Overlap-Transcript-Assignment.txt

awk '{if ($10 == "TSS" || $10 == "Proximal" || $10 == "Exon") print $0}' $rpeaks \
    | awk 'FNR==NR {x[$4];next} ($1 in x)' - Overlap-Transcript-Assignment.txt > tmp.candidate
awk 'FNR==NR {x[$4];next} ($2 in x)' $geneDir/Only-One-Exon.txt tmp.candidate | \
    awk '{print $1 "\t" "one-exon"}' | sort -u > tmp.candidate-list
awk '{if ($3-$2 > 500) print $0}' $geneDir/Exon-One-Basic.bed | \
    awk 'FNR==NR {x[$4];next} ($2 in x)' - tmp.candidate | \
    awk 'FNR==NR {x[$4];next} !($2 in x)' $geneDir/Only-One-Exon.txt - |  \
         awk '{print $1 "\t" "over-500bp"}' | sort -u >> tmp.candidate-list

sort -k1,1 tmp.candidate-list > Qualified-Candidate-TSS.txt

echo "Filtering results ..."
python $scriptDir/filter-all-summaries.py Cat-All-Summaries.txt \
    Overlap-Gene-Assignment.txt $geneDir/Transcript-Gene-Table.txt \
    Qualified-Candidate-TSS.txt > Filtered-All-Summaries.txt

sort -u gene-lists.txt > Read-Linked-Gene-Assignment.txt
rm gene-lists.txt
sort -u transcript-list.txt > Read-Linked-Transcript-Assignment.txt
rm transcript-list.txt

grep TSS $rpeaks | bedtools intersect -s -b - -a $geneDir/TSS.Basic.bed | \
    awk '{print $(NF)}' | sort -u | wc -l | awk '{print "Overlap gene coverage: " $1}'
awk '{if ($3 == "Candidate") print $1}' Filtered-All-Summaries.txt | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' - Overlap-Gene-Assignment.txt | \
    awk '{print $2}' - Read-Linked-Gene-Assignment.txt | sort -u | wc -l | \
     awk '{print "Linked gene coverage: " $1}'


awk -F "\t" '{if ($5 == "Verified GENCODE TSS ") print $0}' Filtered-All-Summaries.txt \
    | awk 'FNR==NR {x[$1];next} ($1 in x)' - Read-Linked-Gene-Assignment.txt  \
    | awk '{print $2}' | sort -u > tmp.x
wc -l tmp.x | awk '{print "\t" "verified GENCODE TSS: " $1}'
awk -F "\t" '{if ($5 == "Verified unannotated TSS ") print $0}' Filtered-All-Summaries.txt \
    | awk 'FNR==NR {x[$1];next} ($1 in x)' - Read-Linked-Gene-Assignment.txt  \
    | awk '{print $2}' | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.x - | sort -u > tmp.y
wc -l tmp.y | awk '{print "\t" "verified unannotated TSS: " $1}'
cat tmp.y >> tmp.x
awk '{if ($3 == "Candidate") print $1}' Filtered-All-Summaries.txt | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' - Overlap-Gene-Assignment.txt  \
    | awk '{print $2}' | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.x - | sort -u | wc -l | \
    awk '{print "\t" "candidate GENCODE TSS: " $1}'


grep TSS $rpeaks | bedtools intersect -s -b - -a $geneDir/TSS.Basic.bed | \
    awk '{print $4}' | sort -u | wc -l | awk '{print "Overlap TSS coverage: " $1}'
awk '{if ($3 == "Candidate") print $1}' Filtered-All-Summaries.txt | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' - Overlap-Transcript-Assignment.txt | \
    awk '{print $2}' - Read-Linked-Transcript-Assignment.txt | sort -u | wc -l | \
    awk '{print "Linked TSS coverage: " $1}'


awk -F "\t" '{if ($5 == "Verified GENCODE TSS ") print $0}' Filtered-All-Summaries.txt \
    | awk 'FNR==NR {x[$1];next} ($1 in x)' - Read-Linked-Transcript-Assignment.txt  \
    | awk '{print $2}' | sort -u > tmp.x
wc -l tmp.x | awk '{print "\t" "verified GENCODE TSS: " $1}'
awk -F "\t" '{if ($5 == "Verified unannotated TSS ") print $0}' Filtered-All-Summaries.txt \
    | awk 'FNR==NR {x[$1];next} ($1 in x)' - Read-Linked-Transcript-Assignment.txt  \
    | awk '{print $2}' | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.x - | sort -u > tmp.y
wc -l tmp.y | awk '{print "\t" "verified unannotated TSS: " $1}'
cat tmp.y >> tmp.x
awk '{if ($3 == "Candidate") print $1}' Filtered-All-Summaries.txt | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' - Overlap-Transcript-Assignment.txt  \
    | awk '{print $2}' | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.x - | sort -u | wc -l | \
    awk '{print "\t" "candidate GENCODE TSS: " $1}'

rm tmp.*
