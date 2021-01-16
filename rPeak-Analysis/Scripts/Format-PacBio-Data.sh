#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

biosample=$1
genomeRef=~/Lab/Reference/Human/hg38/chromInfo.txt

if [ $biosample == "K562" ]
then
    accessions=("ENCFF546DOT" "ENCFF709YES") 
elif [ $biosample == "GM12878" ]
then
    accessions=("ENCFF247TLH" "ENCFF431IOE" "ENCFF520MMC" "ENCFF626GWM")
fi

for accession in ${accessions[@]}
do
    echo $accession

    if [ ! -f $accession.bam ]
    then
        wget https://www.encodeproject.org/files/$accession/@@download/$accession.bam
    fi
    
    samtools sort -@ 16 $accession.bam -o $accession.tmp.bam
done

echo "Merging bams ..."
samtools merge -@ 16 merge.bam *.tmp.bam

echo "Sorting and indexing bam ..."
samtools sort merge.bam -o $biosample-Merge.bam
samtools index $biosample-Merge.bam

echo "Separating strands ..."
samtools view -b -F 16 $biosample-Merge.bam -o pos.bam
samtools view -b -f 16 $biosample-Merge.bam -o neg.bam

echo "Converting to bigWigs ..."
~/bin/bamToBigWig pos.bam ~/Lab/Reference/Human/hg38/chromInfo.txt $biosample.full.plus.bigWig
~/bin/bamToBigWig neg.bam ~/Lab/Reference/Human/hg38/chromInfo.txt $biosample.full.minus.bigWig

###Extracting 5'ends

echo "Extracting 5'ends ..."
bedtools bamtobed -i pos.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_'$biosample'" \
    "\t" $5 "\t" $6}' > p.bed
bedtools bamtobed -i neg.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"_'$biosample'" \
    "\t" $5 "\t" $6}' > n.bed

awk '{print $1 "\t" $2 "\t" $2+1 "\t" $4 "\t" $5}' p.bed > tp5.bed
awk '{print $1 "\t" $3-1 "\t" $3 "\t" $4 "\t" $5}' n.bed > tn5.bed

echo "Converting to bedGraphs ..."
bedtools genomecov -bg -i tp5.bed -g $genomeRef | sort -k1,1 -k2,2n > p5.bg
bedtools genomecov -bg -i tn5.bed -g $genomeRef | sort -k1,1 -k2,2n > n5.bg

echo "Converting to bigWigs ..."
~/bin/bedGraphToBigWig p5.bg $genomeRef $biosample.5ends.plus.bigWig
~/bin/bedGraphToBigWig n5.bg $genomeRef $biosample.5ends.minus.bigWig

mv tp5.bed $biosample.5ends.plus.bed
mv tn5.bed $biosample.5ends.minus.bed
mv p.bed $biosample.full.plus.bed
mv n.bed $biosample.full.minus.bed

rm *tmp.bam pos.bam neg.bam p5.bg n5.bg merge.bam
cp $biosample.*.bigWig /home/moorej3/Public-HTML/TSS-Annotations/Track-Hub/hg38/data
