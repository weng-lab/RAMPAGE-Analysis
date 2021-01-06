#!/bin/bash

#Jill E Moore
#Weng Lab
#UMass Medical School
#January 2021

file=RAMPAGE-RNA-Match-List.txt

peakDir=~/Lab/ENCODE/RAMPAGE/Peaks
signalDir=/data/projects/encode/data
dataDir=/data/zusers/zhangx/projects/rampage/0_rampage_peak/

k=$(wc -l $file | awk '{print $1}')
rm -f tmp
for j in `seq 1 1 $k`;
do

    fields=$(awk '{if (NR == '$j') print NF}' $file) 
    exp=$(awk '{if (NR == '$j') print $1}' $file)
    reads=$(awk '{if (NR == '$j') print $(NF-1)}' $file)
    expRNA=$(awk '{if (NR == '$j') print $3}' $file)
    
    echo $exp
    
    awk '{if ($6 == "+") print $1 "\t" $2 "\t" $3 "\t" $4}' \
        $peakDir/$exp-RPM2-70.bed | grep -v "EBV" > plus.bed
    awk '{if ($6 == "-") print $1 "\t" $2 "\t" $3 "\t" $4}' \
        $peakDir/$exp-RPM2-70.bed | grep -v "EBV" > minus.bed
    
    if [ "$fields" -eq 9 ] #when two replicates are available
    then
        echo "9 fields"
        bw1P=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$4".bigWig"}' $file)
        bw2P=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$5".bigWig"}' $file)
        bw1N=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$6".bigWig"}' $file)
        bw2N=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$7".bigWig"}' $file)   
    
        ~/bin/bigWigAverageOverBed $bw1P plus.bed plus1.out
        ~/bin/bigWigAverageOverBed $bw2P plus.bed plus2.out 
        ~/bin/bigWigAverageOverBed $bw1N minus.bed minus1.out
        ~/bin/bigWigAverageOverBed $bw2N minus.bed minus2.out

        paste plus1.out plus2.out | awk '{print $1 "\t" ($5 + $11)/'$reads'*1000000}' \
            > tmp.rna-summary
        paste minus1.out minus2.out | awk '{print $1 "\t" ($5 + $11)/'$reads'*1000000}' \
            >> tmp.rna-summary
    
        sort -k1,1 tmp.rna-summary > tmp.sort
        mv tmp.sort tmp.rna-summary
        
    else #when only one replicate is available
        echo "7 fields"
        bw1P=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$4".bigWig"}' $file)
        bw1N=$(awk '{if (NR == '$j') print "'$signalDir'/"$3"/"$5".bigWig"}' $file)

        ~/bin/bigWigAverageOverBed $bw1P plus.bed plus1.out
        ~/bin/bigWigAverageOverBed $bw1N minus.bed minus1.out
        
        awk '{print $1 "\t" $5/'$reads'*1000000}' plus1.out minus1.out \
            | sort -k1,1 > tmp.rna-summary
    fi
        
    bwP=$dataDir/bw/$exp"_rampage_plus.bw"
    bwN=$dataDir/bw/$exp"_rampage_minus.bw"
    reads=$(awk '{if (NR == '$j') print $2}' $file)
    
    ~/bin/bigWigAverageOverBed $bwP plus.bed plus.out
    ~/bin/bigWigAverageOverBed $bwN minus.bed minus.out
    
    awk '{print $1 "\t" $5/'$reads'*1000000}' plus.out minus.out \
        | sort -k1,1 > tmp.rampage-summary

    ~/Projects/RAMPAGE/Determine-Genomic-Context.sh $peakDir/$exp-RPM2-70.bed
    sort -k4,4 genomic-context-orientation | grep -v "EBV" > tmp.peaks
    paste tmp.peaks tmp.rna-summary tmp.rampage-summary \
        | awk '{print $4 "\t" $10 "\t" $13 "\t" $15}' > $exp.RAMPAGE-RNA.Summary.txt
    awk '{if ( log($4)/log(10) > log($3)/log(10) +1) print $0}' \
        $exp.RAMPAGE-RNA.Summary.txt > tmp.filter
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.filter $peakDir/$exp-RPM2-70.bed \
        > $exp-RPM2-70.RNA-Filtered.bed
    
done
