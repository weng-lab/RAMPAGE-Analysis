
bio=$1
workingDir=~/Lab/ENCODE/RAMPAGE/ORF
scriptDir=~/GitHub/RAMPAGE-Analysis/rPeak-Analysis/Scripts

annotations=Filtered-All-Summaries.txt
#group="Verified GENCODE TSS"
#group="TSS unannotated transcript"
#group="Verified unannotated TSS"
#group="Candidate GENCODE TSS"
#group="Local transcription"
rm -f ORF-Summary.$bio.txt

groups=( "Verified_GENCODE_TSS" "TSS_unannotated_transcript" "Verified_unannotated_TSS" "Candidate_GENCODE_TSS" "Local_transcription") 
for groupKey in ${groups[@]}
do
echo $groupKey
group=$(echo $groupKey | awk '{gsub(/_/," ");print}')
sam=$bio-Merge.sam

cd $workingDir

awk -F "\t" -v var="$group" '{if ($5 == var) print $0}' $annotations | \
    awk 'FNR==NR {x[$1];next} ($4 in x)' - ../hg38-rPeaks.bed > tmp.bed
awk '{if ($6 == "+") print $0}' tmp.bed | \
    bedtools intersect -wo -a stdin -b ../PacBio/$bio.5ends.plus.bed > tmp.plus-intersect
awk '{if ($6 == "-") print $0}' tmp.bed | \
    bedtools intersect -wo -a stdin -b ../PacBio/$bio.5ends.minus.bed > tmp.minus-intersect

rm -f $groupKey.ORF-Summary.txt
modes=(plus minus)
for mode in ${modes[@]}
do
    awk '{print $4}' tmp.$mode-intersect | sort -u > tmp.$mode-list
    l=$(wc -l tmp.$mode-list | awk '{print $1}')
    echo $mode $l
    jobid=$(sbatch --nodes 1 --array=1-$l --mem=10G --time=04:00:00 \
        --output=/home/moorej3/Job-Logs/jobid_%A_%a.output \
        --error=/home/moorej3/Job-Logs/jobid_%A_%a.error \
        $scriptDir/Batch-Scan-ORFs.sh $mode tmp.$mode-list tmp.$mode-intersect \
        $sam $bio | awk '{print $NF}')
        
    sleep 20
    list=100
    while [ $list -gt 1 ]
    do
        list=$(squeue -j $jobid | wc -l | awk '{print $1}')
        echo -e "jobs still running: $list"
        sleep 10
    done
    
    cat ~/Job-Logs/jobid_$jobid_*.output | awk '{print $0 "\t" "'$groupKey'"}' >> ORF-Summary.$bio.txt
    rm ~/Job-Logs/jobid_$jobid_*
done
done
