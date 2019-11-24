# RAMPAGE-Analysis

## Fetch unique reads and remove PCR duplicates
```
rm_pcr.py --uniq --output rampage_peak_folder --thread 10 RAMPAGE.bam
```

## Call peaks with F-seq
```
call_peak.py rampage_peak_folder
```
