# Weng Lab RAMPAGE rPeak Pipeline

---

## Step 0 - Calling RAMPAGE peaks


## Step 1 - Filtering peaks

```
./1_Peak-Filtering.sh
```

Requires:
* `Determine-Genomic-Context.sh`
* `RAMPAGE-RNA-Match-List.txt`

## Step 2 - Annotating rPeaks

```
./2_rPeak-Annotation.sh
```

Requires:
* `pick-best-peak.py`
* `filter-rPeaks.py`
* `accession-rPeaks.py`
