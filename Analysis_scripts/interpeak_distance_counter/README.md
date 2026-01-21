# interpeak_distance_counter.py

This script computes **adjacent peak distance distributions** from a BED-like peak file, optionally stratified
by **chromatin state**, while supporting score thresholding by:

- a **single percentile**
- a **sweep of percentiles**
- a **target number of peaks**

It also supports configurable inter-peak distance bounds via:

- `--min-distance` (default: 50 bp)
- `--max-distance` (default: 1000 bp)

---

## What the script does

Given a set of genomic peaks:

1. Loads all peaks from an input BED-like file.
2. Optionally assigns each peak a chromatin state using an interval overlap query against a chromatin-state BED file.
3. Filters peaks based on a score threshold (percentile, sweep, or target peaks).
4. For each chromosome:
   - scans **adjacent kept peaks**
   - counts the distance **only if both peaks share the same chromatin state label**
   - counts distances only if `min_distance <= distance <= max_distance`
5. Writes a `*_nuc_dis.txt` report containing:
   - raw counts
   - smoothed counts
   - normalized (%) values
   - mode/median/mean summaries
   - duplicate position diagnostics

---

## Requirements

Python packages used:

- `numpy`
- `scipy`
- `intervaltree`

Example install:

```bash
mamba install -c conda-forge numpy scipy intervaltree
# or
pip install numpy scipy intervaltree
```

---

## Inputs

### 1) Peak file (required)

BED-like, whitespace-delimited.

Minimum required columns:

- `col0 = chrom`
- `col1 = start`
- `col2 = end`

Score filtering:

- default `--score-column 3` (4th column)
- must be numeric (float OK)

Peak position:

- default: midpoint `(start + end) // 2`
- or set `--position-column` to use a column containing the peak position directly

Example peak lines:

```
chr1  100  267  42.7
chr1  300  467  10.2
```

---

### 2) Chromatin state file (optional)

BED-like, whitespace-delimited.

Minimum required columns:

- `col0 = chrom`
- `col1 = start`
- `col2 = end`
- `col3 = state_label`

Example chromatin state lines:

```
chr1  0     10000  13_Het
chr1  10000 20000  1_Active_Promoter
```

If you omit this file, peaks are labeled as a single state: `All`.

---

## Chromatin state assignment

For each peak position `pos`, the script finds all chromatin-state intervals overlapping that position.

Rules:

- if **one** state overlaps → `state_label = that state`
- if **multiple** states overlap → `state_label = state1|state2|...` (sorted)
- if **none** overlap → `state_label = NA`

### Optional: combine euchromatin states

Using:

```bash
--combine-euchromatin
```

Any overlap that includes state labels beginning with:

`1_ 2_ 4_ 5_ 6_ 7_ 9_ 10_ 11_`

is collapsed to:

`Euchromatin`

---

## Filtering modes

Choose one of the following.

### A) Single percentile mode (default)

Threshold is chosen as:

- `threshold = np.percentile(scores, P)`

Peaks are kept if:

- `score >= threshold`

Command:

```bash
./script.py peaks.bed states.bed --score-percentile 25
```

Output filename convention:

- `<peakfile>_<statefile>_scorepct<P>_nuc_dis.txt`
- if no chromatin state file is used:
  - `<peakfile>_whole_genome_scorepct<P>_nuc_dis.txt`

> Note: `P` is the percentile **requested**, not the score threshold value.

---

### B) Percentile sweep mode

Generates many outputs over a percentile range.

Command (0 to 99 by 1):

```bash
./script.py peaks.bed states.bed \
  --pct-range --pct-lower 0 --pct-upper 99 --pct-step 1
```

Produces:

- `..._scorepct0_nuc_dis.txt`
- `..._scorepct1_nuc_dis.txt`
- ...
- `..._scorepct99_nuc_dis.txt`

---

### C) Target-peaks mode (keep ~N peaks)

Chooses an “effective percentile” so that approximately `N` peaks are retained.

Definitions:

- `eff_percentile = 100 * (1 - target_peaks / n_used)`
- `threshold = np.percentile(scores, eff_percentile)`

Command:

```bash
./script.py peaks.bed states.bed --target-peaks 7850000
```

Produces:

- `..._scorepct<effective_percentile>_nuc_dis.txt`

---

## Optional output: write filtered BED files

Add:

```bash
--write-filtered-bed
```

Behavior:

- Single percentile / target-peaks:
  - writes one filtered file
- Sweep mode:
  - writes one per percentile

Filtered BED filename format:

- `<input>_scorepct<P><ext>`

Example:

```bash
./script.py peaks.bed states.bed \
  --pct-range --pct-lower 0 --pct-upper 5 --pct-step 1 \
  --write-filtered-bed
```

---

## Distance bounds (NEW)

By default the script counts distances only in:

- 50–1000 bp

You can change this using:

- `--min-distance <int>`
- `--max-distance <int>`

Example (restrict to 120–250 bp):

```bash
./script.py peaks.bed states.bed \
  --score-percentile 10 \
  --min-distance 120 --max-distance 250
```

Validation rules:

- both must be `>= 0`
- `--max-distance` must be `>= --min-distance`

---

## Distance counting rules (critical)

Distances are computed per chromosome from peaks sorted by position, but only **kept peaks** (passing the score threshold) are considered.

A distance is counted only when:

1. Both peaks are kept (`score >= threshold`)
2. They are **adjacent among kept peaks**
3. They have **exactly the same** `state_label`
4. Distance satisfies: `min_distance <= distance <= max_distance`

Therefore, this measures **within-state adjacent spacing**, not generic adjacent spacing.

---

## Output file format (`*_nuc_dis.txt`)

The output file contains multiple sections.

### 1) Header lines

Includes score-filtering metadata:

- score column used
- requested percentile
- target peaks (if any)
- effective percentile
- chosen threshold
- usable scored lines
- total scanned lines

And distance-filtering metadata:

- `min_distance`
- `max_distance`

---

### 2) Chromosome distances by chromatin state

For each chromosome and each state:

- Mode distance (from smoothed counts)
- Median distance (weighted, raw)
- Mean distance (weighted, raw)

Then a tab-delimited table with columns:

1. chrom  
2. state  
3. distance (bp)  
4. count (raw)  
5. smoothed_count (Savitzky–Golay; rounded; non-negative)  
6. normalized_percent = `count / total_for_(chrom,state) * 100`  
7. smoothed_normalized_percent  

---

### 3) Chromosome distances (all states)

Same as above but pooling all states per chromosome.

---

### 4) Whole genome distances by chromatin state

Same as above but pooling all chromosomes per state.

---

### 5) Whole genome distances (all states)

A single pooled distribution.

---

### 6) Duplicate position report

Lists peak positions that occur more than once among kept peaks:

- `chrom   position   states...`

---

## Common usage recipes

### Whole-genome mode (no chromatin states), single run

```bash
./script.py peaks.bed --score-percentile 50
```

---

### ChromHMM stratified sweep 0–99

```bash
./script.py peaks.bed wgEncodeBroadHmmGm12878HMM.bed \
  --pct-range --pct-lower 0 --pct-upper 99 --pct-step 1
```

---

### Keep ~7.85 million peaks (e.g. WPS-comparable) + write filtered BED

```bash
./script.py peaks.bed wgEncodeBroadHmmGm12878HMM.bed \
  --target-peaks 7850000 --write-filtered-bed
```

---

### Setting score and position columns

Example: score in column 5, peak position in column 7

```bash
./script.py peaks.bed states.bed \
  --score-column 5 --position-column 7 --score-percentile 10
```

---
