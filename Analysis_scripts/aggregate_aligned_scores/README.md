# BigWig window alignment + aggregation (meta-profile)

This mini-pipeline has two steps:

1. **Extract strand-aware per-base score windows** from a BigWig around each feature in an alignment BED/BED-like file  
   → produces a `*.txt` scores matrix (1 line per feature).
2. **Aggregate** those aligned windows into a single **mean meta-profile**  
   → produces a `*_agg.tsv` file (relative position vs mean score).

---

## 1) Extract aligned windows from a BigWig (shell script)

### Script
**`extract_bigwig_windows.sh`**

### Purpose
For each row in an alignment BED/BED-like file:

- choose an **alignment point** (default = midpoint between `start` and `end`)
- extract a fixed window `[center-window_half, center+window_half)` from the BigWig
- expand it to **per-base values**, filling missing bases with `nan`
- if `strand == "-"`, **reverse** the vector so all rows are oriented consistently
- write one space-separated line of scores per feature

### Requirements
- `bigWigToBedGraph` in `$PATH` (UCSC tools)
- `awk`
- bash

### Usage
```bash
chmod +x extract_bigwig_windows.sh

./extract_bigwig_windows.sh \
  -w BH01_chrAll_dup1_WPS_window120_lower120_upper180_WPS.bw \
  -a BH01_chrAll_dup1_WPS_window120_lower120_upper180_nucleosome_regions_20k_rand_peaks.bed
```

### Output
- Default output filename: `<alignment_file_basename>_BW.txt`
- Format: one line per feature, space-separated numeric values (`nan` where missing)

### Common options
#### Use a custom alignment point column
If your BED-like file already contains an explicit alignment coordinate (e.g. column 7),
you can align to that column instead of midpoint:

```bash
./extract_bigwig_windows.sh -w track.bw -a peaks.bed --point-col 7
```

#### Change window size
Window length is `2 * window_half + 1` bases.

Default: `--window-half 3000` → output length = `6001`

```bash
./extract_bigwig_windows.sh -w track.bw -a peaks.bed --window-half 2000
```

#### Column control
Defaults are BED-like:

- chrom = column 1  
- start = column 2  
- end   = column 3  
- strand = column 6  

If your file differs:

```bash
./extract_bigwig_windows.sh -w track.bw -a peaks.bed --strand-col 4
```

---

## 2) Aggregate aligned windows (Python script)

### Script
**`aggregate_aligned_scores.py`**

### Purpose
Reads the score matrix produced above and computes the **mean score at each aligned position**.

Optional QC filters skip bad rows:

- skip NaN/Inf lines
- skip lines with ≥ N consecutive zeros
- skip lines containing values above a max threshold
- skip lines of inconsistent length

### Requirements
- Python 3
- `numpy`
- `tqdm`

Install example (mamba/conda):
```bash
mamba install numpy tqdm
```

### Usage
```bash
python3 aggregate_aligned_scores.py \
  --scores BH01_chrAll_dup1_WPS_window120_lower120_upper180_nucleosome_regions_20k_rand_peaks_BW.txt
```

### Output
Default: `<scores_basename>_agg.tsv`

TSV columns:
- `relative_position` (0 = center)
- `score` (mean across all kept regions)

### Common options
```bash
python3 aggregate_aligned_scores.py \
  --scores scores.txt \
  --zero_thresh 0 \
  --max_score 300 \
  --output meta_profile.tsv
```

---

## Full example pipeline

```bash
# 1) Extract strand-oriented windows around each feature
./extract_bigwig_windows.sh \
  -w BH01_chrAll_dup1_WPS_window120_lower120_upper180_WPS.bw \
  -a BH01_chrAll_dup1_WPS_window120_lower120_upper180_nucleosome_regions_20k_rand_peaks.bed \
  -o BH01_WPS_aligned_windows.txt

# 2) Aggregate into mean meta-profile
python3 aggregate_aligned_scores.py \
  --scores BH01_WPS_aligned_windows.txt \
  --output BH01_WPS_meta_profile.tsv
```

---

## Notes

- Coordinates extracted from BigWig use **0-based, half-open intervals**: `[start, end)`.
- When strand is missing / not `+` or `-`, it is randomly assigned so that orientations remain balanced.

