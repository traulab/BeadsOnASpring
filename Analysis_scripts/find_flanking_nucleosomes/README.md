# find_flanking_nucleosomes.py

This script finds the **+N** and **−N** nucleosome peaks relative to each TF binding site and writes one output TSV per offset.

## What it does

- Loads nucleosome peaks from a BED (`.bed` or `.bed.gz`), grouped by chromosome and sorted by **peak center**.
- For each TF binding site in one or more TF BED files (`.bed` or `.bed.gz`):
  - Determines the TF **reference position** (“TF center”):
    - By default uses the **midpoint** between BED start/end (`(start+end)//2`)
    - If `--tf_pos_col` is set, uses the integer position from that column instead
  - Determines TF **strand**:
    - Uses strand from TF BED column `--tf_strand_col` (default: **6**, BED6 standard)
    - If missing/invalid, strand is **inferred** from the nearest nucleosome peak center
  - Retrieves:
    - `plusN`: the Nth nucleosome peak in the **plus direction** relative to the TF strand
    - `minusN`: the Nth nucleosome peak in the **minus direction** relative to the TF strand

### Direction definition

For a TF on:

- **'+' strand**
  - `plusN` = downstream (higher genomic coordinate)
  - `minusN` = upstream (lower genomic coordinate)

- **'−' strand**
  - `plusN` = upstream (lower genomic coordinate)
  - `minusN` = downstream (higher genomic coordinate)

## Inputs

### TF BED(s)
Positional BED with at least 3 columns:

```
chrom  start  end
```

Optional strand column:
- If present at `--tf_strand_col` (default: **6**, 1-based) and contains `+` or `-`, it is used.
- Otherwise the script infers strand from the nearest nucleosome center.

Optional TF position column:
- If `--tf_pos_col` is set, the TF reference position is taken from that 1-based column.
- Otherwise TF position is the midpoint of BED start/end.

You can pass one or more files and/or glob patterns:
- `TF1.bed TF2.bed.gz data/*.bed`

### Nucleosome BED
BED with at least 3 columns:

```
chrom  start  end
```

Optional score column:
- Defaults to **column 4** (1-based), set with `--nuc_score_col`.
- If missing, score is recorded as `NA`.

Supports `.bed` or `.bed.gz`.

## Outputs

For each TF BED and each offset `k` in `1..N`, the script writes:

- `<tf_prefix>_plus<k>.tsv`
- `<tf_prefix>_minus<k>.tsv`

By default outputs are written to the current directory; change with `--out_dir`.

Each TSV line is:

```
chrom    nuc_start    nuc_end    strand    nuc_score
```

## Usage

Show help:

```bash
python tf_peak_to_nucleosome_peaks.py -h
```

Typical run (find ±100 peaks):

```bash
python tf_peak_to_nucleosome_peaks.py \
  data/CTCF_sites.bed.gz data/RAD21_sites.bed \
  -n peaks/nucleosome_peaks.bed.gz \
  --num 100 \
  --out_dir out
```

### BED6 strand (standard) — no need to set anything
If your TF BED is BED6 with strand in column 6 (standard), the defaults already work:

```bash
python tf_peak_to_nucleosome_peaks.py TF_sites.bed -n nucleosome_peaks.bed --num 50
```

### Using an explicit TF position column
If you have a single-base motif position in column 7 (example), use:

```bash
python tf_peak_to_nucleosome_peaks.py TF_sites.bed \
  -n nucleosome_peaks.bed \
  --num 50 \
  --tf_pos_col 7
```

### Strand in a different column
If your strand is in column 4 instead of 6:

```bash
python tf_peak_to_nucleosome_peaks.py TF_sites.bed \
  -n nucleosome_peaks.bed \
  --num 50 \
  --tf_strand_col 4
```

### Nucleosome score not in column 4
If nucleosome score is in column 5:

```bash
python tf_peak_to_nucleosome_peaks.py TF_sites.bed \
  -n nucleosome_peaks.bed \
  --nuc_score_col 5
```

