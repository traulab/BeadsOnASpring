# fragment_size_counter.py

This script assigns **paired-end BAM fragments** to **chromatin-state BED intervals** and computes **fragment length distributions per state**, optionally producing overlay **density plots**.

## What it does

- Reads a coordinate-sorted paired-end **BAM** (with index).
- Reads a **BED** (or BED-like) file containing chromatin-state intervals.
- For each properly paired fragment, computes:
  - `fragment_start = min(start/end positions of both mates)`
  - `fragment_end   = max(start/end positions of both mates)`
  - `fragment_mid   = (fragment_start + fragment_end) // 2`
  - `fragment_length = fragment_end - fragment_start`
- Assigns the fragment to the **BED interval that contains `fragment_mid`** (midpoint-based overlap).
- Aggregates counts of `fragment_length` per chromatin state.
- Writes one output file per state as `length<TAB>count`.
- Optionally generates overlaid density plots over a chosen length window.

## Inputs

### 1) BAM
- Must be **paired-end** and **coordinate sorted** (recommended).
- Must have an index (`.bai` or `.csi`) so `pysam.fetch()` works.
- Contig naming can be `chr1` or `1` style — the script detects BAM style and harmonizes BED contig names.

### 2) BED (chromatin states)
Expected fields (configurable via CLI args):
- `chrom` (col 1)
- `start` (col 2, 0-based)
- `end` (col 3, half-open)
- `state label` (default **col 4**, set with `--bed_state_col`)
- Optional `RGB` (default **col 9**, set with `--bed_rgb_col`) formatted `R,G,B` in 0–255.

The script assumes per-chromosome intervals are **sorted** and effectively **non-overlapping** for fast linear scanning.

## Outputs

All outputs go to `--output_dir` (default: `./<bam_basename_without_ext>/`).

### Per-state TSVs
One file per state:

```
<bam>_<state>_<label>.txt
```

where:
- `<bam>` is the BAM basename without `.bam`
- `<state>` is the state label sanitized for filenames
- `<label>` is either:
  - `autosomes1-22` (default mode), or
  - the chromosome name (if `--chromosome` is used)

Each TSV is:

```
fragment_length_bp<TAB>count
```

### Optional plots (`--plot`)
Creates:
1) **Main overlay plot** across all states:
   - Default name: `<bam>_<label>_dens_<xmin>-<xmax>.png`
2) **Focus overlay plot** including `euchromatin` and states matching `--focus_prefixes`:
   - Default name: `<bam>_<label>_dens_<xmin>-<xmax>_focus.png`

Density is **normalized within the plotting window** `[plot_xmin, plot_xmax]` *per state*:
`dens(L) = count(L) / sum_{plot_xmin..plot_xmax} count(L)`.

## Synthetic categories

The script can create pooled synthetic categories:

### `euchromatin`
- Built by pooling states whose names start with any prefix in `--euchromatin_prefixes`
- Default prefixes match your previous ChromHMM conventions.

### `whole-genome`
- Built by summing counts across all states (excluding synthetic categories to avoid double counting).

## Installation

### Python dependencies
- `pysam`
- `tqdm`
- `matplotlib` (only needed if using `--plot`)

Example (conda/mamba):

```bash
mamba install -c bioconda pysam
mamba install -c conda-forge tqdm matplotlib
```

## Usage

### Help
```bash
python fragment_size_counter_ChromHMM.py -h
```

### Default: pool autosomes 1–22
```bash
python fragment_size_counter_ChromHMM.py \
  -b sample.bam \
  -d wgEncodeBroadHmmGm12878HMM.bed
```

### Single chromosome
```bash
python fragment_size_counter_ChromHMM.py \
  -b sample.bam \
  -d wgEncodeBroadHmmGm12878HMM.bed \
  --chromosome chr19
```

### Create plots over a custom range
```bash
python fragment_size_counter_ChromHMM.py \
  -b sample.bam \
  -d wgEncodeBroadHmmGm12878HMM.bed \
  --plot \
  --plot_xmin 0 \
  --plot_xmax 500
```

### If your BED has state and/or RGB in different columns
```bash
python fragment_size_counter_ChromHMM.py \
  -b sample.bam \
  -d states.bed \
  --bed_state_col 4 \
  --bed_rgb_col 9
```

### Customize synthetic categories / overrides
```bash
python fragment_size_counter_ChromHMM.py \
  -b sample.bam \
  -d wgEncodeBroadHmmGm12878HMM.bed \
  --euchromatin_prefixes "1_,2_,4_,5_,6_,7_,9_,10_,11_" \
  --grey_override_prefixes "13_" \
  --euchromatin_color_rgb "128,0,128" \
  --grey_override_color_rgb "128,128,128" \
  --whole_genome_color_rgb "0,0,0"
```

## Important notes / limitations

### Chunk-local pairing
Reads are paired **within each chunk only**. If mates fall into different chunks, they may not be paired/processed.

- This is usually fine for **coordinate-sorted BAMs** with typical insert sizes, especially with a large `--chunk_size`.
- If you have very large inserts or unusual pairing patterns, consider increasing `--chunk_size`.

### Interval assumptions
Assignment uses a linear scan (`state_idx`) across sorted intervals. For correctness, per-chromosome BED intervals should be **sorted** and ideally **non-overlapping**.

### Assignment coordinate
Assignment is based on the **fragment midpoint**.
Choose this versus start/end assignment based on what best matches your biological question.

## License
Add your preferred license here.
