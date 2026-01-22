# WPS with Nucleosome + Breakpoint Peak Calling (Kircher-style)

Tools for generating **WPS (Windowed Protection Score)** tracks from paired-end BAM files, applying **rolling-median baseline correction**, optional **Savitzky–Golay smoothing**, and calling:

- **Positive peaks** → **nucleosome regions**
- **Negative peaks (troughs)** → **breakpoint peaks** (called by flipping the sign and re-using the same peak caller)

This implementation is designed to match **Martin Kircher’s 2015 WPS behavior** as closely as possible while operating on standard BAM inputs (0-based, half-open coordinates via `pysam`).

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Concept: WPS scoring](#concept-wps-scoring)
  - [Kircher-style overlap semantics](#kircher-style-overlap-semantics)
  - [Fragment size filtering](#fragment-size-filtering)
- [Baseline correction and smoothing](#baseline-correction-and-smoothing)
- [Peak calling (Kircher-style)](#peak-calling-kircher-style)
- [Script usage](#script-usage)
  - [Basic examples](#basic-examples)
  - [Arguments](#arguments)
  - [How contigs / regions are processed](#how-contigs--regions-are-processed)
- [Outputs](#outputs)
  - [1) Combined per-base scores bedGraph](#1-combined-per-base-scores-bedgraph)
  - [2) Nucleosome regions (positive peaks)](#2-nucleosome-regions-positive-peaks)
  - [3) Breakpoint peaks (negative peaks)](#3-breakpoint-peaks-negative-peaks)
- [Converting to bigWig / bigBed](#converting-to-bigwig--bigbed)
- [Conda / Mamba environment setup](#conda--mamba-environment-setup)
- [Notes / gotchas](#notes--gotchas)

---

## Overview

Pipeline (high level):

1. Reads **paired-end fragments** from one or more BAM files across contigs or user-specified regions.
2. Filters read-pairs:
   - unmapped/mate-unmapped removed
   - same-strand pairs removed
   - duplicate/qcfail removed
   - soft/hard-clipped or padded CIGAR removed
3. Collapses duplicate fragments by coordinate tuple `(chrom, frag_start, frag_end)` to `--max-duplicates`.
4. Optionally subsamples fragments (`--subsample`).
5. For each fragment whose length is within `[--frag-lower, --frag-upper]`:
   - Adds a **Kircher-equivalent WPS kernel** across the fragment footprint.
   - Builds a **coverage** track (+1 for each covered base).
   - Builds a **dyad** track (+1 at fragment center).
6. Computes **rolling-median baseline** (default 1000 bp) and subtracts it.
7. Optionally smooths WPS / baselined WPS (Savitzky–Golay; default 21/2).
8. Calls peaks/regions on the **baselined, smoothed** track:
   - **positive runs** → nucleosome regions
   - **negative troughs** (via sign flip) → breakpoint peaks
9. Writes:
   - `*_combined_scores.bedGraph` (multi-column per-base track)
   - `*_nucleosome_regions.bed` (Kircher-style BED8 calls)
   - `*_breakpoint_peaks.bed` (Kircher-style BED8 calls)

---

## Requirements

### Python
- Python 3

### Required Python libraries
- `numpy`
- `scipy`
- `pysam`
- `tqdm`

### Additional tools (optional but recommended)
- `samtools` (inspect BAMs, indexing)
- UCSC tools:
  - `bedGraphToBigWig`
  - `bedToBigBed` *(for peak call tracks)*

---

## Concept: WPS scoring

Windowed Protection Score (WPS) is a fragmentomics signal intended to quantify **nucleosome protection**:

- For each genome position, define a symmetric protection window (e.g. 120 bp).
- Fragments that fully span the window contribute positively.
- Fragments that truncate within the window contribute negatively.

This script implements WPS **via a per-fragment kernel sum**, which is computationally efficient and reproduces Kircher’s effective overlap semantics.

### Kircher-style overlap semantics

Kircher’s original implementation uses bx-python interval logic, which behaves like **half-open intervals**.
When combined with his 1-based inclusive read coordinates, the effective behavior is:

- protection window length behaves like `2*half` (e.g. 120 bp), not 121 bp
- effective fragment is 1 bp shorter on the right for span/truncate classification

This script reproduces that behavior by precomputing a **Kircher-exact kernel** for each fragment length.

### Fragment size filtering

Only fragments in the configured length range contribute to WPS:

```
--frag-lower 120
--frag-upper 180
```

Coverage and dyad tracks are computed using true fragment coordinates.

---

## Baseline correction and smoothing

WPS is baseline-shifted using a rolling median:

- Default baseline window: **1000 bp**
- Baseline is computed on the chosen WPS track (raw or smoothed; configurable in code version)
- Subtracted to generate a median-centered track:
  - `wps_baselined = wps - rolling_median(wps, 1000)`

Smoothing uses Savitzky–Golay:

- Default `--sg-window 21`
- Default `--sg-order 2`

The peak caller uses:

- `wps_baselined_smoothed`

---

## Peak calling (Kircher-style)

Peak calling replicates Kircher’s `evaluateValues()` logic:

1. Build **positive runs** where `ivalue > 0`
2. Merge runs if separated by gaps ≤ `--peak-merge-gap` (default 5 bp), filling the gap with zeros
3. Evaluate merged runs:
   - Accept region lengths in `[--peak-minlen, --peak-maxlen]` (default 50–150 bp)
   - For long runs:
     - if region length in `[--peak-maxlen, 3*--peak-maxlen]` (150–450 bp), split by **median thresholding** into candidate windows
     - reject if region length > 450 bp (default)
4. **Score threshold:** report only if the **maximum value inside a called window** exceeds:
   - `--peak-varicutoff` (default 5)

This cutoff corresponds to Kircher’s `variCutoff`.

---

## Script usage

### Basic examples

#### Whole genome (all contigs in BAM header)
```bash
python3 wps_with_nucleosome_peak_calling.py \
  -b sample.bam
```

#### Restrict to a contig
```bash
python3 wps_with_nucleosome_peak_calling.py \
  -b sample.bam \
  -c chr12
```

#### Restrict to a genomic interval
Coordinates are interpreted as **0-based** with **end-exclusive** semantics:
```bash
python3 wps_with_nucleosome_peak_calling.py \
  -b sample.bam \
  -c chr12:52621135-52641135
```

#### Match Kircher-like defaults (typical)
```bash
python3 wps_with_nucleosome_peak_calling.py \
  -b SRR452466.bowtie.sorted.bam \
  --protection 120 \
  --frag-lower 120 --frag-upper 180 \
  -c chr12:52621135-52641135
```

#### Multiple BAMs (pooled fragments)
```bash
python3 wps_with_nucleosome_peak_calling.py \
  -b sample1.bam sample2.bam \
  -c chr12
```

---

## Arguments

### Core inputs

| Argument | Meaning |
|---|---|
| `-b/--bamfiles` | One or more paired-end BAM files |
| `-o/--out_prefix` | Output prefix (default derived from BAM basenames + region info) |
| `-c/--contigs` | One or more contigs, optionally ranges `chr:start-end` |

### WPS scoring parameters

| Argument | Meaning |
|---|---|
| `--protection` | Protection window size (bp), default `120` |
| `--frag-lower` | Minimum fragment length contributing to WPS |
| `--frag-upper` | Maximum fragment length contributing to WPS |
| `--max-duplicates` | Max allowed fragments with identical `(chrom,start,end)` |
| `--subsample` | Keep each fragment with probability p (e.g. 0.5) |

### Chunking / padding parameters

| Argument | Meaning |
|---|---|
| `--chunk-bp` | Chunk size per contig (default 100,000 bp) |
| `--overlap-bp` | Overlap padding (default 1,000 bp) |

### Baseline / smoothing parameters

| Argument | Meaning |
|---|---|
| `--baseline-window` | Rolling median window (default 1000 bp) |
| `--sg-window` | Savitzky–Golay smoothing window (odd, default 21) |
| `--sg-order` | Savitzky–Golay polynomial order (default 2) |

### Peak calling parameters (Kircher defaults)

| Argument | Meaning |
|---|---|
| `--peak-minlen` | Minimum length for candidate regions/windows (default 50 bp) |
| `--peak-maxlen` | Maximum length for reported windows (default 150 bp) |
| `--peak-maxregion` | Reject positive regions longer than this (default 450 bp) |
| `--peak-merge-gap` | Merge positive runs across gaps ≤ this (default 5 bp) |
| `--peak-varicutoff` | Minimum max score in window to report (default 5.0) |

---

## How contigs / regions are processed

To scale to genome-wide BAMs efficiently, the script processes contigs in **windows**:

- window length: **`--chunk-bp`** (default 100,000 bp)
- overlap padding: **`--overlap-bp`** (default 1,000 bp on each side)

For each chunk, scoring is computed on the **adjusted** window including overlaps. Output writing trims back to the **original non-overhang** interval to avoid duplicated output from overlaps.

---

## Outputs

All outputs are prefixed with:

```
<out_prefix>_prot<PROTECTION>_lower<LOWER>_upper<UPPER>_maxdup<MAXDUP>
```

### 1) Combined per-base scores bedGraph

**File:**
```
<out_prefix>_combined_scores.bedGraph
```

**Format:**
1. chrom
2. start
3. end
4. coverage
5. dyad
6. wps
7. wps_smoothed
8. wps_baselined
9. wps_baselined_smoothed

Example:
```
chr12   52621135  52621136  41  6  12.0  11.8  2.4  2.1
```

---

### 2) Nucleosome regions (positive peaks)

**File:**
```
<out_prefix>_nucleosome_regions.bed
```

This is a **BED8** file formatted like Kircher’s output:

Columns:
1. `chrom`
2. `start` (BED 0-based)
3. `end` (BED end-exclusive)
4. `name` (chr:start-end in 1-based display style)
5. `score` (max score in called window; reported as integer)
6. `strand` (`+`)
7. `thickStart` (peak midpoint start)
8. `thickEnd` (peak midpoint end)

Example:
```
chr12   51748712   51748745   12:51748713-51748745   6   +   51748728   51748729
```

---

### 3) Breakpoint peaks (negative peaks)

**File:**
```
<out_prefix>_breakpoint_peaks.bed
```

Same BED8 structure as nucleosome peaks, but called on the sign-flipped track to identify troughs.

---

## Converting to bigWig / bigBed

### Extract a single column from combined bedGraph

Example: baselined+smoothed WPS (last column):
```bash
in="<prefix>_combined_scores.bedGraph"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$9}' "$in" > wps_baselined_smoothed.bedGraph
```

### bedGraph → bigWig
```bash
bedGraphToBigWig wps_baselined_smoothed.bedGraph chrom.sizes wps_baselined_smoothed.bw
```

### Peaks → bigBed

Peaks are already BED8 and can be converted directly:
```bash
bedToBigBed <prefix>_nucleosome_regions.bed chrom.sizes <prefix>_nucleosome_regions.bb
bedToBigBed <prefix>_breakpoint_peaks.bed chrom.sizes <prefix>_breakpoint_peaks.bb
```

---

## Conda / Mamba environment setup

Example minimal environment:

```bash
mamba create -n wps_env -c conda-forge python=3.11 numpy scipy pysam tqdm
conda activate wps_env
```

---

## Notes

- Input BAMs must be indexed (`.bai` in same directory).
- The `-c contig:start-end` interval uses BED-style semantics:
  - start is 0-based
  - end is exclusive
- Duplicate filtering is coordinate-based:
  ```
  (chrom, fragment_start, fragment_end)
  ```
- Peak calling follows Kircher’s constraints:
  - report windows 50–150 bp
  - allow splitting of 150–450 bp regions by median-threshold windows
  - reject >450 bp regions
  - require max score > 5 to report

