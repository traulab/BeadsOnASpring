# chromatin_state_shuffler.py

This script shuffles ChromHMM (or any **BED4**) segments **within each chromosome** while preserving:

- the **multiset of segment lengths** per chromosome  
- the **multiset of state/compartment labels** per chromosome  
- the **(length, label) pairing** (pairs are shuffled as intact tuples)

It writes a new BED file where shuffled segments are placed at candidate “anchor” positions sampled from the original annotated regions.

---

## What it does (high level)

For each chromosome:

1. **Read segments** from the input BED and collect `(segment_length_bp, state_label)` tuples into `D_compartments[chrom]`.
2. Build a list of candidate **anchor starts** every `--anchor-step` bp inside all **non-skip-label** segments (default skip label: `NA`) into `D_regions[chrom]`.

   Example: if a segment starts at 1000 and `--anchor-step=200`, anchors are `1000, 1200, 1400, ...`.

3. **Shuffle** the `(length, label)` list for that chromosome.
4. Walk anchor starts in increasing order and place each shuffled segment at the first anchor that **does not overlap** the previously placed segment.
5. **Write** placed segments to an output BED4 file.

---

## Input

### Required format: BED4 (tab-delimited)

Each line must have at least 4 columns:

1. `chrom`
2. `start` (0-based)
3. `end` (exclusive)
4. `state_label` (e.g., ChromHMM state)

Example:

```tsv
chr1	0	10000	1_TssA
chr1	10000	20000	13_Het
```

### Special handling

- Lines where `state_label == --skip-label` are **ignored** (default: `NA`)
- Segments with `end <= start` are skipped defensively
- Comment lines beginning with `#` and blank lines are ignored

---

## Output

A BED4 file specified by `--output-bed` (default: `wgEncodeBroadHmmGm12878HMM_shuffled.bed`).

Each output line is:

```tsv
chrom	shuffled_start	shuffled_end	state_label
```

Where:

- `shuffled_end = shuffled_start + segment_length_bp`
- `state_label` stays paired with the original segment length it came with

---

## Configuration

The script is configured via **command-line arguments** (defaults shown):

- `-i, --input-bed`  
  Default: `wgEncodeBroadHmmGm12878HMM.bed`

- `-o, --output-bed`  
  Default: `wgEncodeBroadHmmGm12878HMM_shuffled.bed`

- `--anchor-step`  
  Default: `200` (bp)  
  Controls the spacing between candidate anchor starts within annotated regions.

- `--skip-label`  
  Default: `NA`  
  Segment label treated as unannotated and skipped.

- `--seed`  
  Default: unset  
  Random seed for reproducible shuffling.

---

## Usage

### Minimal usage (defaults)

```bash
python shuffle_chromhmm.py
```

### Typical usage

```bash
python shuffle_chromhmm.py \
  --input-bed wgEncodeBroadHmmGm12878HMM.bed \
  --output-bed wgEncodeBroadHmmGm12878HMM_shuffled.bed \
  --anchor-step 200 \
  --seed 1
```

---

## Reproducibility (randomness)

If you want the same output across runs, pass a seed:

```bash
--seed 1
```

(Any integer seed is fine.)

---

## Important notes / limitations

### 1) Not all segments may be placed

Placement depends on anchor availability and non-overlap constraints. If anchors are sparse
(e.g., large `--anchor-step`) or segments are long, the script may run out of valid anchors
before all shuffled segments are placed.

**Tip:** Reduce `--anchor-step` (e.g., 50 or 100) to increase placement opportunities.

### 2) Placement restricted to originally annotated regions

Anchors are generated only from regions whose label is **not** equal to `--skip-label` (default: `NA`).
Shuffled segments can only begin at these anchors.

### 3) Output segments can extend beyond the original annotated interval that produced an anchor

Anchors are drawn from within an original segment, but the **shuffled segment length** may be longer
than the remaining space to that segment’s end. The script does **not** currently check whether
`seg_end` remains within:

- the original annotated segment,
- any annotated region at all, or
- chromosome bounds.

If this matters downstream, consider adding a constraint to ensure `seg_end` stays within
an allowed region and/or a chrom-size file.

### 4) Chromosome order

Chromosomes are processed in the internal dictionary order. If you want a specific ordering
(e.g., chr1..chr22, chrX, chrY), sort chromosome keys before iterating.

---

## Suggested improvements (optional)

If you plan to reuse this often, consider adding:

- A summary report with per-chromosome counts:
  - segments read
  - anchors generated
  - segments placed
  - segments dropped
- Bounds checks against chromosome sizes (`chrom.sizes`)
- Optional “place-within-original-interval” constraint

---

