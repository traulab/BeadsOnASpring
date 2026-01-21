#!/usr/bin/env python3
"""
@author Andrew D Johnston

Given one or more TF binding-site BED files, this script finds the +/- Nth nucleosome
peak relative to each TF site and outputs per-offset TSV files.

Key ideas
- Nucleosome peaks are loaded once (BED or BED.gz), grouped by chromosome and sorted by
  peak center for fast binary search (bisect).
- For each TF site, we compute its reference position ("TF center") and choose an orientation ("strand"):
    * TF center defaults to midpoint between BED start/end, but can be overridden with --tf_pos_col.
    * Strand uses a column from the TF BED if present (configurable; BED6 default col 6).
      Otherwise, infer an orientation by looking at the nearest nucleosome center.
- "+N" and "-N" are defined relative to the TF strand:
    * On '+' strand: plus = downstream (increasing genomic coordinate),
                   minus = upstream (decreasing coordinate).
    * On '-' strand: plus = upstream (decreasing coordinate),
                   minus = downstream (increasing coordinate).

Outputs
- For each TF BED and each offset k in [1..N], write:
    <tf_prefix>_plus<k>.tsv
    <tf_prefix>_minus<k>.tsv

Each output line:
    chrom  nuc_start  nuc_end  strand  nuc_score

Notes / assumptions
- This is a purely positional lookup based on nucleosome peak centers; it does not enforce
  that the nucleosome is within any particular window from the TF site.
- Nucleosome peak BED is expected to have at least 3 columns; score column is configurable.
"""

import argparse
from bisect import bisect_left
import gzip
import glob
import os
import sys
from typing import Optional

from tqdm import tqdm


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Find +N and -N nucleosome peaks for each TF binding site."
    )
    p.add_argument(
        "tf_beds",
        nargs="+",
        help="TF BED file(s) or glob pattern(s), e.g. *_hg19.bed or data/*.bed.gz",
    )
    p.add_argument(
        "-n",
        "--nuc_bed",
        required=True,
        help="Nucleosome peak BED file (.bed or .bed.gz).",
    )
    p.add_argument(
        "--num",
        type=int,
        default=100,
        help="Number of nucleosome peaks to find in each direction (default: 100).",
    )
    p.add_argument(
        "--tf_strand_col",
        type=int,
        default=6,
        help="1-based column in TF BED that contains strand (+/-). BED6 strand is column 6 (default: 6). "
             "If missing/invalid, strand is inferred.",
    )
    p.add_argument(
        "--tf_pos_col",
        type=int,
        default=None,
        help="Optional 1-based column containing a single-base TF position (e.g. motif center). "
             "If set, TF center is taken from this column. Default: use midpoint of BED start/end.",
    )
    p.add_argument(
        "--nuc_score_col",
        type=int,
        default=4,
        help="1-based column in nucleosome BED that contains peak score (default: 4). If missing, score is 'NA'.",
    )
    p.add_argument(
        "--out_dir",
        default=".",
        help="Directory to write output TSVs (default: current directory).",
    )
    return p.parse_args()


# ----------------------------
# I/O helpers
# ----------------------------
def smart_open(path: str):
    """Open a text file; use gzip if filename ends with .gz."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def expand_globs(patterns: list[str]) -> list[str]:
    """
    Expand glob patterns.

    - If a pattern matches files, include those.
    - If it doesn't match but is an existing path, include it as-is.
    - Otherwise warn.
    """
    out: list[str] = []
    for pat in patterns:
        matches = glob.glob(pat)
        if matches:
            out.extend(matches)
        else:
            if os.path.exists(pat):
                out.append(pat)
            else:
                print(f"[WARN] No files matched pattern: {pat}", file=sys.stderr)
    return out


# ----------------------------
# Nucleosome loading
# ----------------------------
def load_nucleosomes(nuc_bed_file: str, score_col_1based: int) -> dict:
    """
    Load nucleosome peaks into a dict keyed by chromosome.

    Returns:
      nuc_dict[chrom] = {
          'centers': [center0, center1, ...]  (sorted)
          'records': [(center, start, end, score), ...] (same order)
      }
    """
    score_idx = score_col_1based - 1
    tmp: dict[str, list[tuple[int, int, int, str]]] = {}

    with smart_open(nuc_bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip().split()
            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            score = parts[score_idx] if len(parts) > score_idx else "NA"
            center = (start + end) // 2

            tmp.setdefault(chrom, []).append((center, start, end, score))

    # Sort by center and build centers list (for bisect) once per chromosome
    out: dict[str, dict[str, list]] = {}
    for chrom, recs in tmp.items():
        recs.sort(key=lambda x: x[0])
        out[chrom] = {"centers": [r[0] for r in recs], "records": recs}

    return out


# ----------------------------
# Core lookup logic
# ----------------------------
def infer_strand(nuc_data: dict, tf_center: int) -> str:
    """
    Infer a strand for a TF site based on the nearest nucleosome center.

    Rule: return '+' if the nearest nucleosome center is left of tf_center, else '-'.
    """
    centers = nuc_data["centers"]
    idx = bisect_left(centers, tf_center)

    if idx == 0:
        nearest = centers[0]
    elif idx == len(centers):
        nearest = centers[-1]
    else:
        before = centers[idx - 1]
        after = centers[idx]
        nearest = before if abs(before - tf_center) <= abs(after - tf_center) else after

    return "+" if nearest < tf_center else "-"


def find_n_peak(nuc_data: dict, tf_center: int, strand: str, offset: int, direction: str):
    """
    Find the +N or -N nucleosome peak relative to tf_center and strand.

    direction: 'plus' or 'minus' (relative to TF strand)
    offset: 1-based index (1 = closest in that direction by center order)
    """
    centers = nuc_data["centers"]
    recs = nuc_data["records"]

    idx = bisect_left(centers, tf_center)

    if strand == "+":
        target_idx = idx + (offset - 1) if direction == "plus" else idx - offset
    elif strand == "-":
        target_idx = idx - offset if direction == "plus" else idx + (offset - 1)
    else:
        return None

    if 0 <= target_idx < len(recs):
        return recs[target_idx]  # (center, start, end, score)
    return None


def tf_strand_from_parts(parts: list[str], strand_col_1based: int) -> Optional[str]:
    """Return '+'/'-' if present and valid in the configured column, else None."""
    idx = strand_col_1based - 1
    if len(parts) > idx and parts[idx] in ("+", "-"):
        return parts[idx]
    return None


def tf_center_from_parts(parts: list[str], pos_col_1based: Optional[int]) -> Optional[int]:
    """
    Return TF reference position.

    - If pos_col_1based is provided, read TF position from that column (must be int).
    - Otherwise return None so caller can use BED midpoint.
    """
    if pos_col_1based is None:
        return None
    idx = pos_col_1based - 1
    if len(parts) <= idx:
        return None
    try:
        return int(parts[idx])
    except ValueError:
        return None


# ----------------------------
# Processing
# ----------------------------
def process_tf_bed(
    tf_bed: str,
    nuc_dict: dict,
    num: int,
    tf_strand_col_1based: int,
    tf_pos_col_1based: Optional[int],
    out_dir: str,
):
    """
    Process one TF BED, writing per-offset TSVs to out_dir.

    For each TF site, we determine:
      - TF center: from --tf_pos_col if set, else midpoint(start,end)
      - strand: from --tf_strand_col if valid, else inferred from nearest nucleosome
    Then we emit the +/-k nucleosome peaks for k=1..num.
    """
    # Pre-allocate per-offset buffers
    buffers = {f"plus{k}": [] for k in range(1, num + 1)}
    buffers.update({f"minus{k}": [] for k in range(1, num + 1)})

    with smart_open(tf_bed) as f:
        for line in tqdm(
            f,
            desc=os.path.basename(tf_bed),
            unit="line",
            dynamic_ncols=True,
        ):
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip().split()
            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            # TF reference position
            tf_center = tf_center_from_parts(parts, tf_pos_col_1based)
            if tf_center is None:
                tf_center = (start + end) // 2

            nuc_data = nuc_dict.get(chrom)
            if nuc_data is None or not nuc_data["records"]:
                continue

            # TF orientation
            strand = tf_strand_from_parts(parts, tf_strand_col_1based)
            if strand is None:
                strand = infer_strand(nuc_data, tf_center)

            for k in range(1, num + 1):
                peak_plus = find_n_peak(nuc_data, tf_center, strand, k, "plus")
                if peak_plus:
                    _, s, e, score = peak_plus
                    buffers[f"plus{k}"].append(f"{chrom}\t{s}\t{e}\t{strand}\t{score}")

                peak_minus = find_n_peak(nuc_data, tf_center, strand, k, "minus")
                if peak_minus:
                    _, s, e, score = peak_minus
                    buffers[f"minus{k}"].append(f"{chrom}\t{s}\t{e}\t{strand}\t{score}")

    # Derive prefix from TF BED filename
    base = os.path.basename(tf_bed)
    prefix = base.replace(".bed.gz", "").replace(".bed", "").replace(".gz", "")
    os.makedirs(out_dir, exist_ok=True)

    for key, lines in buffers.items():
        out_path = os.path.join(out_dir, f"{prefix}_{key}.tsv")
        with open(out_path, "w") as out_f:
            if lines:
                out_f.write("\n".join(lines) + "\n")


def main():
    args = parse_args()

    if args.num <= 0:
        print("[ERROR] --num must be > 0", file=sys.stderr)
        sys.exit(1)

    tf_beds = expand_globs(args.tf_beds)
    if not tf_beds:
        print("[ERROR] No TF BED files found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loading nucleosomes: {args.nuc_bed}")
    nuc_dict = load_nucleosomes(args.nuc_bed, score_col_1based=args.nuc_score_col)

    for tf_bed in tf_beds:
        print(f"[INFO] Processing TF BED: {tf_bed}")
        process_tf_bed(
            tf_bed=tf_bed,
            nuc_dict=nuc_dict,
            num=args.num,
            tf_strand_col_1based=args.tf_strand_col,
            tf_pos_col_1based=args.tf_pos_col,
            out_dir=args.out_dir,
        )


if __name__ == "__main__":
    main()
