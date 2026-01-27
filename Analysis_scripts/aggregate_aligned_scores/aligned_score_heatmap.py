#!/usr/bin/env python3
"""
@author Andrew D Johnston

Plot a heatmap of per-peak nucleosome score vectors (e.g. NPS/WPS-like tracks)
that are already aligned such that the centre position corresponds to x=0.

Input format:
- A plain text file where each line corresponds to one region/peak.
- Each line contains space-separated numeric values representing the score track
  across a fixed-width window centred on the alignment point.

What the script does:
1) Loads score vectors line-by-line and filters out bad rows:
   - non-numeric or non-finite values
   - rows containing >= --zero_thresh consecutive zeros
   - rows containing any value > --max_score
2) Optionally limits the number of rows:
   - --subsample-mode first: stop reading after --max_lines valid rows (fast)
   - --subsample-mode random: load all valid rows then randomly subsample to --max_lines
     (use --seed for reproducibility)
3) Optionally crops the x-axis window to the centre fraction using --breadth
   (always centred on x=0).
4) Optionally sorts rows prior to plotting:
   - center: descending score at x=0
   - rise_after_min: distance from the local minimum (left of centre) to the first
     subsequent zero-crossing, tie-broken by centre score
   - absmean: descending mean absolute score across the window
   - max: descending maximum score across the window
   - unsorted: keep input order
5) Produces two output figures:
   - Heatmap image showing sorted regions (rows) vs genomic position (columns),
     using a diverging colormap (seismic) with optional --vmin/--vmax.
   - Mean profile plot showing the average score at each position across all rows.

Outputs:
- <basename>_heatmap.png (unless --output specified)
- <basename>_heatmap_mean.svg (mean profile)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot heatmap of centred and sorted FFT-aligned nucleosome scores."
    )
    parser.add_argument(
        "--scores",
        required=True,
        help="Input score file (one line = one peak, space-separated values).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PNG file for the heatmap. Defaults to <basename>_heatmap.png",
    )
    parser.add_argument(
        "--zero_thresh",
        type=int,
        default=5,
        help="Consecutive zeros threshold for line skipping.",
    )
    parser.add_argument(
        "--max_score",
        type=float,
        default=300.0,
        help="Max score allowed per value.",
    )

    parser.add_argument(
        "--max_lines",
        type=int,
        default=None,
        help="Maximum number of valid lines to include. "
             "If --subsample-mode=first, reading stops after this many valid rows. "
             "If --subsample-mode=random, all valid rows are loaded then randomly subsampled to this many.",
    )
    parser.add_argument(
        "--subsample-mode",
        choices=["first", "random"],
        default="first",
        help="How to apply --max_lines: "
             "'first' = stop reading after max_lines valid rows (fast); "
             "'random' = load all valid rows then randomly subsample to max_lines.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed used when --subsample-mode=random (for reproducibility).",
    )

    parser.add_argument(
        "--breadth",
        type=float,
        default=1.0,
        help="Fraction of scores to plot from centre (e.g. 0.5 = middle 50%%, +1 adjusted).",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Minimum score for color scale (e.g. -200). Values below will be colored the same.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Maximum score for color scale (e.g. 200). Values above will be colored the same.",
    )
    parser.add_argument(
        "--sort-mode",
        choices=["center", "rise_after_min", "absmean", "max", "unsorted"],
        default="center",
        help="Sorting mode: "
             "'center' = score at x=0 (default), "
             "'rise_after_min' = distance from local min to zero-crossing (tie-breaker: center score), "
             "'absmean' = mean(abs(row)), "
             "'max' = max(row), "
             "'unsorted' = original order.",
    )
    parser.add_argument(
        "--label-replace",
        default="+1 NPS peak",
        help="Replacement text for the default '+1 NPS peak' label.",
    )
    parser.add_argument(
        "--mean_ylim",
        type=float,
        default=None,
        help="Symmetric y-axis limit for mean profile plot (e.g., 40 sets ylim to [-40, 40]).",
    )
    return parser.parse_args()


def has_consecutive_zeros(arr, threshold):
    if threshold <= 0:
        return False
    return np.any(
        np.convolve(arr == 0.0, np.ones(threshold, dtype=int), mode="valid") == threshold
    )


def load_valid_scores(file_path, zero_thresh, max_score,
                      max_lines=None, subsample_mode="first", seed=None):
    valid_rows = []
    skipped = 0
    rng = np.random.default_rng(seed) if seed is not None else None

    with open(file_path) as f:
        for line_num, line in enumerate(
            tqdm(f, desc="Loading", unit="line", dynamic_ncols=True), 1
        ):
            if subsample_mode == "first" and max_lines and len(valid_rows) >= max_lines:
                break

            parts = line.strip().split()
            if not parts:
                skipped += 1
                continue

            try:
                values = np.array([float(x) for x in parts], dtype=np.float64)
            except ValueError:
                skipped += 1
                continue

            if (not np.all(np.isfinite(values))
                    or has_consecutive_zeros(values, zero_thresh)
                    or np.any(values > max_score)):
                skipped += 1
                continue

            valid_rows.append(values)

    if not valid_rows:
        raise RuntimeError("No valid lines found.")

    lengths = {row.shape[0] for row in valid_rows}
    if len(lengths) != 1:
        raise RuntimeError(
            f"Inconsistent row lengths after filtering: {sorted(lengths)}. "
            "All rows must be the same length."
        )

    print(f"[INFO] Loaded {len(valid_rows)} valid rows, skipped {skipped}")

    if subsample_mode == "random" and max_lines is not None and len(valid_rows) > max_lines:
        if rng is None:
            rng = np.random.default_rng()
        idx = rng.choice(len(valid_rows), size=max_lines, replace=False)
        idx.sort()
        valid_rows = [valid_rows[i] for i in idx]
        print(f"[INFO] Randomly subsampled to {len(valid_rows)} rows (seed={seed})")

    return np.vstack(valid_rows)


def compute_sort_key_rise_after_min(row, center_index):
    left = row[: center_index + 1]
    min_index = int(np.argmin(left))
    if row[min_index] >= 0:
        return np.inf

    for i in range(min_index + 1, len(row)):
        if row[i] >= 0:
            return i - min_index

    return np.inf


def plot_heatmap(score_matrix, output_file, breadth=1.0,
                 vmin=None, vmax=None, sort_mode="center",
                 label_replace="+1 NPS peak", mean_ylim=None):
    plt.rcParams.update(
        {
            "axes.titlesize": 24,
            "axes.labelsize": 16,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 16,
            "font.size": 16,
        }
    )

    total_cols = score_matrix.shape[1]
    if not (0.0 < breadth <= 1.0):
        raise ValueError("Breadth must be between 0 and 1.")

    if breadth < 1.0:
        keep_cols = int(total_cols * breadth) + 1
        keep_cols = min(keep_cols, total_cols)
        offset = (total_cols - keep_cols) // 2
        score_matrix = score_matrix[:, offset: offset + keep_cols]
        print(f"[INFO] Cropped to middle {keep_cols} columns (offset={offset})")
    else:
        print(f"[INFO] Using full width ({total_cols} columns)")

    centre_index = score_matrix.shape[1] // 2

    if sort_mode == "center":
        sort_order = np.argsort(-score_matrix[:, centre_index])
        score_matrix = score_matrix[sort_order]
        print("[INFO] Sorted by descending score at center (x=0).")

    elif sort_mode == "rise_after_min":
        distances = np.array(
            [compute_sort_key_rise_after_min(row, centre_index) for row in score_matrix],
            dtype=float,
        )
        center_vals = score_matrix[:, centre_index]

        for idx, (d, cv) in enumerate(zip(distances, center_vals)):
            print(f"[DIST] Row {idx:4d} -> distance={d}, center_score={cv:.3f}")

        sort_order = np.lexsort((-center_vals, distances))
        score_matrix = score_matrix[sort_order]
        print("[INFO] Sorted by distance from local min to zero-crossing, tie-broken by center score.")

    elif sort_mode == "absmean":
        mean_abs = np.mean(np.abs(score_matrix), axis=1)
        sort_order = np.argsort(-mean_abs)
        score_matrix = score_matrix[sort_order]
        print("[INFO] Sorted by descending mean(abs(row)).")

    elif sort_mode == "max":
        row_max = np.max(score_matrix, axis=1)
        sort_order = np.argsort(-row_max)
        score_matrix = score_matrix[sort_order]
        print("[INFO] Sorted by descending max(row)).")

    elif sort_mode == "unsorted":
        print("[INFO] No sorting applied (original order retained).")

    num_cols = score_matrix.shape[1]
    half_len = num_cols // 2
    x_vals = np.arange(-half_len, -half_len + num_cols)

    auto_abs = float(np.max(np.abs(score_matrix))) if score_matrix.size else 1.0
    if auto_abs == 0:
        auto_abs = 1.0

    vmin = -auto_abs if vmin is None else vmin
    vmax = auto_abs if vmax is None else vmax

    print(f"[INFO] Final matrix shape: {score_matrix.shape}")
    print(f"[INFO] Using colormap range: vmin={vmin}, vmax={vmax}")
    print(f"[INFO] Score value range: [{score_matrix.min():.2f}, {score_matrix.max():.2f}]")

    plt.figure(figsize=(15, 5))
    im = plt.imshow(
        score_matrix,
        aspect="auto",
        cmap="seismic",
        vmin=vmin,
        vmax=vmax,
        extent=[x_vals[0], x_vals[-1], 0, score_matrix.shape[0]],
    )

    cbar = plt.colorbar(im)
    cbar.set_label("Nucleosome Protection Score (NPS)", fontsize=16)

    plt.xlabel(f"{label_replace} (bp)")
    plt.ylabel("Region index (sorted)")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"[INFO] Heatmap saved to {output_file}")

    mean_profile = np.mean(score_matrix, axis=0)

    plt.figure(figsize=(15, 5))
    plt.plot(x_vals, mean_profile, color="black", linewidth=1.5)
    plt.xlabel(f"{label_replace} (bp)")
    plt.ylabel("Mean NPS")
    plt.xlim(x_vals[0], x_vals[-1])
    if mean_ylim is not None:
        plt.ylim(-mean_ylim, mean_ylim)
    plt.tight_layout()
    mean_outfile = output_file.replace(".png", "_mean.svg")
    plt.savefig(mean_outfile, dpi=300)
    print(f"[INFO] Mean profile saved to {mean_outfile}")


def main():
    args = parse_args()

    if args.output is None:
        base = os.path.basename(args.scores).rsplit(".", 1)[0]
        args.output = f"{base}_heatmap.png"
        print(f"[INFO] Output file not specified. Using default: {args.output}")

    score_matrix = load_valid_scores(
        args.scores,
        args.zero_thresh,
        args.max_score,
        max_lines=args.max_lines,
        subsample_mode=args.subsample_mode,
        seed=args.seed,
    )

    plot_heatmap(
        score_matrix,
        args.output,
        breadth=args.breadth,
        vmin=args.vmin,
        vmax=args.vmax,
        sort_mode=args.sort_mode,
        label_replace=args.label_replace,
        mean_ylim=args.mean_ylim,
    )


if __name__ == "__main__":
    main()
