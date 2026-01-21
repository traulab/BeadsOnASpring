#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Andrew D Johnston

Shuffle ChromHMM segments within each chromosome while preserving:
  - The set of segment lengths per chromosome
  - The set of compartment/state labels per chromosome
  - (Length,label) pairing stays intact (we shuffle tuples)

How it works (high level):
1) Read the input BED. For each chromosome:
   - Collect all (segment_length, state_label) tuples into D_compartments[chrom]
   - Create a list of candidate "start anchors" every --anchor-step bp inside all non-NA segments
     (these are potential positions to place shuffled segments) into D_regions[chrom]
2) For each chromosome:
   - Shuffle the list of (segment_length, state_label)
   - Walk the candidate anchors in increasing order and place the next shuffled segment
     as long as it doesn't overlap the previously placed segment.
3) Write the resulting shuffled segments to an output BED.

Notes:
- Input is assumed BED4: chrom, start, end, label
- Rows with label == --skip-label (default: NA) are ignored
"""

from __future__ import annotations

import argparse
import sys
import random
from collections import defaultdict


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Shuffle ChromHMM BED segments within each chromosome using anchor starts."
    )
    p.add_argument(
        "-i", "--input-bed",
        default="wgEncodeBroadHmmGm12878HMM.bed",
        help="Input BED4 file (default: %(default)s)"
    )
    p.add_argument(
        "-o", "--output-bed",
        default="wgEncodeBroadHmmGm12878HMM_shuffled.bed",
        help="Output BED file (default: %(default)s)"
    )
    p.add_argument(
        "--anchor-step",
        type=int,
        default=200,
        help="Spacing (bp) between candidate anchor starts within annotated regions (default: %(default)s)"
    )
    p.add_argument(
        "--skip-label",
        default="NA",
        help="Label to skip (unannotated regions) (default: %(default)s)"
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible shuffling (default: none)"
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.anchor_step <= 0:
        raise ValueError(f"--anchor-step must be > 0 (got {args.anchor_step})")

    if args.seed is not None:
        random.seed(args.seed)

    # For each chromosome, holds a list of (segment_length_bp, state_label) tuples
    D_compartments: dict[str, list[tuple[int, str]]] = defaultdict(list)

    # For each chromosome, holds candidate anchor starts (chrom, start_pos)
    D_regions: dict[str, list[int]] = defaultdict(list)

    # ----------------------------
    # Read input BED and build per-chromosome structures
    # ----------------------------
    with open(args.input_bed, "r") as f:
        for lineno, raw_line in enumerate(f, start=1):
            s = raw_line.rstrip("\n")
            if not s or s.startswith("#"):
                continue

            fields = s.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{args.input_bed}: line {lineno} has <4 columns; expected BED4 (chrom,start,end,label)"
                )

            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise ValueError(f"{args.input_bed}: line {lineno} start/end not integer: {s}")

            comp = fields[3]

            if comp == args.skip_label:
                continue

            seg_len = end - start
            if seg_len <= 0:
                # Malformed or empty interval; skip defensively
                continue

            # Keep (len,label) pairing intact, but shuffle ordering later
            D_compartments[chrom].append((seg_len, comp))

            # Candidate anchors every anchor_step bp, starting at segment start
            # seg_len=1000, step=200 -> anchors at start+0, +200, +400, +600, +800
            n_anchors = seg_len // args.anchor_step
            if n_anchors <= 0:
                continue

            for x in range(n_anchors):
                D_regions[chrom].append(start + args.anchor_step * x)

    # ----------------------------
    # Place shuffled segments and write output
    # ----------------------------
    with open(args.output_bed, "w") as out:
        for chrom in D_regions:
            if chrom not in D_compartments or not D_compartments[chrom]:
                continue

            anchors = D_regions[chrom]
            if not anchors:
                continue

            # Ensure anchors are in increasing order (they should be, but keep it robust)
            anchors.sort()

            # Shuffle the (len,label) pairs
            segs = D_compartments[chrom]
            random.shuffle(segs)

            prev_end = anchors[0]

            # Walk anchors; place next segment at first anchor that doesn't overlap previous
            seg_idx = 0
            for anchor_start in anchors:
                if anchor_start < prev_end:
                    continue
                if seg_idx >= len(segs):
                    break

                seg_len, state_label = segs[seg_idx]
                seg_end = anchor_start + seg_len

                out.write(f"{chrom}\t{anchor_start}\t{seg_end}\t{state_label}\n")

                prev_end = seg_end
                seg_idx += 1

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        # Allows piping to head, etc.
        raise SystemExit(141)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(1)



