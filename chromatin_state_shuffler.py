# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:04:05 2019

@author: Andrew D Johnston
"""

# -*- coding: utf-8 -*-
"""
Shuffle ChromHMM segments within each chromosome while preserving:
  - The set of segment lengths per chromosome
  - The set of compartment/state labels per chromosome
  - (Length,label) pairing stays intact (we shuffle tuples)

How it works (high level):
1) Read the input BED. For each chromosome:
   - Collect all (segment_length, state_label) tuples into D_compartments[chrom]
   - Create a list of candidate "start anchors" every 200 bp inside all non-NA segments
     (these are potential positions to place shuffled segments) into D_regions[chrom]
2) For each chromosome:
   - Shuffle the list of (segment_length, state_label)
   - Walk the candidate anchors in increasing order and place the next shuffled segment
     as long as it doesn't overlap the previously placed segment.
3) Write the resulting shuffled segments to an output BED.

"""

from random import shuffle

# For each chromosome, holds a list of (segment_length_bp, state_label) tuples
D_compartments = {}

# For each chromosome, holds a list of candidate anchor starts (chrom, start_pos)
# Starts are sampled every 200 bp within each non-NA segment.
D_regions = {}

INPUT_BED = "wgEncodeBroadHmmGm12878HMM.bed"
OUTPUT_BED = "wgEncodeBroadHmmGm12878HMM_shuffled.bed"
ANCHOR_STEP = 200  # bp spacing between candidate anchor starts

# ----------------------------
# Read input BED and build per-chromosome structures
# ----------------------------
with open(INPUT_BED, "r") as F_chromatin:
    for raw_line in F_chromatin:
        fields = raw_line.rstrip("\n").split("\t")

        # Expect BED4: chrom, start, end, state/compartment
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        comp = fields[3]

        # Skip unannotated regions
        if comp == "NA":
            continue

        seg_len = end - start
        if seg_len <= 0:
            # Defensive: skip malformed/empty intervals
            continue

        # Add this segment's (length,label) to the per-chromosome pool that we'll shuffle
        if chrom in D_compartments:
            D_compartments[chrom].append((seg_len, comp))
        else:
            D_compartments[chrom] = [(seg_len, comp)]

        # Create candidate anchor start positions every ANCHOR_STEP bp inside this segment
        # e.g. seg_len=1000 -> anchors at start+0, start+200, start+400, start+600, start+800
        n_anchors = seg_len // ANCHOR_STEP
        if n_anchors <= 0:
            continue

        if chrom not in D_regions:
            D_regions[chrom] = []

        for x in range(n_anchors):
            anchor_start = start + ANCHOR_STEP * x
            D_regions[chrom].append((chrom, anchor_start))

# ----------------------------
# Place shuffled segments and write output
# ----------------------------
with open(OUTPUT_BED, "w") as F_out:
    for chrom in D_regions:
        # If we somehow have anchors but no segments, skip
        if chrom not in D_compartments or not D_compartments[chrom]:
            continue

        # Randomize the order of (length,label) pairs for this chromosome
        shuffle(D_compartments[chrom])

        # Track the end coordinate of the most recently written segment
        # Initialize with the first anchor start so the first placement is allowed.
        prev_end = D_regions[chrom][0][1]

        # Walk anchors in order; place the next shuffled segment at the first
        # anchor that does not overlap the previous placed segment.
        for region in D_regions[chrom]:
            anchor_start = region[1]

            # Only place if this anchor is at/after the end of the last segment
            if anchor_start < prev_end:
                continue

            # If we've placed all shuffled segments, stop for this chromosome
            if not D_compartments[chrom]:
                break

            seg_len, state_label = D_compartments[chrom][0]
            seg_end = anchor_start + seg_len

            # Write BED line: chrom, start, end, label
            F_out.write(f"{chrom}\t{anchor_start}\t{seg_end}\t{state_label}\n")

            # Update end of last placed segment and consume the used (len,label)
            prev_end = seg_end
            D_compartments[chrom].pop(0)

                    

