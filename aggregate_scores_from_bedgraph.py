#!/usr/bin/env python3

import argparse
import pyBigWig
import numpy as np
import glob
import os
from tqdm import tqdm

def parse_regions(region_file):
    """Reads the region BED-like file, loads all regions (first 4 columns)."""
    regions = []
    with open(region_file, 'r') as f:
        for line in f:
            cols = line.strip().split()[:4]  # Use only the first 4 columns
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"

            strand = cols[3] if len(cols) > 3 else "NA"  # Assume + if column 4 is missing
            center = (start + end) // 2
            regions.append((chrom, start, end, center, strand))
    
    if regions:
        print(f"First region in the BED file: {regions[0]}")
    else:
        print("No regions found in the BED file.")
    
    return regions

def process_bigwig_and_aggregate_scores(bigwig_file, regions, up, down):
    """Processes each region and aggregates average scores at each relative position."""
    total_scores = np.zeros(up + down)
    region_count = 0

    with pyBigWig.open(bigwig_file) as bw:
        # Initialize the tqdm progress bar
        for region in tqdm(regions, desc=f"Processing {os.path.basename(bigwig_file)}", unit="region"):
            chrom, start, end, center, strand = region
            if chrom not in bw.chroms():
                continue  # Skip if the chromosome is not present in the BigWig file

            query_start = max(0, center - up)
            query_end = center + down + 1

            scores = bw.values(chrom, query_start, query_end, numpy=True)

            if scores is None or all(np.isnan(scores)):
                continue  # Skip if no valid scores are found for the region

            # Check for runs of 0s with length >= 5
            zero_run_length = 0
            has_long_zero_run = False

            for score in scores:
                if np.isnan(score) or score == 0:
                    zero_run_length += 1
                else:
                    zero_run_length = 0
                
                if zero_run_length >= 5:
                    has_long_zero_run = True
                    break

            if has_long_zero_run:
                continue  # Skip regions with long runs of 0s

            region_count += 1
            for i, score in enumerate(scores):
                if np.isnan(score):
                    continue  # Skip NaNs (missing data)

                relative_pos = i

                # Flip the position if the strand is negative
                if strand == "-" or (strand == "NA" and np.random.random() < 0.5):
                    relative_pos = up + down - relative_pos - 1

                if 0 <= relative_pos < (up + down):
                    total_scores[relative_pos] += score

    return total_scores, region_count

def write_output(average_scores, output_file, up, down):
    """Writes the average scores at each relative position to the output file."""
    with open(output_file, 'w') as f:
        for i, score in enumerate(average_scores):
            pos = i - up  # Position relative to the center
            f.write(f"{pos}\t{score:.4f}\n")
    print(f"Output written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process multiple BigWig files and aggregate average scores over regions.")
    parser.add_argument("--bigwig_pattern", required=True, help="Glob pattern for the BigWig files (e.g., 'snyder_male_samples_*.wig').")
    parser.add_argument("--regions", required=True, help="Path to the BED-like file containing regions.")
    parser.add_argument("--up", type=int, required=True, help="Number of base pairs upstream to align.")
    parser.add_argument("--down", type=int, required=True, help="Number of base pairs downstream to align.")
    parser.add_argument("--output", required=True, help="Output file for aggregated average scores.")

    args = parser.parse_args()

    # Parse the regions from the BED-like file
    regions = parse_regions(args.regions)

    # Gather all BigWig files matching the pattern
    bigwig_files = glob.glob(args.bigwig_pattern)
    
    if not bigwig_files:
        print("No BigWig files found matching the pattern.")
        return

    # Aggregate scores across all BigWig files
    aggregate_scores = np.zeros(args.up + args.down)
    total_region_count = 0

    for bigwig_file in bigwig_files:
        print(f"Processing {bigwig_file}...")
        scores, region_count = process_bigwig_and_aggregate_scores(bigwig_file, regions, args.up, args.down)
        aggregate_scores += scores
        total_region_count += region_count

    # Calculate final average scores
    if total_region_count > 0:
        final_average_scores = aggregate_scores / total_region_count
    else:
        final_average_scores = np.zeros_like(aggregate_scores)

    # Write the final output
    write_output(final_average_scores, args.output, args.up, args.down)

    print(f"{args.output} count: {total_region_count}")

if __name__ == "__main__":
    main()
