#!/usr/bin/env python3
"""
@author Andrew D Johnston 
@author Fiach Antaw

Fragmentomics scoring + peak calling pipeline.

What this script does (high level):
1) Reads paired-end fragments from one or more BAMs in a region (or whole genome).
2) Filters duplicates (same fragment coords) to max N and optionally subsamples.
3) For each fragment length in a specified range, adds a precomputed "PNS-like" score
   distribution across the fragment (or a padded version if shorter than mode length).
4) Also computes simple coverage and a dyad count (fragment center) track.
5) Smooths the PNS track (Savitzky–Golay).
6) Calls:
   - positive peaks ("nucleosome regions") on smoothed PNS
   - negative peaks ("breakpoint peaks") by flipping the PNS sign and re-calling
7) Writes:
   - one combined bedGraph with multiple score tracks per base
   - a BED-like file of nucleosome regions with peak prominence + coverage metrics
   - a BED-like file of breakpoint peaks with analogous metrics
"""

import sys
from tqdm import tqdm
import argparse
import pysam
import numpy as np
import os
from scipy.signal import savgol_filter
from collections import defaultdict
import random


def generate_paired_reads(bamfile, contig=None, start=None, end=None,
                          max_duplicates=0, subsample=None):
    """
    Iterate over reads in a BAM region and yield properly paired reads as (fwd, rev).

    Key behaviors:
    - Uses query_name to match mates (not relying on BAM being name-sorted).
    - Skips unmapped or mate-unmapped reads.
    - Skips read-pairs mapped to the same strand (invalid for proper PE).
    - Collapses duplicates: the same fragment coordinate tuple is allowed up to
      max_duplicates occurrences.
    - Optional subsample in [0..1]: keep each fragment with probability=subsample.
    """
    _unpaired_reads = dict()                 # store first mate until second appears
    fragment_counts = defaultdict(int)       # duplicate counter by fragment coords
    read_count = 0
    subsample_removed_count = 0

    # NOTE: read_total is computed but not used (left from earlier progress logic)
    try:
        read_total = sum(1 for _ in bamfile.fetch(contig, start, end, multiple_iterators=True))
    except:
        return

    for read in bamfile.fetch(contig, start, end, multiple_iterators=True):

        # Skip unmapped or mate-unmapped reads
        if read.is_unmapped or read.mate_is_unmapped:
            continue

        # Ensure coordinates exist (some pathological records can lack these)
        if read.reference_end is None or read.next_reference_start is None:
            continue

        name = read.query_name

        # First time we see this query_name, store it and continue.
        # Second time, we have the mate and can process the pair.
        if name not in _unpaired_reads:
            _unpaired_reads[name] = read
            continue

        mate = _unpaired_reads[name]
        del _unpaired_reads[name]

        # Skip if both alignments are on same strand: improper / weird mapping
        if read.is_reverse == mate.is_reverse:
            continue

        # Define the fragment by genomic coordinates spanning both mates.
        # This is used to identify duplicate fragments.
        fragment_key = (
            read.reference_name,
            min(read.reference_start, mate.reference_start),
            max(read.reference_end, mate.reference_end),
        )

        # Enforce max_duplicates: allow only N fragments with identical coords
        if fragment_counts[fragment_key] > max_duplicates:
            continue
        fragment_counts[fragment_key] += 1

        # Optional subsampling: keep with probability=subsample
        if subsample is not None and random.random() > subsample:
            subsample_removed_count += 1
            continue

        # Yield in (forward-strand read, reverse-strand read) order for consistency
        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read

        read_count += 1


def generate_fragment_ranges(bamfile, contig, start, end, max_duplicates, subsample):
    """
    Convert paired reads into fragment genomic intervals (frag_start, frag_end).

    Important:
    - Enforces a consistent ordering and sanity checks to avoid inverted coords.
    - Restricts fragments to those whose forward read starts within [start, end).
      (Then coverage later further requires frag_end <= end.)
    """
    for r_fwd, r_rev in generate_paired_reads(bamfile, contig, start, end, max_duplicates, subsample):
        # Ensure r_fwd is the forward-strand read (defensive)
        if r_fwd.is_reverse:
            r_fwd, r_rev = r_rev, r_fwd

        # Sanity checks to prevent weird/inverted pairs
        if r_fwd.reference_start > r_rev.reference_start:
            continue
        if r_rev.reference_end < r_fwd.reference_end:
            continue

        # Only accept fragments whose *start* is inside the region
        # (This avoids emitting fragments starting outside but overlapping in.)
        if r_fwd.reference_start < start or r_fwd.reference_start >= end:
            continue

        # Fragment range is from leftmost start to rightmost end (end is exclusive)
        yield r_fwd.reference_start, r_rev.reference_end


def precompute_distributions(pns_frag_range, mode_DNA_length):
    """
    Precompute per-fragment-length score “kernels” used to add PNS signal.

    For each fragment length:
    - Build a triangle score distribution over the first mode_DNA_length bases, 
    - Mirror it from the fragment end
    - Sum the start and end scores
    - Center by subtracting the mean so each kernel sums ~0.
    """
    distributions = {}

    for fragment_length in pns_frag_range:

        # If fragment is shorter than the mode length, expand the kernel length
        # to accomodate the mode length extending beyond the fragment length at both ends
        if fragment_length < mode_DNA_length:
            total_length = mode_DNA_length + (mode_DNA_length - fragment_length)
        else:
            total_length = fragment_length

        # This midpoint logic is based on the MODE length, not total_length:
        # we create a triangle that rises to 1 at the midpoint of the mode-length window
        # and falls back to 0 at the end of the mode-length window.
        midpoint = (mode_DNA_length - 1) // 2
        second_half_start = midpoint + 1

        scores = np.zeros(total_length)

        # Build the "start" triangle only across the first mode_DNA_length positions.
        # Past that (if total_length > mode_DNA_length), scores remain 0 here.
        for i in range(total_length):
            if i <= midpoint:
                scores[i] = i / midpoint
            elif i <= mode_DNA_length - 1:
                scores[i] = 1 - (i - second_half_start) / midpoint

        # Mirror the distribution from the fragment end and sum:
        end_scores = scores[::-1]
        combined_scores = scores + end_scores

        # Mean-center so the kernel has ~0 mean (prvents baseline drift)
        midpoint_val = np.mean(combined_scores)
        centered_scores = [x - midpoint_val for x in combined_scores]

        distributions[fragment_length] = centered_scores

    return distributions


def score_contig(bamfiles, contig, start, end, mode_DNA_length, pns_frag_range,
                 max_duplicates, distributions, subsample):
    """
    Score a genomic interval [start, end) on one contig.

    Outputs 4 tracks (arrays length = end-start):
    - coverage: +1 for each fragment (in range) for each base it covers
    - dyad: +1 at the fragment center position
    - pns: adds the length-specific precomputed kernel across the fragment footprint
    - pns_smoothed: Savitzky–Golay filtered version of pns

    Coordinates & indexing:
    - Arrays are 0-based relative to 'start' (the adjusted_start for the region).
    - Later code converts these back to absolute genome coordinates when writing outputs.
    """
    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=int)
    dyad = np.zeros(ref_len, dtype=int)
    pns = np.zeros(ref_len, dtype=float)

    filter_frag_count = 0
    total_frag_count = 0

    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(
            bamfile, contig, start, end, max_duplicates, subsample
        ):
            frag_length = frag_end - frag_start
            total_frag_count += 1

            # Kernel for this fragment length (only defined for frag lengths in pns_frag_range)
            fragment_scores = distributions.get(frag_length)

            # Add PNS-like kernel if fragment length in allowed range
            if frag_length in pns_frag_range and fragment_scores:

                # For short fragments, place the (expanded) kernel centered on the
                # fragment center; for long fragments, place it over the fragment.
                if frag_length < mode_DNA_length:
                    total_length = mode_DNA_length + (mode_DNA_length - frag_length)
                    fragment_center = frag_start + (frag_length // 2) - start
                    start_pos = fragment_center - (total_length // 2)
                    end_pos = start_pos + total_length
                else:
                    start_pos = frag_start - start
                    end_pos = frag_end - start

                # Clamp kernel to array bounds (region edges / overlap overhangs)
                if start_pos < 0:
                    fragment_scores = fragment_scores[-start_pos:]
                    start_pos = 0
                if end_pos > ref_len:
                    fragment_scores = fragment_scores[:ref_len - start_pos]

                # Add into pns track
                if 0 <= start_pos < ref_len:
                    pns[start_pos:start_pos + len(fragment_scores)] += fragment_scores[:end_pos - start_pos]

            # Coverage and dyad are only computed for fragments entirely contained
            # in the scoring region; also restricted to the allowed length range.
            if frag_start >= start and frag_end <= end:
                if frag_length in pns_frag_range:
                    coverage[frag_start - start:frag_end - start] += 1
                    filter_frag_count += 1

                    fragment_center = frag_start + (frag_length // 2) - start
                    if 0 <= fragment_center < ref_len:
                        dyad[fragment_center] += 1

    # Smooth PNS to reduce noise before peak calling
    window_size = 21
    polyorder = 2
    pns_smoothed = savgol_filter(pns, window_size, polyorder)

    # Package as tuples: (contig, start_of_array_in_genome_coords, array)
    # Note: start here is the *adjusted_start* used for this region.
    scores = {
        'coverage': [(contig, start, coverage)],
        'pns_smoothed': [(contig, start, pns_smoothed)],
        'pns': [(contig, start, pns)],
        'dyad': [(contig, start, dyad)],
    }

    return scores, pns_frag_range


def find_peaks_and_regions(scores, original_start, min_length=50, max_neg_run=5):
    """
    Segment the 1D score array into alternating positive and negative regions, then call peaks.

    Logic:
    - First, find "positive regions": consecutive scores > 0, allowing brief dips (<= max_neg_run)
      before closing the region. Regions shorter than min_length are discarded.
    - Between successive positive regions, find the most negative point => negative peak.
    - Within each positive region, find the max point => positive peak.
    - Also define a "region midpoint" peak: midpoint of region boundaries (not max position).
      (This is stored as adjusted_positive_peaks, used as a PNS peak coordinate.)

    Coordinate handling:
    - Input 'scores' uses array indices [0..len-1].
    - This function converts peaks/regions to absolute genome coordinates by adding original_start.
      (In calling code, original_start is the region's adjusted_start in genome coords.)
    """
    positive_regions = []
    current_region = None
    searching_for_positive = True  # start in positive scanning mode

    for i in range(len(scores)):
        score = scores[i]

        if searching_for_positive:
            if score > 0:
                if current_region is None:
                    current_region = [i, i]
                else:
                    current_region[1] = i
                neg_count = 0
            else:
                # We were in a positive region but now hit <=0.
                # Allow up to max_neg_run consecutive non-positive positions before closing.
                if current_region:
                    neg_count += 1
                    if neg_count == 1:
                        last_positive_end = i - 1

                    if neg_count >= max_neg_run:
                        if last_positive_end - current_region[0] + 1 >= min_length:
                            positive_regions.append([current_region[0], last_positive_end])
                        current_region = None
                        searching_for_positive = False
        else:
            # Searching for negative region (we don't store negative regions; just skip until positive resumes)
            if score <= 0:
                if current_region is None:
                    current_region = [i, i]
                else:
                    current_region[1] = i
            else:
                current_region = [i, i]
                searching_for_positive = True
                neg_count = 0

    # If we ended while inside a positive region, add it if long enough
    if current_region:
        if searching_for_positive and current_region[1] - current_region[0] + 1 >= min_length:
            positive_regions.append(current_region)

    # Negative peaks: for each gap between positive regions, take argmin in the inter-region segment
    negative_peaks = []
    negative_peak_scores = []
    for i in range(1, len(positive_regions)):
        prev_end = positive_regions[i - 1][1]
        next_start = positive_regions[i][0]
        inter_region_scores = scores[prev_end + 1:next_start]
        if len(inter_region_scores) > 0:
            most_negative_index = np.argmin(inter_region_scores) + prev_end + 1
            most_negative_score = scores[most_negative_index]
            negative_peaks.append(most_negative_index)
            negative_peak_scores.append(most_negative_score)

    # Positive peaks: max within each positive region
    positive_peaks = []
    positive_peak_scores = []
    for region in positive_regions:
        region_scores = scores[region[0]:region[1] + 1]
        peak_index = np.argmax(region_scores) + region[0]
        positive_peaks.append(peak_index)
        positive_peak_scores.append(scores[peak_index])

    # Convert positive regions to absolute genomic coords
    positive_peak_regions = [(region[0] + original_start, region[1] + original_start)
                             for region in positive_regions]

    # "Adjusted" peak is region midpoint in absolute coords (stable center-of-region peak)
    adjusted_positive_peaks = [(region[0] + region[1]) // 2 + original_start
                               for region in positive_regions]

    # Convert peak indices to absolute coords
    positive_peaks = [p + original_start for p in positive_peaks]
    negative_peaks = [p + original_start for p in negative_peaks]

    return (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (adjusted_positive_peaks, positive_peak_regions),
    )


def write_bedgraph(scores, contigs, out_prefix, first_region=False):
    """
    Write a multi-column bedGraph-like file with per-base scores for multiple tracks.

    Output format per line:
        chrom  start  end  coverage  pns_smoothed  pns  dyad
    (exact order depends on scores.keys() iteration order)

    Important coordinate note:
    - scores arrays are relative to the region's adjusted_start.
    - contigs contains (original_start, original_end) boundaries for the *non-overhang* region.
      This function only writes positions within [original_start, original_end), trimming overlaps.
    """
    mode = 'w' if first_region else 'a'
    original_start, original_end = contigs[0]
    filename = f"{out_prefix}_combined_scores.bedGraph"

    with open(filename, mode) as f:
        # Use the first score type to establish (contig, start, length)
        first_score_type = list(scores.keys())[0]
        contig_data = scores[first_score_type]

        for contig, start, score_array in contig_data:
            if not contig.startswith("chr"):
                contig = f"chr{contig}"

            for i in range(len(score_array)):
                position = start + i  # absolute genome coordinate
                if original_start <= position < original_end:
                    line = [contig, str(position), str(position + 1)]

                    # Append each score type’s value at this position
                    for stype in scores.keys():
                        stype_score_array = scores[stype][0][2]
                        if i < len(stype_score_array):
                            # Write ints cleanly, floats to 2 decimals
                            if stype_score_array[i] == int(stype_score_array[i]):
                                line.append(f"{int(stype_score_array[i])}")
                            else:
                                line.append(f"{stype_score_array[i]:.2f}")
                        else:
                            line.append("0")

                    f.write("\t".join(line) + "\n")


def write_nucleosome_peaks(peaks, contigs, out_prefix,
                          first_region=False, flip_scores=False):
    """
    Write called peaks and their nucleosome/breakpoint regions into a BED-like file.

    Each output row contains:
    chrom, region_start, region_end,
    prominence_score,
    adjusted_peak (midpoint),
    upstream_score, upstream_peak_pos,
    downstream_score, downstream_peak_pos,
    peak_score, peak_pos,
    max_coverage, max_coverage_pos

    Prominence is computed as:
        peak_score - mean(upstream_flank_score, downstream_flank_score)
    where flank scores come from nearest negative peaks on either side.
    For flipped score calling (breakpoints), we flip sign so "prominence" is still positive-ish.
    """
    mode = 'w' if first_region else 'a'
    original_start, original_end = contigs[0]

    nucleosome_filename = f"{out_prefix}.bed"
    with open(nucleosome_filename, mode) as f:
        for (contig, original_start), peak_data in peaks.items():
            if not contig.startswith("chr"):
                contig = f"chr{contig}"

            num_positive_peaks = len(peak_data['adjusted_peaks'])
            num_negative_peaks = len(peak_data['negative_peaks'])

            for i in range(num_positive_peaks):
                adjusted_peak = peak_data['adjusted_peaks'][i]
                nucleosome_region_start = peak_data['nucleosome_regions'][i][0]
                nucleosome_region_end = peak_data['nucleosome_regions'][i][1]
                peak = peak_data['positive_peaks'][i]

                # Only write peaks whose max-position lies inside the non-overhang region
                if not (original_start <= peak < original_end):
                    continue

                # Find closest upstream and downstream negative peaks around this positive peak
                upstream_index = None
                downstream_index = None

                for j in range(num_negative_peaks):
                    if peak_data['negative_peaks'][j] < peak:
                        upstream_index = j
                    else:
                        break

                for j in range(num_negative_peaks):
                    if peak_data['negative_peaks'][j] > peak:
                        downstream_index = j
                        break

                # If no flanking negative peak exists on a side, fall back to the positive peak
                # (so prominence becomes 0 on that side)
                if upstream_index is not None:
                    upstream_negative_peak = peak_data['negative_peaks'][upstream_index]
                    upstream_score = peak_data['negative_peak_scores'][upstream_index]
                else:
                    upstream_negative_peak = peak
                    upstream_score = peak_data['positive_peak_scores'][i]

                if downstream_index is not None:
                    downstream_negative_peak = peak_data['negative_peaks'][downstream_index]
                    downstream_score = peak_data['negative_peak_scores'][downstream_index]
                else:
                    downstream_negative_peak = peak
                    downstream_score = peak_data['positive_peak_scores'][i]

                # For the "breakpoint" call track, scores were negated to turn troughs into peaks.
                # flip_scores flips sign back for reporting so trough depth is represented consistently.
                if flip_scores:
                    upstream_score *= -1
                    downstream_score *= -1
                    peak_score = -1 * peak_data['positive_peak_scores'][i]
                else:
                    peak_score = peak_data['positive_peak_scores'][i]

                average_flanking_score = np.mean([upstream_score, downstream_score])
                prominence_score = peak_score - average_flanking_score

                max_coverage = peak_data['max_coverages'][i]
                max_position = peak_data['max_positions'][i]

                f.write(
                    f"{contig}\t{nucleosome_region_start}\t{nucleosome_region_end}\t"
                    f"{prominence_score:.2f}\t{adjusted_peak}\t"
                    f"{upstream_score:.2f}\t{upstream_negative_peak}\t"
                    f"{downstream_score:.2f}\t{downstream_negative_peak}\t"
                    f"{peak_score:.2f}\t{peak}\t"
                    f"{max_coverage}\t{max_position}\n"
                )


def split_into_regions(contig, start, end, max_length=100000, overlap=1000):
    """
    Split a long contig interval into windows for memory/time efficiency.

    Each region includes:
      (contig,
       adjusted_start, adjusted_end,   # includes 'overlap' padding on both sides
       original_start, original_end)   # the "true" non-overhang window to keep outputs for

    The overlap padding ensures that smoothing and peak-calling near edges of a window
    still has context from neighboring bases; outputs are then trimmed back to original_*
    in write_bedgraph and write_nucleosome_peaks.
    """
    regions = []
    current_start = start
    while current_start < end:
        adjusted_start = max(current_start - overlap, start - overlap, 0)
        adjusted_end = min(current_start + max_length + overlap, end + overlap)
        original_start = current_start
        original_end = min(current_start + max_length, end)
        regions.append((contig, adjusted_start, adjusted_end, original_start, original_end))
        current_start = original_end
    return regions


def call_and_write_peaks(
    scores,
    coverage_scores,
    adjusted_start,   # array index 0 corresponds to this absolute coordinate
    original_start,
    original_end,
    contig,
    out_prefix,
    first_region,
    peak_type_label,
    flip_scores,
):
    """
    Peak calling + output writing wrapper for one score track over one window.

    Inputs:
    - scores: 1D numpy array for the window (relative to adjusted_start)
    - coverage_scores: 1D int array matching scores length (relative to adjusted_start)
    - adjusted_start: absolute coordinate of index 0 in these arrays
    - original_start/original_end: boundaries of the non-overlap part of this window

    Process:
    1) Call positive regions/peaks and inter-region negative peaks via find_peaks_and_regions().
       That function returns absolute genomic coordinates (because we pass adjusted_start).
    2) For each called region, compute max coverage and its position within the region.
    3) Write peaks/regions to output BED-like file (append after first region).
    """
    positive_peaks, negative_peaks, adjusted_peaks = find_peaks_and_regions(scores, adjusted_start, 50, 5)

    max_coverages = []
    max_positions = []
    arr_len = coverage_scores.shape[0]

    # adjusted_peaks[1] contains region intervals in absolute coords (start,end)
    for start, end in adjusted_peaks[1]:
        # Convert absolute coords back to window-relative indices
        region_start_idx = max(0, start - adjusted_start)
        region_end_idx = min(arr_len - 1, end - adjusted_start)

        if region_end_idx < region_start_idx:
            max_coverages.append(0)
            max_positions.append(0)
            continue

        # Slice is inclusive, so add 1 to end index
        region_coverage = coverage_scores[region_start_idx:region_end_idx + 1]

        if region_coverage.size > 0:
            local_argmax = int(np.argmax(region_coverage))
            max_coverages.append(int(region_coverage[local_argmax]))
            max_positions.append(region_start_idx + local_argmax + adjusted_start)
        else:
            max_coverages.append(0)
            max_positions.append(0)

    peaks = {
        (contig, original_start): {
            'positive_peaks': positive_peaks[0],
            'positive_peak_scores': positive_peaks[1],
            'negative_peaks': negative_peaks[0],
            'negative_peak_scores': negative_peaks[1],
            'adjusted_peaks': adjusted_peaks[0],
            'nucleosome_regions': adjusted_peaks[1],
            'max_coverages': max_coverages,
            'max_positions': max_positions,
        }
    }

    write_nucleosome_peaks(
        peaks,
        [(original_start, original_end)],
        out_prefix + peak_type_label,
        first_region,
        flip_scores,
    )


def require_bam_indexes(bam_paths, parser=None):
    """
    Ensure each BAM has an accompanying BAI index in the same directory.
    Exits with an error if any index is missing.

    Accepts either:
      - <bam>.bai
      - <bam_basename>.bai   (e.g. sample.bai for sample.bam)

    Args:
        bam_paths (list[str]): BAM file paths.
        parser (argparse.ArgumentParser | None): If provided, uses parser.error()
            for consistent CLI error formatting; otherwise raises FileNotFoundError.
    """
    missing = []

    for bam in bam_paths:
        bam = os.path.abspath(bam)
        bam_dir = os.path.dirname(bam)

        # Two common index naming conventions
        idx1 = bam + ".bai"  # e.g. sample.bam.bai
        idx2 = os.path.join(bam_dir, os.path.splitext(os.path.basename(bam))[0] + ".bai")  # e.g. sample.bai

        if not (os.path.exists(idx1) or os.path.exists(idx2)):
            missing.append((bam, idx1, idx2))

    if missing:
        msg_lines = ["ERROR: Missing BAM index (.bai) for the following BAM file(s).",
                     "Each BAM input must have its .bai index in the SAME directory as the BAM:",
                     ""]
        for bam, idx1, idx2 in missing:
            msg_lines.append(f"  BAM: {bam}")
            msg_lines.append(f"    expected: {idx1}")
            msg_lines.append(f"         or : {idx2}")
            msg_lines.append("")
        msg_lines.append("Create indexes with:")
        msg_lines.append("  samtools index <file.bam>")
        msg = "\n".join(msg_lines)

        if parser is not None:
            parser.error(msg)
        else:
            raise FileNotFoundError(msg)


def main():
    parser = argparse.ArgumentParser(description="Score fragmentomics data.")
    parser.add_argument('-b', '--bamfiles', nargs='+', help='BAM file(s) to process')
    parser.add_argument('-o', '--out_prefix', help='prefix for output files (default: based on BAM names and contigs)')
    parser.add_argument('-c', '--contigs', nargs='+',
                        help='limit to contig(s) and optional range, e.g. "2:100000-200000"')
    parser.add_argument('--mode-length', type=int, default=167, help='Mode fragment length (used in kernel geometry)')
    parser.add_argument('--frag-lower', type=int, default=127, help='Lower fragment length to include')
    parser.add_argument('--frag-upper', type=int, default=207, help='Upper fragment length to include')
    parser.add_argument('--max-duplicates', type=int, default=0, help='Max allowed duplicate fragments with same coords')
    parser.add_argument('--subsample', type=float, default=None,
        help='Subsampling proportion (e.g., 0.5 to subsample 50%% of the reads)')

    args = parser.parse_args()

    require_bam_indexes(args.bamfiles, parser=parser)

    # Default output prefix: join BAM basenames, optionally append contig if single contig
    if not args.out_prefix:
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in args.bamfiles]
        args.out_prefix = f"{'_'.join(bam_basenames)}"
        if args.contigs and len(args.contigs) == 1:
            safe_contig = args.contigs[0].replace(":", "_")
            args.out_prefix = f"{args.out_prefix}_{safe_contig}"

    # Add mode and frag range to make outputs self-describing
    args.out_prefix = f"{args.out_prefix}_mode{args.mode_length}_lower{args.frag_lower}_upper{args.frag_upper}"

    # Open BAMs
    bamfiles = []
    for bamfile_path in args.bamfiles:
        try:
            bamfile = pysam.AlignmentFile(bamfile_path, "rb")
            bamfiles.append(bamfile)
        except FileNotFoundError:
            parser.error(f"Unable to open bamfile {bamfile_path} (file not found)")
            return -1
        except Exception as e:
            parser.error(f"Unable to open bamfile {bamfile_path}: {str(e)}")
            return -1

    # Build list of regions to process:
    # - If user provides contigs with ranges: split into 100kb windows with overlaps
    # - If user provides contigs without ranges: whole contig split into windows
    # - If nothing provided: process all references in BAM[0]
    contigs = []
    if args.contigs:
        for contig_range in args.contigs:
            if ':' in contig_range:
                contig, positions = contig_range.split(':')
                start, end = map(int, positions.split('-'))
                contigs.extend(split_into_regions(contig, start, end))
            else:
                contig = contig_range
                start, end = 0, bamfiles[0].get_reference_length(contig)
                contigs.extend(split_into_regions(contig, start, end))
    else:
        for contig in bamfiles[0].references:
            start, end = 0, bamfiles[0].get_reference_length(contig)
            contigs.extend(split_into_regions(contig, start, end))

    # Precompute scoring kernels for each fragment length in [frag_lower, frag_upper]
    pns_frag_range = range(args.frag_lower, args.frag_upper + 1)
    distributions = precompute_distributions(pns_frag_range, args.mode_length)

    # Remove old outputs (append mode is used across windows)
    combined_bedgraph = f"{args.out_prefix}_combined_scores.bedGraph"
    nucleosome_bed = f"{args.out_prefix}_nucleosome_regions.bed"
    breakpoint_bed = f"{args.out_prefix}_breakpoint_peaks.bed"
    for fname in [combined_bedgraph, nucleosome_bed, breakpoint_bed]:
        if os.path.exists(fname):
            os.remove(fname)

    first_region = True
    for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc='Scoring contigs'):

        # Compute tracks on the overlapped window [adjusted_start, adjusted_end)
        scores, pns_frag_range = score_contig(
            bamfiles,
            contig,
            adjusted_start,
            adjusted_end,
            args.mode_length,
            pns_frag_range,
            args.max_duplicates,
            distributions,
            args.subsample
        )

        # Extract arrays (relative to adjusted_start)
        pns_smoothed_scores = scores['pns_smoothed'][0][2]
        coverage_scores = scores['coverage'][0][2]

        # Call positive peaks => nucleosome regions
        call_and_write_peaks(
            scores=pns_smoothed_scores,
            coverage_scores=coverage_scores,
            adjusted_start=adjusted_start,
            original_start=original_start,
            original_end=original_end,
            contig=contig,
            out_prefix=args.out_prefix,
            first_region=first_region,
            peak_type_label="_nucleosome_regions",
            flip_scores=False,
        )

        # Call negative peaks by flipping sign => breakpoint peaks
        flipped_scores = -1 * pns_smoothed_scores
        call_and_write_peaks(
            scores=flipped_scores,
            coverage_scores=coverage_scores,
            adjusted_start=adjusted_start,
            original_start=original_start,
            original_end=original_end,
            contig=contig,
            out_prefix=args.out_prefix,
            first_region=first_region,
            peak_type_label="_breakpoint_peaks",
            flip_scores=True,
        )

        # Write per-base scores, trimmed to the non-overlap region
        write_bedgraph(scores, [(original_start, original_end)], args.out_prefix, first_region)

        first_region = False

    return 0


if __name__ == '__main__':
    sys.exit(main())
