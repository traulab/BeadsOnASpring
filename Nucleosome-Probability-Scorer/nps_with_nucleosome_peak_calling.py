import sys
from tqdm import tqdm
import argparse
import pysam
import numpy as np
import os
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from collections import defaultdict
import random 

def generate_paired_reads(bamfile, contig=None, start=None, end=None, max_duplicates=50, subsample=None):
    """Generator taking paired-end data from a sam/bam file and yielding both pairs as (fwd, rev),
       allowing a specified number of duplicate fragments and applying optional subsampling."""
    _unpaired_reads = dict()
    fragment_counts = defaultdict(int)
    read_count = 0
    subsample_removed_count = 0
    try:
        read_total = sum(1 for _ in bamfile.fetch(contig, start, end, multiple_iterators=True))
    except:
        return
        
    duplicates_per_100k = read_total/100000 * max_duplicates

    for read in bamfile.fetch(contig, start, end, multiple_iterators=True):

        if read.is_unmapped or read.mate_is_unmapped:
            continue  # Skip unmapped reads

        # Ensure that both reference_end and mate.reference_end are valid
        if read.reference_end is None or read.next_reference_start is None:
            continue  # Skip reads with missing coordinates

        name = read.query_name

        # Pair up the reads by query name
        if name not in _unpaired_reads:
            _unpaired_reads[name] = read
            continue

        # Have found mate, yield both
        mate = _unpaired_reads[name]
        del _unpaired_reads[name]

        # Skip if both pairs are mapped to the same strand
        if read.is_reverse == mate.is_reverse:
            continue

        # Define the fragment by its contig and start/end positions
        fragment_key = (read.reference_name, min(read.reference_start, mate.reference_start), 
                        max(read.reference_end, mate.reference_end))

        # Track duplicate counts for this fragment
        if fragment_counts[fragment_key] >= max(1, duplicates_per_100k):
            continue

        fragment_counts[fragment_key] += 1

        # Apply subsampling if specified
        if subsample is not None and random.random() > subsample:
            subsample_removed_count += 1
            continue  # Skip this read if the random number is greater than the subsample value

        # Ensure the pairs are in the right order: (fwd, rev)
        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read

        read_count += 1

    # print(f"Processing {read_count}, removed {subsample_removed_count} reads due to {str(subsample)} subsampling")

def generate_fragment_ranges(bamfile, contig, start, end, max_duplicates, subsample):
    """Generate ranges covered by sequenced DNA fragments with optional subsampling."""
    for r_fwd, r_rev in generate_paired_reads(bamfile, contig, start, end, max_duplicates, subsample):
        if r_fwd.is_reverse:
            r_fwd, r_rev = r_rev, r_fwd  # Ensure r_fwd is always the read that maps to the forward strand

        if r_fwd.reference_start > r_rev.reference_start:
            continue
        if r_rev.reference_end < r_fwd.reference_end:
            continue
        if r_fwd.reference_start < start or r_fwd.reference_start >= end:
            continue
        yield r_fwd.reference_start, r_rev.reference_end

def normalize_scores(scores, window_size):
    """
    Normalizes the scores using a sliding window of specified size.
    """
    half_window = window_size // 2
    smoothed_scores = np.zeros_like(scores)
    total_length = len(scores)
    
    for i in range(total_length):
        window_start = max(i - half_window, 0)
        window_end = min(i + half_window + 1, total_length)
        window_sum = np.sum(scores[window_start:window_end])
        window_count = window_end - window_start
        if window_count > 0:
            smoothed_scores[i] = scores[i] - (window_sum / window_count)
            
    return smoothed_scores

def precompute_distributions(nps_frag_range, mode_DNA_length):
    """
    Precompute fragment score distributions for each fragment length in the range
    and store them in a dictionary.
    """
    distributions = {}
    
    for fragment_length in nps_frag_range:
        # Create the distribution for the current fragment length
        if fragment_length < mode_DNA_length:
            total_length = mode_DNA_length + (mode_DNA_length - fragment_length)
        else:
            total_length = fragment_length

        # Calculate dynamic midpoint based on mode_DNA_length
        midpoint = (mode_DNA_length - 1) // 2
        second_half_start = midpoint + 1

        scores = np.zeros(total_length)

        for i in range(total_length):
            if i <= midpoint:
                scores[i] = i / midpoint  # Scale up the score to peak at 1
            elif i <= mode_DNA_length - 1:
                scores[i] = 1 - (i - second_half_start) / midpoint  # Adjust to descend back to 0

        # Second distribution: flip the first distribution
        end_scores = scores[::-1]  # Reverse the whole scores array

        # Sum the two distributions
        combined_scores = scores + end_scores

        # Center the distribution around y-axis
        midpoint_val = np.mean(combined_scores)
        centered_scores = [x - midpoint_val for x in combined_scores]

        # Store the centered scores in the dictionary
        distributions[fragment_length] = centered_scores

    return distributions

def score_contig(bamfiles, contig, start, end, mode_DNA_length, nps_frag_range, max_duplicates, distributions, subsample):
    """
    Applies a range of scoring algorithms to a given reference contig, returning the results.
    Also returns fragment start positions and lengths for further plotting.
    """
    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=int)
    dyad = np.zeros(ref_len, dtype=int)
    nps = np.zeros(ref_len, dtype=float)
    fragment_starts = []
    fragment_lengths = []
    filter_frag_count = 0
    total_frag_count = 0

    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(bamfile, contig, start, end, max_duplicates, subsample):

            frag_length = frag_end - frag_start
            fragment_starts.append(frag_start)  # Collect the fragment start position
            fragment_lengths.append(frag_length)  # Collect the fragment length
            total_frag_count += 1

            # Retrieve precomputed fragment_scores from the dictionary
            fragment_scores = distributions.get(frag_length)

            if frag_length in nps_frag_range and fragment_scores:
                if frag_length < mode_DNA_length:
                    total_length = mode_DNA_length + (mode_DNA_length - frag_length)
                    fragment_center = frag_start + (frag_length // 2) - start
                    start_pos = fragment_center - (total_length // 2)
                    end_pos = start_pos + total_length
                else:
                    start_pos = frag_start - start
                    end_pos = frag_end - start

                if start_pos < 0:
                    fragment_scores = fragment_scores[-start_pos:]
                    start_pos = 0
                if end_pos > ref_len:
                    fragment_scores = fragment_scores[:ref_len - start_pos]

                if 0 <= start_pos < ref_len:
                    nps[start_pos:start_pos + len(fragment_scores)] += fragment_scores[:end_pos - start_pos]

            # Calculate coverage
            if frag_start >= start and frag_end <= end:
                if frag_length in nps_frag_range:
                    coverage[frag_start - start:frag_end - start] += 1
                    filter_frag_count += 1

                    # Add +1 to the center position of the fragment for dyad scoring
                    fragment_center = frag_start + (frag_length // 2) - start
                    if 0 <= fragment_center < ref_len:
                        dyad[fragment_center] += 1

    # print(filter_frag_count)
    # print(total_frag_count)

    # Apply Savitzky-Golay filter
    window_size = 21
    polyorder = 2
    nps_smoothed = savgol_filter(nps, window_size, polyorder)

    # Return scores as tuples in the format (contig, start, score_array)
    scores = {
        'coverage': [(contig, start, coverage)],
        'nps_smoothed': [(contig, start, nps_smoothed)],
        'nps': [(contig, start, nps)],
        'dyad': [(contig, start, dyad)],
    }

    return scores, nps_frag_range

def find_peaks_and_regions(scores, original_start, min_length=50, max_neg_run=5):
    positive_regions = []
    current_region = None
    searching_for_positive = True  # Start by looking for positive regions

    for i in range(len(scores)):
        score = scores[i]

        if searching_for_positive:
            if score > 0:
                if current_region is None:
                    current_region = [i, i]  # Start a new positive region
                else:
                    current_region[1] = i  # Extend the current positive region
                neg_count = 0  # Reset the negative score count
            else:
                if current_region:
                    neg_count += 1
                    if neg_count == 1:
                        last_positive_end = i - 1  # Mark the point where the score first dips to zero or below

                    if neg_count >= max_neg_run:
                        # End the current positive region at the last positive score
                        if last_positive_end - current_region[0] + 1 >= min_length:
                            positive_regions.append([current_region[0], last_positive_end])
                        # Transition to negative region search
                        current_region = None
                        searching_for_positive = False  # Now look for a negative region
        else:
            if score <= 0:
                if current_region is None:
                    current_region = [i, i]  # Start a new negative region
                else:
                    current_region[1] = i  # Extend the current negative region
            else:
                # Start a new positive region
                current_region = [i, i]
                searching_for_positive = True  # Switch back to looking for a positive region
                neg_count = 0

    # Handle the last region if it didn't end with a negative peak
    if current_region:
        if searching_for_positive and current_region[1] - current_region[0] + 1 >= min_length:
            positive_regions.append(current_region)

    negative_peaks = []
    negative_peak_scores = []
    for i in range(1, len(positive_regions)):
        prev_end = positive_regions[i-1][1]
        next_start = positive_regions[i][0]
        inter_region_scores = scores[prev_end + 1:next_start]
        if len(inter_region_scores) > 0:
            most_negative_index = np.argmin(inter_region_scores) + prev_end + 1
            most_negative_score = scores[most_negative_index]
            negative_peaks.append(most_negative_index)
            negative_peak_scores.append(most_negative_score)

    positive_peaks = []
    positive_peak_scores = []
    for region in positive_regions:
        region_scores = scores[region[0]:region[1] + 1]
        peak_index = np.argmax(region_scores) + region[0]
        positive_peaks.append(peak_index)
        positive_peak_scores.append(scores[peak_index])

    # Create adjusted peaks and regions
    positive_peak_regions = [(region[0] + original_start, region[1] + original_start) for region in positive_regions]
    adjusted_positive_peaks = [(region[0] + region[1]) // 2 + original_start for region in positive_regions]

    positive_peaks = [p + original_start for p in positive_peaks]
    negative_peaks = [p + original_start for p in negative_peaks]

    return (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (adjusted_positive_peaks, positive_peak_regions),
    )

def write_bedgraph(scores, contigs, out_prefix, first_region=False):
    # Set the mode for writing files: 'w' for the first region (overwrite), 'a' for subsequent regions (append)
    mode = 'w' if first_region else 'a'

    # Unpack the original start and end from the contigs tuple
    original_start, original_end = contigs[0]

    # Initialize the bedGraph filename for combined scores
    filename = f"{out_prefix}_combined_scores.bedGraph"
    
    # Prepare to write data
    with open(filename, mode) as f:
        # Use the first score type to get the position data, assuming all score types have the same positions
        first_score_type = list(scores.keys())[0]
        contig_data = scores[first_score_type]

        for contig, start, score_array in contig_data:
            # Add "chr" prefix if not present
            if not contig.startswith("chr"):
                contig = f"chr{contig}"

            for i in range(len(score_array)):
                position = start + i
                # Ensure the position falls within the original range, excluding the overhangs
                if original_start <= position < original_end:
                    # Start by writing the contig, start, and end columns
                    line = [contig, str(position), str(position + 1)]
                    
                    # Add the scores for all score types at this position
                    for stype in scores.keys():
                        # Get the score for the current position from each score type
                        stype_score_array = scores[stype][0][2]
                        if i < len(stype_score_array):
                            if stype_score_array[i] == int(stype_score_array[i]):
                                line.append(f"{int(stype_score_array[i])}")
                            else:
                                line.append(f"{stype_score_array[i]:.2f}")
                        else:
                            line.append("0")  # Fill with zeros if no score

                    # Write the final line to the file, including all score types
                    f.write("\t".join(line) + "\n")

def write_nucleosome_peaks(peaks, contigs, out_prefix, first_region=False, flip_scores=False):
    # Set the mode for writing files: 'w' for the first region (overwrite), 'a' for subsequent regions (append)
    mode = 'w' if first_region else 'a'

    # Unpack the original start and end from the contigs tuple
    original_start, original_end = contigs[0]

    # Write nucleosome regions directly to a file
    nucleosome_filename = f"{out_prefix}.bed"
    with open(nucleosome_filename, mode) as f:
        for (contig, original_start), peak_data in peaks.items():
            # Add "chr" prefix if not present
            if not contig.startswith("chr"):
                contig = f"chr{contig}"

            num_positive_peaks = len(peak_data['adjusted_peaks'])
            num_negative_peaks = len(peak_data['negative_peaks'])

            for i in range(num_positive_peaks):
                adjusted_peak = peak_data['adjusted_peaks'][i]
                nucleosome_region_start = peak_data['nucleosome_regions'][i][0]
                nucleosome_region_end = peak_data['nucleosome_regions'][i][1]
                peak = peak_data['positive_peaks'][i]

                # Ensure the peak falls within the original range, excluding the overhangs
                if not (original_start <= peak < original_end):
                    continue

                # Find the closest upstream and downstream negative peaks
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

                # Safely retrieve the scores for these peaks
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

                # Flip scores if requested
                if flip_scores:
                    upstream_score *= -1
                    downstream_score *= -1
                    peak_score = -1 * peak_data['positive_peak_scores'][i]
                else:
                    peak_score = peak_data['positive_peak_scores'][i]

                # Use the strongest negative peak (most negative score) to calculate prominence
                average_flanking_score = np.mean([upstream_score, downstream_score])
                prominence_score = peak_score - average_flanking_score

                # Retrieve max coverage and position for the current peak
                max_coverage = peak_data['max_coverages'][i]
                max_position = peak_data['max_positions'][i]

                # Write to the output file
                f.write(f"{contig}\t{nucleosome_region_start}\t{nucleosome_region_end}\t"
                        f"{prominence_score:.2f}\t{adjusted_peak}\t"
                        f"{upstream_score:.2f}\t{upstream_negative_peak}\t"
                        f"{downstream_score:.2f}\t{downstream_negative_peak}\t"
                        f"{peak_score:.2f}\t{peak}\t"
                        f"{max_coverage}\t{max_position}\n")

def split_into_regions(contig, start, end, max_length=100000, overlap=1000):
    regions = []
    current_start = start
    while current_start < end:
        adjusted_start = max(current_start - overlap, start - overlap, 0)
        adjusted_end = min(current_start + max_length + overlap, end + overlap)
        original_start = current_start
        original_end = min(current_start + max_length, end)
        regions.append((contig, adjusted_start, adjusted_end, original_start, original_end))
        current_start = original_end  # Move start forward considering the original end, not the overlap
    return regions

def call_and_write_peaks(
    scores,
    coverage_scores,
    adjusted_start,
    original_start,
    original_end,
    contig,
    out_prefix,
    first_region,
    peak_type_label,
    flip_scores,
):
    """
    Calls peaks on the given scores, calculates coverage-based features, and writes results to BED.
    """
    # Call peaks using the same detection logic
    positive_peaks, negative_peaks, adjusted_peaks = find_peaks_and_regions(scores, adjusted_start, 50, 5)

    # Compute peak coverage features
    max_coverages = []
    max_positions = []
    for start, end in adjusted_peaks[1]:
        region_start = start - original_start
        region_end = end - original_start
        region_coverage = coverage_scores[region_start:region_end]

        if region_coverage.size > 0:
            max_coverages.append(np.max(region_coverage))
            max_positions.append(np.argmax(region_coverage) + region_start + original_start)
        else:
            max_coverages.append(0)
            max_positions.append(0)

    # Package peak data
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

    # Write peaks to BED file
    write_nucleosome_peaks(
        peaks,
        [(original_start, original_end)],
        out_prefix + peak_type_label,
        first_region,
        flip_scores,
    )

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Score fragmentomics data.")
    parser.add_argument('-b', '--bamfiles', nargs='+', help='BAM file(s) to process')
    parser.add_argument('-o', '--out_prefix', help='prefix to apply to output files (default: generated based on BAM file names and contigs)')
    parser.add_argument('-c', '--contigs', nargs='+', help='limit calculation to specified contig(s) and optional range (e.g. "2:100000-200000")')
    parser.add_argument('--mode-length', type=int, default=166, help='Mode fragment length (default: 166)')
    parser.add_argument('--frag-lower', type=int, default=126, help='Lower limit for fragment size (default: 126)')
    parser.add_argument('--frag-upper', type=int, default=186, help='Upper limit for fragment size (default: 186)')
    parser.add_argument('--max-duplicates', type=int, default=50, help='Maximum allowed duplicates per 100k fragments within 100k bases (default: 50)')
    parser.add_argument('--no_plot', action='store_true', help='Disable plotting of nps_smooth scores')
    parser.add_argument('--subsample', type=float, default=None, help='Subsampling proportion (e.g., 0.5 to subsample 50% of the reads)')

    args = parser.parse_args()

    # Generate default out_prefix if not provided
    if not args.out_prefix:
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in args.bamfiles]
        args.out_prefix = f"{'_'.join(bam_basenames)}"
        
        # Check if contigs are specified and if there's exactly one contig
        if args.contigs and len(args.contigs) == 1:
            single_contig = args.contigs[0]
            # Add the single contig to the output prefix
            args.out_prefix = f"{args.out_prefix}_{single_contig}"

    args.out_prefix = f"{args.out_prefix}_mode{str(args.mode_length)}_lower{str(args.frag_lower)}_upper{str(args.frag_upper)}"

    # If --bamfiles is specified, follow the original BAM workflow
    if args.bamfiles:
        # Open BAM files
        bamfiles = []
        for bamfile_path in args.bamfiles:
            try:
                bamfile = pysam.AlignmentFile(bamfile_path, "rb")
                bamfiles.append(bamfile)
            except FileNotFoundError as e:
                parser.error(f"Unable to open bamfile {bamfile_path} (file not found)")
                return -1
            except Exception as e:
                parser.error(f"Unable to open bamfile {bamfile_path}: {str(e)}")
                return -1

        # Split contigs based on specified ranges or all contigs
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

        # Precompute distributions for fragment lengths in the range
        nps_frag_range = range(args.frag_lower, args.frag_upper + 1)
        distributions = precompute_distributions(nps_frag_range, args.mode_length)

        first_region = True
        for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc='Scoring contigs'):
            # Score the contig to get nps and coverage arrays
            scores, nps_frag_range = score_contig(
                bamfiles, contig, adjusted_start, adjusted_end, args.mode_length, nps_frag_range, args.max_duplicates, distributions, args.subsample
            )
            
            # Extract nps smoothed scores and coverage
            nps_smoothed_scores = scores['nps_smoothed'][0][2]
            coverage_scores = scores['coverage'][0][2]

            call_and_write_peaks(
                scores=nps_smoothed_scores,
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

            flipped_scores = -1 * nps_smoothed_scores

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

            # Write the scores to the combined bedGraph
            write_bedgraph(scores, [(original_start, original_end)], args.out_prefix, first_region)

            first_region = False  # Set to False after processing the first region

        return 0

if __name__ == '__main__':
    sys.exit(main())
