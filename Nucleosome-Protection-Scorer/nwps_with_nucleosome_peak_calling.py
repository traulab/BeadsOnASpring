import sys
from tqdm import tqdm
import argparse
import pysam
import numpy as np
import os
import gc
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
import pyBigWig

def generate_paired_reads(bamfile, contig=None, start=None, end=None):
    """Generator taking paired-end data from a sam/bam file and yielding both pairs as (fwd, rev)"""
    _unpaired_reads = dict()
    for read in bamfile.fetch(contig, start, end, multiple_iterators=True):
        # Skip reads that aren't part of a fully-mapped pair
        if not read.is_proper_pair:
            continue
        
        name = read.query_name
        if name not in _unpaired_reads:
            _unpaired_reads[name] = read
            continue
        
        # Have found mate, yield both
        mate = _unpaired_reads[name]
        del _unpaired_reads[name]
        
        # Check if both pairs are mapped to the same strand
        if read.is_reverse == mate.is_reverse:
            continue
        
        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read

def generate_fragment_ranges(bamfile, contig, start, end):
    """Generate ranges covered by sequenced DNA fragments"""
    for r_fwd, r_rev in generate_paired_reads(bamfile, contig, start, end):
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

def create_distribution(fragment_length, nwps_frag_range, mode_DNA_length):
    # Ensure the fragment length falls within the nwps_frag_range
    if fragment_length not in nwps_frag_range:
        return None  # Return None for fragment lengths outside the scoring range

    if fragment_length < mode_DNA_length:
        total_length = mode_DNA_length + (mode_DNA_length - fragment_length)
    else:
        total_length = fragment_length

    # Calculate dynamic midpoint based on mode_DNA_length
    midpoint = (mode_DNA_length - 1) // 2  # equivalent to 82 for mode_DNA_length = 166
    second_half_start = midpoint + 1  # equivalent to 83 for mode_DNA_length = 166

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

    # Shifting the distribution to between 1 and -1
    combined_scores = [x - 1 for x in combined_scores]

    # Center the distribution around y-axis
    max_val = max(combined_scores)
    min_val = min(combined_scores)
    midpoint_val = (max_val + min_val) / 2
    centered_scores = [x - midpoint_val for x in combined_scores]

    return centered_scores


def score_contig(bamfiles, contig, start, end, mode_DNA_length, nwps_frag_range):
    """
    Applies a range of scoring algorithms to a given reference contig, returning the results.
    Also returns fragment start positions and lengths for further plotting.
    """

    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=int)
    nwps = np.zeros(ref_len, dtype=float)
    fragment_starts = []
    fragment_lengths = []

    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(bamfile, contig, start, end):

            frag_length = frag_end - frag_start
            fragment_starts.append(frag_start)  # Collect the fragment start position
            fragment_lengths.append(frag_length)  # Collect the fragment length

            # Calculate nwps using combined linear distributions
            fragment_scores = create_distribution(frag_length, nwps_frag_range, mode_DNA_length)

            if frag_length in nwps_frag_range:
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
                    nwps[start_pos:start_pos + len(fragment_scores)] += fragment_scores[:end_pos - start_pos]

            # Calculate coverage
            if frag_start >= start and frag_end <= end:
                if frag_length in nwps_frag_range:
                    coverage[frag_start - start:frag_end - start] += 1

    # Calculate average coverage and normalize other scores by this average
    avg_coverage = np.mean(coverage)
    if avg_coverage > 0:
        nwps_cov_adjusted = nwps / avg_coverage

    # Apply Savitzky-Golay filter
    window_size = 21  # should be odd, larger than polyorder, and appropriate to your data size
    polyorder = 2
    nwps_smoothed = savgol_filter(nwps, window_size, polyorder)

    # Return scores as tuples in the format (contig, start, score_array)
    scores = {
        'coverage': [(contig, start, coverage)],
        'nwps_smoothed': [(contig, start, nwps_smoothed)],
        'nwps': [(contig, start, nwps)],
    }

    return scores, fragment_starts, fragment_lengths, nwps_frag_range

def find_peaks_and_regions(scores, original_start, original_end, min_length=42, max_neg_run=42):
    positive_regions = []
    current_region = None
    stored_positive_start = None
    searching_for_positive = True  # Start by looking for positive regions

    for i in range(len(scores)):
        score = scores[i]

        if searching_for_positive:
            if score > 0:
                if current_region is None:
                    stored_positive_start = i  # Track where the score first becomes positive
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

    # Calculate the negative peaks between nucleosome regions
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

    positive_peak_regions = [(region[0] + original_start, region[1] + original_start) for region in positive_regions]
    adjusted_positive_peaks = [(region[0] + region[1]) // 2 + original_start for region in positive_regions]

    positive_peaks = [p + original_start for p in positive_peaks]
    negative_peaks = [p + original_start for p in negative_peaks]

    return (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (adjusted_positive_peaks, positive_peak_regions),
    )

def write_bedgraph(scores, peaks, contigs, out_prefix, nwps_frag_range, first_region=False, fragment_starts=None, fragment_lengths=None):
    mode = 'w' if first_region else 'a'  # Overwrite for the first region, append for subsequent regions

    # Unpack the original start and end from the contigs tuple
    contig_name, original_start, original_end = contigs[0]

    # Add "chr" prefix if not present
    # if not contig_name.startswith("chr"):
    #     contig_name = f"chr{contig_name}"

    # Initialize the bedGraph filename for combined scores
    filename = f"{out_prefix}.combined_scores.bedGraph"
    
    # Open the file and write the headers if it's the first region
    if first_region:
        with open(filename, 'w') as f:
            # headers = ['chrom', 'start', 'end'] + list(scores.keys())  # Column names for contig, start, end, and all score types
            f.write("")

    # Prepare to write data
    with open(filename, 'a') as f:
        # Use the first score type to get the position data, assuming all score types have the same positions
        first_score_type = list(scores.keys())[0]
        contig_data = scores[first_score_type]

        for contig, start, score_array in contig_data:
            # Add "chr" prefix if not present
            # if not contig.startswith("chr"):
            #     contig = f"chr{contig}"

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

def write_bedgraph_and_peaks(scores, peaks, contigs, out_prefix, nwps_frag_range, first_region=False, fragment_starts=None, fragment_lengths=None):
    # Set the mode for writing files: 'w' for the first region (overwrite), 'a' for subsequent regions (append)
    mode = 'w' if first_region else 'a'

    # Unpack the original start and end from the contigs tuple
    contig_name, original_start, original_end = contigs[0]

    # Add "chr" prefix if not present
    # if not contig_name.startswith("chr"):
    #     contig_name = f"chr{contig_name}"

    # Write the scores to the combined bedGraph
    write_bedgraph(scores, peaks, contigs, out_prefix, nwps_frag_range, first_region, fragment_starts, fragment_lengths)

    # Write nucleosome regions directly to a file
    nucleosome_filename = f"{out_prefix}_nucleosome_regions.bed"
    with open(nucleosome_filename, mode) as f:
        for (contig, original_start), peak_data in peaks.items():
            # Add "chr" prefix if not present
            # if not contig.startswith("chr"):
            #     contig = f"chr{contig}"

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
                    upstream_score = peak_data['positive_peak_scores'][i]  # Use the positive peak score if no upstream negative peak is found

                if downstream_index is not None:
                    downstream_negative_peak = peak_data['negative_peaks'][downstream_index]
                    downstream_score = peak_data['negative_peak_scores'][downstream_index]
                else:
                    downstream_negative_peak = peak
                    downstream_score = peak_data['positive_peak_scores'][i]  # Use the positive peak score if no downstream negative peak is found

                # Use the strongest negative peak (most negative score) to calculate prominence
                average_flanking_score = np.mean([upstream_score, downstream_score])
                prominence_score = peak_data['positive_peak_scores'][i] - average_flanking_score

                # Write to the output file in the required format
                f.write(f"{contig}\t{nucleosome_region_start}\t{nucleosome_region_end}\t"
                        f"{prominence_score:.2f}\t{adjusted_peak}\t"
                        f"{upstream_score:.2f}\t{upstream_negative_peak}\t"
                        f"{downstream_score:.2f}\t{downstream_negative_peak}\t"
                        f"{peak_data['positive_peak_scores'][i]:.2f}\t{peak}\n")

def plot_nwps_scores_directly(start, score_array, peak_data, out_prefix, contig_name, original_start, original_end, nwps_frag_range, score_type, fragment_starts=None, fragment_lengths=None):
    """
    Plots the nwps_smoothed_scores with positive and negative peaks directly from the input data,
    ensuring the plot matches exactly what is written to the bedGraph.
    Also includes plotting individual fragment distributions with vertical stacking, sorted by fragment start positions.
    """
    num_scores = len(score_array)
    plot_width = max(11.25, num_scores / 500)  # Calculate plot width based on the number of scores

    # Adjust height to accommodate fragment plots if fragment_starts is not None
    plot_height = 5 + len(fragment_starts) * 0.3 if fragment_starts else 5
    plot_height = 11.25

    x = np.arange(start, start + num_scores)
    
    # Create the plot with the calculated dimensions
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))
    ax.plot(x, score_array, label='NWPS Smoothed Scores', zorder=1)
    
    # Plot positive peaks as green dots
    # ax.scatter(peak_data['positive_peaks'], peak_data['positive_peak_scores'],
    #            color='green', label='Positive Peak', zorder=5)

    # Plot negative peaks as red dots
    # ax.scatter(peak_data['negative_peaks'], peak_data['negative_peak_scores'],
    #            color='red', label='Negative Peak', zorder=5)
    
    # Add a horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Set the x-limits to the original contig region
    ax.set_xlim(original_start, original_end)

    # Add labels and title
    ax.set_xlabel('Position')
    ax.set_ylabel('NWPS Score')
    ax.set_title(f'NWPS Scores for {contig_name}:{original_start}-{original_end}')
    ax.legend()

    all_y_values = list(score_array) + peak_data['positive_peak_scores'] + peak_data['negative_peak_scores']

    plt.tight_layout()

    plot_filename = f"{out_prefix}.{contig_name}_{original_start}-{original_end}_{score_type}"
    plt.savefig(f"{plot_filename}.svg")
    plt.close()

    print(f"Plot saved to {plot_filename}")

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

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Score fragmentomics data.",
        epilog="""By default, this script will process all contigs in
            the supplied BAM file. To limit processing to specific
            contigs, use the --contig argument.

            Note that the supplied file must be indexed (i.e.
            must have an associated .bai file).
            """
    )
    parser.add_argument('-b', '--bamfiles', nargs='+', required=True, help='BAM file(s) to process')
    parser.add_argument('-o', '--out_prefix', help='prefix to apply to output files (default: generated based on BAM file names and contigs)')
    parser.add_argument('-c', '--contigs', nargs='+', help='limit calculation to specified contig(s) and optional range (e.g. "2:100000-200000")')
    parser.add_argument('--mode-length', type=int, default=166, help='Mode fragment length (default: 166)')
    parser.add_argument('--frag-lower', type=int, default=126, help='Lower limit for fragment size (default: 126)')
    parser.add_argument('--frag-upper', type=int, default=186, help='Upper limit for fragment size (default: 186)')
    parser.add_argument('--no_plot', action='store_true', help='Disable plotting of nwps_smooth scores')
    args = parser.parse_args()

    # Generate default out_prefix if not provided
    if not args.out_prefix:
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in args.bamfiles]
        bam_part = "_".join(bam_basenames)
        if args.contigs:
            contig_part = "_".join(["{}-{}-{}".format(c.split(':')[0], *c.split(':')[1].split('-')) if ':' in c else c for c in args.contigs])
        else:
            contig_part = "all_contigs"
        args.out_prefix = f"{bam_part}_{contig_part}_mode{args.mode_length}_lower{args.frag_lower}_upper{args.frag_upper}"

    # Print the out_prefix for debugging
    print(f"Using output prefix: {args.out_prefix}")

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

    # Generate list of ranges to calculate over
    contigs = []
    if args.contigs is not None:
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
        # No ranges specified, default to all contigs
        for contig in bamfiles[0].references:
            start, end = 0, bamfiles[0].get_reference_length(contig)
            contigs.extend(split_into_regions(contig, start, end))

    # Print debug information about the ranges
    print("Contigs and ranges to process:")
    for contig, adjusted_start, adjusted_end, original_start, original_end in contigs:
        print(f"Contig: {contig}, Start: {original_start}, End: {original_end}")

    nwps_frag_range = range(args.frag_lower, args.frag_upper + 1)

    # Score the contigs and accumulate peaks
    first_region = True
    for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc='Scoring contigs'):
        scores, fragment_starts, fragment_lengths, nwps_frag_range = score_contig(
            bamfiles, contig, adjusted_start, adjusted_end, args.mode_length, nwps_frag_range
        )
        
        nwps_smoothed_scores = scores['nwps_smoothed'][0][2]  # Extract the score array

        # Get the peaks and positive peak regions while excluding overhang regions
        positive_peaks, negative_peaks, adjusted_peaks = find_peaks_and_regions(nwps_smoothed_scores, adjusted_start, adjusted_end, args.mode_length//4, args.mode_length//16)
        
        # Structure peaks data as a dictionary
        peaks = {
            (contig, original_start): {
                'positive_peaks': positive_peaks[0],
                'positive_peak_scores': positive_peaks[1],
                'negative_peaks': negative_peaks[0],
                'negative_peak_scores': negative_peaks[1],
                'adjusted_peaks': adjusted_peaks[0],
                'nucleosome_regions': adjusted_peaks[1],
            }
        }

        # Write the individual outputs (coverage, nwps_smoothed, peaks, positive peak regions)
        write_bedgraph_and_peaks(
            scores, 
            peaks, 
            [(contig, original_start, original_end)], 
            args.out_prefix, 
            nwps_frag_range, 
            first_region, 
            fragment_starts, 
            fragment_lengths, 
        )
                   
        first_region = False  # Set to False after processing the first region

        # Cleanup to free memory
        del peaks
        del scores
        del fragment_starts
        del fragment_lengths
        gc.collect()  # Manually trigger garbage collection
    
    return 0

if __name__ == '__main__':
    sys.exit(main())