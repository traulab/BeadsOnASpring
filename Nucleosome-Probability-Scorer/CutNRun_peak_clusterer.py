import pandas as pd
import numpy as np
import argparse
import os
from tqdm import tqdm

# Function to load BED-style files into a DataFrame
def load_bed_file(filename):
    print(f"---- DEBUG: Loading file {filename} ----")
    columns = ['contig', 'nucleosome_region_start', 'nucleosome_region_end', 'prominence_score',
               'center_position', 'upstream_negative_score', 'upstream_negative_peak_pos',
               'downstream_negative_score', 'downstream_negative_peak_pos', 'positive_peak_score', 'positive_peak_pos']
    df = pd.read_csv(filename, sep="\t", names=columns)
    print(f"---- DEBUG: File {filename} loaded ----")
    return df

# Binary search to find the first relevant control peak that might overlap
def binary_search_control(control_peaks, treatment_start):
    lo, hi = 0, len(control_peaks) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if control_peaks.iloc[mid]['nucleosome_region_end'] < treatment_start:
            lo = mid + 1
        else:
            hi = mid - 1
    return lo - 20 # This will point to the first potentially overlapping control peak

# Binary search to find the first relevant control cluster that might overlap
def binary_search_cluster(clusters, treatment_start):
    lo, hi = 0, len(clusters) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        # Access the first peak's 'nucleosome_region_end' in the cluster
        cluster_mid_peak = clusters[mid][0][0][1]  # clusters[mid] is the whole cluster, [0][0][1] accesses the first peak's pandas Series
        if cluster_mid_peak['nucleosome_region_end'] < treatment_start:
            lo = mid + 1
        else:
            hi = mid - 1
    return lo - 20  # This will point to the first potentially overlapping control cluster

def classify_peaks(df, threshold, control_df=None, control_threshold=None, control_clusters=None):
    # Classify treatment peaks based on threshold
    df['above_threshold'] = df['prominence_score'] >= threshold

    if control_df is not None and control_clusters is not None:
        # Filter control peaks to only include above-threshold control peaks
        control_df['above_threshold'] = control_df['prominence_score'] >= control_threshold
        above_threshold_control_df = control_df[control_df['above_threshold']]

        # Sort both treatment and control DataFrames for efficient searching
        df = df.sort_values(by=['contig', 'nucleosome_region_start'])
        above_threshold_control_df = above_threshold_control_df.sort_values(by=['contig', 'nucleosome_region_start'])

        # Process each contig in treatment data
        for contig in df['contig'].unique():
            print(f"Processing contig: {contig}")
            treatment_peaks = df[df['contig'] == contig]
            control_peaks = above_threshold_control_df[above_threshold_control_df['contig'] == contig]

            # Filter relevant clusters based on the current contig
            relevant_clusters = []
            for cluster in control_clusters:
                # Each cluster contains a list of tuples in the first element
                cluster_peaks = cluster[0]  # This is the list of ('A' or 'B', pandas Series) tuples
                
                # Access the first peak's Series
                first_peak_in_cluster = cluster_peaks[0][1]  # Access the pandas Series from the first tuple
                
                # Ensure it's a pandas Series and then check the contig
                if isinstance(first_peak_in_cluster, pd.Series):
                    if first_peak_in_cluster['contig'] == contig:  # Correctly access the contig
                        relevant_clusters.append(cluster)
                else:
                    # Print the structure in case of unexpected type
                    print(f"DEBUG: Unexpected structure - {type(first_peak_in_cluster)}")

            print(f"Relevant clusters for contig {contig}: {len(relevant_clusters)}")

            control_len = len(control_peaks)
            cluster_len = len(relevant_clusters)

            # Continue with the binary search and filtering as before
            for treatment_idx, treatment_peak in tqdm(treatment_peaks.iterrows(), total=len(treatment_peaks), desc=f"Processing contig {contig}"):
                if not treatment_peak['above_threshold']:
                    continue  # Only process above-threshold treatment peaks

                # Binary search to find the first relevant control peak
                control_idx = binary_search_control(control_peaks, treatment_peak['nucleosome_region_start'])

                # Check control peaks for overlap
                for i in range(control_idx, control_len):
                    control_peak = control_peaks.iloc[i]

                    # If control peak starts after the treatment peak ends, stop
                    if control_peak['nucleosome_region_start'] > treatment_peak['nucleosome_region_end']:
                        break

                    # Check if there's a true overlap between the current treatment and control peak
                    if (treatment_peak['nucleosome_region_start'] <= control_peak['nucleosome_region_end'] and
                        treatment_peak['nucleosome_region_start'] >= control_peak['nucleosome_region_start']):
                        df.at[treatment_idx, 'above_threshold'] = False
                        break
                    elif (treatment_peak['nucleosome_region_end'] <= control_peak['nucleosome_region_end'] and
                        treatment_peak['nucleosome_region_end'] >= control_peak['nucleosome_region_start']):
                        df.at[treatment_idx, 'above_threshold'] = False
                        break

                # Binary search for relevant control cluster
                cluster_idx = binary_search_cluster(relevant_clusters, treatment_peak['nucleosome_region_start'])

                if cluster_len == 0:
                    continue

                # Check if treatment peaks fall between the start and end of the control cluster
                for j in range(cluster_idx, cluster_len):
                    control_cluster = relevant_clusters[j]

                    # Find the overall start and end of the cluster (from the first peak's start to the last peak's end)
                    control_start = control_cluster[0][0][1]['nucleosome_region_start']  # Accessing the first peak
                    control_end = control_cluster[0][-1][1]['nucleosome_region_end']  # Accessing the last peak

                    if control_start > treatment_peak['nucleosome_region_end']:
                        break

                    # If the current cluster starts after the treatment peak ends, stop
                    if control_start > treatment_peak['nucleosome_region_end']:
                        break

                    # Filter out treatment peaks that overlap (even partially) with the control cluster
                    if (treatment_peak['nucleosome_region_start'] >= control_start and 
                        treatment_peak['nucleosome_region_start'] <= control_end):
                        df.at[treatment_idx, 'above_threshold'] = False
                        break  # No need to check further clusters if filtered
                    elif (treatment_peak['nucleosome_region_end'] <= control_end and 
                        treatment_peak['nucleosome_region_end'] >= control_start):
                        df.at[treatment_idx, 'above_threshold'] = False
                        break  # No need to check further clusters if filtered

    return df

# Function to find clusters based on the sliding window approach using integer scoring (+10, -2)
def find_clusters(df, consecutive):
    print(f"---- DEBUG: Finding clusters ----")
    clusters = []
    current_cluster = []
    current_score = 0  # Use integer-based scoring: +10 for above-threshold, -2 for below-threshold
    peak_count_since_last_above_threshold = 0
    valid_cluster = False

    for i, row in df.iterrows():
        if row['above_threshold']:
            # Start or continue a cluster with an above-threshold peak
            current_cluster.append(('A', row))  # A = above-threshold
            current_score += 10  # Add +10 for each above-threshold peak
            valid_cluster = True  # A valid cluster starts with an above-threshold peak
            peak_count_since_last_above_threshold = 0  # Reset the count since the last above-threshold peak
        else:
            # Handle below-threshold peaks, but only if the cluster has already started
            if valid_cluster:
                current_cluster.append(('B', row))  # B = below-threshold
                peak_count_since_last_above_threshold += 1

        # Close the cluster as soon as 5 consecutive below-threshold peaks are reached
        if peak_count_since_last_above_threshold == consecutive:
            # Remove the last 5 below-threshold peaks and restore the score
            for _ in range(consecutive):
                if current_cluster and current_cluster[-1][0] == 'B':
                    current_cluster.pop()

            # Finalize and add the current cluster if it has more than 1 peak
            if valid_cluster and len(current_cluster) > 1:
                max_prominence = max([peak[1]['prominence_score'] for peak in current_cluster])

                final_score = int(current_score / 10)

                if final_score >= 1:
                    clusters.append((current_cluster, final_score, max_prominence))

            # Reset cluster tracking as soon as 5 below-threshold peaks are found
            current_cluster = []
            current_score = 0
            valid_cluster = False
            peak_count_since_last_above_threshold = 0

    # Final cluster check for any open cluster at the end
    if len(current_cluster) > 1 and valid_cluster and current_score >= 10:
        final_score = int(current_score / 10)
        max_prominence = max([peak[1]['prominence_score'] for peak in current_cluster])
        clusters.append((current_cluster, final_score, max_prominence))

    print(f"---- DEBUG: Clustering completed ----")
    return clusters

# Updated function to write output files (including full details in the BED files for treatment and control clusters)
def write_output_files(clusters, control_clusters, treatment_df, control_df, output_prefix, multiplier):
    threshold_suffix = f"_{multiplier}SD_threshold"  # Create the suffix with the multiplier

    print(f"---- DEBUG: Writing output files ----")
    
    # Treatment cluster output (clusters with full details)
    with open(f"{output_prefix}_treatment_clusters{threshold_suffix}.txt", 'w') as f1:
        f1.write("chromosome\tstart\tend\tcluster_size\tpeak_count\tmax_prominence\tmost_prominent_peak_center\n")
        for cluster, score, max_prominence in clusters:
            start = cluster[0][1]['nucleosome_region_start']
            end = cluster[-1][1]['nucleosome_region_end']
            size = len(cluster)

            # Find the center position of the most prominent peak
            most_prominent_peak = max(cluster, key=lambda x: x[1]['prominence_score'])
            most_prominent_peak_center = most_prominent_peak[1]['center_position']

            f1.write(f"{cluster[0][1]['contig']}\t{start}\t{end}\t{size}\t{score}\t{max_prominence}\t{most_prominent_peak_center}\n")

    # Treatment cluster BED output (including all original columns)
    with open(f"{output_prefix}_treatment_cluster_peaks{threshold_suffix}.bed", 'w') as f1_bed:
        for cluster, _, _ in clusters:
            for peak in cluster:
                row = peak[1]
                f1_bed.write("\t".join(map(str, row)) + "\n")
    
    # Control cluster output (clusters with full details)
    with open(f"{output_prefix}_control_clusters{threshold_suffix}.txt", 'w') as f2:
        f2.write("chromosome\tstart\tend\tcluster_size\tpeak_count\tmax_prominence\tmost_prominent_peak_center\n")
        for cluster, score, max_prominence in control_clusters:
            start = cluster[0][1]['nucleosome_region_start']
            end = cluster[-1][1]['nucleosome_region_end']
            size = len(cluster)

            # Find the center position of the most prominent peak
            most_prominent_peak = max(cluster, key=lambda x: x[1]['prominence_score'])
            most_prominent_peak_center = most_prominent_peak[1]['center_position']

            f2.write(f"{cluster[0][1]['contig']}\t{start}\t{end}\t{size}\t{score}\t{max_prominence}\t{most_prominent_peak_center}\n")

    # Control cluster BED output (including all original columns)
    with open(f"{output_prefix}_control_cluster_peaks{threshold_suffix}.bed", 'w') as f2_bed:
        for cluster, _, _ in control_clusters:
            for peak in cluster:
                row = peak[1]
                f2_bed.write("\t".join(map(str, row)) + "\n")

    # Output filtered treatment peaks (with full details)
    treatment_filtered = treatment_df[treatment_df['above_threshold'] == True]
    treatment_filtered.to_csv(f"{output_prefix}_filtered_treatment_peaks{threshold_suffix}.bed", sep="\t", index=False, header=False)

    # Output filtered control peaks (with full details)
    control_filtered = control_df[control_df['above_threshold'] == True]
    control_filtered.to_csv(f"{output_prefix}_filtered_control_peaks{threshold_suffix}.bed", sep="\t", index=False, header=False)

    print(f"---- DEBUG: Output files written ----")

# Main function to handle separate treatment and control threshold calculations and clustering
def main():
    parser = argparse.ArgumentParser(description="Cluster analysis on BED-style peak files.")
    parser.add_argument("treatment_file", help="Input BED-style file for treatment.")
    parser.add_argument("control_file", help="Input BED-style file for negative control.")
    parser.add_argument("--multiplier", type=float, default=1, help="Multiplier for standard deviation threshold (default=1).")
    parser.add_argument("--output_prefix", type=str, default=None, help="Prefix for output files.")
    
    args = parser.parse_args()

    # Set output prefix based on treatment file if not provided
    if args.output_prefix is None:
        args.output_prefix = os.path.splitext(os.path.basename(args.treatment_file))[0]

    print(f"---- DEBUG: Starting script ----")

    # Load the input files
    treatment_df = load_bed_file(args.treatment_file)
    control_df = load_bed_file(args.control_file)

    # Calculate thresholds for treatment and control
    treatment_threshold = treatment_df['prominence_score'].mean() + args.multiplier * treatment_df['prominence_score'].std()
    control_threshold = control_df['prominence_score'].mean() + args.multiplier * control_df['prominence_score'].std()

    print(f"Treatment mean: {treatment_df['prominence_score'].mean()}, threshold: {treatment_threshold}")
    print(f"Control mean: {control_df['prominence_score'].mean()}, threshold: {control_threshold}")

    # Classify control peaks and find clusters without filtering itself
    control_classified = classify_peaks(control_df, control_threshold)  # No control_df and control_clusters arguments here
    control_clusters = find_clusters(control_classified, 3)

    # Classify treatment peaks, including filtering for overlaps with control clusters
    treatment_classified = classify_peaks(treatment_df, treatment_threshold, control_df, control_threshold, control_clusters)

    # Find treatment clusters
    treatment_clusters = find_clusters(treatment_classified, 5)

    # Write output files
    write_output_files(treatment_clusters, control_clusters, treatment_classified, control_df, args.output_prefix, args.multiplier)

    print(f"---- DEBUG: Script completed ----")

if __name__ == "__main__":
    main()