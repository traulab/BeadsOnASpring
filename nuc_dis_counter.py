
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import argparse
from collections import defaultdict
from intervaltree import IntervalTree
import numpy as np
from scipy.signal import savgol_filter
from scipy import stats
from scipy.stats import mode
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def smooth_and_find_mode(distance_dict, window_length=35, polyorder=3):
    """
    Applies Savitzky-Golay smoothing to the count values associated with each distance and finds the mode of the smoothed data.
    Also calculates the median and mean of the unsmoothed data.
    
    Args:
    - distance_dict: Dictionary with distances as keys and counts as values.
    - window_length: Window length for Savitzky-Golay smoothing.
    - polyorder: Polynomial order for Savitzky-Golay smoothing.
    
    Returns:
    - mode_dist: The mode of the smoothed data.
    - median_dist: The median of the unsmoothed data.
    - mean_dist: The mean of the unsmoothed data.
    - smoothed_dict: Dictionary with distances and smoothed counts.
    """
    # Extract distances and counts as arrays
    distances = sorted(distance_dict.keys())
    counts = [distance_dict[dist] for dist in distances]

    # Create an expanded array of distances based on counts for median and mean calculation
    expanded_distances = []
    for dist, count in distance_dict.items():
        expanded_distances.extend([dist] * count)

    # Calculate the median and mean of the unsmoothed data
    median_dist = np.median(expanded_distances) if expanded_distances else None
    mean_dist = np.mean(expanded_distances) if expanded_distances else None

    # Ensure we have enough data points for smoothing
    if len(distances) > window_length:
        # Apply Savitzky-Golay smoothing to the counts (not distances)
        try:
            smoothed_counts = savgol_filter(counts, window_length=window_length, polyorder=polyorder)
        except ValueError as e:
            print(f"Error during smoothing: {e}")
            smoothed_counts = counts  # If smoothing fails, return unsmoothed counts
    else:
        print(f"Warning: Not enough data points for smoothing (data length: {len(distances)}).")
        smoothed_counts = counts  # If insufficient data, don't smooth

    # Create a smoothed dictionary for counts
    smoothed_dict = defaultdict(int)
    for dist, smoothed_count in zip(distances, smoothed_counts):
        smoothed_dict[dist] = max(0, int(smoothed_count))  # Ensure non-negative counts

    # Calculate the mode of the smoothed counts
    mode_index = np.argmax(smoothed_counts)  # Find the index of the maximum smoothed count
    mode_value = distances[mode_index]  # Get the corresponding distance

    return mode_value, median_dist, mean_dist, smoothed_dict

def smooth_normalized_values(normalized_values, window_length=35, polyorder=3):
    """
    Smooths the normalized values using Savitzky-Golay filter.
    
    Args:
    - normalized_values: List of normalized values to smooth.
    - window_length: Window length for Savitzky-Golay smoothing.
    - polyorder: Polynomial order for Savitzky-Golay smoothing.
    
    Returns:
    - List of smoothed normalized values.
    """
    if len(normalized_values) > window_length:
        try:
            smoothed_normalized = savgol_filter(normalized_values, window_length=window_length, polyorder=polyorder)
        except ValueError as e:
            print(f"Error during smoothing normalized values: {e}")
            smoothed_normalized = normalized_values  # If smoothing fails, return unsmoothed values
    else:
        smoothed_normalized = normalized_values  # If insufficient data, don't smooth

    return smoothed_normalized


def process_bed_file(input_file, chromatin_state_file, output_file, combine_euchromatin=False):
    chrom_distances = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # {chromosome: {state: {distance: count}}}
    chrom_all_distances = defaultdict(lambda: defaultdict(int))  # {chromosome: {distance: count}}
    genome_distances = defaultdict(lambda: defaultdict(int))  # {state: {distance: count}}
    genome_all_distances = defaultdict(int)  # {distance: count}
    duplicate_positions = defaultdict(lambda: defaultdict(list))  # {chromosome: {position: [states]}}

    euchromatin_states = ("1_", "2_", "4_", "5_", "6_", "7_", "9_", "10_", "11_")

    # Read chromatin state BED file into interval trees for each chromosome
    chromatin_state_trees = defaultdict(IntervalTree)  # {chromosome: IntervalTree}
    with open(chromatin_state_file, 'r') as csfile:
        for line in csfile:
            cols = line.strip().split()
            chrom = f"{cols[0]}"
            start = int(cols[1])
            end = int(cols[2])
            state = cols[3]  # Chromatin state is in column 4
            chromatin_state_trees[chrom].addi(start, end, state)
            if any(state.startswith(euchromatin_state) for euchromatin_state in euchromatin_states):
                chromatin_state_trees[chrom].addi(start, end, "Euchromatin")

    # Track the last nucleosome position and chromatin state for each chromosome
    prev_position_by_chrom = defaultdict(lambda: None)  # Track previous position for each chromosome
    prev_state_by_chrom = defaultdict(lambda: None)  # Track previous chromatin state for each chromosome
    
    with open(input_file, 'r') as infile:
        
        for line in infile:
            cols = line.strip().split()
            #chrom = f"chr{cols[0]}"
            chrom = f"{cols[0]}"
            position = int(cols[4])  # Nucleosome position is in column 5
            
            # Find overlap with chromatin states using interval trees
            overlaps = chromatin_state_trees[chrom][position]
            if overlaps:
                chromatin_state_match = ','.join(sorted(set([ov.data for ov in overlaps])))
            else:
                chromatin_state_match = "None"
            
            if combine_euchromatin and any(chromatin_state_match.startswith(eu_state) for eu_state in euchromatin_states):
                chromatin_state_match = "Euchromatin"

            # Check for duplicate positions within the same chromosome
            if prev_position_by_chrom[chrom] is not None and prev_position_by_chrom[chrom] == position:
                duplicate_positions[chrom][position].append(chromatin_state_match)
            
            # Calculate distance only if both positions share the same chromatin state
            prev_position = prev_position_by_chrom.get(chrom)
            prev_state = prev_state_by_chrom.get(chrom)
            if prev_position is not None and prev_position != position and prev_state == chromatin_state_match:
                distance = abs(position - prev_position)  # Ensure positive distance
                
                if distance > 1000 or distance < 50:
                    prev_position_by_chrom[chrom] = position
                    prev_state_by_chrom[chrom] = chromatin_state_match
                    continue  # Don't update prev_position and prev_state here
                
                # Update chromosome distances by chromatin state
                chrom_distances[chrom][chromatin_state_match][distance] += 1
                
                # Update chromosome distances (state: "All")
                chrom_all_distances[chrom][distance] += 1
                
                # Update whole genome distances by chromatin state
                genome_distances[chromatin_state_match][distance] += 1
                
                # Update whole genome distances (state: "All")
                genome_all_distances[distance] += 1

            prev_position_by_chrom[chrom] = position
            prev_state_by_chrom[chrom] = chromatin_state_match

    # Main output file
    with open(output_file, 'w') as outfile:
        # Chromosome Distances (broken down by chromatin state)
        outfile.write("Chromosome Distances by Chromatin State:\n")
        for chrom in chrom_distances:
            outfile.write(f"\n{chrom}\n")  # Empty line before each chromosome
            for state in chrom_distances[chrom]:
                mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(chrom_distances[chrom][state])
                if mode_dist is None:
                    print(f"Skipping smoothing for {chrom}, state: {state} due to errors.")
                    continue
                
                # Calculate total counts for normalization
                total_count = sum(chrom_distances[chrom][state].values())
                
                # Prepare normalized and smoothed normalized counts
                normalized_values = []
                for dist in sorted(chrom_distances[chrom][state].keys()):
                    count = chrom_distances[chrom][state][dist]
                    normalized_count = round((count / total_count) * 100, 3) if total_count > 0 else 0
                    normalized_values.append(normalized_count)
                
                smoothed_normalized_values = smooth_normalized_values(normalized_values)

                # Output the mode, median, and smoothed normalized counts
                outfile.write(f"\nState: {state} (Mode Distance: {mode_dist}, Median Distance: {median_dist}, Mean Distance: {mean_dist})\n")
                for i, dist in enumerate(sorted(chrom_distances[chrom][state].keys())):
                    smoothed_count = smoothed_counts.get(dist, 0)
                    count = chrom_distances[chrom][state][dist]
                    normalized_count = normalized_values[i]
                    smoothed_normalized_count = smoothed_normalized_values[i]
                    outfile.write(f"{chrom}\t{state}\t{dist}\t{count}\t{smoothed_count}\t{normalized_count:.3f}\t{smoothed_normalized_count:.3f}\n")
            
        # Chromosome Distances (not broken down by chromatin state, state = "All")
        outfile.write("\nChromosome Distances (All States):\n")
        for chrom in chrom_all_distances:
            mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(chrom_all_distances[chrom])
            if mode_dist is None:
                print(f"Skipping smoothing for {chrom} (All States) due to errors.")
                continue
            
            # Calculate total counts for normalization
            total_count = sum(chrom_all_distances[chrom].values())

            # Prepare normalized and smoothed normalized counts
            normalized_values = []
            for dist in sorted(chrom_all_distances[chrom].keys()):
                count = chrom_all_distances[chrom][dist]
                normalized_count = round((count / total_count) * 100, 3) if total_count > 0 else 0
                normalized_values.append(normalized_count)
            
            smoothed_normalized_values = smooth_normalized_values(normalized_values)

            # Output the mode, median, and smoothed normalized counts
            outfile.write(f"\n{chrom} (Mode Distance: {mode_dist}, Median Distance: {median_dist}, Mean Distance: {mean_dist})\n")
            for i, dist in enumerate(sorted(chrom_all_distances[chrom].keys())):
                smoothed_count = smoothed_counts.get(dist, 0)
                count = chrom_all_distances[chrom][dist]
                normalized_count = normalized_values[i]
                smoothed_normalized_count = smoothed_normalized_values[i]
                outfile.write(f"{chrom}\tAll\t{dist}\t{count}\t{smoothed_count}\t{normalized_count:.3f}\t{smoothed_normalized_count:.3f}\n")

        # Whole-Genome Distances by Chromatin State
        outfile.write("\nWhole Genome Distances by Chromatin State:\n")
        for state in genome_distances:
            mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(genome_distances[state])
            if mode_dist is None:
                print(f"Skipping smoothing for whole genome, state: {state} due to errors.")
                continue
            
            # Calculate total counts for normalization
            total_count = sum(genome_distances[state].values())

            # Prepare normalized and smoothed normalized counts
            normalized_values = []
            for dist in sorted(genome_distances[state].keys()):
                count = genome_distances[state][dist]
                normalized_count = round((count / total_count) * 100, 3) if total_count > 0 else 0
                normalized_values.append(normalized_count)
            
            smoothed_normalized_values = smooth_normalized_values(normalized_values)

            # Output the mode, median, and smoothed normalized counts
            outfile.write(f"\nState: {state} (Mode Distance: {mode_dist}, Median Distance: {median_dist}, Mean Distance: {mean_dist})\n")
            for i, dist in enumerate(sorted(genome_distances[state].keys())):
                smoothed_count = smoothed_counts.get(dist, 0)
                count = genome_distances[state][dist]
                normalized_count = normalized_values[i]
                smoothed_normalized_count = smoothed_normalized_values[i]
                outfile.write(f"whole-genome\t{state}\t{dist}\t{count}\t{smoothed_count}\t{normalized_count:.3f}\t{smoothed_normalized_count:.3f}\n")

        # Whole-Genome Distances (state: "All")
        mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(genome_all_distances)
        if mode_dist is None:
            print(f"Skipping smoothing for whole genome (All States) due to errors.")
        else:
            # Calculate total counts for normalization
            total_count = sum(genome_all_distances.values())

            # Prepare normalized and smoothed normalized counts
            normalized_values = []
            for dist in sorted(genome_all_distances.keys()):
                count = genome_all_distances[dist]
                normalized_count = round((count / total_count) * 100, 3) if total_count > 0 else 0
                normalized_values.append(normalized_count)
            
            smoothed_normalized_values = smooth_normalized_values(normalized_values)

            # Output the mode, median, and smoothed normalized counts
            outfile.write(f"\nWhole Genome (All States) (Mode Distance: {mode_dist}, Median Distance: {median_dist}, Mean Distance: {mean_dist})\n")
            for i, dist in enumerate(sorted(genome_all_distances.keys())):
                smoothed_count = smoothed_counts.get(dist, 0)
                count = genome_all_distances[dist]
                normalized_count = normalized_values[i]
                smoothed_normalized_count = smoothed_normalized_values[i]
                outfile.write(f"whole-genome\tAll\t{dist}\t{count}\t{smoothed_count}\t{normalized_count:.3f}\t{smoothed_normalized_count:.3f}\n")

        # Output duplicate positions with their chromatin states
        outfile.write("\nDuplicate Positions with Chromatin States:\n")
        for chrom in duplicate_positions:
            for position in duplicate_positions[chrom]:
                states = ', '.join(duplicate_positions[chrom][position])
                outfile.write(f"{chrom}\t{position}\t{states}\n")

        # Prepare a dictionary to store distances per chromosome and chromatin state
        distance_data = defaultdict(lambda: defaultdict(list))  # {chromosome: {state: [distances]}}

        # Fill the dictionary with expanded distances
        for chrom in chrom_distances:
            for state in chrom_distances[chrom]:
                for dist, count in chrom_distances[chrom][state].items():
                    expanded_distances = [dist] * count  # Create a list of distances repeated 'count' times
                    distance_data[chrom][state].extend(expanded_distances)

        # Convert the dictionary to a DataFrame
        df_list = []
        for chrom, states in distance_data.items():
            max_len = max(len(distances) for distances in states.values())  # Find the maximum number of distances for any state
            row_data = {"Chromosome": chrom}  # Start with the chromosome column
            for state, distances in states.items():
                # Pad distances with None to ensure all columns have the same length
                row_data[state] = distances + [None] * (max_len - len(distances))
            df_list.append(pd.DataFrame(row_data))  # Append DataFrame for this chromosome

        # Concatenate all chromosomes into a single DataFrame
        df = pd.concat(df_list, ignore_index=True)

        # Write the DataFrame to a CSV file
        dis_output_file = output_file.replace('.txt', '_distances.csv')
        df.to_csv(dis_output_file, index=False)

    chromosome_sizes = {
        'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
        'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
        'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
        'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'whole-genome': sum([
            249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
            141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 
            81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560
        ])
    }

    # Sort the chromosomes by genomic size (largest to smallest)
    chromosome_order = sorted(chromosome_sizes, key=chromosome_sizes.get, reverse=True)

    # Read the CSV file with distances
    df = pd.read_csv(dis_output_file)

    # Select all autosomal chromosomes (without filtering them out yet)
    autosomal_df = df[df['Chromosome'].str.startswith('chr') & ~df['Chromosome'].isin(['chrX', 'chrY'])]

    # Append a new entry for "whole-genome" by combining all autosomal distances
    autosomal_df.loc[:, 'Chromosome'] = 'whole-genome'

    # Combine the dataframes
    combined_df = pd.concat([df, autosomal_df], ignore_index=True)

    chromatin_states_filtered = ['Euchromatin'] + [col for col in df.columns if col.startswith('13_')]

    # Only select columns that exist in the DataFrame
    available_columns = [col for col in chromatin_states_filtered if col in combined_df.columns]

    filtered_combined_df = combined_df[['Chromosome'] + available_columns]

    # Melt the dataframe to long format for easier boxplot creation
    melted_df = filtered_combined_df.melt(id_vars=['Chromosome'], var_name='Chromatin State', value_name='Distance')

    # Drop rows where Distance is NaN (because some states have fewer distances)
    melted_df = melted_df.dropna(subset=['Distance'])

    # Ensure chromosomes are ordered by their size for plotting
    melted_df['Chromosome'] = pd.Categorical(melted_df['Chromosome'], categories=chromosome_order, ordered=True)

    # Plot 1: Grouped by Chromosome
    plot_df = melted_df[melted_df['Chromosome'].isin(['chrX', 'chr7', 'whole-genome'])]

    plt.figure(figsize=(12, 6))
    sns.boxplot(x='Chromosome', y='Distance', hue='Chromatin State', data=melted_df, showfliers=False, width=0.6)

    plt.title('Comparison of Euchromatin and 13_ States Across Chromosomes (Ordered by Size)')
    plt.xlabel('Chromosome')
    plt.ylabel('Distance')
    plt.xticks(rotation=45)
    plt.legend(title='Chromatin State')
    plt.tight_layout()
    plt.savefig(output_file.replace('.txt', '_by_chromosome.svg'))

    # Plot 2: Grouped by Chromatin State
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='Chromatin State', y='Distance', hue='Chromosome', data=melted_df, showfliers=False, width=0.6)

    plt.title('Comparison of Chromosomes Grouped by Chromatin State (Ordered by Size)')
    plt.xlabel('Chromatin State')
    plt.ylabel('Distance')
    plt.xticks(rotation=45)
    plt.legend(title='Chromosome')
    plt.tight_layout()
    plt.savefig(output_file.replace('.txt', '_by_state.svg'))

    # Perform Mann-Whitney U tests for all chromosomes, comparing "Euchromatin" with "13_" states
    mannwhitney_results = []

    # Perform Mann-Whitney U tests on the full melted_df (no chromosome filtering)
    for chrom in melted_df['Chromosome'].unique():  # Perform for all chromosomes
        # Get distances for "Euchromatin" in the current chromosome
        euchromatin_distances = melted_df[(melted_df['Chromosome'] == chrom) & (melted_df['Chromatin State'] == 'Euchromatin')]['Distance']
        
        # Get distances for "13_" states in the current chromosome
        for state in [s for s in chromatin_states_filtered if s.startswith('13_')]:
            state_distances = melted_df[(melted_df['Chromosome'] == chrom) & (melted_df['Chromatin State'] == state)]['Distance']
            
            # Check if there are enough distances to compare
            if len(euchromatin_distances) > 1 and len(state_distances) > 1:
                # Perform a Mann-Whitney U test (non-parametric test)
                u_stat, p_value = stats.mannwhitneyu(euchromatin_distances, state_distances, alternative='two-sided')
                
                # Store the result in a list
                mannwhitney_results.append({
                    'Chromosome': chrom,
                    'State': state,
                    'U-statistic': u_stat,
                    'p-value': p_value
                })

    # Convert the Mann-Whitney U test results into a DataFrame for easy viewing
    mannwhitney_df = pd.DataFrame(mannwhitney_results)

    # Display the Mann-Whitney U test results
    print(mannwhitney_df)

    # Optionally, save the results to a CSV file
    mannwhitney_output_file = output_file.replace('.txt', '_mannwhitney_results_filtered.csv')
    mannwhitney_df.to_csv(mannwhitney_output_file, index=False)

    # Additional output file for statistics (mode, median, mean)
    stats_output_file = output_file.replace('.txt', '_stats.txt')
    with open(stats_output_file, 'w') as stats_outfile:
        stats_outfile.write("Chromosome\tState\tMode (Smoothed)\tMedian (Smoothed)\tMean (Smoothed)\tMode (Unsmoothed)\tMedian (Unsmoothed)\tMean (Unsmoothed)\n")
        
        # Write statistics for chromosomal data
        for chrom in chrom_distances:
            for state in chrom_distances[chrom]:
                mode_smooth, median_smooth, mean_smooth, _ = smooth_and_find_mode(chrom_distances[chrom][state])
                mode_unsmooth, median_unsmooth, mean_unsmooth, _ = smooth_and_find_mode(chrom_distances[chrom][state], window_length=1, polyorder=0)
                stats_outfile.write(f"{chrom}\t{state}\t{mode_smooth}\t{median_smooth}\t{mean_smooth}\t{mode_unsmooth}\t{median_unsmooth}\t{mean_unsmooth}\n")
        
        # Write statistics for genome-wide data
        stats_outfile.write("\nWhole Genome\n")
        for state in genome_distances:
            mode_smooth, median_smooth, mean_smooth, _ = smooth_and_find_mode(genome_distances[state])
            mode_unsmooth, median_unsmooth, mean_unsmooth, _ = smooth_and_find_mode(genome_distances[state], window_length=1, polyorder=0)
            stats_outfile.write(f"whole-genome\t{state}\t{mode_smooth}\t{median_smooth}\t{mean_smooth}\t{mode_unsmooth}\t{median_unsmooth}\t{mean_unsmooth}\n")

        # Write genome-wide statistics for "All"
        mode_smooth, median_smooth, mean_smooth, _ = smooth_and_find_mode(genome_all_distances)
        mode_unsmooth, median_unsmooth, mean_unsmooth, _ = smooth_and_find_mode(genome_all_distances, window_length=1, polyorder=0)
        stats_outfile.write(f"whole-genome\tAll\t{mode_smooth}\t{median_smooth}\t{mean_smooth}\t{mode_unsmooth}\t{median_unsmooth}\t{mean_unsmooth}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate nucleosome distances and counts by chromosome, chromatin state, and whole genome, and report duplicate positions.")
    parser.add_argument("input", help="Path to the input BED-like file.")
    parser.add_argument("chromatin_state", help="Path to the chromatin state BED file.")
    parser.add_argument("output", help="Path to the output file.")
    parser.add_argument("--combine-euchromatin", action="store_true", help="Combine chromatin states starting with 1_, 2_, 4_, 5_, 6_, 7_, 9_, 10_, and 11_ into 'Euchromatin'.")

    args = parser.parse_args()

    process_bed_file(args.input, args.chromatin_state, args.output, args.combine_euchromatin)
