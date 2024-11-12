import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, norm
from collections import defaultdict
from intervaltree import IntervalTree
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time

# Set matplotlib to use the Agg backend for environments without display
plt.switch_backend('Agg')

# Function to annotate chromatin states
def annotate_chromatin_states(nucleosome_data, chromatin_state_file):
    chromatin_state_trees = defaultdict(IntervalTree)

    # Start timing the reading of the chromatin state file
    start_time = time.time()

    # Read the chromatin state BED file and store the intervals in IntervalTree for each chromosome
    print(f"Reading chromatin state file: {chromatin_state_file}")
    with open(chromatin_state_file, 'r') as csfile:
        for i, line in enumerate(csfile):
            cols = line.strip().split()
            chrom = f"{cols[0]}"
            start = int(cols[1])
            end = int(cols[2])
            state = cols[3]
            chromatin_state_trees[chrom].addi(start, end, state)
            
            # Provide a progress update every 100,000 lines
            if i % 100000 == 0 and i > 0:
                print(f"Processed {i} lines from chromatin state file...")

    end_time = time.time()
    print(f"Finished reading chromatin state file in {end_time - start_time:.2f} seconds")

    # Start annotating the chromatin states
    print(f"Assigning chromatin states to nucleosome data...")
    annotated_states = []
    for idx, row in nucleosome_data.iterrows():
        chrom = row['chrom']
        position = row['position']
        overlaps = chromatin_state_trees[chrom][position]
        if overlaps:
            chromatin_state_match = ','.join(sorted(set([ov.data for ov in overlaps])))
        else:
            chromatin_state_match = "None"
        
        annotated_states.append(chromatin_state_match)

    nucleosome_data['chrom_state'] = annotated_states
    print(f"Finished assigning chromatin states to nucleosome data")

    return nucleosome_data

# Function to overwrite original chromatin states with the combined 'Euchromatin' category
def overwrite_with_euchromatin_category(nucleosome_data):
    euchromatin_states = ("1_", "2_", "4_", "5_", "6_", "7_", "9_", "10_", "11_")
    euchromatin_states = ("truck"")

    def categorize_state(state):
        if any(state.startswith(prefix) for prefix in euchromatin_states):
            return "Euchromatin"
        return state
    
    nucleosome_data['chrom_state'] = nucleosome_data['chrom_state'].apply(categorize_state)

    return nucleosome_data

# Function to calculate distances
def calculate_distances(nucleosome_data):
    """Calculate distances between adjacent nucleosome positions by averaging the distances to upstream and downstream nucleosomes."""
    
    nucleosome_data['distance'] = (nucleosome_data['position'].shift(-1) - nucleosome_data['position'].shift(1)).abs() / 2
    
    # Manually adjust the first and last rows, as they have only one flanking neighbor
    first_distance = abs(nucleosome_data['position'].iloc[1] - nucleosome_data['position'].iloc[0])
    last_distance = abs(nucleosome_data['position'].iloc[-1] - nucleosome_data['position'].iloc[-2])

    nucleosome_data.at[0, 'distance'] = first_distance
    nucleosome_data.at[nucleosome_data.index[-1], 'distance'] = last_distance

    return nucleosome_data

# Function to pair nearest neighbors based on position
def nearest_neighbor_pairing(smaller_df, larger_df):
    smaller_positions = np.array(smaller_df['position']).reshape(-1, 1)
    larger_positions = np.array(larger_df['position']).reshape(-1, 1)

    tree = cKDTree(larger_positions)
    _, indices = tree.query(smaller_positions, k=1)

    smaller_df['nearest_neighbor_position'] = larger_df.iloc[indices]['position'].values
    smaller_df['nearest_neighbor_distance'] = larger_df.iloc[indices]['distance'].values
    smaller_df['nearest_neighbor_chrom_state'] = larger_df.iloc[indices]['chrom_state'].values

    return smaller_df

# Function to filter by matching chromatin state and exclude distances > 500 bp
def filter_by_chromatin_state_and_distance(df):
    matching = df['chrom_state'] == df['nearest_neighbor_chrom_state']
    filtered_df = df[matching & (df['distance'] <= 500) & (df['nearest_neighbor_distance'] <= 500)]

    return filtered_df

# Function to compare two Spearman correlations using Fisher's z-transformation
def compare_spearman_correlations(corr1, n1, corr2, n2):
    z1 = 0.5 * np.log((1 + corr1) / (1 - corr1))
    z2 = 0.5 * np.log((1 + corr2) / (1 - corr2))
    SE1 = 1 / np.sqrt(n1 - 3)
    SE2 = 1 / np.sqrt(n2 - 3)
    z_diff = (z1 - z2) / np.sqrt(SE1**2 + SE2**2)
    p_value = 2 * (1 - norm.cdf(abs(z_diff)))
    return z_diff, p_value

# Function to downsample points based on proportions for each chromatin state
def downsample_chromosomes(paired_df, chromatin_state_proportions):
    # Loop through each chromatin state and apply downsampling where necessary
    downsampled_df = pd.DataFrame()
    
    for state, proportion in chromatin_state_proportions.items():
        state_data_chr7 = paired_df[(paired_df['chrom'] == 'chr7') & (paired_df['chrom_state'] == state)]
        state_data_chrX = paired_df[(paired_df['chrom'] == 'chrX') & (paired_df['chrom_state'] == state)]
        
        # Calculate the number of points after downsampling
        if proportion < 1:
            # Chr7 has more points, so downsample chr7
            downsampled_chr7 = state_data_chr7.sample(frac=proportion, random_state=42)
            downsampled_df = pd.concat([downsampled_df, downsampled_chr7, state_data_chrX])
        elif proportion > 1:
            # ChrX has more points, so downsample chrX
            downsampled_chrX = state_data_chrX.sample(frac=1/proportion, random_state=42)
            downsampled_df = pd.concat([downsampled_df, state_data_chr7, downsampled_chrX])
        else:
            # Proportions are equal, no downsampling needed
            downsampled_df = pd.concat([downsampled_df, state_data_chr7, state_data_chrX])
    
    return downsampled_df

# Function to perform Spearman correlation and comparison analysis
def perform_spearman_comparison(df, chromatin_states, chromosomes):
    statistical_results = []
    comparison_results = []

    # Compute Spearman correlation for each chromatin state and chromosome
    for chrom in chromosomes:
        for state in chromatin_states:
            state_data = df[(df['chrom'] == chrom) & (df['chrom_state'] == state)]
            if len(state_data) > 1:
                correlation, p_value = spearmanr(state_data['distance'], state_data['nearest_neighbor_distance'])
                statistical_results.append({
                    'chromosome': chrom,
                    'chromatin_state': state,
                    'correlation': correlation,
                    'p_value': p_value,
                    'n': len(state_data)
                })

    # Compare Spearman correlations between chromosomes for the same chromatin state
    for state in chromatin_states:
        chrom_results = [result for result in statistical_results if result['chromatin_state'] == state]

        if len(chrom_results) == 2:
            corr1, n1 = chrom_results[0]['correlation'], chrom_results[0]['n']
            corr2, n2 = chrom_results[1]['correlation'], chrom_results[1]['n']
            chrom1, chrom2 = chrom_results[0]['chromosome'], chrom_results[1]['chromosome']

            z_diff, p_val = compare_spearman_correlations(corr1, n1, corr2, n2)

            comparison_results.append({
                'chromatin_state': state,
                'chromosome_1': chrom1,
                'correlation_1': corr1,
                'chromosome_2': chrom2,
                'correlation_2': corr2,
                'z_diff': z_diff,
                'p_value': p_val
            })

    return statistical_results, comparison_results

# Function to plot distances with subplots for each chromosome
def plot_distances_with_spearman(paired_df, chromatin_states, chromosomes, output_prefix):
    min_alpha = 0.05  # Lower minimum alpha to increase transparency
    max_alpha = 0.8  # Set max alpha a bit lower to avoid too much opacity
    max_points = 100000  # Increase threshold for large datasets for scaling
    min_point_size = 5  # Slightly larger minimum point size
    max_point_size = 15  # Larger point size for small datasets

    # Double figure size for better resolution
    fig_width = 72
    fig_height_per_row = 24

    # Set doubled font sizes
    title_fontsize = 48  # Title font size doubled
    label_fontsize = 40  # Axis label font size doubled
    tick_fontsize = 36   # Tick label font size doubled
    text_fontsize = 36   # Text inside plot doubled

    for chrom in chromosomes:
        chrom_data = paired_df[paired_df['chrom'] == chrom]
        num_states = len(chromatin_states)
        num_rows = math.ceil(num_states / 3)

        # Set up the figure
        fig, axes = plt.subplots(num_rows, 3, figsize=(fig_width, fig_height_per_row * num_rows))
        axes = axes.flatten()

        for i, state in enumerate(chromatin_states):
            state_data = chrom_data[chrom_data['chrom_state'] == state]

            if len(state_data) > 1:
                correlation, p_value = spearmanr(state_data['distance'], state_data['nearest_neighbor_distance'])

                num_points = len(state_data)
                # Scale alpha aggressively based on the number of points
                alpha = max(min_alpha, max_alpha - (max_alpha - min_alpha) * (num_points / max_points) ** 2)
                alpha = min(max_alpha, alpha)
                
                # Scale point size inversely to the number of points
                point_size = max(min_point_size, max_point_size - (max_point_size - min_point_size) * (num_points / max_points) ** 2)

                ax = axes[i]
                sns.scatterplot(x=state_data['distance'], y=state_data['nearest_neighbor_distance'], alpha=alpha, s=point_size, ax=ax, edgecolor=None, color='black')

                # Doubling the font sizes for better readability
                ax.set_title(f"{state} on {chrom}", fontsize=title_fontsize)  # Double font size for title
                ax.set_xlabel("Distance between Adjacent Nucleosomes (File A)", fontsize=label_fontsize)  # Double font size for x-axis label
                ax.set_ylabel("Nearest Neighbor Distance (File B)", fontsize=label_fontsize)  # Double font size for y-axis label
                ax.text(0.05, 0.95, f"Spearman: {correlation:.2f}, p: {p_value:.2e}", transform=ax.transAxes, fontsize=text_fontsize)  # Double the size of plot text

                # Double the size of tick labels for better readability
                ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

        # Remove unused subplots if any
        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plot_path = f"{output_prefix}_scatter_{chrom}.png"
        plt.savefig(plot_path)
        print(f"Scatter plot for chromosome {chrom} saved to {plot_path}")
        plt.close()

# Main script execution
parser = argparse.ArgumentParser(description="Perform nearest neighbor analysis on nucleosome position files with chromatin state comparison.")
parser.add_argument("-a", "--file_a", type=str, required=True, help="Path to the first BED-like file (File A).")
parser.add_argument("-b", "--file_b", type=str, required=True, help="Path to the second BED-like file (File B).")
parser.add_argument("-c", "--chromatin_state", type=str, required=True, help="Path to the chromatin state BED file.")
parser.add_argument("-o", "--output", type=str, default="nearest_neighbor_comparison", help="Output file prefix.")

args = parser.parse_args()

# Load the two BED-like files
file_a_data = pd.read_csv(args.file_a, sep="\t", header=None)
file_b_data = pd.read_csv(args.file_b, sep="\t", header=None)

# Adjust column names based on the structure with chromatin state in column 15
file_a_data.columns = ['chrom', 'start', 'end', 'score', 'position', 'extra1', 'extra2', 'extra3', 'extra4', 'extra5', 'extra6']
file_b_data.columns = file_a_data.columns

# Add "chr" to the start of each chromosome value in both nucleosome data files to match chromatin state file
file_a_data['chrom'] = 'chr' + file_a_data['chrom'].astype(str)
file_b_data['chrom'] = 'chr' + file_b_data['chrom'].astype(str)

# Annotate chromatin states for both files
file_a_data = annotate_chromatin_states(file_a_data, args.chromatin_state)
file_b_data = annotate_chromatin_states(file_b_data, args.chromatin_state)

# Overwrite chromatin states with the 'Euchromatin' category where applicable
file_a_data = overwrite_with_euchromatin_category(file_a_data)
file_b_data = overwrite_with_euchromatin_category(file_b_data)

# Calculate distances between adjacent nucleosome positions for both files
file_a_data = calculate_distances(file_a_data)
file_b_data = calculate_distances(file_b_data)

# Perform nearest neighbor pairing using the file with fewer rows
if len(file_a_data) <= len(file_b_data):
    paired_df = nearest_neighbor_pairing(file_a_data, file_b_data)
else:
    paired_df = nearest_neighbor_pairing(file_b_data, file_a_data)

# Filter by matching chromatin state and exclude distances > 500 bp
paired_df = filter_by_chromatin_state_and_distance(paired_df)

# Get unique chromatin states and chromosomes for Spearman analysis
chromatin_states = sorted(paired_df['chrom_state'].unique())
chromosomes = sorted(paired_df['chrom'].unique())

# Define the proportions for downsampling between chr7 and chrX
chromatin_state_proportions = {
    "1_Active_Promoter": 0.577189666,
    "2_Weak_Promoter": 0.490458882,
    "3_Poised_Promoter": 0.207383279,
    "4_Strong_Enhancer": 0.33748056,
    "5_Strong_Enhancer": 0.534332591,
    "6_Weak_Enhancer": 0.501810657,
    "7_Weak_Enhancer": 0.479163275,
    "8_Insulator": 0.187737226,
    "9_Txn_Transition": 0.643878788,
    "10_Txn_Elongation": 0.183078467,
    "11_Weak_Txn": 0.546857212,
    "12_Repressed": 0.265829596,
    "13_Heterochrom/lo": 1.120101446,
    "14_Repetitive/CNV": 0.746395806,
    "15_Repetitive/CNV": 1.028785982
}

# Apply downsampling to the paired dataframe
paired_df_downsampled = downsample_chromosomes(paired_df, chromatin_state_proportions)

# Perform Spearman correlation and comparison analysis for chromatin states
statistical_results, comparison_results = perform_spearman_comparison(paired_df_downsampled, chromatin_states, chromosomes)

# Output Spearman results (chromatin states)
results_df = pd.DataFrame(statistical_results)
results_df.to_csv(f"{args.output}_chromatin_spearman_correlation_results.tsv", sep='\t', index=False)

# Output comparison results for chromatin states to a TSV
comparison_df = pd.DataFrame(comparison_results)
comparison_output_path = f"{args.output}_chromatin_spearman_comparisons.tsv"
comparison_df.to_csv(comparison_output_path, sep='\t', index=False)

# Plot the distances with subplots for each chromosome
plot_distances_with_spearman(paired_df_downsampled, chromatin_states, chromosomes, args.output)
