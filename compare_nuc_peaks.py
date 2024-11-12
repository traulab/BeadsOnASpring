import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import argparse
import matplotlib.pyplot as plt

# Set matplotlib to use the Agg backend for environments without display
plt.switch_backend('Agg')

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Perform nearest neighbor analysis on nucleosome position files.")
parser.add_argument("-a", "--file_a", type=str, required=True, help="Path to the first BED-like file (File A).")
parser.add_argument("-b", "--file_b", type=str, required=True, help="Path to the second BED-like file (File B).")
parser.add_argument("-o", "--output", type=str, default="nearest_neighbor_comparison", help="Output file prefix.")
parser.add_argument("--score_percentile", type=int, default=0, help="Percentile threshold for scores. Only paired peaks above this percentile will be used.")

# Parse the arguments
args = parser.parse_args()

# Load the two BED-like files
file_a_data = pd.read_csv(args.file_a, sep="\t", header=None)
file_b_data = pd.read_csv(args.file_b, sep="\t", header=None)

# Adjust column names based on the structure with chromatin state in column 15
file_a_data.columns = ['chrom', 'start', 'end', 'score', 'position', 'extra1', 'extra2', 'extra3', 'extra4', 'extra5',
                       'extra6', 'b_chrom', 'b_start', 'b_end', 'chrom_state', 'b_score', 'b_strand', 
                       'b_thickstart', 'b_thickend', 'b_rgb']
file_b_data.columns = file_a_data.columns

# Normalize the scores by (score / mean * 100)
file_a_data['normalized_score'] = (file_a_data['score'] / file_a_data['score'].mean()) * 100
file_b_data['normalized_score'] = (file_b_data['score'] / file_b_data['score'].mean()) * 100

# Get the unique chromosomes
chromosomes = sorted(set(file_a_data['chrom'].unique()) & set(file_b_data['chrom'].unique()))

# Initialize data structures to hold all data for plotting
all_matched_data = pd.DataFrame()

# Loop through each chromosome and perform analysis
for chrom in chromosomes:
    print(f"Processing chromosome {chrom}")

    # Subset data by chromosome
    file_a_chrom = file_a_data[file_a_data['chrom'] == chrom]
    file_b_chrom = file_b_data[file_b_data['chrom'] == chrom]

    # Use scipy's cKDTree to find nearest neighbors (pair peaks)
    ref_positions = np.array(file_a_chrom['position']).reshape(-1, 1)
    comp_positions = np.array(file_b_chrom['position']).reshape(-1, 1)

    tree = cKDTree(comp_positions)
    distances, indices = tree.query(ref_positions, k=1)

    # Get the nearest neighbors' positions, normalized scores, and chromatin states from the comparison file
    nearest_neighbors = file_b_chrom.iloc[indices]

    # Calculate the signed distance and difference in normalized scores
    file_a_chrom['nearest_distance'] = nearest_neighbors['position'].values - file_a_chrom['position']
    file_a_chrom['nearest_position'] = nearest_neighbors['position'].values
    file_a_chrom['nearest_normalized_score'] = nearest_neighbors['normalized_score'].values
    file_a_chrom['score_diff'] = file_a_chrom['normalized_score'] - file_a_chrom['nearest_normalized_score']
    file_a_chrom['nearest_chromatin_state'] = nearest_neighbors['chrom_state'].values

    # Filter by matching chromatin states between the two files
    matched_data = file_a_chrom[file_a_chrom['chrom_state'] == file_a_chrom['nearest_chromatin_state']]

    # Apply score percentile filtering after pairing
    if args.score_percentile > 0:
        score_thresh_a = np.percentile(matched_data['normalized_score'], args.score_percentile)
        score_thresh_b = np.percentile(matched_data['nearest_normalized_score'], args.score_percentile)
        matched_data = matched_data[(matched_data['normalized_score'] > score_thresh_a) & (matched_data['nearest_normalized_score'] > score_thresh_b)]

    matched_data.loc[:, 'chrom'] = chrom  # Add chromosome information

    # Append the matched data to the full dataset
    all_matched_data = pd.concat([all_matched_data, matched_data], ignore_index=True)

# Get the list of unique chromatin states
chromatin_states = sorted(all_matched_data['chrom_state'].unique())

# Create subplots for boxplots and histograms
fig_boxplot, ax_boxplot = plt.subplots(len(chromatin_states), 1, figsize=(12, 6 * len(chromatin_states)))
fig_hist, ax_hist = plt.subplots(len(chromatin_states), 1, figsize=(12, 6 * len(chromatin_states)))

# Plot the boxplots and histograms for each chromatin state
for i, state in enumerate(chromatin_states):
    state_data = all_matched_data[all_matched_data['chrom_state'] == state]
    chromosomes = sorted(state_data['chrom'].unique())
    
    distance_data = [state_data[state_data['chrom'] == chrom]['nearest_distance'] for chrom in chromosomes]
    score_diff_data = [state_data[state_data['chrom'] == chrom]['score_diff'] for chrom in chromosomes]
    
    # Boxplot for distances and score differences
    ax_boxplot[i].boxplot(distance_data, positions=np.arange(len(chromosomes)), widths=0.4, patch_artist=True, 
                          boxprops=dict(facecolor="blue", alpha=0.6), showfliers=False, tick_labels=chromosomes)
    ax_boxplot[i].boxplot(score_diff_data, positions=np.arange(len(chromosomes)) + 0.5, widths=0.4, patch_artist=True, 
                          boxprops=dict(facecolor="green", alpha=0.6), showfliers=False)
    ax_boxplot[i].set_xticks(np.arange(len(chromosomes)) + 0.25)
    ax_boxplot[i].set_xticklabels(chromosomes)
    ax_boxplot[i].set_title(f'Comparison for Chromatin State: {state}')
    ax_boxplot[i].set_ylabel('Distance (blue) / Score Difference (green)')

    # Histograms for distances and score differences
    for chrom in chromosomes:
        chrom_data = state_data[state_data['chrom'] == chrom]
        ax_hist[i].hist(chrom_data['nearest_distance'], bins=np.arange(-150, 151, 1), density=True, alpha=0.5, 
                        label=f'Chrom {chrom} Distance', histtype='step', linewidth=2)
        ax_hist[i].hist(chrom_data['score_diff'], bins=np.arange(-150, 151, 1), density=True, alpha=0.5, 
                        label=f'Chrom {chrom} Score Diff', histtype='step', linestyle='dotted', linewidth=2)
    ax_hist[i].set_xlim([-150, 150])
    ax_hist[i].set_title(f'Histograms for Chromatin State: {state}')
    ax_hist[i].set_xlabel('Distance / Score Difference')
    ax_hist[i].set_ylabel('Density')
    ax_hist[i].legend(loc='upper right')

# Save the combined boxplot and histogram figures
plt.tight_layout()
fig_boxplot.savefig(f"{args.output}_boxplot_comparison.png")
fig_hist.savefig(f"{args.output}_histogram_comparison.png")

# Display a message after completion
print(f"Chromatin state comparison complete. Results saved with prefix {args.output}.")
