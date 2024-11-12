import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
from intervaltree import IntervalTree
from scipy.stats import kruskal
import seaborn as sns

# Set matplotlib to use the Agg backend for environments without display
plt.switch_backend('Agg')

def annotate_chromatin_states(nucleosome_data, chromatin_state_file):
    chromatin_state_trees = defaultdict(IntervalTree)

    # Read the chromatin state BED file and store the intervals in IntervalTree for each chromosome
    with open(chromatin_state_file, 'r') as csfile:
        for line in csfile:
            cols = line.strip().split()
            chrom = f"{cols[0]}"
            start = int(cols[1])
            end = int(cols[2])
            state = cols[3]
            chromatin_state_trees[chrom].addi(start, end, state)

    # Annotate chromatin states for each nucleosome position
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
    return nucleosome_data

def filter_chromatin_states(data, exclude_states):
    """Exclude chromatin states that start with specific values or are 'None'."""
    filtered_data = data[~data['chrom_state'].str.startswith(tuple(exclude_states))]
    filtered_data = filtered_data[filtered_data['chrom_state'] != "None"]
    return filtered_data

def create_euchromatin_category(data):
    """Create a new category 'Euchromatin' by extracting relevant states and aggregating them."""
    euchromatin_states = ("1_", "2_", "4_", "5_", "6_", "7_", "9_", "10_", "11_")
    
    # Extract distances for chromatin states that match Euchromatin categories
    euchromatin_data = data[data['chrom_state'].str.startswith(euchromatin_states)].copy()
    
    # Assign 'Euchromatin' to these rows
    euchromatin_data['chrom_state'] = 'Euchromatin'
    
    return euchromatin_data

def sort_chromatin_states(chromatin_states):
    """Sort chromatin states by numeric prefix, with 'Euchromatin' at the end."""
    # Extract the numeric prefix (before the underscore)
    sorted_states = sorted(
        [state for state in chromatin_states if state != 'Euchromatin'], 
        key=lambda x: int(x.split('_')[0])
    )
    # Add 'Euchromatin' at the end
    sorted_states.append('Euchromatin')
    return sorted_states

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Perform nearest neighbor analysis on nucleosome position files with chromatin state comparison.")
parser.add_argument("-a", "--file_a", type=str, required=True, help="Path to the first BED-like file (File A).")
parser.add_argument("-b", "--file_b", type=str, required=True, help="Path to the second BED-like file (File B).")
parser.add_argument("-c", "--chromatin_state", type=str, required=True, help="Path to the chromatin state BED file.")
parser.add_argument("--file_c", type=str, help="Path to the alternative BED-like file (File C) for chrX.")
parser.add_argument("-o", "--output", type=str, default="nearest_neighbor_comparison", help="Output file prefix.")
parser.add_argument("--score_percentile", type=int, default=0, help="Percentile threshold for scores. Only paired peaks above this percentile will be used.")

# Parse the arguments
args = parser.parse_args()

# Load the two BED-like files
file_a_data = pd.read_csv(args.file_a, sep="\t", header=None)
file_b_data = pd.read_csv(args.file_b, sep="\t", header=None)

# Adjust column names based on the structure with chromatin state in column 15
file_a_data.columns = ['chrom', 'start', 'end', 'score', 'position', 'extra1', 'extra2', 'extra3', 'extra4', 'extra5', 'extra6']
file_b_data.columns = file_a_data.columns

# Add "chr" to the start of each chromosome value
file_a_data['chrom'] = 'chr' + file_a_data['chrom'].astype(str)
file_b_data['chrom'] = 'chr' + file_b_data['chrom'].astype(str)

# If file_c is provided, load it for chromosome X
if args.file_c:
    file_c_data = pd.read_csv(args.file_c, sep="\t", header=None)
    file_c_data.columns = file_a_data.columns
    file_c_data['chrom'] = 'chr' + file_c_data['chrom'].astype(str)
else:
    file_c_data = None

# Annotate chromatin states for both files
file_a_data = annotate_chromatin_states(file_a_data, args.chromatin_state)
file_b_data = annotate_chromatin_states(file_b_data, args.chromatin_state)
if file_c_data is not None:
    file_c_data = annotate_chromatin_states(file_c_data, args.chromatin_state)

# Normalize the scores by (score / mean * 100)
file_a_data['normalized_score'] = (file_a_data['score'] / file_a_data['score'].mean()) * 100
file_b_data['normalized_score'] = (file_b_data['score'] / file_b_data['score'].mean()) * 100
if file_c_data is not None:
    file_c_data['normalized_score'] = (file_c_data['score'] / file_c_data['score'].mean()) * 100

# Get the unique chromosomes
chromosomes = sorted(set(file_a_data['chrom'].unique()) & set(file_b_data['chrom'].unique()))

# Initialize data structures to hold all data for plotting
all_matched_data = pd.DataFrame()

# Loop through each chromosome and perform analysis
for chrom in chromosomes:
    print(f"Processing chromosome {chrom}")

    # Subset data by chromosome
    file_a_chrom = file_a_data[file_a_data['chrom'] == chrom].copy()

    # Use file_c data for chrX if provided
    if chrom == 'chrX' and file_c_data is not None:
        file_b_chrom = file_c_data[file_c_data['chrom'] == chrom].copy()
    else:
        file_b_chrom = file_b_data[file_b_data['chrom'] == chrom].copy()

    # Use scipy's cKDTree to find nearest neighbors (pair peaks)
    ref_positions = np.array(file_a_chrom['position']).reshape(-1, 1)
    comp_positions = np.array(file_b_chrom['position']).reshape(-1, 1)

    tree = cKDTree(comp_positions)
    distances, indices = tree.query(ref_positions, k=1)

    # Get the nearest neighbors' positions, normalized scores, and chromatin states from the comparison file
    nearest_neighbors = file_b_chrom.iloc[indices]

    # Calculate the absolute distance and difference in normalized scores
    file_a_chrom.loc[:, 'nearest_distance'] = np.abs(nearest_neighbors['position'].values - file_a_chrom['position'].values)
    file_a_chrom.loc[:, 'nearest_position'] = nearest_neighbors['position'].values
    file_a_chrom.loc[:, 'nearest_normalized_score'] = nearest_neighbors['normalized_score'].values
    file_a_chrom.loc[:, 'score_diff'] = np.abs(file_a_chrom['normalized_score'].values - file_a_chrom['nearest_normalized_score'].values)
    file_a_chrom.loc[:, 'nearest_chromatin_state'] = nearest_neighbors['chrom_state'].values

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

# Exclude chromatin states starting with 14_, 15_, and None
exclude_states = ['14_', '15_']
filtered_data = filter_chromatin_states(all_matched_data, exclude_states)

# Create the 'Euchromatin' category and add it to the filtered data
euchromatin_data = create_euchromatin_category(filtered_data)
filtered_data = pd.concat([filtered_data, euchromatin_data])

# Sort chromatin states by numeric prefix and add 'Euchromatin' at the end
chromatin_states = sort_chromatin_states(filtered_data['chrom_state'].unique())

# Perform statistical analysis for each chromatin state
statistical_results = []

# Loop through each chromatin state
for state in chromatin_states:
    state_data = filtered_data[filtered_data['chrom_state'] == state]
    chromosomes = sorted(state_data['chrom'].unique())
    
    # Collect the distance data for each chromosome and compute medians
    medians = {}
    distance_data_per_chrom = []
    
    for chrom in chromosomes:
        chrom_data = state_data[state_data['chrom'] == chrom]['nearest_distance'].values
        distance_data_per_chrom.append(chrom_data)
        medians[chrom] = np.median(chrom_data)
    
    # Perform Kruskal-Wallis or Mann-Whitney U test depending on the number of chromosomes
    if len(chromosomes) > 2:
        stat, p_value = kruskal(*distance_data_per_chrom)
        test_name = 'Kruskal-Wallis'
    elif len(chromosomes) == 2:
        stat, p_value = mannwhitneyu(distance_data_per_chrom[0], distance_data_per_chrom[1])
        test_name = 'Mann-Whitney U'
    else:
        stat = None
        p_value = None
        test_name = 'None'
    
    # Store the results
    statistical_results.append({
        'chromatin_state': state,
        'test': test_name,
        'stat': stat,
        'p_value': p_value,
        **{f'median_{chrom}': medians[chrom] for chrom in chromosomes}
    })

# Output statistical results to TSV
stat_results_df = pd.DataFrame(statistical_results)
output_tsv_path = f"{args.output}_statistical_results.tsv"
stat_results_df.to_csv(output_tsv_path, sep='\t', index=False)

print(f"Statistical results saved to {output_tsv_path}")

# Prepare data for the combined plot
combined_distances = []
combined_scores = []
chromatin_labels = []
chromosome_labels = []

# Loop through each chromatin state and chromosome, collecting the data
for state in chromatin_states:
    state_data = filtered_data[filtered_data['chrom_state'] == state]
    for chrom in sorted(state_data['chrom'].unique()):
        distances = state_data[state_data['chrom'] == chrom]['nearest_distance']
        scores = state_data[state_data['chrom'] == chrom]['score_diff']
        combined_distances.extend(distances)  # Flatten distances
        combined_scores.extend(scores)  # Flatten scores
        chromatin_labels.extend([state] * len(distances))  # Repeat chromatin state for each distance
        chromosome_labels.extend([chrom] * len(distances))  # Repeat chromosome for each distance

# Ensure arrays have the same length
assert len(combined_distances) == len(chromatin_labels) == len(chromosome_labels), "Array lengths are not equal!"
assert len(combined_scores) == len(chromatin_labels), "Array lengths for scores are not equal!"

# Create DataFrames for plotting
distance_plot_df = pd.DataFrame({
    'Distance': combined_distances,
    'Chromatin State': chromatin_labels,
    'Chromosome': chromosome_labels
})

score_plot_df = pd.DataFrame({
    'Score Difference': combined_scores,
    'Chromatin State': chromatin_labels,
    'Chromosome': chromosome_labels
})

# Create the distance boxplot without outliers
fig_combined_distance, ax_combined_distance = plt.subplots(figsize=(12, 8))
sns.boxplot(x='Chromatin State', y='Distance', hue='Chromosome', data=distance_plot_df, palette="Set3", ax=ax_combined_distance, showfliers=False)

# Customize plot labels and appearance
ax_combined_distance.set_title('Chromatin State and Chromosome Distance Comparison')
ax_combined_distance.set_xlabel('Chromatin State')
ax_combined_distance.set_ylabel('Distance')
ax_combined_distance.legend(title='Chromosome')

# Save the combined distance boxplot
plt.tight_layout()
fig_combined_distance.savefig(f"{args.output}_chromatin_chromosome_distance_boxplot.svg")

# Create the score difference boxplot without outliers
fig_combined_score, ax_combined_score = plt.subplots(figsize=(12, 8))
sns.boxplot(x='Chromatin State', y='Score Difference', hue='Chromosome', data=score_plot_df, palette="Set3", ax=ax_combined_score, showfliers=False)

# Customize plot labels and appearance
ax_combined_score.set_title('Chromatin State and Chromosome Score Difference Comparison')
ax_combined_score.set_xlabel('Chromatin State')
ax_combined_score.set_ylabel('Score Difference')
ax_combined_score.legend(title='Chromosome')

# Save the combined score boxplot
plt.tight_layout()
fig_combined_score.savefig(f"{args.output}_chromatin_chromosome_score_boxplot.svg")

# Display messages after completion
print(f"Combined distance boxplot saved as {args.output}_chromatin_chromosome_distance_boxplot.svg")
print(f"Combined score boxplot saved as {args.output}_chromatin_chromosome_score_boxplot.svg")
