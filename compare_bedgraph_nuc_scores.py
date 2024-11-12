import os
import pandas as pd
import numpy as np
import gzip
from tqdm import tqdm
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import argparse
import random
from matplotlib.colors import Normalize, ListedColormap

# Set matplotlib to use the Agg backend for environments without display
plt.switch_backend('Agg')

# Define the lengths of hg19 chromosomes
hg19_chromosome_lengths = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566
}

# Limit number of points to read
MAX_POINTS = 10000000

# Function to count the total number of lines in a bedGraph file
def count_lines(file_path):
    with open_file(file_path) as f:
        return sum(1 for _ in f)

# Function to open a file, handling gzipped files as well as plain text files
def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r') 

# Function to process bedGraph files line by line with progress bar and apply random downsampling
def process_bedgraph_with_downsampling(male_file, female_file, total_lines):
    accumulated_data = []

    print(f"Processing bedGraph files line by line with potential downsampling")

    # Calculate downsampling proportion
    proportion = 1.0
    if total_lines > MAX_POINTS:
        proportion = MAX_POINTS / total_lines
        print(f"Applying downsampling with proportion: {proportion:.4f} (based on {total_lines} total lines)")

    with open_file(male_file) as male_f, open_file(female_file) as female_f:
        with tqdm(total=total_lines, desc="Processing bedGraph lines", unit="lines") as pbar:
            for male_line, female_line in zip(male_f, female_f):
                # Apply random sampling based on the proportion
                if random.random() <= proportion:
                    male_cols = male_line.strip().split()
                    male_chrom = 'chr' + male_cols[0]
                    male_start = int(male_cols[1])
                    male_score = float(male_cols[5])

                    female_cols = female_line.strip().split()
                    female_chrom = 'chr' + female_cols[0]
                    female_start = int(female_cols[1])
                    female_score = float(female_cols[5])

                    if male_score == 0 and female_score == 0:
                        continue

                    if male_chrom == female_chrom and male_start == female_start:
                        accumulated_data.append({
                            'chrom': male_chrom,
                            'start': male_start,
                            'male_score': male_score,
                            'female_score': female_score
                        })

                pbar.update(1)

    combined_data = pd.DataFrame(accumulated_data)

    threshold_score = 500
    combined_data = combined_data[(combined_data['male_score'] <= threshold_score) & (combined_data['female_score'] <= threshold_score)]
    combined_data = combined_data[(combined_data['male_score'] >= -threshold_score) & (combined_data['female_score'] >= -threshold_score)]

    return combined_data

# Function to standardize (z-score normalize) scores based on mean and standard deviation
def standardize_scores(df):
    """
    Standardize (z-score) the male and female scores in the dataframe.
    
    Args:
    - df: The dataframe containing 'male_score' and 'female_score'.
    
    Returns:
    - df: The dataframe with standardized scores.
    """
    df['male_score_z'] = (df['male_score'] - df['male_score'].mean()) / df['male_score'].std()
    df['female_score_z'] = (df['female_score'] - df['female_score'].mean()) / df['female_score'].std()
    
    print("Male score mean and standard deviation (before standardization):", df['male_score'].mean(), df['male_score'].std())
    print("Female score mean and standard deviation (before standardization):", df['female_score'].mean(), df['female_score'].std())
    
    return df

# Custom lightened colormap function for yellow-to-red transition
def lighten_color_map(cmap_name, min_value=0.2, max_value=1.0):
    cmap = plt.get_cmap(cmap_name)
    new_cmap = cmap(np.linspace(min_value, max_value, cmap.N))
    return ListedColormap(new_cmap)

def plot_scores(state_data, chromatin_state, output_prefix, correlation, p_value, cmap_name):
    if len(state_data) > 0:

        plt.figure(figsize=(8, 8))
        plt.rcParams.update({'font.size': 16})

        # Calculate the maximum absolute value from both male and female scores
        max_abs_value = max(state_data['male_score_z'].abs().max(), state_data['female_score_z'].abs().max())

        # Create the hexbin plot with the provided colormap
        hb = plt.hexbin(state_data['male_score_z'], state_data['female_score_z'], gridsize=600, 
                        cmap=lighten_color_map(cmap_name, 0.3, 1.0), bins='log')

        # Normalize counts to use them for adjusting alpha
        counts = hb.get_array()
        norm = Normalize(vmin=counts.min(), vmax=counts.max())
        alphas = norm(counts)

        # Get the colors and modify alpha in the RGBA colors
        colors = plt.get_cmap(cmap_name)(norm(counts))  # Get RGBA values from the custom colormap
        for i, alpha in enumerate(alphas):
            colors[i, -1] = alpha  # Adjust the alpha value

        # Update the hexbin collection with the new colors (including alpha)
        hb.set_facecolor(colors)

        plt.title(f'Chromatin State: {chromatin_state}', fontsize=18)
        plt.xlabel('Male Score (z-score)', fontsize=16)
        plt.ylabel('Female Score (z-score)', fontsize=16)

        max_abs_value = 10

        # Set the x and y limits to be the same, based on the max absolute value
        plt.xlim(-max_abs_value, max_abs_value)
        plt.ylim(-max_abs_value, max_abs_value)

        # Annotate with correlation and p-value
        plt.text(0.05, 0.95, f'Correlation: {correlation:.2f}\nP-value: {p_value:.2e}', fontsize=14, ha='left', va='top', transform=plt.gca().transAxes)

        plot_file = f"{output_prefix}_chromatin_state_{chromatin_state}_male_vs_female_density.png"
        plt.savefig(plot_file, dpi=300)
        print(f"Density plot with adjusted alpha and custom colormap saved to {plot_file}")

        plt.close()
    else:
        print(f"No data available for chromatin state {chromatin_state}, skipping plot.")

# Main script execution
parser = argparse.ArgumentParser(description="Compare male and female bedGraph scores and correlate by chromatin state.")
parser.add_argument("-m", "--male_file", type=str, required=True, help="Path to the male bedGraph file.")
parser.add_argument("-f", "--female_file", type=str, required=True, help="Path to the female bedGraph file.")
parser.add_argument("-s", "--chromatin_state", type=str, required=True, help="The chromatin state being used.")
parser.add_argument("-o", "--output", type=str, default="male_vs_female_comparison", help="Output file prefix.")
parser.add_argument("--chromosome", type=str, required=True, help="Specify the chromosome (e.g., chr7)")
parser.add_argument("--colormap", type=str, default="YlOrRd_r", help="Specify a colormap for the hexbin plot (default: YlOrRd_r).")

args = parser.parse_args()

# Count the total lines in the male file (since the number of lines is the same for both male and female)
total_lines = count_lines(args.male_file)
print(f"Total number of lines in the bedGraph file: {total_lines}")

# Process the bedGraph files with potential downsampling
merged_data = process_bedgraph_with_downsampling(args.male_file, args.female_file, total_lines)

# Standardize the male and female scores (z-score normalization)
merged_data = standardize_scores(merged_data)

# Perform Spearman correlation analysis for the provided chromatin state
if len(merged_data) > 1:
    correlation, p_value = pearsonr(merged_data['male_score_z'], merged_data['female_score_z'])
    spearman_results = [{
        'chromatin_state': args.chromatin_state,
        'correlation': correlation,
        'p_value': p_value,
        'n': len(merged_data)
    }]
    print(f"Chromatin State: {args.chromatin_state}, Correlation: {correlation}, p-value: {p_value}")

    # Plot the scores with the provided colormap
    plot_scores(merged_data, args.chromatin_state, args.output, correlation, p_value, args.colormap)
else:
    print(f"No data available for chromatin state {args.chromatin_state}, skipping correlation.")

# Save correlation results
output_file = f"{args.output}_spearman_correlation_results_chromatin_state_{args.chromatin_state}.tsv"
spearman_results_df = pd.DataFrame(spearman_results)
spearman_results_df.to_csv(output_file, sep='\t', index=False)
print(f"Spearman correlation results saved to {output_file}")