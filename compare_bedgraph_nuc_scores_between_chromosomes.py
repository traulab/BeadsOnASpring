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

# Function to open a file, handling gzipped files as well as plain text files
def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r') 

# Function to count total lines in a bedGraph file for proportion calculation
def count_lines(file_path):
    with open_file(file_path) as f:
        return sum(1 for _ in f)

# Function to process bedGraph files for two chromosomes with downsampling and a loading bar
def process_bedgraph_for_chromosome_comparison(chrom1_file, chrom2_file, chrom1_name, chrom2_name):
    accumulated_data = []

    print(f"Processing bedGraph files for {chrom1_name} and {chrom2_name} with potential downsampling")

    # Count total lines in one of the files (assumes equal number of lines)
    total_lines = count_lines(chrom1_file)
    
    # Calculate downsampling proportion
    proportion = 1.0
    if total_lines > MAX_POINTS:
        proportion = MAX_POINTS / total_lines
        print(f"Applying downsampling with proportion: {proportion:.4f} (based on {total_lines} total lines)")

    with open_file(chrom1_file) as chrom1_f, open_file(chrom2_file) as chrom2_f:
        with tqdm(total=total_lines, desc="Processing bedGraph lines", unit="lines") as pbar:
            for line1, line2 in zip(chrom1_f, chrom2_f):
                # Apply random sampling based on the proportion
                if random.random() <= proportion:
                    cols1 = line1.strip().split()
                    cols2 = line2.strip().split()

                    chrom1_start = int(cols1[1])
                    chrom2_start = int(cols2[1])
                    chrom1_score = float(cols1[5])
                    chrom2_score = float(cols2[5])

                    if chrom1_score == 0 and chrom2_score == 0:
                        continue

                    # Ensure that start positions match
                    if chrom1_start == chrom2_start:
                        accumulated_data.append({
                            'start': chrom1_start,
                            'chrom1_score': chrom1_score,
                            'chrom2_score': chrom2_score
                        })

                pbar.update(1)

    combined_data = pd.DataFrame(accumulated_data)

    # Truncate to the length of the shorter chromosome if needed
    max_length = min(hg19_chromosome_lengths[chrom1_name], hg19_chromosome_lengths[chrom2_name])
    combined_data = combined_data[combined_data['start'] < max_length]

    # Filter extreme scores
    threshold_score = 500
    combined_data = combined_data[(combined_data['chrom1_score'] <= threshold_score) & (combined_data['chrom2_score'] <= threshold_score)]
    combined_data = combined_data[(combined_data['chrom1_score'] >= -threshold_score) & (combined_data['chrom2_score'] >= -threshold_score)]

    return combined_data

# Function to standardize (z-score normalize) scores based on mean and standard deviation
def standardize_scores(df):
    df['chrom1_score_z'] = (df['chrom1_score'] - df['chrom1_score'].mean()) / df['chrom1_score'].std()
    df['chrom2_score_z'] = (df['chrom2_score'] - df['chrom2_score'].mean()) / df['chrom2_score'].std()
    
    print("Chromosome 1 score mean and standard deviation (before standardization):", df['chrom1_score'].mean(), df['chrom1_score'].std())
    print("Chromosome 2 score mean and standard deviation (before standardization):", df['chrom2_score'].mean(), df['chrom2_score'].std())
    
    return df

# Custom lightened colormap function for yellow-to-red transition
def lighten_color_map(cmap_name, min_value=0.2, max_value=1.0):
    cmap = plt.get_cmap(cmap_name)
    new_cmap = cmap(np.linspace(min_value, max_value, cmap.N))
    return ListedColormap(new_cmap)

def plot_scores(data, output_prefix, correlation, p_value, cmap_name):
    if len(data) > 0:
        plt.figure(figsize=(8, 8))
        plt.rcParams.update({'font.size': 16})

        max_abs_value = max(data['chrom1_score_z'].abs().max(), data['chrom2_score_z'].abs().max())
        hb = plt.hexbin(data['chrom1_score_z'], data['chrom2_score_z'], gridsize=600, cmap=lighten_color_map(cmap_name, 0.3, 1.0), bins='log')

        counts = hb.get_array()
        norm = Normalize(vmin=counts.min(), vmax=counts.max())
        alphas = norm(counts)

        colors = plt.get_cmap(cmap_name)(norm(counts))
        for i, alpha in enumerate(alphas):
            colors[i, -1] = alpha
        hb.set_facecolor(colors)

        plt.title('Chromosome Comparison', fontsize=18)
        plt.xlabel('Chromosome 1 Score (z-score)', fontsize=16)
        plt.ylabel('Chromosome 2 Score (z-score)', fontsize=16)

        max_abs_value = 20
        plt.xlim(-max_abs_value, max_abs_value)
        plt.ylim(-max_abs_value, max_abs_value)

        plt.text(0.05, 0.95, f'Correlation: {correlation:.2f}\nP-value: {p_value:.2e}', fontsize=14, ha='left', va='top', transform=plt.gca().transAxes)

        plot_file = f"{output_prefix}_chr_comparison_density.png"
        plt.savefig(plot_file, dpi=300)
        print(f"Density plot with adjusted alpha and custom colormap saved to {plot_file}")

        plt.close()
    else:
        print("No data available for chromosome comparison, skipping plot.")

# Main script execution
parser = argparse.ArgumentParser(description="Compare two chromosome bedGraph scores and correlate.")
parser.add_argument("--chrom1_file", type=str, required=True, help="Path to the first chromosome bedGraph file.")
parser.add_argument("--chrom2_file", type=str, required=True, help="Path to the second chromosome bedGraph file.")
parser.add_argument("--chrom1_name", type=str, required=True, help="Name of the first chromosome (e.g., chrX).")
parser.add_argument("--chrom2_name", type=str, required=True, help="Name of the second chromosome (e.g., chr7).")
parser.add_argument("-o", "--output", type=str, default="chr_comparison", help="Output file prefix.")
parser.add_argument("--colormap", type=str, default="YlOrRd_r", help="Specify a colormap for the hexbin plot (default: YlOrRd_r).")

args = parser.parse_args()

# Process the bedGraph files for comparison
merged_data = process_bedgraph_for_chromosome_comparison(args.chrom1_file, args.chrom2_file, args.chrom1_name, args.chrom2_name)

# Standardize the scores (z-score normalization)
merged_data = standardize_scores(merged_data)

# Perform Pearson correlation analysis
if len(merged_data) > 1:
    correlation, p_value = pearsonr(merged_data['chrom1_score_z'], merged_data['chrom2_score_z'])
    print(f"Correlation: {correlation}, p-value: {p_value}")

    # Plot the scores
    plot_scores(merged_data, args.output, correlation, p_value, args.colormap)
else:
    print("No data available for correlation.")

# Save correlation results
output_file = f"{args.output}_correlation_results.tsv"
correlation_results = pd.DataFrame([{
    'chromosome1': args.chrom1_name,
    'chromosome2': args.chrom2_name,
    'correlation': correlation,
    'p_value': p_value,
    'n': len(merged_data)
}])
correlation_results.to_csv(output_file, sep='\t', index=False)
print(f"Correlation results saved to {output_file}")
