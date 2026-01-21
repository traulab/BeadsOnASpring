## Table of Contents

- [PNS with Nucleosome Peak Calling](#PNS-with-nucleosome-peak-calling)
- [Requirements](#requirements)
- [PNS: Probabilistic Nucleosome Score](#PNS-probabilistic-nucleosome-score)
  - [Key Concepts](#key-concepts)
- [`PNS_with_nucleosome_peak_calling.py`](#PNS_with_nucleosome_peak_callingpy)
  1. [Running on Individual Contigs](#1-running-on-individual-contigs)
  2. [Output Files](#2-output-files)
     - [Combined Scores (bedGraph) File Format](#a-combined-scores-bedgraph-file-format)
     - [Nucleosome Regions (bed) File Format](#b-nucleosome-regions-bed-file-format)
- [Setting up the Conda/Mamba Environment](#setting-up-the-condamamba-environment)
  1. [Installing Conda](#1-installing-conda)
  2. [Installing Mamba via Conda](#2-installing-mamba-via-conda)
  3. [Creating the `PNS` Environment](#3-creating-the-PNS-environment)
  4. [Activating the `PNS` Environment](#4-activating-the-PNS-environment)
  5. [Verify Installation](#5-verify-installation)
  6. [Deactivating the Environment](#6-deactivating-the-environment)

---

# PNS with Nucleosome Peak Calling

This repository provides tools for calling nucleosome peaks from BAM files using nucleosome positioning analysis. Primary script: 

`PNS_with_nucleosome_peak_calling.py`: Calculates Probabilistic Nucleosome Score (PNS) and identifies nucleosome peak regions.

## Requirements

- **Python 3**
- **Required Python libraries**: `pandas`, `numpy`, `scipy`, `pysam`, `matplotlib`, `tqdm`, `pyBigWig`
- **Additional tools**: `samtools`, `bedGraphToBigWig`

---

## PNS: Probabilistic Nucleosome Score

The Probabilistic Nucleosome Score (PNS) measures the degree of confidence in nucleosome positioning based on the length and distribution of paired-end DNA fragments. Each fragment is assigned a probability distribution based on the mode fragment length of the sample.

![PNS Example](IGV_PNS_example.png)

### Key Concepts:

1. **Mode Fragment Length**: This is the most frequent fragment length in the sample, representing the most probable size of a nucleosome-protected fragment. Mode-sized fragments have the highest confidence for determining nucleosome positioning.
2. **Dyad Likelihood Distributions**: Each fragment is assigned a dyad likelihood distribution reflecting the confidence in nucleosome positioning. Mode-sized fragments have the greatest confidence of nucleosome positioning at the center (+1) and breakpoints at the start and end (-1). Non-mode fragments are given broader distributions to reflect reduced confidence.

---

## `PNS_with_nucleosome_peak_calling.py`

### 1. Running on Individual Contigs

To analyze individual contigs, use the following command:

```bash
python3 PNS_with_nucleosome_peak_calling.py -b <your_bam_file.bam> -c <chromosome> --mode-length <int> --frag-lower <int> --frag-upper <int>
```

### 2. Output Files

#### a. Combined Scores (`bedGraph`) File Format

This file contains the scores for each contig and is formatted as follows:

- **Columns**:
  1. `chrom` (e.g., `chr1`)
  2. `start` (e.g., `1000`)
  3. `end` (e.g., `1001`)
  4. `coverage` (e.g., `50`)
  5. `PNS_smoothed` (e.g., `3.75`)
  6. `PNS` (e.g., `3.8`)
  7. `Dyad count` (e.g., `11`)

To convert to a bedGraph with a single PNS score per base:
```bash
awk 'BEGIN { OFS = "\t" } { print $1, $2, $3, $6 }' combined_scores.bedGraph > PNS.bedGraph
```

To convert to bigWig:
```bash
bedGraphToBigWig PNS.bedGraph chrom.sizes PNS.bw
```

chrom.sizes file can be generated from the BAM used to generate the scores bedGraph, or using the bedGraph as below:

```bash
awk '{if ($2+0 > max[$1]) max[$1] = $3} END {for (c in max) print c, max[c]}' PNS.bedGraph > chrom.sizes
```

#### b. Nucleosome Regions (`bed-like`) File Format

This file contains nucleosome regions identified in the analysis. The format includes:

- **Columns**:
  1. `contig`
  2. `nucleosome_region_start`
  3. `nucleosome_region_end`
  4. `prominence_score`
  5. `center_position`
  6. `upstream_negative_peak_score`
  7. `upstream_negative_peak_pos`
  8. `downstream_negative_peak_score`
  9. `downstream_negative_peak_pos`
  10. `positive_peak_score`
  11. `positive_peak_pos`
  12. `maximum_seq_coverage`
  13. `maximum_seq_coverage_pos`

To convert these nucleosome protection peak calls to standard BED format that can then be converted to a bigBed:

```bash
infile="nucleosome_regions.bed"
outfile="${infile%.bed}_converted.bed"

awk 'BEGIN { OFS="\t" }
{
    chr = $1
    start = $2
    end = $3
    name = chr ":" (start+1) "-" end
    summit = $5
    score = int($4 + 0.5)  # manual round
    if (score < 0) score = -score
    if (score > 1000) score = 1000
    strand = "+"
    thickEnd = summit + 1

    print chr, start, end, name, score, strand, summit, thickEnd
}' "$infile" > "$outfile"
```

---

## Setting up the Conda/Mamba Environment

### 1. Installing Conda

First, download and install Miniconda. Use the following commands to download the installer and run it:

```bash
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the instructions to complete the installation and initialize Miniconda

# Once the installation is complete, activate the changes to your shell
source ~/.bashrc
```

### 2. Installing Mamba via Conda

Mamba is a fast, efficient package manager that enhances Conda's functionality. It uses Conda as its underlying system but significantly speeds up package installation, dependency resolution, and environment creation.

Install Mamba using the following Conda command:

```bash
# Install mamba in the base environment
conda install mamba -n base -c conda-forge
```

### 3. Creating the `PNS` Environment

Once Mamba is installed, you can create a new environment using the `PNS_with_dependencies.yaml` file, which contains all the necessary dependencies for this project.

```bash
# Create the PNS environment from the YAML file
mamba env create -f PNS_with_dependencies.yaml
```

### 4. Activating the `PNS` Environment

After creating the environment, activate it with the following command:

```bash
# Activate the PNS environment
conda activate PNS_env
```

### 5. Verify Installation

To verify that all dependencies are installed correctly, you can list all the installed packages in the `PNS` environment:

```bash
conda list
```

This will display all the packages and their versions in your environment.

### 6. Deactivating the Environment

When you are done working, you can deactivate the environment by running:

```bash
conda deactivate
```
