# wavelength_cosine_mixture.py

This script generates an aggregate cosine signal from
a distribution of wavelengths.

It simulates how a spectrum of wavelengths interferes when summed, producing:

-   an averaged cosine waveform
-   a positive-peak envelope
-   a wavelength distribution (probabilities + integerized counts)
-   a heatmap showing each individual cosine instance sorted by
    wavelength

------------------------------------------------------------------------

## Core idea

1.  Choose a wavelength distribution (skewed Cauchy or Gaussian).

2.  Evaluate the distribution on integer wavelengths.

3.  Convert probabilities into integer counts that sum to
    `--total-counts`.

4.  For each wavelength w, generate:

    cos( 2π \* (\|x\| - offset) / w )

5.  Sum all cosine components and divide by total count.

6.  Extract an envelope by joining positive local maxima.

7.  Render each wave instance as a heatmap (downsampled if too large).

------------------------------------------------------------------------

## Outputs

All outputs share a common prefix.

### TSV files

-   \*\_wavelength_distribution.tsv\
    wavelength + probability

-   \*\_wavelength_counts.tsv\
    wavelength + count

-   \*\_avg_coswave.tsv\
    x + averaged cosine

-   \*\_pospeak_envelope.tsv\
    x + envelope

### Plots

-   \*\_avg_coswave.svg\
    Average cosine + envelope

-   \*\_instance_waves_heatmap.png\
    One row per cosine instance (sorted by wavelength)

------------------------------------------------------------------------

## Example

Cauchy:

python wavelength_cosine_mixture.py --a 0 --mode 192 --scale 12 --xmin 1 --xmax 2000 \
--total-counts 41465 --heatmap-max-rows 10000

Gaussian:

python wavelength_cosine_mixture.py --dist gaussian --mean 182 --sigma 8.5 --xmin 1 \
--xmax 2000 --peak-offset 130 --total-counts 41465 --heatmap-max-rows 10000

------------------------------------------------------------------------

## Command-line arguments

### Distribution selection

-   `--dist {skewcauchy,gaussian}`\
    Which wavelength distribution to use.

### Skewed Cauchy parameters

-   `--mode`\
    Location of the peak (acts as the mode / center).

-   `--a`\
    Skew parameter.\
    Positive = right-skewed, negative = left-skewed, zero = symmetric.

-   `--scale`\
    Width of the distribution (\>0). Larger values give broader
    wavelength spread.

### Gaussian parameters

-   `--mean`\
    Mean wavelength (overrides `--mode` if provided).

-   `--sigma`\
    Standard deviation (\>0). Controls spread.

### Wavelength grid

-   `--xmin`\
    Minimum integer wavelength considered.

-   `--xmax`\
    Maximum integer wavelength considered.

-   `--total-counts`\
    Total number of cosine instances generated (after discretization).

### Wave geometry

-   `--wave-xmin`\
    Minimum x position for cosine evaluation.

-   `--wave-xmax`\
    Maximum x position for cosine evaluation.

-   `--step`\
    Step size for x (default 1).

-   `--peak-offset`\
    Forces all cosine maxima to align at ±offset via `( |x| - offset )`.

### Envelope

-   `--env-min-distance`\
    Minimum spacing (in x samples) between detected envelope peaks.

-   `--env-min-height`\
    Minimum peak height for envelope detection (default 0).

### Heatmap

-   `--heatmap-max-rows`\
    Maximum number of rows in instance heatmap.\
    If `total-counts` exceeds this, rows are downsampled.

-   `--heatmap-seed`\
    RNG seed for reproducible heatmap downsampling.

### Output

-   `--outprefix`\
    Prefix for all output files.\
    If omitted, a prefix is automatically generated from parameters.

------------------------------------------------------------------------

## Conceptual framing

This script performs:

**wavelength spectrum → cosine basis → interference pattern**

Similar in spirit to:

-   Fourier mixtures
-   nucleosome phasing superposition
-   cfDNA fragment periodicity interference

