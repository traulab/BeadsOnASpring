#!/usr/bin/env python3
"""
@author Andrew D Johnston

This script builds a *mixture of cosine waves* whose wavelengths are drawn from a
user-defined wavelength distribution. It was designed for cases where you want a
distribution of wavelengths (e.g. centered near 186) and then want to see what the
aggregate “interference pattern” looks like when those cosine components are summed.

Core idea
---------
1) Choose a wavelength distribution on integer wavelengths [xmin..xmax]
   - skewcauchy: Azzalini skew-Cauchy (shape parameter a sets skew)
   - gaussian:   Normal distribution (mean, sigma)
2) Convert the distribution on the integer grid into:
   (a) a discrete probability distribution P(wavelength)
   (b) integer counts C(wavelength) that sum to total_counts
3) Generate a mirrored-offset cosine wave for each wavelength:
       cos( 2π * (|x| - peak_offset) / wavelength )
   Using |x| makes the pattern symmetric around 0 and causes all cosine peaks to
   align at x = +peak_offset and x = -peak_offset (when peak_offset > 0).
4) Compute:
   - the *average* aggregate signal:
         avg(x) = ( Σ_w C(w) * cos(...) ) / total_counts
   - a simple “positive-peak envelope” by finding local maxima above 0 and
     linearly connecting them.
5) Output:
   - wavelength_distribution.tsv  : wavelength <tab> probability (no header)
   - wavelength_counts.tsv        : wavelength <tab> integer count (no header; count>0 only)
   - avg_coswave.tsv              : x <tab> avg(x) (no header)
   - pospeak_envelope.tsv         : x <tab> envelope(x) (no header; only defined x)
   - avg_coswave.svg              : plot of avg(x) + envelope, with y fixed to [-1, 1]
   - instance_waves_heatmap.png   : heatmap where each row is one wave instance
                                   (one row per “copy” implied by counts), sorted by
                                   wavelength low→high. If total_counts is huge, the
                                   heatmap is downsampled to heatmap_max_rows.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import find_peaks


# -----------------------------
# Distribution definitions
# -----------------------------
def skew_cauchy_pdf(x: np.ndarray, a: float, loc: float, scale: float) -> np.ndarray:
    """
    Azzalini skew-Cauchy PDF evaluated at x.

    z = (x - loc) / scale
    PDF = (2/(pi*scale)) * 1/(1+z^2) * Phi(a*z)

    Parameters
    ----------
    x : array-like
        Points to evaluate the PDF at.
    a : float
        Skew parameter. a>0 gives right-skew, a<0 gives left-skew. a=0 is symmetric Cauchy.
    loc : float
        Location parameter (controls where the peak sits; you are using this as "mode").
    scale : float
        Width parameter (>0). Larger = broader.
    """
    z = (x - loc) / scale
    return (2.0 / (np.pi * scale)) * (1.0 / (1.0 + z**2)) * norm.cdf(a * z)


def gaussian_pdf(x: np.ndarray, mean: float, sigma: float) -> np.ndarray:
    """Gaussian (normal) PDF evaluated at x."""
    return norm.pdf(x, loc=mean, scale=sigma)


# -----------------------------
# Discretization: PDF -> integer counts
# -----------------------------
def counts_from_pdf(pdf: np.ndarray, total_counts: int) -> np.ndarray:
    """
    Convert a nonnegative "pdf-on-grid" array into integer counts that sum to total_counts.

    Method:
    - Normalize pdf to probabilities
    - expected = probs * total_counts
    - counts = floor(expected)
    - distribute the leftover remainder to the largest fractional parts

    This creates a deterministic integerization with minimal distortion.

    Returns
    -------
    counts : np.ndarray[int]
        Same length as pdf, sums to total_counts.
    """
    pdf = np.asarray(pdf, dtype=float)
    pdf[pdf < 0] = 0.0

    s = pdf.sum()
    if s <= 0:
        raise RuntimeError("PDF sum is zero across the requested range.")
    probs = pdf / s
    expected = probs * float(total_counts)

    counts = np.floor(expected).astype(int)
    remainder = total_counts - counts.sum()

    if remainder > 0:
        frac = expected - counts
        order = np.argsort(frac)[::-1]        # biggest fractional parts first
        counts[order[:remainder]] += 1
    elif remainder < 0:
        # Unlikely, but keep safe
        frac = expected - counts
        order = np.argsort(frac)             # smallest fractional parts first
        counts[order[:(-remainder)]] -= 1

    if counts.sum() != total_counts:
        raise RuntimeError("Internal rounding error: counts do not sum to total_counts.")
    return counts


# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate a wavelength distribution (skewCauchy or Gaussian), discretize to counts, "
            "sum mirrored-offset cosine waves, output SVG+TSVs, positive-peak envelope, and a "
            "PNG heatmap with one row per wave instance (optionally downsampled)."
        )
    )

    # Distribution choice
    p.add_argument(
        "--dist",
        choices=["skewcauchy", "gaussian"],
        default="skewcauchy",
        help="Distribution to use for wavelength model (default: skewcauchy)",
    )

    # Location-ish params
    p.add_argument(
        "--mode",
        type=float,
        default=186.0,
        help="For skewcauchy: loc (peak location). For gaussian: used as mean if --mean not set.",
    )
    p.add_argument(
        "--mean",
        type=float,
        default=None,
        help="Gaussian mean (overrides --mode for gaussian).",
    )

    # SkewCauchy params
    p.add_argument(
        "--a",
        type=float,
        default=5.0,
        help="Skew parameter for skewcauchy (a>0 right-skew, a<0 left-skew).",
    )
    p.add_argument(
        "--scale",
        type=float,
        default=50.0,
        help="Scale (>0) for skewcauchy (width).",
    )

    # Gaussian params
    p.add_argument(
        "--sigma",
        type=float,
        default=20.0,
        help="Gaussian std dev (>0).",
    )

    # Wavelength grid and total instances
    p.add_argument("--xmin", type=int, default=1, help="Min integer wavelength.")
    p.add_argument("--xmax", type=int, default=2000, help="Max integer wavelength.")
    p.add_argument(
        "--total-counts",
        type=int,
        default=100000,
        help="Total wave instances implied by discretized counts (sum of counts).",
    )

    # x-axis for waveform generation
    p.add_argument("--wave-xmin", type=int, default=-2500, help="Waveform x min.")
    p.add_argument("--wave-xmax", type=int, default=2500, help="Waveform x max.")
    p.add_argument("--step", type=int, default=1, help="Waveform x step.")

    # Peak alignment (mirrored about 0)
    p.add_argument(
        "--peak-offset",
        type=float,
        default=0.0,
        help=(
            "Mirror-symmetric shift so cosine maxima align at x=+offset and x=-offset. "
            "Implemented via (|x| - offset)."
        ),
    )

    # Envelope controls
    p.add_argument(
        "--env-min-distance",
        type=int,
        default=10,
        help="Minimum distance (in x samples) between envelope peaks.",
    )
    p.add_argument(
        "--env-min-height",
        type=float,
        default=0.0,
        help="Minimum peak height for envelope (default 0.0 → peaks above 0).",
    )

    # Heatmap controls
    p.add_argument(
        "--heatmap-max-rows",
        type=int,
        default=20000,
        help=(
            "Max rows to render in instance heatmap. If total_counts exceeds this, heatmap rows "
            "are downsampled by sampling wavelengths proportional to counts."
        ),
    )
    p.add_argument(
        "--heatmap-seed",
        type=int,
        default=1,
        help="RNG seed for heatmap downsampling.",
    )

    # Output naming
    p.add_argument(
        "--outprefix",
        default=None,
        help="Output prefix (default derived from parameters).",
    )
    return p.parse_args()


def main():
    args = parse_args()

    # -----------------------------
    # Basic argument validation
    # -----------------------------
    if args.xmin > args.xmax:
        raise ValueError("--xmin must be <= --xmax")
    if args.total_counts <= 0:
        raise ValueError("--total-counts must be > 0")
    if args.step <= 0:
        raise ValueError("--step must be > 0")
    if args.peak_offset < 0:
        raise ValueError("--peak-offset must be >= 0")
    if args.heatmap_max_rows <= 0:
        raise ValueError("--heatmap-max-rows must be > 0")

    # -----------------------------
    # 1) Build wavelength grid (integers)
    # -----------------------------
    wl = np.arange(args.xmin, args.xmax + 1, dtype=int)   # integer wavelengths
    wl_f = wl.astype(float)                               # float for PDF evaluation

    # -----------------------------
    # 2) Evaluate chosen distribution on the wavelength grid
    # -----------------------------
    if args.dist == "skewcauchy":
        if args.scale <= 0:
            raise ValueError("--scale must be > 0")
        pdf = skew_cauchy_pdf(wl_f, a=args.a, loc=args.mode, scale=args.scale)
        title_dist = f"skewCauchy a={args.a:g}, mode={args.mode:g}, scale={args.scale:g}"
        prefix = args.outprefix or (
            f"skewcauchy_a{args.a:g}_mode{args.mode:g}_scale{args.scale:g}_"
            f"xmin{args.xmin}_xmax{args.xmax}_N{args.total_counts}_off{args.peak_offset:g}"
        )
    else:
        if args.sigma <= 0:
            raise ValueError("--sigma must be > 0")
        mean = args.mean if args.mean is not None else args.mode
        pdf = gaussian_pdf(wl_f, mean=mean, sigma=args.sigma)
        title_dist = f"Gaussian mean={mean:g}, sigma={args.sigma:g}"
        prefix = args.outprefix or (
            f"gaussian_mean{mean:g}_sigma{args.sigma:g}_"
            f"xmin{args.xmin}_xmax{args.xmax}_N{args.total_counts}_off{args.peak_offset:g}"
        )

    # Ensure nonnegative and normalize to discrete probabilities over integer wl
    pdf_nonneg = np.clip(pdf, 0, None)
    if pdf_nonneg.sum() <= 0:
        raise RuntimeError("PDF evaluated to zero across wavelength range; check params/range.")
    probs = pdf_nonneg / pdf_nonneg.sum()

    # -----------------------------
    # 3) Convert probabilities to integer counts that sum to total_counts
    # -----------------------------
    counts = counts_from_pdf(pdf_nonneg, args.total_counts).astype(int)

    # -----------------------------
    # 4) Define outputs (no headers by design, Excel-friendly)
    # -----------------------------
    out_dist_tsv     = f"{prefix}_wavelength_distribution.tsv"   # wavelength\tprobability
    out_counts_tsv   = f"{prefix}_wavelength_counts.tsv"         # wavelength\tcount (only count>0)
    out_wave_tsv     = f"{prefix}_avg_coswave.tsv"               # x\tavg_signal
    out_env_tsv      = f"{prefix}_pospeak_envelope.tsv"          # x\tenvelope (only defined region)
    out_svg          = f"{prefix}_avg_coswave.svg"               # average + envelope plot
    out_heatmap_png  = f"{prefix}_instance_waves_heatmap.png"    # instance heatmap

    # -----------------------------
    # 5) Write wavelength distribution TSV (probabilities)
    # -----------------------------
    with open(out_dist_tsv, "w") as f:
        for w, p in zip(wl, probs):
            if p > 0:
                f.write(f"{w}\t{p}\n")

    # -----------------------------
    # 6) Write wavelength counts TSV (integer counts)
    # -----------------------------
    with open(out_counts_tsv, "w") as f:
        for w, c in zip(wl, counts):
            if c > 0:
                f.write(f"{w}\t{c}\n")

    # -----------------------------
    # 7) Build the aggregate average cosine waveform
    # -----------------------------
    x = np.arange(args.wave_xmin, args.wave_xmax + 1, args.step, dtype=float)

    # Mirrored shift: x_eff = |x| - peak_offset
    # Ensures maxima at x=+offset and x=-offset for all wavelengths
    x_eff = np.abs(x) - float(args.peak_offset)

    # Weighted sum of cosines; then divide by total_counts to get an average
    summed = np.zeros_like(x, dtype=float)
    nz = np.nonzero(counts)[0]
    for idx in nz:
        w = float(wl[idx])
        c = float(counts[idx])
        summed += c * np.cos(2.0 * np.pi * x_eff / w)

    avg = summed / float(args.total_counts)

    # Write average waveform TSV
    with open(out_wave_tsv, "w") as f:
        for xi, yi in zip(x.astype(int), avg):
            f.write(f"{xi}\t{yi}\n")

    # -----------------------------
    # 8) Positive-peak envelope (connect peaks above 0)
    # -----------------------------
    # Find local maxima in avg() that exceed env_min_height (default 0.0) and are spaced apart
    peaks, _ = find_peaks(avg, height=args.env_min_height, distance=args.env_min_distance)

    # Envelope will be NaN except between first and last detected peak,
    # where we linearly interpolate between successive peak points
    env = np.full_like(avg, np.nan, dtype=float)
    if peaks.size >= 2:
        xp = x[peaks]
        yp = avg[peaks]
        mask = (x >= xp[0]) & (x <= xp[-1])
        env[mask] = np.interp(x[mask], xp, yp)

    # Write envelope TSV (only finite)
    with open(out_env_tsv, "w") as f:
        for xi, yi in zip(x.astype(int), env):
            if np.isfinite(yi):
                f.write(f"{xi}\t{yi}\n")

    # Plot average waveform + envelope as SVG; fix y-axis to [-1, 1]
    plt.figure(figsize=(12, 5))
    plt.plot(x, avg, lw=1.5, label="average cosine")
    if peaks.size >= 2:
        plt.plot(x, env, lw=1.2, label="envelope (pos peaks)")
    plt.xlabel("x")
    plt.ylabel("Average signal")
    plt.title(f"Average cosine + envelope ({title_dist}; peak offset ±{args.peak_offset:g})")
    plt.xlim(args.wave_xmin, args.wave_xmax)
    plt.ylim(-1, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_svg)
    plt.close()

    # -----------------------------
    # 9) Instance heatmap (one row per instance), sorted by wavelength
    # -----------------------------
    # WARNING: total_counts rows can be huge, so we downsample rows if needed.
    total_instances = int(counts.sum())  # should equal args.total_counts
    max_rows = int(args.heatmap_max_rows)
    rng = np.random.default_rng(args.heatmap_seed)

    if total_instances <= max_rows:
        # Exact expansion: one row per instance
        wl_instances = np.repeat(wl, counts)
        downsampled = False
    else:
        # Downsample: sample wavelengths with probability proportional to counts
        p_w = counts.astype(float) / counts.sum()
        wl_instances = rng.choice(wl, size=max_rows, replace=True, p=p_w)
        downsampled = True

    # Sort so the heatmap is ordered low→high wavelength
    wl_instances.sort()
    used_rows = wl_instances.size

    # Build heatmap matrix (rows=instances, cols=x)
    # float32 reduces memory (and is plenty for cos values)
    heat = np.empty((used_rows, x.size), dtype=np.float32)
    for i, w in enumerate(wl_instances):
        heat[i, :] = np.cos(2.0 * np.pi * x_eff / float(w))

    # Plot heatmap using the same colormap used in your nucleosome heatmap script: "seismic"
    plt.figure(figsize=(12, 8))
    im = plt.imshow(
        heat,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
        cmap="seismic",     # red-white-blue diverging (same as your NPS heatmap script)
        vmin=-1,
        vmax=1,
        extent=[x[0], x[-1], 0, used_rows],
    )
    cbar = plt.colorbar(im)
    cbar.set_label("Cosine amplitude", fontsize=12)

    plt.xlabel("x")
    plt.ylabel("Wave instance index (sorted by wavelength)")
    title = "Instance cosine waves heatmap"
    if downsampled:
        title += f" (downsampled to {used_rows} rows from {total_instances})"
    else:
        title += f" ({used_rows} rows)"
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_heatmap_png, dpi=300)
    plt.close()

    # -----------------------------
    # 10) Summary
    # -----------------------------
    print("Wrote:")
    print(f"  {out_dist_tsv}        (wavelength\\tprobability; no header)")
    print(f"  {out_counts_tsv}      (wavelength\\tcount; no header; count>0 only)")
    print(f"  {out_wave_tsv}        (x\\tavg_signal; no header)")
    print(f"  {out_env_tsv}         (x\\tenvelope; no header; only defined region)")
    print(f"  {out_svg}")
    print(f"  {out_heatmap_png}     (rows={'downsampled' if downsampled else 'all'}; sorted by wavelength)")

if __name__ == "__main__":
    main()
