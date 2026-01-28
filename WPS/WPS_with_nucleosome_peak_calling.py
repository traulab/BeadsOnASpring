#!/usr/bin/env python3
"""
@author Andrew D Johnston

Kircher-style WPS scoring + median-centering + Kircher-closely-matched peak calling

This script computes per-base tracks from paired-end BAM fragments:
  - coverage: fragment overlap depth (0-based, half-open)
  - dyad: fragment midpoint counts
  - wps: Kircher-equivalent WPS kernel sum
  - wps_smoothed: Savitzky–Golay smoothed WPS
  - mWPS: median-centered WPS (wps - rolling_median(wps))
  - sm_mWPS: Kircher "--smoother" analogue:
        sm_mWPS = wps_smoothed - rolling_median(wps)
    (i.e., smooth first, subtract the median of the RAW WPS window)

Peak calling (matched to Kircher 2015):
  - Operates on a per-bp track (sm_mWPS by default).
  - Builds positive runs where ivalue > 0, merging across gaps <= 5 bp (default),
    filling gaps with zeros.
  - For each merged run, applies Kircher evaluateValues logic:
      * minlength <= L <= maxlength: choose the contiguous >= median window with max SUM
      * maxlength <= L <= 3*maxlength: emit each contiguous >= median window
        whose length is within [minlength, maxlength]
      * otherwise: emit nothing
  - Uses Kircher's median(values) behavior for region median-thresholding (upper-middle).
  - Uses Kircher's continousWindows() behavior (including its "drop-first-after-gap").
  - Uses 1-based internal ipos coordinates for calls, and outputs BED like Kircher:
      chr{chrom}  (cstart-1)  cend  name=chrom:cstart-cend  score=cval  +  thickStart=cmiddle-1  thickEnd=cmiddle
  - Midpoint uses "round-half-up" to emulate Kircher/Python2 rounding.

Coordinate conventions:
  - Internal scoring arrays are indexed by 0-based genome coordinates.
  - Peak calling uses 1-based ipos to match Kircher’s output semantics.
  - BED output is 0-based start, end-exclusive, matching Kircher’s printed columns.
"""

import sys
import argparse
import os
import random
import math
from collections import defaultdict

import numpy as np
import pysam
from tqdm import tqdm
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import savgol_filter


# ----------------------------
# CIGAR filter (Kircher-style)
# ----------------------------
def is_softclipped_or_padded(cigartuples):
    """Return True if CIGAR has S/H/P ops (soft clip, hard clip, pad)."""
    if not cigartuples:
        return False
    for op, _ln in cigartuples:
        if op in (4, 5, 6):
            return True
    return False


# ----------------------------
# Paired-end fragment generator
# ----------------------------
def generate_paired_reads(bamfile, contig=None, start=None, end=None, max_duplicates=0, subsample=None):
    """
    Yield paired reads as (fwd, rev) with:
      - unmapped/mate-unmapped filtered
      - same-strand pairs filtered
      - duplicate/qcfail filtered
      - soft/hard clip/pad filtered
      - duplicate fragments limited to max_duplicates by (contig, frag_start, frag_end)
      - optional subsampling fraction in [0,1]
    """
    unpaired = {}
    frag_counts = defaultdict(int)

    try:
        it = bamfile.fetch(contig, start, end, multiple_iterators=True)
    except Exception:
        return

    for read in it:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_duplicate or read.is_qcfail:
            continue
        if is_softclipped_or_padded(read.cigartuples):
            continue
        if read.reference_end is None or read.next_reference_start is None:
            continue

        qn = read.query_name
        if qn not in unpaired:
            unpaired[qn] = read
            continue

        mate = unpaired.pop(qn)

        if mate.is_unmapped or mate.is_duplicate or mate.is_qcfail:
            continue
        if is_softclipped_or_padded(mate.cigartuples):
            continue

        if read.is_reverse == mate.is_reverse:
            continue

        if subsample is not None and random.random() > subsample:
            continue

        frag_start = min(read.reference_start, mate.reference_start)
        frag_end = max(read.reference_end, mate.reference_end)
        frag_contig = read.reference_name

        if frag_end <= frag_start:
            continue

        key = (frag_contig, frag_start, frag_end)
        if frag_counts[key] > max_duplicates:
            continue
        frag_counts[key] += 1

        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read


def generate_fragment_ranges(bamfile, contig, start, end, max_duplicates, subsample):
    """Yield fragments as (frag_start, frag_end) in 0-based half-open coords, requiring overlap with [start,end)."""
    for r_fwd, r_rev in generate_paired_reads(bamfile, contig, start, end, max_duplicates, subsample):
        frag_start = min(r_fwd.reference_start, r_rev.reference_start)
        frag_end = max(r_fwd.reference_end, r_rev.reference_end)
        if frag_end <= frag_start:
            continue
        if frag_end <= start or frag_start >= end:
            continue
        yield frag_start, frag_end


# ----------------------------
# Kircher-equivalent WPS kernel (effective bx semantics)
# ----------------------------
def wps_kernel_kircher_exact(L_true: int, protection: int = 120) -> np.ndarray:
    """
    Kernel matching Kircher's effective computation with bx intervals.

    For protection=120 => half=60, effective window length behaves like 120 (2*half),
    and the interval semantics make the effective overlap/span logic equivalent to:

      total_len = L_true + 2*half - 2
      flank     = 2*half - 1
      mid       = L_true - 2*half

    If mid <= 0 => all -1.
    """
    half = protection // 2
    if L_true <= 0:
        return np.array([], dtype=np.int8)

    total_len = L_true + 2 * half - 2
    if total_len <= 0:
        return np.array([], dtype=np.int8)

    flank = 2 * half - 1
    mid = L_true - 2 * half

    if mid <= 0:
        return np.full(total_len, -1, dtype=np.int8)

    k = np.empty(total_len, dtype=np.int8)
    k[:flank] = -1
    k[flank:flank + mid] = +1
    k[flank + mid:] = -1
    return k


def precompute_distributions_kircher_exact(wps_frag_range, protection: int = 120):
    """Precompute kernels for each TRUE fragment length."""
    return {int(L): wps_kernel_kircher_exact(int(L), protection=protection) for L in wps_frag_range}


# ----------------------------
# Rolling median baseline (window=1000)
# ----------------------------
def rolling_median(x: np.ndarray, window: int = 1000) -> np.ndarray:
    """
    Rolling median for even window sizes.
    Places the median at the RIGHT-middle index (half = window//2).
    Returns NaN where the full window doesn't fit.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < window:
        return np.full(n, np.nan, dtype=float)

    wins = sliding_window_view(x, window_shape=window)
    med = np.median(wins, axis=1)  # length n-window+1

    out = np.full(n, np.nan, dtype=float)
    half = window // 2
    out[half: n - half + 1] = med
    return out


# ----------------------------
# Scoring
# ----------------------------
def score_contig(
    bamfiles,
    contig,
    start,
    end,
    protection,
    wps_frag_range,
    max_duplicates,
    distributions,
    subsample,
    baseline_window=1000,
    sg_window=21,
    sg_order=2,
):
    """
    Compute tracks over [start,end):
      - coverage: TRUE fragment overlap
      - dyad: TRUE fragment midpoint
      - wps: Kircher-equivalent kernel sum
      - wps_smoothed: Savitzky–Golay on raw wps
      - mWPS: wps - rolling median of raw wps
      - sm_mWPS: wps_smoothed - rolling median of raw wps   (Kircher --smoother behavior)
    """
    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=np.int32)
    dyad = np.zeros(ref_len, dtype=np.int32)
    wps = np.zeros(ref_len, dtype=np.float64)

    half = protection // 2

    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(bamfile, contig, start, end, max_duplicates, subsample):
            L_true = frag_end - frag_start
            if L_true <= 0:
                continue

            cov_s = max(frag_start, start)
            cov_e = min(frag_end, end)
            if cov_e > cov_s:
                coverage[cov_s - start: cov_e - start] += 1

            frag_center = frag_start + (L_true - 1) // 2
            if start <= frag_center < end:
                dyad[frag_center - start] += 1

            if L_true not in wps_frag_range:
                continue

            kernel = distributions.get(L_true, None)
            if kernel is None or kernel.size == 0:
                continue

            # Effective placement matching Kircher’s bx overlap domain.
            kernel_start_genome = frag_start - half + 1

            k_s = kernel_start_genome - start
            k_e = k_s + kernel.size

            arr_s = max(k_s, 0)
            arr_e = min(k_e, ref_len)
            if arr_e <= arr_s:
                continue

            ker_s = arr_s - k_s
            ker_e = ker_s + (arr_e - arr_s)

            wps[arr_s:arr_e] += kernel[ker_s:ker_e].astype(np.float64)

    # Smooth raw WPS (Kircher smoother uses a SG-like smoother on raw WPS values)
    if ref_len >= sg_window and sg_window % 2 == 1:
        wps_smoothed = savgol_filter(wps, sg_window, sg_order)
    else:
        wps_smoothed = wps.copy()

    # Baseline is the rolling median of RAW WPS
    baseline = rolling_median(wps, window=baseline_window)
    if np.isnan(baseline).any():
        valid = np.where(~np.isnan(baseline))[0]
        if valid.size > 0:
            first, last = valid[0], valid[-1]
            baseline[:first] = baseline[first]
            baseline[last + 1:] = baseline[last]
        else:
            baseline[:] = 0.0

    mWPS = wps - baseline
    sm_mWPS = wps_smoothed - baseline

    return {
        "coverage": [(contig, start, coverage)],
        "dyad": [(contig, start, dyad)],
        "wps": [(contig, start, wps)],
        "wps_smoothed": [(contig, start, wps_smoothed)],
        "mWPS": [(contig, start, mWPS)],
        "sm_mWPS": [(contig, start, sm_mWPS)],
    }


# ----------------------------
# Kircher-matched peak calling
# ----------------------------
def kircher_median(values):
    """
    Match Kircher's median(values) (NOT numpy median):
      helper.sort()
      if round(lVal*0.5) == int(lVal*0.5): return helper[int(lVal*0.5)]
      else: average of the two indexed values

    For even length, this returns the UPPER middle element.
    """
    helper = list(values)
    helper.sort()
    lVal = len(helper)
    point = 0.5
    if round(lVal * point) == int(lVal * point):
        return float(helper[int(lVal * point)])
    else:
        return float((helper[int(round(lVal * point))] + helper[int(lVal * point)]) * 0.5)


def kircher_continous_windows(region_pairs):
    """
    Faithful to Kircher's continousWindows(region):
    when a gap occurs, it appends the current window and RESETS state,
    but does NOT start a new window with the current (pos,val) point.
    """
    res = []
    cstart, cend = None, None
    cmax = None
    csum = 0.0
    for pos, val in region_pairs:
        if cmax is None:
            cend = pos
            cstart = pos
            cmax = val
            csum = val
        else:
            if cend + 1 == pos:
                cend = pos
                csum += val
                if cmax < val:
                    cmax = val
            else:
                res.append((csum, cstart, cend, cmax))
                cstart, cend = None, None
                cmax = None
                csum = 0.0
                # Intentionally do NOT initialize new window with (pos,val).
    if cmax is not None:
        res.append((csum, cstart, cend, cmax))
    return res


def round_half_up(x):
    """Emulate Python2-style rounding for .5 (round-half-up)."""
    return int(math.floor(float(x) + 0.5))


def evaluate_values_kircher(chrom_nochr, start1, end1, values, report=True, minlength=50, maxlength=150, vari_cutoff=5.0):
    """
    Closely match Kircher evaluateValues(chrom,start,end,values,report).

    Inputs:
      - chrom_nochr: "1","2",...,"X","Y"
      - start1/end1: 1-based inclusive positions (ipos semantics)
      - values: list of ivalue values aligned to positions start1..end1
    Returns:
      list of BED rows tuples:
        (chrom_with_chr, bed_start0, bed_end0, name, score_int, "+", thickStart0, thickEnd0)
    """
    if start1 is None or end1 is None or not values:
        return []

    L = len(values)
    calls = []

    if maxlength >= L >= minlength:
        cMed = kircher_median(values)
        region_pairs = [(pos, val) for pos, val in zip(range(start1, end1 + 1), values) if val >= cMed]
        windows = kircher_continous_windows(region_pairs)
        if not windows:
            return []
        windows.sort()
        _score_sum, cstart, cend, cval = windows[-1]
        if report and (cval > vari_cutoff):
            cmiddle = cstart + round_half_up((cend - cstart) * 0.5)
            chrom_out = f"chr{chrom_nochr}"
            bed_start = cstart - 1
            bed_end = cend
            name = f"{chrom_nochr}:{cstart}-{cend}"
            score_int = int(round(float(cval)))
            thick_start = cmiddle - 1
            thick_end = cmiddle
            calls.append((chrom_out, bed_start, bed_end, name, score_int, "+", thick_start, thick_end))
        return calls

    elif (3 * maxlength) >= L >= maxlength:
        cMed = kircher_median(values)
        region_pairs = [(pos, val) for pos, val in zip(range(start1, end1 + 1), values) if val >= cMed]
        windows = kircher_continous_windows(region_pairs)
        for _score_sum, cstart, cend, cval in windows:
            seg_len = (cend - cstart + 1)
            if maxlength >= seg_len >= minlength:
                if report and (cval > vari_cutoff):
                    cmiddle = cstart + round_half_up((cend - cstart) * 0.5)
                    chrom_out = f"chr{chrom_nochr}"
                    bed_start = cstart - 1
                    bed_end = cend
                    name = f"{chrom_nochr}:{cstart}-{cend}"
                    score_int = int(round(float(cval)))
                    thick_start = cmiddle - 1
                    thick_end = cmiddle
                    calls.append((chrom_out, bed_start, bed_end, name, score_int, "+", thick_start, thick_end))
        return calls

    return []


def call_peaks_kircher_matched(
    contig,
    adjusted_start0,
    track,
    merge_gap_bp=5,
    minlength=50,
    maxlength=150,
    vari_cutoff=5.0,
):
    """
    Kircher-matched positive-run builder + evaluateValues.

    Uses ipos (1-based) positions:
      ipos = adjusted_start0 + i + 1

    Merge rule:
      merge if ipos <= cend + merge_gap_bp  (step=1)
      fill missing positions with zeros.
    """
    # Chrom name for filtering (Kircher uses "1".. "22","X","Y")
    chrom_nochr = contig.replace("chr", "")
    allowed = {str(i) for i in range(1, 23)} | {"X", "Y"}
    report = chrom_nochr in allowed

    calls = []
    cstart = None   # 1-based inclusive
    cend = None     # 1-based inclusive
    clist = []

    n = int(track.size)
    for i in range(n):
        ivalue = float(track[i])
        ipos = int(adjusted_start0 + i + 1)  # 1-based

        if ivalue > 0:
            if (cend is not None) and (ipos <= (cend + merge_gap_bp)):
                while (cend + 1) < ipos:
                    cend += 1
                    clist.append(0.0)
                clist.append(ivalue)
                cend = ipos
            else:
                if cstart is not None and clist and report:
                    calls.extend(
                        evaluate_values_kircher(
                            chrom_nochr, cstart, cend, clist,
                            report=report,
                            minlength=minlength,
                            maxlength=maxlength,
                            vari_cutoff=vari_cutoff,
                        )
                    )
                clist = [ivalue]
                cstart = ipos
                cend = ipos

    if cstart is not None and clist and report:
        calls.extend(
            evaluate_values_kircher(
                chrom_nochr, cstart, cend, clist,
                report=report,
                minlength=minlength,
                maxlength=maxlength,
                vari_cutoff=vari_cutoff,
            )
        )

    return calls


# ----------------------------
# Output writers
# ----------------------------
def write_bedgraph(scores, contigs, out_prefix, first_region=False):
    mode = "w" if first_region else "a"
    original_start, original_end = contigs[0]
    filename = f"{out_prefix}_combined_scores.bedGraph"

    with open(filename, mode) as f:
        first_score_type = list(scores.keys())[0]
        contig_data = scores[first_score_type]

        for contig, start, score_array in contig_data:
            contig_out = contig if contig.startswith("chr") else f"chr{contig}"

            for i in range(len(score_array)):
                position = start + i
                if original_start <= position < original_end:
                    line = [contig_out, str(position), str(position + 1)]
                    for stype in scores.keys():
                        arr = scores[stype][0][2]
                        val = arr[i]
                        if isinstance(val, (np.integer, int)) or (isinstance(val, (np.floating, float)) and float(val).is_integer()):
                            line.append(str(int(val)))
                        else:
                            line.append(f"{float(val):.2f}")
                    f.write("\t".join(line) + "\n")


def write_bed_rows(rows, path, mode):
    with open(path, mode) as f:
        for chrom, start, end, name, score, strand, thick_start, thick_end in rows:
            f.write(
                f"{chrom}\t{int(start)}\t{int(end)}\t{name}\t{int(score)}\t{strand}\t{int(thick_start)}\t{int(thick_end)}\n"
            )


def split_into_regions(contig, start, end, contig_len, max_length=100000, overlap=1000):
    """
    Split [start,end) into core chunks of length <= max_length.
    Each chunk is scored on an expanded interval with +/- overlap padding,
    clipped to [0, contig_len).

    Returns tuples:
      (contig, adjusted_start, adjusted_end, original_start, original_end)
    """
    regions = []
    current_start = start
    while current_start < end:
        original_start = current_start
        original_end = min(current_start + max_length, end)

        adjusted_start = max(original_start - overlap, 0)
        adjusted_end = min(original_end + overlap, contig_len)

        regions.append((contig, adjusted_start, adjusted_end, original_start, original_end))
        current_start = original_end
    return regions


# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser(description="Kircher-style WPS scoring + median-centering + Kircher-matched peak calling.")
    parser.add_argument("-b", "--bamfiles", nargs="+", required=True, help="BAM file(s) to process")
    parser.add_argument("-o", "--out_prefix", default=None, help="Output prefix")
    parser.add_argument(
        "-c", "--contigs", nargs="+", default=None,
        help='Limit to contig(s) and optional range, e.g. "12:51730340-52039340" or "12"'
    )
    parser.add_argument("--protection", type=int, default=120, help="Protection window (bp), default 120.")
    parser.add_argument("--frag-lower", type=int, default=127, help="Lower fragment size to include in WPS")
    parser.add_argument("--frag-upper", type=int, default=207, help="Upper fragment size to include in WPS")
    parser.add_argument("--max-duplicates", type=int, default=0, help="Maximum allowed duplicate fragments (same coords)")
    parser.add_argument("--subsample", type=float, default=None, help="Subsampling proportion (e.g. 0.5 keeps ~50%)")
    parser.add_argument("--chunk-bp", type=int, default=100000, help="Chunk size per contig")
    parser.add_argument("--overlap-bp", type=int, default=1000, help="Overlap padding for edge-safe scoring")
    parser.add_argument("--baseline-window", type=int, default=1000, help="Rolling median window for baseline subtraction")
    parser.add_argument("--sg-window", type=int, default=21, help="Savitzky-Golay window (odd)")
    parser.add_argument("--sg-order", type=int, default=2, help="Savitzky-Golay polynomial order")

    # Peak calling controls (Kircher defaults)
    parser.add_argument("--peak-minlen", type=int, default=50, help="Minimum length (bp) for candidate windows")
    parser.add_argument("--peak-maxlen", type=int, default=150, help="Maximum length (bp) for reported windows")
    parser.add_argument("--peak-maxregion", type=int, default=450, help="Reject merged regions longer than this (bp)")
    parser.add_argument("--peak-merge-gap", type=int, default=5, help="Merge positive runs if gap <= this many bp")
    parser.add_argument("--peak-varicutoff", type=float, default=5.0, help="Minimum max score to report a peak window")

    args = parser.parse_args()

    bamfiles = []
    for p in args.bamfiles:
        try:
            bamfiles.append(pysam.AlignmentFile(p, "rb"))
        except Exception as e:
            parser.error(f"Unable to open BAM {p}: {e}")
            return 2

    if not args.out_prefix:
        bnames = [os.path.splitext(os.path.basename(b))[0] for b in args.bamfiles]
        args.out_prefix = "_".join(bnames)
        if args.contigs and len(args.contigs) == 1:
            args.out_prefix = f"{args.out_prefix}_{args.contigs[0]}"

    args.out_prefix = (
        f"{args.out_prefix}_prot{args.protection}"
        f"_lower{args.frag_lower}_upper{args.frag_upper}"
        f"_maxdup{args.max_duplicates}"
    )

    contigs = []
    if args.contigs:
        for contig_range in args.contigs:
            if ':' in contig_range:
                contig, positions = contig_range.split(':')
                start, end = map(int, positions.split('-'))
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(split_into_regions(
                    contig, start, end, contig_len,
                    max_length=args.chunk_bp,
                    overlap=args.overlap_bp
                ))
            else:
                contig = contig_range
                start, end = 0, bamfiles[0].get_reference_length(contig)
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(split_into_regions(
                    contig, start, end, contig_len,
                    max_length=args.chunk_bp,
                    overlap=args.overlap_bp
                ))
    else:
        for contig in bamfiles[0].references:
            start, end = 0, bamfiles[0].get_reference_length(contig)
            contig_len = bamfiles[0].get_reference_length(contig)
            contigs.extend(split_into_regions(
                contig, start, end, contig_len,
                max_length=args.chunk_bp,
                overlap=args.overlap_bp
            ))

    wps_frag_range = set(range(args.frag_lower, args.frag_upper + 1))
    distributions = precompute_distributions_kircher_exact(wps_frag_range, protection=args.protection)

    combined_bedgraph = f"{args.out_prefix}_combined_scores.bedGraph"
    nuc_bed = f"{args.out_prefix}_nucleosome_regions.bed"
    brk_bed = f"{args.out_prefix}_breakpoint_peaks.bed"
    for fn in (combined_bedgraph, nuc_bed, brk_bed):
        if os.path.exists(fn):
            os.remove(fn)

    first_region = True

    for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc="Scoring contigs"):
        scores = score_contig(
            bamfiles=bamfiles,
            contig=contig,
            start=adjusted_start,
            end=adjusted_end,
            protection=args.protection,
            wps_frag_range=wps_frag_range,
            max_duplicates=args.max_duplicates,
            distributions=distributions,
            subsample=args.subsample,
            baseline_window=args.baseline_window,
            sg_window=args.sg_window,
            sg_order=args.sg_order,
        )

        # Kircher --smoother analogue:
        #   ivalue = smoother(WPS) - median(WPS window)
        track = scores["sm_mWPS"][0][2]

        nuc_rows = call_peaks_kircher_matched(
            contig=contig,
            adjusted_start0=adjusted_start,
            track=track,
            merge_gap_bp=args.peak_merge_gap,
            minlength=args.peak_minlen,
            maxlength=args.peak_maxlen,
            vari_cutoff=args.peak_varicutoff,
        )
        brk_rows = call_peaks_kircher_matched(
            contig=contig,
            adjusted_start0=adjusted_start,
            track=(-1.0 * track),
            merge_gap_bp=args.peak_merge_gap,
            minlength=args.peak_minlen,
            maxlength=args.peak_maxlen,
            vari_cutoff=args.peak_varicutoff,
        )

        # Enforce maximum merged-region rejection (Kircher uses 3*maxlength)
        if args.peak_maxregion != 3 * args.peak_maxlen:
            nuc_rows = [r for r in nuc_rows if (r[2] - r[1]) <= args.peak_maxregion]
            brk_rows = [r for r in brk_rows if (r[2] - r[1]) <= args.peak_maxregion]

        # Keep only calls overlapping the core non-overlap interval
        def keep_core(rows):
            out = []
            for chrom, s, e, name, score, strand, ts, te in rows:
                # rows are BED: [s,e) 0-based
                if e <= original_start or s >= original_end:
                    continue
                out.append((chrom, s, e, name, score, strand, ts, te))
            return out

        nuc_rows = keep_core(nuc_rows)
        brk_rows = keep_core(brk_rows)

        write_bed_rows(nuc_rows, nuc_bed, mode=("w" if first_region else "a"))
        write_bed_rows(brk_rows, brk_bed, mode=("w" if first_region else "a"))

        write_bedgraph(scores, [(original_start, original_end)], args.out_prefix, first_region)
        first_region = False

    for b in bamfiles:
        try:
            b.close()
        except Exception:
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
