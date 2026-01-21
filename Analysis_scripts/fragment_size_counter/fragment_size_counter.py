#!/usr/bin/env python3
"""
@author Andrew D Johnston

Assign paired-end BAM fragments to chromatin-state BED intervals and compute fragment-length
distributions per state.

Key behavior
- Fragments are assigned by the *fragment midpoint* (integer midpoint between fragment start and end)
  into the BED interval that contains that midpoint coordinate.
- Processing is done in coordinate chunks per chromosome to limit memory.
- Outputs one TSV per state:  <bam>_<state>_<label>.txt  with:  length<TAB>count
- Optional overlay density plots in a user-defined length window.

BED assumptions (configurable via args)
- chrom: required
- start,end: required
- state label column: default col4 (1-based)
- RGB column: default col9 (1-based) formatted "R,G,B" (0–255), optional

Notes / constraints
- Pairing is performed *within each chunk only*. If mates land in different chunks, they may not pair.
  (This is acceptable for coordinate-sorted BAMs with typical insert sizes when chunk_size is large.)
- The BED intervals for a chromosome are assumed to be sorted and non-overlapping for the fast scan.
"""

import argparse
import os
import re
from collections import defaultdict

import pysam
from tqdm import tqdm

# ---------- Default plotting colors (matplotlib RGB tuples in 0..1) ----------
DEFAULT_EUCHROMATIN_COLOR = (0.5, 0.0, 0.5)  # purple
DEFAULT_MID_GREY = (0.5, 0.5, 0.5)           # override for selected prefixes
DEFAULT_WHOLE_GENOME_COLOR = (0.0, 0.0, 0.0) # black


# ----------------------------
# Small helpers
# ----------------------------
def sanitize_filename(text: str) -> str:
    """Make a safe filename component from a state label."""
    return re.sub(r'[\\/:*?"<>|]', "_", text)


def detect_bam_chr_style(bamfile: pysam.AlignmentFile) -> str:
    """Infer whether BAM contigs are like 'chr1' or '1'."""
    refs = list(bamfile.references)
    if "chr1" in refs:
        return "chr"
    if "1" in refs:
        return "nochr"
    chr_like = sum(1 for r in refs if r.startswith("chr"))
    return "chr" if chr_like >= (len(refs) / 2) else "nochr"


def normalize_contig(contig: str, target_style: str) -> str:
    """Convert a contig name to match BAM style ('chr' vs 'nochr')."""
    if target_style == "chr":
        return contig if contig.startswith("chr") else f"chr{contig}"
    return contig[3:] if contig.startswith("chr") else contig


def parse_rgb_to_mpl(rgb: str):
    """Convert 'R,G,B' (0-255) to matplotlib tuple (r,g,b) in 0..1. Return None if invalid."""
    try:
        parts = [p.strip() for p in rgb.split(",")]
        if len(parts) != 3:
            return None
        r, g, b = (int(parts[0]), int(parts[1]), int(parts[2]))
        r = max(0, min(255, r))
        g = max(0, min(255, g))
        b = max(0, min(255, b))
        return (r / 255.0, g / 255.0, b / 255.0)
    except Exception:
        return None


# ----------------------------
# BED parsing / harmonization
# ----------------------------
def parse_bed_file_with_colors(
    bed_file: str,
    state_col_1based: int = 4,
    rgb_col_1based: int = 9,
):
    """
    Parse BED-like file into:
      chromatin_states: dict[chrom] -> list[(start, end, state)]
      state_colors: dict[state] -> matplotlib RGB tuple (0..1), if provided

    Columns are 1-based (like typical CLI tools). RGB is optional.
    """
    state_idx = state_col_1based - 1
    rgb_idx = rgb_col_1based - 1

    chromatin_states = defaultdict(list)
    state_colors = {}

    with open(bed_file, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            fields = s.split()
            if len(fields) < 3:
                raise ValueError(f"BED line has <3 columns: {line.rstrip()}")
            if len(fields) <= state_idx:
                raise ValueError(
                    f"BED line does not have state column {state_col_1based}: {line.rstrip()}"
                )

            chrom = fields[0]
            start, end = int(fields[1]), int(fields[2])
            state = fields[state_idx]
            chromatin_states[chrom].append((start, end, state))

            if len(fields) > rgb_idx:
                mpl = parse_rgb_to_mpl(fields[rgb_idx])
                if mpl is not None and state not in state_colors:
                    state_colors[state] = mpl

    # Sort intervals per chromosome for fast linear scan assignment
    for chrom in chromatin_states:
        chromatin_states[chrom].sort(key=lambda x: (x[0], x[1]))

    return chromatin_states, state_colors


def harmonize_bed_to_bam(chromatin_states: dict, bam_style: str) -> dict:
    """Normalize BED contig keys to match BAM contig naming style."""
    out = defaultdict(list)
    for chrom, intervals in chromatin_states.items():
        norm = normalize_contig(chrom, bam_style)
        out[norm].extend(intervals)
    for chrom in out:
        out[chrom].sort(key=lambda x: (x[0], x[1]))
    return out


def autosomes_1_22_present(bamfile: pysam.AlignmentFile, bam_style: str) -> list[str]:
    """Return autosomes 1-22 in BAM style, filtered to those present in the BAM."""
    chroms = [normalize_contig(str(i), bam_style) for i in range(1, 23)]
    present = set(bamfile.references)
    return [c for c in chroms if c in present]


# ----------------------------
# BAM read pairing and assignment
# ----------------------------
def generate_paired_reads_from_list(reads):
    """
    Pair reads by query_name within the given list.

    This is chunk-local pairing: mates that never appear in the same chunk won't be paired.
    """
    unpaired = {}
    for read in reads:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        name = read.query_name
        other = unpaired.pop(name, None)
        if other is None:
            unpaired[name] = read
            continue

        # Skip if orientations are not a proper FR pair (coarse filter)
        if read.is_reverse == other.is_reverse:
            continue

        # Yield in (read1_forward, read2_reverse) order
        if not read.is_reverse:
            yield read, other
        else:
            yield other, read


def process_chunk(reads_in_chunk, current_states, length_counts_by_state, state_idx):
    """
    For a list of reads from a chromosome chunk:
      - pair reads by name
      - compute fragment_start/fragment_end and fragment length
      - assign the fragment to the BED interval containing fragment_midpoint
    Returns updated state_idx for efficient scanning through sorted intervals.
    """
    paired = list(generate_paired_reads_from_list(reads_in_chunk))
    paired.sort(key=lambda x: min(x[0].reference_start, x[1].reference_start))

    for read1, read2 in paired:
        fragment_start = min(
            read1.reference_start, read1.reference_end,
            read2.reference_start, read2.reference_end,
        )
        fragment_end = max(
            read1.reference_start, read1.reference_end,
            read2.reference_start, read2.reference_end,
        )
        fragment_length = fragment_end - fragment_start
        if fragment_length <= 0:
            continue

        # Use midpoint for assignment (integer bp coordinate)
        fragment_mid = (fragment_start + fragment_end) // 2

        # Advance interval pointer until we find the interval whose end is > fragment_mid
        while state_idx < len(current_states) and fragment_mid >= current_states[state_idx][1]:
            state_idx += 1

        if state_idx < len(current_states):
            istart, iend, state = current_states[state_idx]
            if istart <= fragment_mid < iend:
                length_counts_by_state[state][fragment_length] += 1

    return state_idx


def process_chromosome_in_chunks(
    bamfile,
    chromatin_states,
    chromosome,
    length_counts_by_state,
    chunk_size=1_000_000,
):
    """Process one chromosome in coordinate chunks and update length_counts_by_state."""
    current_states = chromatin_states.get(chromosome, [])
    if not current_states:
        print(f"WARNING: No BED intervals found for {chromosome}. Skipping.")
        return 0

    chrom_length = bamfile.get_reference_length(chromosome)
    state_idx = 0
    assigned = 0

    for start in tqdm(range(0, chrom_length, chunk_size), desc=f"Processing {chromosome}", unit="chunk"):
        end = min(start + chunk_size, chrom_length)

        reads_in_chunk = []
        for read in bamfile.fetch(chromosome, start, end):
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            reads_in_chunk.append(read)

        before = sum(sum(d.values()) for d in length_counts_by_state.values())
        state_idx = process_chunk(reads_in_chunk, current_states, length_counts_by_state, state_idx)
        after = sum(sum(d.values()) for d in length_counts_by_state.values())
        assigned += (after - before)

    return assigned


# ----------------------------
# Synthetic categories + color overrides
# ----------------------------
def add_euchromatin_category(
    length_counts_by_state,
    state_colors,
    euchromatin_prefixes,
    euchromatin_name,
    euchromatin_color,
):
    """Create a pooled 'euchromatin' category by summing states with matching prefixes."""
    eu_counts = defaultdict(int)
    for state, counts in list(length_counts_by_state.items()):
        if any(state.startswith(p) for p in euchromatin_prefixes):
            for L, c in counts.items():
                eu_counts[L] += c

    if eu_counts:
        length_counts_by_state[euchromatin_name] = eu_counts
        if state_colors is not None:
            state_colors[euchromatin_name] = euchromatin_color


def add_whole_genome_category(
    length_counts_by_state,
    state_colors,
    whole_name,
    whole_color,
    exclude_states,
):
    """Create 'whole-genome' as the sum of all (non-excluded) state counts."""
    wg_counts = defaultdict(int)
    for state, counts in list(length_counts_by_state.items()):
        if state in exclude_states:
            continue
        for L, c in counts.items():
            wg_counts[L] += c

    if wg_counts:
        length_counts_by_state[whole_name] = wg_counts
        if state_colors is not None:
            state_colors[whole_name] = whole_color


def apply_color_overrides(
    state_colors,
    length_counts_by_state,
    euchromatin_name,
    euchromatin_color,
    whole_name,
    whole_color,
    grey_prefixes,
    grey_color,
):
    """
    Force final plotting colors:
      - euchromatin: fixed color
      - whole-genome: fixed color
      - states with any prefix in grey_prefixes: fixed mid-grey
    """
    state_colors[euchromatin_name] = euchromatin_color
    state_colors[whole_name] = whole_color
    for state in length_counts_by_state.keys():
        if any(state.startswith(p) for p in grey_prefixes):
            state_colors[state] = grey_color


# ----------------------------
# Output + plotting
# ----------------------------
def write_output(length_counts_by_state, bamfile_path, output_dir, label):
    """Write one output per state: <bam>_<state>_<label>.txt with length<TAB>count."""
    bam_base = os.path.basename(bamfile_path)
    bam_base_noext = bam_base[:-4] if bam_base.endswith(".bam") else os.path.splitext(bam_base)[0]

    for state, counts in length_counts_by_state.items():
        if not counts:
            continue
        sanitized_state = sanitize_filename(state)
        out_path = os.path.join(output_dir, f"{bam_base_noext}_{sanitized_state}_{label}.txt")
        with open(out_path, "w") as f:
            for length in sorted(counts):
                f.write(f"{length}\t{counts[length]}\n")


def plot_density_overlay(
    length_counts_by_state,
    state_colors,
    out_png,
    x_min,
    x_max,
    title=None,
    max_states=None,
    include_states=None,
):
    """
    Overlaid density curves for each state.

    Density is normalized *within [x_min, x_max]* for each state:
      dens(L) = count(L) / sum_{x_min..x_max} count(L)

    Ordering: states with the highest peak density are plotted first.
    """
    import matplotlib.pyplot as plt

    xs = list(range(x_min, x_max + 1))

    items = []
    for state, counts in length_counts_by_state.items():
        if include_states is not None and state not in include_states:
            continue

        ys = [counts.get(L, 0) for L in xs]
        total = sum(ys)
        if total == 0:
            continue

        dens = [y / total for y in ys]
        items.append((max(dens), total, state, dens))

    items.sort(key=lambda t: (t[0], t[1]), reverse=True)
    if max_states is not None:
        items = items[:max_states]

    plt.figure(figsize=(11, 6))
    for max_dens, total, state, dens in items:
        color = (state_colors or {}).get(state)
        label = f"{state} (n={total:,}, max={max_dens:.4g})"
        if color is None:
            plt.plot(xs, dens, label=label)
        else:
            plt.plot(xs, dens, label=label, color=color)

    plt.xlim(x_min, x_max)
    plt.xlabel("Fragment length (bp)")
    plt.ylabel(f"Density (normalized within {x_min}–{x_max} bp)")
    if title:
        plt.title(title)

    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def build_focus_subset(length_counts_by_state, focus_prefixes, include_names):
    """
    Return a set of states to include in the focus plot.

    - Any state whose name starts with a prefix in focus_prefixes
    - Any state exactly matching a name in include_names
    """
    keep = set()
    for state in length_counts_by_state.keys():
        if state in include_names or any(state.startswith(p) for p in focus_prefixes):
            keep.add(state)
    return keep


# ----------------------------
# CLI / main
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Assign paired-end BAM fragments to chromatin states (BED) and compute fragment length distributions."
    )
    p.add_argument("-b", "--bamfile", required=True, help="Path to BAM file (coordinate-sorted recommended).")
    p.add_argument("-d", "--bedfile", required=True, help="Path to BED file containing chromatin-state intervals.")

    # BED format controls
    p.add_argument("--bed_state_col", type=int, default=4, help="1-based BED column for state label (default: 4).")
    p.add_argument("--bed_rgb_col", type=int, default=9, help="1-based BED column for RGB 'R,G,B' (default: 9).")

    # Output controls
    p.add_argument(
        "-o",
        "--output_dir",
        default=None,
        help="Output directory (default: ./<bam_basename_without_ext>/).",
    )

    # What to process
    p.add_argument(
        "--chromosome",
        default=None,
        help="Chromosome to process (e.g. '1' or 'chr1'). If omitted, pools autosomes 1–22 (default).",
    )
    p.add_argument("--chunk_size", type=int, default=1_000_000, help="Chunk size in bp (default: 1,000,000).")

    # Synthetic categories
    p.add_argument(
        "--euchromatin_prefixes",
        default="1_,2_,4_,5_,6_,7_,9_,10_,11_",
        help="Comma-separated prefixes to pool into euchromatin (default matches your prior hardcode).",
    )
    p.add_argument("--euchromatin_name", default="euchromatin", help="Name for pooled euchromatin category.")
    p.add_argument(
        "--whole_genome_name",
        default="whole-genome",
        help="Name for pooled whole-genome category.",
    )

    # Color overrides
    p.add_argument(
        "--euchromatin_color_rgb",
        default="128,0,128",
        help="Euchromatin color as 'R,G,B' (0–255). Default: 128,0,128 (purple).",
    )
    p.add_argument(
        "--whole_genome_color_rgb",
        default="0,0,0",
        help="Whole-genome color as 'R,G,B' (0–255). Default: 0,0,0 (black).",
    )
    p.add_argument(
        "--grey_override_prefixes",
        default="13_",
        help="Comma-separated prefixes to force mid-grey color (default: '13_').",
    )
    p.add_argument(
        "--grey_override_color_rgb",
        default="128,128,128",
        help="Grey override color as 'R,G,B' (0–255). Default: 128,128,128.",
    )

    # Plotting
    p.add_argument(
        "--plot",
        action="store_true",
        help="Make density plots over a length window. Writes main plot + focus plot.",
    )
    p.add_argument("--plot_xmin", type=int, default=0, help="Plot x-axis min length (default: 0).")
    p.add_argument("--plot_xmax", type=int, default=500, help="Plot x-axis max length (default: 500).")
    p.add_argument(
        "--plot_out",
        default=None,
        help="Main plot PNG path (default: <output_dir>/<bam>_<label>_dens_<xmin>-<xmax>.png).",
    )
    p.add_argument(
        "--plot_max_states",
        type=int,
        default=None,
        help="If set, only plot the top N states by peak density (main plot only).",
    )
    p.add_argument(
        "--focus_prefixes",
        default="8_,13_",
        help="Comma-separated prefixes for the focus plot (default: '8_,13_').",
    )

    args = p.parse_args()

    # Determine default output directory if not provided
    bam_base = os.path.basename(args.bamfile)
    bam_base_noext = bam_base[:-4] if bam_base.endswith(".bam") else os.path.splitext(bam_base)[0]
    if args.output_dir is None:
        args.output_dir = os.path.join(os.getcwd(), bam_base_noext)
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse colors from args
    euch_color = parse_rgb_to_mpl(args.euchromatin_color_rgb) or DEFAULT_EUCHROMATIN_COLOR
    wg_color = parse_rgb_to_mpl(args.whole_genome_color_rgb) or DEFAULT_WHOLE_GENOME_COLOR
    grey_color = parse_rgb_to_mpl(args.grey_override_color_rgb) or DEFAULT_MID_GREY

    euch_prefixes = tuple([p.strip() for p in args.euchromatin_prefixes.split(",") if p.strip()])
    grey_prefixes = tuple([p.strip() for p in args.grey_override_prefixes.split(",") if p.strip()])
    focus_prefixes = tuple([p.strip() for p in args.focus_prefixes.split(",") if p.strip()])

    # Load BED intervals and optional colors
    chromatin_states_raw, state_colors = parse_bed_file_with_colors(
        args.bedfile,
        state_col_1based=args.bed_state_col,
        rgb_col_1based=args.bed_rgb_col,
    )

    # Open BAM, harmonize contig naming between BED and BAM
    bamfile = pysam.AlignmentFile(args.bamfile, "rb")
    bam_style = detect_bam_chr_style(bamfile)
    chromatin_states = harmonize_bed_to_bam(chromatin_states_raw, bam_style)

    # Combined counts across all processed chromosomes
    length_counts_by_state = defaultdict(lambda: defaultdict(int))

    # Decide which chromosome(s) to process
    if args.chromosome is None:
        chroms_to_process = autosomes_1_22_present(bamfile, bam_style)
        if not chroms_to_process:
            raise SystemExit("ERROR: None of autosomes 1–22 were found in BAM references.")
        label = "autosomes1-22"
        print("No --chromosome provided; pooling autosomes 1–22:")
        print("  " + ", ".join(chroms_to_process))
    else:
        chrom = normalize_contig(args.chromosome, bam_style)
        if chrom not in bamfile.references:
            example = bamfile.references[0] if bamfile.references else "(none)"
            raise SystemExit(
                f"ERROR: Chromosome '{chrom}' not found in BAM references.\n"
                f"  You requested: {args.chromosome}\n"
                f"  BAM style: {bam_style}\n"
                f"  Example BAM contig: {example}\n"
            )
        chroms_to_process = [chrom]
        label = chrom

    # Process BAM by chromosome, in chunks
    for chrom in chroms_to_process:
        print(f"\nProcessing chromosome: {chrom} (BAM style: {bam_style})")
        assigned = process_chromosome_in_chunks(
            bamfile=bamfile,
            chromatin_states=chromatin_states,
            chromosome=chrom,
            length_counts_by_state=length_counts_by_state,
            chunk_size=args.chunk_size,
        )
        print(f"  Assigned fragments (this chrom): {assigned:,}")

    bamfile.close()

    # Add synthetic categories (pooled)
    add_euchromatin_category(
        length_counts_by_state,
        state_colors=state_colors,
        euchromatin_prefixes=euch_prefixes,
        euchromatin_name=args.euchromatin_name,
        euchromatin_color=euch_color,
    )
    add_whole_genome_category(
        length_counts_by_state,
        state_colors=state_colors,
        whole_name=args.whole_genome_name,
        whole_color=wg_color,
        exclude_states=(args.euchromatin_name, args.whole_genome_name),
    )

    # Enforce final color overrides (including grey prefixes)
    apply_color_overrides(
        state_colors=state_colors,
        length_counts_by_state=length_counts_by_state,
        euchromatin_name=args.euchromatin_name,
        euchromatin_color=euch_color,
        whole_name=args.whole_genome_name,
        whole_color=wg_color,
        grey_prefixes=grey_prefixes,
        grey_color=grey_color,
    )

    # Write outputs
    write_output(length_counts_by_state, args.bamfile, args.output_dir, label)

    # Report (exclude synthetic to avoid confusing "total assigned")
    assigned_total = sum(sum(d.values()) for d in length_counts_by_state.values())
    nstates = sum(1 for d in length_counts_by_state.values() if sum(d.values()) > 0)
    print(f"\nDone. Wrote {nstates} state files (including any synthetic categories).")
    print(f"Output directory: {args.output_dir}")

    # Optional plotting
    if args.plot:
        xmin, xmax = args.plot_xmin, args.plot_xmax
        out_png = args.plot_out or os.path.join(
            args.output_dir, f"{bam_base_noext}_{label}_dens_{xmin}-{xmax}.png"
        )
        plot_density_overlay(
            length_counts_by_state=length_counts_by_state,
            state_colors=state_colors,
            out_png=out_png,
            x_min=xmin,
            x_max=xmax,
            title=f"{bam_base_noext} {label}: fragment length density by chromatin state",
            max_states=args.plot_max_states,
            include_states=None,
        )
        print(f"Wrote plot: {out_png}")

        # Focus plot (e.g., euchromatin + selected prefixes)
        focus_states = build_focus_subset(
            length_counts_by_state,
            focus_prefixes=focus_prefixes,
            include_names={args.euchromatin_name},
        )
        out_png2 = os.path.join(
            args.output_dir, f"{bam_base_noext}_{label}_dens_{xmin}-{xmax}_focus.png"
        )
        plot_density_overlay(
            length_counts_by_state=length_counts_by_state,
            state_colors=state_colors,
            out_png=out_png2,
            x_min=xmin,
            x_max=xmax,
            title=f"{bam_base_noext} {label}: focus states ({args.euchromatin_name} + {','.join(focus_prefixes)})",
            max_states=None,
            include_states=focus_states,
        )
        print(f"Wrote plot: {out_png2}")


if __name__ == "__main__":
    main()
