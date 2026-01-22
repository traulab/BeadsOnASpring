#!/usr/bin/env bash
set -euo pipefail

# extract_bigwig_windows.sh
#
# For each record in an "alignment" BED-like file, extract a fixed window from a BigWig
# around an alignment point, expand to per-base values (nan where missing),
# optionally reverse-complement the vector based on strand, and write one line per record.
#
# Default alignment point: midpoint between start & end (BED columns 2&3, 0-based, half-open).
# Optional alignment point: take value from a specified column (1-based column index).
#
# Requirements: bigWigToBedGraph, awk

usage() {
  cat <<'EOF'
Usage:
  extract_bigwig_windows.sh -w BIGWIG -a ALIGNMENT_FILE [options]

Required:
  -w, --bigwig FILE           Input BigWig
  -a, --alignment-file FILE   Input BED/BED-like alignment file (tab-delimited)

Options:
  -o, --out FILE              Output file (default: <alignment_basename>_BW.txt)
  --window-half N             Half-window size in bp (default: 3001)
                              Output length = 2*N (default: 6002 bases)
  --chrom-col N               1-based column for chrom (default: 1)
  --start-col N               1-based column for start (default: 2)
  --end-col N                 1-based column for end (default: 3)
  --strand-col N              1-based column for strand (default: 6). If missing/., random +/-.
  --point-col N               1-based column for alignment point (absolute genome coordinate).
                              If set, this overrides midpoint calculation.
                              (Value should be an integer base coordinate in the same coordinate
                               system you want to pass to bigWigToBedGraph -start/-end.)
  --skip-header               Skip the first line of the alignment file
  -h, --help                  Show this help

Example (midpoint):
  extract_bigwig_windows.sh -w track.bw -a peaks.bed

Example (use column 7 as point):
  extract_bigwig_windows.sh -w track.bw -a peaks.bed --point-col 7

Notes:
- bigWigToBedGraph expects 0-based half-open coordinates for -start/-end.
- If strand is "-", the output vector is reversed.
EOF
}

# Defaults
BIGWIG=""
ALIGNMENT_FILE=""
OUT=""
WINDOW_HALF=3000
CHROM_COL=1
START_COL=2
END_COL=3
STRAND_COL=6
POINT_COL=0
SKIP_HEADER=0

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|--bigwig) BIGWIG="$2"; shift 2 ;;
    -a|--alignment-file) ALIGNMENT_FILE="$2"; shift 2 ;;
    -o|--out) OUT="$2"; shift 2 ;;
    --window-half) WINDOW_HALF="$2"; shift 2 ;;
    --chrom-col) CHROM_COL="$2"; shift 2 ;;
    --start-col) START_COL="$2"; shift 2 ;;
    --end-col) END_COL="$2"; shift 2 ;;
    --strand-col) STRAND_COL="$2"; shift 2 ;;
    --point-col) POINT_COL="$2"; shift 2 ;;
    --skip-header) SKIP_HEADER=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

# Validate
if [[ -z "$BIGWIG" || -z "$ALIGNMENT_FILE" ]]; then
  echo "ERROR: --bigwig and --alignment-file are required." >&2
  usage
  exit 2
fi
if [[ ! -f "$BIGWIG" ]]; then
  echo "ERROR: BigWig not found: $BIGWIG" >&2
  exit 2
fi
if [[ ! -f "$ALIGNMENT_FILE" ]]; then
  echo "ERROR: Alignment file not found: $ALIGNMENT_FILE" >&2
  exit 2
fi
if ! command -v bigWigToBedGraph >/dev/null 2>&1; then
  echo "ERROR: bigWigToBedGraph not found in PATH." >&2
  exit 2
fi
if ! command -v awk >/dev/null 2>&1; then
  echo "ERROR: awk not found in PATH." >&2
  exit 2
fi

# Default output name
if [[ -z "$OUT" ]]; then
  PREFIX="$(basename "$ALIGNMENT_FILE")"
  PREFIX="${PREFIX%.bed}"
  PREFIX="${PREFIX%.BED}"
  OUT="${PREFIX}_BW.txt"
fi

WINDOW_LEN=$(( 2 * WINDOW_HALF + 1))

# Process
# We let awk parse columns robustly and choose alignment point.
# Strand: if not +/-, randomize using awk's srand/rand so it's per-line stochastic.
awk \
  -v CHROM_COL="$CHROM_COL" \
  -v START_COL="$START_COL" \
  -v END_COL="$END_COL" \
  -v STRAND_COL="$STRAND_COL" \
  -v POINT_COL="$POINT_COL" \
  -v WINDOW_HALF="$WINDOW_HALF" \
  -v WINDOW_LEN="$WINDOW_LEN" \
  -v BIGWIG="$BIGWIG" \
  -v SKIP_HEADER="$SKIP_HEADER" \
'BEGIN{
  OFS="\t";
  srand();
}
function get_field(arr, n,   v){
  v = (n <= length(arr) ? arr[n] : "");
  return v;
}
NR==1 && SKIP_HEADER==1 { next }
# Skip blank/comment lines
/^[[:space:]]*$/ { next }
$1 ~ /^#/ { next }
{
  # Split on tabs/spaces
  n = split($0, f, /[[:space:]]+/)

  chr = get_field(f, CHROM_COL)
  s0  = get_field(f, START_COL)
  e0  = get_field(f, END_COL)
  strand = get_field(f, STRAND_COL)

  if (chr == "" || s0 == "" || e0 == "") next
  start = int(s0)
  end   = int(e0)

  # Choose alignment point
  if (POINT_COL > 0) {
    p0 = get_field(f, POINT_COL)
    if (p0 == "" || p0 ~ /[^0-9-]/) next
    center = int(p0)
  } else {
    center = int((start + end) / 2)
  }

  start_window = center - WINDOW_HALF
  end_window   = center + WINDOW_HALF + 1

  # Normalize/randomize strand
  if (strand != "+" && strand != "-") {
    if (rand() < 0.5) strand = "-"
    else strand = "+"
  }

  # Run bigWigToBedGraph and expand to per-base vector of length WINDOW_LEN.
  cmd = "bigWigToBedGraph \"" BIGWIG "\" stdout -chrom=\"" chr "\" -start=" start_window " -end=" end_window

  # init
  delete scores
  last = start_window
  i = 0

  while ((cmd | getline line) > 0) {
    m = split(line, b, "\t")
    if (m < 4) continue
    bstart = int(b[2]); bend = int(b[3]); val = b[4]

    # fill gaps
    while (last < bstart && i < WINDOW_LEN) {
      scores[i++] = "nan"
      last++
    }
    # fill covered
    for (j = bstart; j < bend && i < WINDOW_LEN; j++) {
      scores[i++] = val
      last++
    }
    if (i >= WINDOW_LEN) break
  }
  close(cmd)

  # fill tail
  while (i < WINDOW_LEN) scores[i++] = "nan"

  # output vector (reverse if strand == "-")
  if (strand == "-") {
    for (j = WINDOW_LEN - 1; j >= 0; j--) {
      printf "%s%s", scores[j], (j==0 ? "\n" : " ")
    }
  } else {
    for (j = 0; j < WINDOW_LEN; j++) {
      printf "%s%s", scores[j], (j==WINDOW_LEN-1 ? "\n" : " ")
    }
  }
}
' "$ALIGNMENT_FILE" > "$OUT"

echo "Wrote: $OUT" >&2
