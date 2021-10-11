#!/usr/bin/env python3

"""
FragScorer scores fragmentomics data.

Scores are calculated in the 'score_range' function, to add a new score
simply add a new array to the dictionary returned by that function and
add code to the for-loop to adjust the score for each fragment.
"""

import sys
from tqdm import tqdm
import argparse
import pysam
import numpy as np
import pyBigWig

def generate_paired_reads(bamfile, contig=None, start=None, end=None):
    """Generator taking paired-end data from a sam/bam file and yielding both pairs as (fwd, rev)"""
    _unpaired_reads = dict()
    for read in bamfile.fetch(contig, start, end, multiple_iterators=True):
        # Skip reads that aren't part of a fully-mapped pair
        if not read.is_proper_pair:
            continue
        
        name = read.query_name
        if name not in _unpaired_reads:
            _unpaired_reads[name] = read
            continue
        
        # Have found mate, yield both
        mate = _unpaired_reads[name]
        del _unpaired_reads[name]
        assert read.is_reverse != mate.is_reverse, "Both pairs mapped to same strand?"
        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read

def generate_fragment_ranges(bamfile, contig, start, end):
    """Generate ranges covered by sequenced DNA fragments"""
    for r_fwd, r_rev in generate_paired_reads(bamfile, contig, start, end):
        assert r_fwd.reference_start <= r_rev.reference_start
        assert r_rev.reference_end >= r_fwd.reference_end
        yield r_fwd.reference_start, r_rev.reference_end

def score_contig(bamfile, contig):
    """
    Applies a range of scoring algorithms to a given reference contig, returning the results.
    """

    # Get total number of fragments from the index, if possible
    num_mapped_reads = None
    for idx_stat in bamfile.get_index_statistics():
        if idx_stat.contig == contig:
            num_mapped_reads = idx_stat.mapped
            break
    assert num_mapped_reads is not None, \
            "contig not present in index stats"
    num_pairs = num_mapped_reads // 2

    ## Initialize score arrays
    ref_len = bamfile.get_reference_length(contig)
    coverage = np.zeros(ref_len, dtype=int)   # Coverage
    l_wps    = np.zeros(ref_len, dtype=int)   # Windowed Protection Score (long fraction)
    s_wps    = np.zeros(ref_len, dtype=int)   # Windowed Protection Score (short fraction)
    l_bds    = np.zeros(ref_len, dtype=float) # Breakpoint Density Score (Blackman window, window = 120)
    sws      = np.zeros(ref_len, dtype=float) # Scaled Window Score (Blackman window)

    ## Store array references in dict (this gets returned later)
    scores = {
            'coverage': coverage,
            'l_wps': l_wps,
            's_wps': s_wps,
            'l_bds': l_bds,
            'sws': sws
    }

    # internal variables
    l_wps_frag_range = range(120, 180)
    l_wps_window = 120
    s_wps_frag_range = range(35, 80)
    s_wps_window = 30
    l_bds_window = np.blackman(120)
    sws_frag_range = range(120, 180)
    sws_windows = [np.blackman(l) for l in sws_frag_range]

    ## Iterate over all reads in contig and calculate scores
    for start, end in tqdm(
            generate_fragment_ranges(bamfile, contig, 0, ref_len),
            total=num_pairs, unit='frags'):

        frag_length = end-start

        # Calculate coverage
        coverage[start:end] += 1

        # Calculate l_wps
        if frag_length in l_wps_frag_range:           
            # Add bonus for coverage
            
            l_wps[start + l_wps_window//2:end - l_wps_window//2 + 1] += 1
            
            # Add penalty for proximity to fragment end
            l_wps[start - 1 - l_wps_window//2:start + l_wps_window//2] -= 1
            l_wps[end + 1 - l_wps_window//2:end + l_wps_window//2 + 2] -= 1
        
        # Calculate s_wps
        if frag_length in s_wps_frag_range:
        #     # Add bonus for coverage
            s_wps[start:end] += 1

        #     # Add penalty for proximity to fragment end
            s_wps[start - s_wps_window//2:start + s_wps_window//2] -= 1
            s_wps[end - s_wps_window//2:end - + s_wps_window//2] -= 1

        # Calculate l_bds
        l_bds[start - len(l_bds_window)//2:start + len(l_bds_window)//2] += l_bds_window
        l_bds[end - len(l_bds_window)//2:end + len(l_bds_window)//2] += l_bds_window

        # Calculate sws
        if frag_length in sws_frag_range:
            w = sws_windows[frag_length - sws_frag_range.start]
            sws[start:end] += w

    return scores

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Score fragmentomics data.",
            epilog="""By default, this script will process all contigs in
                   the supplied BAM file. To limit processing to specific
                   contigs, use the --contig argument.

                   Note that the supplied file must be indexed (i.e.
                   must have an associated .bai file).
                   """
            )
    parser.add_argument('bamfile', help='bam file to process')
    parser.add_argument('out_prefix', help='prefix to apply to output files')
    parser.add_argument('-c', '--contig',
            help='limit calculation to specified contig(s)',
            dest='contigs', action='append', nargs='*')
    args = parser.parse_args()

    # Open SAM file
    try:
        bamfile = pysam.AlignmentFile(args.bamfile, "rb")
    except FileNotFoundError as e:
        parser.error("Unable to open bamfile (file not found)")
        # Note: parser.error should exit, so this return never gets
        # reached in normal circumstances.
        return -1
    except e:
        parser.error("Unable to open bamfile: {}".format(str(e)))
        return -1

    # Generate list of ranges to calculate over
    contigs = []
    if args.contigs is not None:
        contigs = [x[0] for x in args.contigs]
    else:
        # No ranges specified, default to all contigs
        for contig in bamfile.references:
            contigs.append(contig)

    # Score the contigs
    results = {}
    for contig in tqdm(contigs, desc='Scoring contigs'):
        results[contig] = score_contig(bamfile, contig)

    # Rearrange data to allow bigwig files to be generated
    scores = {}
    for contig in results:
        contig_scores = results[contig]
        for score_type in contig_scores:
            if score_type not in scores:
                scores[score_type] = []
            scores[score_type].append((contig, contig_scores[score_type]))

    for score_type in tqdm(scores, desc="Writing files"):
        # Open file
        filename = args.out_prefix + '.{}.bw'.format(score_type)
        bw = pyBigWig.open(filename, 'w')
        # Create header
        header = []
        for contig, contig_length in zip(bamfile.references, bamfile.lengths):
            if contig in results:
                header.append((contig, contig_length))
        bw.addHeader(header)
        
        # Write data to file
        for contig, s in scores[score_type]:
            bw.addEntries(contig, 1, values=[float(x) for x in s], span=1, step=1)
        bw.close()

    return 0

if __name__ == '__main__':
    sys.exit(main())
