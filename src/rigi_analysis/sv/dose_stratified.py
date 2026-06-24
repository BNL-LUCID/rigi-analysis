#!/usr/bin/env python
"""Dose-stratified INV-DBS analysis - VERSION 4 (with Temporal Concordance).

Key changes from v3:
- REQUIRES temporal concordance: INV timepoint must match DBS pattern
- Only pairs where INV and DBS appear at the same timepoint are counted
- This ensures we're measuring true mutagenic coupling, not random co-location

Deduplication strategy:
1. DBS deduplication: Within dose bin only (not across all doses)
2. SV deduplication: Within dose×time only (preserves dose-response signal)
3. Dose binning: By dose RATE (A-C = Low, D-E = High)

Usage:
    python dose_stratified_inv_dbs_v4.py \
        --annotsv-dir path/to/annotsv \
        --mutation-dir path/to/mutations \
        --output-dir dose_analysis_v4 \
        --inv-size 0 \
        --window 10 \
        --mega-threshold 50000000
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

# Dose rates in mGy/hr
DOSE_RATES = {
    'A': 0.001,
    'B': 0.01,
    'C': 0.1,
    'D': 1.0,
    'E': 2.0
}

# Timepoints in hours
TIMEPOINTS = {
    'W1': 168,
    'W2': 336,
    'W3': 504
}

# Dose bins by dose RATE (not cumulative dose)
DOSE_BINS = {
    'Low': ['A', 'B', 'C'],
    'High': ['D', 'E']
}

# Radiation-induced temporal patterns
RADIATION_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']

# Map timepoint to pattern position
TIMEPOINT_TO_PATTERN_POS = {'W1': 1, 'W2': 2, 'W3': 3}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def calculate_cumulative_dose(dose_letter, timepoint):
    """Calculate cumulative dose in mGy."""
    dose_rate = DOSE_RATES.get(dose_letter, 0)
    hours = TIMEPOINTS.get(timepoint, 0)
    return dose_rate * hours


def get_dose_bin(dose_letter):
    """Get dose bin category based on dose RATE."""
    for bin_name, doses in DOSE_BINS.items():
        if dose_letter in doses:
            return bin_name
    return None


def find_column(df, candidates):
    """Find the first matching column name from candidates."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def normalize_chromosome(chrom_series):
    """Normalize chromosome names to 'chrN' format."""
    chrom_series = chrom_series.astype(str)
    return chrom_series.apply(
        lambda x: x if str(x).startswith('chr') else f'chr{x}'
    )


def timepoint_matches_pattern(timepoint, pattern):
    """Check if a timepoint (W1/W2/W3) is present in a DBS pattern.

    Pattern format: '0T00' where positions 1,2,3 = W1,W2,W3
    'T' or 'C' = present, '0' = absent
    """
    if len(pattern) < 4:
        return False
    pos = TIMEPOINT_TO_PATTERN_POS.get(timepoint, 0)
    if pos > 0:
        return pattern[pos] != '0'
    return False


# =============================================================================
# SV LOADING
# =============================================================================

def load_annotsv_files(annotsv_dir, chromosomes=None):
    """Load AnnotSV files for all doses and timepoints.

    DEDUPLICATION: Within each dose×time combination only.
    Same SV at different doses = KEPT (biological signal).
    """
    print("=" * 70)
    print("LOADING STRUCTURAL VARIANTS FROM ANNOTSV FILES")
    print("=" * 70)

    if chromosomes:
        print(f"Filtering to chromosomes: {', '.join(chromosomes)}")
    else:
        print("Analyzing genome-wide (all chromosomes)")

    annotsv_path = Path(annotsv_dir)

    # Find AnnotSV files
    patterns = ['d0_vs_d[A-E]_W[1-3]_annotated.tsv', 'd0_vs_d[A-E]_W[1-3]*.tsv']
    all_files = set()
    for pattern in patterns:
        all_files.update(annotsv_path.glob(pattern))

    # Also try without glob patterns
    for d in ['A', 'B', 'C', 'D', 'E']:
        for w in ['W1', 'W2', 'W3']:
            all_files.update(annotsv_path.glob(f'd0_vs_d{d}_{w}*.tsv'))

    all_files = sorted(list(all_files))

    if not all_files:
        print(f"ERROR: No AnnotSV files found in {annotsv_dir}")
        return pd.DataFrame()

    print(f"Found {len(all_files)} AnnotSV files\n")

    all_svs = []

    for file_path in all_files:
        filename = file_path.name

        # Parse dose and timepoint
        dose = None
        timepoint = None

        for d in ['A', 'B', 'C', 'D', 'E']:
            if f'_vs_d{d}_' in filename:
                dose = d
                break

        for w in ['W1', 'W2', 'W3']:
            if f'_{w}_' in filename or f'_{w}.' in filename:
                timepoint = w
                break

        if dose is None or timepoint is None:
            print(f"  WARNING: Could not parse dose/timepoint from: {filename}")
            continue

        cumulative_dose = calculate_cumulative_dose(dose, timepoint)
        dose_bin = get_dose_bin(dose)

        print(f"  {filename}")
        print(f"    Dose: {dose} ({DOSE_RATES[dose]} mGy/hr), Time: {timepoint}, Bin: {dose_bin}")

        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            raw_count = len(df)

            # Filter to full annotations only
            annotsv_type_col = find_column(df, ['AnnotSV_type', 'Annotation_mode'])
            if annotsv_type_col:
                df = df[df[annotsv_type_col] == 'full'].copy()
                print(f"    Rows: {raw_count:,} → {len(df):,} (after 'full' filter)")
            else:
                print(f"    Rows: {raw_count:,} (WARNING: no AnnotSV_type column)")

            if len(df) == 0:
                continue

            # Normalize chromosome
            chr_col = find_column(df, ['SV_chrom', 'Chromosome', 'Chrom', '#Chromosome'])
            if chr_col is None:
                print("    WARNING: No chromosome column found")
                continue

            df['Chromosome'] = normalize_chromosome(df[chr_col])

            # Filter to target chromosomes
            if chromosomes:
                df = df[df['Chromosome'].isin(chromosomes)].copy()
                if len(df) == 0:
                    continue

            # Add dose metadata
            df['Dose'] = dose
            df['Timepoint'] = timepoint
            df['Dose_Time'] = f"{dose}_{timepoint}"
            df['Cumulative_Dose_mGy'] = cumulative_dose
            df['Dose_Bin'] = dose_bin

            print(f"    Loaded: {len(df):,} SVs")
            all_svs.append(df)

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    if not all_svs:
        print("\nERROR: No SV data loaded!")
        return pd.DataFrame()

    combined = pd.concat(all_svs, ignore_index=True)
    print(f"\n{'─' * 70}")
    print(f"Combined SVs from all files: {len(combined):,}")

    # Find and standardize coordinate columns
    start_col = find_column(combined, ['SV_start', 'Start', 'start'])
    end_col = find_column(combined, ['SV_end', 'End', 'end'])
    type_col = find_column(combined, ['SV_type', 'Type', 'SV_Type'])

    if not all([start_col, end_col, type_col]):
        print("ERROR: Missing required columns (start/end/type)")
        return pd.DataFrame()

    combined = combined.rename(columns={start_col: 'SV_Start', end_col: 'SV_End', type_col: 'SV_Type'})

    # Deduplication within dose×time
    print("\nDeduplicating SVs WITHIN each dose×time combination...")

    def create_sv_key(row):
        chrom = str(row['Chromosome'])
        start = int(row['SV_Start'])
        end = int(row['SV_End'])
        sv_type = str(row['SV_Type'])
        dose_time = str(row['Dose_Time'])
        coord1, coord2 = sorted([start, end])
        return f"{chrom}:{coord1}-{coord2}:{sv_type}:{dose_time}"

    combined['SV_Key'] = combined.apply(create_sv_key, axis=1)

    before_dedup = len(combined)
    combined = combined.drop_duplicates(subset='SV_Key', keep='first')
    after_dedup = len(combined)

    print(f"  Before: {before_dedup:,} → After: {after_dedup:,}")
    print(f"  Removed: {before_dedup - after_dedup:,} technical duplicates")

    # Summary
    print("\nSV counts by dose bin:")
    for bin_name in ['Low', 'High']:
        count = (combined['Dose_Bin'] == bin_name).sum()
        print(f"  {bin_name}: {count:,}")

    return combined


# =============================================================================
# INVERSION FILTERING
# =============================================================================

def filter_inversions(sv_df, min_size=0, max_size=None, mega_threshold=50_000_000):
    """Filter SVs to inversions and add Is_Mega flag."""
    print(f"\n{'=' * 70}")
    print("FILTERING TO INVERSIONS")
    print("=" * 70)

    inv_df = sv_df[sv_df['SV_Type'] == 'INV'].copy()
    print(f"Total inversions: {len(inv_df):,}")

    # Calculate size
    inv_df['SV_Length'] = abs(inv_df['SV_End'] - inv_df['SV_Start'])

    # Add Is_Mega flag
    inv_df['Is_Mega'] = inv_df['SV_Length'] >= mega_threshold
    mega_count = inv_df['Is_Mega'].sum()
    small_count = (~inv_df['Is_Mega']).sum()
    print(f"  Mega (≥{mega_threshold/1e6:.0f}Mb): {mega_count:,}")
    print(f"  Small (<{mega_threshold/1e6:.0f}Mb): {small_count:,}")

    # Apply size filters if specified
    if min_size > 0:
        inv_df = inv_df[inv_df['SV_Length'] >= min_size].copy()
        print(f"After min size filter (≥{min_size:,}bp): {len(inv_df):,}")

    if max_size is not None:
        inv_df = inv_df[inv_df['SV_Length'] <= max_size].copy()
        print(f"After max size filter (≤{max_size:,}bp): {len(inv_df):,}")

    # Summary by dose bin
    print("\nInversions by dose bin:")
    for bin_name in ['Low', 'High']:
        bin_inv = inv_df[inv_df['Dose_Bin'] == bin_name]
        mega = bin_inv['Is_Mega'].sum()
        small = (~bin_inv['Is_Mega']).sum()
        print(f"  {bin_name}: {len(bin_inv):,} total ({mega:,} mega, {small:,} small)")

    return inv_df


# =============================================================================
# DBS LOADING
# =============================================================================

def load_dbs_mutations(mutation_dir, chromosomes=None):
    """Load DBS mutations with dose and Pattern information.

    DEDUPLICATION: Within each dose only (not across doses).
    Same DBS at different doses = KEPT as separate events.
    """
    print(f"\n{'=' * 70}")
    print("LOADING DBS MUTATIONS")
    print("=" * 70)

    if chromosomes:
        print(f"Filtering to chromosomes: {', '.join(chromosomes)}")
    else:
        print("Loading genome-wide")

    mutation_path = Path(mutation_dir)
    dbs_files = list(mutation_path.glob('DBS/DBS_dose_*_merged.csv'))

    if not dbs_files:
        dbs_files = list(mutation_path.glob('**/DBS_dose_*_merged.csv'))

    if not dbs_files:
        print(f"ERROR: No DBS files found in {mutation_dir}")
        return pd.DataFrame()

    print(f"Found {len(dbs_files)} DBS files\n")

    all_dbs = []

    # Extract dose label from filenames like:
    #   DBS_dose_A_merged.csv   -> "A"
    #   DBS_dose_dA_merged.csv  -> "A"  (current convention from merge_annotation.py)
    dose_re = re.compile(r'dose_d?([A-E])_merged', re.IGNORECASE)

    for file_path in sorted(dbs_files):
        filename = file_path.name

        m = dose_re.search(filename)
        if not m:
            print(f"  WARNING: Could not parse dose from: {filename}")
            continue
        dose = m.group(1).upper()

        dose_bin = get_dose_bin(dose)
        print(f"  {filename} → Dose {dose} ({dose_bin})")

        try:
            df = pd.read_csv(file_path, low_memory=False)
            raw_count = len(df)

            # Filter out control
            if 'Sample' in df.columns:
                df = df[~df['Sample'].str.contains('d0', case=False, na=False)].copy()

            # Filter to radiation patterns
            if 'Pattern' in df.columns:
                df = df[df['Pattern'].isin(RADIATION_PATTERNS)].copy()
                print(f"    Filtered to radiation patterns: {raw_count:,} → {len(df):,}")
            else:
                print("    WARNING: No Pattern column found!")
                continue

            if len(df) == 0:
                continue

            # Normalize chromosome
            chr_col = find_column(df, ['Chromosome', 'Chrom', 'Chr'])
            if chr_col is None:
                continue

            df['Chromosome'] = normalize_chromosome(df[chr_col])

            # Filter chromosomes
            if chromosomes:
                df = df[df['Chromosome'].isin(chromosomes)].copy()
                if len(df) == 0:
                    continue

            # Get position
            pos_col = find_column(df, ['Start', 'Position', 'Pos'])
            if pos_col is None:
                continue

            df = df.rename(columns={pos_col: 'Position'})

            # Add dose metadata
            df['Dose'] = dose
            df['Dose_Bin'] = dose_bin

            print(f"    Loaded: {len(df):,} radiation DBS")
            all_dbs.append(df)

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    if not all_dbs:
        print("\nERROR: No DBS data loaded!")
        return pd.DataFrame()

    combined = pd.concat(all_dbs, ignore_index=True)
    print(f"\n{'─' * 70}")
    print(f"Combined DBS from all files: {len(combined):,}")

    # Deduplication within dose (not dose bin)
    print("\nDeduplicating DBS WITHIN each dose...")

    combined['DBS_Key'] = (
        combined['Chromosome'] + ':' +
        combined['Position'].astype(str) + ':' +
        combined['Dose']
    )

    before_dedup = len(combined)
    combined = combined.drop_duplicates(subset='DBS_Key', keep='first')
    after_dedup = len(combined)

    print(f"  Before: {before_dedup:,} → After: {after_dedup:,}")
    print(f"  Removed: {before_dedup - after_dedup:,} duplicates within same dose")

    # Summary
    print("\nDBS counts by dose bin:")
    for bin_name in ['Low', 'High']:
        count = (combined['Dose_Bin'] == bin_name).sum()
        print(f"  {bin_name}: {count:,}")

    # Pattern distribution
    print("\nDBS by pattern:")
    for pattern in ['0T00', '00T0', '000T']:
        count = (combined['Pattern'] == pattern).sum()
        print(f"  {pattern}: {count:,}")
    multi = combined[combined['Pattern'].isin(['0TT0', '00TT', '0T0T', '0TTT'])]
    print(f"  Multi-timepoint: {len(multi):,}")

    return combined


# =============================================================================
# INV-DBS PAIRING WITH TEMPORAL CONCORDANCE
# =============================================================================

def find_inv_dbs_pairs(inv_df, dbs_df, window_size=10):
    """Find INV-DBS co-occurrence with TEMPORAL CONCORDANCE requirement.

    CRITICAL: Only pairs where:
    1. Same dose (INV dose == DBS dose)
    2. Same timepoint (INV timepoint present in DBS pattern)
    3. Spatial proximity (DBS within window of INV breakpoint)

    This ensures we're measuring true mutagenic coupling from the same repair event.
    """
    print(f"\n{'=' * 70}")
    print("FINDING INV-DBS PAIRS WITH TEMPORAL CONCORDANCE")
    print(f"Window: ±{window_size}bp")
    print("=" * 70)

    if 'Pattern' not in dbs_df.columns:
        print("ERROR: DBS data missing Pattern column - cannot check temporal concordance!")
        return {}

    results = {}

    for dose_bin in ['Low', 'High']:
        print(f"\n{'─' * 70}")
        print(f"{dose_bin.upper()} DOSE (doses {', '.join(DOSE_BINS[dose_bin])})")
        print('─' * 70)

        bin_invs = inv_df[inv_df['Dose_Bin'] == dose_bin].copy()
        bin_dbs = dbs_df[dbs_df['Dose_Bin'] == dose_bin].copy()

        mega_inv = bin_invs['Is_Mega'].sum()
        small_inv = (~bin_invs['Is_Mega']).sum()
        print(f"Inversions: {len(bin_invs):,} ({mega_inv:,} mega, {small_inv:,} small)")
        print(f"DBS mutations: {len(bin_dbs):,}")

        if len(bin_invs) == 0 or len(bin_dbs) == 0:
            print("  → Skipping (insufficient data)")
            results[dose_bin] = pd.DataFrame()
            continue

        concordant_pairs = []
        discordant_pairs = []

        # Group by chromosome for efficiency
        chroms_inv = set(bin_invs['Chromosome'].unique())
        chroms_dbs = set(bin_dbs['Chromosome'].unique())
        chroms_both = sorted(chroms_inv & chroms_dbs)

        print(f"Chromosomes with both: {len(chroms_both)}")

        for chrom in chroms_both:
            chrom_invs = bin_invs[bin_invs['Chromosome'] == chrom]
            chrom_dbs = bin_dbs[bin_dbs['Chromosome'] == chrom]

            for _, inv in chrom_invs.iterrows():
                sv_start = inv['SV_Start']
                sv_end = inv['SV_End']
                inv_dose = inv['Dose']
                inv_timepoint = inv['Timepoint']

                # Filter DBS to SAME DOSE
                dose_matched_dbs = chrom_dbs[chrom_dbs['Dose'] == inv_dose]

                if len(dose_matched_dbs) == 0:
                    continue

                dbs_positions = dose_matched_dbs['Position'].values
                dbs_patterns = dose_matched_dbs['Pattern'].values

                # Find DBS near breakpoints
                near_start = (
                    (dbs_positions >= sv_start - window_size) &
                    (dbs_positions <= sv_start + window_size)
                )

                near_end = (
                    (dbs_positions >= sv_end - window_size) &
                    (dbs_positions <= sv_end + window_size)
                )

                # Count CONCORDANT DBS (same timepoint)
                dbs_at_start_concordant = 0
                dbs_at_end_concordant = 0
                dbs_at_start_total = 0
                dbs_at_end_total = 0

                # Check start breakpoint
                start_indices = np.where(near_start)[0]
                for idx in start_indices:
                    dbs_at_start_total += 1
                    if timepoint_matches_pattern(inv_timepoint, dbs_patterns[idx]):
                        dbs_at_start_concordant += 1

                # Check end breakpoint
                end_indices = np.where(near_end)[0]
                for idx in end_indices:
                    dbs_at_end_total += 1
                    if timepoint_matches_pattern(inv_timepoint, dbs_patterns[idx]):
                        dbs_at_end_concordant += 1

                # Only count CONCORDANT pairs
                dbs_count_concordant = dbs_at_start_concordant + dbs_at_end_concordant
                dbs_count_total = dbs_at_start_total + dbs_at_end_total

                if dbs_count_concordant > 0:
                    concordant_pairs.append({
                        'Chromosome': chrom,
                        'SV_Start': sv_start,
                        'SV_End': sv_end,
                        'SV_Length': inv['SV_Length'],
                        'Is_Mega': inv['Is_Mega'],
                        'Dose': inv_dose,
                        'Timepoint': inv_timepoint,
                        'Cumulative_Dose_mGy': inv['Cumulative_Dose_mGy'],
                        'Dose_Bin': dose_bin,
                        'DBS_Count': dbs_count_concordant,
                        'DBS_at_Start': dbs_at_start_concordant,
                        'DBS_at_End': dbs_at_end_concordant,
                        'Concordant': True
                    })

                # Also track discordant for comparison (optional output)
                dbs_count_discordant = dbs_count_total - dbs_count_concordant
                if dbs_count_discordant > 0:
                    discordant_pairs.append({
                        'Chromosome': chrom,
                        'SV_Start': sv_start,
                        'SV_End': sv_end,
                        'SV_Length': inv['SV_Length'],
                        'Is_Mega': inv['Is_Mega'],
                        'Dose': inv_dose,
                        'Timepoint': inv_timepoint,
                        'Cumulative_Dose_mGy': inv['Cumulative_Dose_mGy'],
                        'Dose_Bin': dose_bin,
                        'DBS_Count': dbs_count_discordant,
                        'DBS_at_Start': dbs_at_start_total - dbs_at_start_concordant,
                        'DBS_at_End': dbs_at_end_total - dbs_at_end_concordant,
                        'Concordant': False
                    })

        # Results - only concordant pairs
        if concordant_pairs:
            result_df = pd.DataFrame(concordant_pairs)
            results[dose_bin] = result_df

            mega_pairs = result_df['Is_Mega'].sum()
            small_pairs = (~result_df['Is_Mega']).sum()

            print("\n  CONCORDANT INV-DBS pairs (same dose + timepoint):")
            print(f"    Total: {len(result_df):,}")
            print(f"    Mega: {mega_pairs:,}")
            print(f"    Small: {small_pairs:,}")
            print(f"    Total DBS at breakpoints: {result_df['DBS_Count'].sum():,}")

            if discordant_pairs:
                disc_df = pd.DataFrame(discordant_pairs)
                print(f"\n  DISCORDANT pairs (same dose, different timepoint): {len(disc_df):,}")
        else:
            print("\n  No CONCORDANT INV-DBS pairs found")
            results[dose_bin] = pd.DataFrame()

            if discordant_pairs:
                disc_df = pd.DataFrame(discordant_pairs)
                print(f"  (Found {len(disc_df):,} discordant pairs - excluded)")

    return results


# =============================================================================
# GENE EXTRACTION
# =============================================================================

def extract_genes(inv_df, inv_dbs_df):
    """Extract genes from inversions that have CONCORDANT DBS at breakpoints."""
    if len(inv_dbs_df) == 0:
        return pd.DataFrame()

    gene_col = find_column(inv_df, ['Gene_name', 'Gene', 'SYMBOL'])
    if gene_col is None:
        return pd.DataFrame()

    genes_list = []

    for _, pair in inv_dbs_df.iterrows():
        matches = inv_df[
            (inv_df['Chromosome'] == pair['Chromosome']) &
            (inv_df['SV_Start'] == pair['SV_Start']) &
            (inv_df['SV_End'] == pair['SV_End']) &
            (inv_df['Dose'] == pair['Dose']) &
            (inv_df['Timepoint'] == pair['Timepoint'])
        ]

        for _, inv in matches.iterrows():
            genes_str = str(inv.get(gene_col, ''))
            if pd.isna(genes_str) or genes_str in ['', 'nan', 'NA']:
                continue

            for sep in [';', '/']:
                genes_str = genes_str.replace(sep, '|')

            for gene in genes_str.split('|'):
                gene = gene.strip()
                if gene and gene not in ['', 'nan', 'NA']:
                    genes_list.append({
                        'Gene': gene,
                        'Chromosome': pair['Chromosome'],
                        'SV_Start': pair['SV_Start'],
                        'SV_End': pair['SV_End'],
                        'SV_Length': pair['SV_Length'],
                        'Is_Mega': pair['Is_Mega'],
                        'Dose': pair['Dose'],
                        'Timepoint': pair['Timepoint'],
                        'Dose_Bin': pair['Dose_Bin'],
                        'DBS_Count': pair['DBS_Count']
                    })

    if not genes_list:
        return pd.DataFrame()

    return pd.DataFrame(genes_list)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Dose-stratified INV-DBS analysis with TEMPORAL CONCORDANCE (v4)"
    )
    parser.add_argument('--annotsv-dir', required=True, help='Directory with AnnotSV files')
    parser.add_argument('--mutation-dir', required=True, help='Directory with mutation files')
    parser.add_argument('--output-dir', default='dose_analysis_v4', help='Output directory')
    parser.add_argument('--inv-size', type=int, default=0, help='Minimum inversion size (bp)')
    parser.add_argument('--max-inv-size', type=int, default=None, help='Maximum inversion size (bp)')
    parser.add_argument('--window', type=int, default=10, help='DBS window around breakpoints (bp)')
    parser.add_argument('--mega-threshold', type=int, default=50_000_000,
                        help='Threshold for Is_Mega flag (bp, default: 50Mb)')
    parser.add_argument('--chromosomes', nargs='+', help='Chromosomes to analyze')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Header
    print("=" * 70)
    print("DOSE-STRATIFIED INV-DBS ANALYSIS WITH TEMPORAL CONCORDANCE (v4)")
    print("=" * 70)
    print("\n*** IMPORTANT: Only CONCORDANT pairs are counted ***")
    print("*** INV and DBS must have SAME dose AND timepoint ***")
    print(f"\nDose bins: Low = {DOSE_BINS['Low']}, High = {DOSE_BINS['High']}")
    print(f"Window size: ±{args.window}bp")
    print(f"Mega threshold: ≥{args.mega_threshold/1e6:.0f}Mb")
    print(f"Output: {args.output_dir}/")

    # Load and process
    sv_df = load_annotsv_files(args.annotsv_dir, chromosomes=args.chromosomes)
    if len(sv_df) == 0:
        return 1

    inv_df = filter_inversions(sv_df, min_size=args.inv_size, max_size=args.max_inv_size,
                                mega_threshold=args.mega_threshold)
    if len(inv_df) == 0:
        return 1

    inv_df.to_csv(output_dir / 'inversions_all.csv', index=False)

    dbs_df = load_dbs_mutations(args.mutation_dir, chromosomes=args.chromosomes)
    if len(dbs_df) == 0:
        return 1

    inv_dbs_results = find_inv_dbs_pairs(inv_df, dbs_df, window_size=args.window)

    # Save results
    print(f"\n{'=' * 70}")
    print("SAVING RESULTS")
    print("=" * 70)

    summary_data = []

    for dose_bin, pairs_df in inv_dbs_results.items():
        bin_lower = dose_bin.lower()

        if len(pairs_df) == 0:
            print(f"\n{dose_bin}: No concordant INV-DBS pairs found")
            summary_data.append({
                'Dose_Bin': dose_bin,
                'INV_with_DBS_total': 0, 'INV_with_DBS_mega': 0, 'INV_with_DBS_small': 0,
                'Unique_genes_total': 0, 'Unique_genes_mega': 0, 'Unique_genes_small': 0
            })
            continue

        # Save pairs
        pairs_file = output_dir / f'inv_dbs_pairs_{bin_lower}.csv'
        pairs_df.to_csv(pairs_file, index=False)
        print(f"\n{dose_bin}:")
        print(f"  Pairs: {pairs_file}")

        # Extract genes
        genes_df = extract_genes(inv_df, pairs_df)

        if len(genes_df) > 0:
            genes_file = output_dir / f'genes_{bin_lower}_dose.csv'
            genes_df.to_csv(genes_file, index=False)
            print(f"  Genes: {genes_file}")

            # Calculate stats
            total_genes = genes_df['Gene'].nunique()
            mega_genes = genes_df[genes_df['Is_Mega']]['Gene'].nunique()
            small_genes = genes_df[~genes_df['Is_Mega']]['Gene'].nunique()

            print(f"  Unique genes: {total_genes:,} total ({mega_genes:,} mega, {small_genes:,} small)")
        else:
            total_genes = mega_genes = small_genes = 0

        mega_inv = pairs_df['Is_Mega'].sum()
        small_inv = (~pairs_df['Is_Mega']).sum()

        summary_data.append({
            'Dose_Bin': dose_bin,
            'INV_with_DBS_total': len(pairs_df),
            'INV_with_DBS_mega': mega_inv,
            'INV_with_DBS_small': small_inv,
            'Unique_genes_total': total_genes,
            'Unique_genes_mega': mega_genes,
            'Unique_genes_small': small_genes
        })

    # Summary
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / 'analysis_summary.csv'
    summary_df.to_csv(summary_file, index=False)

    print(f"\n{'=' * 70}")
    print("SUMMARY (CONCORDANT PAIRS ONLY)")
    print("=" * 70)
    print(summary_df.to_string(index=False))

    print(f"\n{'─' * 70}")
    print("NOTE: Only INV-DBS pairs with SAME dose AND timepoint are counted.")
    print("This ensures we're measuring true mutagenic coupling from the same repair event.")
    print('─' * 70)

    print(f"\nAll results saved to: {args.output_dir}/")
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    exit(main())
