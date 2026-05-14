#!/usr/bin/env python3
"""
Concordance Definition Comparison
===================================
Compares two concordance definitions on existing INV-DBS pairs:

  Loose (current): INV timepoint is present anywhere in DBS pattern
                   e.g. W2 INV + 00TT DBS → concordant

  Strict: DBS must be single-timepoint AND match INV timepoint exactly
          e.g. W2 INV + 00TT DBS → discordant
               W2 INV + 00T0 DBS → concordant

Usage:
    python check_concordance_definitions.py \
        --pairs figure_temporal_dynamics_v3_pairs.csv

Or re-derive pairs from raw data (same args as original script):
    python check_concordance_definitions.py \
        --annotsv-dir ./AnnotSV_files \
        --dbs-data ./merged_mutation_data_files/DBS \
        --window 10
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import re

# Single-timepoint radiation patterns only
SINGLE_TIMEPOINT_PATTERNS = {'0T00': 'W1', '00T0': 'W2', '000T': 'W3'}
RADIATION_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']
TIMEPOINT_TO_PATTERN_POS = {'W1': 1, 'W2': 2, 'W3': 3}


def timepoint_in_pattern(timepoint, pattern):
    """Loose: INV timepoint is present anywhere in DBS pattern."""
    if len(pattern) < 4:
        return False
    pos = TIMEPOINT_TO_PATTERN_POS.get(timepoint, 0)
    return pos > 0 and pattern[pos] != '0'


def timepoint_exact_match(timepoint, pattern):
    """Strict: DBS is single-timepoint AND matches INV timepoint exactly."""
    expected = {'W1': '0T00', 'W2': '00T0', 'W3': '000T'}.get(timepoint)
    return pattern == expected


def evaluate_concordance(pairs_df):
    """Apply both definitions and print comparison."""

    pairs_df = pairs_df.copy()
    pairs_df['Concordant_Loose'] = pairs_df.apply(
        lambda r: timepoint_in_pattern(r['INV_Timepoint'], r['DBS_Pattern']), axis=1
    )
    pairs_df['Concordant_Strict'] = pairs_df.apply(
        lambda r: timepoint_exact_match(r['INV_Timepoint'], r['DBS_Pattern']), axis=1
    )

    n = len(pairs_df)

    print("\n" + "=" * 70)
    print("CONCORDANCE DEFINITION COMPARISON")
    print("=" * 70)
    print(f"\nTotal pairs: {n}")

    # Show which pairs differ between definitions
    differs = pairs_df[pairs_df['Concordant_Loose'] != pairs_df['Concordant_Strict']]
    print(f"\nPairs where definitions disagree: {len(differs)}")
    if len(differs) > 0:
        print("\n  DBS patterns causing disagreement:")
        print(differs[['INV_Timepoint', 'DBS_Pattern',
                        'Concordant_Loose', 'Concordant_Strict']].to_string(index=False))

    print("\n" + "-" * 70)
    print(f"{'Metric':<35} {'Loose':>10} {'Strict':>10}")
    print("-" * 70)

    # Overall
    loose_overall = pairs_df['Concordant_Loose'].mean() * 100
    strict_overall = pairs_df['Concordant_Strict'].mean() * 100
    print(f"{'Overall concordance':<35} {loose_overall:>9.1f}% {strict_overall:>9.1f}%")

    # By timepoint
    for tp in ['W1', 'W2', 'W3']:
        tp_pairs = pairs_df[pairs_df['INV_Timepoint'] == tp]
        if len(tp_pairs) > 0:
            loose = tp_pairs['Concordant_Loose'].mean() * 100
            strict = tp_pairs['Concordant_Strict'].mean() * 100
            print(f"  {tp} (n={len(tp_pairs):<3}){'':<25} {loose:>9.1f}% {strict:>9.1f}%")

    # By size
    def get_size_cat(length):
        if length < 10_000:      return 'Small (<10kb)'
        elif length < 1_000_000: return 'Medium (10kb-1Mb)'
        elif length < 50_000_000:return 'Large (1-50Mb)'
        else:                    return 'Mega (≥50Mb)'

    pairs_df['Size_Cat'] = pairs_df['INV_Length'].apply(get_size_cat)
    print()
    for size in ['Small (<10kb)', 'Medium (10kb-1Mb)', 'Large (1-50Mb)', 'Mega (≥50Mb)']:
        s = pairs_df[pairs_df['Size_Cat'] == size]
        if len(s) > 0:
            loose = s['Concordant_Loose'].mean() * 100
            strict = s['Concordant_Strict'].mean() * 100
            print(f"  {size:<33} {loose:>9.1f}% {strict:>9.1f}%")

    print("-" * 70)
    print(f"\nRandom expectation (null): 33.3%")

    return pairs_df


# =============================================================================
# OPTIONAL: re-derive pairs from raw data (copied from original script)
# =============================================================================

def find_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None

def normalize_chromosome(chrom_series):
    chrom_series = chrom_series.astype(str)
    return chrom_series.apply(lambda x: x if x.startswith('chr') else f'chr{x}')

def load_annotsv_inversions(annotsv_dir, mega_threshold=50_000_000):
    annotsv_path = Path(annotsv_dir)
    all_files = sorted(list(set(
        list(annotsv_path.glob('d0_vs_d*_W*_annotated.tsv')) +
        list(annotsv_path.glob('d0_vs_d*_W*.tsv'))
    )))
    if not all_files:
        raise FileNotFoundError(f"No AnnotSV files found in {annotsv_dir}")

    all_invs = []
    for file_path in all_files:
        filename = file_path.name
        dose = next((d for d in 'ABCDE' if f'_vs_d{d}_' in filename), None)
        timepoint = next((w for w in ['W1','W2','W3'] if f'_{w}_' in filename or f'_{w}.' in filename), None)
        if not dose or not timepoint:
            continue
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        at_col = find_column(df, ['AnnotSV_type', 'Annotation_mode'])
        if at_col:
            df = df[df[at_col] == 'full']
        type_col = find_column(df, ['SV_type', 'Type', 'SV_Type'])
        if type_col:
            df = df[df[type_col] == 'INV']
        if len(df) == 0:
            continue
        chr_col = find_column(df, ['SV_chrom', 'Chromosome', 'Chrom', '#Chromosome'])
        start_col = find_column(df, ['SV_start', 'Start', 'start'])
        end_col = find_column(df, ['SV_end', 'End', 'end'])
        df['Chromosome'] = normalize_chromosome(df[chr_col])
        df = df.rename(columns={start_col: 'SV_Start', end_col: 'SV_End'})
        df['SV_Length'] = abs(df['SV_End'] - df['SV_Start'])
        df['Is_Mega'] = df['SV_Length'] >= mega_threshold
        df['Dose'] = dose
        df['Timepoint'] = timepoint
        all_invs.append(df[['Chromosome','SV_Start','SV_End','SV_Length','Is_Mega','Dose','Timepoint']])

    combined = pd.concat(all_invs, ignore_index=True)
    combined['INV_Key'] = (combined['Chromosome'] + ':' + combined['SV_Start'].astype(str) +
                           '-' + combined['SV_End'].astype(str) + ':' +
                           combined['Dose'] + ':' + combined['Timepoint'])
    return combined.drop_duplicates(subset='INV_Key', keep='first')


def load_dbs_mutations(mutation_dir):
    mutation_path = Path(mutation_dir)
    dbs_files = sorted(list(mutation_path.glob('DBS_dose_*_merged.csv')))
    if not dbs_files:
        raise FileNotFoundError(f"No DBS files found in {mutation_dir}")
    # Extract dose label from filenames like:
    #   DBS_dose_A_merged.csv   -> "A"
    #   DBS_dose_dA_merged.csv  -> "A"  (current convention from merge_annotation.py)
    dose_re = re.compile(r'dose_d?([A-E])_merged', re.IGNORECASE)

    all_dbs = []
    for file_path in dbs_files:
        m = dose_re.search(file_path.name)
        if not m:
            print(f"  skip {file_path.name}: no dose_[d]X_merged pattern found")
            continue
        dose = m.group(1).upper()
        df = pd.read_csv(file_path, low_memory=False)
        if 'Sample' in df.columns:
            df = df[~df['Sample'].str.contains('d0', case=False, na=False)]
        if 'Pattern' in df.columns:
            df = df[df['Pattern'].isin(RADIATION_PATTERNS)]
        chr_col = find_column(df, ['Chromosome', 'Chrom', 'Chr'])
        pos_col = find_column(df, ['Start', 'Position', 'Pos'])
        df['Chromosome'] = normalize_chromosome(df[chr_col])
        df = df.rename(columns={pos_col: 'Position'})
        df['Dose'] = dose
        all_dbs.append(df)
    return pd.concat(all_dbs, ignore_index=True)


def derive_pairs(annotsv_dir, dbs_dir, window=10):
    print("Loading inversions and DBS from raw data...")
    inv_df = load_annotsv_inversions(annotsv_dir)
    dbs_df = load_dbs_mutations(dbs_dir)
    pairs = []
    for chrom in sorted(set(inv_df['Chromosome']) & set(dbs_df['Chromosome'])):
        chrom_invs = inv_df[inv_df['Chromosome'] == chrom]
        chrom_dbs = dbs_df[dbs_df['Chromosome'] == chrom]
        for _, inv in chrom_invs.iterrows():
            dose_dbs = chrom_dbs[chrom_dbs['Dose'] == inv['Dose']]
            if len(dose_dbs) == 0:
                continue
            pos = dose_dbs['Position'].values
            for bp, bp_type in [(inv['SV_Start'], 'Start'), (inv['SV_End'], 'End')]:
                mask = (pos >= bp - window) & (pos <= bp + window)
                for i in np.where(mask)[0]:
                    pairs.append({
                        'INV_Timepoint': inv['Timepoint'],
                        'INV_Length': inv['SV_Length'],
                        'Is_Mega': inv['Is_Mega'],
                        'DBS_Pattern': dose_dbs['Pattern'].values[i],
                        'Breakpoint': bp_type
                    })
    pairs_df = pd.DataFrame(pairs)
    return pairs_df.drop_duplicates(keep='first')


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', help='Path to existing pairs CSV')
    parser.add_argument('--annotsv-dir', help='AnnotSV directory (if re-deriving pairs)')
    parser.add_argument('--dbs-data', help='DBS directory (if re-deriving pairs)')
    parser.add_argument('--window', type=int, default=10)
    args = parser.parse_args()

    if args.pairs:
        print(f"Loading pairs from: {args.pairs}")
        pairs_df = pd.read_csv(args.pairs)
    elif args.annotsv_dir and args.dbs_data:
        pairs_df = derive_pairs(args.annotsv_dir, args.dbs_data, args.window)
    else:
        parser.error("Provide either --pairs or both --annotsv-dir and --dbs-data")

    result_df = evaluate_concordance(pairs_df)
    out = (Path(args.pairs).stem + '_concordance_comparison.csv'
           if args.pairs else 'concordance_comparison.csv')
    result_df.to_csv(out, index=False)
    print(f"\n✓ Saved: {out}")


if __name__ == "__main__":
    main()