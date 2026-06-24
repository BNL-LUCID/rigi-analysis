#!/usr/bin/env python3
"""
SV Type-Specific Enrichment Decay Analysis
==========================================
Calculate mutation enrichment at multiple window sizes for each SV type separately.

This completes the puzzle:
- Figure 4A shows DBS decay for ALL SVs combined
- This analysis shows decay curves for INV, TRA, DUP, DEL separately
- Tests whether INV maintains strong enrichment at larger windows

Usage:
    python sv_type_specific_decay.py \
        --sv-catalog sv_temporal_catalog.csv \
        --mutations merged_mutation_data_files/DBS \
        --output sv_type_decay_analysis \
        --windows 10,25,50,100
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

try:
    import dask.dataframe as dd
    from dask.diagnostics import ProgressBar
    DASK_AVAILABLE = True
    print("✓ Dask available - using parallel processing")
except ImportError:
    DASK_AVAILABLE = False
    print("⚠ Dask not available - falling back to pandas (will be slow)")

try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False

# Constants
RADIATION_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '0T0T', '00TT', '0TTT']
CONTROL_PATTERNS = ['0C00', '00C0', '000C', '0CC0', '0C0C', '00CC', '0CCC']

CHROM_SIZES = {
    '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
    '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
    '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
    '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
    '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
    '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
}

GENOME_SIZE = sum(CHROM_SIZES.values())
VALID_CHROMS = set(CHROM_SIZES.keys())


def load_sv_catalog(sv_file):
    """Load SV catalog and extract breakpoints by type."""
    print(f"\nLoading SV catalog: {sv_file}")
    
    df = pd.read_csv(sv_file, low_memory=False)
    print(f"  Total SVs: {len(df):,}")
    
    # Filter to radiation SVs with valid chromosomes
    df = df[df['Pattern'].isin(RADIATION_PATTERNS)].copy()
    df['Chrom'] = df['Chrom'].astype(str).str.replace('chr', '', regex=False)
    df = df[df['Chrom'].isin(VALID_CHROMS)]
    
    print(f"  Radiation SVs in valid chroms: {len(df):,}")
    
    # Count by SV type
    sv_counts = df['SV_Type'].value_counts()
    print("\n  SVs by type:")
    for sv_type, count in sv_counts.items():
        print(f"    {sv_type}: {count:,}")
    
    return df


def extract_breakpoints_by_type(sv_df):
    """Extract breakpoint coordinates for each SV type."""
    
    breakpoints_by_type = {}
    
    for sv_type in ['INV', 'TRA', 'DUP', 'DEL', 'INS']:
        sv_subset = sv_df[sv_df['SV_Type'] == sv_type].copy()
        
        if len(sv_subset) == 0:
            continue
        
        breakpoints = []
        
        for _, row in sv_subset.iterrows():
            chrom = str(row['Chrom'])
            
            # Add start breakpoint
            if pd.notna(row.get('Start')):
                breakpoints.append({
                    'Chrom': chrom,
                    'Pos': int(row['Start']),
                    'SV_Type': sv_type
                })
            
            # Add end breakpoint
            if pd.notna(row.get('End')):
                breakpoints.append({
                    'Chrom': chrom,
                    'Pos': int(row['End']),
                    'SV_Type': sv_type
                })
        
        bp_df = pd.DataFrame(breakpoints)
        bp_df = bp_df.drop_duplicates()
        
        breakpoints_by_type[sv_type] = bp_df
        print(f"  {sv_type}: {len(bp_df):,} unique breakpoints")
    
    return breakpoints_by_type


def load_mutations(mut_dir):
    """Load mutation data using dask for efficiency."""
    print(f"\nLoading mutations from: {mut_dir}")

    # Non-recursive on purpose: pass a single mutation type's dir
    # (e.g. merged_data/DBS/), not an aggregate root. Mixing types
    # dilutes signal -- SNV's much larger genomic spread drowns out
    # type-specific clustering near SV breakpoints.
    mut_files = sorted(Path(mut_dir).glob("*.csv"))
    print(f"  Found {len(mut_files)} CSV files")

    if len(mut_files) == 0:
        raise FileNotFoundError(f"No CSV files found in {mut_dir}")
    
    # Read first file to detect column names
    sample_df = pd.read_csv(mut_files[0], nrows=5)
    print(f"  Detected columns: {sample_df.columns.tolist()}")
    
    # Find chromosome column (could be 'Chrom', 'Chromosome', '#CHROM', etc.)
    chrom_col = None
    pos_col = None
    pattern_col = None
    
    for col in sample_df.columns:
        col_lower = col.lower()
        if 'chrom' in col_lower and chrom_col is None:
            chrom_col = col
        if ('pos' in col_lower or 'start' in col_lower) and pos_col is None:
            pos_col = col
        if 'pattern' in col_lower and pattern_col is None:
            pattern_col = col
    
    if chrom_col is None:
        raise ValueError(f"Could not find chromosome column in {mut_files[0]}")
    if pos_col is None:
        raise ValueError(f"Could not find position column in {mut_files[0]}")
    if pattern_col is None:
        raise ValueError(f"Could not find pattern column in {mut_files[0]}")
    
    print(f"  Using columns: Chrom='{chrom_col}', Pos='{pos_col}', Pattern='{pattern_col}'")
    
    if DASK_AVAILABLE:
        print("  Using Dask for efficient loading...")
        
        # Only read the columns we actually need. This avoids dtype-inference
        # mismatches across files (e.g. Indel_Type is NaN/float in DBS files
        # but string in ID files, which crashes dd.read_csv if loaded).
        dtype_dict = {
            chrom_col: 'object',
            pattern_col: 'object',
        }

        # Pass the explicit file list rather than a non-recursive glob string,
        # so subdir layouts (merged_data/{TYPE}/*.csv) are handled too.
        ddf = dd.read_csv(
            [str(f) for f in mut_files],
            usecols=[chrom_col, pos_col, pattern_col],
            dtype=dtype_dict,
            blocksize='64MB',
            assume_missing=True,
        )
        
        # Rename to standard columns
        ddf = ddf.rename(columns={chrom_col: 'Chrom', pos_col: 'Pos', pattern_col: 'Pattern'})
        
        # Filter and clean
        ddf['Chrom'] = ddf['Chrom'].astype(str).str.replace('chr', '', regex=False)
        
        # Filter to valid chromosomes
        ddf = ddf[ddf['Chrom'].isin(list(VALID_CHROMS))]
        
        # Filter to radiation and control patterns
        rad_mask = ddf['Pattern'].isin(RADIATION_PATTERNS)
        ctrl_mask = ddf['Pattern'].isin(CONTROL_PATTERNS)
        ddf = ddf[rad_mask | ctrl_mask]
        
        # Keep only needed columns
        ddf = ddf[['Chrom', 'Pos', 'Pattern']]
        
        # Compute to pandas (only filtered data)
        print("  Computing filtered dataframe...")
        with ProgressBar():
            mut_df = ddf.compute()
        
    else:
        print("  Using pandas (this may be slow for large files)...")
        dfs = []
        for f in mut_files:
            print(f"    Loading {f.name}...")
            df = pd.read_csv(
                f,
                usecols=[chrom_col, pos_col, pattern_col],
                low_memory=False,
            )

            # Rename columns
            df = df.rename(columns={chrom_col: 'Chrom', pos_col: 'Pos', pattern_col: 'Pattern'})
            dfs.append(df[['Chrom', 'Pos', 'Pattern']])
        
        mut_df = pd.concat(dfs, ignore_index=True)
        
        # Filter and clean
        mut_df['Chrom'] = mut_df['Chrom'].astype(str).str.replace('chr', '', regex=False)
        mut_df = mut_df[mut_df['Chrom'].isin(VALID_CHROMS)]
        
        # Filter to radiation and control patterns
        rad_mask = mut_df['Pattern'].isin(RADIATION_PATTERNS)
        ctrl_mask = mut_df['Pattern'].isin(CONTROL_PATTERNS)
        mut_df = mut_df[rad_mask | ctrl_mask].copy()
    
    print(f"  Total mutations: {len(mut_df):,}")
    
    rad_total = mut_df['Pattern'].isin(RADIATION_PATTERNS).sum()
    ctrl_total = mut_df['Pattern'].isin(CONTROL_PATTERNS).sum()
    
    print(f"    Radiation: {rad_total:,}")
    print(f"    Control: {ctrl_total:,}")
    
    return mut_df


def count_mutations_near_breakpoints_chrom(chrom, breakpoints_df, mutations_df, windows):
    """Count mutations near breakpoints for one chromosome - optimized version."""
    
    # Filter to this chromosome
    bp_chrom = breakpoints_df[breakpoints_df['Chrom'] == chrom].copy()
    mut_chrom = mutations_df[mutations_df['Chrom'] == chrom].copy()
    
    if len(bp_chrom) == 0 or len(mut_chrom) == 0:
        return []
    
    # Sort for efficiency
    bp_chrom = bp_chrom.sort_values('Pos')
    mut_chrom = mut_chrom.sort_values('Pos')
    
    bp_positions = bp_chrom['Pos'].values
    mut_positions = mut_chrom['Pos'].values
    mut_patterns = mut_chrom['Pattern'].values
    
    results = []
    
    # For each window size
    for window in windows:
        rad_total = 0
        ctrl_total = 0
        
        # For each breakpoint
        for bp_pos in bp_positions:
            window_start = bp_pos - window
            window_end = bp_pos + window
            
            # Use searchsorted for fast range query
            start_idx = np.searchsorted(mut_positions, window_start, side='left')
            end_idx = np.searchsorted(mut_positions, window_end, side='right')
            
            if start_idx < end_idx:
                patterns_in_window = mut_patterns[start_idx:end_idx]
                
                rad_count = sum(p in RADIATION_PATTERNS for p in patterns_in_window)
                ctrl_count = sum(p in CONTROL_PATTERNS for p in patterns_in_window)
                
                rad_total += rad_count
                ctrl_total += ctrl_count
        
        if rad_total > 0 or ctrl_total > 0:
            results.append({
                'Chrom': chrom,
                'Window': window,
                'Rad_Count': rad_total,
                'Ctrl_Count': ctrl_total
            })
    
    return results


def analyze_sv_type(sv_type, breakpoints_df, mutations_df, windows, n_jobs=-1):
    """Analyze one SV type across all windows."""
    
    print(f"\n  Analyzing {sv_type}...")
    print(f"    Breakpoints: {len(breakpoints_df):,}")
    
    # Get chromosomes with both breakpoints and mutations
    chroms_with_bp = set(breakpoints_df['Chrom'].unique())
    chroms_with_mut = set(mutations_df['Chrom'].unique())
    chroms_to_process = sorted(chroms_with_bp & chroms_with_mut)
    
    print(f"    Chromosomes: {len(chroms_to_process)}")
    
    # Process chromosomes in parallel
    if JOBLIB_AVAILABLE and n_jobs != 1:
        results = Parallel(n_jobs=n_jobs)(
            delayed(count_mutations_near_breakpoints_chrom)(
                chrom, breakpoints_df, mutations_df, windows
            ) for chrom in chroms_to_process
        )
        results = [item for sublist in results for item in sublist]
    else:
        results = []
        for chrom in chroms_to_process:
            chrom_results = count_mutations_near_breakpoints_chrom(
                chrom, breakpoints_df, mutations_df, windows
            )
            results.extend(chrom_results)
    
    # Aggregate results
    if len(results) == 0:
        return None
    
    results_df = pd.DataFrame(results)
    
    # Sum by window
    summary = results_df.groupby('Window').agg({
        'Rad_Count': 'sum',
        'Ctrl_Count': 'sum'
    }).reset_index()
    
    summary['SV_Type'] = sv_type
    summary['N_Breakpoints'] = len(breakpoints_df)
    
    # Calculate enrichment
    rad_total = mutations_df['Pattern'].isin(RADIATION_PATTERNS).sum()
    ctrl_total = mutations_df['Pattern'].isin(CONTROL_PATTERNS).sum()
    
    rad_density = rad_total / GENOME_SIZE
    ctrl_density = ctrl_total / GENOME_SIZE
    
    for idx, row in summary.iterrows():
        window = row['Window']
        n_bp = row['N_Breakpoints']
        window_area = 2 * window * n_bp
        
        rad_expected = rad_density * window_area
        ctrl_expected = ctrl_density * window_area
        
        summary.at[idx, 'Rad_Expected'] = rad_expected
        summary.at[idx, 'Ctrl_Expected'] = ctrl_expected
        summary.at[idx, 'Rad_Enrichment'] = row['Rad_Count'] / rad_expected if rad_expected > 0 else 0
        summary.at[idx, 'Ctrl_Enrichment'] = row['Ctrl_Count'] / ctrl_expected if ctrl_expected > 0 else 0
        summary.at[idx, 'Enrichment_Ratio'] = summary.at[idx, 'Rad_Enrichment'] / summary.at[idx, 'Ctrl_Enrichment'] if summary.at[idx, 'Ctrl_Enrichment'] > 0 else 0
    
    print(f"    Results:")
    for _, row in summary.iterrows():
        print(f"      {int(row['Window'])}bp: Rad={row['Rad_Enrichment']:.2f}x, Ctrl={row['Ctrl_Enrichment']:.2f}x, Ratio={row['Enrichment_Ratio']:.2f}x")
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description='SV type-specific enrichment decay analysis'
    )
    parser.add_argument('--sv-catalog', required=True, help='SV temporal catalog CSV')
    parser.add_argument('--mutations', required=True, help='Directory with mutation CSV files')
    parser.add_argument('--output', default='sv_type_decay_analysis', help='Output directory')
    parser.add_argument('--windows', default='10,25,50,100', help='Comma-separated window sizes')
    parser.add_argument('--n-jobs', type=int, default=-1, help='Parallel jobs')
    
    args = parser.parse_args()
    
    # Parse windows
    windows = [int(w) for w in args.windows.split(',')]
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("SV TYPE-SPECIFIC ENRICHMENT DECAY ANALYSIS")
    print("="*80)
    print(f"Windows: {windows} bp")
    
    # Load data
    sv_df = load_sv_catalog(args.sv_catalog)
    breakpoints_by_type = extract_breakpoints_by_type(sv_df)
    mutations_df = load_mutations(args.mutations)
    
    # Analyze each SV type
    print("\n" + "="*80)
    print("ANALYZING EACH SV TYPE")
    print("="*80)
    
    all_results = []
    
    for sv_type in ['INV', 'TRA', 'DUP', 'DEL']:
        if sv_type not in breakpoints_by_type:
            print(f"\n  Skipping {sv_type} (no breakpoints)")
            continue
        
        result = analyze_sv_type(
            sv_type,
            breakpoints_by_type[sv_type],
            mutations_df,
            windows,
            args.n_jobs
        )
        
        if result is not None:
            all_results.append(result)
    
    # Combine and save
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        output_file = output_dir / 'sv_type_enrichment_decay.csv'
        combined_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved: {output_file}")
        
        # Print summary
        print("\n" + "="*80)
        print("SUMMARY: ENRICHMENT DECAY BY SV TYPE")
        print("="*80)
        
        for sv_type in ['INV', 'TRA', 'DUP', 'DEL']:
            subset = combined_df[combined_df['SV_Type'] == sv_type]
            if len(subset) > 0:
                print(f"\n{sv_type}:")
                print(f"  Window   Rad_Enrich   Ctrl_Enrich   Ratio")
                print(f"  -------  -----------  ------------  ------")
                for _, row in subset.iterrows():
                    print(f"  {int(row['Window']):>5}bp  {row['Rad_Enrichment']:>10.2f}x  {row['Ctrl_Enrichment']:>11.2f}x  {row['Enrichment_Ratio']:>5.2f}x")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()