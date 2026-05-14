#!/usr/bin/env python3
"""
SV-Mutation Positional Correlation Analysis (Dask Version)
===========================================================
Test whether mutations cluster near SV breakpoints using Dask for large files.

Hypothesis: Radiation → DSBs → SVs → Error-prone repair → Nearby mutations

Optimizations:
- Dask for lazy evaluation and chunked processing
- Early filtering to radiation patterns only
- Chromosome-parallel processing
- Memory-efficient joins

Usage:
    python sv_mutation_correlation_dask.py \
        --sv-catalog sv_temporal_catalog.csv \
        --mutations ./merged_mutation_data_files/ \
        --output sv_mutation_correlation \
        --windows 1000,5000,10000,50000 \
        --plot
"""

import pandas as pd
import numpy as np
import time
from pathlib import Path
import argparse
from scipy import stats
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Try to import dask
try:
    import dask.dataframe as dd
    from dask.diagnostics import ProgressBar
    DASK_AVAILABLE = True
    print("✓ Dask available - using parallel processing")
except ImportError:
    DASK_AVAILABLE = False
    print("⚠ Dask not available - falling back to pandas")

# Try to import joblib for parallel processing
try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
    print("✓ Joblib available - using parallel chromosome processing")
except ImportError:
    JOBLIB_AVAILABLE = False
    print("⚠ Joblib not available - sequential processing")

# Try to import joblib for parallel processing
try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
    print("✓ Joblib available - using parallel chromosome processing")
except ImportError:
    JOBLIB_AVAILABLE = False
    print("⚠ Joblib not available - sequential processing")


# =============================================================================
# CONSTANTS
# =============================================================================

# Pattern definitions - ONLY radiation patterns for filtering
RADIATION_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '0T0T', '00TT', '0TTT']
CONTROL_PATTERNS = ['0C00', '00C0', '000C', '0CC0', '0C0C', '00CC', '0CCC']

# Pattern to timepoint mapping
PATTERN_TIMEPOINTS = {
    '0T00': ['W1'], '0C00': ['W1'],
    '00T0': ['W2'], '00C0': ['W2'],
    '000T': ['W3'], '000C': ['W3'],
    '0TT0': ['W1', 'W2'], '0CC0': ['W1', 'W2'],
    '0T0T': ['W1', 'W3'], '0C0C': ['W1', 'W3'],
    '00TT': ['W2', 'W3'], '00CC': ['W2', 'W3'],
    '0TTT': ['W1', 'W2', 'W3'], '0CCC': ['W1', 'W2', 'W3']
}

# Default window sizes (bp)
DEFAULT_WINDOWS = [1000, 5000, 10000, 50000]

# Chromosome sizes (hg38)
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


# =============================================================================
# DATA LOADING WITH DASK
# =============================================================================

def load_sv_catalog(filepath):
    """Load SV catalog with breakpoint positions."""
    print(f"Loading SV catalog: {filepath}")
    df = pd.read_csv(filepath)
    
    # Standardize column names
    col_map = {
        'SV_chrom': 'Chrom', 'Chrom': 'Chrom', 'Chr': 'Chrom',
        'SV_start': 'Start', 'Start': 'Start',
        'SV_end': 'End', 'End': 'End',
        'SV_type': 'SV_Type', 'SV_Type': 'SV_Type',
        'Pattern': 'Pattern', 'Dose': 'Dose'
    }
    
    for old, new in col_map.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})
    
    # Clean chromosome names
    if 'Chrom' in df.columns:
        df['Chrom'] = df['Chrom'].astype(str).str.replace('chr', '')
        df = df[df['Chrom'].isin(VALID_CHROMS)]
    
    # Filter to radiation patterns ONLY
    if 'Pattern' in df.columns:
        df = df[df['Pattern'].isin(RADIATION_PATTERNS)].copy()
    
    print(f"  Loaded {len(df):,} radiation SVs")
    
    return df


def load_mutations_dask(filepath):
    """
    Load mutation files with Dask for memory efficiency.
    Uses usecols to only load needed columns and ProgressBar for tracking.
    """
    filepath = Path(filepath)
    
    all_patterns = set(RADIATION_PATTERNS + CONTROL_PATTERNS)
    
    if filepath.is_dir():
        print(f"Loading mutation files from directory: {filepath}")

        # Non-recursive on purpose: pass a single mutation type's dir
        # (e.g. merged_data/DBS/), not an aggregate root. Mixing types
        # dilutes signal — SNV's much larger genomic spread drowns out
        # type-specific clustering near SV breakpoints.
        csv_files = sorted(filepath.glob('*.csv'))
        print(f"  Found {len(csv_files)} CSV files")
        
        if DASK_AVAILABLE:
            # Use Dask with usecols for efficiency
            dfs = []
            
            for f in csv_files:
                print(f"  Loading {f.name}...")
                
                try:
                    # Load columns we need - include End for indels
                    ddf = dd.read_csv(
                        str(f),
                        usecols=['Chromosome', 'Start', 'End', 'Pattern'],
                        dtype={'Chromosome': str, 'Start': float, 'End': float, 'Pattern': str},
                        blocksize='64MB',
                        assume_missing=True
                    )
                    
                    # Filter to valid patterns
                    ddf = ddf[ddf['Pattern'].isin(all_patterns)]
                    
                    # Clean chromosome
                    ddf['Chrom'] = ddf['Chromosome'].str.replace('chr', '')
                    ddf = ddf[ddf['Chrom'].isin(VALID_CHROMS)]
                    
                    dfs.append(ddf[['Chrom', 'Start', 'End', 'Pattern']])
                    
                except Exception as e:
                    print(f"    ERROR: {e}")
            
            if dfs:
                print(f"\n  Computing combined dataframe...")
                combined_ddf = dd.concat(dfs)
                
                with ProgressBar():
                    df = combined_ddf.compute()
                
                # Convert positions to int
                df['Start'] = pd.to_numeric(df['Start'], errors='coerce').fillna(0).astype(int)
                df['End'] = pd.to_numeric(df['End'], errors='coerce').fillna(0).astype(int)
                
                # For SNVs where End might be missing or same as Start
                df.loc[df['End'] == 0, 'End'] = df.loc[df['End'] == 0, 'Start'] + 1
                
                print(f"\n  Total filtered mutations: {len(df):,}")
            else:
                raise ValueError("No mutation files loaded")
        
        else:
            # Pandas fallback - file by file
            all_dfs = []
            total_raw = 0
            total_filtered = 0
            
            for i, f in enumerate(csv_files):
                print(f"  [{i+1}/{len(csv_files)}] {f.name}...", end=' ', flush=True)
                
                try:
                    chunk_df = pd.read_csv(
                        f,
                        usecols=['Chromosome', 'Start', 'End', 'Pattern'],
                        dtype={'Chromosome': str, 'Pattern': str},
                        low_memory=False
                    )
                    raw_count = len(chunk_df)
                    total_raw += raw_count
                    
                    chunk_df = chunk_df[chunk_df['Pattern'].isin(all_patterns)]
                    chunk_df['Chrom'] = chunk_df['Chromosome'].str.replace('chr', '')
                    chunk_df = chunk_df[chunk_df['Chrom'].isin(VALID_CHROMS)]
                    chunk_df['Start'] = pd.to_numeric(chunk_df['Start'], errors='coerce').fillna(0).astype(int)
                    chunk_df['End'] = pd.to_numeric(chunk_df['End'], errors='coerce').fillna(0).astype(int)
                    chunk_df.loc[chunk_df['End'] == 0, 'End'] = chunk_df.loc[chunk_df['End'] == 0, 'Start'] + 1
                    
                    filtered_count = len(chunk_df)
                    total_filtered += filtered_count
                    
                    all_dfs.append(chunk_df[['Chrom', 'Start', 'End', 'Pattern']])
                    print(f"{raw_count:,} → {filtered_count:,}")
                    
                except Exception as e:
                    print(f"ERROR: {e}")
            
            df = pd.concat(all_dfs, ignore_index=True)
            print(f"\n  Total filtered mutations: {len(df):,}")
    
    else:
        # Single file
        print(f"Loading single mutation file: {filepath}")
        
        if DASK_AVAILABLE:
            ddf = dd.read_csv(
                str(filepath),
                usecols=['Chromosome', 'Start', 'End', 'Pattern'],
                dtype={'Chromosome': str, 'Start': float, 'End': float, 'Pattern': str},
                blocksize='64MB'
            )
            
            ddf = ddf[ddf['Pattern'].isin(all_patterns)]
            ddf['Chrom'] = ddf['Chromosome'].str.replace('chr', '')
            ddf = ddf[ddf['Chrom'].isin(VALID_CHROMS)]
            
            with ProgressBar():
                df = ddf[['Chrom', 'Start', 'End', 'Pattern']].compute()
            
            df['Start'] = pd.to_numeric(df['Start'], errors='coerce').fillna(0).astype(int)
            df['End'] = pd.to_numeric(df['End'], errors='coerce').fillna(0).astype(int)
            df.loc[df['End'] == 0, 'End'] = df.loc[df['End'] == 0, 'Start'] + 1
        else:
            df = pd.read_csv(
                filepath,
                usecols=['Chromosome', 'Start', 'End', 'Pattern'],
                dtype={'Chromosome': str, 'Pattern': str},
                low_memory=False
            )
            df['Chrom'] = df['Chromosome'].str.replace('chr', '')
            df = df[df['Pattern'].isin(all_patterns)]
            df = df[df['Chrom'].isin(VALID_CHROMS)]
            df['Start'] = pd.to_numeric(df['Start'], errors='coerce').fillna(0).astype(int)
            df['End'] = pd.to_numeric(df['End'], errors='coerce').fillna(0).astype(int)
            df.loc[df['End'] == 0, 'End'] = df.loc[df['End'] == 0, 'Start'] + 1
            df = df[['Chrom', 'Start', 'End', 'Pattern']]
    
    # Summary
    rad_count = df['Pattern'].isin(RADIATION_PATTERNS).sum()
    ctrl_count = df['Pattern'].isin(CONTROL_PATTERNS).sum()
    print(f"    Radiation (T): {rad_count:,}")
    print(f"    Control (C):   {ctrl_count:,}")
    
    return df


# =============================================================================
# BREAKPOINT EXTRACTION
# =============================================================================

def extract_breakpoints(sv_df):
    """Extract individual breakpoints from SVs."""
    print("\nExtracting SV breakpoints...")
    
    breakpoints = []
    
    for _, row in sv_df.iterrows():
        chrom = str(row.get('Chrom', ''))
        start = int(row.get('Start', 0))
        end = int(row.get('End', 0))
        sv_type = row.get('SV_Type', 'UNK')
        pattern = row.get('Pattern', '')
        dose = row.get('Dose', '')
        
        # Breakpoint 1 (start)
        breakpoints.append({
            'Chrom': chrom,
            'Position': start,
            'SV_Type': sv_type,
            'Pattern': pattern,
            'Dose': dose
        })
        
        # Breakpoint 2 (end) - only if different
        if end != start and end > 0:
            breakpoints.append({
                'Chrom': chrom,
                'Position': end,
                'SV_Type': sv_type,
                'Pattern': pattern,
                'Dose': dose
            })
    
    bp_df = pd.DataFrame(breakpoints)
    
    # Remove duplicates
    bp_df = bp_df.drop_duplicates(subset=['Chrom', 'Position', 'Pattern'])
    
    print(f"  Unique breakpoints: {len(bp_df):,}")
    
    return bp_df


# =============================================================================
# EFFICIENT CHROMOSOME-PARALLEL ANALYSIS
# =============================================================================

def analyze_chromosome_optimized(chrom, bp_chrom, mut_chrom, window_sizes):
    """
    Analyze mutations near breakpoints for a single chromosome.
    OPTIMIZED: Uses sorted arrays and binary search for efficiency.
    """
    results = []
    
    if len(bp_chrom) == 0 or len(mut_chrom) == 0:
        return results
    
    # Sort mutations by Start for efficient searching
    mut_sorted = mut_chrom.sort_values('Start').reset_index(drop=True)
    mut_starts = mut_sorted['Start'].values.astype(np.int64)
    mut_ends = mut_sorted['End'].values.astype(np.int64)
    mut_patterns = mut_sorted['Pattern'].values
    
    # Pre-compute pattern classification (once for all)
    is_rad = np.array([p in RADIATION_PATTERNS for p in mut_patterns])
    is_ctrl = np.array([p in CONTROL_PATTERNS for p in mut_patterns])
    
    # Get breakpoint data
    bp_positions = bp_chrom['Position'].values.astype(np.int64)
    bp_sv_types = bp_chrom['SV_Type'].values
    bp_sv_patterns = bp_chrom['Pattern'].values
    
    # Use largest window to find max search range
    max_window = max(window_sizes)
    
    for i in range(len(bp_positions)):
        bp_pos = bp_positions[i]
        sv_type = bp_sv_types[i]
        sv_pattern = bp_sv_patterns[i]
        
        # Find candidate mutations using binary search on largest window
        # Mutations that could overlap: Start <= bp_pos + max_window AND End >= bp_pos - max_window
        # Use Start for left bound: find first mutation where Start > bp_pos + max_window
        right_idx = np.searchsorted(mut_starts, bp_pos + max_window, side='right')
        
        # Quick skip if no candidates
        if right_idx == 0:
            continue
        
        # Get candidate slice
        cand_starts = mut_starts[:right_idx]
        cand_ends = mut_ends[:right_idx]
        cand_patterns = mut_patterns[:right_idx]
        cand_is_rad = is_rad[:right_idx]
        cand_is_ctrl = is_ctrl[:right_idx]
        
        # Pre-filter: only keep mutations where End >= bp_pos - max_window
        valid_mask = cand_ends >= (bp_pos - max_window)
        
        if not valid_mask.any():
            continue
        
        # Apply filter
        cand_starts = cand_starts[valid_mask]
        cand_ends = cand_ends[valid_mask]
        cand_patterns = cand_patterns[valid_mask]
        cand_is_rad = cand_is_rad[valid_mask]
        cand_is_ctrl = cand_is_ctrl[valid_mask]
        
        sv_timepoints = set(PATTERN_TIMEPOINTS.get(sv_pattern, []))
        
        for window in window_sizes:
            win_start = bp_pos - window
            win_end = bp_pos + window
            
            # Check overlap: NOT (end < win_start OR start > win_end)
            overlaps = ~((cand_ends < win_start) | (cand_starts > win_end))
            
            if not overlaps.any():
                continue
            
            # Count
            rad_count = (overlaps & cand_is_rad).sum()
            ctrl_count = (overlaps & cand_is_ctrl).sum()
            
            # Concordance and distance for radiation
            concordant = 0
            mean_dist = np.nan
            min_dist = np.nan
            
            if rad_count > 0:
                rad_mask = overlaps & cand_is_rad
                rad_patterns = cand_patterns[rad_mask]
                
                for pat in rad_patterns:
                    mut_timepoints = set(PATTERN_TIMEPOINTS.get(pat, []))
                    if sv_timepoints & mut_timepoints:
                        concordant += 1
                
                # Distances
                rad_starts = cand_starts[rad_mask]
                rad_ends = cand_ends[rad_mask]
                
                distances = np.where(
                    (rad_starts <= bp_pos) & (bp_pos <= rad_ends),
                    0,
                    np.minimum(np.abs(rad_starts - bp_pos), np.abs(rad_ends - bp_pos))
                )
                mean_dist = distances.mean()
                min_dist = distances.min()
            
            results.append({
                'Chrom': chrom,
                'BP_Position': int(bp_pos),
                'SV_Type': sv_type,
                'SV_Pattern': sv_pattern,
                'Window': window,
                'Rad_Count': int(rad_count),
                'Ctrl_Count': int(ctrl_count),
                'Concordant': int(concordant),
                'Mean_Distance': mean_dist,
                'Min_Distance': min_dist
            })
    
    return results


def count_mutations_parallel(breakpoints_df, mutations_df, window_sizes, n_jobs=-1):
    """
    Count mutations near breakpoints using parallel processing across chromosomes.
    
    Args:
        n_jobs: Number of parallel jobs (-1 = all cores)
    """
    print("\nCounting mutations near breakpoints...")
    print(f"  Window sizes: {window_sizes}")
    
    # Group by chromosome
    bp_by_chrom = {c: g for c, g in breakpoints_df.groupby('Chrom')}
    mut_by_chrom = {c: g for c, g in mutations_df.groupby('Chrom')}
    
    all_chroms = sorted(set(bp_by_chrom.keys()) & set(mut_by_chrom.keys()), 
                        key=lambda x: int(x) if x.isdigit() else 100)
    
    print(f"  Chromosomes with both SVs and mutations: {len(all_chroms)}")
    
    # Count totals
    total_bp = sum(len(bp_by_chrom[c]) for c in all_chroms)
    total_mut = sum(len(mut_by_chrom[c]) for c in all_chroms)
    
    print(f"  Total breakpoints: {total_bp:,}")
    print(f"  Total mutations:   {total_mut:,}")
    
    # Show per-chromosome stats
    print(f"\n  Per-chromosome breakdown:")
    for chrom in all_chroms:
        n_bp = len(bp_by_chrom[chrom])
        n_mut = len(mut_by_chrom[chrom])
        print(f"    Chr {chrom:>2}: {n_bp:>6,} BP × {n_mut:>10,} mutations")
    
    if JOBLIB_AVAILABLE and len(all_chroms) > 1:
        print(f"\n  Processing {len(all_chroms)} chromosomes in PARALLEL (n_jobs={n_jobs})...")
        start_time = time.time()
        
        # Parallel processing with joblib
        results_list = Parallel(n_jobs=n_jobs, verbose=10, backend='loky')(
            delayed(analyze_chromosome_optimized)(
                chrom, 
                bp_by_chrom[chrom], 
                mut_by_chrom[chrom], 
                window_sizes
            )
            for chrom in all_chroms
        )
        
        elapsed = time.time() - start_time
        print(f"\n  Parallel processing completed in {elapsed:.1f}s")
        
        # Flatten results
        all_results = []
        for chrom_results in results_list:
            all_results.extend(chrom_results)
    
    else:
        # Sequential processing with progress
        print(f"\n  Processing chromosomes SEQUENTIALLY...")
        all_results = []
        
        for chrom in all_chroms:
            bp_chrom = bp_by_chrom[chrom]
            mut_chrom = mut_by_chrom[chrom]
            
            n_bp = len(bp_chrom)
            n_mut = len(mut_chrom)
            
            print(f"    Chr {chrom:>2}: {n_bp:,} BP × {n_mut:,} mut...", end=' ', flush=True)
            
            start_time = time.time()
            chrom_results = analyze_chromosome_optimized(chrom, bp_chrom, mut_chrom, window_sizes)
            all_results.extend(chrom_results)
            elapsed = time.time() - start_time
            
            if chrom_results:
                rad_found = sum(r['Rad_Count'] for r in chrom_results if r['Window'] == min(window_sizes))
                print(f"→ {rad_found:,} T mutations ({elapsed:.1f}s)")
            else:
                print(f"→ 0 ({elapsed:.1f}s)")
    
    results_df = pd.DataFrame(all_results)
    print(f"\n  Total results: {len(results_df):,}")
    
    return results_df


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_background_density(mutations_df):
    """Calculate genome-wide background mutation density."""
    print("\nCalculating background mutation density...")
    
    total_mutations = len(mutations_df)
    
    rad_total = mutations_df['Pattern'].isin(RADIATION_PATTERNS).sum()
    ctrl_total = mutations_df['Pattern'].isin(CONTROL_PATTERNS).sum()
    
    rad_density = rad_total / GENOME_SIZE
    ctrl_density = ctrl_total / GENOME_SIZE
    
    print(f"  Radiation mutations: {rad_total:,} ({rad_density:.2e}/bp)")
    print(f"  Control mutations:   {ctrl_total:,} ({ctrl_density:.2e}/bp)")
    
    return {
        'radiation': rad_density,
        'control': ctrl_density,
        'rad_total': rad_total,
        'ctrl_total': ctrl_total
    }


def analyze_enrichment(results_df, background, n_breakpoints, output_dir):
    """Analyze mutation enrichment near SV breakpoints."""
    print("\n" + "=" * 70)
    print("MUTATION ENRICHMENT NEAR SV BREAKPOINTS")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    print(f"\n  Unique breakpoints: {n_breakpoints:,}")
    
    print(f"\n  {'Window':<10} {'Class':<12} {'Observed':<15} {'Expected':<15} {'Enrichment':<12} {'p-value'}")
    print("  " + "-" * 80)
    
    enrichment_results = []
    
    for window in sorted(results_df['Window'].unique()):
        window_data = results_df[results_df['Window'] == window]
        
        # Window area = 2 * window * n_breakpoints
        window_area = 2 * window * n_breakpoints
        
        for mut_class, col, density_key in [('Radiation', 'Rad_Count', 'radiation'), 
                                             ('Control', 'Ctrl_Count', 'control')]:
            observed = window_data[col].sum()
            expected = background[density_key] * window_area
            
            enrichment = observed / expected if expected > 0 else 0
            
            # Poisson test
            if expected > 0:
                pval = 1 - stats.poisson.cdf(int(observed) - 1, expected)
            else:
                pval = 1.0
            
            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
            
            print(f"  {window:<10} {mut_class:<12} {observed:<15,} {expected:<15,.0f} {enrichment:<12.2f} {pval:.2e} {sig}")
            
            enrichment_results.append({
                'Window': window,
                'Mutation_Class': mut_class,
                'Observed': observed,
                'Expected': expected,
                'Enrichment': enrichment,
                'P_Value': pval
            })
    
    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df.to_csv(output_dir / 'enrichment_by_window.csv', index=False)
    
    return enrichment_df


def analyze_pattern_concordance(results_df, output_dir):
    """Analyze temporal pattern concordance."""
    print("\n" + "=" * 70)
    print("TEMPORAL PATTERN CONCORDANCE")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    # Use smallest window
    smallest_window = results_df['Window'].min()
    window_data = results_df[results_df['Window'] == smallest_window]
    
    total_rad = window_data['Rad_Count'].sum()
    concordant = window_data['Concordant'].sum()
    discordant = total_rad - concordant
    
    print(f"\n  Radiation mutations near radiation SV breakpoints:")
    print(f"    Total:      {total_rad:,}")
    print(f"    Concordant: {concordant:,} ({100*concordant/total_rad:.1f}%)")
    print(f"    Discordant: {discordant:,} ({100*discordant/total_rad:.1f}%)")
    
    # Expected ~33% if random
    expected_conc = total_rad / 3
    
    if total_rad > 0:
        chi2, pval = stats.chisquare([concordant, discordant], [expected_conc, total_rad - expected_conc])
        print(f"\n  Chi-square vs random (33%):")
        print(f"    χ² = {chi2:.2f}, p = {pval:.2e}")
        
        if pval < 0.05 and concordant > expected_conc:
            print(f"    → SIGNIFICANT: Higher concordance than expected")
            print(f"    → Supports: Same repair event causes both SV and mutations")
    
    # By SV pattern
    print(f"\n  Concordance by SV Pattern:")
    pattern_conc = window_data.groupby('SV_Pattern').agg({
        'Rad_Count': 'sum',
        'Concordant': 'sum'
    }).reset_index()
    
    pattern_conc['Concordance_Pct'] = 100 * pattern_conc['Concordant'] / pattern_conc['Rad_Count']
    pattern_conc = pattern_conc.sort_values('Rad_Count', ascending=False)
    
    print(f"  {'Pattern':<12} {'Rad Muts':<12} {'Concordant':<12} {'%'}")
    print("  " + "-" * 45)
    
    for _, row in pattern_conc.iterrows():
        print(f"  {row['SV_Pattern']:<12} {row['Rad_Count']:<12,.0f} {row['Concordant']:<12,.0f} {row['Concordance_Pct']:.1f}")
    
    pattern_conc.to_csv(output_dir / 'pattern_concordance.csv', index=False)
    
    return concordant, discordant


def analyze_by_sv_type(results_df, output_dir):
    """Analyze by SV type."""
    print("\n" + "=" * 70)
    print("SV TYPE-SPECIFIC ANALYSIS")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    smallest_window = results_df['Window'].min()
    window_data = results_df[results_df['Window'] == smallest_window]
    
    print(f"\n  Mutations near breakpoints by SV type (window={smallest_window}bp):")
    print(f"  {'SV Type':<10} {'Rad Muts':<15} {'Ctrl Muts':<15} {'Ratio':<10} {'Conc %'}")
    print("  " + "-" * 60)
    
    type_stats = window_data.groupby('SV_Type').agg({
        'Rad_Count': 'sum',
        'Ctrl_Count': 'sum',
        'Concordant': 'sum'
    }).reset_index()
    
    type_stats['Ratio'] = type_stats['Rad_Count'] / type_stats['Ctrl_Count'].replace(0, np.nan)
    type_stats['Concordance_Pct'] = 100 * type_stats['Concordant'] / type_stats['Rad_Count'].replace(0, np.nan)
    
    type_stats = type_stats.sort_values('Rad_Count', ascending=False)
    
    for _, row in type_stats.iterrows():
        print(f"  {row['SV_Type']:<10} {row['Rad_Count']:<15,.0f} {row['Ctrl_Count']:<15,.0f} {row['Ratio']:<10.2f} {row['Concordance_Pct']:.1f}")
    
    type_stats.to_csv(output_dir / 'sv_type_mutations.csv', index=False)
    
    return type_stats


def analyze_distance_distribution(results_df, output_dir):
    """Analyze mean distance of mutations from breakpoints."""
    print("\n" + "=" * 70)
    print("DISTANCE ANALYSIS")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    smallest_window = results_df['Window'].min()
    window_data = results_df[results_df['Window'] == smallest_window]
    
    # Filter to breakpoints with mutations
    with_muts = window_data[window_data['Rad_Count'] > 0]
    
    mean_dist = with_muts['Mean_Distance'].mean()
    median_dist = with_muts['Mean_Distance'].median()
    min_overall = with_muts['Min_Distance'].min()
    
    print(f"\n  Distance of radiation mutations from breakpoints:")
    print(f"    Mean of means:   {mean_dist:.0f} bp")
    print(f"    Median of means: {median_dist:.0f} bp")
    print(f"    Closest overall: {min_overall:.0f} bp")
    
    # Distribution
    dist_bins = [0, 100, 250, 500, 750, 1000]
    
    print(f"\n  Breakpoints by closest mutation distance:")
    for i in range(len(dist_bins) - 1):
        low, high = dist_bins[i], dist_bins[i+1]
        count = ((with_muts['Min_Distance'] >= low) & (with_muts['Min_Distance'] < high)).sum()
        pct = 100 * count / len(with_muts) if len(with_muts) > 0 else 0
        print(f"    {low}-{high}bp: {count:,} ({pct:.1f}%)")


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_plots(results_df, enrichment_df, output_dir):
    """Create visualization plots."""
    import matplotlib.pyplot as plt
    
    output_dir = Path(output_dir)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Enrichment by window
    ax = axes[0, 0]
    for mut_class, color in [('Radiation', '#e74c3c'), ('Control', '#3498db')]:
        data = enrichment_df[enrichment_df['Mutation_Class'] == mut_class]
        ax.plot(data['Window'], data['Enrichment'], 'o-', color=color, 
                label=mut_class, markersize=10, linewidth=2)
    
    ax.axhline(y=1, color='black', linestyle='--', linewidth=2, label='Expected')
    ax.set_xlabel('Window Size (bp)')
    ax.set_ylabel('Enrichment (Observed/Expected)')
    ax.set_title('Mutation Enrichment Near SV Breakpoints', fontweight='bold')
    ax.legend()
    ax.set_xscale('log')
    
    # 2. Radiation vs Control counts
    ax = axes[0, 1]
    smallest_window = results_df['Window'].min()
    window_data = results_df[results_df['Window'] == smallest_window]
    
    rad_total = window_data['Rad_Count'].sum()
    ctrl_total = window_data['Ctrl_Count'].sum()
    
    ax.bar(['Radiation\n(T patterns)', 'Control\n(C patterns)'], 
           [rad_total, ctrl_total], color=['#e74c3c', '#3498db'], alpha=0.7)
    ax.set_ylabel('Mutations near breakpoints')
    ax.set_title(f'Mutations within {smallest_window}bp of Breakpoints', fontweight='bold')
    
    for i, v in enumerate([rad_total, ctrl_total]):
        ax.text(i, v + v*0.02, f'{v:,}', ha='center', fontsize=12)
    
    # 3. Pattern concordance
    ax = axes[1, 0]
    concordant = window_data['Concordant'].sum()
    discordant = window_data['Rad_Count'].sum() - concordant
    
    ax.pie([concordant, discordant], labels=['Concordant\n(same timepoint)', 'Discordant'],
           colors=['#27ae60', '#e74c3c'], autopct='%1.1f%%',
           explode=(0.05, 0), startangle=90)
    ax.set_title('Temporal Pattern Concordance\n(Radiation mutations near Radiation SVs)', 
                 fontweight='bold')
    
    # 4. By SV type
    ax = axes[1, 1]
    type_data = window_data.groupby('SV_Type').agg({
        'Rad_Count': 'sum',
        'Ctrl_Count': 'sum'
    }).reset_index()
    type_data = type_data.sort_values('Rad_Count', ascending=False)
    
    x = range(len(type_data))
    width = 0.35
    
    ax.bar([i - width/2 for i in x], type_data['Rad_Count'], 
           width, label='Radiation', color='#e74c3c', alpha=0.7)
    ax.bar([i + width/2 for i in x], type_data['Ctrl_Count'],
           width, label='Control', color='#3498db', alpha=0.7)
    
    ax.set_xticks(x)
    ax.set_xticklabels(type_data['SV_Type'])
    ax.set_xlabel('SV Type')
    ax.set_ylabel('Mutation Count')
    ax.set_title('Mutations Near Breakpoints by SV Type', fontweight='bold')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sv_mutation_correlation.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n  Saved: {output_dir / 'sv_mutation_correlation.png'}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='SV-Mutation Positional Correlation (Optimized)')
    parser.add_argument('--sv-catalog', '-s', required=True, help='SV catalog with breakpoints')
    parser.add_argument('--mutations', '-m', required=True, help='Mutation directory or file')
    parser.add_argument('--output', '-o', default='sv_mutation_correlation', help='Output directory')
    parser.add_argument('--windows', '-w', default='1000,5000,10000,50000',
                        help='Window sizes in bp (comma-separated)')
    parser.add_argument('--n-jobs', '-j', type=int, default=-1,
                        help='Number of parallel jobs (-1 = all cores)')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    window_sizes = [int(w) for w in args.windows.split(',')]
    
    print("=" * 70)
    print("SV-MUTATION POSITIONAL CORRELATION (OPTIMIZED)")
    print("=" * 70)
    print("\nHypothesis: Radiation → SVs → Error-prone repair → Nearby mutations")
    print(f"Window sizes: {window_sizes}")
    print(f"Parallel jobs: {args.n_jobs}")
    
    # Load SV data (small, use pandas)
    sv_df = load_sv_catalog(args.sv_catalog)
    
    # Load mutations (large, use dask with filtering)
    mut_df = load_mutations_dask(args.mutations)
    
    # Extract breakpoints
    breakpoints_df = extract_breakpoints(sv_df)
    n_breakpoints = len(breakpoints_df)
    
    # Calculate background density
    background = calculate_background_density(mut_df)
    
    # Count mutations near breakpoints (PARALLEL)
    results_df = count_mutations_parallel(breakpoints_df, mut_df, window_sizes, n_jobs=args.n_jobs)
    
    # Save intermediate results
    results_df.to_csv(output_dir / 'breakpoint_mutation_counts.csv', index=False)
    
    # Analyze enrichment
    enrichment_df = analyze_enrichment(results_df, background, n_breakpoints, output_dir)
    
    # Analyze concordance
    analyze_pattern_concordance(results_df, output_dir)
    
    # Analyze by SV type
    analyze_by_sv_type(results_df, output_dir)
    
    # Analyze distance
    analyze_distance_distribution(results_df, output_dir)
    
    # Create plots
    if args.plot:
        print("\nGenerating plots...")
        create_plots(results_df, enrichment_df, output_dir)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    smallest_window = min(window_sizes)
    sw_rad = enrichment_df[(enrichment_df['Window'] == smallest_window) & 
                           (enrichment_df['Mutation_Class'] == 'Radiation')]
    sw_ctrl = enrichment_df[(enrichment_df['Window'] == smallest_window) & 
                            (enrichment_df['Mutation_Class'] == 'Control')]
    
    if len(sw_rad) > 0:
        rad_enrich = sw_rad['Enrichment'].values[0]
        print(f"\n  Radiation mutation enrichment ({smallest_window}bp): {rad_enrich:.2f}x")
    
    if len(sw_ctrl) > 0:
        ctrl_enrich = sw_ctrl['Enrichment'].values[0]
        print(f"  Control mutation enrichment ({smallest_window}bp):   {ctrl_enrich:.2f}x")
    
    if len(sw_rad) > 0 and rad_enrich > 1.5:
        print(f"\n  ✓ SUPPORTS HYPOTHESIS:")
        print(f"    Radiation mutations cluster near SV breakpoints")
        print(f"    → Error-prone repair during SV formation causes mutations")
    elif len(sw_rad) > 0 and rad_enrich <= 1:
        print(f"\n  ✗ DOES NOT SUPPORT HYPOTHESIS:")
        print(f"    No enrichment of mutations near breakpoints")
    
    print(f"\n  Output directory: {output_dir}")


if __name__ == "__main__":
    main()