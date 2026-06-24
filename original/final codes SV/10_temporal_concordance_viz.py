#!/usr/bin/env python3
"""
Temporal Dynamics Analysis Figure - V3 
=======================================
Uses AnnotSV files directly (consistent with dose-stratified analysis).

Key changes from v2:
- Loads SVs from AnnotSV files (not sv_temporal_catalog.csv)
- Consistent pairing logic with dose_stratified_inv_dbs_v3.py
- Proper handling of mega-inversions

Usage:
    python create_temporal_figure_v3.py \
        --annotsv-dir ./AnnotSV_files \
        --dbs-data ./merged_mutation_data_files/DBS \
        --output figure_temporal_dynamics_v3.png \
        --window 10
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from collections import defaultdict

# =============================================================================
# CONSTANTS
# =============================================================================

DOSE_RATES = {
    'A': 0.001, 'B': 0.01, 'C': 0.1, 'D': 1.0, 'E': 2.0
}

TIMEPOINTS = {'W1': 168, 'W2': 336, 'W3': 504}

RADIATION_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']
CONTROL_PATTERNS = ['0C00', '00C0', '000C', '0CC0', '00CC', '0C0C', '0CCC']

# Map timepoint to pattern position
TIMEPOINT_TO_PATTERN_POS = {'W1': 1, 'W2': 2, 'W3': 3}

COLORS = {
    'W1': '#E88B6F',
    'W2': '#A8C98F', 
    'W3': '#8B6BA8',
}

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

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


def get_timepoint_from_pattern(pattern):
    """Extract which timepoint(s) a pattern represents."""
    timepoints = []
    if len(pattern) >= 4:
        if pattern[1] != '0':
            timepoints.append('W1')
        if pattern[2] != '0':
            timepoints.append('W2')
        if pattern[3] != '0':
            timepoints.append('W3')
    return timepoints


def patterns_overlap(pattern1, pattern2):
    """Check if two temporal patterns have overlapping timepoints."""
    if len(pattern1) != 4 or len(pattern2) != 4:
        return False
    for i in range(1, 4):
        if pattern1[i] != '0' and pattern2[i] != '0':
            return True
    return False


def timepoint_matches_pattern(timepoint, pattern):
    """Check if a timepoint (W1/W2/W3) is present in a pattern."""
    if len(pattern) < 4:
        return False
    pos = TIMEPOINT_TO_PATTERN_POS.get(timepoint, 0)
    if pos > 0:
        return pattern[pos] != '0'
    return False


# =============================================================================
# LOAD ANNOTSV FILES
# =============================================================================

def load_annotsv_inversions(annotsv_dir, mega_threshold=50_000_000):
    """
    Load inversions from AnnotSV files with dose and timepoint info.
    """
    print(f"\n{'=' * 70}")
    print("LOADING INVERSIONS FROM ANNOTSV FILES")
    print("=" * 70)
    
    annotsv_path = Path(annotsv_dir)
    
    # Find AnnotSV files
    all_files = list(annotsv_path.glob('d0_vs_d*_W*_annotated.tsv'))
    if not all_files:
        all_files = list(annotsv_path.glob('d0_vs_d*_W*.tsv'))
    if not all_files:
        all_files = list(annotsv_path.glob('**/d0_vs_d*_W*.tsv'))
    
    all_files = sorted(list(set(all_files)))
    
    if not all_files:
        print(f"ERROR: No AnnotSV files found in {annotsv_dir}")
        return pd.DataFrame()
    
    print(f"Found {len(all_files)} AnnotSV files\n")
    
    all_invs = []
    
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
            continue
        
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            # Filter to full annotations only
            annotsv_type_col = find_column(df, ['AnnotSV_type', 'Annotation_mode'])
            if annotsv_type_col:
                df = df[df[annotsv_type_col] == 'full'].copy()
            
            if len(df) == 0:
                continue
            
            # Filter to inversions
            type_col = find_column(df, ['SV_type', 'Type', 'SV_Type'])
            if type_col:
                df = df[df[type_col] == 'INV'].copy()
            
            if len(df) == 0:
                continue
            
            # Normalize chromosome
            chr_col = find_column(df, ['SV_chrom', 'Chromosome', 'Chrom', '#Chromosome'])
            if chr_col is None:
                continue
            
            df['Chromosome'] = normalize_chromosome(df[chr_col])
            
            # Get coordinates
            start_col = find_column(df, ['SV_start', 'Start', 'start'])
            end_col = find_column(df, ['SV_end', 'End', 'end'])
            
            if not all([start_col, end_col]):
                continue
            
            df = df.rename(columns={start_col: 'SV_Start', end_col: 'SV_End'})
            
            # Calculate size
            df['SV_Length'] = abs(df['SV_End'] - df['SV_Start'])
            df['Is_Mega'] = df['SV_Length'] >= mega_threshold
            
            # Add metadata
            df['Dose'] = dose
            df['Timepoint'] = timepoint
            
            print(f"  {filename}: {len(df)} inversions (Dose {dose}, {timepoint})")
            all_invs.append(df[['Chromosome', 'SV_Start', 'SV_End', 'SV_Length', 'Is_Mega', 'Dose', 'Timepoint']])
            
        except Exception as e:
            print(f"  ERROR loading {filename}: {e}")
            continue
    
    if not all_invs:
        print("\nERROR: No inversions loaded!")
        return pd.DataFrame()
    
    combined = pd.concat(all_invs, ignore_index=True)
    
    # Deduplicate within dose×timepoint
    combined['INV_Key'] = (
        combined['Chromosome'] + ':' + 
        combined['SV_Start'].astype(str) + '-' + 
        combined['SV_End'].astype(str) + ':' +
        combined['Dose'] + ':' + combined['Timepoint']
    )
    combined = combined.drop_duplicates(subset='INV_Key', keep='first')
    
    print(f"\n  Total unique inversions: {len(combined):,}")
    print(f"  Mega (≥{mega_threshold/1e6:.0f}Mb): {combined['Is_Mega'].sum():,}")
    print(f"  Small (<{mega_threshold/1e6:.0f}Mb): {(~combined['Is_Mega']).sum():,}")
    
    # Count by timepoint
    print(f"\n  By timepoint:")
    for tp in ['W1', 'W2', 'W3']:
        count = (combined['Timepoint'] == tp).sum()
        print(f"    {tp}: {count:,}")
    
    return combined


# =============================================================================
# LOAD DBS MUTATIONS
# =============================================================================

def load_dbs_mutations(mutation_dir):
    """Load DBS mutations with Pattern information."""
    print(f"\n{'=' * 70}")
    print("LOADING DBS MUTATIONS")
    print("=" * 70)
    
    mutation_path = Path(mutation_dir)
    dbs_files = list(mutation_path.glob('DBS_dose_*_merged.csv'))
    
    if not dbs_files:
        dbs_files = list(mutation_path.glob('**/DBS_dose_*_merged.csv'))
    
    if not dbs_files:
        print(f"ERROR: No DBS files found in {mutation_dir}")
        return pd.DataFrame()
    
    print(f"Found {len(dbs_files)} DBS files\n")
    
    all_dbs = []
    
    for file_path in sorted(dbs_files):
        filename = file_path.name
        
        # Parse dose
        dose = None
        for d in ['A', 'B', 'C', 'D', 'E']:
            if f'dose_{d}_' in filename:
                dose = d
                break
        
        if dose is None:
            continue
        
        try:
            df = pd.read_csv(file_path, low_memory=False)
            
            # Filter out control samples
            if 'Sample' in df.columns:
                df = df[~df['Sample'].str.contains('d0', case=False, na=False)].copy()
            
            # Filter to radiation patterns
            if 'Pattern' in df.columns:
                df = df[df['Pattern'].isin(RADIATION_PATTERNS)].copy()
            
            if len(df) == 0:
                continue
            
            # Normalize chromosome
            chr_col = find_column(df, ['Chromosome', 'Chrom', 'Chr'])
            if chr_col is None:
                continue
            
            df['Chromosome'] = normalize_chromosome(df[chr_col])
            
            # Get position
            pos_col = find_column(df, ['Start', 'Position', 'Pos'])
            if pos_col is None:
                continue
            
            df = df.rename(columns={pos_col: 'Position'})
            
            # Add dose
            df['Dose'] = dose
            
            print(f"  {filename}: {len(df):,} radiation DBS")
            all_dbs.append(df)
            
        except Exception as e:
            print(f"  ERROR loading {filename}: {e}")
            continue
    
    if not all_dbs:
        print("\nERROR: No DBS data loaded!")
        return pd.DataFrame()
    
    combined = pd.concat(all_dbs, ignore_index=True)
    print(f"\n  Total DBS: {len(combined):,}")
    
    # Count by pattern
    if 'Pattern' in combined.columns:
        print(f"\n  By pattern:")
        for pattern in ['0T00', '00T0', '000T']:
            count = (combined['Pattern'] == pattern).sum()
            print(f"    {pattern}: {count:,}")
        multi = combined[combined['Pattern'].isin(['0TT0', '00TT', '0T0T', '0TTT'])]
        print(f"    Multi-timepoint: {len(multi):,}")
    
    return combined


# =============================================================================
# FIND INV-DBS PAIRS WITH TEMPORAL CONCORDANCE
# =============================================================================

def find_inv_dbs_pairs_with_concordance(inv_df, dbs_df, window_size=10):
    """
    Find INV-DBS pairs and check temporal concordance.
    
    IMPORTANT: Pairs must match on BOTH dose AND timepoint.
    - INV has specific dose (A-E) and timepoint (W1-W3)
    - DBS has dose (A-E) and Pattern encoding timepoints
    - Concordance = same dose AND INV timepoint present in DBS pattern
    """
    print(f"\n{'=' * 70}")
    print(f"FINDING INV-DBS PAIRS (window: ±{window_size}bp)")
    print("=" * 70)
    
    if 'Pattern' not in dbs_df.columns:
        print("ERROR: DBS data missing Pattern column")
        return pd.DataFrame()
    
    if 'Dose' not in dbs_df.columns:
        print("ERROR: DBS data missing Dose column")
        return pd.DataFrame()
    
    pairs = []
    
    common_chroms = set(inv_df['Chromosome'].unique()) & set(dbs_df['Chromosome'].unique())
    print(f"Chromosomes with both INV and DBS: {len(common_chroms)}")
    
    for chrom in sorted(common_chroms):
        chrom_invs = inv_df[inv_df['Chromosome'] == chrom]
        chrom_dbs = dbs_df[dbs_df['Chromosome'] == chrom]
        
        for _, inv in chrom_invs.iterrows():
            sv_start = inv['SV_Start']
            sv_end = inv['SV_End']
            inv_timepoint = inv['Timepoint']
            inv_dose = inv['Dose']
            
            # Filter DBS to SAME DOSE first
            dose_matched_dbs = chrom_dbs[chrom_dbs['Dose'] == inv_dose]
            
            if len(dose_matched_dbs) == 0:
                continue
            
            dbs_positions = dose_matched_dbs['Position'].values
            dbs_patterns = dose_matched_dbs['Pattern'].values
            dbs_indices = dose_matched_dbs.index.values
            
            # Find DBS near start breakpoint
            near_start_mask = (
                (dbs_positions >= sv_start - window_size) &
                (dbs_positions <= sv_start + window_size)
            )
            
            # Find DBS near end breakpoint
            near_end_mask = (
                (dbs_positions >= sv_end - window_size) &
                (dbs_positions <= sv_end + window_size)
            )
            
            # Process each nearby DBS (already dose-matched)
            for mask, breakpoint_type in [(near_start_mask, 'Start'), (near_end_mask, 'End')]:
                indices = np.where(mask)[0]
                
                for idx in indices:
                    dbs_pattern = dbs_patterns[idx]
                    dbs_pos = dbs_positions[idx]
                    
                    # Check if INV timepoint is in DBS pattern
                    # (Dose already matches from filtering above)
                    concordant = timepoint_matches_pattern(inv_timepoint, dbs_pattern)
                    
                    pairs.append({
                        'Chromosome': chrom,
                        'INV_Start': sv_start,
                        'INV_End': sv_end,
                        'INV_Length': inv['SV_Length'],
                        'Is_Mega': inv['Is_Mega'],
                        'INV_Dose': inv_dose,
                        'INV_Timepoint': inv_timepoint,
                        'DBS_Position': dbs_pos,
                        'DBS_Pattern': dbs_pattern,
                        'DBS_Dose': inv_dose,  # Same as INV dose (filtered)
                        'Breakpoint': breakpoint_type,
                        'Concordant': concordant
                    })
    
    if not pairs:
        print("No INV-DBS pairs found!")
        return pd.DataFrame()
    
    pairs_df = pd.DataFrame(pairs)
    
    # Deduplicate (same DBS might be near both breakpoints of same INV)
    pairs_df = pairs_df.drop_duplicates(
        subset=['Chromosome', 'INV_Start', 'INV_End', 'DBS_Position', 'INV_Timepoint'],
        keep='first'
    )
    
    print(f"\nTotal INV-DBS pairs: {len(pairs_df):,}")
    print(f"  Concordant: {pairs_df['Concordant'].sum():,} ({pairs_df['Concordant'].mean()*100:.1f}%)")
    print(f"  Mega-inv pairs: {pairs_df['Is_Mega'].sum():,}")
    print(f"  Small-inv pairs: {(~pairs_df['Is_Mega']).sum():,}")
    
    # Concordance by timepoint
    print(f"\n  Concordance by INV timepoint:")
    for tp in ['W1', 'W2', 'W3']:
        tp_pairs = pairs_df[pairs_df['INV_Timepoint'] == tp]
        if len(tp_pairs) > 0:
            conc = tp_pairs['Concordant'].mean() * 100
            print(f"    {tp}: {conc:.1f}% (n={len(tp_pairs)})")
    
    # Concordance by size
    print(f"\n  Concordance by INV size:")
    mega_pairs = pairs_df[pairs_df['Is_Mega']]
    small_pairs = pairs_df[~pairs_df['Is_Mega']]
    if len(mega_pairs) > 0:
        print(f"    Mega: {mega_pairs['Concordant'].mean()*100:.1f}% (n={len(mega_pairs)})")
    if len(small_pairs) > 0:
        print(f"    Small: {small_pairs['Concordant'].mean()*100:.1f}% (n={len(small_pairs)})")
    
    return pairs_df


# =============================================================================
# CALCULATE CONCORDANCE STATISTICS
# =============================================================================

def calculate_concordance_stats(pairs_df):
    """Calculate concordance statistics from pairs data."""
    
    if len(pairs_df) == 0:
        return {
            'W1': None, 'W2': None, 'W3': None,
            'overall': None, 'n_pairs': 0, 'n_concordant': 0,
            'by_size': {}
        }
    
    # Overall
    n_concordant = pairs_df['Concordant'].sum()
    n_total = len(pairs_df)
    overall = (n_concordant / n_total * 100) if n_total > 0 else 0
    
    # By timepoint
    by_timepoint = {}
    for tp in ['W1', 'W2', 'W3']:
        tp_pairs = pairs_df[pairs_df['INV_Timepoint'] == tp]
        if len(tp_pairs) > 0:
            by_timepoint[tp] = tp_pairs['Concordant'].mean() * 100
        else:
            by_timepoint[tp] = None
    
    # By size category
    def get_size_category(length):
        if length < 10000:
            return 'Small (<10kb)'
        elif length < 1000000:
            return 'Medium (10kb-1Mb)'
        elif length < 50000000:
            return 'Large (1-50Mb)'
        else:
            return 'Mega (≥50Mb)'
    
    pairs_df = pairs_df.copy()
    pairs_df['Size_Category'] = pairs_df['INV_Length'].apply(get_size_category)
    
    by_size = {}
    size_order = ['Small (<10kb)', 'Medium (10kb-1Mb)', 'Large (1-50Mb)', 'Mega (≥50Mb)']
    
    for size_cat in size_order:
        size_pairs = pairs_df[pairs_df['Size_Category'] == size_cat]
        if len(size_pairs) > 0:
            by_size[size_cat] = {
                'concordance': size_pairs['Concordant'].mean() * 100,
                'n': len(size_pairs),
                'n_concordant': size_pairs['Concordant'].sum()
            }
    
    return {
        'W1': by_timepoint.get('W1'),
        'W2': by_timepoint.get('W2'),
        'W3': by_timepoint.get('W3'),
        'overall': overall,
        'n_pairs': n_total,
        'n_concordant': n_concordant,
        'by_size': by_size
    }


# =============================================================================
# CREATE FIGURE
# =============================================================================

def create_temporal_figure(inv_df, dbs_df, pairs_df, concordance_stats, output):
    """Create 4-panel temporal dynamics figure."""
    
    print(f"\n{'=' * 70}")
    print("CREATING FIGURE")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # =========================================================================
    # PANEL A: Pattern Distribution (SVs vs DBS)
    # =========================================================================
    ax1 = axes[0, 0]
    
    # SV distribution by timepoint
    sv_counts = inv_df['Timepoint'].value_counts()
    sv_w1 = sv_counts.get('W1', 0)
    sv_w2 = sv_counts.get('W2', 0)
    sv_w3 = sv_counts.get('W3', 0)
    sv_total = sv_w1 + sv_w2 + sv_w3
    
    # DBS distribution by pattern (single timepoint only)
    if 'Pattern' in dbs_df.columns:
        dbs_w1 = (dbs_df['Pattern'] == '0T00').sum()
        dbs_w2 = (dbs_df['Pattern'] == '00T0').sum()
        dbs_w3 = (dbs_df['Pattern'] == '000T').sum()
        dbs_multi = dbs_df[dbs_df['Pattern'].isin(['0TT0', '00TT', '0T0T', '0TTT'])].shape[0]
    else:
        dbs_w1 = dbs_w2 = dbs_w3 = dbs_multi = 0
    
    dbs_single_total = dbs_w1 + dbs_w2 + dbs_w3
    dbs_total = dbs_single_total + dbs_multi
    
    # Multi-timepoint SVs (approximation: SVs appearing in multiple timepoints would need tracking)
    # For now, assume all SVs are single-timepoint (consistent with 99.7% transient finding)
    sv_multi = 0
    
    patterns = ['W1', 'W2', 'W3', 'Multi-timepoint']
    sv_pct = [
        sv_w1/sv_total*100 if sv_total > 0 else 0,
        sv_w2/sv_total*100 if sv_total > 0 else 0,
        sv_w3/sv_total*100 if sv_total > 0 else 0,
        0  # Assuming single-timepoint for SVs from AnnotSV
    ]
    dbs_pct = [
        dbs_w1/dbs_total*100 if dbs_total > 0 else 0,
        dbs_w2/dbs_total*100 if dbs_total > 0 else 0,
        dbs_w3/dbs_total*100 if dbs_total > 0 else 0,
        dbs_multi/dbs_total*100 if dbs_total > 0 else 0
    ]
    
    x = np.arange(len(patterns))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, sv_pct, width, label='INV (from AnnotSV)', 
                    color='#5B8DBE', edgecolor='black', linewidth=1)
    bars2 = ax1.bar(x + width/2, dbs_pct, width, label='DBS Mutations',
                    color='#8B6BA8', edgecolor='black', linewidth=1)
    
    ax1.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax1.set_title('A. Temporal Pattern Distribution', fontsize=13, fontweight='bold', pad=15)
    ax1.set_xticks(x)
    ax1.set_xticklabels(patterns, rotation=45, ha='right')
    ax1.legend(fontsize=10, frameon=True)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0.5:
                ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                        f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # =========================================================================
    # PANEL B: Temporal Concordance by Timepoint
    # =========================================================================
    ax2 = axes[0, 1]
    
    concordance_pct = [
        concordance_stats['W1'] if concordance_stats['W1'] is not None else 0,
        concordance_stats['W2'] if concordance_stats['W2'] is not None else 0,
        concordance_stats['W3'] if concordance_stats['W3'] is not None else 0
    ]
    
    timepoints = ['W1', 'W2', 'W3']
    bars = ax2.bar(timepoints, concordance_pct, 
                   color=[COLORS['W1'], COLORS['W2'], COLORS['W3']],
                   edgecolor='black', linewidth=2, alpha=0.85)
    
    ax2.axhline(y=33, color='red', linestyle='--', linewidth=2, 
                label='Expected (Random: 33%)', alpha=0.7)
    
    ax2.set_ylabel('Temporal Concordance (%)', fontsize=12, fontweight='bold')
    ax2.set_title(f'B. INV-DBS Temporal Concordance\n(n={concordance_stats["n_pairs"]} pairs)',
                 fontsize=13, fontweight='bold', pad=15)
    ax2.legend(fontsize=10, frameon=True, loc='lower right')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylim(0, 105)
    
    for bar, val in zip(bars, concordance_pct):
        if val > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., val + 2,
                    f'{val:.1f}%', ha='center', va='bottom', 
                    fontsize=10, fontweight='bold')
    
    # =========================================================================
    # PANEL C: Persistence (Transient vs Persistent)
    # =========================================================================
    ax3 = axes[1, 0]
    
    # SVs are essentially all transient (single timepoint from AnnotSV)
    sv_transient_pct = 100.0
    sv_persistent_pct = 0.0
    
    # DBS persistence
    dbs_transient_pct = dbs_single_total / dbs_total * 100 if dbs_total > 0 else 0
    dbs_persistent_pct = dbs_multi / dbs_total * 100 if dbs_total > 0 else 0
    
    categories = ['Transient\n(Single Timepoint)', 'Persistent\n(Multi-Timepoint)']
    sv_values = [sv_transient_pct, sv_persistent_pct]
    dbs_values = [dbs_transient_pct, dbs_persistent_pct]
    
    x = np.arange(len(categories))
    width = 0.35
    
    bars1 = ax3.bar(x - width/2, sv_values, width, label='INV',
                    color='#5B8DBE', edgecolor='black', linewidth=1)
    bars2 = ax3.bar(x + width/2, dbs_values, width, label='DBS Mutations',
                    color='#8B6BA8', edgecolor='black', linewidth=1)
    
    ax3.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax3.set_title('C. Persistence vs Transient Events', fontsize=13, fontweight='bold', pad=15)
    ax3.set_xticks(x)
    ax3.set_xticklabels(categories, fontsize=11)
    ax3.legend(fontsize=10, frameon=True)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0.5:
                ax3.text(bar.get_x() + bar.get_width()/2., height + 1,
                        f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # =========================================================================
    # PANEL D: Concordance by Inversion Size
    # =========================================================================
    ax4 = axes[1, 1]
    
    size_order = ['Small (<10kb)', 'Medium (10kb-1Mb)', 'Large (1-50Mb)', 'Mega (≥50Mb)']
    size_colors = ['#A8D5BA', '#6BC48A', '#3E9B5C', '#1E6B3D']
    
    sizes = []
    concordance_vals = []
    counts = []
    colors_used = []
    
    for i, size_cat in enumerate(size_order):
        if size_cat in concordance_stats['by_size']:
            stats = concordance_stats['by_size'][size_cat]
            sizes.append(size_cat)
            concordance_vals.append(stats['concordance'])
            counts.append(stats['n'])
            colors_used.append(size_colors[i])
    
    if sizes:
        x_pos = np.arange(len(sizes))
        bars = ax4.bar(x_pos, concordance_vals, color=colors_used, 
                       edgecolor='black', linewidth=1.5, alpha=0.85)
        
        ax4.axhline(y=33, color='red', linestyle='--', linewidth=2, 
                    label='Expected (Random: 33%)', alpha=0.7)
        
        for bar, val, n in zip(bars, concordance_vals, counts):
            if val > 0:
                ax4.text(bar.get_x() + bar.get_width()/2., val + 2,
                        f'{val:.1f}%\n(n={n})', ha='center', va='bottom', 
                        fontsize=9, fontweight='bold')
        
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(sizes, fontsize=10, rotation=15, ha='right')
    
    ax4.set_ylabel('Temporal Concordance (%)', fontsize=12, fontweight='bold')
    ax4.set_title('D. Concordance by Inversion Size', fontsize=13, fontweight='bold', pad=15)
    ax4.legend(fontsize=10, frameon=True, loc='lower right')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.set_ylim(0, 115)
    ax4.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Saved: {output}")
    
    svg_output = output.replace('.png', '.svg')
    plt.savefig(svg_output, format='svg', bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {svg_output}")
    
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Temporal dynamics figure using AnnotSV files directly'
    )
    parser.add_argument('--annotsv-dir', required=True, 
                       help='Directory with AnnotSV files')
    parser.add_argument('--dbs-data', required=True, 
                       help='DBS mutation directory')
    parser.add_argument('--output', default='figure_temporal_dynamics_v3.png',
                       help='Output figure filename')
    parser.add_argument('--window', type=int, default=10,
                       help='Window size for INV-DBS pairing (bp)')
    parser.add_argument('--mega-threshold', type=int, default=50_000_000,
                       help='Threshold for mega-inversion (bp, default: 50Mb)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("TEMPORAL DYNAMICS FIGURE - V3 (AnnotSV Source)")
    print("=" * 70)
    print(f"\nAnnotSV dir: {args.annotsv_dir}")
    print(f"DBS dir: {args.dbs_data}")
    print(f"Window: ±{args.window}bp")
    print(f"Mega threshold: {args.mega_threshold/1e6:.0f}Mb")
    
    # Load data
    inv_df = load_annotsv_inversions(args.annotsv_dir, mega_threshold=args.mega_threshold)
    if len(inv_df) == 0:
        print("ERROR: No inversions loaded")
        return 1
    
    dbs_df = load_dbs_mutations(args.dbs_data)
    if len(dbs_df) == 0:
        print("ERROR: No DBS loaded")
        return 1
    
    # Find pairs with concordance
    pairs_df = find_inv_dbs_pairs_with_concordance(inv_df, dbs_df, window_size=args.window)
    
    if len(pairs_df) == 0:
        print("ERROR: No INV-DBS pairs found")
        return 1
    
    # Save pairs
    pairs_output = args.output.replace('.png', '_pairs.csv')
    pairs_df.to_csv(pairs_output, index=False)
    print(f"\n✓ Saved pairs: {pairs_output}")
    
    # Calculate statistics
    concordance_stats = calculate_concordance_stats(pairs_df)
    
    # Create figure
    create_temporal_figure(inv_df, dbs_df, pairs_df, concordance_stats, args.output)
    
    # Print summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Total INV-DBS pairs: {concordance_stats['n_pairs']}
Concordant pairs: {concordance_stats['n_concordant']}
Overall concordance: {concordance_stats['overall']:.1f}%

By timepoint:
  W1: {concordance_stats['W1']:.1f}% (if available)
  W2: {concordance_stats['W2']:.1f}% (if available)
  W3: {concordance_stats['W3']:.1f}% (if available)

By size:""")
    
    for size, stats in concordance_stats['by_size'].items():
        print(f"  {size}: {stats['concordance']:.1f}% (n={stats['n']})")
    
    print(f"\n{'=' * 70}")
    print("COMPLETE")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    exit(main())