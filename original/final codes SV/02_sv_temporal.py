#!/usr/bin/env python3
"""
SV Temporal Pattern Tracker
============================
Track structural variants across timepoints (W1, W2, W3) to assign temporal
patterns similar to point mutations (0T00, 0TTT, 000T, etc.)

SVs are matched across timepoints by:
- Same chromosome
- Same SV type
- Breakpoints within tolerance (default ±500bp)

Output patterns:
- 0T00: SV present only at W1
- 00T0: SV present only at W2
- 000T: SV present only at W3
- 0TT0: SV present at W1 and W2, gone by W3
- 0TTT: SV present at W1, persisted through W3
- 00TT: SV appeared at W2, persisted to W3
- 0T0T: SV at W1, absent W2, reappeared W3 (intermittent)

Note: First position is always '0' since we're comparing to d0 (baseline)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import re
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

CHROM_ORDER = [str(i) for i in range(1, 23)] + ['X', 'Y']

# Breakpoint matching tolerance
BREAKPOINT_TOLERANCE = 500  # bp

# Valid temporal patterns (first char always 0 = not in baseline)
VALID_PATTERNS = [
    '0T00',  # W1 only
    '00T0',  # W2 only
    '000T',  # W3 only
    '0TT0',  # W1+W2
    '0T0T',  # W1+W3 (intermittent)
    '00TT',  # W2+W3
    '0TTT',  # W1+W2+W3 (persistent)
]

PATTERN_CATEGORIES = {
    'transient_early': ['0T00'],
    'transient_mid': ['00T0'],
    'transient_late': ['000T'],
    'persistent_full': ['0TTT'],
    'persistent_late': ['00TT'],
    'declining': ['0TT0'],
    'intermittent': ['0T0T'],
}


# =============================================================================
# SV LOADING
# =============================================================================

def load_sv_file(filepath):
    """Load and deduplicate AnnotSV file, preserving gene annotations."""
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    if 'Annotation_mode' in df.columns:
        df = df[df['Annotation_mode'] == 'full'].copy()
    # Find columns
    chrom_col = next((c for c in ['SV_chrom', 'Chromosome', 'Chr', '#CHROM'] if c in df.columns), None)
    start_col = next((c for c in ['SV_start', 'Start', 'POS'] if c in df.columns), None)
    end_col = next((c for c in ['SV_end', 'End', 'END'] if c in df.columns), None)
    type_col = next((c for c in ['SV_type', 'SVTYPE', 'Type'] if c in df.columns), None)
    gene_col = next((c for c in ['Gene_name', 'Gene name', 'GENE'] if c in df.columns), None)
    pli_col = next((c for c in ['GnomAD_pLI', 'pLI', 'ExAC_pLI'] if c in df.columns), None)
    
    if not all([chrom_col, start_col, end_col]):
        raise ValueError(f"Missing required columns in {filepath}")
    
    df = df.rename(columns={chrom_col: 'Chrom', start_col: 'Start', end_col: 'End'})
    if type_col:
        df = df.rename(columns={type_col: 'SV_Type'})
    else:
        df['SV_Type'] = 'Unknown'
    
    if gene_col:
        df = df.rename(columns={gene_col: 'Gene'})
    else:
        df['Gene'] = None
    
    if pli_col:
        df = df.rename(columns={pli_col: 'pLI'})
    else:
        df['pLI'] = None
    
    # Clean chromosome
    df['Chrom'] = df['Chrom'].astype(str).str.replace('chr', '')
    df = df[df['Chrom'].isin(CHROM_ORDER)].copy()
    
    # Deduplicate by AnnotSV_ID or coordinates
    if 'AnnotSV_ID' in df.columns:
        df['SV_ID'] = df['AnnotSV_ID'].str.replace(r'_\d+$', '', regex=True)
        df = df.drop_duplicates(subset=['SV_ID'], keep='first')
    else:
        df = df.drop_duplicates(subset=['Chrom', 'Start', 'End', 'SV_Type'], keep='first')
    
    # Calculate length
    df['SV_Length'] = abs(df['End'] - df['Start'])
    
    # Create unique key for matching
    df['Match_Key'] = df['Chrom'] + '_' + df['SV_Type'] + '_' + df['Start'].astype(str) + '_' + df['End'].astype(str)
    
    return df[['Chrom', 'Start', 'End', 'SV_Type', 'SV_Length', 'Match_Key', 'Gene', 'pLI']].reset_index(drop=True)


# =============================================================================
# SV MATCHING ACROSS TIMEPOINTS
# =============================================================================

def match_sv(sv_row, other_df, tolerance=BREAKPOINT_TOLERANCE):
    """
    Check if an SV matches any SV in another timepoint's data.
    
    Matching criteria:
    - Same chromosome
    - Same SV type
    - Start positions within tolerance
    - End positions within tolerance
    """
    chrom = sv_row['Chrom']
    sv_type = sv_row['SV_Type']
    start = sv_row['Start']
    end = sv_row['End']
    
    # Filter to same chrom and type
    candidates = other_df[
        (other_df['Chrom'] == chrom) & 
        (other_df['SV_Type'] == sv_type)
    ]
    
    if len(candidates) == 0:
        return False
    
    # Check breakpoint proximity
    start_match = (candidates['Start'] >= start - tolerance) & (candidates['Start'] <= start + tolerance)
    end_match = (candidates['End'] >= end - tolerance) & (candidates['End'] <= end + tolerance)
    
    return (start_match & end_match).any()


def build_sv_catalog(sv_w1, sv_w2, sv_w3, tolerance=BREAKPOINT_TOLERANCE):
    """
    Build a unified catalog of all SVs across timepoints with temporal patterns.
    
    Strategy:
    1. Start with all unique SVs from all timepoints
    2. For each SV, check presence in W1, W2, W3
    3. Assign temporal pattern
    """
    print("  Building unified SV catalog...")
    
    # Combine all SVs
    all_svs = []
    
    # Add W1 SVs
    if len(sv_w1) > 0:
        w1_copy = sv_w1.copy()
        w1_copy['Source'] = 'W1'
        all_svs.append(w1_copy)
    
    # Add W2 SVs
    if len(sv_w2) > 0:
        w2_copy = sv_w2.copy()
        w2_copy['Source'] = 'W2'
        all_svs.append(w2_copy)
    
    # Add W3 SVs
    if len(sv_w3) > 0:
        w3_copy = sv_w3.copy()
        w3_copy['Source'] = 'W3'
        all_svs.append(w3_copy)
    
    if not all_svs:
        return pd.DataFrame()
    
    combined = pd.concat(all_svs, ignore_index=True)
    
    # Deduplicate across timepoints (keep first occurrence)
    # Group by approximate coordinates
    combined['Chrom_Type'] = combined['Chrom'] + '_' + combined['SV_Type']
    combined = combined.sort_values(['Chrom_Type', 'Start', 'End'])
    
    # Build unique SV list by merging similar SVs
    unique_svs = []
    processed = set()
    
    print(f"  Processing {len(combined)} total SV entries...")
    
    for idx, row in combined.iterrows():
        # Create a signature for this SV
        sig = (row['Chrom'], row['SV_Type'], row['Start'], row['End'])
        
        # Check if we've already processed a similar SV
        skip = False
        for p_sig in processed:
            if (p_sig[0] == sig[0] and  # same chrom
                p_sig[1] == sig[1] and  # same type
                abs(p_sig[2] - sig[2]) <= tolerance and  # start within tolerance
                abs(p_sig[3] - sig[3]) <= tolerance):    # end within tolerance
                skip = True
                break
        
        if skip:
            continue
        
        processed.add(sig)
        
        # Check presence in each timepoint
        in_w1 = match_sv(row, sv_w1, tolerance) if len(sv_w1) > 0 else False
        in_w2 = match_sv(row, sv_w2, tolerance) if len(sv_w2) > 0 else False
        in_w3 = match_sv(row, sv_w3, tolerance) if len(sv_w3) > 0 else False
        
        # Build pattern (first char is always 0 = not in baseline)
        pattern = '0'
        pattern += 'T' if in_w1 else '0'
        pattern += 'T' if in_w2 else '0'
        pattern += 'T' if in_w3 else '0'
        
        unique_svs.append({
            'Chrom': row['Chrom'],
            'Start': row['Start'],
            'End': row['End'],
            'SV_Type': row['SV_Type'],
            'SV_Length': row['SV_Length'],
            'Gene': row.get('Gene', None),
            'pLI': row.get('pLI', None),
            'In_W1': in_w1,
            'In_W2': in_w2,
            'In_W3': in_w3,
            'Pattern': pattern
        })
    
    result_df = pd.DataFrame(unique_svs)
    print(f"  Identified {len(result_df)} unique SVs with temporal patterns")
    
    return result_df


# =============================================================================
# PATTERN ANALYSIS
# =============================================================================

def categorize_pattern(pattern):
    """Categorize a pattern into biological meaning."""
    for category, patterns in PATTERN_CATEGORIES.items():
        if pattern in patterns:
            return category
    return 'other'


def summarize_patterns(sv_catalog):
    """Summarize temporal pattern distribution."""
    if len(sv_catalog) == 0:
        return {}
    
    pattern_counts = sv_catalog['Pattern'].value_counts().to_dict()
    
    # Add category summaries
    sv_catalog['Category'] = sv_catalog['Pattern'].apply(categorize_pattern)
    category_counts = sv_catalog['Category'].value_counts().to_dict()
    
    # Calculate percentages
    total = len(sv_catalog)
    pattern_pct = {k: round(100 * v / total, 1) for k, v in pattern_counts.items()}
    
    return {
        'total_unique_svs': total,
        'pattern_counts': pattern_counts,
        'pattern_percentages': pattern_pct,
        'category_counts': category_counts
    }


def analyze_by_sv_type(sv_catalog):
    """Analyze patterns broken down by SV type."""
    if len(sv_catalog) == 0:
        return pd.DataFrame()
    
    results = []
    for sv_type in sv_catalog['SV_Type'].unique():
        type_df = sv_catalog[sv_catalog['SV_Type'] == sv_type]
        pattern_counts = type_df['Pattern'].value_counts()
        
        for pattern, count in pattern_counts.items():
            results.append({
                'SV_Type': sv_type,
                'Pattern': pattern,
                'Count': count,
                'Pct_of_Type': round(100 * count / len(type_df), 1)
            })
    
    return pd.DataFrame(results)


# =============================================================================
# BATCH ANALYSIS
# =============================================================================

def analyze_dose(sv_dir, dose, tolerance=BREAKPOINT_TOLERANCE):
    """Analyze temporal patterns for a single dose across W1, W2, W3."""
    sv_dir = Path(sv_dir)
    
    # Load SV files for each timepoint
    sv_data = {}
    for tp in ['W1', 'W2', 'W3']:
        # Try different filename patterns
        patterns = [
            f'd0_vs_d{dose}_{tp}_annotated.tsv',
            f'd0_vs_d{dose}_{tp}.tsv',
            f'd{dose}_{tp}_annotated.tsv'
        ]
        
        found = False
        for pattern in patterns:
            filepath = sv_dir / pattern
            if filepath.exists():
                sv_data[tp] = load_sv_file(filepath)
                print(f"    {tp}: {len(sv_data[tp])} SVs")
                found = True
                break
        
        if not found:
            print(f"    {tp}: No file found")
            sv_data[tp] = pd.DataFrame(columns=['Chrom', 'Start', 'End', 'SV_Type', 'SV_Length', 'Match_Key'])
    
    # Build temporal catalog
    sv_catalog = build_sv_catalog(
        sv_data.get('W1', pd.DataFrame()),
        sv_data.get('W2', pd.DataFrame()),
        sv_data.get('W3', pd.DataFrame()),
        tolerance
    )
    
    if len(sv_catalog) > 0:
        sv_catalog['Dose'] = f'd{dose}'
    
    return sv_catalog


def batch_analysis(sv_dir, output_dir, tolerance=BREAKPOINT_TOLERANCE):
    """Run temporal pattern analysis for all doses."""
    sv_dir = Path(sv_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("SV TEMPORAL PATTERN ANALYSIS")
    print("=" * 60)
    print(f"Breakpoint tolerance: ±{tolerance} bp")
    
    all_catalogs = []
    all_summaries = []
    
    for dose in ['A', 'B', 'C', 'D', 'E']:
        print(f"\n--- DOSE {dose} ---")
        
        sv_catalog = analyze_dose(sv_dir, dose, tolerance)
        
        if len(sv_catalog) > 0:
            all_catalogs.append(sv_catalog)
            
            # Summarize
            summary = summarize_patterns(sv_catalog)
            summary['dose'] = f'd{dose}'
            all_summaries.append(summary)
            
            print(f"\n  Pattern distribution:")
            for pattern, count in sorted(summary['pattern_counts'].items()):
                pct = summary['pattern_percentages'][pattern]
                print(f"    {pattern}: {count:,} ({pct}%)")
    
    # Combine results
    if all_catalogs:
        combined_catalog = pd.concat(all_catalogs, ignore_index=True)
        combined_catalog.to_csv(output_dir / 'sv_temporal_catalog.csv', index=False)
        
        # Create summary table
        summary_rows = []
        for dose in ['A', 'B', 'C', 'D', 'E']:
            dose_df = combined_catalog[combined_catalog['Dose'] == f'd{dose}']
            if len(dose_df) == 0:
                continue
            
            row = {'Dose': f'd{dose}', 'Total_SVs': len(dose_df)}
            for pattern in VALID_PATTERNS:
                row[pattern] = (dose_df['Pattern'] == pattern).sum()
            summary_rows.append(row)
        
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(output_dir / 'sv_temporal_summary.csv', index=False)
        
        # By SV type analysis
        type_analysis = analyze_by_sv_type(combined_catalog)
        type_analysis.to_csv(output_dir / 'sv_temporal_by_type.csv', index=False)
        
        print("\n" + "=" * 60)
        print("OVERALL SUMMARY")
        print("=" * 60)
        print(summary_df.to_string(index=False))
        
        print(f"\nResults saved to: {output_dir}")
        
        return combined_catalog, summary_df
    
    return pd.DataFrame(), pd.DataFrame()


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_temporal_patterns(summary_df, catalog_df, output_dir):
    """Create visualizations of temporal patterns."""
    import matplotlib.pyplot as plt
    
    output_dir = Path(output_dir)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Pattern distribution by dose (stacked bar)
    ax = axes[0, 0]
    if len(summary_df) > 0:
        pattern_cols = [c for c in summary_df.columns if c.startswith('0')]
        if pattern_cols:
            summary_df.set_index('Dose')[pattern_cols].plot(
                kind='bar', stacked=True, ax=ax, 
                colormap='Set2'
            )
            ax.set_xlabel('Dose')
            ax.set_ylabel('Number of SVs')
            ax.set_title('SV Temporal Patterns by Dose', fontweight='bold')
            ax.legend(title='Pattern', bbox_to_anchor=(1.02, 1), loc='upper left')
            ax.tick_params(axis='x', rotation=0)
    
    # 2. Pattern percentages heatmap
    ax = axes[0, 1]
    if len(summary_df) > 0 and len(pattern_cols) > 0:
        pct_df = summary_df.set_index('Dose')[pattern_cols].copy()
        pct_df = pct_df.div(pct_df.sum(axis=1), axis=0) * 100
        
        im = ax.imshow(pct_df.values, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(len(pct_df.columns)))
        ax.set_xticklabels(pct_df.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(pct_df.index)))
        ax.set_yticklabels(pct_df.index)
        
        for i in range(len(pct_df.index)):
            for j in range(len(pct_df.columns)):
                ax.text(j, i, f'{pct_df.values[i, j]:.1f}%', 
                       ha='center', va='center', fontsize=8)
        
        ax.set_title('Pattern Distribution (% per dose)', fontweight='bold')
        plt.colorbar(im, ax=ax, label='Percentage')
    
    # 3. SV type vs pattern
    ax = axes[1, 0]
    if len(catalog_df) > 0:
        type_pattern = catalog_df.groupby(['SV_Type', 'Pattern']).size().unstack(fill_value=0)
        type_pattern.plot(kind='bar', ax=ax, colormap='Set2')
        ax.set_xlabel('SV Type')
        ax.set_ylabel('Count')
        ax.set_title('Temporal Patterns by SV Type', fontweight='bold')
        ax.legend(title='Pattern', bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=45)
    
    # 4. Persistent vs Transient summary
    ax = axes[1, 1]
    if len(catalog_df) > 0:
        # Categorize
        def get_category(pattern):
            if pattern in ['0TTT']:
                return 'Persistent (all weeks)'
            elif pattern in ['00TT', '0TT0']:
                return 'Persistent (2 weeks)'
            elif pattern in ['0T00', '00T0', '000T']:
                return 'Transient (1 week)'
            elif pattern in ['0T0T']:
                return 'Intermittent'
            else:
                return 'Other'
        
        catalog_df['Category'] = catalog_df['Pattern'].apply(get_category)
        cat_by_dose = catalog_df.groupby(['Dose', 'Category']).size().unstack(fill_value=0)
        
        cat_by_dose.plot(kind='bar', ax=ax, colormap='Set1')
        ax.set_xlabel('Dose')
        ax.set_ylabel('Number of SVs')
        ax.set_title('SV Persistence Categories by Dose', fontweight='bold')
        ax.legend(title='Category', bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=0)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sv_temporal_patterns.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to: {output_dir}")


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Track SV temporal patterns across timepoints')
    parser.add_argument('sv_dir', help='Directory with AnnotSV files')
    parser.add_argument('--output', '-o', default='sv_temporal_results', help='Output directory')
    parser.add_argument('--tolerance', '-t', type=int, default=1000, 
                        help='Breakpoint matching tolerance in bp (default: 1000)')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    
    args = parser.parse_args()
    
    catalog_df, summary_df = batch_analysis(args.sv_dir, args.output, args.tolerance)
    
    if args.plot and len(catalog_df) > 0:
        plot_temporal_patterns(summary_df, catalog_df, args.output)