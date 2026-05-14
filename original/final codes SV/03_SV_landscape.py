#!/usr/bin/env python3
"""
Figure 3: SV Landscape - Revised with dose and temporal analysis
================================================================
Panel A: Total SV counts by type
Panel B: SV counts by dose (A→E)
Panel C: SV counts by timepoint (W1→W3)
Panel D: Heatmap - Dose × Timepoint

Usage:
    python figure3_sv_landscape_v2.py \
        --annotsv-dir ./SV_files/annoted_passed/ \
        --output figure3_sv_landscape_v2.png
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re

# Professional colors
COLORS = {
    'INV': '#2E5A88', 'TRA': '#5B8DBE', 'DEL': '#8C8C8C',
    'DUP': '#404040', 'INS': '#D4D4D4', 'BND': '#5B8DBE'
}

SV_ORDER = ['INV', 'TRA', 'DEL', 'DUP', 'INS']

# Dose labels and order
DOSE_ORDER = ['A', 'B', 'C', 'D', 'E']
DOSE_LABELS = {
    'A': '0.36',
    'B': '0.20', 
    'C': '0.40',
    'D': '1.47',
    'E': '2.62'
}

# Timepoint order
WEEK_ORDER = ['W1', 'W2', 'W3']

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11


def parse_filename(filename):
    """Extract dose and week from filename like d0_vs_dA_W1_annotated.tsv"""
    match = re.search(r'd0_vs_d([A-E])_W(\d)', filename)
    if match:
        dose = match.group(1)
        week = f'W{match.group(2)}'
        return dose, week
    return None, None


def load_all_data(annotsv_dir):
    """Load all AnnotSV files and return structured data."""
    print("\nLoading AnnotSV files...")
    
    annotsv_path = Path(annotsv_dir)
    files = list(annotsv_path.glob('*.tsv'))
    print(f"  Found {len(files)} files")
    
    all_data = []
    
    for f in files:
        dose, week = parse_filename(f.name)
        if dose is None:
            print(f"  Skipping {f.name} - couldn't parse dose/week")
            continue
        
        try:
            df = pd.read_csv(f, sep='\t', low_memory=False)
            
            # Filter to full annotation mode for unique SV counting
            if 'Annotation_mode' in df.columns:
                # For counting unique SVs, use 'full' rows (one per SV)
                df_full = df[df['Annotation_mode'] == 'full'].copy()
            else:
                df_full = df.copy()
            
            # Map BND to TRA
            if 'SV_type' in df_full.columns:
                df_full.loc[df_full['SV_type'] == 'BND', 'SV_type'] = 'TRA'
            
            # Count unique SVs by AnnotSV_ID within this file
            if 'AnnotSV_ID' in df_full.columns:
                unique_svs = df_full.drop_duplicates(subset=['AnnotSV_ID'])
            else:
                unique_svs = df_full
            
            # Add metadata
            unique_svs = unique_svs.copy()
            unique_svs['Dose'] = dose
            unique_svs['Week'] = week
            unique_svs['File'] = f.name
            
            all_data.append(unique_svs)
            
            sv_counts = unique_svs['SV_type'].value_counts()
            print(f"  {f.name}: {len(unique_svs)} unique SVs")
            
        except Exception as e:
            print(f"  Error loading {f.name}: {e}")
    
    if not all_data:
        return None
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\n  Total rows: {len(combined):,}")
    
    return combined


def calculate_statistics(df):
    """Calculate various statistics from the data."""
    stats = {}
    
    # Total counts by SV type
    stats['total_by_type'] = df.groupby('SV_type').size().reindex(SV_ORDER, fill_value=0)
    
    # Counts by dose and SV type
    stats['by_dose'] = df.groupby(['Dose', 'SV_type']).size().unstack(fill_value=0)
    stats['by_dose'] = stats['by_dose'].reindex(DOSE_ORDER)
    stats['by_dose'] = stats['by_dose'].reindex(columns=[s for s in SV_ORDER if s in stats['by_dose'].columns], fill_value=0)
    
    # Counts by week and SV type
    stats['by_week'] = df.groupby(['Week', 'SV_type']).size().unstack(fill_value=0)
    stats['by_week'] = stats['by_week'].reindex(WEEK_ORDER)
    stats['by_week'] = stats['by_week'].reindex(columns=[s for s in SV_ORDER if s in stats['by_week'].columns], fill_value=0)
    
    # Heatmap: Dose × Week (total SVs)
    stats['dose_week_total'] = df.groupby(['Dose', 'Week']).size().unstack(fill_value=0)
    stats['dose_week_total'] = stats['dose_week_total'].reindex(DOSE_ORDER)
    stats['dose_week_total'] = stats['dose_week_total'].reindex(columns=WEEK_ORDER, fill_value=0)
    
    # Heatmap: Dose × Week for each SV type
    stats['dose_week_by_type'] = {}
    for sv in SV_ORDER:
        sv_df = df[df['SV_type'] == sv]
        if len(sv_df) > 0:
            hm = sv_df.groupby(['Dose', 'Week']).size().unstack(fill_value=0)
            hm = hm.reindex(DOSE_ORDER, fill_value=0)
            hm = hm.reindex(columns=WEEK_ORDER, fill_value=0)
            stats['dose_week_by_type'][sv] = hm
    
    return stats


def create_figure(stats, output):
    """Create the final figure."""
    print("\nCreating figure...")
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # =========================================================================
    # Panel A: Total SV counts by type (horizontal bar)
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    
    total_counts = stats['total_by_type']
    sv_order_sorted = total_counts.sort_values().index.tolist()
    colors = [COLORS.get(sv, '#888888') for sv in sv_order_sorted]
    
    bars = ax1.barh(range(len(sv_order_sorted)), [total_counts[sv] for sv in sv_order_sorted],
                    color=colors, edgecolor='black', alpha=0.85, linewidth=1.5)
    
    ax1.set_yticks(range(len(sv_order_sorted)))
    ax1.set_yticklabels(sv_order_sorted, fontsize=12, fontweight='bold')
    ax1.set_xlabel('Unique SV Events', fontsize=12, fontweight='bold')
    ax1.set_title('A. Total SV Counts by Type', fontsize=14, fontweight='bold', pad=15, loc='left')
    ax1.set_xscale('log')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(axis='x', alpha=0.3)
    
    # Add count labels
    for i, sv in enumerate(sv_order_sorted):
        val = total_counts[sv]
        if val > 0:
            ax1.text(val * 1.3, i, f'{int(val):,}', ha='left', va='center', 
                    fontsize=10, fontweight='bold')
    
    # =========================================================================
    # Panel B: SV counts by dose
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])
    
    by_dose = stats['by_dose']
    x = np.arange(len(DOSE_ORDER))
    width = 0.15
    
    # Plot bars for each SV type
    sv_types_to_plot = [sv for sv in SV_ORDER if sv in by_dose.columns and by_dose[sv].sum() > 0]
    n_types = len(sv_types_to_plot)
    
    for i, sv in enumerate(sv_types_to_plot):
        offset = (i - n_types/2 + 0.5) * width
        counts = [by_dose.loc[d, sv] if d in by_dose.index else 0 for d in DOSE_ORDER]
        ax2.bar(x + offset, counts, width, label=sv, color=COLORS.get(sv, '#888888'),
               edgecolor='black', linewidth=0.5, alpha=0.85)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{d}\n({DOSE_LABELS[d]})' for d in DOSE_ORDER], fontsize=10)
    ax2.set_xlabel('Dose Level (mGy/hr)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('SV Count', fontsize=12, fontweight='bold')
    ax2.set_title('B. SV Counts by Dose', fontsize=14, fontweight='bold', pad=15, loc='left')
    ax2.legend(loc='upper left', fontsize=9, frameon=True)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.3)
    
    # =========================================================================
    # Panel C: SV counts by timepoint (grouped bars)
    # =========================================================================
    ax3 = fig.add_subplot(gs[1, 0])
    
    by_week = stats['by_week']
    x = np.arange(len(WEEK_ORDER))
    width = 0.15
    
    # Grouped bar chart
    for i, sv in enumerate(sv_types_to_plot):
        offset = (i - n_types/2 + 0.5) * width
        counts = [by_week.loc[w, sv] if w in by_week.index else 0 for w in WEEK_ORDER]
        ax3.bar(x + offset, counts, width, label=sv, color=COLORS.get(sv, '#888888'),
               edgecolor='black', linewidth=0.5, alpha=0.85)
    
    ax3.set_xticks(x)
    ax3.set_xticklabels(['Week 1', 'Week 2', 'Week 3'], fontsize=11, fontweight='bold')
    ax3.set_xlabel('Timepoint', fontsize=12, fontweight='bold')
    ax3.set_ylabel('SV Count', fontsize=12, fontweight='bold')
    ax3.set_title('C. SV Counts by Timepoint', fontsize=14, fontweight='bold', pad=15, loc='left')
    ax3.legend(loc='upper left', fontsize=9, frameon=True)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.grid(axis='y', alpha=0.3)
    
    # =========================================================================
    # Panel D: Heatmap - Dose × Timepoint
    # =========================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    
    heatmap_data = stats['dose_week_total']
    
    # Create annotation with values
    annot = heatmap_data.values.astype(int)
    
    sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='Blues',
                ax=ax4, linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'SV Count'})
    
    ax4.set_yticklabels([f'{d} ({DOSE_LABELS[d]} mGy/hr)' for d in DOSE_ORDER], 
                        rotation=0, fontsize=10)
    ax4.set_xticklabels(['Week 1', 'Week 2', 'Week 3'], fontsize=10)
    ax4.set_xlabel('Timepoint', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Dose Level', fontsize=12, fontweight='bold')
    ax4.set_title('D. SV Burden: Dose × Timepoint', fontsize=14, fontweight='bold', pad=15, loc='left')
    
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.savefig(output.replace('.png', '.svg'), bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output}")
    print(f"  Saved: {output.replace('.png', '.pdf')}")
    print(f"  Saved: {output.replace('.png', '.svg')}")


def print_summary(stats):
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    print("\nTotal SV counts by type:")
    for sv in SV_ORDER:
        count = stats['total_by_type'].get(sv, 0)
        print(f"  {sv}: {count:,}")
    print(f"  TOTAL: {stats['total_by_type'].sum():,}")
    
    print("\nSV counts by dose:")
    print(stats['by_dose'].to_string())
    print(f"\nTotals: {stats['by_dose'].sum(axis=1).to_dict()}")
    
    print("\nSV counts by timepoint:")
    print(stats['by_week'].to_string())
    print(f"\nTotals: {stats['by_week'].sum(axis=1).to_dict()}")
    
    print("\nDose × Timepoint heatmap:")
    print(stats['dose_week_total'].to_string())


def main():
    parser = argparse.ArgumentParser(
        description='Figure 3: SV Landscape with dose and temporal analysis'
    )
    parser.add_argument('--annotsv-dir', required=True,
                       help='Directory with AnnotSV files')
    parser.add_argument('--output', default='figure3_sv_landscape_v2.png',
                       help='Output figure filename')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("FIGURE 3: SV LANDSCAPE (DOSE & TEMPORAL)")
    print("=" * 70)
    
    # Load data
    df = load_all_data(args.annotsv_dir)
    
    if df is None or len(df) == 0:
        print("ERROR: No data loaded")
        return 1
    
    # Calculate statistics
    stats = calculate_statistics(df)
    
    # Print summary
    print_summary(stats)
    
    # Create figure
    create_figure(stats, args.output)
    
    # Save CSV summaries
    output_dir = Path(args.output).parent
    stats['by_dose'].to_csv(output_dir / 'sv_counts_by_dose.csv')
    stats['by_week'].to_csv(output_dir / 'sv_counts_by_week.csv')
    stats['dose_week_total'].to_csv(output_dir / 'sv_counts_dose_week.csv')
    print(f"\nSaved CSV files to {output_dir}")
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    exit(main())