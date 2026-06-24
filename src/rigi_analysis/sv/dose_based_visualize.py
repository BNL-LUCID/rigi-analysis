#!/usr/bin/env python3
"""Dose Response Figure: 3-Panel Layout (Updated for Concordant Pairs).
====================================================================
A. Genome-wide INV-DBS distribution
B. Summary bar chart (INV-DBS pairs by dose)
C. Functional category distribution (6 condensed categories)

Usage:
    python dose_response_figure_v4.py \
        --inv-dbs-low inv_dbs_pairs_low.csv \
        --inv-dbs-high inv_dbs_pairs_high.csv \
        --categorized-genes categorized_genes_v4.csv \
        --output dose_response_figure_v4
"""

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

# =============================================================================
# CONFIGURATION
# =============================================================================

COLORS = {
    'low': '#4A90D9',      # Blue
    'high': '#D94A4A',     # Red
    'other': '#808080',    # Gray
}

# Category order for plotting (6 condensed categories)
CATEGORY_ORDER = [
    'Signal Transduction',
    'Gene Expression',
    'Cell Structure & Adhesion',
    'Cell Cycle & DNA Damage',
    'Development & Differentiation',
    'Metabolism & Other'
]

# Short labels for x-axis (single line for rotated display)
CATEGORY_SHORT = {
    'Signal Transduction': 'Signal Transduction',
    'Gene Expression': 'Gene Expression',
    'Cell Structure & Adhesion': 'Cell Structure & Adhesion',
    'Cell Cycle & DNA Damage': 'Cell Cycle & DNA Damage',
    'Development & Differentiation': 'Development & Diff.',
    'Metabolism & Other': 'Metabolism & Other'
}


# =============================================================================
# PANEL A: GENOME-WIDE VIEW
# =============================================================================

def create_genome_view(ax, inv_low_df, inv_high_df):
    """Create genome-wide view of concordant INV-DBS events.
    No mega-inversions (none have concordant DBS).
    """
    # Chromosome sizes (GRCh38, in Mb)
    chrom_sizes = {
        'chr1': 248.9, 'chr2': 242.2, 'chr3': 198.3, 'chr4': 190.2,
        'chr5': 181.5, 'chr6': 170.8, 'chr7': 159.3, 'chr8': 145.1,
        'chr9': 138.4, 'chr10': 133.8, 'chr11': 135.1, 'chr12': 133.3,
        'chr13': 114.4, 'chr14': 107.0, 'chr15': 102.0, 'chr16': 90.3,
        'chr17': 83.3, 'chr18': 80.4, 'chr19': 58.6, 'chr20': 64.4,
        'chr21': 46.7, 'chr22': 50.8, 'chrX': 156.0, 'chrY': 57.2
    }

    chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    # Calculate cumulative positions
    cumulative = {}
    pos = 0
    for chrom in chroms:
        cumulative[chrom] = pos
        pos += chrom_sizes.get(chrom, 100)
    total_length = pos

    # Plot chromosome background
    for i, chrom in enumerate(chroms):
        start = cumulative[chrom]
        end = start + chrom_sizes.get(chrom, 100)
        color = '#F0F0F0' if i % 2 == 0 else '#E0E0E0'
        ax.axvspan(start, end, alpha=0.5, color=color)

        # Label
        mid = (start + end) / 2
        label = chrom.replace('chr', '')
        ax.text(mid, 1.15, label, ha='center', va='bottom', fontsize=8, fontweight='bold')

    def get_unique_svs(df):
        """Get unique SVs with deduplication."""
        if df is None or len(df) == 0:
            return []

        tolerance = 100
        unique_svs = []
        df_sorted = df.sort_values(['Chromosome', 'SV_Start']).reset_index(drop=True)

        for _, row in df_sorted.iterrows():
            chrom = str(row['Chromosome'])
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            start = row['SV_Start']
            end = row['SV_End']
            length = row.get('SV_Length', abs(end - start))

            # Check for match with existing
            matched = False
            for usv in unique_svs:
                if (usv['Chromosome'] == chrom and
                    abs(usv['SV_Start'] - start) <= tolerance and
                    abs(usv['SV_End'] - end) <= tolerance):
                    usv['count'] += 1
                    matched = True
                    break

            if not matched:
                unique_svs.append({
                    'Chromosome': chrom,
                    'SV_Start': start,
                    'SV_End': end,
                    'SV_Length': length,
                    'count': 1
                })

        return unique_svs

    def plot_inversions(svs, y_center, track_height, color):
        """Plot inversions as vertical lines."""
        for sv in svs:
            chrom = sv['Chromosome']
            if chrom not in cumulative:
                continue

            mid_pos = cumulative[chrom] + (sv['SV_Start'] + sv['SV_End']) / 2 / 1e6
            sv_length = sv['SV_Length']

            # Line width based on size
            if sv_length > 0:
                base_width = 1.5 + (np.log10(max(sv_length, 100)) - 2) * 0.5
                base_width = max(1.5, min(4, base_width))
            else:
                base_width = 1.5

            # Scale by recurrence
            line_width = base_width * (1 + 0.3 * min(sv['count'] - 1, 2))
            line_height = track_height * 0.7

            # For large inversions (>1Mb), draw as rectangle
            if sv_length > 1_000_000:
                start = cumulative[chrom] + sv['SV_Start'] / 1e6
                end = cumulative[chrom] + sv['SV_End'] / 1e6
                rect = mpatches.Rectangle(
                    (start, y_center - line_height/2), end - start, line_height,
                    linewidth=2, edgecolor=color, facecolor=color, alpha=0.4
                )
                ax.add_patch(rect)
            else:
                ax.vlines(mid_pos, y_center - line_height/2, y_center + line_height/2,
                         color=color, linewidth=line_width, alpha=0.85)

        return len(svs)

    # Get unique SVs
    svs_low = get_unique_svs(inv_low_df)
    svs_high = get_unique_svs(inv_high_df)

    # Plot tracks
    track_height = 0.28
    n_low = plot_inversions(svs_low, 0.7, track_height, COLORS['low'])
    n_high = plot_inversions(svs_high, 0.3, track_height, COLORS['high'])

    print(f"  Panel A: Low dose - {n_low} unique SVs")
    print(f"  Panel A: High dose - {n_high} unique SVs")

    # Add track backgrounds
    ax.axhspan(0.7 - track_height/2, 0.7 + track_height/2, alpha=0.15, color=COLORS['low'])
    ax.axhspan(0.3 - track_height/2, 0.3 + track_height/2, alpha=0.15, color=COLORS['high'])

    # Styling
    ax.set_xlim(0, total_length)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0.3, 0.7])
    ax.set_yticklabels(['High Dose\n(D,E)', 'Low Dose\n(A,B,C)'], fontsize=11, fontweight='bold')
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('A. Genome-Wide Distribution of Temporally Concordant INV-DBS Pairs',
                 fontsize=14, fontweight='bold', pad=12)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['low'], label='Low dose', alpha=0.8),
        mpatches.Patch(facecolor=COLORS['high'], label='High dose', alpha=0.8),
        plt.Line2D([0], [0], color='gray', linewidth=1.5, label='Small inversion (<1Mb)'),
        mpatches.Patch(facecolor='gray', edgecolor='gray', alpha=0.4, label='Large inversion (≥1Mb)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9, frameon=True,
              fancybox=True, framealpha=0.95)


# =============================================================================
# PANEL B: SUMMARY STATISTICS
# =============================================================================

def create_summary_panel(ax, n_pairs_low, n_pairs_high, n_genes_low, n_genes_high):
    """Create summary bar chart of concordant INV-DBS pairs."""
    categories = ['Low Dose\n(A,B,C)', 'High Dose\n(D,E)']
    values = [n_pairs_low, n_pairs_high]
    colors = [COLORS['low'], COLORS['high']]

    bars = ax.bar(categories, values, color=colors, edgecolor='black', linewidth=2, width=0.5)

    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + max(values)*0.03,
                f'{val}', ha='center', va='bottom', fontsize=18, fontweight='bold')

    ax.set_ylabel('Concordant INV-DBS Pairs', fontsize=13, fontweight='bold')
    ax.set_title('B. INV-DBS Co-occurrence by Dose', fontsize=14, fontweight='bold', pad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, max(values) * 1.3)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=10)

    # Add gene counts annotation
    legend_text = (f'Protein-coding genes:\n'
                   f'  Low: {n_genes_low}\n'
                   f'  High: {n_genes_high}')
    ax.text(0.98, 0.98, legend_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', alpha=0.95))


# =============================================================================
# PANEL C: FUNCTIONAL CATEGORY BAR PLOT
# =============================================================================

def create_category_panel(ax, categorized_csv):
    """Create side-by-side bar plot of functional categories by dose.
    Uses 6 condensed categories.
    """
    df = pd.read_csv(categorized_csv)

    # Count genes per category and dose
    counts = df.groupby(['Functional_Category', 'Dose']).size().unstack(fill_value=0)

    # Ensure both Low and High columns exist
    if 'Low' not in counts.columns:
        counts['Low'] = 0
    if 'High' not in counts.columns:
        counts['High'] = 0

    # Reorder categories
    categories = [c for c in CATEGORY_ORDER if c in counts.index]
    counts = counts.reindex(categories)
    counts = counts.fillna(0)

    # Setup bar positions
    x = np.arange(len(counts))
    width = 0.35

    # Create bars
    bars_low = ax.bar(x - width/2, counts['Low'], width, label='Low Dose (A,B,C)',
                      color=COLORS['low'], edgecolor='black', linewidth=1.2)
    bars_high = ax.bar(x + width/2, counts['High'], width, label='High Dose (D,E)',
                       color=COLORS['high'], edgecolor='black', linewidth=1.2)

    # Add value labels on bars
    def add_labels(bars, values):
        for bar, val in zip(bars, values):
            if val > 0:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.3,
                       f'{int(val)}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    add_labels(bars_low, counts['Low'])
    add_labels(bars_high, counts['High'])

    # Styling
    ax.set_ylabel('Number of Genes', fontsize=13, fontweight='bold')
    ax.set_title('C. Functional Categories of Affected Genes', fontsize=14, fontweight='bold', pad=12)

    # X-axis labels - rotated to avoid overlap
    labels = [CATEGORY_SHORT.get(c, c) for c in counts.index]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9, ha='right', rotation=30)

    # Y-axis
    max_val = max(counts['Low'].max(), counts['High'].max())
    ax.set_ylim(0, max_val * 1.2)
    ax.tick_params(axis='y', labelsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    ax.legend(loc='upper right', fontsize=11, frameon=True, framealpha=0.95)

    # Add totals annotation
    total_low = int(counts['Low'].sum())
    total_high = int(counts['High'].sum())
    ax.text(0.02, 0.98, f'Total genes: Low={total_low}, High={total_high}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9))


# =============================================================================
# MAIN FIGURE
# =============================================================================

def create_figure(args):
    """Create the 3-panel figure."""
    print("=" * 70)
    print("DOSE RESPONSE FIGURE - CONCORDANT PAIRS (v4)")
    print("=" * 70)

    # Load INV-DBS pairs
    print("\nLoading concordant INV-DBS pairs...")
    inv_low_df = pd.read_csv(args.inv_dbs_low) if args.inv_dbs_low else None
    inv_high_df = pd.read_csv(args.inv_dbs_high) if args.inv_dbs_high else None

    n_pairs_low = len(inv_low_df) if inv_low_df is not None else 0
    n_pairs_high = len(inv_high_df) if inv_high_df is not None else 0

    print(f"  Low dose pairs: {n_pairs_low}")
    print(f"  High dose pairs: {n_pairs_high}")

    # Load categorized genes
    print("\nLoading categorized genes...")
    cat_df = pd.read_csv(args.categorized_genes)
    n_genes_low = len(cat_df[cat_df['Dose'] == 'Low'])
    n_genes_high = len(cat_df[cat_df['Dose'] == 'High'])
    print(f"  Low dose genes: {n_genes_low}")
    print(f"  High dose genes: {n_genes_high}")

    # Create figure
    print("\nCreating figure...")
    fig = plt.figure(figsize=(16, 11))

    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.85, 1], width_ratios=[0.9, 1.1],
                  hspace=0.35, wspace=0.25)

    # Panel A: Genome view (spans full width)
    ax_genome = fig.add_subplot(gs[0, :])
    create_genome_view(ax_genome, inv_low_df, inv_high_df)

    # Panel B: Summary stats
    ax_summary = fig.add_subplot(gs[1, 0])
    create_summary_panel(ax_summary, n_pairs_low, n_pairs_high, n_genes_low, n_genes_high)

    # Panel C: Functional categories
    ax_categories = fig.add_subplot(gs[1, 1])
    create_category_panel(ax_categories, args.categorized_genes)

    # Save outputs
    output_path = Path(args.output)

    # Save PNG
    png_path = output_path.with_suffix('.png')
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
    print(f"\n✓ Saved: {png_path}")

    # Save SVG
    svg_path = output_path.with_suffix('.svg')
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white', pad_inches=0.3)
    print(f"✓ Saved: {svg_path}")

    plt.close()

    print("\n" + "=" * 70)
    print("FIGURE COMPLETE")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description='Create dose response figure with concordant pairs')
    parser.add_argument('--inv-dbs-low', required=True, help='Low dose INV-DBS pairs CSV')
    parser.add_argument('--inv-dbs-high', required=True, help='High dose INV-DBS pairs CSV')
    parser.add_argument('--categorized-genes', required=True, help='Categorized genes CSV')
    parser.add_argument('--output', default='dose_response_figure_v4', help='Output filename (without extension)')

    args = parser.parse_args()
    create_figure(args)


if __name__ == "__main__":
    main()
