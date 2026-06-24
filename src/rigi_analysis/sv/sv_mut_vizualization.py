#!/usr/bin/env python3
"""Figure 4: Inversion-Specific DBS Enrichment - Four Lines of Evidence.
====================================================================
CORRECTED VERSION - Uses actual PASS-filtered results

Unified story showing INV-DBS coupling through:
A. Distance-dependent coupling (enrichment decay)
B. Type-specific coupling (SV type enrichment)
C. Size-dependent coupling (mega-inversion threshold)
D. Gene-level enrichment (mega-inversion effect)

Usage:
    python create_figure4_corrected.py \
        --correlation-dir sv_correlation_PASS \
        --size-analysis inv_size_PASS \
        --output figure4_inv_dbs_unified.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Professional color scheme
COLORS = {
    'Radiation': '#D4896A',
    'Control': '#7FA6C9',
    'Expected': '#2C2C2C',
    'INV': '#2E5A88',
    'TRA': '#5B8DBE',
    'DEL': '#8C8C8C',
    'DUP': '#404040',
    'INS': '#D4D4D4',
    'Mega': '#E74C3C',
    'Other': '#95A5A6'
}

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = [
    'Arial', 'Liberation Sans', 'DejaVu Sans', 'sans-serif']
plt.rcParams['font.size'] = 11


def load_enrichment_data(corr_dir):
    """Load DBS enrichment by window from overall analysis."""
    # Try multiple possible paths - prioritize DBS_overall
    possible_files = [
        Path(corr_dir) / 'DBS_overall' / 'enrichment_by_window.csv',
        Path(corr_dir) / 'enrichment_by_window.csv',
        Path(corr_dir) / 'DBS' / 'enrichment_by_window.csv',
    ]

    for enrichment_file in possible_files:
        if enrichment_file.exists():
            print(f"  Loading overall enrichment: {enrichment_file}")
            return pd.read_csv(enrichment_file)

    print("  WARNING: enrichment_by_window.csv not found")
    return None


def load_sv_type_data(sv_type_dir):
    """Load SV type-specific DBS data from sv_type_specific_decay.py output."""
    possible_files = [
        Path(sv_type_dir) / 'DBS' / 'sv_type_enrichment_decay.csv',
        Path(sv_type_dir) / 'sv_type_enrichment_decay.csv',
        Path(sv_type_dir) / 'sv_type_enrichment_decay_DBS.csv',
    ]

    for file_path in possible_files:
        if file_path.exists():
            print(f"  Loading SV-type enrichment: {file_path}")
            df = pd.read_csv(file_path)

            # Filter to 10bp window if multiple windows present
            if 'Window' in df.columns:
                df = df[df['Window'] == 10].copy()

            # Use Enrichment_Ratio as the key metric
            if 'Enrichment_Ratio' in df.columns:
                df['Ratio'] = df['Enrichment_Ratio']
            elif 'Rad_Ctrl_Ratio' in df.columns:
                df['Ratio'] = df['Rad_Ctrl_Ratio']
            elif 'Rad_Count' in df.columns and 'Ctrl_Count' in df.columns:
                df['Ratio'] = df['Rad_Count'] / df['Ctrl_Count'].replace(0, np.nan)

            return df

    print("  WARNING: sv_type_enrichment_decay.csv not found")
    return None


def load_size_analysis(size_dir):
    """Load inversion size analysis."""
    possible_files = [
        Path(size_dir) / 'size_class_results.csv',
        Path(size_dir) / 'size_class_summary.csv',
    ]

    for size_file in possible_files:
        if size_file.exists():
            print(f"  Loading size analysis: {size_file}")
            return pd.read_csv(size_file)

    print("  WARNING: size_class_results.csv not found")
    return None


def create_figure4(corr_dir, sv_type_dir, size_dir, output):
    """Create clean 4-panel Figure 4."""
    print("\n" + "="*80)
    print("CREATING FIGURE 4: INV-DBS ENRICHMENT (4 PANELS)")
    print("="*80)

    # Load data
    print("\nLoading data...")
    enrichment_df = load_enrichment_data(corr_dir)
    sv_type_df = load_sv_type_data(sv_type_dir)
    size_df = load_size_analysis(size_dir)

    # Create 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.subplots_adjust(hspace=0.35, wspace=0.30)

    # =========================================================================
    # PANEL A: Distance-Dependent Coupling (Overall enrichment decay)
    # =========================================================================
    ax1 = axes[0, 0]

    if enrichment_df is not None and len(enrichment_df) > 0:
        rad_data = enrichment_df[enrichment_df['Mutation_Class'] == 'Radiation'].sort_values('Window')
        ctrl_data = enrichment_df[enrichment_df['Mutation_Class'] == 'Control'].sort_values('Window')

        if len(rad_data) > 0:
            ax1.plot(rad_data['Window'], rad_data['Enrichment'],
                    'o-', color=COLORS['Radiation'], label='Radiation',
                    markersize=10, linewidth=2.5, markeredgecolor='black', markeredgewidth=0.5)

        if len(ctrl_data) > 0:
            ax1.plot(ctrl_data['Window'], ctrl_data['Enrichment'],
                    'o-', color=COLORS['Control'], label='Control',
                    markersize=10, linewidth=2.5, markeredgecolor='black', markeredgewidth=0.5)

        ax1.axhline(y=1.0, color=COLORS['Expected'], linestyle='--',
                   linewidth=2, label='Expected', alpha=0.7)

        # Annotate first point (no hardcoded value)
        if len(rad_data) > 0:
            first_enrich = rad_data.iloc[0]['Enrichment']
            ax1.text(rad_data.iloc[0]['Window'], first_enrich + 0.05,
                    f'{first_enrich:.2f}×', ha='center', va='bottom',
                    fontsize=10, fontweight='bold', color=COLORS['Radiation'])

        ax1.legend(frameon=True, fontsize=10, loc='upper right')
        ax1.set_xscale('log')

        if len(rad_data) > 0:
            ax1.set_ylim(0.4, max(rad_data['Enrichment'].max() * 1.15, 1.8))
    else:
        ax1.text(0.5, 0.5, 'Data not available', ha='center', va='center',
                fontsize=12, transform=ax1.transAxes)

    ax1.set_xlabel('Window Size (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Enrichment (Observed/Expected)', fontsize=12, fontweight='bold')
    ax1.set_title('A. Distance-Dependent Coupling',
                 fontsize=13, fontweight='bold', pad=15)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # =========================================================================
    # PANEL B: Type-Specific Coupling (SV type enrichment at 10bp)
    # =========================================================================
    ax2 = axes[0, 1]

    if sv_type_df is not None and 'Ratio' in sv_type_df.columns:
        main_types = ['DEL', 'DUP', 'TRA', 'INV']
        sv_plot = sv_type_df[sv_type_df['SV_Type'].isin(main_types)].copy()
        sv_plot = sv_plot.sort_values('Ratio')

        if len(sv_plot) > 0:
            colors_sv = [COLORS.get(sv, '#888') for sv in sv_plot['SV_Type']]

            y_pos = np.arange(len(sv_plot))
            bars = ax2.barh(y_pos, sv_plot['Ratio'],
                           color=colors_sv, edgecolor='black', linewidth=1.5, alpha=0.85)

            # Highlight INV bar (highest enrichment)
            if 'INV' in sv_plot['SV_Type'].values:
                inv_idx = list(sv_plot['SV_Type']).index('INV')
                bars[inv_idx].set_edgecolor('red')
                bars[inv_idx].set_linewidth(3)
                bars[inv_idx].set_alpha(1.0)

            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(sv_plot['SV_Type'], fontsize=12, fontweight='bold')
            ax2.axvline(x=1.0, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)

            # Add value labels (no hardcoded values)
            for i, (bar, val) in enumerate(zip(bars, sv_plot['Ratio'])):
                if pd.notna(val):
                    ax2.text(val + 0.15, i, f'{val:.2f}×',
                            va='center', ha='left', fontsize=10, fontweight='bold')

            # Set x-axis limit based on data
            max_ratio = sv_plot['Ratio'].max()
            ax2.set_xlim(0, max_ratio * 1.3)
    else:
        ax2.text(0.5, 0.5, 'Data not available', ha='center', va='center',
                fontsize=12, transform=ax2.transAxes)

    ax2.set_xlabel('DBS Enrichment (Rad/Ctrl Ratio)', fontsize=12, fontweight='bold')
    ax2.set_title('B. Type-Specific Coupling (10bp)',
                 fontsize=13, fontweight='bold', pad=15)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # =========================================================================
    # PANEL C: Size-Dependent Coupling
    # =========================================================================
    ax3 = axes[1, 0]

    if size_df is not None and len(size_df) > 0:
        size_order = ['Tiny', 'Small', 'Medium', 'Large', 'Very Large', 'Huge', 'Mega']
        size_plot = size_df[size_df['Size_Class'].isin(size_order)].copy()

        if len(size_plot) > 0:
            size_plot['Size_Class'] = pd.Categorical(size_plot['Size_Class'],
                                                      categories=size_order, ordered=True)
            size_plot = size_plot.sort_values('Size_Class')

            # Color mega differently
            colors_size = [COLORS['Mega'] if sc == 'Mega' else COLORS['Other']
                          for sc in size_plot['Size_Class']]

            x_pos = np.arange(len(size_plot))
            bars = ax3.bar(x_pos, size_plot['Percent_With_DBS'],
                          color=colors_size, edgecolor='black', linewidth=1.5, alpha=0.85)

            # Highlight mega bar
            if 'Mega' in size_plot['Size_Class'].values:
                mega_idx = list(size_plot['Size_Class']).index('Mega')
                bars[mega_idx].set_edgecolor('darkred')
                bars[mega_idx].set_linewidth(3)

            ax3.set_xticks(x_pos)

            size_ranges = {
                'Tiny': '<1kb', 'Small': '1-10kb', 'Medium': '10-100kb',
                'Large': '0.1-1Mb', 'Very Large': '1-10Mb', 'Huge': '10-50Mb',
                'Mega': '≥50Mb'
            }

            size_labels = [size_ranges.get(sc, str(sc)) for sc in size_plot['Size_Class']]
            ax3.set_xticklabels(size_labels, fontsize=9, rotation=45, ha='right')

            # Add value labels
            for i, (bar, val) in enumerate(zip(bars, size_plot['Percent_With_DBS'])):
                if val > 0:
                    ax3.text(i, val + 0.5, f'{val:.1f}%',
                            ha='center', va='bottom', fontsize=9, fontweight='bold')
    else:
        ax3.text(0.5, 0.5, 'Data not available', ha='center', va='center',
                fontsize=12, transform=ax3.transAxes)

    ax3.set_ylabel('% Genes with DBS', fontsize=12, fontweight='bold')
    ax3.set_title('C. Size-Dependent Coupling',
                 fontsize=13, fontweight='bold', pad=15)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # =========================================================================
    # PANEL D: Gene-Level Enrichment (Mega vs Smaller)
    # =========================================================================
    ax4 = axes[1, 1]

    if size_df is not None and len(size_df) > 0 and 'Size_Class' in size_df.columns:
        mega_data = size_df[size_df['Size_Class'] == 'Mega']
        non_mega = size_df[size_df['Size_Class'] != 'Mega']

        if len(mega_data) > 0 and len(non_mega) > 0:
            pct_mega = mega_data['Percent_With_DBS'].values[0]

            # Calculate non-mega percentage properly
            if 'N_Genes' in size_df.columns and 'N_With_DBS' in size_df.columns:
                total_genes = non_mega['N_Genes'].sum()
                total_with_dbs = non_mega['N_With_DBS'].sum()
                pct_non_mega = (total_with_dbs / total_genes * 100) if total_genes > 0 else 0
                n_mega = int(mega_data['N_Genes'].values[0])
                n_non_mega = int(total_genes)
            else:
                # Fallback: weighted average
                pct_non_mega = non_mega['Percent_With_DBS'].mean()
                n_mega = 'N/A'
                n_non_mega = 'N/A'

            enrichment = pct_mega / pct_non_mega if pct_non_mega > 0 else 0

            categories = ['Mega-Inversions\n(≥50 Mb)', 'Smaller Inversions\n(<50 Mb)']
            values = [pct_mega, pct_non_mega]
            colors = [COLORS['Mega'], COLORS['Other']]

            bars = ax4.bar(categories, values, color=colors,
                          edgecolor='black', linewidth=2, alpha=0.85, width=0.5)

            # Add value labels with sample sizes
            for i, (bar, val) in enumerate(zip(bars, values)):
                n_val = n_mega if i == 0 else n_non_mega
                label = f'{val:.1f}%\n(n={n_val:,})' if isinstance(n_val, int) else f'{val:.1f}%'
                ax4.text(i, val + 0.5, label,
                        ha='center', va='bottom', fontsize=11, fontweight='bold')

            # Add enrichment annotation (calculate p-value text from data if available)
            if 'N_Genes' in size_df.columns:
                # Try to get p-value from mega_test_results.csv
                mega_test_file = Path(size_dir) / 'mega_test_results.csv'
                if mega_test_file.exists():
                    mega_test = pd.read_csv(mega_test_file)
                    if 'P_Value' in mega_test.columns:
                        p_val = mega_test['P_Value'].values[0]
                        p_text = f'p = {p_val:.2e}' if p_val > 1e-300 else 'p < 1e-300'
                    else:
                        p_text = 'p < 0.001'
                else:
                    p_text = 'p < 0.001'
            else:
                p_text = ''

            stats_text = f'{enrichment:.2f}× Enrichment'
            if p_text:
                stats_text += f'\n{p_text}'

            ax4.text(0.5, max(values) * 0.75, stats_text,
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                             edgecolor='orange', linewidth=2, alpha=0.9))

            ax4.set_ylim(0, max(values) * 1.35)
    else:
        ax4.text(0.5, 0.5, 'Data not available', ha='center', va='center',
                fontsize=12, transform=ax4.transAxes)

    ax4.set_ylabel('% Genes with DBS', fontsize=12, fontweight='bold')
    ax4.set_title('D. Gene-Level Enrichment',
                 fontsize=13, fontweight='bold', pad=15)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax4.set_axisbelow(True)

    # Save
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"\n✓ Saved: {output}")
    print(f"✓ Saved: {output.replace('.png', '.pdf')}")
    plt.close()

    # Print summary
    print("\n" + "="*80)
    print("FIGURE 4 SUMMARY")
    print("="*80)

    if sv_type_df is not None and 'Ratio' in sv_type_df.columns:
        inv_ratio = sv_type_df[sv_type_df['SV_Type'] == 'INV']['Ratio'].values
        if len(inv_ratio) > 0:
            print(f"\nPanel B - INV-DBS enrichment: {inv_ratio[0]:.2f}×")

    if size_df is not None and 'Size_Class' in size_df.columns:
        mega_data = size_df[size_df['Size_Class'] == 'Mega']
        if len(mega_data) > 0:
            print(f"Panel C/D - Mega-inversion DBS: {mega_data['Percent_With_DBS'].values[0]:.1f}%")

    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Figure 4: INV-DBS Enrichment (corrected version)'
    )
    parser.add_argument('--correlation-dir', required=True,
                       help='Directory with overall correlation results (for Panel A)')
    parser.add_argument('--sv-type-dir', required=True,
                       help='Directory with SV-type specific decay results (for Panel B - 7.13x)')
    parser.add_argument('--size-analysis', required=True,
                       help='Directory with inversion size analysis (for Panels C & D)')
    parser.add_argument('--output', default='figure4_inv_dbs.png',
                       help='Output filename')

    args = parser.parse_args()

    create_figure4(args.correlation_dir, args.sv_type_dir, args.size_analysis, args.output)


if __name__ == "__main__":
    main()
