#!/usr/bin/env python3
"""
Direct Inversion Size Analysis from AnnotSV + DBS
=================================================
Analyze INV-DBS enrichment by inversion size WITHOUT gene features.

Directly counts:
1. INV sizes from AnnotSV
2. DBS within 10bp of INV breakpoints
3. Enrichment by size class

Usage:
    python inversion_size_direct_analysis.py \
        --annotsv-dir path/to/annotsv \
        --dbs-dir path/to/DBS \
        --window 10 \
        --output inversion_size_direct.png
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from scipy import stats
from collections import defaultdict
import glob
from fast_dbs_search import fast_dbs_count

SIZE_CLASSES = [
    ('Tiny', 0, 1_000),
    ('Small', 1_000, 10_000),
    ('Medium', 10_000, 100_000),
    ('Large', 100_000, 1_000_000),
    ('Very Large', 1_000_000, 10_000_000),
    ('Huge', 10_000_000, 50_000_000),
    ('Mega', 50_000_000, float('inf'))
]

COLORS = {
    'With_DBS': '#E88B6F',
    'Without_DBS': '#5B8DBE'
}


def load_inv_from_annotsv(annotsv_dir):
    """Load INV events with sizes and breakpoints."""
    print("\nLoading INV data from AnnotSV...")
    
    files = list(Path(annotsv_dir).glob('*.tsv'))
    print(f"  Found {len(files)} files")
    
    all_inv = []
    
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', low_memory=False)
            
            # Keep only full annotations
            if 'Annotation_mode' in df.columns:
                df = df[df['Annotation_mode'] == 'full'].copy()
            
            # Filter to INV
            if 'SV_type' in df.columns:
                df = df[df['SV_type'] == 'INV'].copy()
            
            if len(df) == 0:
                continue
            
            # Extract and process data (vectorized approach)
            df['Chrom'] = df['SV_chrom'].astype(str).str.replace('chr', '')
            df['Start'] = df['SV_start']
            df['End'] = df['SV_end']
            df['Size'] = df['SV_length'].abs()
            df['Genes'] = df['Gene_name'].astype(str).str.split(';')
            
            # Explode genes (creates one row per gene)
            df = df.explode('Genes')
            df['Gene'] = df['Genes'].str.strip()
            
            # Filter valid genes
            df = df[df['Gene'].notna()]
            df = df[df['Gene'] != '']
            df = df[df['Gene'] != 'nan']
            
            # Select final columns
            inv_data = df[['Gene', 'Chrom', 'Start', 'End', 'Size']].copy()
            inv_data['Breakpoint_Start'] = inv_data['Start']
            inv_data['Breakpoint_End'] = inv_data['End']
            
            if len(inv_data) > 0:
                all_inv.append(inv_data)
                
        except Exception as e:
            print(f"  Warning: {f.name}: {e}")
    
    if all_inv:
        result = pd.concat(all_inv, ignore_index=True)
        
        # CRITICAL: Deduplicate - same INV appears across multiple files
        print(f"  Before deduplication: {len(result):,} gene-INV pairs")
        
        result = result.drop_duplicates(subset=['Gene', 'Chrom', 'Start', 'End'])
        
        print(f"  After deduplication: {len(result):,} unique gene-INV pairs")
        print(f"  Unique genes: {result['Gene'].nunique():,}")
        
        return result
    
    return pd.DataFrame()


def load_dbs_mutations(dbs_dir):
    """Load DBS mutations."""
    print("\nLoading DBS mutations...")
    
    dbs_files = glob.glob(f"{dbs_dir}/*.csv")
    print(f"  Found {len(dbs_files)} files")
    
    all_dbs = []
    
    for f in dbs_files:
        try:
            df = pd.read_csv(f, low_memory=False)
            
            # Keep only essential columns
            if 'Chromosome' in df.columns and 'Start' in df.columns:
                df = df[['Chromosome', 'Start']].copy()
                df['Chromosome'] = df['Chromosome'].astype(str).str.replace('chr', '')
                all_dbs.append(df)
        except:
            pass
    
    if all_dbs:
        combined = pd.concat(all_dbs, ignore_index=True)
        print(f"  Total DBS: {len(combined):,}")
        return combined
    
    return pd.DataFrame()


def find_dbs_near_inv(inv_df, dbs_df, window=10):
    """Find DBS within window of INV breakpoints."""
    print(f"\nFinding DBS within {window}bp of INV breakpoints...")
    
    # Create breakpoint list (vectorized - MUCH faster)
    # Start breakpoints
    bp_start = inv_df[['Gene', 'Chrom', 'Start', 'Size']].copy()
    bp_start.columns = ['Gene', 'Chrom', 'Pos', 'Size']
    
    # End breakpoints
    bp_end = inv_df[['Gene', 'Chrom', 'End', 'Size']].copy()
    bp_end.columns = ['Gene', 'Chrom', 'Pos', 'Size']
    
    # Combine
    bp_df = pd.concat([bp_start, bp_end], ignore_index=True)
    
    # Store gene sizes (vectorized)
    gene_sizes = inv_df.set_index('Gene')['Size'].to_dict()
    
    # Fast parallel DBS counting
    gene_dbs = fast_dbs_count(bp_df, dbs_df, window=window)
    
    # Create result DataFrame
    results = []
    for gene in gene_sizes.keys():
        dbs_count = gene_dbs.get(gene, 0)
        results.append({
            'Gene': gene,
            'Size': gene_sizes[gene],
            'DBS_Count': dbs_count,
            'Has_DBS': 1 if dbs_count > 0 else 0
        })
    
    result_df = pd.DataFrame(results)
    
    print(f"  Genes analyzed: {len(result_df):,}")
    print(f"  Genes with DBS: {result_df['Has_DBS'].sum():,}")
    
    return result_df


def classify_by_size(df):
    """Classify genes by max INV size."""
    print("\nClassifying by size...")
    
    # Get max size per gene
    gene_max_size = df.groupby('Gene')['Size'].max().reset_index()
    gene_max_size.columns = ['Gene', 'Max_Size']
    
    # Merge with DBS info
    gene_dbs = df.groupby('Gene').agg({
        'DBS_Count': 'sum',
        'Has_DBS': 'max'
    }).reset_index()
    
    merged = gene_max_size.merge(gene_dbs, on='Gene')
    
    # Classify
    def get_size_class(size):
        for name, min_s, max_s in SIZE_CLASSES:
            if min_s <= size < max_s:
                return name
        return 'Unknown'
    
    merged['Size_Class'] = merged['Max_Size'].apply(get_size_class)
    
    return merged


def analyze_by_size_class(df):
    """Analyze DBS enrichment by size class."""
    print("\nAnalyzing by size class...")
    
    results = []
    
    for name, min_size, max_size in SIZE_CLASSES:
        class_genes = df[df['Size_Class'] == name]
        
        if len(class_genes) == 0:
            continue
        
        n_genes = len(class_genes)
        n_with_dbs = class_genes['Has_DBS'].sum()
        pct_with_dbs = 100 * n_with_dbs / n_genes
        
        mean_dbs = class_genes['DBS_Count'].mean()
        total_dbs = class_genes['DBS_Count'].sum()
        
        results.append({
            'Size_Class': name,
            'Min_Size_Mb': min_size / 1e6,
            'Max_Size_Mb': max_size / 1e6 if max_size != float('inf') else np.inf,
            'N_Genes': n_genes,
            'N_With_DBS': n_with_dbs,
            'Percent_With_DBS': pct_with_dbs,
            'Total_DBS': total_dbs,
            'Mean_DBS_Per_Gene': mean_dbs
        })
        
        print(f"\n  {name}:")
        print(f"    Genes: {n_genes:,}")
        print(f"    With DBS: {n_with_dbs:,} ({pct_with_dbs:.1f}%)")
        print(f"    Mean DBS/gene: {mean_dbs:.2f}")
    
    return pd.DataFrame(results)


def test_mega_vs_small(df, threshold=50_000_000):
    """Statistical test: mega vs smaller inversions."""
    print(f"\n{'='*80}")
    print(f"MEGA-INVERSION TEST (threshold = {threshold/1e6:.0f} Mb)")
    print(f"{'='*80}")
    
    mega = df[df['Max_Size'] >= threshold]
    small = df[df['Max_Size'] < threshold]
    
    mega_dbs = mega['Has_DBS'].sum()
    small_dbs = small['Has_DBS'].sum()
    
    mega_pct = 100 * mega_dbs / len(mega) if len(mega) > 0 else 0
    small_pct = 100 * small_dbs / len(small) if len(small) > 0 else 0
    
    print(f"\nMega-inversions (≥{threshold/1e6:.0f} Mb):")
    print(f"  N genes: {len(mega):,}")
    print(f"  With DBS: {mega_dbs:,} ({mega_pct:.1f}%)")
    
    print(f"\nSmaller inversions (<{threshold/1e6:.0f} Mb):")
    print(f"  N genes: {len(small):,}")
    print(f"  With DBS: {small_dbs:,} ({small_pct:.1f}%)")
    
    # Fisher's exact test
    if len(mega) > 0 and len(small) > 0:
        contingency = [
            [mega_dbs, len(mega) - mega_dbs],
            [small_dbs, len(small) - small_dbs]
        ]
        
        odds_ratio, p_value = stats.fisher_exact(contingency)
        
        print(f"\nFisher's exact test:")
        print(f"  Odds ratio: {odds_ratio:.2f}×")
        print(f"  P-value: {p_value:.2e}")
        
        if p_value < 0.05:
            if odds_ratio > 1:
                print(f"  ✓ Mega-inversions SIGNIFICANTLY more enriched for DBS")
            else:
                print(f"  ✓ Smaller inversions SIGNIFICANTLY more enriched for DBS")
        else:
            print(f"  ✗ No significant difference")
        
        # Effect size
        enrichment = mega_pct / small_pct if small_pct > 0 else 0
        print(f"  Enrichment: {enrichment:.2f}× (mega vs small)")
        
        return {
            'Threshold_Mb': threshold / 1e6,
            'Mega_N': len(mega),
            'Mega_Pct': mega_pct,
            'Small_N': len(small),
            'Small_Pct': small_pct,
            'Odds_Ratio': odds_ratio,
            'P_Value': p_value,
            'Enrichment': enrichment
        }
    
    return None


def create_visualization(size_results, gene_df, mega_stats, output):
    """Create 4-panel visualization."""
    print("\nCreating visualization...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # =========================================================================
    # Panel A: DBS % by Size Class
    # =========================================================================
    ax = axes[0, 0]
    
    # Filter to classes with data
    plot_data = size_results[size_results['N_Genes'] > 0].copy()
    
    x = np.arange(len(plot_data))
    colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(plot_data)))
    
    bars = ax.bar(x, plot_data['Percent_With_DBS'],
                  color=colors, edgecolor='black', linewidth=1.5, alpha=0.85)
    
    ax.set_xticks(x)
    ax.set_xticklabels(plot_data['Size_Class'], rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('% Genes with INV-DBS', fontsize=12, fontweight='bold')
    ax.set_title('A. DBS Enrichment by Inversion Size', fontsize=13, fontweight='bold', pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)
    
    # Add percentages
    for i, (bar, pct) in enumerate(zip(bars, plot_data['Percent_With_DBS'])):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
               f'{pct:.1f}%', ha='center', va='bottom',
               fontsize=9, fontweight='bold')
    
    # =========================================================================
    # Panel B: Gene Counts by Size Class
    # =========================================================================
    ax = axes[0, 1]
    
    bars = ax.bar(x, plot_data['N_Genes'],
                  color=colors, edgecolor='black', linewidth=1.5, alpha=0.85)
    
    ax.set_xticks(x)
    ax.set_xticklabels(plot_data['Size_Class'], rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('B. Gene Distribution by Size', fontsize=13, fontweight='bold', pad=15)
    ax.set_yscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, which='both')
    
    # =========================================================================
    # Panel C: Size Distribution (With vs Without DBS)
    # =========================================================================
    ax = axes[1, 0]
    
    with_dbs = gene_df[gene_df['Has_DBS'] == 1]['Max_Size'].values / 1e6
    without_dbs = gene_df[gene_df['Has_DBS'] == 0]['Max_Size'].values / 1e6
    
    ax.hist(without_dbs, bins=50, range=(0, 200), alpha=0.5,
           label=f'Without DBS (n={len(without_dbs):,})',
           color=COLORS['Without_DBS'], edgecolor='black', linewidth=0.5)
    
    ax.hist(with_dbs, bins=50, range=(0, 200), alpha=0.7,
           label=f'With DBS (n={len(with_dbs):,})',
           color=COLORS['With_DBS'], edgecolor='black', linewidth=0.5)
    
    ax.axvline(50, color='red', linestyle='--', linewidth=2, alpha=0.7,
              label='Mega threshold (50 Mb)')
    
    ax.set_xlabel('Max Inversion Size (Mb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title('C. Size Distribution: With vs Without DBS', fontsize=13, fontweight='bold', pad=15)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)
    
    # =========================================================================
    # Panel D: Mega vs Small Comparison
    # =========================================================================
    ax = axes[1, 1]
    
    if mega_stats:
        categories = ['Mega\n(≥50 Mb)', 'Smaller\n(<50 Mb)']
        percentages = [mega_stats['Mega_Pct'], mega_stats['Small_Pct']]
        counts = [mega_stats['Mega_N'], mega_stats['Small_N']]
        
        bars = ax.bar(categories, percentages,
                     color=[COLORS['With_DBS'], COLORS['Without_DBS']],
                     edgecolor='black', linewidth=2, alpha=0.85)
        
        ax.set_ylabel('% Genes with DBS', fontsize=12, fontweight='bold')
        ax.set_title('D. Mega-Inversion Enrichment', fontsize=13, fontweight='bold', pad=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3)
        
        # Add values
        for i, (bar, pct, n) in enumerate(zip(bars, percentages, counts)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   f'{pct:.1f}%\n(n={n:,})',
                   ha='center', va='bottom',
                   fontsize=10, fontweight='bold')
        
        # Add statistics
        stats_text = (f"Enrichment: {mega_stats['Enrichment']:.2f}×\n"
                     f"OR: {mega_stats['Odds_Ratio']:.2f}\n"
                     f"p = {mega_stats['P_Value']:.2e}")
        
        ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✓ Saved: {output}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Direct inversion size analysis from AnnotSV + DBS'
    )
    parser.add_argument('--annotsv-dir', required=True,
                       help='Directory with AnnotSV TSV files')
    parser.add_argument('--dbs-dir', required=True,
                       help='Directory with DBS CSV files')
    parser.add_argument('--window', type=int, default=10,
                       help='Window size for pairing (bp)')
    parser.add_argument('--output', default='inversion_size_direct.png',
                       help='Output figure filename')
    parser.add_argument('--mega-threshold', type=int, default=50_000_000,
                       help='Threshold for mega-inversions (bp)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("DIRECT INVERSION SIZE ANALYSIS")
    print("="*80)
    
    # Load INV data
    inv_df = load_inv_from_annotsv(args.annotsv_dir)
    if len(inv_df) == 0:
        print("ERROR: No INV data loaded")
        return 1
    
    # Load DBS data
    dbs_df = load_dbs_mutations(args.dbs_dir)
    if len(dbs_df) == 0:
        print("ERROR: No DBS data loaded")
        return 1
    
    # Find DBS near INV
    paired_df = find_dbs_near_inv(inv_df, dbs_df, args.window)
    
    # Classify by size
    gene_df = classify_by_size(paired_df)
    
    # Analyze by size class
    size_results = analyze_by_size_class(gene_df)
    
    # Test mega vs small
    mega_stats = test_mega_vs_small(gene_df, args.mega_threshold)
    
    # Save results
    output_dir = Path(args.output).parent
    if output_dir == Path('.'):
        output_dir = Path('.')
    
    size_results.to_csv(output_dir / 'size_class_results.csv', index=False)
    gene_df.to_csv(output_dir / 'genes_by_size.csv', index=False)
    
    if mega_stats:
        pd.DataFrame([mega_stats]).to_csv(output_dir / 'mega_test_results.csv', index=False)
    
    print(f"\n✓ Saved: {output_dir / 'size_class_results.csv'}")
    print(f"✓ Saved: {output_dir / 'genes_by_size.csv'}")
    
    # Create visualization
    create_visualization(size_results, gene_df, mega_stats, args.output)
    
    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)
    
    return 0


if __name__ == "__main__":
    exit(main())