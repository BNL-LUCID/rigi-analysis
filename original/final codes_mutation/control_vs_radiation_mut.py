#!/usr/bin/env python3
"""
Control vs Radiation Mutation Pattern Analysis
===============================================
Compare mutation patterns between control (C) and radiation (T) conditions.

Pattern structure: [W0][W1][W2][W3]
- 0 = Absent in both
- C = Present in Control only
- B = Present in Both
- T = Present in Treatment (radiation) only

Key comparisons:
- 0T** vs 0C** patterns (radiation-induced vs control-specific)
- Gene overlap between persistent patterns
- Pathway enrichment differences

Usage:
    python control_vs_radiation_patterns.py \
        --catalog mutation_gene_catalog.csv \
        --output control_comparison \
        --plot
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# PATTERN DEFINITIONS
# =============================================================================

# Radiation-specific patterns (not in baseline, appears in radiation)
RADIATION_PATTERNS = {
    'transient': ['0T00', '00T0', '000T'],
    'persistent': ['0TT0', '0T0T', '00TT', '0TTT']
}

# Control-specific patterns (not in baseline, appears in control only)
CONTROL_PATTERNS = {
    'transient': ['0C00', '00C0', '000C'],
    'persistent': ['0CC0', '0C0C', '00CC', '0CCC']
}

# Pattern pairs for direct comparison
PATTERN_PAIRS = [
    ('0T00', '0C00', 'Early-cleared'),
    ('00T0', '00C0', 'Mid-transient'),
    ('000T', '000C', 'Late-appearing'),
    ('0TT0', '0CC0', 'Early-mid persistent'),
    ('0T0T', '0C0C', 'Intermittent'),
    ('00TT', '00CC', 'Late-persistent'),
    ('0TTT', '0CCC', 'Fully-persistent')
]

DOSES = ['dA', 'dB', 'dC', 'dD', 'dE']
DOSE_RATES = {'dA': 0.001, 'dB': 0.01, 'dC': 0.1, 'dD': 1.0, 'dE': 2.0}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_mutation_catalog(filepath):
    """Load mutation catalog with all patterns."""
    print(f"Loading: {filepath}")
    df = pd.read_csv(filepath)
    print(f"  Total records: {len(df):,}")
    print(f"  Columns: {df.columns.tolist()}")
    return df


# =============================================================================
# PATTERN COMPARISON
# =============================================================================

def compare_pattern_distributions(df, output_dir):
    """Compare radiation vs control pattern distributions."""
    print("\n" + "=" * 70)
    print("PATTERN DISTRIBUTION COMPARISON")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    # Get pattern counts
    if 'Count' in df.columns:
        pattern_counts = df.groupby('Pattern')['Count'].sum()
    else:
        pattern_counts = df['Pattern'].value_counts()
    
    # Direct pattern pair comparison
    print("\n  Direct Pattern Pair Comparison:")
    print(f"  {'Radiation':<12} {'Control':<12} {'Rad Count':<15} {'Ctrl Count':<15} {'Ratio':<10} {'Label'}")
    print("  " + "-" * 85)
    
    comparison_results = []
    
    for rad_pat, ctrl_pat, label in PATTERN_PAIRS:
        rad_count = pattern_counts.get(rad_pat, 0)
        ctrl_count = pattern_counts.get(ctrl_pat, 0)
        
        if ctrl_count > 0:
            ratio = rad_count / ctrl_count
        else:
            ratio = float('inf') if rad_count > 0 else 0
        
        print(f"  {rad_pat:<12} {ctrl_pat:<12} {rad_count:<15,} {ctrl_count:<15,} {ratio:<10.2f} {label}")
        
        comparison_results.append({
            'Radiation_Pattern': rad_pat,
            'Control_Pattern': ctrl_pat,
            'Label': label,
            'Radiation_Count': rad_count,
            'Control_Count': ctrl_count,
            'Ratio_Rad_Ctrl': ratio
        })
    
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df.to_csv(output_dir / 'pattern_pair_comparison.csv', index=False)
    
    # Total radiation-specific vs control-specific
    print("\n  Aggregate Comparison:")
    
    rad_transient = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['transient'])
    rad_persistent = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['persistent'])
    ctrl_transient = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['transient'])
    ctrl_persistent = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['persistent'])
    
    print(f"    Radiation transient (0T00+00T0+000T):  {rad_transient:,}")
    print(f"    Control transient (0C00+00C0+000C):    {ctrl_transient:,}")
    print(f"    Ratio: {rad_transient/ctrl_transient:.3f}")
    
    print(f"\n    Radiation persistent (0TT0+0T0T+00TT+0TTT): {rad_persistent:,}")
    print(f"    Control persistent (0CC0+0C0C+00CC+0CCC):   {ctrl_persistent:,}")
    print(f"    Ratio: {rad_persistent/ctrl_persistent:.3f}")
    
    # Chi-square test: Is radiation enriched for persistent patterns?
    print("\n  Chi-Square Test (Transient vs Persistent × Radiation vs Control):")
    contingency = [
        [rad_transient, rad_persistent],
        [ctrl_transient, ctrl_persistent]
    ]
    chi2, pval, dof, expected = chi2_contingency(contingency)
    print(f"    χ² = {chi2:.2f}, p = {pval:.2e}")
    
    rad_persist_rate = rad_persistent / (rad_transient + rad_persistent) * 100
    ctrl_persist_rate = ctrl_persistent / (ctrl_transient + ctrl_persistent) * 100
    print(f"    Radiation persistence rate: {rad_persist_rate:.2f}%")
    print(f"    Control persistence rate:   {ctrl_persist_rate:.2f}%")
    
    if pval < 0.05:
        if rad_persist_rate > ctrl_persist_rate:
            print("    → SIGNIFICANT: Radiation has HIGHER persistence rate")
        else:
            print("    → SIGNIFICANT: Control has HIGHER persistence rate")
    else:
        print("    → Not significant: Similar persistence rates")
    
    # Save aggregate
    aggregate = {
        'Category': ['Radiation_Transient', 'Radiation_Persistent', 
                     'Control_Transient', 'Control_Persistent'],
        'Count': [rad_transient, rad_persistent, ctrl_transient, ctrl_persistent],
        'Persistence_Rate': [0, rad_persist_rate, 0, ctrl_persist_rate]
    }
    pd.DataFrame(aggregate).to_csv(output_dir / 'aggregate_comparison.csv', index=False)
    
    return comparison_df


def compare_genes_by_pattern(df, output_dir):
    """
    Compare genes hit by radiation vs control persistent patterns.
    
    NOTE: Patterns are POSITIONAL - a gene can have BOTH 0TTT mutations at one position
    AND 0CCC mutations at another position. These are independent events.
    Gene lists should NOT subtract overlaps - each pattern's genes are valid for enrichment.
    """
    print("\n" + "=" * 70)
    print("GENE OVERLAP ANALYSIS")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    # Get genes for each persistent pattern
    rad_persistent_patterns = RADIATION_PATTERNS['persistent']
    ctrl_persistent_patterns = CONTROL_PATTERNS['persistent']
    
    rad_genes = set(df[df['Pattern'].isin(rad_persistent_patterns)]['Gene'].dropna().unique())
    ctrl_genes = set(df[df['Pattern'].isin(ctrl_persistent_patterns)]['Gene'].dropna().unique())
    
    print(f"\n  Radiation persistent genes: {len(rad_genes):,}")
    print(f"  Control persistent genes:   {len(ctrl_genes):,}")
    
    # Overlap statistics (for comparison, NOT for subtraction)
    overlap = rad_genes & ctrl_genes
    
    print(f"\n  Genes with BOTH rad and ctrl persistent mutations: {len(overlap):,}")
    print(f"    (at DIFFERENT positions - these are independent events)")
    
    # Jaccard similarity
    jaccard = len(overlap) / len(rad_genes | ctrl_genes) if len(rad_genes | ctrl_genes) > 0 else 0
    print(f"\n  Jaccard similarity: {jaccard:.3f}")
    
    # Save COMPLETE gene lists (no subtraction - patterns are positional)
    pd.DataFrame({'Gene': list(rad_genes)}).to_csv(
        output_dir / 'genes_radiation_persistent.txt', index=False, header=False)
    pd.DataFrame({'Gene': list(ctrl_genes)}).to_csv(
        output_dir / 'genes_control_persistent.txt', index=False, header=False)
    
    # Compare 0TTT vs 0CCC specifically
    print("\n  Fully Persistent Comparison (0TTT vs 0CCC):")
    
    genes_0TTT = set(df[df['Pattern'] == '0TTT']['Gene'].dropna().unique())
    genes_0CCC = set(df[df['Pattern'] == '0CCC']['Gene'].dropna().unique())
    
    print(f"    0TTT genes: {len(genes_0TTT):,}")
    print(f"    0CCC genes: {len(genes_0CCC):,}")
    
    overlap_full = genes_0TTT & genes_0CCC
    print(f"    Genes with BOTH patterns: {len(overlap_full):,}")
    print(f"    Jaccard: {len(overlap_full) / len(genes_0TTT | genes_0CCC):.3f}")
    
    # Save COMPLETE gene lists for pathway enrichment
    # These are the correct files - ALL genes with each pattern
    pd.DataFrame({'Gene': list(genes_0TTT)}).to_csv(
        output_dir / 'genes_0TTT.txt', index=False, header=False)
    pd.DataFrame({'Gene': list(genes_0CCC)}).to_csv(
        output_dir / 'genes_0CCC.txt', index=False, header=False)
    
    print(f"\n  Saved gene lists for pathway analysis:")
    print(f"    - genes_0TTT.txt ({len(genes_0TTT):,} genes) - ALL radiation persistent")
    print(f"    - genes_0CCC.txt ({len(genes_0CCC):,} genes) - ALL control persistent")
    print(f"\n  NOTE: Use these COMPLETE lists for enrichment (not subtracted)")
    
    return {
        'rad_genes': rad_genes,
        'ctrl_genes': ctrl_genes,
        'overlap': overlap,
        'genes_0TTT': genes_0TTT,
        'genes_0CCC': genes_0CCC
    }


def compare_by_dose(df, output_dir):
    """Compare radiation vs control patterns by dose."""
    print("\n" + "=" * 70)
    print("DOSE-SPECIFIC COMPARISON")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    if 'Dose' not in df.columns:
        print("  No Dose column found")
        return None
    
    print("\n  Radiation vs Control Mutations by Dose:")
    print(f"  {'Dose':<8} {'Rad Trans':<12} {'Rad Pers':<12} {'Ctrl Trans':<12} {'Ctrl Pers':<12} {'Rad/Ctrl':<10}")
    print("  " + "-" * 70)
    
    results = []
    
    for dose in DOSES:
        dose_df = df[df['Dose'] == dose]
        
        if 'Count' in dose_df.columns:
            pattern_counts = dose_df.groupby('Pattern')['Count'].sum()
        else:
            pattern_counts = dose_df['Pattern'].value_counts()
        
        rad_trans = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['transient'])
        rad_pers = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['persistent'])
        ctrl_trans = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['transient'])
        ctrl_pers = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['persistent'])
        
        rad_total = rad_trans + rad_pers
        ctrl_total = ctrl_trans + ctrl_pers
        ratio = rad_total / ctrl_total if ctrl_total > 0 else 0
        
        print(f"  {dose:<8} {rad_trans:<12,} {rad_pers:<12,} {ctrl_trans:<12,} {ctrl_pers:<12,} {ratio:<10.2f}")
        
        results.append({
            'Dose': dose,
            'Dose_Rate': DOSE_RATES[dose],
            'Rad_Transient': rad_trans,
            'Rad_Persistent': rad_pers,
            'Ctrl_Transient': ctrl_trans,
            'Ctrl_Persistent': ctrl_pers,
            'Rad_Total': rad_total,
            'Ctrl_Total': ctrl_total,
            'Ratio_Rad_Ctrl': ratio
        })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_dir / 'dose_comparison.csv', index=False)
    
    # Dose-response for ratio
    from scipy.stats import spearmanr
    log_doses = [np.log10(DOSE_RATES[d]) for d in results_df['Dose']]
    ratios = results_df['Ratio_Rad_Ctrl'].values
    
    rho, pval = spearmanr(log_doses, ratios)
    print(f"\n  Dose-Response (Rad/Ctrl ratio vs log Dose):")
    print(f"    Spearman ρ = {rho:.3f}, p = {pval:.3f}")
    
    if pval < 0.05:
        if rho > 0:
            print("    → SIGNIFICANT: Higher doses have MORE radiation-specific mutations")
        else:
            print("    → SIGNIFICANT: Higher doses have FEWER radiation-specific mutations")
    
    return results_df


def compare_mutation_types(df, output_dir):
    """Compare mutation types in radiation vs control patterns."""
    print("\n" + "=" * 70)
    print("MUTATION TYPE COMPARISON")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    if 'Mutation_Type' not in df.columns:
        print("  No Mutation_Type column found")
        return None
    
    # Radiation patterns
    rad_df = df[df['Pattern'].isin(RADIATION_PATTERNS['transient'] + RADIATION_PATTERNS['persistent'])]
    ctrl_df = df[df['Pattern'].isin(CONTROL_PATTERNS['transient'] + CONTROL_PATTERNS['persistent'])]
    
    if 'Count' in df.columns:
        rad_types = rad_df.groupby('Mutation_Type')['Count'].sum()
        ctrl_types = ctrl_df.groupby('Mutation_Type')['Count'].sum()
    else:
        rad_types = rad_df['Mutation_Type'].value_counts()
        ctrl_types = ctrl_df['Mutation_Type'].value_counts()
    
    print("\n  Mutation Type Distribution:")
    print(f"  {'Type':<12} {'Radiation':<15} {'Rad %':<10} {'Control':<15} {'Ctrl %':<10} {'Ratio'}")
    print("  " + "-" * 75)
    
    all_types = set(rad_types.index) | set(ctrl_types.index)
    rad_total = rad_types.sum()
    ctrl_total = ctrl_types.sum()
    
    results = []
    for mut_type in sorted(all_types):
        rad_count = rad_types.get(mut_type, 0)
        ctrl_count = ctrl_types.get(mut_type, 0)
        rad_pct = 100 * rad_count / rad_total if rad_total > 0 else 0
        ctrl_pct = 100 * ctrl_count / ctrl_total if ctrl_total > 0 else 0
        ratio = rad_pct / ctrl_pct if ctrl_pct > 0 else 0
        
        print(f"  {mut_type:<12} {rad_count:<15,} {rad_pct:<10.1f} {ctrl_count:<15,} {ctrl_pct:<10.1f} {ratio:<.2f}")
        
        results.append({
            'Mutation_Type': mut_type,
            'Radiation_Count': rad_count,
            'Radiation_Pct': rad_pct,
            'Control_Count': ctrl_count,
            'Control_Pct': ctrl_pct,
            'Ratio': ratio
        })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_dir / 'mutation_type_comparison.csv', index=False)
    
    # Chi-square test
    print("\n  Chi-Square Test (Mutation Type × Radiation vs Control):")
    contingency = pd.DataFrame({
        'Radiation': rad_types,
        'Control': ctrl_types
    }).fillna(0)
    
    chi2, pval, dof, expected = chi2_contingency(contingency)
    print(f"    χ² = {chi2:.2f}, p = {pval:.2e}")
    
    if pval < 0.05:
        print("    → SIGNIFICANT: Mutation type distribution differs between radiation and control")
    
    return results_df


def analyze_late_crisis(df, output_dir):
    """
    Compare late-appearing patterns (000T vs 000C) - the crisis pattern.
    
    NOTE: Patterns are POSITIONAL - no subtraction of gene lists.
    """
    print("\n" + "=" * 70)
    print("LATE CRISIS ANALYSIS (000T vs 000C)")
    print("=" * 70)
    
    output_dir = Path(output_dir)
    
    # Get counts
    if 'Count' in df.columns:
        pattern_counts = df.groupby('Pattern')['Count'].sum()
    else:
        pattern_counts = df['Pattern'].value_counts()
    
    count_000T = pattern_counts.get('000T', 0)
    count_000C = pattern_counts.get('000C', 0)
    
    print(f"\n  000T (radiation late-appearing): {count_000T:,}")
    print(f"  000C (control late-appearing):   {count_000C:,}")
    print(f"  Ratio (000T/000C): {count_000T/count_000C:.3f}")
    
    # By dose
    if 'Dose' in df.columns:
        print("\n  Late Crisis by Dose:")
        print(f"  {'Dose':<8} {'000T':<15} {'000C':<15} {'Ratio':<10}")
        print("  " + "-" * 50)
        
        for dose in DOSES:
            dose_df = df[df['Dose'] == dose]
            if 'Count' in dose_df.columns:
                pc = dose_df.groupby('Pattern')['Count'].sum()
            else:
                pc = dose_df['Pattern'].value_counts()
            
            t_count = pc.get('000T', 0)
            c_count = pc.get('000C', 0)
            ratio = t_count / c_count if c_count > 0 else 0
            
            print(f"  {dose:<8} {t_count:<15,} {c_count:<15,} {ratio:<10.2f}")
    
    # Gene comparison
    genes_000T = set(df[df['Pattern'] == '000T']['Gene'].dropna().unique())
    genes_000C = set(df[df['Pattern'] == '000C']['Gene'].dropna().unique())
    
    overlap = genes_000T & genes_000C
    
    print(f"\n  Gene Comparison:")
    print(f"    000T genes: {len(genes_000T):,}")
    print(f"    000C genes: {len(genes_000C):,}")
    print(f"    Genes with BOTH patterns: {len(overlap):,}")
    print(f"    Jaccard: {len(overlap)/len(genes_000T | genes_000C):.3f}")
    
    # Save COMPLETE gene lists (no subtraction - patterns are positional)
    pd.DataFrame({'Gene': list(genes_000T)}).to_csv(
        output_dir / 'genes_000T_late_crisis.txt', index=False, header=False)
    pd.DataFrame({'Gene': list(genes_000C)}).to_csv(
        output_dir / 'genes_000C_late_control.txt', index=False, header=False)
    
    print(f"\n  Saved: genes_000T_late_crisis.txt ({len(genes_000T):,} genes)")
    print(f"  Saved: genes_000C_late_control.txt ({len(genes_000C):,} genes)")
    
    return {
        'genes_000T': genes_000T,
        'genes_000C': genes_000C,
        'overlap': overlap
    }


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_plots(df, comparison_df, output_dir):
    """Create comparison plots."""
    import matplotlib.pyplot as plt
    
    output_dir = Path(output_dir)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Pattern pair comparison
    ax = axes[0, 0]
    if comparison_df is not None:
        x = range(len(comparison_df))
        width = 0.35
        
        ax.bar([i - width/2 for i in x], comparison_df['Radiation_Count'] / 1e6, 
               width, label='Radiation', color='#e74c3c')
        ax.bar([i + width/2 for i in x], comparison_df['Control_Count'] / 1e6,
               width, label='Control', color='#3498db')
        
        ax.set_xticks(x)
        ax.set_xticklabels(comparison_df['Label'], rotation=45, ha='right')
        ax.set_ylabel('Count (millions)')
        ax.set_title('Radiation vs Control Pattern Counts', fontweight='bold')
        ax.legend()
    
    # 2. Ratio by pattern
    ax = axes[0, 1]
    if comparison_df is not None:
        ratios = comparison_df['Ratio_Rad_Ctrl'].values
        colors = ['#e74c3c' if r > 1 else '#3498db' for r in ratios]
        
        ax.barh(range(len(comparison_df)), ratios, color=colors)
        ax.set_yticks(range(len(comparison_df)))
        ax.set_yticklabels(comparison_df['Label'])
        ax.axvline(x=1, color='black', linestyle='--', linewidth=2)
        ax.set_xlabel('Ratio (Radiation / Control)')
        ax.set_title('Radiation vs Control Ratio', fontweight='bold')
    
    # 3. Persistence comparison
    ax = axes[1, 0]
    if 'Count' in df.columns:
        pattern_counts = df.groupby('Pattern')['Count'].sum()
    else:
        pattern_counts = df['Pattern'].value_counts()
    
    rad_trans = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['transient'])
    rad_pers = sum(pattern_counts.get(p, 0) for p in RADIATION_PATTERNS['persistent'])
    ctrl_trans = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['transient'])
    ctrl_pers = sum(pattern_counts.get(p, 0) for p in CONTROL_PATTERNS['persistent'])
    
    categories = ['Radiation\nTransient', 'Radiation\nPersistent', 
                  'Control\nTransient', 'Control\nPersistent']
    values = [rad_trans/1e6, rad_pers/1e6, ctrl_trans/1e6, ctrl_pers/1e6]
    colors = ['#e74c3c', '#c0392b', '#3498db', '#2980b9']
    
    ax.bar(categories, values, color=colors)
    ax.set_ylabel('Count (millions)')
    ax.set_title('Transient vs Persistent Mutations', fontweight='bold')
    
    # 4. Dose comparison
    ax = axes[1, 1]
    if 'Dose' in df.columns:
        dose_data = []
        for dose in DOSES:
            dose_df = df[df['Dose'] == dose]
            if 'Count' in dose_df.columns:
                pc = dose_df.groupby('Pattern')['Count'].sum()
            else:
                pc = dose_df['Pattern'].value_counts()
            
            rad_total = sum(pc.get(p, 0) for p in RADIATION_PATTERNS['transient'] + RADIATION_PATTERNS['persistent'])
            ctrl_total = sum(pc.get(p, 0) for p in CONTROL_PATTERNS['transient'] + CONTROL_PATTERNS['persistent'])
            dose_data.append((dose, rad_total, ctrl_total))
        
        x = range(len(DOSES))
        width = 0.35
        
        ax.bar([i - width/2 for i in x], [d[1]/1e6 for d in dose_data], 
               width, label='Radiation', color='#e74c3c')
        ax.bar([i + width/2 for i in x], [d[2]/1e6 for d in dose_data],
               width, label='Control', color='#3498db')
        
        ax.set_xticks(x)
        ax.set_xticklabels([f"{d}\n({DOSE_RATES[d]})" for d in DOSES])
        ax.set_ylabel('Count (millions)')
        ax.set_xlabel('Dose (mGy/hr)')
        ax.set_title('Mutations by Dose', fontweight='bold')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'control_vs_radiation_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n  Saved: {output_dir / 'control_vs_radiation_comparison.png'}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Control vs Radiation Pattern Analysis')
    parser.add_argument('--catalog', '-c', required=True, help='Mutation gene catalog')
    parser.add_argument('--output', '-o', default='control_comparison', help='Output directory')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("CONTROL vs RADIATION MUTATION PATTERN ANALYSIS")
    print("=" * 70)
    
    # Load data
    df = load_mutation_catalog(args.catalog)
    
    # Run analyses
    comparison_df = compare_pattern_distributions(df, output_dir)
    gene_results = compare_genes_by_pattern(df, output_dir)
    dose_results = compare_by_dose(df, output_dir)
    type_results = compare_mutation_types(df, output_dir)
    crisis_results = analyze_late_crisis(df, output_dir)
    
    # Create plots
    if args.plot:
        print("\nGenerating plots...")
        create_plots(df, comparison_df, output_dir)
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutput files saved to: {output_dir}")
    print("\nGene lists for pathway analysis:")
    print("  - genes_0TTT_only.txt (radiation-specific persistent)")
    print("  - genes_0CCC_only.txt (control-specific persistent)")
    print("  - genes_000T_only_late_crisis.txt (radiation-specific late)")
    print("  - genes_000C_only_late_control.txt (control-specific late)")


if __name__ == "__main__":
    main()