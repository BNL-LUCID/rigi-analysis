#!/usr/bin/env python3
"""
Categorize genes into condensed functional groups (6 categories).

Usage:
    python categorize_genes_condensed.py \
        --input annotated_gene.csv \
        --output categorized_genes_v4.csv
"""

import pandas as pd
import argparse

# =============================================================================
# CONDENSED CATEGORIES (6 instead of 11)
# =============================================================================

CATEGORIES = {
    'Signal Transduction': [
        'MAPK', 'Ras', 'kinase', 'signaling', 'signal transduction', 
        'receptor', 'GTPase', 'phosphorylation', 'PI3K', 'Wnt',
        'TGF', 'Notch', 'growth factor', 'GPCR', 'phosphatase',
        'tyrosine kinase', 'serine/threonine', 'G protein', 'cascade'
    ],
    'Cell Structure & Adhesion': [
        'actin', 'microtubule', 'adhesion', 'cytoskeleton', 'migration',
        'integrin', 'junction', 'motility', 'kinesin', 'myosin',
        'tubulin', 'focal adhesion', 'extracellular matrix', 'cadherin',
        'cell polarity', 'transport', 'vesicle', 'endocytosis', 'trafficking',
        'membrane', 'localization', 'secretion'
    ],
    'Gene Expression': [
        'chromatin', 'histone', 'transcription', 'DNA binding', 
        'RNA polymerase', 'nucleosome', 'epigenetic', 'methylation',
        'acetylation', 'zinc finger', 'transcription factor',
        'translation', 'mRNA', 'ribosom', 'RNA processing', 'splicing',
        'RNA binding', 'spliceosome', 'polycomb'
    ],
    'Cell Cycle & DNA Damage': [
        'cell cycle', 'mitotic', 'mitosis', 'DNA repair', 'DNA damage',
        'checkpoint', 'cytokinesis', 'chromosome', 'segregation',
        'cell division', 'apoptosis', 'cohesin', 'G1/S', 'G2/M',
        'centrosome', 'spindle', 'kinetochore', 'cell death',
        'ubiquitin', 'proteasome', 'autophagy'
    ],
    'Development & Differentiation': [
        'development', 'morphogenesis', 'differentiation', 'embryonic',
        'organogenesis', 'patterning', 'stem cell', 'neurogenesis',
        'angiogenesis', 'neuron', 'synapse', 'synaptic', 'axon',
        'dendrite', 'neural', 'brain', 'nervous system',
        'immune', 'T cell', 'B cell', 'cytokine'
    ],
    'Metabolism & Other': [
        'metaboli', 'biosynthesis', 'catabolic', 'glucose', 'lipid',
        'mitochondri', 'ATP', 'oxidative', 'respiration', 'enzyme',
        'ceramide', 'sphingolipid', 'glycosaminoglycan'
    ]
}

# Priority order for tie-breaking (if gene matches multiple categories)
CATEGORY_PRIORITY = [
    'Cell Cycle & DNA Damage',
    'Signal Transduction', 
    'Gene Expression',
    'Development & Differentiation',
    'Cell Structure & Adhesion',
    'Metabolism & Other'
]


def categorize_gene(row):
    """Assign gene to functional category based on keywords."""
    
    # Combine text fields for searching
    text_fields = []
    for col in ['GO_Terms', 'Pathways', 'Function', 'Full_Name']:
        if col in row.index and pd.notna(row[col]):
            text_fields.append(str(row[col]).lower())
    
    combined_text = ' '.join(text_fields)
    
    if not combined_text.strip():
        return 'Metabolism & Other'
    
    # Score each category
    scores = {}
    for category, keywords in CATEGORIES.items():
        score = 0
        for keyword in keywords:
            if keyword.lower() in combined_text:
                score += 1
        scores[category] = score
    
    # Get max score
    max_score = max(scores.values())
    
    if max_score == 0:
        return 'Metabolism & Other'
    
    # If tie, use priority order
    top_categories = [cat for cat, score in scores.items() if score == max_score]
    
    if len(top_categories) == 1:
        return top_categories[0]
    else:
        # Return highest priority among tied categories
        for cat in CATEGORY_PRIORITY:
            if cat in top_categories:
                return cat
        return top_categories[0]


def main():
    parser = argparse.ArgumentParser(
        description='Categorize genes into condensed functional groups'
    )
    parser.add_argument('--input', required=True,
                        help='Input annotated genes CSV')
    parser.add_argument('--output', required=True,
                        help='Output CSV with categories')
    
    args = parser.parse_args()
    
    # Load data
    print("=" * 60)
    print("GENE CATEGORIZATION (CONDENSED)")
    print("=" * 60)
    
    df = pd.read_csv(args.input)
    print(f"\nLoaded {len(df)} genes from {args.input}")
    
    # Filter out pseudogenes and non-coding
    def is_real_gene(gene):
        gene = str(gene)
        if 'pseudogene' in gene.lower():
            return False
        if gene.endswith('P1') or gene.endswith('P2') or gene.endswith('P5') or gene.endswith('P6') or gene.endswith('P9'):
            # Check if it's a pseudogene pattern like GOLGA2P6
            if any(c.isdigit() for c in gene[:-2]):
                return False
        return True
    
    # Keep track of excluded genes
    excluded = df[~df['Gene'].apply(is_real_gene)]['Gene'].tolist()
    if excluded:
        print(f"\nExcluding {len(excluded)} pseudogenes/non-coding:")
        for g in excluded:
            print(f"  - {g}")
    
    df = df[df['Gene'].apply(is_real_gene)].copy()
    print(f"\nAnalyzing {len(df)} protein-coding genes")
    
    # Categorize each gene
    df['Functional_Category'] = df.apply(categorize_gene, axis=1)
    
    # Print summary
    print("\n" + "=" * 60)
    print("CATEGORY DISTRIBUTION")
    print("=" * 60)
    
    # Overall counts
    print("\nOverall:")
    for cat in CATEGORY_PRIORITY:
        count = (df['Functional_Category'] == cat).sum()
        if count > 0:
            print(f"  {cat}: {count}")
    
    # By dose
    print("\n" + "-" * 60)
    print("BY DOSE")
    print("-" * 60)
    
    for dose in ['High', 'Low']:
        dose_df = df[df['Dose'] == dose]
        print(f"\n{dose} Dose ({len(dose_df)} genes):")
        for cat in CATEGORY_PRIORITY:
            count = (dose_df['Functional_Category'] == cat).sum()
            if count > 0:
                genes = dose_df[dose_df['Functional_Category'] == cat]['Gene'].tolist()
                print(f"  {cat}: {count}")
                for g in genes[:5]:  # Show first 5
                    pli = dose_df[dose_df['Gene'] == g]['pLI'].values[0]
                    pli_str = f" (pLI={pli:.2f})" if pd.notna(pli) and pli >= 0.5 else ""
                    print(f"    - {g}{pli_str}")
                if len(genes) > 5:
                    print(f"    ... and {len(genes) - 5} more")
    
    # Save
    df.to_csv(args.output, index=False)
    print(f"\n{'=' * 60}")
    print(f"Saved {len(df)} categorized genes to {args.output}")
    print("=" * 60)
    
    # Print cross-tabulation
    print("\n" + "=" * 60)
    print("CROSS-TABULATION (for figure)")
    print("=" * 60)
    
    cross_tab = pd.crosstab(df['Functional_Category'], df['Dose'])
    cross_tab = cross_tab.reindex(CATEGORY_PRIORITY)
    cross_tab = cross_tab.fillna(0).astype(int)
    print(cross_tab.to_string())


if __name__ == '__main__':
    main()