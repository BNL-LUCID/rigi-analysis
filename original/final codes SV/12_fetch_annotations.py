#!/usr/bin/env python3
"""
Fetch gene functional annotations from MyGene.info API.
Updated version - annotates ALL protein-coding genes (no pLI filter).

This script programmatically queries MyGene.info to get:
- Gene full name
- Function summary
- GO biological process terms
- Pathway involvement

Run this LOCALLY (not in Claude environment) as external APIs are blocked there.

Usage:
    python fetch_gene_annotations_v2.py \
        --genes-high genes_high_dose.csv \
        --genes-low genes_low_dose.csv \
        --gnomad gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
        --output annotated_genes_v4.csv
"""

import argparse
import pandas as pd
import requests
import time
from pathlib import Path


def is_protein_coding(gene):
    """Check if gene is likely protein-coding based on naming conventions."""
    gene = str(gene)
    if gene.startswith(('LOC', 'LINC', 'MIR', 'SNORD', 'SNORA')):
        return False
    if '-AS' in gene or '-DT' in gene or '-IT' in gene:
        return False
    if gene.startswith('RNA') and gene[3:].isdigit():
        return False
    return True


def load_gnomad_pli(gnomad_path):
    """Load pLI scores from gnomAD file."""
    print(f"Loading gnomAD pLI scores from {gnomad_path}...")
    
    try:
        # Read the bgz file
        gnomad = pd.read_csv(gnomad_path, sep='\t', compression='gzip',
                            usecols=['gene', 'pLI'])
        
        # Take max pLI per gene (multiple transcripts)
        gnomad_pli = gnomad.groupby('gene')['pLI'].max().reset_index()
        print(f"  Loaded pLI for {len(gnomad_pli)} genes")
        return gnomad_pli
    except Exception as e:
        print(f"  Warning: Could not load gnomAD file: {e}")
        return None


def query_mygene(gene_symbols, batch_size=100):
    """
    Query MyGene.info for multiple genes at once.
    
    Args:
        gene_symbols: List of gene symbols
        batch_size: Number of genes per API call
    
    Returns:
        Dictionary mapping gene symbol to annotation
    """
    base_url = "https://mygene.info/v3/query"
    
    all_results = {}
    
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i+batch_size]
        query = ",".join(batch)
        
        params = {
            "q": query,
            "scopes": "symbol",
            "species": "human",
            "fields": "symbol,name,summary,go.BP,pathway.kegg,pathway.reactome",
            "size": batch_size,
        }
        
        try:
            response = requests.post(base_url, data=params, timeout=30)
            response.raise_for_status()
            results = response.json()
            
            for hit in results:
                if isinstance(hit, dict) and 'symbol' in hit:
                    symbol = hit['symbol']
                    all_results[symbol] = {
                        'Full_Name': hit.get('name', ''),
                        'Summary': hit.get('summary', ''),
                        'GO_Terms': extract_go_terms(hit.get('go', {})),
                        'Pathways': extract_pathways(hit),
                    }
            
            print(f"  Fetched batch {i//batch_size + 1}/{(len(gene_symbols)-1)//batch_size + 1}: {len(batch)} genes")
            time.sleep(0.3)  # Rate limiting
            
        except Exception as e:
            print(f"  Error fetching batch {i//batch_size + 1}: {e}")
    
    return all_results


def extract_go_terms(go_data):
    """Extract GO biological process terms."""
    if not go_data:
        return ""
    
    bp_terms = go_data.get('BP', [])
    if isinstance(bp_terms, dict):
        bp_terms = [bp_terms]
    
    terms = []
    for bp in bp_terms[:10]:  # Top 10 terms
        if isinstance(bp, dict) and 'term' in bp:
            terms.append(bp['term'])
    
    return "; ".join(terms)


def extract_pathways(hit):
    """Extract pathway information."""
    pathways = []
    
    # KEGG pathways
    kegg = hit.get('pathway', {}).get('kegg', [])
    if isinstance(kegg, dict):
        kegg = [kegg]
    for p in kegg[:5]:
        if isinstance(p, dict) and 'name' in p:
            pathways.append(p['name'])
    
    # Reactome pathways  
    reactome = hit.get('pathway', {}).get('reactome', [])
    if isinstance(reactome, dict):
        reactome = [reactome]
    for p in reactome[:5]:
        if isinstance(p, dict) and 'name' in p:
            pathways.append(p['name'])
    
    return "; ".join(pathways[:8])


def main():
    parser = argparse.ArgumentParser(
        description='Fetch gene annotations from MyGene.info for all protein-coding genes'
    )
    parser.add_argument('--genes-high', required=True, 
                        help='High dose genes CSV')
    parser.add_argument('--genes-low', required=True,
                        help='Low dose genes CSV')
    parser.add_argument('--gnomad', default=None,
                        help='gnomAD pLI file (optional, .bgz or .csv)')
    parser.add_argument('--output', default='annotated_genes_v4.csv',
                        help='Output CSV')
    
    args = parser.parse_args()
    
    # Load gene files
    print("=" * 60)
    print("LOADING GENE FILES")
    print("=" * 60)
    
    high_df = pd.read_csv(args.genes_high)
    low_df = pd.read_csv(args.genes_low)
    
    print(f"High dose: {len(high_df)} rows, {high_df['Gene'].nunique()} unique genes")
    print(f"Low dose: {len(low_df)} rows, {low_df['Gene'].nunique()} unique genes")
    
    # Combine and deduplicate
    all_genes = pd.concat([
        high_df[['Gene']].drop_duplicates().assign(Dose='High'),
        low_df[['Gene']].drop_duplicates().assign(Dose='Low')
    ])
    
    # Handle genes appearing in both doses (keep High)
    all_genes = all_genes.drop_duplicates(subset='Gene', keep='first')
    
    print(f"\nCombined unique genes: {len(all_genes)}")
    
    # Filter to protein-coding
    all_genes['is_pc'] = all_genes['Gene'].apply(is_protein_coding)
    pc_genes = all_genes[all_genes['is_pc']].drop(columns=['is_pc']).copy()
    
    print(f"Protein-coding genes: {len(pc_genes)}")
    print(f"  High dose: {(pc_genes['Dose'] == 'High').sum()}")
    print(f"  Low dose: {(pc_genes['Dose'] == 'Low').sum()}")
    
    # Load gnomAD pLI if provided
    gnomad_pli = None
    if args.gnomad:
        gnomad_pli = load_gnomad_pli(args.gnomad)
    
    # Query MyGene.info
    print("\n" + "=" * 60)
    print("FETCHING ANNOTATIONS FROM MYGENE.INFO")
    print("=" * 60)
    
    gene_list = pc_genes['Gene'].tolist()
    annotations = query_mygene(gene_list)
    print(f"\nRetrieved annotations for {len(annotations)} genes")
    
    # Build result dataframe
    print("\n" + "=" * 60)
    print("BUILDING OUTPUT")
    print("=" * 60)
    
    results = []
    for _, row in pc_genes.iterrows():
        gene = row['Gene']
        dose = row['Dose']
        
        annot = annotations.get(gene, {})
        
        # Truncate summary
        summary = annot.get('Summary', '')
        if len(summary) > 300:
            summary = summary[:300] + '...'
        
        result = {
            'Gene': gene,
            'Dose': dose,
            'Full_Name': annot.get('Full_Name', ''),
            'Function': summary,
            'GO_Terms': annot.get('GO_Terms', ''),
            'Pathways': annot.get('Pathways', ''),
        }
        
        # Add pLI if available
        if gnomad_pli is not None:
            pli_match = gnomad_pli[gnomad_pli['gene'] == gene]
            if len(pli_match) > 0:
                result['pLI'] = pli_match['pLI'].values[0]
            else:
                result['pLI'] = None
        
        results.append(result)
    
    result_df = pd.DataFrame(results)
    
    # Sort by dose, then pLI (if available)
    if 'pLI' in result_df.columns:
        result_df = result_df.sort_values(
            by=['Dose', 'pLI'], 
            ascending=[True, False],
            na_position='last'
        )
    else:
        result_df = result_df.sort_values(by=['Dose', 'Gene'])
    
    # Save
    result_df.to_csv(args.output, index=False)
    print(f"\nSaved {len(result_df)} annotated genes to {args.output}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    print(f"\nBy dose:")
    print(f"  High: {(result_df['Dose'] == 'High').sum()} genes")
    print(f"  Low: {(result_df['Dose'] == 'Low').sum()} genes")
    
    if 'pLI' in result_df.columns:
        high_pli = result_df[result_df['pLI'] >= 0.9]
        print(f"\nHigh pLI (≥0.9): {len(high_pli)} genes")
        print(f"  High dose: {(high_pli['Dose'] == 'High').sum()}")
        print(f"  Low dose: {(high_pli['Dose'] == 'Low').sum()}")
    
    print(f"\nGenes with GO annotations: {(result_df['GO_Terms'] != '').sum()}")
    print(f"Genes with pathway annotations: {(result_df['Pathways'] != '').sum()}")
    
    # Show sample
    print("\n" + "-" * 60)
    print("SAMPLE OUTPUT (first 10 genes)")
    print("-" * 60)
    for _, row in result_df.head(10).iterrows():
        pli_str = f" pLI={row['pLI']:.2f}" if 'pLI' in row and pd.notna(row['pLI']) else ""
        print(f"\n{row['Gene']} ({row['Dose']}){pli_str}")
        if row['Full_Name']:
            print(f"  Name: {row['Full_Name'][:60]}...")
        if row['GO_Terms']:
            print(f"  GO: {row['GO_Terms'][:80]}...")


if __name__ == "__main__":
    main()