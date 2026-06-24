#!/usr/bin/env python
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import argparse
from pathlib import Path
import re

# Define pattern groups
def categorize_pattern(pattern):
    """
    Categorize patterns into biologically meaningful groups with descriptive labels.
    
    Using only the specified radiation patterns.
    """
    # Specific radiation patterns as requested
    radiation_patterns = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']
    
    if pattern in radiation_patterns:
        return 'Radiation-specific'
    # Control-specific mutations
    elif pattern in ['C000', '0C00', '00C0', '000C']:
        return 'Control-specific'
    # Baseline mutations (present in both control and radiation)
    elif pattern == 'CBBB':
        return 'Baseline (CBBB)'
    # Other patterns
    else:
        return 'Other patterns'
def _locate_combined_files(annotated_dir, pattern_dir, mutation_type, dose):
    """
    Locate single combined-chromosomes annotated and pattern files.

    Returns (annotated_path, pattern_path) or (None, None) if combined-mode
    inputs aren't present (caller falls back to per-chromosome mode).
    """
    annotated_candidates = [
        os.path.join(annotated_dir, f"all_{mutation_type}_annotated.csv"),
        os.path.join(annotated_dir, mutation_type, f"all_{mutation_type}_annotated.csv"),
        os.path.join(annotated_dir, f"all_{mutation_type}_annotated2.csv"),
        os.path.join(annotated_dir, mutation_type, f"all_{mutation_type}_annotated2.csv"),
        # MNS-style fallback (Step 4 default output name when not renamed)
        os.path.join(annotated_dir, mutation_type, "annotated_mutations.csv"),
    ]
    annotated_path = next((p for p in annotated_candidates if os.path.exists(p)), None)
    if annotated_path is None:
        # Glob fallback: any single non-per-chromosome file (incl. type subdir)
        glob_patterns = [
            os.path.join(annotated_dir, f"*_{mutation_type}_annotated*.csv"),
            os.path.join(annotated_dir, mutation_type, f"*_{mutation_type}_annotated*.csv"),
        ]
        candidates = []
        for gp in glob_patterns:
            candidates.extend(glob.glob(gp))
        candidates = [
            p for p in candidates
            if not re.match(r'(\d+|[XY])_', os.path.basename(p))
        ]
        if len(candidates) == 1:
            annotated_path = candidates[0]

    pattern_candidates = [
        os.path.join(pattern_dir, f"dose_{dose}", f"mutation_annotations_dose_{dose}_{mutation_type}.csv"),
        os.path.join(pattern_dir, mutation_type, f"dose_{dose}", f"mutation_annotations_dose_{dose}_{mutation_type}.csv"),
    ]
    pattern_path = next((p for p in pattern_candidates if os.path.exists(p)), None)

    return annotated_path, pattern_path


def _merge_combined(annotated_path, pattern_path, mutation_type, dose, output_dir=None):
    """Combined-mode merge: one annotated file + one pattern file → merged DataFrame."""
    print(f"Loading annotated file: {annotated_path}")
    annotated_df = pd.read_csv(annotated_path)
    print(f"  {len(annotated_df):,} rows, columns: {annotated_df.columns.tolist()}")

    print(f"Loading pattern file: {pattern_path}")
    pattern_df = pd.read_csv(pattern_path)
    print(f"  {len(pattern_df):,} rows, columns: {pattern_df.columns.tolist()}")

    if 'Pattern' not in pattern_df.columns:
        raise ValueError(f"'Pattern' column missing from {pattern_path}")
    if 'MutationID' not in annotated_df.columns:
        raise ValueError(f"'MutationID' column missing from {annotated_path}")

    # Pattern_Group derived from Pattern (downstream code expects this column)
    pattern_df['Pattern_Group'] = pattern_df['Pattern'].apply(categorize_pattern)

    # Build join keys (chr_pos_ref_alt). Annotated MutationID has trailing
    # _sample_timepoint, so trim to the first 4 components.
    annotated_df['MatchKey'] = annotated_df['MutationID'].apply(
        lambda x: '_'.join(str(x).split('_')[:4]) if pd.notna(x) else x
    )
    if 'PermanentMutationID' in pattern_df.columns:
        pattern_df['MatchKey'] = pattern_df['PermanentMutationID']
    else:
        pattern_df['MatchKey'] = pattern_df['MutationID']

    sample_anno = annotated_df['MatchKey'].iloc[0] if not annotated_df.empty else None
    sample_pat = pattern_df['MatchKey'].iloc[0] if not pattern_df.empty else None
    print(f"Sample annotated MatchKey: {sample_anno}")
    print(f"Sample pattern   MatchKey: {sample_pat}")

    pattern_cols = ['MatchKey', 'Pattern', 'Pattern_Group']
    for extra in ('W0', 'W1', 'W2', 'W3', 'Category'):
        if extra in pattern_df.columns:
            pattern_cols.append(extra)

    merged = pd.merge(
        annotated_df,
        pattern_df[pattern_cols],
        on='MatchKey',
        how='inner',
    )
    print(f"Merged: {len(merged):,} rows "
          f"(annotated={len(annotated_df):,}, pattern={len(pattern_df):,})")

    if merged.empty:
        print("Warning: merge returned 0 rows. Sample keys:")
        print(f"  annotated: {annotated_df['MatchKey'].head(3).tolist()}")
        print(f"  pattern:   {pattern_df['MatchKey'].head(3).tolist()}")

    # Save to output_dir/merged_data/{type}_dose_{dose}_merged.csv (matches what
    # main() looks for on subsequent runs). If output_dir wasn't passed, fall
    # back to a project-relative location for backward compat.
    if output_dir:
        save_dir = os.path.join(output_dir, "merged_data")
    else:
        save_dir = os.path.join("merged_data", mutation_type)
    os.makedirs(save_dir, exist_ok=True)
    output_file = os.path.join(save_dir, f"{mutation_type}_dose_{dose}_merged.csv")
    merged.to_csv(output_file, index=False)
    print(f"Saved merged data to {output_file}")
    return merged


def load_and_merge_data(annotated_dir, pattern_dir, mutation_type, dose, chromosomes=None, output_dir=None):
    """
    Load and merge annotated data with pattern data for a specific mutation type and dose.

    Tries combined-file layout first (single all-chromosomes file per side).
    Falls back to per-chromosome layout if combined files aren't present.
    """
    if chromosomes is None:
        annotated_path, pattern_path = _locate_combined_files(
            annotated_dir, pattern_dir, mutation_type, dose
        )
        if annotated_path and pattern_path:
            print(f"Combined-file mode: {annotated_path} + {pattern_path}")
            return _merge_combined(annotated_path, pattern_path, mutation_type, dose, output_dir=output_dir)
        print("Combined-file inputs not found, falling back to per-chromosome mode")

    merged_data = []

    # If no chromosomes specified, detect from available files
    if not chromosomes:
        # Find available chromosomes from pattern files
        pattern_path = f"{pattern_dir}/{mutation_type}/dose_{dose}/mutation_annotations_dose_{dose}_chr_*.csv"
        print(f"Looking for pattern files matching: {pattern_path}")
        
        pattern_files = glob.glob(pattern_path)
        if not pattern_files:
            print(f"No pattern files found for path: {pattern_path}")
            # Try alternate path format
            alt_pattern_path = f"{pattern_dir}/dose_{dose}/mutation_annotations_dose_{dose}_chr_*.csv"
            print(f"Trying alternate path: {alt_pattern_path}")
            pattern_files = glob.glob(alt_pattern_path)
            if not pattern_files:
                print(f"No pattern files found for alternate path either")
                return pd.DataFrame()
        
        print(f"Found {len(pattern_files)} pattern files")
        chromosomes = []
        for f in pattern_files:
            match = re.search(r'chr_([^.]+)', os.path.basename(f))
            if match:
                chromosomes.append(match.group(1))
    
    print(f"Processing {mutation_type} mutations for dose {dose} on chromosomes: {chromosomes}")
    
    # Check if annotated directory exists
    if not os.path.exists(annotated_dir):
        print(f"Error: Annotated directory {annotated_dir} does not exist")
        return pd.DataFrame()
    
    # Find all annotated files for this mutation type
    annotated_pattern = f"*_{mutation_type}_annotated.csv"  # Added * to match annotated2.csv files
    print(f"Searching for annotated files matching: {annotated_pattern}")
    all_annotated = glob.glob(f"{annotated_dir}/{annotated_pattern}", recursive=True)
    print(f"Found {len(all_annotated)} annotated files for {mutation_type}")
    
    # Create a mapping from chromosome to annotated file
    chr_to_file = {}
    for file_path in all_annotated:
        filename = os.path.basename(file_path)
        # Extract chromosome from filename pattern like "2_ID_annotated.csv" or "2_ID_annotated2.csv"
        match = re.match(r'(\d+|[XY])_.*_annotated.*\.csv', filename)
        if match:
            chr_name = match.group(1)
            # Prefer annotated2.csv if available (has more columns)
            if chr_name in chr_to_file:
                if "annotated" in filename:
                    chr_to_file[chr_name] = file_path
            else:
                chr_to_file[chr_name] = file_path
            print(f"Mapped chromosome {chr_name} to file {file_path}")
    
    for chrom in chromosomes:
        # Check if we have an annotated file for this chromosome
        if chrom in chr_to_file:
            annotated_file = chr_to_file[chrom]
        else:
            print(f"Warning: No annotated file found for chromosome {chrom}")
            continue
            
        print(f"Loading annotated file: {annotated_file}")
        try:
            annotated_df = pd.read_csv(annotated_file)
            print(f"Loaded {len(annotated_df)} rows from annotated file")
        except Exception as e:
            print(f"Error loading annotated file {annotated_file}: {e}")
            continue
        
        # Try different possible paths for pattern files
        pattern_paths = [
            f"{pattern_dir}/{mutation_type}/dose_{dose}/mutation_annotations_dose_{dose}_chr_{chrom}.csv",
            f"{pattern_dir}/dose_{dose}/mutation_annotations_dose_{dose}_chr_{chrom}.csv"
        ]
        
        pattern_file = None
        for path in pattern_paths:
            if os.path.exists(path):
                pattern_file = path
                break
                
        if not pattern_file:
            print(f"Warning: No matching pattern file found for chromosome {chrom}, tried paths: {pattern_paths}")
            continue
            
        print(f"Loading pattern file: {pattern_file}")
        try:
            pattern_df = pd.read_csv(pattern_file)
            print(f"Loaded {len(pattern_df)} rows from pattern file")
        except Exception as e:
            print(f"Error loading pattern file {pattern_file}: {e}")
            continue
        
        # Print column names for debugging
        print(f"Annotated file columns: {annotated_df.columns.tolist()}")
        print(f"Pattern file columns: {pattern_df.columns.tolist()}")
        
        # Check if necessary columns exist
        if 'Pattern' not in pattern_df.columns:
            print(f"Error: 'Pattern' column missing from pattern file {pattern_file}")
            continue
            
        if 'MutationID' not in annotated_df.columns:
            print(f"Error: 'MutationID' column missing from annotated file {annotated_file}")
            continue
            
        # Ensure we have a PermanentMutationID column
        if 'PermanentMutationID' not in pattern_df.columns:
            print(f"Error: 'PermanentMutationID' column missing from pattern file {pattern_file}")
            # Try to find an equivalent column
            id_columns = [col for col in pattern_df.columns if 'id' in col.lower() or 'mutation' in col.lower()]
            if id_columns:
                print(f"Found potential ID columns: {id_columns}, using {id_columns[0]}")
                pattern_df['PermanentMutationID'] = pattern_df[id_columns[0]]
            else:
                continue
        
        # Add pattern group
        pattern_df['Pattern_Group'] = pattern_df['Pattern'].apply(categorize_pattern)
        
        # Show a sample of pattern data for debugging
        print(f"Sample patterns and groups:\n{pattern_df[['Pattern', 'Pattern_Group']].value_counts().head(5)}")
        
        # ---- Create matching keys that will align between the two dataframes ----
        
        # Get sample IDs to understand format
        if not annotated_df.empty:
            sample_anno_id = annotated_df['MutationID'].iloc[0]
            print(f"Sample MutationID format: {sample_anno_id}")
        else:
            print("Warning: Annotated dataframe is empty")
            continue
            
        if not pattern_df.empty:
            sample_pattern_id = pattern_df['PermanentMutationID'].iloc[0]
            print(f"Sample PermanentMutationID format: {sample_pattern_id}")
        else:
            print("Warning: Pattern dataframe is empty")
            continue
            
        # Extract common components for matching based on the observed formats
        # For annotated format like "Y_3028275_TG_CA_d0_W0"
        # And pattern format like "Y_10071536_TG_CT"
        # We need to create matching keys with just chromosome, position, ref, alt
        
        # For annotated dataframe, create a new match key by extracting components
        annotated_df['MatchKey'] = annotated_df['MutationID'].apply(
            lambda x: '_'.join(x.split('_')[:4]) if isinstance(x, str) and '_' in x else x
        )
        print(f"Created MatchKey from MutationID by keeping only first 4 components")
        print(f"Sample original MutationID: {sample_anno_id}")
        print(f"Sample MatchKey: {annotated_df['MatchKey'].iloc[0]}")
        
        # For pattern dataframe, check if we have all 4 components
        components = sample_pattern_id.split('_')
        if len(components) >= 4:
            # The PermanentMutationID already has the format we need
            pattern_df['MatchKey'] = pattern_df['PermanentMutationID']
            print("Pattern PermanentMutationID already has the right format for matching")
        else:
            # In case it has a different format
            print("Pattern PermanentMutationID has unexpected format, trying to build matching key")
            if all(col in pattern_df.columns for col in ['Chromosome', 'Start', 'Ref', 'Alt']):
                pattern_df['MatchKey'] = pattern_df['Chromosome'] + '_' + \
                                         pattern_df['Start'].astype(str) + '_' + \
                                         pattern_df['Ref'] + '_' + \
                                         pattern_df['Alt']
                print(f"Built MatchKey from individual columns")
            else:
                print(f"Cannot create matching key: missing necessary columns in pattern file")
                continue
                
        # Quick check if the ref/alt bases might need to be reversed in one of the datasets
        # This is a simplified check - might need more sophisticated matching
        merge_key_annotated = 'MatchKey'  # Default value
        
        if len(merged_data) == 0:  # Only check on first iteration
            # Get a few examples of MatchKey from both dataframes
            anno_keys = annotated_df['MatchKey'].head(10).tolist()
            pattern_keys = pattern_df['MatchKey'].head(10).tolist()
            
            # Check for any exact matches
            exact_matches = set(anno_keys).intersection(set(pattern_keys))
            if exact_matches:
                print(f"Found {len(exact_matches)} exact matches between keys")
            else:
                print("No exact matches between keys, checking for reversed ref/alt")
                
                # Check if we need to swap ref/alt in one dataset
                # Example: Convert "Y_123_TG_CA" to "Y_123_CA_TG"
                sample_pairs = []
                for a_key in anno_keys[:5]:
                    parts = a_key.split('_')
                    if len(parts) >= 4:
                        reversed_key = f"{parts[0]}_{parts[1]}_{parts[3]}_{parts[2]}"
                        if reversed_key in pattern_keys:
                            sample_pairs.append((a_key, reversed_key))
                
                if sample_pairs:
                    print(f"Found {len(sample_pairs)} matches with reversed ref/alt!")
                    print(f"Example: {sample_pairs[0][0]} matches with {sample_pairs[0][1]}")
                    print("Creating reversed match key")
                    
                    # Create reversed key for annotated dataframe
                    parts = annotated_df['MatchKey'].str.split('_', expand=True)
                    if len(parts.columns) >= 4:
                        annotated_df['ReversedMatchKey'] = parts[0] + '_' + parts[1] + '_' + parts[3] + '_' + parts[2]
                        print("Using ReversedMatchKey for matching")
                        merge_key_annotated = 'ReversedMatchKey'
        
        # Attempt the merge
        print(f"Merging on {merge_key_annotated} to MatchKey")
        try:
            merged = pd.merge(
                annotated_df, 
                pattern_df[['MatchKey', 'Pattern', 'Pattern_Group']], 
                left_on=merge_key_annotated, 
                right_on='MatchKey',
                how='inner'
            )
            
            print(f"Merged data has {len(merged)} rows")
            
            if len(merged) == 0:
                print(f"Warning: Merge resulted in empty DataFrame for chromosome {chrom}")
                
                # Show samples for debugging
                print(f"Sample annotated match keys: {annotated_df[merge_key_annotated].head(5).tolist()}")
                print(f"Sample pattern match keys: {pattern_df['MatchKey'].head(5).tolist()}")
                
                # Try a more flexible approach - using string contains (this is just for debugging)
                print("Checking for partial string matches...")
                found = 0
                for i, anno_id in enumerate(annotated_df[merge_key_annotated].head(20)):
                    for j, pat_id in enumerate(pattern_df['MatchKey'].head(20)):
                        # Check if key components match partially
                        anno_parts = str(anno_id).split('_')
                        pat_parts = str(pat_id).split('_')
                        
                        # Check if chromosome and ref/alt match
                        if len(anno_parts) >= 4 and len(pat_parts) >= 4:
                            if anno_parts[0] == pat_parts[0] and anno_parts[2] == pat_parts[2] and anno_parts[3] == pat_parts[3]:
                                print(f"Potential position mismatch: {anno_id} vs {pat_id}")
                                found += 1
                        
                        if found >= 5:
                            break
                    if found >= 5:
                        break
                        
                if found == 0:
                    print("No partial matches found in sample IDs")
            else:
                merged_data.append(merged)
                
        except Exception as e:
            print(f"Error merging data for chromosome {chrom}: {e}")
    
    if merged_data:
        result = pd.concat(merged_data, ignore_index=True)
        print(f"Final merged data has {len(result)} rows")
        
        # Save the merged data to a CSV file for future use
        output_dir = os.path.join("merged_data", mutation_type)
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        result.to_csv(output_file, index=False)
        print(f"Saved merged data to {output_file}")
        
        return result
    else:
        print(f"No data found for {mutation_type}, dose {dose}")
        return pd.DataFrame()

def analyze_genomic_distribution(merged_df, output_dir, mutation_type, dose):
    """
    Analyze the genomic distribution of mutation patterns and generate visualizations.
    
    Args:
        merged_df: DataFrame with merged annotation and pattern data
        output_dir: Directory to save output files
        mutation_type: Type of mutation being analyzed
        dose: Dose level being analyzed
    """
    if merged_df.empty:
        print(f"No data to analyze for {mutation_type}, dose {dose}")
        return
        
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze specific radiation patterns
    analyze_specific_patterns(merged_df, output_dir, mutation_type, dose)
    
    # 1. Distribution of pattern groups across gene locations
    if 'Gene_Location' in merged_df.columns:
        plt.figure(figsize=(12, 8))
        
        # Count by gene location and pattern group
        location_counts = pd.crosstab(
            merged_df['Gene_Location'], 
            merged_df['Pattern_Group'],
            normalize='columns'
        ) * 100  # Convert to percentage
        
        # Plot
        location_counts.plot(kind='bar', stacked=False)
        plt.title(f'Distribution of Mutation Patterns by Gene Location\n{mutation_type} - Dose {dose}')
        plt.ylabel('Percentage')
        plt.xlabel('Gene Location')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_gene_location_distribution.png", dpi=300)
        plt.close()
        
        # Save data
        location_counts.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_gene_location_distribution.csv")
    
    # 2. Distribution of pattern groups across feature types
    if 'Feature_Type' in merged_df.columns:
        # Get top 10 most common feature types
        top_features = merged_df['Feature_Type'].value_counts().nlargest(10).index.tolist()
        filtered_df = merged_df[merged_df['Feature_Type'].isin(top_features)]
        
        plt.figure(figsize=(14, 8))
        
        # Count by feature type and pattern group
        feature_counts = pd.crosstab(
            filtered_df['Feature_Type'], 
            filtered_df['Pattern_Group'],
            normalize='columns'
        ) * 100  # Convert to percentage
        
        # Plot
        feature_counts.plot(kind='bar', stacked=False)
        plt.title(f'Distribution of Mutation Patterns by Feature Type\n{mutation_type} - Dose {dose}')
        plt.ylabel('Percentage')
        plt.xlabel('Feature Type')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_feature_type_distribution.png", dpi=300)
        plt.close()
        
        # Save data
        feature_counts.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_feature_type_distribution.csv")
    
    # 3. Strand bias analysis (if relevant columns exist)
    if 'Gene_Strand' in merged_df.columns and 'Mutation_Strand' in merged_df.columns:
        # Define transcription strand (same or opposite)
        merged_df['Transcription_Strand'] = 'Non-genic'
        
        # For genic regions
        genic_mask = merged_df['Gene_Location'] != 'Intergenic'
        merged_df.loc[genic_mask & (merged_df['Gene_Strand'] == merged_df['Mutation_Strand']), 'Transcription_Strand'] = 'Same'
        merged_df.loc[genic_mask & (merged_df['Gene_Strand'] != merged_df['Mutation_Strand']), 'Transcription_Strand'] = 'Opposite'
        
        plt.figure(figsize=(10, 6))
        
        # Count by transcription strand and pattern group
        strand_counts = pd.crosstab(
            merged_df['Transcription_Strand'], 
            merged_df['Pattern_Group'],
            normalize='columns'
        ) * 100  # Convert to percentage
        
        # Plot
        strand_counts.plot(kind='bar')
        plt.title(f'Strand Bias Analysis by Mutation Pattern\n{mutation_type} - Dose {dose}')
        plt.ylabel('Percentage')
        plt.xlabel('Transcription Strand Relation')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_strand_bias.png", dpi=300)
        plt.close()
        
        # Save data
        strand_counts.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_strand_bias.csv")
    
    # 4. Most affected genes with pattern-specific breakdown
    if 'Gene_Name' in merged_df.columns:
        # Remove intergenic regions
        genic_df = merged_df[merged_df['Gene_Location'] != 'Intergenic']
        
        # Count mutations per gene and pattern group
        gene_pattern_counts = pd.crosstab(genic_df['Gene_Name'], genic_df['Pattern_Group'])
        
        # Get top 20 most affected genes
        gene_total_counts = gene_pattern_counts.sum(axis=1).sort_values(ascending=False)
        top_genes = gene_total_counts.head(20).index.tolist()
        
        # Filter for top genes
        top_gene_counts = gene_pattern_counts.loc[top_genes]
        
        # Calculate percentage of each pattern within each gene
        gene_percentages = top_gene_counts.div(top_gene_counts.sum(axis=1), axis=0) * 100
        
        # Create two visualizations: 
        # 1. Raw counts
        plt.figure(figsize=(12, 10))
        ax = top_gene_counts.plot(kind='barh', stacked=True)
        plt.title(f'Top 20 Genes with Most Mutations\n{mutation_type} - Dose {dose}')
        plt.xlabel('Number of Mutations')
        plt.ylabel('Gene')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add total count labels
        for i, gene in enumerate(top_genes):
            total = gene_total_counts[gene]
            ax.text(total + (total * 0.01), i, f' {total:,}', va='center')
            
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_top_genes_counts.png", dpi=300)
        plt.close()
        
        # 2. Percentage breakdown
        plt.figure(figsize=(12, 10))
        gene_percentages.plot(kind='barh', stacked=True)
        plt.title(f'Pattern Distribution in Top 20 Genes (Percentage)\n{mutation_type} - Dose {dose}')
        plt.xlabel('Percentage of Mutations')
        plt.ylabel('Gene')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_top_genes_percentage.png", dpi=300)
        plt.close()
        
        # Save data
        top_gene_counts.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_top_genes_counts.csv")
        gene_percentages.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_top_genes_percentage.csv")
        
        # 3. Separate analysis for radiation-specific vs. control-specific vs. baseline patterns in top genes
        radiation_patterns = [col for col in gene_pattern_counts.columns if 'Radiation' in col]
        control_patterns = [col for col in gene_pattern_counts.columns if 'Control' in col]
        baseline_patterns = [col for col in gene_pattern_counts.columns if 'Baseline' in col]
        
        # Calculate sums for each category
        if radiation_patterns and control_patterns and baseline_patterns:
            gene_pattern_sums = pd.DataFrame({
                'Radiation': gene_pattern_counts[radiation_patterns].sum(axis=1),
                'Control': gene_pattern_counts[control_patterns].sum(axis=1),
                'Baseline': gene_pattern_counts[baseline_patterns].sum(axis=1)
            })
            
            # Get top genes for each category
            top_radiation = gene_pattern_sums['Radiation'].nlargest(10).index.tolist()
            top_control = gene_pattern_sums['Control'].nlargest(10).index.tolist()
            top_baseline = gene_pattern_sums['Baseline'].nlargest(10).index.tolist()
            
            # Plot top genes for each category
            plt.figure(figsize=(12, 8))
            gene_pattern_sums.loc[top_radiation]['Radiation'].sort_values(ascending=True).plot(kind='barh')
            plt.title(f'Top 10 Genes Most Affected by Radiation-Specific Patterns\n{mutation_type} - Dose {dose}')
            plt.xlabel('Number of Radiation-Specific Mutations')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_top_radiation_genes.png", dpi=300)
            plt.close()
            
            plt.figure(figsize=(12, 8))
            gene_pattern_sums.loc[top_control]['Control'].sort_values(ascending=True).plot(kind='barh')
            plt.title(f'Top 10 Genes Most Affected by Control-Specific Patterns\n{mutation_type} - Dose {dose}')
            plt.xlabel('Number of Control-Specific Mutations')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_top_control_genes.png", dpi=300)
            plt.close()
            
            plt.figure(figsize=(12, 8))
            gene_pattern_sums.loc[top_baseline]['Baseline'].sort_values(ascending=True).plot(kind='barh')
            plt.title(f'Top 10 Genes Most Affected by Baseline Patterns\n{mutation_type} - Dose {dose}')
            plt.xlabel('Number of Baseline Mutations')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_top_baseline_genes.png", dpi=300)
            plt.close()
            
            # Save the data
            gene_pattern_sums.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_category_sums.csv")
    
    # 5. Heatmap of individual patterns distribution within specific genomic contexts
    if 'Gene_Location' in merged_df.columns and 'Pattern' in merged_df.columns:
        plt.figure(figsize=(12, 10))
        
        # Get the specific radiation patterns
        specific_radiation_patterns = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']
        pattern_df = merged_df[merged_df['Pattern'].isin(specific_radiation_patterns)]
        
        if not pattern_df.empty:
            # Create the crosstab
            pattern_location = pd.crosstab(
                pattern_df['Pattern'], 
                pattern_df['Gene_Location'],
                normalize='index'  # Normalize by pattern
            ) * 100  # Convert to percentage
            
            # Plot heatmap
            sns.heatmap(pattern_location, annot=True, fmt='.1f', cmap='viridis', cbar_kws={'label': 'Percentage'})
            plt.title(f'Distribution of Radiation Patterns Across Genomic Regions\n{mutation_type} - Dose {dose}')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_location_heatmap.png", dpi=300)
            plt.close()
            
            # Save data
            pattern_location.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_location_heatmap.csv")
    
    # 6. For InDels, analyze additional annotations if present
    if mutation_type == 'ID':
        # Check for InDel-specific columns
        indel_columns = [col for col in merged_df.columns if col.startswith('Indel_')]
        
        if indel_columns:
            print(f"Found InDel-specific columns: {indel_columns}")
            
            # Create a directory for InDel-specific analyses
            indel_dir = os.path.join(output_dir, 'indel_analysis')
            os.makedirs(indel_dir, exist_ok=True)
            
            # Analyze each InDel property
            for col in indel_columns:
                if col in ['Indel_Type', 'Indel_Mechanism']:
                    # Categorical columns
                    if merged_df[col].nunique() <= 15:  # Only if not too many categories
                        # Distribution by pattern group
                        plt.figure(figsize=(12, 8))
                        property_pattern = pd.crosstab(
                            merged_df[col],
                            merged_df['Pattern_Group'],
                            normalize='index'
                        ) * 100
                        
                        property_pattern.plot(kind='bar', stacked=False)
                        plt.title(f'{col} Distribution by Pattern Group\n{mutation_type} - Dose {dose}')
                        plt.ylabel('Percentage')
                        plt.xlabel(col)
                        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
                        plt.tight_layout()
                        plt.savefig(f"{indel_dir}/{mutation_type}_dose_{dose}_{col}_by_pattern.png", dpi=300)
                        plt.close()
                        
                        # Save data
                        property_pattern.to_csv(f"{indel_dir}/{mutation_type}_dose_{dose}_{col}_by_pattern.csv")
                
                elif col in ['Indel_Size', 'Repeat_Length']:
                    # Numeric columns
                    # Create size categories for better visualization
                    if col == 'Indel_Size':
                        bins = [-100, -10, -5, -2, -1, 1, 2, 5, 10, 100]
                        labels = ['<-10', '-10 to -5', '-5 to -2', '-2 to -1', '-1 to 1', '1 to 2', '2 to 5', '5 to 10', '>10']
                    else:  # Repeat_Length
                        bins = [0, 1, 2, 3, 4, 5, 10, 20, 100]
                        labels = ['0', '1', '2', '3', '4', '5-9', '10-19', '20+']
                    
                    merged_df[f'{col}_Category'] = pd.cut(merged_df[col], bins=bins, labels=labels)
                    
                    # Distribution by pattern group
                    plt.figure(figsize=(12, 8))
                    property_pattern = pd.crosstab(
                        merged_df[f'{col}_Category'],
                        merged_df['Pattern_Group'],
                        normalize='columns'
                    ) * 100
                    
                    property_pattern.plot(kind='bar', stacked=False)
                    plt.title(f'{col} Distribution by Pattern Group\n{mutation_type} - Dose {dose}')
                    plt.ylabel('Percentage')
                    plt.xlabel(col)
                    plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.tight_layout()
                    plt.savefig(f"{indel_dir}/{mutation_type}_dose_{dose}_{col}_by_pattern.png", dpi=300)
                    plt.close()
                    
                    # Save data
                    property_pattern.to_csv(f"{indel_dir}/{mutation_type}_dose_{dose}_{col}_by_pattern.csv")
    
    # 7. Summary statistics
    summary = {
        'Total_Mutations': len(merged_df),
        'Pattern_Group_Counts': merged_df['Pattern_Group'].value_counts().to_dict(),
        'Pattern_Group_Percentages': (merged_df['Pattern_Group'].value_counts(normalize=True) * 100).to_dict(),
        'Top_Patterns': merged_df['Pattern'].value_counts().head(10).to_dict()
    }
    
    # Save summary as CSV
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_summary.csv", index=False)
    
    print(f"Analysis completed for {mutation_type}, dose {dose}")


def compare_doses(annotated_dir, pattern_dir, mutation_type, output_dir, doses=None):
    """
    Compare mutation patterns across different doses for a specific mutation type.
    
    Args:
        annotated_dir: Directory containing annotated CSV files
        pattern_dir: Directory containing pattern CSV files
        mutation_type: One of 'SNV', 'DBS', 'MNS', 'ID'
        output_dir: Directory to save output files
        doses: List of doses to compare (default: A-E)
    """
    if doses is None:
        doses = ['A', 'B', 'C', 'D', 'E']
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to the merged data directory - use the project-level merged data directory
    merged_data_dir = os.path.join("merged_data", mutation_type)
    
    # Check if combined file already exists
    combined_file = os.path.join(output_dir, f"{mutation_type}_all_doses_merged.csv")
    if os.path.exists(combined_file):
        print(f"Loading existing combined dose data from {combined_file}")
        try:
            # Instead of loading the entire file, we'll analyze it in chunks
            print("Processing existing combined file in chunks...")
            process_combined_data_in_chunks(combined_file, output_dir, mutation_type, doses)
            return
        except Exception as e:
            print(f"Error processing existing combined file: {e}. Will regenerate.")
    
    # If we don't have a combined file or couldn't process it,
    # we'll process each dose individually and combine the results
    
    # 1. First, compute basic statistics for each dose separately
    dose_stats = {}
    for dose in doses:
        merged_file = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        
        if not os.path.exists(merged_file):
            print(f"Warning: No merged file found for dose {dose}. Skipping.")
            continue
        
        print(f"Processing dose {dose} statistics...")
        try:
            # Process the file in chunks to collect statistics
            dose_stats[dose] = compute_dose_statistics(merged_file, dose)
        except Exception as e:
            print(f"Error processing dose {dose}: {e}. Skipping.")
    
    if not dose_stats:
        print(f"No data found for {mutation_type} across doses")
        return
    
    # 2. Generate pattern group distribution plot
    try:
        generate_pattern_group_plot(dose_stats, output_dir, mutation_type)
    except Exception as e:
        print(f"Error generating pattern group distribution plot: {e}")
    
    # 3. Process each dose individually for genomic context analysis
    for dose in doses:
        merged_file = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        
        if not os.path.exists(merged_file):
            continue
        
        print(f"Analyzing genomic contexts for dose {dose}...")
        try:
            analyze_genomic_contexts(merged_file, dose, doses, output_dir, mutation_type)
        except Exception as e:
            print(f"Error analyzing genomic contexts for dose {dose}: {e}")
        
        # Force garbage collection after each dose
        gc.collect()
    
    # 4. Gene-specific analysis
    try:
        analyze_gene_specific_patterns(merged_data_dir, doses, output_dir, mutation_type)
    except Exception as e:
        print(f"Error in gene-specific analysis: {e}")
    
    # 5. For InDels, do specialized analysis
    if mutation_type == 'ID':
        try:
            analyze_indel_properties(merged_data_dir, doses, output_dir, mutation_type)
        except Exception as e:
            print(f"Error in InDel-specific analysis: {e}")
    
    print(f"Dose comparison completed for {mutation_type}")

def compute_dose_statistics(file_path, dose):
    """
    Compute basic statistics for a dose by processing the file in chunks.
    """
    pattern_group_counts = defaultdict(int)
    total_mutations = 0
    
    # Process in chunks
    for chunk in pd.read_csv(file_path, usecols=['Pattern', 'Pattern_Group'], chunksize=50000):
        # Count pattern groups
        chunk_counts = chunk['Pattern_Group'].value_counts()
        for group, count in chunk_counts.items():
            pattern_group_counts[group] += count
        
        # Add to total mutations
        total_mutations += len(chunk)
        
        # Free memory
        del chunk
        gc.collect()
    
    return {
        'total': total_mutations,
        'pattern_groups': pattern_group_counts
    }

def generate_pattern_group_plot(dose_stats, output_dir, mutation_type):
    """
    Generate pattern group distribution plot from dose statistics.
    """
    # Create DataFrame for plotting
    data = []
    for dose, stats in dose_stats.items():
        total = stats['total']
        for group, count in stats['pattern_groups'].items():
            percentage = (count / total) * 100
            data.append({'Dose': dose, 'Pattern_Group': group, 'Percentage': percentage})
    
    if not data:
        print("No data available for pattern group distribution plot")
        return
    
    df = pd.DataFrame(data)
    
    # Pivot for plotting
    pivot_df = df.pivot(index='Dose', columns='Pattern_Group', values='Percentage')
    
    # Plot
    plt.figure(figsize=(10, 6))
    pivot_df.plot(kind='bar', stacked=True)
    plt.title(f'Distribution of Mutation Pattern Groups Across Doses\n{mutation_type}')
    plt.ylabel('Percentage')
    plt.xlabel('Dose')
    plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{mutation_type}_dose_comparison_pattern_groups.png", dpi=300)
    plt.close()
    
    # Save data
    pivot_df.to_csv(f"{output_dir}/{mutation_type}_dose_comparison_pattern_groups.csv")
    
    # Free memory
    del df
    del pivot_df
    gc.collect()

def analyze_genomic_contexts(file_path, dose, all_doses, output_dir, mutation_type):
    """
    Analyze genomic contexts for a specific dose.
    """
    # We'll just check if the file has the needed column
    sample = pd.read_csv(file_path, nrows=1)
    if 'Gene_Location' not in sample.columns:
        print(f"Gene_Location column not found in {file_path}, skipping genomic context analysis")
        return
    
    # Process the file in chunks
    location_counts = defaultdict(lambda: defaultdict(int))
    pattern_groups = set()
    
    for chunk in pd.read_csv(file_path, usecols=['Gene_Location', 'Pattern_Group'], chunksize=50000):
        # Count by gene location and pattern group
        for _, row in chunk.iterrows():
            location = row['Gene_Location']
            pattern_group = row['Pattern_Group']
            location_counts[pattern_group][location] += 1
            pattern_groups.add(pattern_group)
        
        # Free memory
        del chunk
        gc.collect()
    
    # Process each pattern group
    for pattern_group in pattern_groups:
        # Skip if too few data points
        total = sum(location_counts[pattern_group].values())
        if total < 10:
            print(f"Skipping {pattern_group} - too few data points ({total})")
            continue
        
        # Convert counts to percentages
        percentages = {}
        for location, count in location_counts[pattern_group].items():
            percentages[location] = (count / total) * 100
        
        # Save data for this dose and pattern group
        data = {'Dose': dose, 'Gene_Location': [], 'Percentage': []}
        for location, percentage in percentages.items():
            data['Gene_Location'].append(location)
            data['Percentage'].append(percentage)
        
        # Save to file for later aggregation
        safe_pattern_name = pattern_group.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '')
        dose_file = f"{output_dir}/{mutation_type}_{safe_pattern_name}_dose_{dose}_genomic.csv"
        pd.DataFrame(data).to_csv(dose_file, index=False)
    
    # Free memory
    del location_counts
    gc.collect()

def analyze_gene_specific_patterns(merged_data_dir, doses, output_dir, mutation_type):
    """
    Analyze gene-specific patterns across doses.
    """
    # First check if any dose has gene information
    has_gene_info = False
    for dose in doses:
        file_path = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        if not os.path.exists(file_path):
            continue
        
        sample = pd.read_csv(file_path, nrows=1)
        if 'Gene_Name' in sample.columns:
            has_gene_info = True
            break
    
    if not has_gene_info:
        print(f"No Gene_Name information found for {mutation_type}, skipping gene-specific analysis")
        return
    
    # Count top genes for each dose
    top_genes_by_dose = {}
    for dose in doses:
        file_path = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        if not os.path.exists(file_path):
            continue
        
        gene_counts = defaultdict(int)
        total_processed = 0
        
        # Read in chunks
        for chunk in pd.read_csv(file_path, usecols=['Gene_Name', 'Gene_Location'], chunksize=50000):
            # Skip intergenic regions
            genic_chunk = chunk[chunk['Gene_Location'] != 'Intergenic']
            
            # Count genes
            chunk_counts = genic_chunk['Gene_Name'].value_counts()
            for gene, count in chunk_counts.items():
                if pd.notna(gene):  # Skip NaN gene names
                    gene_counts[gene] += count
            
            total_processed += len(chunk)
            
            # Free memory
            del chunk
            gc.collect()
        
        # Get top genes
        gene_count_series = pd.Series(gene_counts)
        top_genes = gene_count_series.nlargest(20).index.tolist()
        top_genes_by_dose[dose] = top_genes
        
        # Save top genes for this dose
        top_gene_data = {'Gene': top_genes, 'Count': [gene_counts[gene] for gene in top_genes]}
        pd.DataFrame(top_gene_data).to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_top_genes.csv", index=False)
    
    # Find genes that appear in the top 20 across multiple doses
    all_top_genes = []
    for genes in top_genes_by_dose.values():
        all_top_genes.extend(genes)
    
    gene_freq = pd.Series(all_top_genes).value_counts()
    multi_dose_genes = gene_freq[gene_freq > 1].index.tolist()
    
    if multi_dose_genes:
        print(f"Found {len(multi_dose_genes)} genes that appear in top 20 across multiple doses")
        
        # For each multi-dose gene, create a separate file with its data
        for gene in multi_dose_genes[:10]:  # Limit to top 10 genes
            gene_data = {'Dose': [], 'Count': []}
            
            for dose in doses:
                file_path = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
                if not os.path.exists(file_path):
                    continue
                
                count = 0
                for chunk in pd.read_csv(file_path, usecols=['Gene_Name'], chunksize=50000):
                    count += len(chunk[chunk['Gene_Name'] == gene])
                    
                    # Free memory
                    del chunk
                    gc.collect()
                
                gene_data['Dose'].append(dose)
                gene_data['Count'].append(count)
            
            # Save gene data
            pd.DataFrame(gene_data).to_csv(f"{output_dir}/{mutation_type}_{gene}_dose_counts.csv", index=False)

def analyze_indel_properties(merged_data_dir, doses, output_dir, mutation_type):
    """
    Analyze InDel-specific properties across doses.
    """
    # Check for InDel properties
    indel_columns = []
    for dose in doses:
        file_path = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
        if not os.path.exists(file_path):
            continue
        
        sample = pd.read_csv(file_path, nrows=1)
        indel_columns = [col for col in sample.columns if col.startswith('Indel_')]
        if indel_columns:
            break
    
    if not indel_columns:
        print(f"No InDel-specific columns found for {mutation_type}")
        return
    
    # For each InDel property, analyze distribution across doses
    for column in indel_columns:
        if column in ['Indel_Type', 'Indel_Mechanism']:
            # Categorical columns - count frequencies
            value_counts = {dose: defaultdict(int) for dose in doses}
            dose_totals = defaultdict(int)
            
            for dose in doses:
                file_path = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
                if not os.path.exists(file_path):
                    continue
                
                # Read in chunks
                for chunk in pd.read_csv(file_path, usecols=[column], chunksize=50000):
                    # Count values
                    chunk_counts = chunk[column].value_counts()
                    for value, count in chunk_counts.items():
                        if pd.notna(value):  # Skip NaN values
                            value_counts[dose][value] += count
                            dose_totals[dose] += count
                    
                    # Free memory
                    del chunk
                    gc.collect()
            
            # Convert to percentages and save
            data = []
            for dose in doses:
                if dose_totals[dose] > 0:
                    for value, count in value_counts[dose].items():
                        percentage = (count / dose_totals[dose]) * 100
                        data.append({'Dose': dose, column: value, 'Percentage': percentage})
            
            if data:
                df = pd.DataFrame(data)
                df.to_csv(f"{output_dir}/{mutation_type}_{column}_by_dose.csv", index=False)
            
            # Free memory
            del value_counts
            del dose_totals
            gc.collect()

def process_combined_data_in_chunks(combined_file, output_dir, mutation_type, doses):
    """
    Process a combined data file in chunks to generate comparison plots.
    """
    # 1. First pass: count pattern groups by dose
    pattern_group_by_dose = defaultdict(lambda: defaultdict(int))
    dose_totals = defaultdict(int)
    
    print("Pass 1: Counting pattern groups by dose...")
    for chunk in pd.read_csv(combined_file, usecols=['Dose', 'Pattern_Group'], chunksize=50000):
        # Count by dose and pattern group
        for _, row in chunk.iterrows():
            dose = row['Dose']
            pattern_group = row['Pattern_Group']
            pattern_group_by_dose[dose][pattern_group] += 1
            dose_totals[dose] += 1
        
        # Free memory
        del chunk
        gc.collect()
    
    # Generate pattern group distribution plot
    data = []
    for dose, group_counts in pattern_group_by_dose.items():
        for group, count in group_counts.items():
            percentage = (count / dose_totals[dose]) * 100
            data.append({'Dose': dose, 'Pattern_Group': group, 'Percentage': percentage})
    
    if data:
        df = pd.DataFrame(data)
        pivot_df = df.pivot(index='Dose', columns='Pattern_Group', values='Percentage')
        
        # Plot
        plt.figure(figsize=(10, 6))
        pivot_df.plot(kind='bar', stacked=True)
        plt.title(f'Distribution of Mutation Pattern Groups Across Doses\n{mutation_type}')
        plt.ylabel('Percentage')
        plt.xlabel('Dose')
        plt.legend(title='Pattern Group', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_comparison_pattern_groups.png", dpi=300)
        plt.close()
        
        # Save data
        pivot_df.to_csv(f"{output_dir}/{mutation_type}_dose_comparison_pattern_groups.csv")
    
    # Free memory
    del pattern_group_by_dose
    del dose_totals
    gc.collect()
    
    # 2. Second pass: check for genomic context data
    print("Pass 2: Checking for genomic context data...")
    sample = pd.read_csv(combined_file, nrows=1)
    if 'Gene_Location' in sample.columns:
        # For memory efficiency, we'll just generate distribution for radiation-specific patterns
        analyze_specific_pattern_genomics(combined_file, output_dir, mutation_type, doses)
    
    # 3. For gene-specific analysis, we'll handle in a separate pass if needed
    if 'Gene_Name' in sample.columns:
        print("Pass 3: Analyzing top genes across doses...")
        analyze_genes_in_chunks(combined_file, output_dir, mutation_type, doses)
    
    # 4. For InDels, handle in a separate pass if needed
    indel_columns = [col for col in sample.columns if col.startswith('Indel_')]
    if indel_columns and mutation_type == 'ID':
        print("Pass 4: Analyzing InDel properties...")
        analyze_indel_properties_in_chunks(combined_file, output_dir, mutation_type, doses, indel_columns)

def analyze_specific_pattern_genomics(combined_file, output_dir, mutation_type, doses):
    """
    Analyze specific pattern distributions across genomic contexts.
    """
    # For radiation patterns specifically
    radiation_pattern_counts = {dose: defaultdict(int) for dose in doses}
    rad_pattern_totals = defaultdict(int)
    
    for chunk in pd.read_csv(combined_file, 
                             usecols=['Dose', 'Pattern_Group', 'Gene_Location'], 
                             chunksize=50000):
        # Only keep radiation-specific patterns
        rad_chunk = chunk[chunk['Pattern_Group'] == 'Radiation-specific']
        
        # Count by dose and location
        for _, row in rad_chunk.iterrows():
            dose = row['Dose']
            location = row['Gene_Location']
            radiation_pattern_counts[dose][location] += 1
            rad_pattern_totals[dose] += 1
        
        # Free memory
        del chunk
        gc.collect()
    
    # Generate plot data
    data = []
    for dose, loc_counts in radiation_pattern_counts.items():
        if rad_pattern_totals[dose] > 0:
            for location, count in loc_counts.items():
                percentage = (count / rad_pattern_totals[dose]) * 100
                data.append({'Dose': dose, 'Gene_Location': location, 'Percentage': percentage})
    
    if data:
        df = pd.DataFrame(data)
        # Only keep top locations
        top_locations = df.groupby('Gene_Location')['Percentage'].sum().nlargest(8).index
        df_filtered = df[df['Gene_Location'].isin(top_locations)]
        
        # Pivot and plot
        pivot_df = df_filtered.pivot(index='Dose', columns='Gene_Location', values='Percentage')
        
        plt.figure(figsize=(12, 8))
        pivot_df.plot(kind='bar', stacked=True)
        plt.title(f'Genomic Distribution of Radiation Patterns Across Doses\n{mutation_type}')
        plt.ylabel('Percentage')
        plt.xlabel('Dose')
        plt.legend(title='Gene Location', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_radiation_patterns_location_by_dose.png", dpi=300)
        plt.close()
        
        # Save data
        pivot_df.to_csv(f"{output_dir}/{mutation_type}_radiation_patterns_location_by_dose.csv")
    
    # Free memory
    del radiation_pattern_counts
    del rad_pattern_totals
    gc.collect()

def analyze_genes_in_chunks(combined_file, output_dir, mutation_type, doses):
    """
    Analyze gene-specific patterns in chunks.
    """
    # Count genes by dose - only for genic regions
    gene_counts = {dose: defaultdict(int) for dose in doses}
    
    for chunk in pd.read_csv(combined_file, 
                            usecols=['Dose', 'Gene_Name', 'Gene_Location'], 
                            chunksize=50000):
        # Filter for genic regions
        genic_chunk = chunk[chunk['Gene_Location'] != 'Intergenic']
        
        # Count genes by dose
        for _, row in genic_chunk.iterrows():
            dose = row['Dose']
            gene = row['Gene_Name']
            if pd.notna(gene):  # Skip NaN gene names
                gene_counts[dose][gene] += 1
        
        # Free memory
        del chunk
        gc.collect()
    
    # Find top genes across all doses
    all_genes = {}
    for dose, counts in gene_counts.items():
        for gene, count in counts.items():
            if gene in all_genes:
                all_genes[gene] += count
            else:
                all_genes[gene] = count
    
    # Get top genes
    top_genes = sorted(all_genes.items(), key=lambda x: x[1], reverse=True)[:20]
    top_gene_names = [gene for gene, _ in top_genes]
    
    # Create gene count data
    data = []
    for dose in doses:
        for gene in top_gene_names:
            count = gene_counts[dose].get(gene, 0)
            if count > 0:  # Only include non-zero counts
                data.append({'Dose': dose, 'Gene': gene, 'Count': count})
    
    if data:
        df = pd.DataFrame(data)
        
        # Save data for future visualization
        df.to_csv(f"{output_dir}/{mutation_type}_top_genes_by_dose.csv", index=False)
        
        # Create a pivot table for easier visualization
        pivot_df = df.pivot(index='Gene', columns='Dose', values='Count')
        pivot_df = pivot_df.fillna(0)
        
        # Save the pivot table
        pivot_df.to_csv(f"{output_dir}/{mutation_type}_top_genes_dose_matrix.csv")
    
    # Free memory
    del gene_counts
    del all_genes
    gc.collect()

def analyze_indel_properties_in_chunks(combined_file, output_dir, mutation_type, doses, indel_columns):
    """
    Analyze InDel-specific properties in chunks.
    """
    # Process categorical columns
    for column in indel_columns:
        if column in ['Indel_Type', 'Indel_Mechanism']:
            # Count values by dose
            value_counts = {dose: defaultdict(int) for dose in doses}
            dose_totals = defaultdict(int)
            
            for chunk in pd.read_csv(combined_file, usecols=['Dose', column], chunksize=50000):
                for _, row in chunk.iterrows():
                    dose = row['Dose']
                    value = row[column]
                    if pd.notna(value):  # Skip NaN values
                        value_counts[dose][value] += 1
                        dose_totals[dose] += 1
                
                # Free memory
                del chunk
                gc.collect()
            
            # Generate plot data
            data = []
            for dose, counts in value_counts.items():
                if dose_totals[dose] > 0:
                    for value, count in counts.items():
                        percentage = (count / dose_totals[dose]) * 100
                        data.append({'Dose': dose, column: value, 'Percentage': percentage})
            
            if data:
                df = pd.DataFrame(data)
                
                # Limit to top categories if there are many
                if df[column].nunique() > 8:
                    top_values = df.groupby(column)['Percentage'].mean().nlargest(8).index
                    df = df[df[column].isin(top_values)]
                
                # Save data
                df.to_csv(f"{output_dir}/{mutation_type}_{column}_by_dose.csv", index=False)
            
            # Free memory
            del value_counts
            del dose_totals
            gc.collect()
def analyze_specific_patterns(merged_df, output_dir, mutation_type, dose):
    """
    Analyze specific radiation patterns (not just pattern groups) to understand
    the distribution of different temporal patterns.
    
    Args:
        merged_df: DataFrame with merged annotation and pattern data
        output_dir: Directory to save output files
        mutation_type: Type of mutation being analyzed
        dose: Dose level being analyzed
    """
    if merged_df.empty or 'Pattern' not in merged_df.columns:
        print(f"Cannot analyze specific patterns: missing required columns")
        return
    
    # Define radiation patterns
    radiation_patterns = ['0T00', '00T0', '000T', '0TT0', '00TT', '0T0T', '0TTT']
    
    # Filter for radiation patterns
    pattern_data = merged_df[merged_df['Pattern'].isin(radiation_patterns)]
    
    if pattern_data.empty:
        print(f"No specific radiation patterns found in the data")
        return
    
    # Count each specific pattern
    pattern_counts = pattern_data['Pattern'].value_counts().reindex(radiation_patterns, fill_value=0)
    total_patterns = pattern_counts.sum()
    pattern_percentages = (pattern_counts / total_patterns * 100).round(2)
    
    # Save the pattern distribution
    pattern_df = pd.DataFrame({
        'Pattern': pattern_counts.index,
        'Count': pattern_counts.values,
        'Percentage': pattern_percentages.values
    })
    pattern_df.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_specific_patterns.csv", index=False)
    
    # Plot the distribution
    plt.figure(figsize=(10, 6))
    ax = pattern_percentages.plot(kind='bar', color='steelblue')
    plt.title(f'Distribution of Specific Radiation Patterns\n{mutation_type} - Dose {dose}')
    plt.ylabel('Percentage of Radiation Patterns')
    plt.xlabel('Pattern')
    
    # Add percentage labels on top of bars
    for i, v in enumerate(pattern_percentages):
        ax.text(i, v + 0.5, f"{v}%", ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_specific_patterns.png", dpi=300)
    plt.close()
    
    # Group patterns by number of timepoints with mutations
    single_timepoint = ['0T00', '00T0', '000T']
    double_timepoint = ['0TT0', '00TT', '0T0T']
    triple_timepoint = ['0TTT']
    
    # Calculate percentage for each group
    single_pct = pattern_data['Pattern'].isin(single_timepoint).sum() / len(pattern_data) * 100
    double_pct = pattern_data['Pattern'].isin(double_timepoint).sum() / len(pattern_data) * 100
    triple_pct = pattern_data['Pattern'].isin(triple_timepoint).sum() / len(pattern_data) * 100
    
    # Plot the distribution by timepoint count
    plt.figure(figsize=(8, 6))
    timepoint_counts = pd.Series({
        'Single Timepoint': single_pct,
        'Double Timepoint': double_pct,
        'Triple Timepoint': triple_pct
    })
    ax = timepoint_counts.plot(kind='bar', color=['lightblue', 'steelblue', 'darkblue'])
    plt.title(f'Distribution of Radiation Patterns by Timepoints Affected\n{mutation_type} - Dose {dose}')
    plt.ylabel('Percentage of Radiation Patterns')
    
    # Add percentage labels on top of bars
    for i, v in enumerate(timepoint_counts):
        ax.text(i, v + 0.5, f"{v:.1f}%", ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_timepoint_distribution.png", dpi=300)
    plt.close()
    
    # If there's enough data, analyze pattern distribution by genomic context
    if 'Gene_Location' in merged_df.columns and len(pattern_data) >= 50:
        # Create a cross-tabulation of patterns by gene location
        location_pattern = pd.crosstab(
            pattern_data['Gene_Location'],
            pattern_data['Pattern'],
            normalize='columns'
        ) * 100
        
        # Plot heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(location_pattern, annot=True, fmt='.1f', cmap='viridis')
        plt.title(f'Distribution of Radiation Patterns by Genomic Location\n{mutation_type} - Dose {dose}')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_location_heatmap.png", dpi=300)
        plt.close()
        
        # Save data
        location_pattern.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_location_distribution.csv")
    
    # For InDels, analyze pattern distribution by InDel properties
    if mutation_type == 'ID' and 'Indel_Type' in merged_df.columns:
        indel_properties = ['Indel_Type', 'Indel_Mechanism', 'Indel_Size']
        
        for property_col in indel_properties:
            if property_col in merged_df.columns:
                if property_col == 'Indel_Size':
                    # Create categories for numeric Indel_Size
                    bins = [-100, -10, -5, -2, -1, 1, 2, 5, 10, 100]
                    labels = ['<-10', '-10 to -5', '-5 to -2', '-2 to -1', '-1 to 1', '1 to 2', '2 to 5', '5 to 10', '>10']
                    pattern_data['Size_Category'] = pd.cut(pattern_data['Indel_Size'], bins=bins, labels=labels)
                    property_values = pattern_data['Size_Category']
                else:
                    property_values = pattern_data[property_col]
                
                # Create a cross-tabulation
                property_pattern = pd.crosstab(
                    property_values,
                    pattern_data['Pattern'],
                    normalize='columns'
                ) * 100
                
                # Plot heatmap
                plt.figure(figsize=(12, 8))
                sns.heatmap(property_pattern, annot=True, fmt='.1f', cmap='viridis')
                plt.title(f'Distribution of Radiation Patterns by {property_col}\n{mutation_type} - Dose {dose}')
                plt.tight_layout()
                plt.savefig(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_{property_col}_heatmap.png", dpi=300)
                plt.close()
                
                # Save data
                property_pattern.to_csv(f"{output_dir}/{mutation_type}_dose_{dose}_pattern_{property_col}_distribution.csv")

    
def main():
    parser = argparse.ArgumentParser(description="Analyze mutation patterns across genomic contexts")
    parser.add_argument('-a', '--annotated-dir', required=True, help='Directory containing annotated CSV files')
    parser.add_argument('-p', '--pattern-dir', required=True, help='Directory containing pattern CSV files')
    parser.add_argument('-o', '--output-dir', default='pattern_analysis_results', help='Output directory')
    parser.add_argument('-m', '--mutation-types', nargs='+', default=['SNV', 'DBS', 'MNS', 'ID'], 
                        help='Mutation types to analyze')
    parser.add_argument('-d', '--doses', nargs='+', default=['A', 'B', 'C', 'D', 'E'], 
                        help='Doses to analyze')
    parser.add_argument('-c', '--chromosomes', nargs='+', help='Specific chromosomes to analyze (optional)')
    parser.add_argument('--compare-doses', action='store_true', help='Compare patterns across doses')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose debugging information')
    
    args = parser.parse_args()
    
    # Print directory information
    print(f"==== Directory Information ====")
    print(f"Annotated directory: {args.annotated_dir}")
    print(f"Pattern directory: {args.pattern_dir}")
    print(f"Output directory: {args.output_dir}")
    
    # Check if directories exist
    if not os.path.exists(args.annotated_dir):
        print(f"ERROR: Annotated directory does not exist: {args.annotated_dir}")
        return
    
    if not os.path.exists(args.pattern_dir):
        print(f"ERROR: Pattern directory does not exist: {args.pattern_dir}")
        return
    
    # List files in annotated directory for debugging (limited depth)
    print("\n==== Annotated Directory Structure ====")
    for root, dirs, files in os.walk(args.annotated_dir, topdown=True):
        level = root.replace(args.annotated_dir, '').count(os.sep)
        if level > 2:  # Limit directory traversal depth
            continue
        indent = ' ' * 4 * level
        print(f"{indent}{os.path.basename(root)}/")
        if level <= 1:  # Only show files up to a certain directory depth
            sub_indent = ' ' * 4 * (level + 1)
            for f in files[:5]:  # Show only first 5 files in each directory
                print(f"{sub_indent}{f}")
            if len(files) > 5:
                print(f"{sub_indent}... ({len(files) - 5} more files)")
    
    # List files in pattern directory for debugging (limited depth)
    print("\n==== Pattern Directory Structure ====")
    for root, dirs, files in os.walk(args.pattern_dir, topdown=True):
        level = root.replace(args.pattern_dir, '').count(os.sep)
        if level > 2:  # Limit directory traversal depth
            continue
        indent = ' ' * 4 * level
        print(f"{indent}{os.path.basename(root)}/")
        if level <= 1:  # Only show files up to a certain directory depth
            sub_indent = ' ' * 4 * (level + 1)
            for f in files[:5]:  # Show only first 5 files in each directory
                print(f"{sub_indent}{f}")
            if len(files) > 5:
                print(f"{sub_indent}... ({len(files) - 5} more files)")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create directory for merged data
    merged_data_dir = os.path.join(args.output_dir, "merged_data")
    os.makedirs(merged_data_dir, exist_ok=True)
    
    # Process each mutation type and dose combination
    for mutation_type in args.mutation_types:
        print(f"\n==== Processing {mutation_type} mutations ====")
        
        # Check for specific mutation type folders
        pattern_mt_path = os.path.join(args.pattern_dir, mutation_type)
        if not os.path.exists(pattern_mt_path):
            print(f"Warning: Pattern directory for {mutation_type} not found at {pattern_mt_path}")
            # Look for alternative paths
            possible_mt_dirs = [d for d in os.listdir(args.pattern_dir) 
                               if os.path.isdir(os.path.join(args.pattern_dir, d)) and mutation_type.lower() in d.lower()]
            if possible_mt_dirs:
                print(f"Found possible alternative directories: {possible_mt_dirs}")
        
        # Individual dose analysis
        for dose in args.doses:
            print(f"\n-- Processing dose {dose} --")
            dose_output_dir = os.path.join(args.output_dir, mutation_type, f"dose_{dose}")
            os.makedirs(dose_output_dir, exist_ok=True)
            
            # Find pattern files for this dose
            dose_pattern_dir = os.path.join(args.pattern_dir, mutation_type, f"dose_{dose}")
            if not os.path.exists(dose_pattern_dir):
                alt_dose_pattern_dir = os.path.join(args.pattern_dir, f"dose_{dose}")
                if os.path.exists(alt_dose_pattern_dir):
                    print(f"Using alternative dose pattern directory: {alt_dose_pattern_dir}")
                else:
                    print(f"Warning: Neither {dose_pattern_dir} nor {alt_dose_pattern_dir} exists")
            
            # Check if merged data already exists
            merged_file = os.path.join(merged_data_dir, f"{mutation_type}_dose_{dose}_merged.csv")
            if os.path.exists(merged_file):
                print(f"Found existing merged data: {merged_file}")
                try:
                    merged_data = pd.read_csv(merged_file)
                    print(f"Loaded {len(merged_data)} rows from existing merged data")
                except Exception as e:
                    print(f"Error loading existing merged data: {e}")
                    merged_data = load_and_merge_data(
                        args.annotated_dir, 
                        args.pattern_dir, 
                        mutation_type, 
                        dose, 
                        args.chromosomes
                    )
            else:
                # Load and analyze data
                merged_data = load_and_merge_data(
                    args.annotated_dir, 
                    args.pattern_dir, 
                    mutation_type, 
                    dose, 
                    args.chromosomes
                )
            
            if not merged_data.empty:
                analyze_genomic_distribution(merged_data, dose_output_dir, mutation_type, dose)
            else:
                print(f"No data available for {mutation_type}, dose {dose}")
        
        # Dose comparison (if requested)
        if args.compare_doses:
            comparison_dir = os.path.join(args.output_dir, mutation_type, "dose_comparison")
            compare_doses(args.annotated_dir, args.pattern_dir, mutation_type, comparison_dir, args.doses)
    
    print("\nAnalysis complete!")

if __name__ == '__main__':
    main()
