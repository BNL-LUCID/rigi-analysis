#!/usr/bin/env python3
"""Mutation Gene Pattern Catalog Generator.
========================================
Aggregates mutation data by Gene + Pattern + Dose + Mutation Type.
Uses dask for efficient processing of large mutation files.

Input:
- Directory with merged mutation CSVs (SNV_dose_*_merged.csv, ID_dose_*_merged.csv, etc.)

Output:
- mutation_gene_catalog.csv: Aggregated counts by Gene, Pattern, Dose, Mutation_Type
- mutation_gene_summary.csv: Summary per gene (total counts, pattern breakdown)

Expected input columns:
- Gene_Name (or Gene): Gene symbol
- Pattern: Temporal pattern (0T00, 00T0, 000T, etc.)
- Dose: Dose level (dA, dB, dC, dD, dE)
- Mutation_Type: SNV, DBS, ID, MNS
- Gene_Location: Exonic, Intronic, etc. (optional)
"""

import argparse
import warnings
from pathlib import Path

import dask.dataframe as dd
import pandas as pd
from dask.diagnostics import ProgressBar

warnings.filterwarnings('ignore')

# Valid patterns and doses
VALID_PATTERNS = ['0T00', '00T0', '000T', '0TT0', '0T0T', '00TT', '0TTT']
VALID_DOSES = ['dA', 'dB', 'dC', 'dD', 'dE']
MUTATION_TYPES = ['SNV', 'DBS', 'ID', 'MNS']


def find_mutation_files(mut_dir):
    """Find all mutation CSV files in directory."""
    mut_dir = Path(mut_dir)

    # Try different patterns
    patterns = [
        '*_merged.csv',
        'SNV_*.csv',
        'ID_*.csv',
        'DBS_*.csv',
        'MNS_*.csv',
        '*mutation*.csv'
    ]

    all_files = []
    for pattern in patterns:
        files = list(mut_dir.glob(pattern))
        all_files.extend(files)

    # Remove duplicates
    all_files = list(set(all_files))

    return sorted(all_files)


def identify_columns(ddf):
    """Identify key columns in the dataframe."""
    cols = ddf.columns.tolist()

    # Find gene column
    gene_col = None
    for possible in ['Gene_Name', 'Gene', 'gene_name', 'GENE', 'gene']:
        if possible in cols:
            gene_col = possible
            break

    # Find pattern column
    pattern_col = None
    for possible in ['Pattern', 'pattern', 'Temporal_Pattern']:
        if possible in cols:
            pattern_col = possible
            break

    # Find dose column
    dose_col = None
    for possible in ['Dose', 'dose', 'DOSE']:
        if possible in cols:
            dose_col = possible
            break

    # Find mutation type column
    mut_type_col = None
    for possible in ['Mutation_Type', 'mutation_type', 'MutationType', 'Type']:
        if possible in cols:
            mut_type_col = possible
            break

    # Find location column
    location_col = None
    for possible in ['Gene_Location', 'Location', 'Feature_Type', 'Consequence']:
        if possible in cols:
            location_col = possible
            break

    return {
        'gene': gene_col,
        'pattern': pattern_col,
        'dose': dose_col,
        'mut_type': mut_type_col,
        'location': location_col
    }


def process_mutation_file(filepath, file_index, total_files):
    """Process a single mutation file and return aggregated counts."""
    filepath = Path(filepath)
    print(f"\n[{file_index}/{total_files}] Processing: {filepath.name}")

    try:
        # Read with dask
        ddf = dd.read_csv(filepath, assume_missing=True, low_memory=False)

        # Identify columns
        col_map = identify_columns(ddf)

        if col_map['gene'] is None:
            print("  ⚠ No gene column found, skipping")
            return None

        print(f"  Columns: gene={col_map['gene']}, pattern={col_map['pattern']}, "
              f"dose={col_map['dose']}, type={col_map['mut_type']}")

        # Select and rename columns
        select_cols = [col_map['gene']]
        rename_map = {col_map['gene']: 'Gene'}

        if col_map['pattern']:
            select_cols.append(col_map['pattern'])
            rename_map[col_map['pattern']] = 'Pattern'

        if col_map['dose']:
            select_cols.append(col_map['dose'])
            rename_map[col_map['dose']] = 'Dose'

        if col_map['mut_type']:
            select_cols.append(col_map['mut_type'])
            rename_map[col_map['mut_type']] = 'Mutation_Type'

        if col_map['location']:
            select_cols.append(col_map['location'])
            rename_map[col_map['location']] = 'Location'

        # Filter to existing columns
        select_cols = [c for c in select_cols if c in ddf.columns]

        ddf_subset = ddf[select_cols].rename(columns=rename_map)

        # Infer mutation type from filename if not in data
        if 'Mutation_Type' not in ddf_subset.columns:
            mut_type = 'SNV'  # default
            fname_upper = filepath.name.upper()
            for mt in MUTATION_TYPES:
                if mt in fname_upper:
                    mut_type = mt
                    break
            ddf_subset['Mutation_Type'] = mut_type
            print(f"  Inferred Mutation_Type from filename: {mut_type}")

        # Infer dose from filename if not in data
        if 'Dose' not in ddf_subset.columns:
            dose = 'dA'  # default
            fname_lower = filepath.name.lower()
            for d in VALID_DOSES:
                if d.lower() in fname_lower or f"dose_{d[-1]}" in fname_lower:
                    dose = d
                    break
            ddf_subset['Dose'] = dose
            print(f"  Inferred Dose from filename: {dose}")

        # Default pattern if missing
        if 'Pattern' not in ddf_subset.columns:
            ddf_subset['Pattern'] = 'Unknown'

        # Filter out invalid genes
        ddf_subset = ddf_subset[ddf_subset['Gene'].notnull()]
        ddf_subset = ddf_subset[ddf_subset['Gene'] != '']
        ddf_subset = ddf_subset[ddf_subset['Gene'] != 'Unknown']

        # Group and count
        group_cols = ['Gene', 'Dose', 'Pattern', 'Mutation_Type']
        if 'Location' in ddf_subset.columns:
            group_cols.append('Location')

        print(f"  Aggregating by: {group_cols}")

        with ProgressBar():
            counts = ddf_subset.groupby(group_cols).size().reset_index()
            counts.columns = group_cols + ['Count']
            result = counts.compute()

        print(f"  ✓ Generated {len(result):,} aggregated records")

        return result

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None


def aggregate_all_mutations(mut_dir, output_dir):
    """Process all mutation files and create aggregated catalog."""
    mut_dir = Path(mut_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MUTATION GENE PATTERN CATALOG GENERATOR")
    print("=" * 70)

    # Find files
    mutation_files = find_mutation_files(mut_dir)
    print(f"\nFound {len(mutation_files)} mutation files in {mut_dir}")

    if not mutation_files:
        print("No mutation files found!")
        return None

    for f in mutation_files[:10]:
        print(f"  - {f.name}")
    if len(mutation_files) > 10:
        print(f"  ... and {len(mutation_files) - 10} more")

    # Process each file
    all_results = []

    for i, filepath in enumerate(mutation_files, 1):
        result = process_mutation_file(filepath, i, len(mutation_files))
        if result is not None and len(result) > 0:
            all_results.append(result)

    if not all_results:
        print("\nNo results generated!")
        return None

    # Combine all results
    print("\n" + "-" * 70)
    print("Combining results...")

    combined = pd.concat(all_results, ignore_index=True)
    print(f"  Total records before aggregation: {len(combined):,}")

    # Re-aggregate in case of overlapping files
    group_cols = ['Gene', 'Dose', 'Pattern', 'Mutation_Type']
    if 'Location' in combined.columns:
        group_cols.append('Location')

    catalog = combined.groupby(group_cols)['Count'].sum().reset_index()
    print(f"  Total records after final aggregation: {len(catalog):,}")

    # Save detailed catalog
    catalog_file = output_dir / 'mutation_gene_catalog.csv'
    catalog.to_csv(catalog_file, index=False)
    print(f"\n✓ Saved: {catalog_file}")

    # Create gene-level summary
    print("\nCreating gene-level summary...")

    gene_summary = create_gene_summary(catalog)
    summary_file = output_dir / 'mutation_gene_summary.csv'
    gene_summary.to_csv(summary_file, index=False)
    print(f"✓ Saved: {summary_file}")

    # Print statistics
    print_statistics(catalog, gene_summary)

    return catalog, gene_summary


def create_gene_summary(catalog):
    """Create gene-level summary with pattern and type breakdowns."""
    summary = catalog.groupby('Gene').agg({
        'Count': 'sum'
    }).reset_index()
    summary.columns = ['Gene', 'Total_Mutations']

    # Add pattern breakdown
    for pattern in VALID_PATTERNS:
        pattern_counts = catalog[catalog['Pattern'] == pattern].groupby('Gene')['Count'].sum()
        summary[f'Mut_{pattern}'] = summary['Gene'].map(pattern_counts).fillna(0).astype(int)

    # Add mutation type breakdown
    for mut_type in MUTATION_TYPES:
        type_counts = catalog[catalog['Mutation_Type'] == mut_type].groupby('Gene')['Count'].sum()
        summary[f'Mut_{mut_type}'] = summary['Gene'].map(type_counts).fillna(0).astype(int)

    # Add dose breakdown
    for dose in VALID_DOSES:
        dose_counts = catalog[catalog['Dose'] == dose].groupby('Gene')['Count'].sum()
        summary[f'Mut_{dose}'] = summary['Gene'].map(dose_counts).fillna(0).astype(int)

    # Add location breakdown if available
    # Note: Raw data uses 'Exon', 'Intron' - searching for these will also match 'Exonic', 'Intronic'
    if 'Location' in catalog.columns:
        for loc, col_name in [('Exon', 'Exonic'), ('Intron', 'Intronic'), ('Intergenic', 'Intergenic')]:
            loc_mask = catalog['Location'].str.contains(loc, case=False, na=False)
            loc_counts = catalog[loc_mask].groupby('Gene')['Count'].sum()
            summary[f'Mut_{col_name}'] = summary['Gene'].map(loc_counts).fillna(0).astype(int)

    # Calculate temporal categories
    transient_cols = [f'Mut_{p}' for p in ['0T00', '00T0', '000T'] if f'Mut_{p}' in summary.columns]
    persistent_cols = [f'Mut_{p}' for p in ['0TT0', '00TT', '0TTT'] if f'Mut_{p}' in summary.columns]

    if transient_cols:
        summary['Mut_Transient'] = summary[transient_cols].sum(axis=1)
    if persistent_cols:
        summary['Mut_Persistent'] = summary[persistent_cols].sum(axis=1)

    # Sort by total mutations
    summary = summary.sort_values('Total_Mutations', ascending=False)

    return summary


def print_statistics(catalog, gene_summary):
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    print(f"\nTotal mutations: {catalog['Count'].sum():,}")
    print(f"Unique genes: {catalog['Gene'].nunique():,}")
    print(f"Catalog records: {len(catalog):,}")

    # Pattern distribution
    print("\nPattern Distribution:")
    pattern_counts = catalog.groupby('Pattern')['Count'].sum().sort_values(ascending=False)
    total = pattern_counts.sum()
    for pattern, count in pattern_counts.items():
        pct = count / total * 100
        print(f"  {pattern}: {count:,} ({pct:.1f}%)")

    # Mutation type distribution
    print("\nMutation Type Distribution:")
    type_counts = catalog.groupby('Mutation_Type')['Count'].sum().sort_values(ascending=False)
    for mut_type, count in type_counts.items():
        pct = count / total * 100
        print(f"  {mut_type}: {count:,} ({pct:.1f}%)")

    # Dose distribution
    print("\nDose Distribution:")
    dose_counts = catalog.groupby('Dose')['Count'].sum().sort_values(ascending=False)
    for dose, count in dose_counts.items():
        pct = count / total * 100
        print(f"  {dose}: {count:,} ({pct:.1f}%)")

    # Top genes
    print("\nTop 20 Genes by Mutation Count:")
    print("-" * 50)
    top_genes = gene_summary.head(20)[['Gene', 'Total_Mutations']]
    for _, row in top_genes.iterrows():
        print(f"  {row['Gene']}: {row['Total_Mutations']:,}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate mutation gene pattern catalog from merged CSV files'
    )
    parser.add_argument(
        '--mut-dir', '-m',
        required=True,
        help='Directory containing merged mutation CSV files'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='mutation_catalog',
        help='Output directory for catalog files'
    )

    args = parser.parse_args()

    catalog, summary = aggregate_all_mutations(args.mut_dir, args.output_dir)

    if catalog is not None:
        print("\n" + "=" * 70)
        print("COMPLETE")
        print("=" * 70)
        print("\nOutput files:")
        print("  - mutation_gene_catalog.csv (detailed by gene/dose/pattern/type)")
        print("  - mutation_gene_summary.csv (gene-level summary)")
        print("\nUse mutation_gene_summary.csv as input to stage1_hotspot_detection.py")


if __name__ == "__main__":
    main()
