#!/usr/bin/env python3
"""Fast DBS Search Module.
=====================
Optimized parallel search for DBS near breakpoints.
Import this in your main analysis scripts.
"""

from collections import defaultdict
from multiprocessing import Pool, cpu_count

import numpy as np


def create_chromosome_index(dbs_df):
    """Create sorted chromosome index for binary search."""
    chrom_index = {}

    for chrom in dbs_df['Chromosome'].unique():
        chrom_dbs = dbs_df[dbs_df['Chromosome'] == chrom]['Start'].values
        chrom_dbs.sort()
        chrom_index[chrom] = chrom_dbs

    return chrom_index


def _count_worker(args):
    """Worker function for parallel processing."""
    bp_chunk, chrom_index, window = args

    gene_dbs = defaultdict(int)

    for _, bp in bp_chunk.iterrows():
        chrom = bp['Chrom']

        if chrom not in chrom_index:
            continue

        pos = bp['Pos']
        gene = bp['Gene']

        chrom_dbs = chrom_index[chrom]

        # Binary search (much faster than filtering)
        left = np.searchsorted(chrom_dbs, pos - window, side='left')
        right = np.searchsorted(chrom_dbs, pos + window, side='right')

        count = right - left

        if count > 0:
            gene_dbs[gene] += count

    return gene_dbs


def fast_dbs_count(breakpoints_df, dbs_df, window=10, n_cores=None, verbose=True):
    """Fast parallel DBS counting.

    Args:
        breakpoints_df: DataFrame with [Gene, Chrom, Pos]
        dbs_df: DataFrame with [Chromosome, Start]
        window: Search window (bp)
        n_cores: CPU cores (None = auto)
        verbose: Print progress

    Returns:
        dict: {gene: dbs_count}
    """
    if n_cores is None:
        n_cores = max(1, cpu_count() - 1)

    if verbose:
        print(f"  Fast DBS search: {len(breakpoints_df):,} breakpoints, {n_cores} cores")

    # Create index
    chrom_index = create_chromosome_index(dbs_df)

    # Split into chunks
    chunk_size = max(1000, len(breakpoints_df) // (n_cores * 4))
    chunks = [
        breakpoints_df.iloc[i:i+chunk_size]
        for i in range(0, len(breakpoints_df), chunk_size)
    ]

    # Parallel processing
    args_list = [(chunk, chrom_index, window) for chunk in chunks]

    with Pool(n_cores) as pool:
        results = pool.map(_count_worker, args_list)

    # Merge
    gene_dbs = defaultdict(int)
    for chunk_result in results:
        for gene, count in chunk_result.items():
            gene_dbs[gene] += count

    if verbose:
        n_with_dbs = sum(1 for c in gene_dbs.values() if c > 0)
        print(f"  ✓ {n_with_dbs:,} genes with DBS")

    return dict(gene_dbs)


def annotate_with_dbs(gene_df, breakpoints_df, dbs_df, window=10, n_cores=None):
    """Add DBS counts to gene DataFrame.

    Args:
        gene_df: DataFrame with 'Gene' column
        breakpoints_df: DataFrame with [Gene, Chrom, Pos]
        dbs_df: DataFrame with [Chromosome, Start]
        window: Search window (bp)
        n_cores: CPU cores

    Returns:
        DataFrame with added [DBS_Count, Has_DBS] columns
    """
    # Count DBS
    gene_dbs = fast_dbs_count(breakpoints_df, dbs_df, window, n_cores)

    # Add to DataFrame
    gene_df['DBS_Count'] = gene_df['Gene'].map(gene_dbs).fillna(0).astype(int)
    gene_df['Has_DBS'] = (gene_df['DBS_Count'] > 0).astype(int)

    return gene_df
