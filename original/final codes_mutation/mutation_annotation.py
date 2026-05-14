#!/usr/bin/env python
"""
Optimized mutation annotation using vectorized operations and efficient lookups.
"""

import os
import pandas as pd
import numpy as np
import logging
import sys
import time
import gc
import argparse
import pickle
import traceback
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    handlers=[
        logging.FileHandler('mutation_annotation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)


class OptimizedMutationAnnotator:
    """Optimized mutation annotator with vectorized operations."""
    
    def __init__(self, annotation_dir: str = "annotations", build: str = "hg38", output_dir: str = "results"):
        """
        Initialize annotator.
        
        Args:
            annotation_dir: Directory containing annotation files
            build: Genome build version
            output_dir: Directory for output files
        """
        self.annotation_dir = os.path.abspath(annotation_dir)
        self.build = build
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.interval_trees_file = os.path.join(self.annotation_dir, f"{build}_interval_trees.pkl")
        self.feature_data_file = os.path.join(self.annotation_dir, f"{build}_feature_data.pkl")
        
        if not os.path.exists(self.interval_trees_file) or not os.path.exists(self.feature_data_file):
            raise FileNotFoundError(
                f"Annotation data not found at:\n"
                f"  {self.interval_trees_file}\n"
                f"  {self.feature_data_file}\n"
                f"Please run annotation_preprocessing.py first."
            )
        
        self.interval_trees = None
        self.feature_data = None
        self.chromosome_cache = {}
        
        self.load_annotation_data()

    def load_annotation_data(self):
        """Load interval trees and feature data."""
        try:
            logger.info(f"Loading interval trees...")
            with open(self.interval_trees_file, 'rb') as f:
                self.interval_trees = pickle.load(f)
            
            logger.info(f"Loading feature data...")
            with open(self.feature_data_file, 'rb') as f:
                self.feature_data = pickle.load(f)
            
            chrom_count = len(self.interval_trees)
            logger.info(f"Loaded annotations for {chrom_count} chromosomes")
            
        except Exception as e:
            logger.error(f"Error loading annotation data: {str(e)}")
            raise

    def get_best_feature_for_position(self, chrom: str, pos: int) -> Optional[Tuple[str, Dict]]:
        """
        Efficiently get best feature for a position using priority system.
        
        Args:
            chrom: Chromosome name
            pos: 0-based position
            
        Returns:
            (feature_type, feature_data) or None
        """
        if chrom not in self.interval_trees:
            return None
        
        best_priority = -1
        best_feature = None
        best_feature_type = None
        
        # Priority: exon(5) > 5utr(4) > 3utr(3) > intron(2) > gene(1) > promoter(0)
        feature_priorities = {
            'exon': 5,
            '5utr': 4,
            '3utr': 3,
            'intron': 2,
            'gene': 1,
            'promoter': 0
        }
        
        # Query each feature type (fast with interval tree)
        for feature_type, tree in self.interval_trees[chrom].items():
            priority = feature_priorities.get(feature_type, -1)
            
            # Skip if lower priority than current best
            if priority <= best_priority:
                continue
            
            # Query tree for overlapping intervals
            overlaps = tree[pos:pos+1]
            
            if overlaps:
                # Get first overlap
                for interval in overlaps:
                    priority_val, feature_id = interval.data
                    
                    best_priority = priority
                    best_feature_id = feature_id
                    best_feature_type = feature_type
                    break
        
        # Retrieve feature data if found
        if best_feature_type and best_feature_id is not None:
            try:
                feature_data = self.feature_data[chrom][best_feature_type][best_feature_id]
                return best_feature_type, feature_data
            except (KeyError, TypeError):
                return None
        
        return None

    def annotate_batch_vectorized(self, batch_df: pd.DataFrame) -> pd.DataFrame:
        """
        Vectorized annotation of mutation batch.
        Much faster than row-by-row processing.
        
        Args:
            batch_df: DataFrame with Chromosome, Start columns
            
        Returns:
            Annotated DataFrame
        """
        # Initialize annotation columns with defaults
        batch_df['Gene_Location'] = 'Intergenic'
        batch_df['Gene_Name'] = 'Unknown'
        batch_df['Gene_Strand'] = 'NA'
        batch_df['Feature_Type'] = 'Intergenic'
        
        # Group by chromosome for efficiency
        for chrom, chrom_group in batch_df.groupby('Chromosome'):
            if chrom not in self.interval_trees:
                continue
            
            # Get indices for this chromosome
            indices = chrom_group.index
            positions = chrom_group['Start'].values
            
            # Annotate all positions for this chromosome efficiently
            annotations = []
            for pos in positions:
                result = self.get_best_feature_for_position(chrom, pos)
                
                if result:
                    feature_type, feature_data = result
                    annotations.append({
                        'location': feature_type,
                        'name': feature_data['name'],
                        'strand': feature_data['strand'],
                        'type': feature_data['type']
                    })
                else:
                    annotations.append({
                        'location': 'Intergenic',
                        'name': 'Unknown',
                        'strand': 'NA',
                        'type': 'Intergenic'
                    })
            
            # Update dataframe using vectorized assignment (FAST!)
            batch_df.loc[indices, 'Gene_Location'] = [a['location'] for a in annotations]
            batch_df.loc[indices, 'Gene_Name'] = [a['name'] for a in annotations]
            batch_df.loc[indices, 'Gene_Strand'] = [a['strand'] for a in annotations]
            batch_df.loc[indices, 'Feature_Type'] = [a['type'] for a in annotations]
        
        return batch_df

    def annotate_mutations(self, mutations_df: pd.DataFrame, 
                          batch_size: int = 50000) -> pd.DataFrame:
        """
        Annotate mutations with genomic features.
        Processes in batches for memory efficiency.
        
        Args:
            mutations_df: DataFrame of mutations
            batch_size: Batch size for processing
            
        Returns:
            Annotated DataFrame
        """
        logger.info(f"Annotating {len(mutations_df)} mutations")
        start_time = time.time()
        
        if len(mutations_df) == 0:
            return mutations_df
        
        result_df = mutations_df.copy()
        result_df['Gene_Location'] = 'Intergenic'
        result_df['Gene_Name'] = 'Unknown'
        result_df['Gene_Strand'] = 'NA'
        result_df['Feature_Type'] = 'Intergenic'

        num_batches = (len(result_df) + batch_size - 1) // batch_size
        logger.info(f"Processing in {num_batches} batches (size: {batch_size})")
        
        for batch_idx in range(num_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(result_df))
            
            batch = result_df.iloc[start_idx:end_idx].copy()
            batch_start = time.time()
            
            batch = self.annotate_batch_vectorized(batch)
            result_df.iloc[start_idx:end_idx] = batch
            
            batch_elapsed = time.time() - batch_start
            mutations_per_sec = len(batch) / max(0.001, batch_elapsed)
            
            logger.info(
                f"Batch {batch_idx+1}/{num_batches}: "
                f"{len(batch)} mutations in {batch_elapsed:.2f}s "
                f"({mutations_per_sec:.1f} mut/sec)"
            )
            
            del batch
            gc.collect()
        
        self._log_annotation_stats(result_df)
        
        total_elapsed = time.time() - start_time
        logger.info(f"Annotation completed in {total_elapsed:.2f}s")
        
        return result_df

    def _log_annotation_stats(self, df: pd.DataFrame):
        """Log annotation statistics."""
        location_counts = df['Gene_Location'].value_counts().to_dict()
        logger.info(f"Annotation summary: {location_counts}")
        
        top_genes = df['Gene_Name'].value_counts().head(5)
        logger.info(f"Top 5 annotated genes: {top_genes.to_dict()}")

    # ========================================================================
    # NEW: Quality & Transcription Analysis (Integrated)
    # ========================================================================
    
    def analyze_quality_transcription(self, annotated_df: pd.DataFrame) -> Dict:
        """
        Analyze quality and transcription status of mutations.
        Automatically detects quality/transcription columns.
        
        Args:
            annotated_df: DataFrame with annotations
            
        Returns:
            Dictionary with quality/transcription analysis results
        """
        logger.info("Analyzing quality and transcription status")
        
        results = {}
        quality_output_dir = os.path.join(self.output_dir, 'quality_transcription_analysis')
        os.makedirs(quality_output_dir, exist_ok=True)
        
        # Detect quality and transcription columns dynamically
        quality_col = None
        transcription_col = None
        
        for col in annotated_df.columns:
            if 'quality' in col.lower():
                quality_col = col
            if 'transcription' in col.lower():
                transcription_col = col
        
        if not quality_col and not transcription_col:
            logger.warning("No quality or transcription columns found - skipping analysis")
            return results
        
        # Quality Analysis
        if quality_col:
            logger.info(f"Analyzing quality column: {quality_col}")
            quality_counts = annotated_df[quality_col].value_counts()
            results['quality'] = quality_counts.to_dict()

        if quality_col and len(quality_counts) > 0:
            # Plot
            fig, ax = plt.subplots(figsize=(10, 6))
            colors = plt.cm.Set2(np.linspace(0, 1, len(quality_counts)))
            quality_counts.plot(kind='bar', ax=ax, color=colors)
            ax.set_title('Mutation Quality Distribution', fontsize=14, fontweight='bold')
            ax.set_xlabel('Quality Annotation', fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            quality_plot = os.path.join(quality_output_dir, 'quality_distribution.png')
            plt.savefig(quality_plot, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"  ✓ Saved: {quality_plot}")
            
            # CSV
            quality_df = quality_counts.reset_index()
            quality_df.columns = ['Quality_Annotation', 'Count']
            quality_df['Percentage'] = (quality_df['Count'] / quality_df['Count'].sum() * 100).round(2)
            quality_csv = os.path.join(quality_output_dir, 'quality_summary.csv')
            quality_df.to_csv(quality_csv, index=False)
            logger.info(f"  ✓ Saved: {quality_csv}")
        
        # Transcription Analysis
        if transcription_col:
            logger.info(f"Analyzing transcription column: {transcription_col}")
            transcription_counts = annotated_df[transcription_col].value_counts()
            results['transcription'] = transcription_counts.to_dict()

        if transcription_col and len(transcription_counts) > 0:
            # Plot
            fig, ax = plt.subplots(figsize=(10, 6))
            colors = plt.cm.Set3(np.linspace(0, 1, len(transcription_counts)))
            transcription_counts.plot(kind='bar', ax=ax, color=colors)
            ax.set_title('Mutation Transcription Status Distribution', fontsize=14, fontweight='bold')
            ax.set_xlabel('Transcription Status', fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            transcription_plot = os.path.join(quality_output_dir, 'transcription_distribution.png')
            plt.savefig(transcription_plot, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"  ✓ Saved: {transcription_plot}")
            
            # CSV
            transcription_df = transcription_counts.reset_index()
            transcription_df.columns = ['Transcription_Status', 'Count']
            transcription_df['Percentage'] = (transcription_df['Count'] / transcription_df['Count'].sum() * 100).round(2)
            transcription_csv = os.path.join(quality_output_dir, 'transcription_summary.csv')
            transcription_df.to_csv(transcription_csv, index=False)
            logger.info(f"  ✓ Saved: {transcription_csv}")
        
        # Cross-tabulation heatmap
        if quality_col and transcription_col:
            logger.info("Creating quality vs transcription heatmap")
            crosstab = pd.crosstab(annotated_df[quality_col], annotated_df[transcription_col])

            if crosstab.empty or crosstab.values.sum() == 0:
                logger.warning(
                    "Quality and transcription annotations are mutually exclusive "
                    "(no row has both populated) - skipping joint heatmap"
                )
                return results

            fig, ax = plt.subplots(figsize=(12, 8))
            sns.heatmap(crosstab, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
                       cbar_kws={'label': 'Count'}, linewidths=0.5)
            ax.set_title('Quality vs Transcription Status', fontsize=14, fontweight='bold')
            ax.set_xlabel('Transcription Status', fontsize=12)
            ax.set_ylabel('Quality Annotation', fontsize=12)
            plt.tight_layout()
            
            crosstab_plot = os.path.join(quality_output_dir, 'quality_vs_transcription_heatmap.png')
            plt.savefig(crosstab_plot, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"  ✓ Saved: {crosstab_plot}")
            
            # CSV
            crosstab_csv = os.path.join(quality_output_dir, 'quality_vs_transcription.csv')
            crosstab.to_csv(crosstab_csv)
            logger.info(f"  ✓ Saved: {crosstab_csv}")
        
        return results

    # ========================================================================
    # NEW: Indel Mechanism Analysis (Integrated)
    # ========================================================================
    
    def analyze_indel_mechanisms(self, annotated_df: pd.DataFrame) -> Dict:
        """
        Analyze indel types, mechanisms, and sizes.
        
        Args:
            annotated_df: DataFrame with annotations
            
        Returns:
            Dictionary with indel analysis results
        """
        logger.info("Analyzing indel mechanisms")

        results = {}

        if 'Mutation_Type' in annotated_df.columns:
            mut_types = set(annotated_df['Mutation_Type'].dropna().unique())
            if mut_types and 'ID' not in mut_types:
                logger.info(
                    f"Mutation_Type is {sorted(mut_types)} - no indels expected, "
                    f"skipping indel analysis"
                )
                return results

        ref_col = 'Ref' if 'Ref' in annotated_df.columns else 'Reference_Allele'
        alt_col = 'Alt' if 'Alt' in annotated_df.columns else 'Alternate_Allele'

        if ref_col not in annotated_df.columns or alt_col not in annotated_df.columns:
            logger.warning(
                f"Reference/alternate allele columns not found "
                f"(looked for Ref/Alt and Reference_Allele/Alternate_Allele) - "
                f"skipping indel analysis"
            )
            return results

        indel_output_dir = os.path.join(self.output_dir, 'indel_mechanism_analysis')
        os.makedirs(indel_output_dir, exist_ok=True)

        # Filter only indels
        indels = annotated_df[
            (annotated_df[ref_col].astype(str).str.len() !=
             annotated_df[alt_col].astype(str).str.len())
        ].copy()

        if len(indels) == 0:
            logger.warning("No indels found in dataset - skipping indel analysis")
            return results

        logger.info(f"Analyzing {len(indels)} indels")

        # Classify indel type
        indels['Indel_Type'] = indels.apply(
            lambda row: 'Deletion' if len(str(row[ref_col])) > len(str(row[alt_col]))
            else 'Insertion',
            axis=1
        )

        # Calculate indel size
        indels['Indel_Size'] = indels.apply(
            lambda row: abs(len(str(row[ref_col])) - len(str(row[alt_col]))),
            axis=1
        )
        
        # Classify frameshift
        indels['Is_Frameshift'] = indels['Indel_Size'] % 3 != 0
        
        # 1. Indel Type Distribution
        indel_type_counts = indels['Indel_Type'].value_counts()
        results['indel_type'] = indel_type_counts.to_dict()
        
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = ['#e74c3c', '#3498db']
        indel_type_counts.plot(kind='bar', ax=ax, color=colors)
        ax.set_title('Indel Type Distribution (Deletion vs Insertion)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Indel Type', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        plt.xticks(rotation=0)
        plt.tight_layout()
        
        indel_type_plot = os.path.join(indel_output_dir, 'indel_type_distribution.png')
        plt.savefig(indel_type_plot, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  ✓ Saved: {indel_type_plot}")
        
        # Indel Type CSV
        indel_type_df = indel_type_counts.reset_index()
        indel_type_df.columns = ['Indel_Type', 'Count']
        indel_type_df['Percentage'] = (indel_type_df['Count'] / indel_type_df['Count'].sum() * 100).round(2)
        indel_type_csv = os.path.join(indel_output_dir, 'indel_type_summary.csv')
        indel_type_df.to_csv(indel_type_csv, index=False)
        logger.info(f"  ✓ Saved: {indel_type_csv}")
        
        # 2. Indel Size Distribution
        fig, ax = plt.subplots(figsize=(10, 6))
        for indel_type in ['Deletion', 'Insertion']:
            subset = indels[indels['Indel_Type'] == indel_type]['Indel_Size']
            if len(subset) > 0:
                ax.hist(subset, bins=30, alpha=0.6, label=indel_type)
        
        ax.set_title('Indel Size Distribution', fontsize=14, fontweight='bold')
        ax.set_xlabel('Size (bp)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_yscale('log')
        ax.legend()
        plt.tight_layout()
        
        indel_size_plot = os.path.join(indel_output_dir, 'indel_size_distribution.png')
        plt.savefig(indel_size_plot, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  ✓ Saved: {indel_size_plot}")
        
        # 3. Location Analysis
        if 'Gene_Location' in annotated_df.columns:
            location_by_type = pd.crosstab(indels['Indel_Type'], indels['Gene_Location'])
            results['location_by_type'] = location_by_type.to_dict()
            
            fig, ax = plt.subplots(figsize=(12, 6))
            location_by_type.plot(kind='bar', ax=ax, stacked=False)
            ax.set_title('Indel Type by Genomic Location', fontsize=14, fontweight='bold')
            ax.set_xlabel('Indel Type', fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            ax.legend(title='Location', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            plt.xticks(rotation=0)
            plt.tight_layout()
            
            location_plot = os.path.join(indel_output_dir, 'indel_by_location.png')
            plt.savefig(location_plot, dpi=300, bbox_inches='tight')
            plt.close()
            logger.info(f"  ✓ Saved: {location_plot}")
            
            # CSV
            location_csv = os.path.join(indel_output_dir, 'indel_by_location.csv')
            location_by_type.to_csv(location_csv)
            logger.info(f"  ✓ Saved: {location_csv}")
        
        # 4. Frameshift Analysis
        frameshift_counts = indels['Is_Frameshift'].value_counts()
        results['frameshift'] = {
            'Frameshift': frameshift_counts.get(True, 0),
            'In-frame': frameshift_counts.get(False, 0)
        }
        
        fig, ax = plt.subplots(figsize=(8, 6))
        frameshift_labels = ['In-frame', 'Frameshift']
        frameshift_values = [
            frameshift_counts.get(False, 0),
            frameshift_counts.get(True, 0)
        ]
        colors_fs = ['#2ecc71', '#e74c3c']
        ax.bar(frameshift_labels, frameshift_values, color=colors_fs)
        ax.set_title('Frameshift vs In-frame Indels', fontsize=14, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12)
        for i, v in enumerate(frameshift_values):
            ax.text(i, v + max(frameshift_values)*0.02, str(v), ha='center', fontweight='bold')
        plt.tight_layout()
        
        frameshift_plot = os.path.join(indel_output_dir, 'frameshift_distribution.png')
        plt.savefig(frameshift_plot, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  ✓ Saved: {frameshift_plot}")
        
        # 5. Comprehensive Summary
        summary = pd.DataFrame({
            'Metric': [
                'Total Indels',
                'Deletions',
                'Insertions',
                'Frameshift Indels',
                'In-frame Indels',
                'Mean Size (bp)',
                'Median Size (bp)',
                'Max Size (bp)',
                'Min Size (bp)'
            ],
            'Count': [
                len(indels),
                len(indels[indels['Indel_Type'] == 'Deletion']),
                len(indels[indels['Indel_Type'] == 'Insertion']),
                frameshift_counts.get(True, 0),
                frameshift_counts.get(False, 0),
                round(indels['Indel_Size'].mean(), 2),
                round(indels['Indel_Size'].median(), 2),
                indels['Indel_Size'].max(),
                indels['Indel_Size'].min()
            ]
        })
        
        summary_csv = os.path.join(indel_output_dir, 'indel_summary.csv')
        summary.to_csv(summary_csv, index=False)
        logger.info(f"  ✓ Saved: {summary_csv}")
        
        results['summary'] = summary.to_dict('records')
        
        return results


def main():
    """Run complete mutation annotation pipeline"""
    parser = argparse.ArgumentParser(description='Annotate mutations with genomic features')
    parser.add_argument('--mutations-df', '-m', type=str, required=True, 
                       help='Pandas DataFrame pickle file with mutations')
    parser.add_argument('--annotation-dir', '-a', type=str, default='annotations',
                       help='Directory containing annotation files')
    parser.add_argument('--build', '-b', type=str, default='hg38',
                       help='Genome build version')
    parser.add_argument('--output-dir', '-o', type=str, default='mutation_annotation_results',
                       help='Output directory')
    
    args = parser.parse_args()
    
    try:
        logger.info("="*60)
        logger.info("MUTATION ANNOTATION PIPELINE")
        logger.info("="*60)
        
        # Initialize annotator
        annotator = OptimizedMutationAnnotator(
            annotation_dir=args.annotation_dir,
            build=args.build,
            output_dir=args.output_dir
        )
        
        # Load mutations
        logger.info(f"Loading mutations from: {args.mutations_df}")
        mutations_df = pd.read_pickle(args.mutations_df)
        logger.info(f"Loaded {len(mutations_df)} mutations")
        
        # Annotate mutations
        logger.info("Annotating mutations with genomic features")
        annotated_df = annotator.annotate_mutations(mutations_df)
        
        # Save annotated mutations
        output_pkl = os.path.join(args.output_dir, 'annotated_mutations.pkl')
        annotated_df.to_pickle(output_pkl)
        logger.info(f"✓ Saved annotated mutations: {output_pkl}")
        
        output_csv = os.path.join(args.output_dir, 'annotated_mutations.csv')
        annotated_df.to_csv(output_csv, index=False)
        logger.info(f"✓ Saved annotated mutations: {output_csv}")
        
        # Quality & Transcription Analysis
        logger.info("\n" + "="*60)
        logger.info("QUALITY & TRANSCRIPTION ANALYSIS")
        logger.info("="*60)
        quality_results = annotator.analyze_quality_transcription(annotated_df)
        
        # Indel Mechanism Analysis
        logger.info("\n" + "="*60)
        logger.info("INDEL MECHANISM ANALYSIS")
        logger.info("="*60)
        indel_results = annotator.analyze_indel_mechanisms(annotated_df)
        
        logger.info("\n" + "="*60)
        logger.info("✓ PIPELINE COMPLETE")
        logger.info("="*60)
        logger.info(f"Results saved to: {args.output_dir}")
        logger.info(f"  - Annotated mutations: annotated_mutations.csv")
        logger.info(f"  - Quality analysis: quality_transcription_analysis/")
        logger.info(f"  - Indel analysis: indel_mechanism_analysis/")
        
        return 0
    
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())