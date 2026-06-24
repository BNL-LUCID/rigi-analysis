#!/usr/bin/env python
"""
Optimized mutation preprocessing with efficient parsing and memory management.
"""

import os
import re
import glob
import pandas as pd
import numpy as np
import logging
import sys
import traceback
import warnings
import time
import argparse
from typing import List, Dict, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    handlers=[
        logging.FileHandler('preprocessing.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

class MutationParser:
    """Efficient mutation file parser with optimized pattern matching."""
    
    def __init__(self):
        """Initialize parser with precompiled patterns."""
        # Precompile all regex patterns once
        self.patterns = {
            'mutation': re.compile(r'\[?([ACGT]+)>([ACGT]+)\]?'),
            'qt_annotation': re.compile(r'^([NQTU]):'),
            'dose': re.compile(r'd([A-E0])', re.IGNORECASE),
            'timepoint': re.compile(r'W(\d+)')
        }
        
        self.feature_priority = {
            'exon': 5,
            '5utr': 4,
            '3utr': 3,
            'intron': 2,
            'gene': 1,
            'promoter': 0
        }

    def extract_dose(self, sample_name: str) -> str:
        """Extract dose from sample name."""
        if sample_name.lower().startswith('d0'):
            return 'Control'
        
        match = self.patterns['dose'].search(sample_name)
        return f"d{match.group(1)}" if match else 'Unknown'

    def extract_timepoint(self, sample_name: str) -> str:
        """Extract timepoint from sample name."""
        match = self.patterns['timepoint'].search(sample_name)
        return f"W{match.group(1)}" if match else 'Unknown'

    def extract_context_info(self, context: str) -> Tuple[Optional[str], Optional[str], str, str]:
        """
        Extract quality, transcription, and mutation info from context in one pass.
        
        Returns:
            (quality_annotation, transcription_annotation, ref, alt)
        """
        quality_annotation = None
        transcription_annotation = None
        ref = None
        alt = None
        
        # Extract quality/transcription annotation at start
        qt_match = self.patterns['qt_annotation'].match(context)
        if qt_match:
            annotation = qt_match.group(1)
            if annotation in ['N', 'Q']:
                quality_annotation = annotation
            elif annotation in ['T', 'U']:
                transcription_annotation = annotation
        
        # Extract mutation
        mut_match = self.patterns['mutation'].search(context)
        if mut_match:
            ref = mut_match.group(1)
            alt = mut_match.group(2)
        
        return quality_annotation, transcription_annotation, ref, alt

    def parse_snv_dbs(self, fields: List[str], chrom: str) -> Optional[Dict]:
        """Parse SNV/DBS mutation efficiently."""
        try:
            if len(fields) < 4:
                return None
            
            sample = fields[0]
            position = int(fields[2])
            context = fields[3]
            
            quality_annot, transc_annot, ref, alt = self.extract_context_info(context)
            
            if not ref or not alt:
                return None
            
            # Determine mutation type by ref/alt length
            if len(ref) == 1 and len(alt) == 1:
                full_mutation_type = 'SNV'
            elif len(ref) == 2 and len(alt) == 2:
                full_mutation_type = 'DBS'
            elif len(ref) >= 3 and len(alt) >= 3 and len(ref) == len(alt):
                full_mutation_type = 'MNS'
            else:
                full_mutation_type = 'OTHER'
            
            strand_info = int(fields[4]) if len(fields) > 4 else 0
            
            return {
                'chromosome': chrom,
                'start': position - 1,  # 0-based
                'end': position + len(ref),
                'sample': sample,
                'ref': ref,
                'alt': alt,
                'context': context,
                'strand': strand_info,
                'dose': self.extract_dose(sample),
                'timepoint': self.extract_timepoint(sample),
                'mutation_type': full_mutation_type,
                'quality_annotation': quality_annot,
                'transcription_annotation': transc_annot,
                'indel_type': None,
                'indel_mechanism': None,
                'indel_size': None,
                'repeat_length': None
            }
        except Exception as e:
            logging.debug(f"Error parsing SNV/DBS: {e}")
            return None

    def parse_id(self, fields: List[str], chrom: str) -> Optional[Dict]:
        """Parse ID (insertion/deletion) mutation."""
        try:
            if len(fields) < 7:
                return None
            
            sample = fields[0]
            position = int(fields[2])
            context_str = fields[3]
            ref = fields[4]
            alt = fields[5]
            strand_info = int(fields[6]) if len(fields) > 6 else 0
            
            # Parse context: "N:5:Del:M:2" or "Q:3:Ins:R:1"
            context_parts = context_str.split(':')
            
            quality_annotation = None
            indel_type = None
            indel_mechanism = None
            indel_size = None
            repeat_length = None
            
            if len(context_parts) >= 1 and context_parts[0] in ['N', 'Q']:
                quality_annotation = context_parts[0]
            
            if len(context_parts) >= 5:
                try:
                    indel_size = int(context_parts[1])
                    indel_type = context_parts[2]  # Del or Ins
                    indel_mechanism = context_parts[3]  # M, C, or R
                    repeat_length = int(context_parts[4])
                except (ValueError, IndexError):
                    pass
            
            if not indel_type:
                return None
            
            return {
                'chromosome': chrom,
                'start': position - 1,
                'end': position + len(ref),
                'sample': sample,
                'ref': ref,
                'alt': alt,
                'context': context_str,
                'strand': strand_info,
                'dose': self.extract_dose(sample),
                'timepoint': self.extract_timepoint(sample),
                'mutation_type': f"ID_{indel_type}",
                'quality_annotation': quality_annotation,
                'transcription_annotation': None,
                'indel_type': indel_type,
                'indel_mechanism': indel_mechanism,
                'indel_size': indel_size,
                'repeat_length': repeat_length
            }
        except Exception as e:
            logging.debug(f"Error parsing ID: {e}")
            return None

    def parse_line(self, line: str, chrom: str, mut_type: str) -> Optional[Dict]:
        """Parse single line based on mutation type."""
        line = line.strip()
        if not line:
            return None
        
        fields = line.split('\t')
        
        if mut_type in ['SNV', 'DBS', 'MNS']:
            return self.parse_snv_dbs(fields, chrom)
        elif mut_type == 'ID':
            return self.parse_id(fields, chrom)
        
        return None

    def read_mutation_file(self, file_path: str, chrom: str, mut_type: str, 
                          chunk_size: int = 10000) -> pd.DataFrame:
        """Read mutation file with streaming to save memory."""
        mutations = []
        
        try:
            with open(file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    parsed = self.parse_line(line, chrom, mut_type)
                    if parsed:
                        mutations.append(parsed)
                    
                    # Log progress every chunk_size lines
                    if line_num % chunk_size == 0:
                        logging.debug(f"Processed {line_num} lines from {file_path}")
        except Exception as e:
            logging.error(f"Error reading {file_path}: {e}")
            logging.error(traceback.format_exc())
        
        if not mutations:
            return pd.DataFrame()
        
        # Create DataFrame efficiently using dict list
        df = pd.DataFrame(mutations)
        
        # Rename columns for consistency
        df.columns = [
            'Chromosome', 'Start', 'End', 'Sample', 'Ref', 'Alt', 'Context',
            'Mutation_Strand', 'Dose', 'Timepoint', 'Mutation_Type',
            'Quality_Annotation', 'Transcription_Annotation',
            'Indel_Type', 'Indel_Mechanism', 'Indel_Size', 'Repeat_Length'
        ]
        
        # Create MutationID using vectorized string operations (much faster)
        df['MutationID'] = (
            df['Chromosome'].astype(str) + '_' +
            df['Start'].astype(str) + '_' +
            df['Ref'] + '_' +
            df['Alt'] + '_' +
            df['Sample']
        )
        
        return df


def read_mutation_files(data_dir: str, chromosome: Optional[str] = None,
                       mut_type: str = 'SNV') -> pd.DataFrame:
    """Read mutation files for specified chromosome and type.

    Searches `data_dir` recursively for `*_seqinfo.txt` so the script works
    regardless of how SigProfilerMatrixGenerator nests its output (the matrix
    generator typically writes
    `<input>/output/<basename>/<TYPE>/.../{chr}_seqinfo.txt`).
    """
    logging.info(f"Starting {mut_type} mutation file reading...")

    # Recursive glob — matches both flat (./SNV/1_seqinfo.txt) and nested
    # (./SNV/CNN96/<sample>/1_seqinfo.txt) layouts.
    if chromosome and chromosome.lower() != 'all':
        file_pattern = os.path.join(data_dir, "**", f"{chromosome}_seqinfo.txt")
    else:
        file_pattern = os.path.join(data_dir, "**", "*_seqinfo.txt")

    file_paths = sorted(glob.glob(file_pattern, recursive=True))
    
    if not file_paths:
        logging.warning(f"No {mut_type} files found: {file_pattern}")
        return pd.DataFrame()
    
    parser = MutationParser()
    all_mutations = []
    
    for file_path in file_paths:
        chrom = os.path.basename(file_path).split('_')[0]
        
        # Skip if specific chromosome requested
        if chromosome and chromosome.lower() != 'all' and chrom != chromosome:
            continue
        
        logging.info(f"Processing {chrom} for {mut_type}...")
        df = parser.read_mutation_file(file_path, chrom, mut_type)
        
        if not df.empty:
            all_mutations.append(df)
    
    # Concatenate all mutations efficiently
    if not all_mutations:
        logging.warning(f"No {mut_type} mutations found")
        return pd.DataFrame()
    
    mutations_df = pd.concat(all_mutations, ignore_index=True)
    
    logging.info(f"Processed {len(mutations_df)} {mut_type} mutations")
    
    # Log annotation statistics
    _log_annotation_stats(mutations_df, mut_type)
    
    return mutations_df


def _log_annotation_stats(df: pd.DataFrame, mut_type: str):
    """Log annotation statistics for QC."""
    if 'Quality_Annotation' in df.columns:
        quality_counts = df['Quality_Annotation'].value_counts().to_dict()
        logging.info(f"Quality annotation counts: {quality_counts}")
    
    if 'Transcription_Annotation' in df.columns:
        transc_counts = df['Transcription_Annotation'].value_counts().to_dict()
        logging.info(f"Transcription annotation counts: {transc_counts}")
    
    if mut_type == 'ID' and 'Indel_Mechanism' in df.columns:
        mechanism_counts = df['Indel_Mechanism'].value_counts().to_dict()
        logging.info(f"Indel mechanism counts: {mechanism_counts}")


def create_annotation_summary(mutations_df: pd.DataFrame, mut_type: str) -> pd.DataFrame:
    """Create annotation summary for QC."""
    summaries = []
    
    if mut_type in ['SNV', 'DBS', 'MNS']:
        # Quality annotation
        if not mutations_df['Quality_Annotation'].isna().all():
            quality = mutations_df['Quality_Annotation'].value_counts().reset_index()
            quality.columns = ['Annotation', 'Count']
            quality['Type'] = 'Quality'
            summaries.append(quality)
        
        # Transcription annotation
        if not mutations_df['Transcription_Annotation'].isna().all():
            transc = mutations_df['Transcription_Annotation'].value_counts().reset_index()
            transc.columns = ['Annotation', 'Count']
            transc['Type'] = 'Transcription'
            summaries.append(transc)
    
    elif mut_type == 'ID':
        # Indel type
        if not mutations_df['Indel_Type'].isna().all():
            indel_type = mutations_df['Indel_Type'].value_counts().reset_index()
            indel_type.columns = ['Annotation', 'Count']
            indel_type['Type'] = 'Indel_Type'
            summaries.append(indel_type)
        
        # Mechanism
        if not mutations_df['Indel_Mechanism'].isna().all():
            mechanism = mutations_df['Indel_Mechanism'].value_counts().reset_index()
            mechanism.columns = ['Annotation', 'Count']
            mechanism['Type'] = 'Mechanism'
            summaries.append(mechanism)
    
    return pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()


def main():
    """Main preprocessing function."""
    parser = argparse.ArgumentParser(
        description='Process mutation files for analysis.',
        epilog=(
            "Example: point --input-dir at SigProfilerMatrixGenerator's output "
            "tree (typically <vcf_dir>/output/<vcf_dir_basename>/) so the "
            "SNV/ DBS/ ID/ MNS/ subfolders it created are picked up directly."
        ),
    )
    parser.add_argument('-i', '--input-dir', type=str, default='.',
                       help=("Directory containing SNV/, DBS/, ID/, MNS/ "
                             "subfolders (default: current directory). "
                             "Each subfolder is searched recursively for "
                             "*_seqinfo.txt files."))
    parser.add_argument('-c', '--chromosome', type=str, default='all',
                       help='Chromosome to process (default: all)')
    parser.add_argument('-o', '--output', type=str, default='processed_data',
                       help='Output directory')
    parser.add_argument('-s', '--summary', type=str, default='summary_data',
                       help='Summary directory')
    parser.add_argument('-t', '--mutation-types', nargs='+',
                       choices=['SNV', 'DBS', 'MNS', 'ID'],
                       default=['SNV', 'DBS', 'MNS', 'ID'],
                       help=("Which mutation types to process. Default: all "
                             "four. Useful when SigProfiler has finished some "
                             "types but not others, or to run multiple "
                             "instances in parallel — e.g. one process for "
                             "SNV DBS MNS, another for ID."))

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    os.makedirs(args.summary, exist_ok=True)

    all_dirs = {
        'SNV': os.path.join(args.input_dir, 'SNV'),
        'DBS': os.path.join(args.input_dir, 'DBS'),
        'MNS': os.path.join(args.input_dir, 'MNS'),
        'ID':  os.path.join(args.input_dir, 'ID'),
    }
    mutation_dirs = {t: all_dirs[t] for t in args.mutation_types}
    logging.info(f"Processing mutation types: {list(mutation_dirs)}")

    for mut_type, data_dir in mutation_dirs.items():
        if not os.path.exists(data_dir):
            logging.warning(f"Directory not found: {data_dir}")
            continue
        
        start_time = time.time()
        logging.info(f"Processing {mut_type}...")
        
        # Read mutations
        mutations_df = read_mutation_files(data_dir, args.chromosome, mut_type)
        
        if mutations_df.empty:
            logging.warning(f"No {mut_type} mutations found")
            continue
        
        # Save mutations
        csv_file = os.path.join(args.output, f"{args.chromosome}_{mut_type}_mutations.csv")
        pkl_file = os.path.join(args.output, f"{args.chromosome}_{mut_type}_mutations.pkl")
        
        mutations_df.to_csv(csv_file, index=False)
        mutations_df.to_pickle(pkl_file)
        
        # Save summary
        summary_df = create_annotation_summary(mutations_df, mut_type)
        if not summary_df.empty:
            summary_file = os.path.join(args.summary, f"{args.chromosome}_{mut_type}_summary.csv")
            summary_df.to_csv(summary_file, index=False)
            logging.info(f"Saved summary to {summary_file}")
        
        elapsed = time.time() - start_time
        logging.info(f"Saved {len(mutations_df)} {mut_type} mutations in {elapsed:.2f}s")
    
    logging.info("Preprocessing completed")


if __name__ == "__main__":
    main()