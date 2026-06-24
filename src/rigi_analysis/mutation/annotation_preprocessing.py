#!/usr/bin/env python
"""Optimized annotation preprocessing with streaming and efficient storage."""

import argparse
import gc
import gzip
import logging
import os
import pickle
import sys
import time
from collections import defaultdict

import pandas as pd
import requests
from intervaltree import IntervalTree

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    handlers=[
        logging.FileHandler('annotation_preprocessing.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

class AnnotationPreprocessor:
    """Efficient annotation preprocessing with streaming and caching."""

    def __init__(self, build="hg38", annotation_dir="annotations"):
        self.build = build
        self.annotation_dir = annotation_dir
        os.makedirs(annotation_dir, exist_ok=True)
        self.logger = logging.getLogger(__name__)

        # Feature priority (higher = higher priority)
        self.feature_priority = {
            'exon': 5,
            '5utr': 4,
            '3utr': 3,
            'intron': 2,
            'gene': 1,
            'promoter': 0
        }
        self.promoter_region = 2000

    def download_ucsc_annotations(self):
        """Download gene annotations from UCSC."""
        refgene_file = os.path.join(self.annotation_dir, "refGene.txt.gz")

        if os.path.exists(refgene_file):
            self.logger.info(f"Using existing file: {refgene_file}")
            return refgene_file

        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{self.build}/database/refGene.txt.gz"
        self.logger.info(f"Downloading from {url}...")

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(refgene_file, 'wb') as f:
                f.write(response.content)

            self.logger.info(f"Downloaded to {refgene_file}")
            return refgene_file
        except Exception as e:
            self.logger.error(f"Download failed: {e}")
            raise

    def stream_refgene_file(self, refgene_file, chunk_size=10000):
        """Stream process refGene file to avoid loading everything into memory."""
        self.logger.info(f"Streaming annotations from {refgene_file}...")

        columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
                  'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']

        with gzip.open(refgene_file, 'rt') as f:
            chunk = []
            for line_num, line in enumerate(f, 1):
                values = line.strip().split('\t')
                chunk.append(dict(zip(columns, values)))

                if len(chunk) >= chunk_size:
                    yield pd.DataFrame(chunk)
                    chunk = []
                    self.logger.info(f"Processed {line_num} transcripts...")

            if chunk:
                yield pd.DataFrame(chunk)

    def process_refgene_streaming(self, refgene_file):
        """Process refGene file with streaming and build trees directly."""
        self.logger.info("Processing annotations with streaming...")

        interval_trees = defaultdict(lambda: defaultdict(IntervalTree))
        feature_data = defaultdict(lambda: defaultdict(dict))

        processed_count = 0

        for chunk_df in self.stream_refgene_file(refgene_file):
            processed_count += len(chunk_df)

            for _, transcript in chunk_df.iterrows():
                chrom = transcript['chrom'].replace('chr', '')
                transcript['name']
                gene_name = transcript['name2']
                strand = transcript['strand']
                tx_start = int(transcript['txStart'])
                tx_end = int(transcript['txEnd'])
                cds_start = int(transcript['cdsStart'])
                cds_end = int(transcript['cdsEnd'])

                # Add gene entry
                self._add_feature(interval_trees, feature_data, chrom,
                                'gene', tx_start, tx_end, gene_name, strand)

                # Process exons, UTRs, introns
                exon_starts = [int(x) for x in transcript['exonStarts'].strip(',').split(',') if x]
                exon_ends = [int(x) for x in transcript['exonEnds'].strip(',').split(',') if x]

                for i, (exon_start, exon_end) in enumerate(zip(exon_starts, exon_ends)):
                    # Add exon
                    self._add_feature(interval_trees, feature_data, chrom,
                                    'exon', exon_start, exon_end, gene_name, strand)

                    # Add UTRs
                    if cds_start > exon_start:
                        utr_end = min(exon_end, cds_start)
                        utr_type = '5utr' if strand == '+' else '3utr'
                        self._add_feature(interval_trees, feature_data, chrom,
                                        utr_type, exon_start, utr_end, gene_name, strand)

                    if cds_end < exon_end:
                        utr_start = max(exon_start, cds_end)
                        utr_type = '3utr' if strand == '+' else '5utr'
                        self._add_feature(interval_trees, feature_data, chrom,
                                        utr_type, utr_start, exon_end, gene_name, strand)

                # Add introns
                for i in range(len(exon_starts) - 1):
                    intron_start = exon_ends[i]
                    intron_end = exon_starts[i + 1]
                    self._add_feature(interval_trees, feature_data, chrom,
                                    'intron', intron_start, intron_end, gene_name, strand)

                # Add promoter
                if strand == '+':
                    prom_start = max(0, tx_start - self.promoter_region)
                    prom_end = tx_start
                else:
                    prom_start = tx_end
                    prom_end = tx_end + self.promoter_region

                self._add_feature(interval_trees, feature_data, chrom,
                                'promoter', prom_start, prom_end, gene_name, strand)

            # Free memory
            gc.collect()

        self.logger.info(f"Processed {processed_count} transcripts")
        return interval_trees, feature_data

    def _add_feature(self, interval_trees, feature_data, chrom, feature_type,
                    start, end, name, strand):
        """Add a feature to trees and data structures."""
        # Fix zero-length intervals
        if start >= end:
            end = start + 1

        priority = self.feature_priority.get(feature_type, 0)
        feature_id = len(feature_data[chrom][feature_type])

        # Store minimal data in feature_data
        feature_data[chrom][feature_type][feature_id] = {
            'name': name,
            'strand': strand,
            'type': feature_type
        }

        # Store interval with priority and ID
        interval_trees[chrom][feature_type][start:end] = (priority, feature_id)

    def save_trees(self, interval_trees, feature_data):
        """Save trees and data using pickle."""
        trees_file = os.path.join(self.annotation_dir, f"{self.build}_interval_trees.pkl")
        data_file = os.path.join(self.annotation_dir, f"{self.build}_feature_data.pkl")

        with open(trees_file, 'wb') as f:
            pickle.dump(dict(interval_trees), f)

        with open(data_file, 'wb') as f:
            pickle.dump(dict(feature_data), f)

        self.logger.info(f"Saved trees to {trees_file}")
        self.logger.info(f"Saved data to {data_file}")

    def run(self, force=False):
        """Run complete preprocessing pipeline."""
        start_time = time.time()

        self.logger.info(f"Starting annotation preprocessing for {self.build}")

        # Download
        refgene_file = self.download_ucsc_annotations()

        # Process with streaming
        interval_trees, feature_data = self.process_refgene_streaming(refgene_file)

        # Save
        self.save_trees(interval_trees, feature_data)

        elapsed = time.time() - start_time
        self.logger.info(f"Completed in {elapsed:.2f} seconds")


def main():
    parser = argparse.ArgumentParser(description='Download and process genomic annotations.')
    parser.add_argument('--build', default="hg38", help='Genome build (default: hg38)')
    parser.add_argument('--annotation-dir', default="annotations", help='Annotation directory')
    parser.add_argument('--force', action='store_true', help='Force reprocessing')

    args = parser.parse_args()

    preprocessor = AnnotationPreprocessor(args.build, args.annotation_dir)
    preprocessor.run(force=args.force)


if __name__ == "__main__":
    main()
