#!/usr/bin/env python3
"""
Filter AnnotSV files to keep only PASS calls.
Creates a new directory with filtered files.
"""

import pandas as pd
from pathlib import Path
import argparse
import os

def filter_annotsv_file(input_path, output_path):
    """
    Filter a single AnnotSV file to keep only PASS calls.
    Returns tuple: (total_rows, pass_rows, total_sv, pass_sv)
    """
    # Read the file
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    
    total_rows = len(df)
    
    # Count unique SVs before filtering (using full annotation rows)
    if 'Annotation_mode' in df.columns:
        total_sv = len(df[df['Annotation_mode'] == 'full'])
    else:
        total_sv = total_rows
    
    # Filter to PASS only
    if 'FILTER' in df.columns:
        df_pass = df[df['FILTER'] == 'PASS'].copy()
    else:
        print(f"  WARNING: No FILTER column found in {input_path.name}")
        df_pass = df.copy()
    
    pass_rows = len(df_pass)
    
    # Count unique SVs after filtering
    if 'Annotation_mode' in df_pass.columns:
        pass_sv = len(df_pass[df_pass['Annotation_mode'] == 'full'])
    else:
        pass_sv = pass_rows
    
    # Save filtered file
    df_pass.to_csv(output_path, sep='\t', index=False)
    
    return total_rows, pass_rows, total_sv, pass_sv


def main():
    parser = argparse.ArgumentParser(description='Filter AnnotSV files to PASS calls only')
    parser.add_argument('--input-dir', required=True, help='Directory containing AnnotSV files')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered files')
    parser.add_argument('--pattern', default='*_annotated.tsv', help='File pattern to match (default: *_annotated.tsv)')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all AnnotSV files
    annotsv_files = list(input_dir.glob(args.pattern))
    
    if not annotsv_files:
        print(f"No files matching '{args.pattern}' found in {input_dir}")
        return
    
    print(f"Found {len(annotsv_files)} AnnotSV files")
    print(f"Output directory: {output_dir}")
    print("=" * 70)
    
    # Track statistics
    stats = []
    
    for input_path in sorted(annotsv_files):
        # Keep original filename since it's in a separate directory
        output_path = output_dir / input_path.name
        
        print(f"\nProcessing: {input_path.name}")
        
        total_rows, pass_rows, total_sv, pass_sv = filter_annotsv_file(input_path, output_path)
        
        pct_rows = 100 * pass_rows / total_rows if total_rows > 0 else 0
        pct_sv = 100 * pass_sv / total_sv if total_sv > 0 else 0
        
        print(f"  Rows: {total_rows:,} → {pass_rows:,} ({pct_rows:.1f}% kept)")
        print(f"  SVs:  {total_sv:,} → {pass_sv:,} ({pct_sv:.1f}% kept)")
        print(f"  Saved: {output_path.name}")
        
        stats.append({
            'file': input_path.name,
            'total_rows': total_rows,
            'pass_rows': pass_rows,
            'total_sv': total_sv,
            'pass_sv': pass_sv,
            'pct_sv_kept': pct_sv
        })
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    total_sv_all = sum(s['total_sv'] for s in stats)
    pass_sv_all = sum(s['pass_sv'] for s in stats)
    
    print(f"\nTotal files processed: {len(stats)}")
    print(f"Total SVs across all files: {total_sv_all:,}")
    print(f"PASS SVs across all files: {pass_sv_all:,}")
    print(f"Overall retention rate: {100 * pass_sv_all / total_sv_all:.1f}%")
    
    # Save summary to file
    summary_path = output_dir / 'filtering_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("AnnotSV PASS Filtering Summary\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Input directory: {input_dir}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Files processed: {len(stats)}\n\n")
        f.write(f"{'File':<50} {'Total SVs':>12} {'PASS SVs':>12} {'% Kept':>10}\n")
        f.write("-" * 84 + "\n")
        for s in stats:
            f.write(f"{s['file']:<50} {s['total_sv']:>12,} {s['pass_sv']:>12,} {s['pct_sv_kept']:>10.1f}%\n")
        f.write("-" * 84 + "\n")
        f.write(f"{'TOTAL':<50} {total_sv_all:>12,} {pass_sv_all:>12,} {100*pass_sv_all/total_sv_all:>10.1f}%\n")
    
    print(f"\nSummary saved to: {summary_path}")
    print(f"\nFiltered files are in: {output_dir}")
    print("\nNext steps:")
    print("  1. Re-run your INV-DBS analysis using the filtered files")
    print("  2. Re-do pathway enrichment on PASS-only inversions")


if __name__ == "__main__":
    main()