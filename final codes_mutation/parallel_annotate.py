#!/usr/bin/env python
"""
Parallel wrapper around mutation_annotation.py.

Splits the input mutations pickle into N roughly-equal chunks (by row count,
preserving chromosome ordering), runs mutation_annotation.py on each chunk
in parallel as a subprocess, then concatenates the per-chunk annotated
outputs into a single combined annotated pkl + csv.

The annotation script itself is not modified — this is purely an outer loop.

Typical use (e.g. SNV on a 72-core box, 24 workers):

    python parallel_annotate.py \
        --input        processed_data/all_SNV_mutations.pkl \
        --annotation-dir annotations \
        --build        hg38 \
        --output       annotated_mutations/SNV/all_SNV_annotated.pkl \
        --workers      24
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd


def split_input(input_pkl: str, work_dir: str, n_chunks: int) -> list:
    print(f"Loading {input_pkl} ...")
    t0 = time.time()
    df = pd.read_pickle(input_pkl)
    print(f"  {len(df):,} rows in {time.time() - t0:.1f}s")

    # Sort by chromosome so each chunk holds contiguous chromosome blocks —
    # makes per-chunk logs easier to skim. Doesn't affect correctness.
    if 'Chromosome' in df.columns:
        df = df.sort_values('Chromosome', kind='stable').reset_index(drop=True)

    chunk_size = (len(df) + n_chunks - 1) // n_chunks
    chunks_dir = os.path.join(work_dir, 'chunks')
    os.makedirs(chunks_dir, exist_ok=True)

    chunk_paths = []
    for i in range(n_chunks):
        start = i * chunk_size
        end = min(start + chunk_size, len(df))
        if start >= len(df):
            break
        chunk_df = df.iloc[start:end].copy()
        chunk_path = os.path.join(chunks_dir, f'chunk_{i:03d}.pkl')
        chunk_df.to_pickle(chunk_path)
        chunk_paths.append(chunk_path)
        chrom_summary = ''
        if 'Chromosome' in chunk_df.columns:
            chrom_summary = f" chroms={chunk_df['Chromosome'].nunique()}"
        print(f"  chunk {i:03d}: {len(chunk_df):,} rows{chrom_summary}")

    return chunk_paths


def run_chunk(chunk_path: str, annotation_dir: str, build: str,
              script: str, work_dir: str) -> dict:
    """Run mutation_annotation.py on a single chunk in a subprocess."""
    chunk_path = os.path.abspath(chunk_path)
    annotation_dir = os.path.abspath(annotation_dir)
    script = os.path.abspath(script)

    chunk_name = os.path.basename(chunk_path).replace('.pkl', '')
    chunk_out_dir = os.path.abspath(os.path.join(work_dir, 'out', chunk_name))
    os.makedirs(chunk_out_dir, exist_ok=True)
    log_path = os.path.join(chunk_out_dir, 'subprocess.log')

    cmd = [
        sys.executable, script,
        '-m', chunk_path,
        '-a', annotation_dir,
        '-b', build,
        '-o', chunk_out_dir,
    ]

    t0 = time.time()
    # cwd=chunk_out_dir isolates mutation_annotation.py's hardcoded
    # FileHandler('mutation_annotation.log') so subprocesses don't race on
    # a shared log file.
    with open(log_path, 'w') as logf:
        result = subprocess.run(
            cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=chunk_out_dir,
        )
    elapsed = time.time() - t0

    annotated_pkl = os.path.join(chunk_out_dir, 'annotated_mutations.pkl')
    return {
        'chunk': chunk_name,
        'returncode': result.returncode,
        'elapsed': elapsed,
        'output_pkl': annotated_pkl if os.path.exists(annotated_pkl) else None,
        'log': log_path,
        'out_dir': chunk_out_dir,
    }


def merge_outputs(chunk_results: list, output_pkl: str, write_csv: bool = True):
    print(f"\nMerging {len(chunk_results)} chunk outputs ...")
    dfs = []
    for r in chunk_results:
        if r['output_pkl'] is None:
            print(f"  MISSING {r['chunk']} (see {r['log']})")
            continue
        df = pd.read_pickle(r['output_pkl'])
        dfs.append(df)
        print(f"  loaded {r['chunk']}: {len(df):,} rows")

    if not dfs:
        print("ERROR: no chunk outputs found, nothing to merge")
        sys.exit(1)

    combined = pd.concat(dfs, ignore_index=True)
    print(f"  combined: {len(combined):,} rows")

    os.makedirs(os.path.dirname(output_pkl) or '.', exist_ok=True)
    combined.to_pickle(output_pkl)
    print(f"Wrote {output_pkl}")
    if write_csv:
        output_csv = os.path.splitext(output_pkl)[0] + '.csv'
        print(f"Writing CSV ({len(combined):,} rows) ...")
        combined.to_csv(output_csv, index=False)
        print(f"Wrote {output_csv}")
    else:
        print("Skipped CSV write (--no-csv). Convert later with pandas if needed.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', '-m', required=True,
                   help='Input mutations pickle (e.g. all_SNV_mutations.pkl)')
    p.add_argument('--annotation-dir', '-a', required=True,
                   help='Directory with hg38_interval_trees.pkl and feature_data.pkl')
    p.add_argument('--build', '-b', default='hg38')
    p.add_argument('--output', '-o', required=True,
                   help='Final combined annotated pkl path. CSV is written alongside.')
    p.add_argument('--workers', '-w', type=int, default=24,
                   help='Parallel workers (default 24). Each loads its own '
                        'copy of the interval-tree annotation data, so RAM '
                        'scales with workers — start at 24 on a 72-core box, '
                        'raise if RAM allows.')
    p.add_argument('--script', default=os.path.join(os.path.dirname(__file__),
                                                    'mutation_annotation.py'),
                   help='Path to mutation_annotation.py (default: alongside this script)')
    p.add_argument('--work-dir', default=None,
                   help='Scratch dir for chunks + per-chunk outputs '
                        '(default: <output_basename>_parallel_work)')
    p.add_argument('--keep-work-dir', action='store_true',
                   help='Keep scratch chunks/outputs after merging (for debugging). '
                        'Default deletes after successful merge.')
    p.add_argument('--no-csv', action='store_true',
                   help='Skip writing the CSV alongside the PKL. Useful for large '
                        'tables where the CSV would be tens of GB.')
    args = p.parse_args()

    # Validate --output is a file path, not a directory. Common mistake to
    # pass a dir like "annotated/" — merge_outputs would crash at write time
    # after hours of work. If it looks dir-shaped, append a default filename.
    if args.output.endswith('/') or args.output.endswith(os.sep) or os.path.isdir(args.output):
        default_name = os.path.basename(args.input).replace('_mutations.pkl', '_annotated.pkl')
        if not default_name.endswith('.pkl'):
            default_name = 'annotated_mutations.pkl'
        new_output = os.path.join(args.output.rstrip('/' + os.sep), default_name)
        print(f"--output looks like a directory; writing to {new_output}")
        args.output = new_output
    if not args.output.endswith('.pkl'):
        sys.exit(f"--output must end in .pkl (got: {args.output})")

    work_dir = args.work_dir or (os.path.splitext(args.output)[0] + '_parallel_work')
    os.makedirs(work_dir, exist_ok=True)
    print(f"Work dir: {work_dir}")
    print(f"Output:   {args.output}")
    print(f"Workers:  {args.workers}")
    print(f"Script:   {args.script}")

    overall_t0 = time.time()

    t0 = time.time()
    chunk_paths = split_input(args.input, work_dir, args.workers)
    print(f"Split done in {time.time() - t0:.1f}s")

    print(f"\nRunning {len(chunk_paths)} chunks across {args.workers} workers ...")
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(run_chunk, cp, args.annotation_dir, args.build,
                      args.script, work_dir): cp
            for cp in chunk_paths
        }
        done = 0
        for f in as_completed(futures):
            r = f.result()
            done += 1
            status = 'OK' if r['returncode'] == 0 else f'FAIL({r["returncode"]})'
            print(f"  [{done}/{len(chunk_paths)}] {r['chunk']} {status} "
                  f"({r['elapsed']:.1f}s)")
            results.append(r)
    print(f"Annotation done in {time.time() - t0:.1f}s")

    failures = [r for r in results if r['returncode'] != 0]
    if failures:
        print(f"\n{len(failures)} chunks failed:")
        for r in failures:
            print(f"  {r['chunk']}: see {r['log']}")
        sys.exit(1)

    merge_outputs(results, args.output, write_csv=not args.no_csv)

    if not args.keep_work_dir:
        print(f"\nCleaning up {work_dir} ...")
        shutil.rmtree(work_dir, ignore_errors=True)

    print(f"\nTotal wall clock: {time.time() - overall_t0:.1f}s")


if __name__ == '__main__':
    main()
