#!/usr/bin/env python
r"""Dose-level parallel wrapper around mutation_pattern_assignment.py.

Each dose's pattern computation is independent (different filtered subset,
different output dir), so we can fan out across processes — one subprocess
per dose. The pattern_assignment script itself is not modified; this helper
just pre-splits the annotated input by dose and launches parallel runs.

The split per-dose input contains:
  - all rows of the control dose(s) (matched against the same regex the
    script uses internally: d0 / D0 / Control)
  - all rows for that one dose

Each subprocess sees only its dose + controls in the input, so the script's
internal dose loop runs once and writes to <output-dir>/<TYPE>/dose_<dose>/.

Typical use (SNV with 5 doses, all in parallel):

    python parallel_pattern_assignment.py \\
        --input        annotated_mutations/SNV/all_SNV_annotated.pkl \\
        --output-dir   pattern_analysis/SNV \\
        --mutation-type SNV \\
        --workers      5
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

CONTROL_PATTERN = re.compile(r'(?:d0|D0|[Cc]ontrol)')


def split_by_dose(input_pkl: str, work_dir: str, mutation_type: str) -> list:
    """Load annotated pickle and write one input pickle per non-control dose
    (each containing controls + that dose). Returns list of (dose, dir) pairs.
    """
    print(f"Loading {input_pkl} ...")
    t0 = time.time()
    df = pd.read_pickle(input_pkl)
    print(f"  {len(df):,} rows in {time.time() - t0:.1f}s")

    if 'Dose' not in df.columns:
        sys.exit("ERROR: input pickle has no 'Dose' column")

    all_doses = df['Dose'].dropna().unique().tolist()
    control_doses = [d for d in all_doses if CONTROL_PATTERN.match(str(d))]
    treated_doses = sorted(set(all_doses) - set(control_doses))

    if not control_doses:
        print(f"  WARNING: no control dose detected (matched against {CONTROL_PATTERN.pattern})")
    print(f"  Controls: {control_doses}")
    print(f"  Treated:  {treated_doses}")

    control_mask = df['Dose'].isin(control_doses)
    control_rows = df[control_mask]
    print(f"  Control rows: {len(control_rows):,}")

    dose_inputs = []
    for dose in treated_doses:
        dose_mask = df['Dose'] == dose
        dose_df = pd.concat([control_rows, df[dose_mask]], ignore_index=True)

        # Each dose gets its own subdir with one pickle named to match the
        # script's input glob: *_<TYPE>_annotated.pkl
        dose_dir = os.path.join(work_dir, f'input_{dose}')
        os.makedirs(dose_dir, exist_ok=True)
        dose_pkl = os.path.join(dose_dir, f'all_{mutation_type}_annotated.pkl')
        dose_df.to_pickle(dose_pkl)
        print(f"  {dose}: {len(dose_df):,} rows (control + dose) -> {dose_pkl}")
        dose_inputs.append((dose, dose_dir))

    return dose_inputs


def run_dose(dose: str, dose_input_dir: str, output_dir: str,
             mutation_type: str, script: str, work_dir: str) -> dict:
    """Run mutation_pattern_assignment.py on a single dose's pre-split input."""
    dose_input_dir = os.path.abspath(dose_input_dir)
    output_dir = os.path.abspath(output_dir)
    script = os.path.abspath(script)

    log_dir = os.path.abspath(os.path.join(work_dir, 'logs'))
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f'{dose}.log')

    # Each subprocess gets its own CWD so any incidental relative-path writes
    # don't collide across workers.
    cwd = os.path.abspath(os.path.join(work_dir, f'cwd_{dose}'))
    os.makedirs(cwd, exist_ok=True)

    cmd = [
        sys.executable, script,
        '-i', dose_input_dir,
        '-o', output_dir,
        '-m', mutation_type,
    ]

    t0 = time.time()
    with open(log_path, 'w') as logf:
        result = subprocess.run(
            cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=cwd,
        )
    elapsed = time.time() - t0

    return {
        'dose': dose,
        'returncode': result.returncode,
        'elapsed': elapsed,
        'log': log_path,
    }


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--input', '-i', required=True,
                   help='Annotated mutations pickle (e.g. all_SNV_annotated.pkl)')
    p.add_argument('--output-dir', '-o', required=True,
                   help='Top-level pattern output dir (e.g. pattern_analysis/SNV). '
                        'Each dose writes to <output-dir>/dose_<dose>/.')
    p.add_argument('--mutation-type', '-m', required=True,
                   choices=['SNV', 'DBS', 'MNS', 'ID'])
    p.add_argument('--workers', '-w', type=int, default=5,
                   help='Parallel workers (default 5, one per dose). Each loads '
                        'its own filtered DataFrame, so RAM scales with workers.')
    p.add_argument('--script', default=os.path.join(
                       os.path.dirname(__file__), 'mutation_pattern_assignment.py'),
                   help='Path to mutation_pattern_assignment.py '
                        '(default: alongside this script)')
    p.add_argument('--work-dir', default=None,
                   help='Scratch dir for per-dose split inputs + logs '
                        '(default: <output-dir>_parallel_work)')
    p.add_argument('--keep-work-dir', action='store_true',
                   help='Keep scratch dir after run (for debugging). '
                        'Default deletes after all doses succeed.')
    args = p.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input not found: {args.input}")
    if not os.path.exists(args.script):
        sys.exit(f"ERROR: script not found: {args.script}")

    work_dir = args.work_dir or (args.output_dir.rstrip('/') + '_parallel_work')
    os.makedirs(work_dir, exist_ok=True)
    print(f"Work dir:   {work_dir}")
    print(f"Output dir: {args.output_dir}")
    print(f"Workers:    {args.workers}")
    print(f"Script:     {args.script}")

    overall_t0 = time.time()

    t0 = time.time()
    dose_inputs = split_by_dose(args.input, work_dir, args.mutation_type)
    print(f"Split done in {time.time() - t0:.1f}s")

    if not dose_inputs:
        sys.exit("ERROR: no treated doses to process")

    print(f"\nRunning {len(dose_inputs)} doses across {args.workers} workers ...")
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(run_dose, dose, dose_dir, args.output_dir,
                      args.mutation_type, args.script, work_dir): dose
            for dose, dose_dir in dose_inputs
        }
        done = 0
        for f in as_completed(futures):
            r = f.result()
            done += 1
            status = 'OK' if r['returncode'] == 0 else f'FAIL({r["returncode"]})'
            print(f"  [{done}/{len(dose_inputs)}] dose={r['dose']} {status} "
                  f"({r['elapsed']:.1f}s)")
            results.append(r)
    print(f"All doses done in {time.time() - t0:.1f}s")

    failures = [r for r in results if r['returncode'] != 0]
    if failures:
        print(f"\n{len(failures)} doses failed:")
        for r in failures:
            print(f"  dose={r['dose']}: see {r['log']}")
        sys.exit(1)

    if not args.keep_work_dir:
        print(f"\nCleaning up {work_dir} ...")
        shutil.rmtree(work_dir, ignore_errors=True)

    print(f"\nTotal wall clock: {time.time() - overall_t0:.1f}s")
    print(f"Pattern outputs at: {args.output_dir}/dose_*/")


if __name__ == '__main__':
    main()
