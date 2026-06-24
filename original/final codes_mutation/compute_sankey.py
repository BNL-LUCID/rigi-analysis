#!/usr/bin/env python
"""
Compute per-dose Sankey trajectory JSONs from annotated mutations.

Unlike the original compute_sankey_flows.py, which emitted pairwise
transitions like "W1_Exposed->W2_Lost: 1234", this version emits full
4-week trajectories like:

    "W0_Absent->W1_Exposed->W2_Lost->W3_Exposed_Recurrent": 1234

This is what the trajectory-aware renderer needs to stack ribbons by full
forward path, eliminating the visual artifact where independently stacked
per-transition ribbons created false visual continuity (a ribbon appearing
to flow W1_Exposed -> ... -> W3_Exposed when actually no mutation followed
that path; it was two unrelated ribbons that happened to align spatially).

State taxonomy (matches manuscript Fig 2 rendering):

    W0: Present, Absent
    W1: Both, Lost, Control, Exposed                (no Recurrent -- W0 is baseline)
    W2: Both, Exposed_Recurrent, Control_Recurrent,
        Exposed, Control, Lost                       (Recurrent if ever-present-before
                                                      AND absent at immediately prior week)
    W3: same as W2, transitioning from W2

Outputs:
  - dose_<dose>/all_chromosomes_trajectories.json -- one per non-control dose
  - combined/all_chromosomes_trajectories.json    -- one for the full dataset

Backward-compatibility files (pairwise transitions, same format as before):
  - dose_<dose>/all_chromosomes_flows.json
  - combined/all_chromosomes_flows.json

Usage:
    python compute_sankey_flows.py \\
        --input      annotated_mutations/DBS/all_DBS_annotated.pkl \\
        --output-dir sankey_flows/DBS
"""

import argparse
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd


CONTROL_PATTERN = re.compile(r'(?:d0|D0|[Cc]ontrol)')


def compute_states_for_subset(df: pd.DataFrame) -> pd.DataFrame:
    """Return DataFrame indexed by PermanentMutationID with W0..W3 state columns.

    A state is "Recurrent" only if the mutation was present in that arm at
    some prior timepoint AND was absent in that arm at the immediately
    preceding timepoint. Mutations that were never present before are
    classified as Exposed/Control (newly appearing), not Recurrent.
    """
    if 'PermanentMutationID' not in df.columns:
        df = df.copy()
        df['PermanentMutationID'] = (
            df['Chromosome'].astype(str) + '_' +
            df['Start'].astype(str) + '_' +
            df['Ref'].astype(str) + '_' +
            df['Alt'].astype(str)
        )

    is_control = (
        df['Dose'].astype(str).str.match(CONTROL_PATTERN).fillna(False)
    )

    flag_df = pd.DataFrame({'PermanentMutationID': df['PermanentMutationID'].values})
    for week in range(4):
        in_week = (df['Timepoint'] == f'W{week}').values
        flag_df[f'control_w{week}'] = in_week & is_control.values
        if week >= 1:
            flag_df[f'treated_w{week}'] = in_week & ~is_control.values

    mut = flag_df.groupby('PermanentMutationID').max()

    c0 = mut['control_w0'].values
    c1, t1 = mut['control_w1'].values, mut['treated_w1'].values
    c2, t2 = mut['control_w2'].values, mut['treated_w2'].values
    c3, t3 = mut['control_w3'].values, mut['treated_w3'].values

    ever_treated_through_w0 = np.zeros(len(mut), dtype=bool)
    ever_control_through_w0 = c0.copy()

    ever_treated_through_w1 = ever_treated_through_w0 | t1
    ever_control_through_w1 = ever_control_through_w0 | c1

    ever_treated_through_w2 = ever_treated_through_w1 | t2
    ever_control_through_w2 = ever_control_through_w1 | c2

    mut['W0'] = np.where(c0, 'W0_Present', 'W0_Absent')

    mut['W1'] = np.select(
        [c1 & t1, c1 & ~t1, ~c1 & t1],
        ['W1_Both', 'W1_Control', 'W1_Exposed'],
        default='W1_Lost',
    )

    treated_recur_at_w2 = ever_treated_through_w1 & ~t1 & ~c2 & t2
    control_recur_at_w2 = ever_control_through_w1 & ~c1 & c2 & ~t2
    mut['W2'] = np.select(
        [
            c2 & t2,
            treated_recur_at_w2,
            control_recur_at_w2,
            ~c2 & t2,
            c2 & ~t2,
        ],
        [
            'W2_Both',
            'W2_Exposed_Recurrent',
            'W2_Control_Recurrent',
            'W2_Exposed',
            'W2_Control',
        ],
        default='W2_Lost',
    )

    treated_recur_at_w3 = ever_treated_through_w2 & ~t2 & ~c3 & t3
    control_recur_at_w3 = ever_control_through_w2 & ~c2 & c3 & ~t3
    mut['W3'] = np.select(
        [
            c3 & t3,
            treated_recur_at_w3,
            control_recur_at_w3,
            ~c3 & t3,
            c3 & ~t3,
        ],
        [
            'W3_Both',
            'W3_Exposed_Recurrent',
            'W3_Control_Recurrent',
            'W3_Exposed',
            'W3_Control',
        ],
        default='W3_Lost',
    )

    return mut[['W0', 'W1', 'W2', 'W3']]


def compute_trajectories(states: pd.DataFrame) -> dict:
    """Aggregate full 4-week paths into a flat dict.

    Keys look like 'W0_Absent->W1_Exposed->W2_Lost->W3_Exposed_Recurrent'.
    Each mutation contributes to exactly one trajectory.
    """
    paths = (
        states['W0'].astype(str) + '->' +
        states['W1'].astype(str) + '->' +
        states['W2'].astype(str) + '->' +
        states['W3'].astype(str)
    )
    counts = paths.value_counts()
    return {str(k): int(v) for k, v in counts.items()}


def trajectories_to_pairwise(trajectories: dict) -> dict:
    """Collapse trajectories to per-transition counts (legacy format).

    Useful for backward-compatibility with consumers expecting the old
    pairwise flows JSON layout.
    """
    flows = {}
    for traj, count in trajectories.items():
        states = traj.split('->')
        for i in range(len(states) - 1):
            key = f'{states[i]}->{states[i+1]}'
            flows[key] = flows.get(key, 0) + count
    return flows


def write_outputs(states: pd.DataFrame, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    trajectories = compute_trajectories(states)
    traj_path = os.path.join(out_dir, 'all_chromosomes_trajectories.json')
    with open(traj_path, 'w') as f:
        json.dump(trajectories, f, indent=2, sort_keys=True)

    flows = trajectories_to_pairwise(trajectories)
    flows_path = os.path.join(out_dir, 'all_chromosomes_flows.json')
    with open(flows_path, 'w') as f:
        json.dump(flows, f, indent=2, sort_keys=True)

    return trajectories, flows, traj_path, flows_path


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--input', '-i', required=True,
                   help='Annotated mutations pickle (Step 4 output, e.g. all_DBS_annotated.pkl)')
    p.add_argument('--output-dir', '-o', required=True,
                   help='Output dir; one dose_<dose>/ subdir written per dose plus combined/')
    args = p.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input not found: {args.input}")

    print(f"Loading {args.input} ...")
    t0 = time.time()
    df = pd.read_pickle(args.input)
    print(f"  {len(df):,} rows in {time.time() - t0:.1f}s")

    for col in ('Chromosome', 'Start', 'Ref', 'Alt', 'Dose', 'Timepoint'):
        if col not in df.columns:
            sys.exit(f"ERROR: input missing required column: {col}")

    all_doses = df['Dose'].dropna().unique().tolist()
    treated_doses = sorted([d for d in all_doses if not CONTROL_PATTERN.match(str(d))])
    control_doses = [d for d in all_doses if CONTROL_PATTERN.match(str(d))]

    if not treated_doses:
        sys.exit("ERROR: no treated doses found (CONTROL_PATTERN matched all)")
    if not control_doses:
        print("WARNING: no control doses detected -- W0 baseline will be empty")

    print(f"Controls: {control_doses}")
    print(f"Treated:  {treated_doses}")

    control_mask = df['Dose'].isin(control_doses)
    n_control_rows = int(control_mask.sum())

    os.makedirs(args.output_dir, exist_ok=True)

    for dose in treated_doses:
        print(f"\n=== Dose {dose} ===")
        t_dose = time.time()
        dose_mask = df['Dose'] == dose
        subset = df[control_mask | dose_mask]
        print(f"  rows: {len(subset):,} "
              f"(controls={n_control_rows:,}, treated={int(dose_mask.sum()):,})")

        states = compute_states_for_subset(subset)
        out_dir = os.path.join(args.output_dir, f'dose_{dose}')
        trajectories, flows, traj_path, flows_path = write_outputs(states, out_dir)

        print(f"  unique mutations: {len(states):,}")
        print(f"  unique trajectories: {len(trajectories)}")
        print(f"  unique pairwise flows: {len(flows)}")
        print(f"  saved: {traj_path}")
        print(f"  saved: {flows_path}")
        print(f"  elapsed: {time.time() - t_dose:.1f}s")

    print(f"\n=== Combined (all doses + controls) ===")
    t_combined = time.time()
    combined_states = compute_states_for_subset(df)
    out_dir = os.path.join(args.output_dir, 'combined')
    trajectories, flows, traj_path, flows_path = write_outputs(combined_states, out_dir)
    print(f"  unique mutations: {len(combined_states):,}")
    print(f"  unique trajectories: {len(trajectories)}")
    print(f"  unique pairwise flows: {len(flows)}")
    print(f"  saved: {traj_path}")
    print(f"  saved: {flows_path}")
    print(f"  elapsed: {time.time() - t_combined:.1f}s")

    print(f"\nAll done in {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()