#!/usr/bin/env python3
"""Repeat Element Analysis - 4-Panel Streamlined Version (v2: any-breakpoint semantics).
====================================================================================
Clean supplementary figure showing:
A. Repeat involvement (pie)         — Both / One / No (unchanged from v1)
B. Repeat categories (bar)          — unchanged from v1
C. Repeat by SV type (bar)          — uses ANY breakpoint (v1 used left only)
D. Control vs Radiation (bar)       — uses ANY breakpoint (v1 used left only)

Difference from v1 (`04_repeat_analysis.py`):
  v1 panels C and D counted only the LEFT breakpoint, systematically
  undercounting SVs whose right breakpoint sits in a repeat. v2 uses
  `Has_Any_Repeat = Has_Left_Repeat | Has_Right_Repeat`, matching the
  manuscript's "either breakpoint in a repeat" convention.

Run both versions on the same input to inspect the difference.
"""

import argparse
import glob

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Professional color scheme (matching Figure 3)
COLORS = {
    'INV': '#2E5A88',
    'TRA': '#5B8DBE',
    'DEL': '#8C8C8C',
    'DUP': '#404040',
    'INS': '#D4D4D4',
    'Control': '#7FA6C9',
    'Radiation': '#D4896A'
}

REPEAT_COLORS = {
    'None': '#B8B8B8',
    'Alu': '#8B6BA8',
    'LINE': '#E88B6F',
    'SINE': '#F0B67F',
    'LTR': '#A8C98F',
    'DNA': '#8FACC9',
    'Other': '#C9B8A8'
}

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = [
    'Arial', 'Liberation Sans', 'DejaVu Sans', 'sans-serif']
plt.rcParams['font.size'] = 11


def categorize_repeat(repeat_type):
    """Categorize repeat into major classes.

    Pattern matching expanded (vs the original v2) to capture variant naming
    used by AnnotSV/RepeatMasker — L3 (LINE), MLT*/THE1* (LTR/MaLR),
    Tigger*/MER* (DNA), FLAM* (Alu ancestor / SINE).

    A new 'Simple' category captures tandem/low-complexity/satellite repeats
    (ALR/Alpha, A-rich, (CA)n, satellites). These are reported separately so
    they can be excluded from "interspersed-repeat" counting in panels C/D.
    """
    if pd.isna(repeat_type) or repeat_type == '' or repeat_type == '.':
        return 'None'

    s = str(repeat_type).lower().strip()

    # ---- Simple / satellite / low-complexity ---------------------------------
    # Counted separately; excluded from Has_Any_Repeat (interspersed) downstream.
    if (
        s.startswith('(')                         # tandem like (CA)n, (T)n
        or 'simple' in s
        or 'low_complexity' in s or 'low complex' in s
        or 'satellite' in s
        or s.startswith('alr/') or 'alpha' in s   # alpha satellite (ALR/Alpha)
        or s.startswith('hsat')                   # human satellite
        or 'a-rich' in s or 't-rich' in s
        or 'g-rich' in s or 'c-rich' in s
    ):
        return 'Simple'

    # ---- Interspersed repeat families ---------------------------------------
    # Alu (SINE family) — includes FLAM_C, the Alu ancestor.
    if 'alu' in s or s.startswith('flam'):
        return 'Alu'

    # LINE family
    if 'line' in s or s.startswith('l1') or s.startswith('l2') or s.startswith('l3'):
        return 'LINE'

    # SINE family (non-Alu)
    if 'sine' in s or 'mir' in s:
        return 'SINE'

    # LTR retrotransposons (ERV; MaLR family includes MLT*, THE1*, MST*)
    if (
        'ltr' in s or 'erv' in s or 'herv' in s
        or s.startswith('mlt') or s.startswith('the1') or s.startswith('mst')
    ):
        return 'LTR'

    # DNA transposons (hAT, Tigger, MER, Charlie, TcMar)
    if (
        'dna' in s or 'hat' in s or 'charlie' in s or 'tcmar' in s
        or s.startswith('tigger') or s.startswith('mer')
    ):
        return 'DNA'

    return 'Other'


def create_sv_key(row):
    """Create unique SV key for deduplication."""
    chrom = row.get('SV_chrom', 'NA')

    # Handle NaN values in coordinates
    start_val = row.get('SV_start', 0)
    end_val = row.get('SV_end', 0)

    # Convert to int, handling NaN
    try:
        start = int(start_val) if pd.notna(start_val) else 0
        end = int(end_val) if pd.notna(end_val) else 0
    except (ValueError, TypeError):
        start = 0
        end = 0

    sv_type = row.get('SV_type', 'NA')

    # Sort coordinates for inversions
    coord_key = f"{min(start,end)}-{max(start,end)}"
    return f"{chrom}:{coord_key}:{sv_type}"


def load_and_process_svs(annotsv_dir, label='Data'):
    """Load AnnotSV files and deduplicate."""
    print(f"\nLoading {label}...")

    files = glob.glob(f"{annotsv_dir}/*.tsv")
    if not files:
        print(f"  No files found in {annotsv_dir}")
        return None

    all_svs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', low_memory=False)
            if 'Annotation_mode' in df.columns:
                df = df[df['Annotation_mode'] == 'full']
            all_svs.append(df)
        except Exception as e:
            print(f"  Error loading {f}: {e}")
            continue

    if not all_svs:
        return None

    combined = pd.concat(all_svs, ignore_index=True)
    print(f"  Loaded {len(combined):,} entries")

    # Deduplicate
    combined['SV_key'] = combined.apply(create_sv_key, axis=1)
    deduped = combined.drop_duplicates(subset=['SV_key'])
    print(f"  After deduplication: {len(deduped):,} unique SVs")

    return deduped


def analyze_repeat_involvement(df):
    """Analyze repeat element involvement."""
    # Find repeat columns
    repeat_cols = [c for c in df.columns if 'repeat' in c.lower()]

    # Try to identify left and right
    left_col = None
    right_col = None

    for col in repeat_cols:
        if 'left' in col.lower():
            left_col = col
        elif 'right' in col.lower():
            right_col = col

    # If not found by left/right, try first two repeat columns
    if not left_col and len(repeat_cols) >= 2:
        left_col = repeat_cols[0]
        right_col = repeat_cols[1]

    if not left_col or not right_col:
        print("  Warning: Could not find repeat columns")
        print(f"  Available columns with 'repeat': {repeat_cols}")
        return None

    print(f"  Using columns: {left_col}, {right_col}")

    # Categorize repeats
    df['Left_Repeat_Cat'] = df[left_col].apply(categorize_repeat)
    df['Right_Repeat_Cat'] = df[right_col].apply(categorize_repeat)

    # Stash the original column names so main() can write the categorization CSV.
    df.attrs['_repeat_left_col'] = left_col
    df.attrs['_repeat_right_col'] = right_col

    # ---- Has_*_Repeat: interspersed-repeat semantics ------------------------
    # We count a breakpoint as "in repeat" only when it sits in an interspersed
    # repeat family (Alu, LINE, SINE, LTR, DNA) or an unclassified 'Other'
    # interspersed element. Tandem/simple/satellite repeats ('Simple') are
    # NOT counted — they don't mediate the SV mechanisms (MMEJ, HR) the
    # manuscript discusses, and including them inflates the percentage above
    # the published values.
    INTERSPERSED = {'Alu', 'LINE', 'SINE', 'LTR', 'DNA', 'Other'}
    df['Has_Left_Repeat']  = df['Left_Repeat_Cat'].isin(INTERSPERSED)
    df['Has_Right_Repeat'] = df['Right_Repeat_Cat'].isin(INTERSPERSED)
    df['Has_Any_Repeat']   = df['Has_Left_Repeat'] | df['Has_Right_Repeat']

    # ---- One-line per-side summary (full breakdown goes to CSV) -------------
    for side, cat_col in [('LEFT', 'Left_Repeat_Cat'),
                          ('RIGHT', 'Right_Repeat_Cat')]:
        vc = df[cat_col].value_counts()
        n_simple = int(vc.get('Simple', 0))
        n_other  = int(vc.get('Other', 0))
        pct_simple = 100 * n_simple / len(df) if len(df) else 0
        pct_other  = 100 * n_other  / len(df) if len(df) else 0
        print(f"  {side}: Simple={n_simple:,} ({pct_simple:.1f}%) | "
              f"Other={n_other:,} ({pct_other:.1f}%) | "
              f"see categorization CSV for full breakdown")

    def classify_involvement(row):
        if row['Has_Left_Repeat'] and row['Has_Right_Repeat']:
            return 'Both in Repeat'
        elif row['Has_Left_Repeat'] or row['Has_Right_Repeat']:
            return 'One in Repeat'
        else:
            return 'No Repeat'

    df['Repeat_Involvement'] = df.apply(classify_involvement, axis=1)

    return df


def write_categorization_csv(datasets, output_csv):
    """Write a long-form CSV capturing every raw repeat-type value seen in any
    input file, the category it was assigned to, and how often it occurred.

    `datasets` is a list of `(source_label, df)` tuples — typically
    `[('Radiation', radiation_df), ('Control', control_df)]`. `df` must have
    been processed by `analyze_repeat_involvement()` so the original
    LEFT/RIGHT repeat column names are available on `df.attrs`.

    Output schema (one row per Source × Side × Raw_Repeat_Type):
        Source, Side, Assigned_Category, Raw_Repeat_Type, Count
    Sorted to put Simple/Other first (so the rows that drove categorization
    decisions are at the top) and then largest counts within each category.
    """
    rows = []
    for source_label, df in datasets:
        if df is None:
            continue
        left_col = df.attrs.get('_repeat_left_col')
        right_col = df.attrs.get('_repeat_right_col')
        for side, raw_col, cat_col in [('LEFT',  left_col,  'Left_Repeat_Cat'),
                                        ('RIGHT', right_col, 'Right_Repeat_Cat')]:
            if raw_col is None:
                continue
            grp = (
                df.groupby([raw_col, cat_col], dropna=False)
                  .size()
                  .reset_index(name='Count')
            )
            grp = grp.rename(columns={raw_col: 'Raw_Repeat_Type',
                                       cat_col: 'Assigned_Category'})
            grp['Source'] = source_label
            grp['Side'] = side
            rows.append(grp)

    if not rows:
        print(f"  WARN: no data to write to {output_csv}")
        return

    out = pd.concat(rows, ignore_index=True)
    out = out[['Source', 'Side', 'Assigned_Category', 'Raw_Repeat_Type', 'Count']]

    # Sort: Source A→Z, Side L→R, Simple/Other first, then Count desc within group
    cat_order = {'Simple': 0, 'Other': 1, 'None': 2,
                 'Alu': 3, 'LINE': 4, 'SINE': 5, 'LTR': 6, 'DNA': 7}
    out['_cat_rank'] = out['Assigned_Category'].map(cat_order).fillna(99)
    out = out.sort_values(
        ['Source', 'Side', '_cat_rank', 'Count'],
        ascending=[True, True, True, False]
    ).drop(columns='_cat_rank')

    out.to_csv(output_csv, index=False)
    print(f"  Saved categorization breakdown: {output_csv}")


def create_4panel_figure(radiation_df, control_df, output):
    """Create streamlined 4-panel figure."""
    print("\nCreating 4-panel figure...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # =========================================================================
    # PANEL A: Repeat Involvement (Pie Chart)
    # =========================================================================
    ax = axes[0, 0]

    involvement = radiation_df['Repeat_Involvement'].value_counts()
    colors_pie = [REPEAT_COLORS['Other'], REPEAT_COLORS['SINE'], REPEAT_COLORS['None']]

    wedges, texts, autotexts = ax.pie(
        involvement.values,
        labels=involvement.index,
        autopct='%1.1f%%',
        colors=colors_pie,
        startangle=90,
        textprops={'fontsize': 11, 'fontweight': 'bold'}
    )

    # Make percentage text white for better contrast
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(12)
        autotext.set_fontweight('bold')

    ax.set_title('A. Repeat Involvement', fontweight='bold', fontsize=13, pad=15)

    # =========================================================================
    # PANEL B: Repeat Categories (Bar Chart)
    # =========================================================================
    ax = axes[0, 1]

    all_repeats = pd.concat([
        radiation_df['Left_Repeat_Cat'],
        radiation_df['Right_Repeat_Cat']
    ])
    repeat_counts = all_repeats[all_repeats != 'None'].value_counts().head(6)

    colors_bar = [REPEAT_COLORS.get(r, '#888') for r in repeat_counts.index]
    bars = ax.barh(range(len(repeat_counts)), repeat_counts.values,
                   color=colors_bar, edgecolor='black', linewidth=1, alpha=0.85)

    ax.set_yticks(range(len(repeat_counts)))
    ax.set_yticklabels(repeat_counts.index, fontsize=11, fontweight='bold')
    ax.set_xlabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('B. Repeat Category', fontweight='bold', fontsize=13, pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, repeat_counts.values)):
        ax.text(val * 1.02, i, f'{int(val):,}',
               va='center', ha='left', fontsize=10, fontweight='bold')

    # =========================================================================
    # PANEL C: Repeat by SV Type (Horizontal Bar)
    # =========================================================================
    ax = axes[1, 0]

    # v2: count any-breakpoint involvement instead of left-only
    sv_repeat = radiation_df.groupby('SV_type')['Has_Any_Repeat'].apply(
        lambda x: (x.sum() / len(x) * 100) if len(x) > 0 else 0
    ).sort_values(ascending=True)

    colors_sv = [COLORS.get(sv, '#888') for sv in sv_repeat.index]
    bars = ax.barh(range(len(sv_repeat)), sv_repeat.values,
                   color=colors_sv, edgecolor='black', linewidth=1, alpha=0.85)

    ax.set_yticks(range(len(sv_repeat)))
    ax.set_yticklabels(sv_repeat.index, fontsize=11, fontweight='bold')
    ax.set_xlabel('% with Repeat Element', fontsize=11, fontweight='bold')
    ax.set_title('C. Repeat by SV Type', fontweight='bold', fontsize=13, pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 110)

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, sv_repeat.values)):
        ax.text(val + 1, i, f'{val:.1f}%',
               va='center', ha='left', fontsize=10, fontweight='bold')

    # =========================================================================
    # PANEL D: Control vs Radiation (Bar Chart with Delta)
    # =========================================================================
    ax = axes[1, 1]

    if control_df is not None:
        # v2: any-breakpoint involvement (not left-only)
        rad_pct = (radiation_df['Has_Any_Repeat'].sum() / len(radiation_df)) * 100
        ctrl_pct = (control_df['Has_Any_Repeat'].sum() / len(control_df)) * 100

        bars = ax.bar(['Control\n(d0)', 'Radiation'],
                     [ctrl_pct, rad_pct],
                     color=[COLORS['Control'], COLORS['Radiation']],
                     edgecolor='black', linewidth=1.5, alpha=0.85,
                     width=0.6)

        ax.set_ylabel('% SVs with Repeat Element', fontsize=11, fontweight='bold')
        ax.set_title('D. Control vs Radiation', fontweight='bold', fontsize=13, pad=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(0, 100)

        # Add delta annotation
        delta = rad_pct - ctrl_pct
        y_max = max(ctrl_pct, rad_pct)
        ax.plot([0, 1], [y_max+8, y_max+8], 'k-', linewidth=2)
        ax.text(0.5, y_max+12, f'Δ = {delta:+.1f}%',
               ha='center', fontsize=11, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                        alpha=0.4, edgecolor='black', linewidth=1))

        # Add values on bars
        for bar, val in zip(bars, [ctrl_pct, rad_pct]):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height/2,
                   f'{val:.1f}%', ha='center', va='center',
                   fontsize=12, fontweight='bold', color='white')
    else:
        ax.text(0.5, 0.5, 'Control data\nnot available',
               ha='center', va='center', fontsize=12, style='italic',
               transform=ax.transAxes)
        ax.set_title('D. Control vs Radiation', fontweight='bold', fontsize=13, pad=15)
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output}")

    # Print summary stats
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print("\nRadiation samples:")
    print(f"  Total SVs: {len(radiation_df):,}")
    involvement = radiation_df['Repeat_Involvement'].value_counts()
    for cat, count in involvement.items():
        pct = 100 * count / len(radiation_df)
        print(f"  {cat}: {count:,} ({pct:.1f}%)")

    if control_df is not None:
        print("\nControl samples:")
        print(f"  Total SVs: {len(control_df):,}")
        ctrl_involvement = control_df['Repeat_Involvement'].value_counts()
        for cat, count in ctrl_involvement.items():
            pct = 100 * count / len(control_df)
            print(f"  {cat}: {count:,} ({pct:.1f}%)")


def resolve_output_path(output_arg: str, default_filename: str = 'figure_s2_repeat_analysis_v2.png') -> str:
    """Accept either a file path (with .png/.pdf/.svg extension) or a directory.

    - If output_arg points at an existing directory, or ends with a path separator,
      or has no extension, treat it as a directory and append `default_filename`.
    - Otherwise treat it as a file path.
    - In either case, ensure the parent directory exists.
    """
    import os

    looks_like_dir = (
        output_arg.endswith(os.sep)
        or output_arg.endswith('/')
        or os.path.isdir(output_arg)
        or os.path.splitext(output_arg)[1] == ''
    )
    if looks_like_dir:
        out_path = os.path.join(output_arg, default_filename)
    else:
        out_path = output_arg

    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    return out_path


def main():
    parser = argparse.ArgumentParser(
        description='Repeat Element Analysis - Streamlined 4-Panel Figure'
    )
    parser.add_argument('--radiation-dir', required=True,
                       help='Directory with radiation AnnotSV files')
    parser.add_argument('--control-dir', required=False,
                       help='Directory with control AnnotSV files')
    parser.add_argument('--output', default='figure_s2_repeat_analysis_v2.png',
                       help=('Output figure path. Accepts either a file path '
                             '(*.png/.pdf/.svg) or a directory — if a directory '
                             'is given, the figure is saved as '
                             'figure_s2_repeat_analysis_v2.png inside it. '
                             'Parent directories are created if missing.'))
    parser.add_argument('--categories-csv', default=None,
                       help=('Output path for the per-raw-type categorization '
                             'breakdown CSV. Defaults to a sibling file named '
                             '<output_basename>_categorization.csv.'))
    args = parser.parse_args()

    output_path = resolve_output_path(args.output)

    print("="*80)
    print("REPEAT ELEMENT ANALYSIS - 4-PANEL VERSION")
    print("="*80)
    print(f"Output: {output_path}")

    # Load radiation data
    radiation_df = load_and_process_svs(args.radiation_dir, 'Radiation')
    if radiation_df is None:
        print("ERROR: Could not load radiation data")
        return

    # Analyze repeats
    radiation_df = analyze_repeat_involvement(radiation_df)
    if radiation_df is None:
        print("ERROR: Could not analyze repeat involvement")
        return

    # Load control if provided
    control_df = None
    if args.control_dir:
        control_df = load_and_process_svs(args.control_dir, 'Control')
        if control_df is not None:
            control_df = analyze_repeat_involvement(control_df)

    # Write categorization breakdown CSV alongside the figure (sibling file).
    import os
    cat_csv = args.categories_csv or (
        os.path.splitext(output_path)[0] + '_categorization.csv'
    )
    write_categorization_csv(
        [('Radiation', radiation_df), ('Control', control_df)],
        cat_csv,
    )

    # Create figure
    create_4panel_figure(radiation_df, control_df, output_path)

    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
