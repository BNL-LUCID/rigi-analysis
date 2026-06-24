#!/usr/bin/env python
r"""Trajectory-aware Sankey renderer.

Unlike the previous renderer, which drew per-transition ribbons with
independent stacking at each node (causing visual artifacts where ribbons
appeared to flow continuously through paths no mutation actually followed),
this version preserves trajectory identity end-to-end.

Each mutation is assigned to one of N unique trajectories (full 4-week
state paths). For each trajectory, all four ribbon segments
(W0->W1, W1->W2, W2->W3) are drawn as a contiguous unit -- the segment at
each transition takes its slot in the source/target node based on the
trajectory's overall layout order, so the ribbon visually traces the
mutation's actual path through all four columns.

Input: a trajectories JSON from compute_sankey_flows_trajectory.py:
    {"W0_Absent->W1_Exposed->W2_Lost->W3_Exposed_Recurrent": 4744, ...}

Usage:
    python render_sankey_trajectory.py \\
        --trajectories-json sankey_flows/DBS/combined/all_chromosomes_trajectories.json \\
        --output            sankey_figures/DBS_combined.png \\
        --title             "Temporal Dynamics for DBS" \\
        --subtitle          "All doses, all chromosomes"
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path

COLUMNS = {
    0: ['W0_Present', 'W0_Absent'],
    1: ['W1_Both', 'W1_Lost', 'W1_Control', 'W1_Exposed'],
    2: ['W2_Both', 'W2_Exposed_Recurrent', 'W2_Control_Recurrent',
        'W2_Exposed', 'W2_Control', 'W2_Lost'],
    3: ['W3_Both', 'W3_Exposed_Recurrent', 'W3_Control_Recurrent',
        'W3_Lost', 'W3_Exposed', 'W3_Control'],
}

NODE_LAYOUT_INDEX = {
    node: (col_idx, pos_idx)
    for col_idx, nodes in COLUMNS.items()
    for pos_idx, node in enumerate(nodes)
}

COLOR_MAP = {
    'Present':           '#8B7355',
    'Absent':            '#A0937D',
    'Both':              '#9890B0',
    'Exposed':           '#C4A5A0',
    'Control':           '#8BA5B5',
    'Lost':              '#B08080',
    'Exposed_Recurrent': '#7DAA9E',
    'Control_Recurrent': '#7DAA9E',
}

LEGEND_ITEMS = [
    ('Both',             '#9890B0'),
    ('Exposed',          '#C4A5A0'),
    ('Control',          '#8BA5B5'),
    ('Lost',             '#B08080'),
    ('Recurrent',        '#7DAA9E'),
    ('Baseline Present', '#8B7355'),
    ('Baseline Absent',  '#A0937D'),
]


def get_node_color(node_name: str) -> str:
    for key in ('Exposed_Recurrent', 'Control_Recurrent',
                'Present', 'Absent', 'Both', 'Exposed', 'Control', 'Lost'):
        if key in node_name:
            return COLOR_MAP[key]
    return '#AAAAAA'


def get_link_color(node_name: str, alpha: float = 0.35) -> tuple:
    base = get_node_color(node_name)
    r = int(base[1:3], 16)
    g = int(base[3:5], 16)
    b = int(base[5:7], 16)
    return (r / 255, g / 255, b / 255, alpha)


def compute_node_sizes_from_trajectories(trajectories: dict) -> dict:
    """For each node, return max(inflow, outflow) summed across trajectories."""
    inflow = defaultdict(int)
    outflow = defaultdict(int)
    for traj, count in trajectories.items():
        states = traj.split('->')
        for i in range(len(states) - 1):
            outflow[states[i]] += count
            inflow[states[i + 1]] += count
    sizes = {}
    for n in set(inflow) | set(outflow):
        sizes[n] = max(inflow[n], outflow[n])
    return sizes


def trajectory_layout_key(traj: str) -> tuple:
    """Sort trajectories by layout position at each week.

    This determines the overall order in which trajectories are stacked.
    Sorting top-down means trajectories ending highest (W3) sit on top,
    then secondary-sorted by W2, W1, W0.
    """
    states = traj.split('->')
    return tuple(NODE_LAYOUT_INDEX.get(s, (99, 99))[1] for s in states)


def draw_sankey_on_axis(ax, trajectories: dict, title: str = '', subtitle: str = '',
                        font_scale: float = 1.0, draw_legend: bool = True):
    """Render the trajectory-aware sankey onto a matplotlib axis."""
    node_sizes = compute_node_sizes_from_trajectories(trajectories)

    ax.set_xlim(-0.1, 5.0)
    ax.set_ylim(-0.22, 1.05)
    ax.axis('off')

    col_x = {0: 0.3, 1: 1.3, 2: 2.5, 3: 3.7}
    node_width = 0.08
    total_height = 0.85
    y_padding = 0.02

    # Lay out nodes
    node_positions = {}
    for col_idx, nodes in COLUMNS.items():
        x = col_x[col_idx]
        col_total = sum(node_sizes.get(n, 0) for n in nodes if node_sizes.get(n, 0) > 0)
        if col_total == 0:
            continue
        y_current = 0.92
        for node in nodes:
            size = node_sizes.get(node, 0)
            if size == 0:
                continue
            height = (size / col_total) * total_height
            node_positions[node] = (x, y_current - height, height)
            y_current -= height + y_padding

    # Draw nodes
    for node, (x, y, h) in node_positions.items():
        color = get_node_color(node)
        rect = mpatches.FancyBboxPatch(
            (x - node_width / 2, y), node_width, h,
            boxstyle='round,pad=0.01,rounding_size=0.01',
            facecolor=color, edgecolor='#666666', linewidth=0.8,
        )
        ax.add_patch(rect)

    # Pre-sort trajectories by layout key. Drawing in this order makes
    # ribbons stack consistently within each node: a trajectory passing
    # through W2_Lost on its way to W3_Both will sit at the same relative
    # position (top) at all transitions touching that path.
    sorted_trajs = sorted(trajectories.items(), key=lambda kv: trajectory_layout_key(kv[0]))

    # Track the current "next ribbon top" at each node, separately for
    # incoming and outgoing sides. Because we draw trajectories in a fixed
    # order, every ribbon segment of trajectory T finds its slot in the
    # source node's outgoing stack and the target node's incoming stack
    # consistently.
    node_out_y = {n: pos[1] + pos[2] for n, pos in node_positions.items()}
    node_in_y = {n: pos[1] + pos[2] for n, pos in node_positions.items()}

    for traj, count in sorted_trajs:
        states = traj.split('->')
        for i in range(len(states) - 1):
            src, tgt = states[i], states[i + 1]
            if src not in node_positions or tgt not in node_positions:
                continue

            src_x, src_y, src_h = node_positions[src]
            tgt_x, tgt_y, tgt_h = node_positions[tgt]

            flow_h_src = (count / node_sizes[src]) * src_h
            flow_h_tgt = (count / node_sizes[tgt]) * tgt_h

            y1_top = node_out_y[src]
            y1_bottom = y1_top - flow_h_src
            y2_top = node_in_y[tgt]
            y2_bottom = y2_top - flow_h_tgt

            node_out_y[src] = y1_bottom
            node_in_y[tgt] = y2_bottom

            x1 = src_x + node_width / 2
            x2 = tgt_x - node_width / 2
            ctrl = (x2 - x1) * 0.4

            verts = [
                (x1, y1_top),
                (x1 + ctrl, y1_top), (x2 - ctrl, y2_top), (x2, y2_top),
                (x2, y2_bottom),
                (x2 - ctrl, y2_bottom), (x1 + ctrl, y1_bottom), (x1, y1_bottom),
                (x1, y1_top),
            ]
            codes = [
                Path.MOVETO,
                Path.CURVE4, Path.CURVE4, Path.CURVE4,
                Path.LINETO,
                Path.CURVE4, Path.CURVE4, Path.CURVE4,
                Path.CLOSEPOLY,
            ]
            path = Path(verts, codes)
            patch = mpatches.PathPatch(
                path, facecolor=get_link_color(src), edgecolor='none', zorder=1,
            )
            ax.add_patch(patch)

    # Node labels
    for node, (x, y, h) in node_positions.items():
        parts = node.split('_', 1)
        state = parts[1].replace('_', ' ') if len(parts) > 1 else ''
        if x < 2:
            label_x = x - node_width / 2 - 0.03
            ha = 'right'
        else:
            label_x = x + node_width / 2 + 0.03
            ha = 'left'
        ax.text(label_x, y + h / 2, state,
                fontsize=14 * font_scale, fontweight='medium',
                ha=ha, va='center', color='#333333')

    # Week labels
    week_labels = ['Week 0\n(Baseline)', 'Week 1', 'Week 2', 'Week 3']
    for col_idx, label in enumerate(week_labels):
        ax.text(col_x[col_idx], -0.05, label,
                fontsize=18 * font_scale, fontweight='bold',
                ha='center', va='top', color='#444444')

    # Title + subtitle
    if title:
        ax.text(2.2, 1.04, title,
                fontsize=28 * font_scale, fontweight='bold',
                ha='center', va='top', color='#333333')
    if subtitle:
        ax.text(2.2, 0.98, subtitle,
                fontsize=18 * font_scale,
                ha='center', va='top', color='#555555')

    if draw_legend:
        legend_y = -0.14
        legend_start_x = 0.2
        legend_spacing = 0.52
        for i, (label, color) in enumerate(LEGEND_ITEMS):
            x_pos = legend_start_x + i * legend_spacing
            rect = mpatches.Rectangle(
                (x_pos, legend_y), 0.05, 0.028,
                facecolor=color, edgecolor='#666666', linewidth=0.5,
            )
            ax.add_patch(rect)
            ax.text(x_pos + 0.06, legend_y + 0.014, label,
                    fontsize=14 * font_scale, va='center', color='#444444')


def render(trajectories: dict, title: str, subtitle: str, output: str, dpi: int = 300):
    """Single-panel render to a standalone PNG."""
    fig, ax = plt.subplots(figsize=(19.2, 10.8), dpi=150)
    draw_sankey_on_axis(ax, trajectories, title=title, subtitle=subtitle,
                        font_scale=1.0, draw_legend=True)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
    plt.savefig(output, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved: {output}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--trajectories-json', '-f', required=True,
                   help='Path to trajectories JSON from compute_sankey_flows_trajectory.py')
    p.add_argument('--output', '-o', required=True,
                   help='Output PNG path')
    p.add_argument('--title', '-t', default='Temporal Mutation Dynamics',
                   help='Main title')
    p.add_argument('--subtitle', '-s', default='',
                   help='Subtitle (e.g. "Dose dA, All Chromosomes")')
    p.add_argument('--dpi', type=int, default=300)
    args = p.parse_args()

    if not os.path.exists(args.trajectories_json):
        sys.exit(f"ERROR: trajectories JSON not found: {args.trajectories_json}")

    with open(args.trajectories_json) as f:
        trajectories = json.load(f)
    print(f"Loaded {len(trajectories)} trajectories from {args.trajectories_json}")
    print(f"  total mutations: {sum(trajectories.values()):,}")

    render(trajectories, args.title, args.subtitle, args.output, args.dpi)


if __name__ == '__main__':
    main()
