"""Full RIGI Analysis Workflow.

Executes both the Mutation and SV pipelines sequentially
using RADICAL tools for distributed task execution.
"""

import argparse
import asyncio
import copy
import json
import os
import sys

from radical.asyncflow import WorkflowEngine
from radical.asyncflow.logging import init_default_logger

from ..utils.datetime_utils import datetime_now

from .backend import get_backend
from .mutation_workflow import MutationPipeline
from .sv_workflow import SVPipeline

CURRENT_DIR = os.getcwd()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Workflow application for Full RIGI Pipeline',
        usage='rigi-analysis-run full_workflow [options]',
    )
    parser.add_argument(
        '-w', '--work-dir', dest='work_dir', type=str, required=False,
        help='workspace for workflow session sandboxes',
    )
    parser.add_argument(
        '-o', '--output-dir', dest='output_dir', type=str, required=False,
        help='directory path for output',
    )
    parser.add_argument(
        '-c', '--config-file', dest='config_file', type=str, required=True,
        help='configuration file with the workflow description',
    )
    parser.add_argument(
        '-t', '--runtime', dest='runtime', type=int, required=False,
        help='requested runtime (min) for the workflow to run',
    )
    return parser.parse_args(sys.argv[1:])


async def run_workflow() -> None:
    """Main asynchronous execution function.

    Runs the mutation pipeline to completion first (producing merged
    CSVs), then feeds those outputs into the SV pipeline.
    """
    args = parse_args()

    config_file = os.path.abspath(args.config_file)
    with open(config_file, 'r') as f:
        cfg = json.load(f)

    cfg['cfg_file'] = config_file

    init_default_logger()

    sandbox_path = os.path.abspath(args.work_dir or CURRENT_DIR)
    run_description = dict(
        **cfg.get('run_description', {}), sandbox=sandbox_path)
    if args.runtime:
        run_description['runtime'] = int(args.runtime)

    backend = await get_backend(run_description)
    flow = await WorkflowEngine.create(backend=backend)

    # General output dir for full workflow
    out_dir = (
        args.output_dir
        or cfg.get(
            'output_dir', os.path.join(CURRENT_DIR, 'out_full'))
    )
    os.makedirs(out_dir, exist_ok=True)

    mut_out = os.path.join(out_dir, 'mutation_output')
    sv_out = os.path.join(out_dir, 'sv_output')

    try:
        # Build mutation pipeline from original config
        mut_pipe = MutationPipeline(
            name='full-mutation-pipe', config=cfg, flow=flow,
            output_dir=mut_out,
        )

        # Deep-copy config for SV pipeline so we don't corrupt the
        # original; then wire the mutation merged dir into it
        sv_cfg = copy.deepcopy(cfg)
        sv_cfg['mutation_merged_dir'] = os.path.join(
            mut_out, 'merged_analysis')
        sv_pipe = SVPipeline(
            name='full-sv-pipe', config=sv_cfg, flow=flow,
            output_dir=sv_out,
        )

        print(f'{datetime_now()} Full Pipeline started')

        # Mutation pipeline must finish first — SV consumes its outputs
        mut_results = await mut_pipe.run()
        print(f'{datetime_now()} Mutation Pipeline finished')

        sv_results = await sv_pipe.run()
        print(f'{datetime_now()} SV Pipeline finished')

        print(f'Full pipeline results: {(mut_results, sv_results)}')
    except Exception:
        raise
    finally:
        await flow.shutdown()


def main() -> None:
    """Entry point."""
    asyncio.run(run_workflow())


if __name__ == '__main__':
    main()
