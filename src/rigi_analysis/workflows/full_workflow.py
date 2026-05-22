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
        description='Workflow application for All RIGI Pipelines',
        usage='rigi-analysis-workflow [options]',
    )
    parser.add_argument(
        '-c', '--config-file', dest='config_file', type=str, required=True,
        help='configuration file with the workflow description',
    )
    parser.add_argument(
        '-o', '--output-dir', dest='output_dir', type=str, required=False,
        help='directory path for workflow output',
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

    init_default_logger()

    backend = await get_backend(cfg.get('run_description', {}))
    flow = await WorkflowEngine.create(backend=backend)

    # General output dir for full workflow
    out_dir = args.output_dir or cfg.get('output_dir', CURRENT_DIR)
    os.makedirs(out_dir, exist_ok=True)

    out_mutation_dir = os.path.join(out_dir, 'out_mutation')
    out_sv_dir = os.path.join(out_dir, 'out_sv')

    try:
        # Build mutation pipeline from original config
        mut_pipe = MutationPipeline(
            name='full-mutation-pipe', 
            config=cfg, 
            flow=flow,
            output_dir=out_mutation_dir,
        )

        # Deep-copy config for SV pipeline so we don't corrupt the
        # original; then wire the mutation merged dir into it
        sv_cfg = copy.deepcopy(cfg)
        sv_pipe = SVPipeline(
            name='full-sv-pipe', 
            config=sv_cfg, 
            flow=flow,
            output_dir=out_sv_dir,
            mutation_output_dir=out_mutation_dir,
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
