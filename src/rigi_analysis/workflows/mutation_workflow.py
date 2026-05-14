"""Mutation Analysis Workflow.

Executes a multi-stage Mutation Analysis pipeline
using RADICAL tools for distributed task execution.
"""

import argparse
import asyncio
import json
import os
import sys
from typing import Optional, Any

from radical.asyncflow import WorkflowEngine
from radical.asyncflow.logging import init_default_logger

from ..utils.datetime_utils import datetime_now

from .backend import get_backend

CURRENT_DIR = os.getcwd()


class MutationPipeline:
    """Mutation Pipeline Execution Class.

    Orchestrates the seven-stage mutation analysis pipeline described
    in ``docs/mutation/README.md``.  Each stage is registered as a
    RADICAL asyncflow executable task; the ``run`` method wires them
    into the correct dependency DAG.
    """

    def __init__(
            self,
            name: str,
            config: dict,
            flow: WorkflowEngine,
            output_dir: Optional[str] = None,
    ):
        """Initialize the mutation pipeline.

        Args:
            name: Unique name for this pipeline instance.
            config: Full configuration dictionary with pipeline
                parameters as top-level keys.
            flow: WorkflowEngine instance from RADICAL tools.
            output_dir: Directory to write results to (overrides config).
        """
        self.flow = flow
        self.name = name

        self.cfg = config
        self.cfg_file = config.get('cfg_file', '')

        self.output_dir = (
            output_dir
            or self.cfg.get('output_dir',
                            os.path.join(CURRENT_DIR, 'out_mutation'))
        )
        os.makedirs(self.output_dir, exist_ok=True)

        # Define directories
        self.vcf_dir = self.cfg.get('vcf_dir', 'vcf_files')
        self.annotations_dir = self.cfg.get(
            'annotations_dir',
            os.path.join(self.output_dir, 'annotations'),
        )
        self.sigprofiler_out = os.path.join(
            self.output_dir, 'sigprofiler_output')
        self.processed_data = os.path.join(
            self.output_dir, 'processed_data')
        self.summary_data = os.path.join(
            self.output_dir, 'summary_data')
        self.annotated_mutations = os.path.join(
            self.output_dir, 'annotated_mutations')
        self.pattern_analysis = os.path.join(
            self.output_dir, 'pattern_analysis')
        self.merged_analysis = os.path.join(
            self.output_dir, 'merged_analysis')
        self.sankey_flows = os.path.join(
            self.output_dir, 'sankey_flows')
        self.sankey_figures = os.path.join(
            self.output_dir, 'sankey_figures')

        self.mutation_types = self.cfg.get(
            'mutation_types', ['SNV', 'DBS', 'ID', 'MNS'])
        self.doses = self.cfg.get(
            'doses', ['dA', 'dB', 'dC', 'dD', 'dE'])
        self.genome_build = self.cfg.get('genome_build', 'hg38')
        self.sigprofiler_ref = self.cfg.get('sigprofiler_ref', 'GRCh38')

        self.register_pipeline_tasks()

    def register_pipeline_tasks(self) -> None:
        """Register all pipeline tasks as asyncflow executables."""
        _task_desc_common = {}

        @self.flow.executable_task
        async def stage_annotation_preprocessing(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 3 — Build interval-tree annotations."""
            return (
                f'rigi-analysis-run annotation_preprocessing'
                f' --build {self.genome_build}'
                f' --annotation-dir {self.annotations_dir}'
            )

        @self.flow.executable_task
        async def stage_sigprofiler(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 1 — Signature Extraction via SigProfiler."""
            return (
                f'rigi-analysis-run sigprofiler'
                f' -i {self.vcf_dir}'
                f' -o {self.sigprofiler_out}'
                f' -r {self.sigprofiler_ref}'
            )

        @self.flow.executable_task
        async def stage_mutation_preprocessing(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 2 — Parse SigProfiler output into per-type pickles."""
            sig_in = os.path.join(
                self.sigprofiler_out, 'output', 'vcf_files')
            return (
                f'rigi-analysis-run mutation_preprocessing'
                f' --input-dir {sig_in}'
                f' --output {self.processed_data}'
                f' --summary {self.summary_data}'
            )

        @self.flow.executable_task
        async def stage_mutation_annotation(
                *args,
                mut_type: str,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 4 — Genomic annotation per mutation type."""
            mut_pkl = os.path.join(
                self.processed_data, f'all_{mut_type}_mutations.pkl')
            out_dir = os.path.join(self.annotated_mutations, mut_type)
            run_cmd = (
                f'rigi-analysis-run mutation_annotation'
                f' -m {mut_pkl} -a {self.annotations_dir}'
                f' -b {self.genome_build} -o {out_dir}'
            )
            mv_cmd = (
                f'mv {out_dir}/annotated_mutations.pkl'
                f' {out_dir}/all_{mut_type}_annotated.pkl'
            )
            return f'{run_cmd} && {mv_cmd}'

        @self.flow.executable_task
        async def stage_mutation_pattern(
                *args,
                mut_type: str,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 5 — Temporal pattern assignment."""
            in_dir = os.path.join(self.annotated_mutations, mut_type)
            out_dir = os.path.join(self.pattern_analysis, mut_type)
            return (
                f'rigi-analysis-run mutation_pattern_assignment'
                f' -i {in_dir} -o {out_dir} -m {mut_type}'
            )

        @self.flow.executable_task
        async def stage_merge_annotation(
                *args,
                mut_type: str,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 6 — Merge annotations with pattern assignments."""
            doses_str = ' '.join(self.doses)
            return (
                f'rigi-analysis-run merge_annotation'
                f' --annotated-dir {self.annotated_mutations}'
                f' --pattern-dir {self.pattern_analysis}'
                f' --output-dir {self.merged_analysis}'
                f' --mutation-types {mut_type}'
                f' --doses {doses_str}'
            )

        @self.flow.executable_task
        async def stage_compute_sankey(
                *args,
                mut_type: str,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 7a — Compute Sankey trajectory flows."""
            in_pkl = os.path.join(
                self.annotated_mutations, mut_type,
                f'all_{mut_type}_annotated.pkl',
            )
            out_dir = os.path.join(self.sankey_flows, mut_type)
            return (
                f'rigi-analysis-run compute_sankey'
                f' --input {in_pkl} --output-dir {out_dir}'
            )

        @self.flow.executable_task
        async def stage_render_sankey(
                *args,
                mut_type: str,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Stage 7b — Render Sankey figures."""
            in_json = os.path.join(
                self.sankey_flows, mut_type,
                'combined', 'all_chromosomes_trajectories.json',
            )
            out_png = os.path.join(
                self.sankey_figures, f'{mut_type}_combined.png')
            return (
                f'rigi-analysis-run render_sankey'
                f' --trajectories-json {in_json}'
                f' --output {out_png}'
                f' --title "Temporal Dynamics for {mut_type}"'
                f' --subtitle "All Doses, All Chromosomes"'
            )

        setattr(self, 'stage_annotation_preprocessing',
                stage_annotation_preprocessing)
        setattr(self, 'stage_sigprofiler', stage_sigprofiler)
        setattr(self, 'stage_mutation_preprocessing',
                stage_mutation_preprocessing)
        setattr(self, 'stage_mutation_annotation',
                stage_mutation_annotation)
        setattr(self, 'stage_mutation_pattern', stage_mutation_pattern)
        setattr(self, 'stage_merge_annotation', stage_merge_annotation)
        setattr(self, 'stage_compute_sankey', stage_compute_sankey)
        setattr(self, 'stage_render_sankey', stage_render_sankey)

    async def run(self) -> tuple[Any]:
        """Execute the mutation pipeline logic asynchronously.

        DAG:
            annotation_preprocessing ──┐
            sigprofiler ──► preprocessing ──┤
                                           ├──► annotation(×N)
                                                  ├──► pattern ──► merge
                                                  └──► compute_sankey ──► render_sankey
        """
        print(f'{datetime_now()} Pipeline {self.name} started')

        # Stage 1 + Stage 3: independent, can run in parallel
        annot_task = self.stage_annotation_preprocessing()
        sig_task = self.stage_sigprofiler()
        await asyncio.gather(annot_task, sig_task)

        # Stage 2: depends on sigprofiler output
        prep_task = self.stage_mutation_preprocessing(sig_task)
        await prep_task

        # Stages 4-7: per mutation type
        merge_tasks = []
        sankey_tasks = []
        for mut_type in self.mutation_types:
            # Stage 4: annotation (depends on preprocessing + annotations)
            ann_task = self.stage_mutation_annotation(
                prep_task, annot_task, mut_type=mut_type)
            # Stage 5: pattern assignment (depends on annotation)
            pat_task = self.stage_mutation_pattern(
                ann_task, mut_type=mut_type)
            # Stage 6: merge (depends on pattern)
            m_task = self.stage_merge_annotation(
                pat_task, mut_type=mut_type)
            merge_tasks.append(m_task)

            # Stage 7a-b: sankey (depends on annotation, parallel to 5-6)
            c_sankey = self.stage_compute_sankey(
                ann_task, mut_type=mut_type)
            r_sankey = self.stage_render_sankey(
                c_sankey, mut_type=mut_type)
            sankey_tasks.append(r_sankey)

        results = await asyncio.gather(*(merge_tasks + sankey_tasks))
        return results


# ----------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Workflow application for Mutation Pipeline',
        usage='rigi-analysis-run mutation_workflow [options]',
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
    """Main asynchronous execution function."""
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

    try:
        p = MutationPipeline(
            name='mutation-pipe', config=cfg, flow=flow,
            output_dir=args.output_dir,
        )
        results = await p.run()
        print(f'Pipeline results: {results}')
    except Exception:
        raise
    finally:
        await flow.shutdown()


def main() -> None:
    """Entry point."""
    asyncio.run(run_workflow())


if __name__ == '__main__':
    main()
