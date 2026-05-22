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

        self.output_dir = (output_dir or
                           self.cfg.get('output_dir',
                                         os.path.join(CURRENT_DIR,
                                                      'out_mutation')))
        os.makedirs(self.output_dir, exist_ok=True)

        self.genome_build = self.cfg.get('genome_build')
        self.sigprofiler_ref = self.cfg.get('sigprofiler_ref')
        self.mutation_types = self.cfg.get('mutation_types', [])
        self.doses = self.cfg.get('doses', [])

        self.annotations_dir = self.cfg.get('annotations_dir')
        assert os.path.exists(self.annotations_dir), \
            'Provided "annotations_dir" does not exist'
        self.vcf_dir = self.cfg.get('vcf_dir')
        assert os.path.exists(self.vcf_dir), \
            'Provided "vcf_dir" does not exist'
        self.sigprofiler_out_dir = os.path.join(
            self.output_dir, 'sigprofiler_output')
        self.processed_data_dir = os.path.join(
            self.output_dir, 'processed_data')
        self.summary_data_dir = os.path.join(
            self.output_dir, 'summary_data')
        self.annotated_mutations_dir = os.path.join(
            self.output_dir, 'annotated_mutations')
        self.pattern_analysis_dir = os.path.join(
            self.output_dir, 'pattern_analysis')
        self.merged_analysis_dir = os.path.join(
            self.output_dir, 'merged_analysis')
        self.sankey_flows_dir = os.path.join(
            self.output_dir, 'sankey_flows')
        self.sankey_figures_dir = os.path.join(
            self.output_dir, 'sankey_figures')

        self.register_pipeline_tasks()

    def register_pipeline_tasks(self) -> None:
        """Register all pipeline tasks as asyncflow executables."""

        @self.flow.executable_task
        async def stage_annotation_preprocessing(*args):
            """Stage 1a — Build interval-tree annotations (one-time setup)."""
            cmd = 'rigi-analysis-run annotation_preprocessing'
            if self.genome_build:
                cmd += f' --build {self.genome_build}'
            cmd += f' --annotation-dir {self.annotations_dir}'
            return cmd

        @self.flow.executable_task
        async def stage_sigprofiler(*args):
            """Stage 1b — Signature Extraction via SigProfiler (per-VCF run)."""
            cmd = ('rigi-analysis-run sigprofiler'
                   f' -i {self.vcf_dir}'
                   f' -o {self.sigprofiler_out_dir}')
            if self.sigprofiler_ref:
                cmd += f' -r {self.sigprofiler_ref}'
            return cmd

        @self.flow.executable_task
        async def stage_mutation_preprocessing(*args):
            """Stage 2 — Parse SigProfiler output into per-type pickles."""
            sig_in = os.path.join(self.vcf_dir, 'output', 'vcf_files')
            return ('rigi-analysis-run mutation_preprocessing'
                    f' --input-dir {sig_in}'
                    f' --output {self.processed_data_dir}'
                    f' --summary {self.summary_data_dir}')

        @self.flow.executable_task
        async def stage_mutation_annotation(*args, mut_type: str):
            """Stage 3 — Genomic annotation per mutation type."""
            mut_pkl = os.path.join(self.processed_data_dir,
                                   f'all_{mut_type}_mutations.pkl')
            out_dir = os.path.join(self.annotated_mutations_dir, mut_type)
            run_cmd = ('rigi-analysis-run mutation_annotation'
                       f' -m {mut_pkl}'
                       f' -a {self.annotations_dir}'
                       f' -o {out_dir}')
            if self.genome_build:
                run_cmd += f' -b {self.genome_build}'
            # Stage 4 expects all_<TYPE>_annotated.pkl naming scheme
            mv_cmd = ('mv'
                      f' {out_dir}/annotated_mutations.pkl'
                      f' {out_dir}/all_{mut_type}_annotated.pkl')
            return f'{run_cmd} && {mv_cmd}'

        @self.flow.executable_task
        async def stage_mutation_pattern(*args, mut_type: str):
            """Stage 4 — Temporal pattern assignment."""
            in_dir = os.path.join(self.annotated_mutations_dir, mut_type)
            out_dir = os.path.join(self.pattern_analysis_dir, mut_type)
            return ('rigi-analysis-run mutation_pattern_assignment'
                    f' -i {in_dir}'
                    f' -o {out_dir}'
                    f' -m {mut_type}')

        @self.flow.executable_task
        async def stage_merge_annotation(*args, mut_type: str):
            """Stage 5 — Merge annotations with pattern assignments."""
            return ('rigi-analysis-run merge_annotation'
                    f' --annotated-dir {self.annotated_mutations_dir}'
                    f' --pattern-dir {self.pattern_analysis_dir}'
                    f' --output-dir {self.merged_analysis_dir}'
                    f' --mutation-types {mut_type}'
                    f' --doses {" ".join(self.doses)}')

        @self.flow.executable_task
        async def stage_compute_sankey(*args, mut_type: str):
            """Stage 6 — Compute Sankey trajectory flows."""
            in_pkl = os.path.join(self.annotated_mutations_dir,
                                  mut_type,
                                  f'all_{mut_type}_annotated.pkl')
            out_dir = os.path.join(self.sankey_flows_dir, mut_type)
            return ('rigi-analysis-run compute_sankey'
                    f' --input {in_pkl}'
                    f' --output-dir {out_dir}')

        @self.flow.executable_task
        async def stage_render_sankey(*args, mut_type: str):
            """Stage 7 — Render Sankey figures."""
            in_json = os.path.join(self.sankey_flows_dir,
                                   mut_type,
                                   'combined',
                                   'all_chromosomes_trajectories.json')
            out_png = os.path.join(self.sankey_figures_dir,
                                   f'{mut_type}_combined.png')
            return ('rigi-analysis-run render_sankey'
                    f' --trajectories-json {in_json}'
                    f' --output {out_png}'
                    f' --title "Temporal Dynamics for {mut_type}"'
                    f' --subtitle "All Doses, All Chromosomes"')

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
            sigprofiler ──► preprocessing ──┐
                                            └──► annotation(×N)
                                                  ├──► pattern ──► merge
                                                  └──► compute_sankey ──► render_sankey
        """
        print(f'{datetime_now()} Pipeline {self.name} started')

        # Stage 1a + Stage 1b: independent, can run in parallel
        annot_task = self.stage_annotation_preprocessing()
        sig_task = self.stage_sigprofiler()
        await asyncio.gather(annot_task, sig_task)

        # Stage 2: depends on sigprofiler output
        prep_task = self.stage_mutation_preprocessing(sig_task)
        await prep_task

        # Stages 3-6: per mutation type
        merge_tasks = []
        sankey_tasks = []
        for mut_type in self.mutation_types:
            # Stage 3: annotation (depends on preprocessing + annotations)
            ann_task = self.stage_mutation_annotation(
                prep_task, annot_task, mut_type=mut_type)
            # Stage 4: pattern assignment (depends on annotation)
            pat_task = self.stage_mutation_pattern(
                ann_task, mut_type=mut_type)
            # Stage 5: merge (depends on pattern)
            m_task = self.stage_merge_annotation(
                pat_task, mut_type=mut_type)
            merge_tasks.append(m_task)

            # Stage 6-7: sankey (depends on annotation, parallel to 4-5)
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
        description='Workflow application for RIGI Mutation Pipeline',
        usage='rigi-analysis-workflow-mutation [options]',
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
    """Main asynchronous execution function."""
    args = parse_args()

    config_file = os.path.abspath(args.config_file)
    with open(config_file, 'r') as f:
        cfg = json.load(f)

    init_default_logger()

    backend = await get_backend(cfg.get('run_description', {}))
    flow = await WorkflowEngine.create(backend=backend)

    try:
        p = MutationPipeline(
            name='mutation-pipe',
            config=cfg,
            flow=flow,
            output_dir=args.output_dir,
        )
        results = await p.run()
        print(f'{datetime_now()} Mutation pipeline results: {results}')
    except Exception:
        raise
    finally:
        await flow.shutdown()


def main() -> None:
    """Entry point."""
    asyncio.run(run_workflow())


if __name__ == '__main__':
    main()
