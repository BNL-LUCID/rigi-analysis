"""Structural Variant Analysis Workflow.

Executes a multi-stage SV Analysis pipeline
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


class SVPipeline:
    """SV Pipeline Execution Class.

    Orchestrates the 14-step SV analysis pipeline described in
    ``docs/sv/README.md``.  Each step is registered as a RADICAL
    asyncflow executable task; the ``run`` method wires them into the
    correct dependency DAG.
    """

    def __init__(
            self,
            name: str,
            config: dict,
            flow: WorkflowEngine,
            output_dir: Optional[str] = None,
    ):
        """Initialize the SV pipeline.

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
                            os.path.join(CURRENT_DIR, 'out_sv'))
        )
        os.makedirs(self.output_dir, exist_ok=True)

        # Directories & inputs
        self.annotsv_dir = self.cfg.get('annotsv_dir', './annotsv')
        self.annotsv_control_dir = self.cfg.get(
            'annotsv_control_dir', './annotsv_control')

        self.mutation_merged_dir = self.cfg.get(
            'mutation_merged_dir', '../mutation/merged')
        self.dbs_mutation_dir = os.path.join(
            self.mutation_merged_dir, 'DBS')
        self.gnomad_file = self.cfg.get(
            'gnomad_file',
            './gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz',
        )

        self.sv_tolerance = str(self.cfg.get('sv_tolerance', '1000'))
        self.windows = str(self.cfg.get('windows', '10,25,50,100'))
        self.single_window = str(self.cfg.get('single_window', '10'))
        self.mega_threshold = str(
            self.cfg.get('mega_threshold', '50000000'))

        # Output paths
        self.annotsv_passed = os.path.join(
            self.output_dir, 'annotsv_passed')
        self.annotsv_passed_control = os.path.join(
            self.output_dir, 'annotsv_passed_control')
        self.sv_temporal = os.path.join(
            self.output_dir, 'sv_temporal')
        self.sv_temporal_catalog = os.path.join(
            self.sv_temporal, 'sv_temporal_catalog.csv')
        self.sv_correlation = os.path.join(
            self.output_dir, 'sv_correlation')
        self.sv_type_decay = os.path.join(
            self.output_dir, 'sv_type_decay')
        self.inv_size = os.path.join(
            self.output_dir, 'inv_size')
        self.dose_stratified = os.path.join(
            self.output_dir, 'dose_stratified')

        self.register_pipeline_tasks()

    def register_pipeline_tasks(self) -> None:
        """Register all pipeline tasks as asyncflow executables."""
        _task_desc_common = {}

        # -- Step 1: PASS filtering -----------------------------------

        @self.flow.executable_task
        async def stage_filter_pass_rad(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 1a — PASS filter on radiation AnnotSV outputs."""
            return (
                f'rigi-analysis-run filter_pass'
                f' -i {self.annotsv_dir}'
                f' -o {self.annotsv_passed}'
            )

        @self.flow.executable_task
        async def stage_filter_pass_ctl(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 1b — PASS filter on control AnnotSV outputs."""
            return (
                f'rigi-analysis-run filter_pass'
                f' -i {self.annotsv_control_dir}'
                f' -o {self.annotsv_passed_control}'
            )

        # -- Step 2-4: Temporal + Landscape + Repeat ------------------

        @self.flow.executable_task
        async def stage_sv_temporal(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 2 — Temporal pattern assignment."""
            return (
                f'rigi-analysis-run sv_temporal'
                f' {self.annotsv_passed}'
                f' --output {self.sv_temporal}'
                f' --tolerance {self.sv_tolerance}'
                f' --plot'
            )

        @self.flow.executable_task
        async def stage_sv_landscape(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 3 — SV landscape figure (Fig 3)."""
            out_png = os.path.join(
                self.output_dir, 'figure3_sv_landscape.png')
            return (
                f'rigi-analysis-run SV_landscape'
                f' --annotsv-dir {self.annotsv_passed}'
                f' --output {out_png}'
            )

        @self.flow.executable_task
        async def stage_repeat_analysis(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 4 — Repeat analysis (Fig S2)."""
            out_png = os.path.join(
                self.output_dir, 'figure_s2_repeat_analysis.png')
            return (
                f'rigi-analysis-run repeat_analysis'
                f' --radiation-dir {self.annotsv_passed}'
                f' --control-dir {self.annotsv_passed_control}'
                f' --output {out_png}'
            )

        # -- Steps 5-8: Correlation branch ----------------------------

        @self.flow.executable_task
        async def stage_sv_mut_correlation(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 5 — Breakpoint-proximal mutation enrichment."""
            return (
                f'rigi-analysis-run sv_mutation_correlation'
                f' --sv-catalog {self.sv_temporal_catalog}'
                f' --mutations {self.mutation_merged_dir}'
                f' --output {self.sv_correlation}'
                f' --windows {self.windows} --plot'
            )

        @self.flow.executable_task
        async def stage_sv_type_decay(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 6 — Distance decay (Fig S3)."""
            return (
                f'rigi-analysis-run sv_type_specific_decay'
                f' --sv-catalog {self.sv_temporal_catalog}'
                f' --mutations {self.dbs_mutation_dir}'
                f' --output {self.sv_type_decay}'
                f' --windows {self.windows}'
            )

        @self.flow.executable_task
        async def stage_inv_size_analysis(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 7 — INV-size × DBS coupling (Fig 4 C, D)."""
            out_png = os.path.join(self.inv_size, 'inv_size.png')
            return (
                f'rigi-analysis-run inversion_size_analysis'
                f' --annotsv-dir {self.annotsv_passed}'
                f' --dbs-dir {self.dbs_mutation_dir}'
                f' --window {self.single_window}'
                f' --mega-threshold {self.mega_threshold}'
                f' --output {out_png}'
            )

        @self.flow.executable_task
        async def stage_sv_mut_viz(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 8 — INV-DBS coupling figure (Fig 4)."""
            out_png = os.path.join(
                self.output_dir, 'figure4_inv_dbs_unified.png')
            return (
                f'rigi-analysis-run sv_mut_vizualization'
                f' --correlation-dir {self.sv_correlation}'
                f' --sv-type-dir {self.sv_type_decay}'
                f' --size-analysis {self.inv_size}'
                f' --output {out_png}'
            )

        # -- Steps 9-10: Temporal concordance branch ------------------

        @self.flow.executable_task
        async def stage_temporal_concordance(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 9 — Temporal concordance (§4.8)."""
            return (
                f'rigi-analysis-run temporal_concordance'
                f' --annotsv-dir {self.annotsv_passed}'
                f' --dbs-data {self.dbs_mutation_dir}'
                f' --window {self.single_window}'
            )

        @self.flow.executable_task
        async def stage_temporal_concordance_viz(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 10 — Concordance figure (Fig S4)."""
            out_png = os.path.join(
                self.output_dir, 'figure_temporal_dynamics.png')
            return (
                f'rigi-analysis-run temporal_concordance_viz'
                f' --annotsv-dir {self.annotsv_passed}'
                f' --dbs-data {self.dbs_mutation_dir}'
                f' --output {out_png}'
                f' --window {self.single_window}'
            )

        # -- Steps 11-14: Dose-stratified branch ----------------------

        @self.flow.executable_task
        async def stage_dose_stratified(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 11 — Dose-stratified INV-DBS pairs (Fig 5B)."""
            return (
                f'rigi-analysis-run dose_stratified'
                f' --annotsv-dir {self.annotsv_passed}'
                f' --mutation-dir {self.mutation_merged_dir}'
                f' --output-dir {self.dose_stratified}'
                f' --window {self.single_window}'
                f' --mega-threshold {self.mega_threshold}'
            )

        @self.flow.executable_task
        async def stage_fetch_annotations(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 12 — Fetch gene annotations (requires network)."""
            high_csv = os.path.join(
                self.dose_stratified, 'genes_high_dose.csv')
            low_csv = os.path.join(
                self.dose_stratified, 'genes_low_dose.csv')
            out_csv = os.path.join(
                self.output_dir, 'annotated_genes.csv')
            return (
                f'rigi-analysis-run fetch_annotations'
                f' --genes-high {high_csv}'
                f' --genes-low {low_csv}'
                f' --gnomad {self.gnomad_file}'
                f' --output {out_csv}'
            )

        @self.flow.executable_task
        async def stage_categorise_genes(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 13 — Categorise genes (6 functional categories)."""
            in_csv = os.path.join(
                self.output_dir, 'annotated_genes.csv')
            out_csv = os.path.join(
                self.output_dir, 'categorized_genes.csv')
            return (
                f'rigi-analysis-run categorise_genes'
                f' --input {in_csv} --output {out_csv}'
            )

        @self.flow.executable_task
        async def stage_dose_visualize(
                *args,
                task_description: dict = _task_desc_common,  # noqa: B006
        ):
            """Step 14 — Dose-response figure (Fig 5)."""
            low_csv = os.path.join(
                self.dose_stratified, 'inv_dbs_pairs_low.csv')
            high_csv = os.path.join(
                self.dose_stratified, 'inv_dbs_pairs_high.csv')
            cat_csv = os.path.join(
                self.output_dir, 'categorized_genes.csv')
            out_prefix = os.path.join(
                self.output_dir, 'figure5_dose_response')
            return (
                f'rigi-analysis-run dose_based_visualize'
                f' --inv-dbs-low {low_csv}'
                f' --inv-dbs-high {high_csv}'
                f' --categorized-genes {cat_csv}'
                f' --output {out_prefix}'
            )

        setattr(self, 'stage_filter_pass_rad', stage_filter_pass_rad)
        setattr(self, 'stage_filter_pass_ctl', stage_filter_pass_ctl)
        setattr(self, 'stage_sv_temporal', stage_sv_temporal)
        setattr(self, 'stage_sv_landscape', stage_sv_landscape)
        setattr(self, 'stage_repeat_analysis', stage_repeat_analysis)
        setattr(self, 'stage_sv_mut_correlation',
                stage_sv_mut_correlation)
        setattr(self, 'stage_sv_type_decay', stage_sv_type_decay)
        setattr(self, 'stage_inv_size_analysis',
                stage_inv_size_analysis)
        setattr(self, 'stage_sv_mut_viz', stage_sv_mut_viz)
        setattr(self, 'stage_temporal_concordance',
                stage_temporal_concordance)
        setattr(self, 'stage_temporal_concordance_viz',
                stage_temporal_concordance_viz)
        setattr(self, 'stage_dose_stratified', stage_dose_stratified)
        setattr(self, 'stage_fetch_annotations',
                stage_fetch_annotations)
        setattr(self, 'stage_categorise_genes', stage_categorise_genes)
        setattr(self, 'stage_dose_visualize', stage_dose_visualize)

    async def run(self) -> tuple[Any]:
        """Execute the SV pipeline logic asynchronously.

        DAG (simplified):
            filter_pass_rad ──┬──► sv_temporal ──┬──► correlation ──┐
                              │                  └──► type_decay ───┤
                              ├──► sv_landscape (terminal)          ├──► sv_mut_viz
                              ├──► inv_size ────────────────────────┘
                              ├──► concordance ──► concordance_viz
                              └──► dose_stratified ──► fetch ──► categorise ──► dose_viz
            filter_pass_ctl ──┴──► repeat_analysis (terminal)
        """
        print(f'{datetime_now()} Pipeline {self.name} started')

        # Step 1: PASS filtering (radiation + control in parallel)
        f_rad = self.stage_filter_pass_rad()
        f_ctl = self.stage_filter_pass_ctl()
        await asyncio.gather(f_rad, f_ctl)

        # Steps 2-4 & 7 & 9: all depend only on f_rad (or f_rad+f_ctl)
        # Launch them all; only await sv_temporal before steps 5-6
        sv_temp = self.stage_sv_temporal(f_rad)
        sv_land = self.stage_sv_landscape(f_rad)           # terminal
        rep_ana = self.stage_repeat_analysis(f_rad, f_ctl)  # terminal
        inv_sz_ana = self.stage_inv_size_analysis(f_rad)    # needed by 8
        tmp_con = self.stage_temporal_concordance(f_rad)    # step 9
        tmp_cn_vz = self.stage_temporal_concordance_viz(tmp_con)  # step 10

        # Steps 11-14: dose-stratified chain (depends only on f_rad)
        dos_str = self.stage_dose_stratified(f_rad)
        fch_ann = self.stage_fetch_annotations(dos_str)
        cat_gns = self.stage_categorise_genes(fch_ann)
        dos_vz = self.stage_dose_visualize(cat_gns)

        # Wait for sv_temporal only to unblock steps 5-6
        await sv_temp

        # Steps 5-6: depend on sv_temporal catalog
        sv_mut_cor = self.stage_sv_mut_correlation(sv_temp)
        sv_typ_dec = self.stage_sv_type_decay(sv_temp)

        # Step 8: depends on steps 5, 6, 7 outputs
        sv_mut_vz = self.stage_sv_mut_viz(
            sv_mut_cor, sv_typ_dec, inv_sz_ana)

        # Gather all terminal tasks
        results = await asyncio.gather(
            sv_land, rep_ana,       # terminal figures
            sv_mut_vz,              # step 8 figure
            tmp_cn_vz,              # step 10 figure
            dos_vz,                 # step 14 figure
        )
        return results


# ------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Workflow application for SV Pipeline',
        usage='rigi-analysis-run sv_workflow [options]',
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
        p = SVPipeline(
            name='sv-pipe', config=cfg, flow=flow,
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
