"""Dry-run tests for the SV Analysis Workflow.

These tests mock the RADICAL stack so no HPC backend is needed.
They verify:
  - All 15 task stages are registered and produce correct CLI commands.
  - The DAG ordering in ``run()`` is correct.
  - Config values (tolerance, windows, threshold) appear in commands.
"""

import asyncio
import copy
import os

import pytest

# Install mocks before any workflow import
from conftest import install_mock_radical, SV_CONFIG, _MockWorkflowEngine

install_mock_radical()

from rigi_analysis.workflows.sv_workflow import SVPipeline  # noqa: E402


@pytest.fixture
def pipeline(tmp_path):
    """Create an SVPipeline wired to a mock WorkflowEngine."""
    cfg = copy.deepcopy(SV_CONFIG)
    flow = _MockWorkflowEngine()
    p = SVPipeline(
        name='test-sv',
        config=cfg,
        flow=flow,
        output_dir=str(tmp_path / 'out'),
    )
    return p


class TestSVPipelineInit:
    """Verify pipeline construction."""

    def test_output_dir_created(self, pipeline, tmp_path):
        assert os.path.isdir(str(tmp_path / 'out'))

    def test_tolerance_is_string(self, pipeline):
        assert pipeline.sv_tolerance == '1000'

    def test_windows_from_config(self, pipeline):
        assert pipeline.windows == '10,25,50,100'

    def test_mega_threshold_is_string(self, pipeline):
        assert pipeline.mega_threshold == '50000000'

    def test_annotsv_dirs(self, pipeline):
        assert pipeline.annotsv_dir == './annotsv'
        assert pipeline.annotsv_control_dir == './annotsv_control'

    def test_dbs_mutation_dir_derived(self, pipeline):
        # By default derived from mutation_merged_dir
        expected = os.path.join(pipeline.mutation_merged_dir, 'DBS')
        assert pipeline.dbs_mutation_dir == expected

    def test_missing_optional_config_keys(self, tmp_path):
        import copy
        from conftest import SV_CONFIG, _MockWorkflowEngine
        
        cfg = copy.deepcopy(SV_CONFIG)
        cfg.pop('windows', None)
        cfg.pop('single_window', None)
        
        flow = _MockWorkflowEngine()
        pipe = SVPipeline(
            name='test-missing-keys',
            config=cfg,
            flow=flow,
            output_dir=str(tmp_path / 'out'),
        )
        assert pipe.windows == '10,25,50,100'
        assert pipe.single_window is None


class TestSVTaskRegistration:
    """Verify all 15 task functions are registered."""

    EXPECTED_STAGES = [
        'stage_filter_pass_rad',
        'stage_filter_pass_ctl',
        'stage_sv_temporal',
        'stage_sv_landscape',
        'stage_repeat_analysis',
        'stage_sv_mut_correlation',
        'stage_sv_type_decay',
        'stage_inv_size_analysis',
        'stage_sv_mut_viz',
        'stage_temporal_concordance',
        'stage_temporal_concordance_viz',
        'stage_dose_stratified',
        'stage_fetch_annotations',
        'stage_categorise_genes',
        'stage_dose_visualize',
    ]

    def test_all_stages_registered(self, pipeline):
        for name in self.EXPECTED_STAGES:
            assert hasattr(pipeline, name), f'Missing stage: {name}'
            assert callable(getattr(pipeline, name))

    def test_stage_count(self):
        assert len(self.EXPECTED_STAGES) == 15


class TestSVCommandStrings:
    """Verify individual task command strings match docs."""

    def _run(self, task):
        """Await a task (which returns an asyncio.Task)."""
        return asyncio.get_event_loop().run_until_complete(task)

    def test_filter_pass_rad(self, pipeline):
        cmd = self._run(pipeline.stage_filter_pass_rad())
        assert 'rigi-analysis-run filter_pass' in cmd
        assert '-i ./annotsv' in cmd

    def test_filter_pass_ctl(self, pipeline):
        cmd = self._run(pipeline.stage_filter_pass_ctl())
        assert 'rigi-analysis-run filter_pass' in cmd
        assert '-i ./annotsv_control' in cmd

    def test_sv_temporal(self, pipeline):
        cmd = self._run(pipeline.stage_sv_temporal())
        assert 'rigi-analysis-run sv_temporal' in cmd
        assert '--tolerance 1000' in cmd
        assert '--plot' in cmd

    def test_sv_landscape(self, pipeline):
        cmd = self._run(pipeline.stage_sv_landscape())
        assert 'rigi-analysis-run SV_landscape' in cmd
        assert 'figure3_sv_landscape.png' in cmd

    def test_repeat_analysis(self, pipeline):
        cmd = self._run(pipeline.stage_repeat_analysis())
        assert 'rigi-analysis-run repeat_analysis' in cmd
        assert '--radiation-dir' in cmd
        assert '--control-dir' in cmd

    def test_sv_mutation_correlation(self, pipeline):
        cmd = self._run(pipeline.stage_sv_mut_correlation())
        assert 'rigi-analysis-run sv_mutation_correlation' in cmd
        assert '--windows 10,25,50,100' in cmd
        assert '--plot' in cmd

    def test_sv_type_decay(self, pipeline):
        cmd = self._run(pipeline.stage_sv_type_decay())
        assert 'rigi-analysis-run sv_type_specific_decay' in cmd
        assert 'DBS' in cmd

    def test_inversion_size(self, pipeline):
        cmd = self._run(pipeline.stage_inv_size_analysis())
        assert 'rigi-analysis-run inversion_size_analysis' in cmd
        assert '--mega-threshold 50000000' in cmd
        assert '--window 10' in cmd

    def test_sv_mut_viz(self, pipeline):
        cmd = self._run(pipeline.stage_sv_mut_viz())
        assert 'rigi-analysis-run sv_mut_vizualization' in cmd
        assert 'figure4_inv_dbs_unified.png' in cmd

    def test_temporal_concordance(self, pipeline):
        cmd = self._run(pipeline.stage_temporal_concordance())
        assert 'rigi-analysis-run temporal_concordance' in cmd
        assert '--window 10' in cmd

    def test_temporal_concordance_viz(self, pipeline):
        cmd = self._run(pipeline.stage_temporal_concordance_viz())
        assert 'rigi-analysis-run temporal_concordance_viz' in cmd
        assert 'figure_temporal_dynamics.png' in cmd

    def test_dose_stratified(self, pipeline):
        cmd = self._run(pipeline.stage_dose_stratified())
        assert 'rigi-analysis-run dose_stratified' in cmd
        assert '--mega-threshold 50000000' in cmd

    def test_fetch_annotations(self, pipeline):
        cmd = self._run(pipeline.stage_fetch_annotations())
        assert 'rigi-analysis-run fetch_annotations' in cmd
        assert '--gnomad' in cmd
        assert 'genes_high_dose.csv' in cmd

    def test_categorise_genes(self, pipeline):
        cmd = self._run(pipeline.stage_categorise_genes())
        assert 'rigi-analysis-run categorise_genes' in cmd
        assert 'annotated_genes.csv' in cmd

    def test_dose_visualize(self, pipeline):
        cmd = self._run(pipeline.stage_dose_visualize())
        assert 'rigi-analysis-run dose_based_visualize' in cmd
        assert 'figure5_dose_response' in cmd
        assert 'categorized_genes.csv' in cmd


class TestSVDAG:
    """Verify full pipeline DAG execution order."""

    def test_run_produces_results(self, pipeline):
        results = asyncio.get_event_loop().run_until_complete(
            pipeline.run()
        )
        assert results is not None
        # 5 terminal tasks: landscape, repeat, sv_mut_viz, concordance_viz, dose_viz
        assert len(results) == 5

    def test_filter_pass_runs_first(self, pipeline):
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands

        first_two = commands[:2]
        assert all('filter_pass' in c for c in first_two)

    def test_sv_temporal_before_correlation(self, pipeline):
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands

        sv_temp_idx = next(
            i for i, c in enumerate(commands) if 'sv_temporal' in c
        )
        corr_idx = next(
            i for i, c in enumerate(commands)
            if 'sv_mutation_correlation' in c
        )
        assert sv_temp_idx < corr_idx

    def test_dose_chain_sequential(self, pipeline):
        """dose_stratified → fetch → categorise → visualize."""
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands

        dose_idx = next(
            i for i, c in enumerate(commands)
            if 'rigi-analysis-run dose_stratified' in c
        )
        fetch_idx = next(
            i for i, c in enumerate(commands) if 'fetch_annotations' in c
        )
        cat_idx = next(
            i for i, c in enumerate(commands) if 'categorise_genes' in c
        )
        viz_idx = next(
            i for i, c in enumerate(commands)
            if 'dose_based_visualize' in c
        )
        assert dose_idx < fetch_idx < cat_idx < viz_idx

    def test_all_14_steps_executed(self, pipeline):
        """All 14 documented CLI commands should appear."""
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands

        expected_commands = [
            'filter_pass',              # ×2
            'sv_temporal',
            'SV_landscape',
            'repeat_analysis',
            'sv_mutation_correlation',
            'sv_type_specific_decay',
            'inversion_size_analysis',
            'sv_mut_vizualization',
            'temporal_concordance ',    # trailing space to avoid matching _viz
            'temporal_concordance_viz',
            'rigi-analysis-run dose_stratified',
            'fetch_annotations',
            'categorise_genes',
            'dose_based_visualize',
        ]
        cmd_text = '\n'.join(commands)
        for needle in expected_commands:
            assert needle in cmd_text, f'Missing command: {needle}'
