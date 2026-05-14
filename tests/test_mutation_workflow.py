"""Dry-run tests for the Mutation Analysis Workflow.

These tests mock the RADICAL stack so no HPC backend is needed.
They verify:
  - All 8 task stages are registered and produce correct CLI commands.
  - The DAG ordering in ``run()`` is correct.
  - Config values are plumbed into command strings.
"""

import asyncio
import copy
import os

import pytest

# Install mocks before any workflow import
from conftest import install_mock_radical, MUTATION_CONFIG, _MockWorkflowEngine

install_mock_radical()

from rigi_analysis.workflows.mutation_workflow import MutationPipeline  # noqa: E402


@pytest.fixture
def pipeline(tmp_path):
    """Create a MutationPipeline wired to a mock WorkflowEngine."""
    cfg = copy.deepcopy(MUTATION_CONFIG)
    flow = _MockWorkflowEngine()
    p = MutationPipeline(
        name='test-mutation',
        config=cfg,
        flow=flow,
        output_dir=str(tmp_path / 'out'),
    )
    return p


class TestMutationPipelineInit:
    """Verify pipeline construction."""

    def test_output_dir_created(self, pipeline, tmp_path):
        assert os.path.isdir(str(tmp_path / 'out'))

    def test_mutation_types_from_config(self, pipeline):
        assert pipeline.mutation_types == ['SNV', 'DBS']

    def test_doses_from_config(self, pipeline):
        assert pipeline.doses == ['dA', 'dB']

    def test_genome_build(self, pipeline):
        assert pipeline.genome_build == 'hg38'


class TestMutationTaskRegistration:
    """Verify all 8 task functions are registered on the pipeline."""

    EXPECTED_STAGES = [
        'stage_annotation_preprocessing',
        'stage_sigprofiler',
        'stage_mutation_preprocessing',
        'stage_mutation_annotation',
        'stage_mutation_pattern',
        'stage_merge_annotation',
        'stage_compute_sankey',
        'stage_render_sankey',
    ]

    def test_all_stages_registered(self, pipeline):
        for name in self.EXPECTED_STAGES:
            assert hasattr(pipeline, name), f'Missing stage: {name}'
            assert callable(getattr(pipeline, name))


class TestMutationCommandStrings:
    """Verify individual task command strings match docs."""

    def _run(self, task):
        """Await a task (which returns an asyncio.Task)."""
        return asyncio.get_event_loop().run_until_complete(task)

    def test_annotation_preprocessing_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_annotation_preprocessing())
        assert 'rigi-analysis-run annotation_preprocessing' in cmd
        assert '--build hg38' in cmd
        assert '--annotation-dir' in cmd

    def test_sigprofiler_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_sigprofiler())
        assert 'rigi-analysis-run sigprofiler' in cmd
        assert '-i /data/vcf_files' in cmd
        assert '-r GRCh38' in cmd

    def test_mutation_preprocessing_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_mutation_preprocessing())
        assert 'rigi-analysis-run mutation_preprocessing' in cmd
        assert '--input-dir' in cmd
        assert 'output/vcf_files' in cmd

    def test_mutation_annotation_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_mutation_annotation(mut_type='SNV'))
        assert 'rigi-analysis-run mutation_annotation' in cmd
        assert 'all_SNV_mutations.pkl' in cmd
        assert '&& mv' in cmd
        assert 'all_SNV_annotated.pkl' in cmd

    def test_mutation_pattern_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_mutation_pattern(mut_type='DBS'))
        assert 'rigi-analysis-run mutation_pattern_assignment' in cmd
        assert '-m DBS' in cmd

    def test_merge_annotation_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_merge_annotation(mut_type='SNV'))
        assert 'rigi-analysis-run merge_annotation' in cmd
        assert '--mutation-types SNV' in cmd
        assert '--doses dA dB' in cmd

    def test_compute_sankey_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_compute_sankey(mut_type='SNV'))
        assert 'rigi-analysis-run compute_sankey' in cmd
        assert 'all_SNV_annotated.pkl' in cmd

    def test_render_sankey_cmd(self, pipeline):
        cmd = self._run(pipeline.stage_render_sankey(mut_type='DBS'))
        assert 'rigi-analysis-run render_sankey' in cmd
        assert 'DBS_combined.png' in cmd
        assert 'Temporal Dynamics for DBS' in cmd


class TestMutationDAG:
    """Verify full pipeline DAG execution order."""

    def test_run_produces_results(self, pipeline):
        results = asyncio.get_event_loop().run_until_complete(
            pipeline.run()
        )
        assert results is not None
        # 2 mutation types × (merge + sankey) = 4 terminal tasks
        assert len(results) == 4

    def test_run_command_order(self, pipeline):
        """Commands should follow the documented stage ordering.

        Expected order for 2 mutation types [SNV, DBS]:
          1. annotation_preprocessing (parallel with sigprofiler)
          2. sigprofiler
          3. mutation_preprocessing
          4. mutation_annotation × 2 (SNV, DBS)
          ...
        """
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands

        # First two must be annotation_preprocessing and sigprofiler
        first_two = set(commands[:2])
        assert any('annotation_preprocessing' in c for c in first_two)
        assert any('sigprofiler' in c for c in first_two)

        # Third must be mutation_preprocessing
        assert 'mutation_preprocessing' in commands[2]

    def test_all_mutation_types_processed(self, pipeline):
        asyncio.get_event_loop().run_until_complete(pipeline.run())
        commands = pipeline.flow.recorder.commands
        cmd_text = '\n'.join(commands)
        for mt in ['SNV', 'DBS']:
            assert f'all_{mt}_mutations.pkl' in cmd_text
            assert f'all_{mt}_annotated.pkl' in cmd_text
