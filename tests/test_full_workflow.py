"""Dry-run tests for the Full RIGI Analysis Workflow.

Verifies:
  - Mutation pipeline runs to completion before SV pipeline starts.
  - Config deep-copy prevents mutation between sub-pipelines.
  - The mutation merged output directory is wired into SV config.
"""

import asyncio
import copy

import pytest

# Install mocks before any workflow import
from conftest import install_mock_radical, FULL_CONFIG, _MockWorkflowEngine

install_mock_radical()

from rigi_analysis.workflows.mutation_workflow import MutationPipeline  # noqa: E402
from rigi_analysis.workflows.sv_workflow import SVPipeline  # noqa: E402


@pytest.fixture
def pipelines(tmp_path):
    """Create both pipelines mimicking full_workflow.run_workflow()."""
    cfg = copy.deepcopy(FULL_CONFIG)
    flow = _MockWorkflowEngine()

    mut_out = str(tmp_path / 'out' / 'mutation_output')
    sv_out = str(tmp_path / 'out' / 'sv_output')

    mut_pipe = MutationPipeline(
        name='full-mutation', config=cfg, flow=flow, output_dir=mut_out,
    )

    # Deep-copy config before mutating and wire the mutation output to SV
    import os
    sv_cfg = copy.deepcopy(cfg)
    sv_cfg['mutation_merged_dir'] = os.path.join(mut_out, 'merged_data')
    os.makedirs(sv_cfg['mutation_merged_dir'], exist_ok=True)
    
    sv_pipe = SVPipeline(
        name='full-sv', config=sv_cfg, flow=flow, output_dir=sv_out,
    )

    return mut_pipe, sv_pipe, flow


class TestFullWorkflowConfig:
    """Verify config isolation between sub-pipelines."""

    def test_original_config_not_mutated(self, tmp_path):
        cfg = copy.deepcopy(FULL_CONFIG)
        original_merged = cfg.get('mutation_merged_dir')
        flow = _MockWorkflowEngine()

        MutationPipeline(
            name='test', config=cfg, flow=flow,
            output_dir=str(tmp_path / 'mut'),
        )
        sv_cfg = copy.deepcopy(cfg)
        sv_cfg['mutation_merged_dir'] = '/wired/path'

        # Original config must be unchanged
        assert cfg.get('mutation_merged_dir') == original_merged


class TestFullWorkflowDAG:
    """Verify end-to-end execution order."""

    def test_mutation_runs_before_sv(self, pipelines):
        mut_pipe, sv_pipe, flow = pipelines

        # Run mutation first
        asyncio.get_event_loop().run_until_complete(mut_pipe.run())
        mut_cmd_count = len(flow.recorder.commands)

        # Then run SV
        asyncio.get_event_loop().run_until_complete(sv_pipe.run())
        all_commands = flow.recorder.commands

        # All mutation commands should come before any SV command
        # (SV starts with filter_pass)
        first_filter = next(
            i for i, c in enumerate(all_commands) if 'filter_pass' in c
        )
        assert first_filter >= mut_cmd_count

    def test_both_pipelines_produce_results(self, pipelines):
        mut_pipe, sv_pipe, _ = pipelines

        mut_results = asyncio.get_event_loop().run_until_complete(
            mut_pipe.run()
        )
        sv_results = asyncio.get_event_loop().run_until_complete(
            sv_pipe.run()
        )

        assert mut_results is not None
        assert sv_results is not None
        assert len(mut_results) > 0
        assert len(sv_results) > 0

    def test_mutation_merged_dir_wired_to_sv_commands(self, pipelines):
        """SV tasks that consume mutation merged CSVs should reference
        the mutation output's merged directory."""
        mut_pipe, sv_pipe, flow = pipelines

        asyncio.get_event_loop().run_until_complete(mut_pipe.run())
        asyncio.get_event_loop().run_until_complete(sv_pipe.run())

        # These are the CLI commands that take --mutations / --mutation-dir
        # pointing to the mutation merged dir.  Match on the CLI verb
        # (the token right after 'rigi-analysis-run') to avoid false
        # positives from path fragments.
        sv_commands = [
            c for c in flow.recorder.commands
            if ('rigi-analysis-run sv_mutation_correlation' in c
                or 'rigi-analysis-run sv_type_specific_decay' in c
                or 'rigi-analysis-run dose_stratified' in c)
        ]
        assert len(sv_commands) > 0, 'No SV commands consuming mutation data found'
        for cmd in sv_commands:
            assert 'merged_data' in cmd, (
                f'Command should reference merged_data: {cmd}'
            )


class TestFullWorkflowUtils:
    """Test shared utilities."""

    def test_datetime_now_format(self):
        from rigi_analysis.utils.datetime_utils import datetime_now
        ts = datetime_now()
        # Should match YYYY-MM-DD HH:MM:SS
        assert len(ts) == 19
        assert ts[4] == '-'
        assert ts[10] == ' '
        assert ts[13] == ':'
