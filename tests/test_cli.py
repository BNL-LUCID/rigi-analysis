"""Tests for Workflow Command-Line Interfaces.

Provides coverage for parse_args, run_workflow, and main
across all three pipeline entry points.
"""

import argparse
import asyncio
import copy
import json
import os
from unittest.mock import patch, mock_open, AsyncMock, MagicMock

import pytest

# Install mocks before importing workflow modules
from conftest import install_mock_radical, FULL_CONFIG, _MockWorkflowEngine

install_mock_radical()

from rigi_analysis.workflows import full_workflow
from rigi_analysis.workflows import mutation_workflow
from rigi_analysis.workflows import sv_workflow


# --- Helpers ---
import builtins
original_open = builtins.open

def smart_open_mock(filename, *args, **kwargs):
    """Only mock json file reads so gettext can read its .mo files."""
    if str(filename).endswith('.json'):
        return mock_open(read_data=json.dumps(FULL_CONFIG))()
    return original_open(filename, *args, **kwargs)


class TestFullWorkflowCLI:
    def test_parse_args_success(self):
        with patch('sys.argv', ['rigi-analysis-workflow', '-c', 'config.json', '-o', '/out']):
            args = full_workflow.parse_args()
            assert args.config_file == 'config.json'
            assert args.output_dir == '/out'

    def test_parse_args_missing_required(self):
        with patch('sys.argv', ['rigi-analysis-workflow']):
            with pytest.raises(SystemExit):
                full_workflow.parse_args()

    @patch('builtins.open', side_effect=smart_open_mock)
    def test_run_workflow_success(self, mock_file, tmp_path):
        out_dir = tmp_path / 'out'
        # Pre-create the expected mutation merged dir so SVPipeline assert passes
        os.makedirs(out_dir / 'out_mutation' / 'merged_data', exist_ok=True)
        with patch('sys.argv', ['rigi-analysis-workflow', '-c', 'dummy.json', '-o', str(out_dir)]):
            asyncio.run(full_workflow.run_workflow())

    @patch('sys.argv', ['rigi-analysis-workflow', '-c', 'dummy.json'])
    @patch('builtins.open', side_effect=smart_open_mock)
    @patch('rigi_analysis.workflows.full_workflow.MutationPipeline.run')
    def test_run_workflow_exception_re_raised(self, mock_mut_run, mock_file):
        mock_mut_run.side_effect = ValueError('Test error')
        with pytest.raises(ValueError, match='Test error'):
            asyncio.run(full_workflow.run_workflow())

    @patch('rigi_analysis.workflows.full_workflow.run_workflow')
    def test_main(self, mock_run_workflow):
        # Simply patch asyncio.run and run_workflow
        with patch('asyncio.run') as mock_asyncio_run:
            full_workflow.main()
            mock_asyncio_run.assert_called_once()


class TestMutationWorkflowCLI:
    def test_parse_args_success(self):
        with patch('sys.argv', ['rigi-analysis-workflow-mutation', '-c', 'config.json', '-o', '/out']):
            args = mutation_workflow.parse_args()
            assert args.config_file == 'config.json'
            assert args.output_dir == '/out'

    @patch('builtins.open', side_effect=smart_open_mock)
    def test_run_workflow_success(self, mock_file, tmp_path):
        out_dir = str(tmp_path / 'out')
        with patch('sys.argv', ['rigi-analysis-workflow-mutation', '-c', 'dummy.json', '-o', out_dir]):
            asyncio.run(mutation_workflow.run_workflow())

    @patch('sys.argv', ['rigi-analysis-workflow-mutation', '-c', 'dummy.json'])
    @patch('builtins.open', side_effect=smart_open_mock)
    @patch('rigi_analysis.workflows.mutation_workflow.MutationPipeline.run')
    def test_run_workflow_exception_re_raised(self, mock_run, mock_file):
        mock_run.side_effect = ValueError('Test mutation error')
        with pytest.raises(ValueError, match='Test mutation error'):
            asyncio.run(mutation_workflow.run_workflow())

    @patch('rigi_analysis.workflows.mutation_workflow.run_workflow')
    def test_main(self, mock_run_workflow):
        with patch('asyncio.run') as mock_asyncio_run:
            mutation_workflow.main()
            mock_asyncio_run.assert_called_once()


class TestSVWorkflowCLI:
    def test_parse_args_success(self):
        with patch('sys.argv', ['rigi-analysis-workflow-sv', '-c', 'config.json', '-o', '/out']):
            args = sv_workflow.parse_args()
            assert args.config_file == 'config.json'
            assert args.output_dir == '/out'

    @patch('builtins.open', side_effect=smart_open_mock)
    def test_run_workflow_success(self, mock_file, tmp_path):
        out_dir = tmp_path / 'out'
        # Pre-create the expected mutation merged dir so SVPipeline assert passes
        os.makedirs(out_dir.parent / 'out_mutation' / 'merged_data', exist_ok=True)
        with patch('sys.argv', ['rigi-analysis-workflow-sv', '-c', 'dummy.json', '-o', str(out_dir)]):
            asyncio.run(sv_workflow.run_workflow())

    @patch('sys.argv', ['rigi-analysis-workflow-sv', '-c', 'dummy.json'])
    @patch('builtins.open', side_effect=smart_open_mock)
    @patch('rigi_analysis.workflows.sv_workflow.SVPipeline.run')
    def test_run_workflow_exception_re_raised(self, mock_run, mock_file):
        mock_run.side_effect = ValueError('Test SV error')
        with pytest.raises(ValueError, match='Test SV error'):
            asyncio.run(sv_workflow.run_workflow())

    @patch('rigi_analysis.workflows.sv_workflow.run_workflow')
    def test_main(self, mock_run_workflow):
        with patch('asyncio.run') as mock_asyncio_run:
            sv_workflow.main()
            mock_asyncio_run.assert_called_once()

