"""Dry-run test infrastructure for RIGI Analysis workflows.

Provides mock replacements for radical.asyncflow and rhapsody so that
the workflow modules can be imported and exercised locally without the
RADICAL stack installed.
"""

import asyncio
import json
import os
import sys
import types
from unittest.mock import MagicMock


# ── Minimal config fixtures ──────────────────────────────────────────

MUTATION_CONFIG = {
    'vcf_dir': '/data/vcf_files',
    'annotations_dir': '/data/annotations',
    'mutation_types': ['SNV', 'DBS'],
    'doses': ['dA', 'dB'],
    'genome_build': 'hg38',
    'sigprofiler_ref': 'GRCh38',
    'run_description': {},
}

SV_CONFIG = {
    'annotsv_dir': '/data/annotsv',
    'annotsv_control_dir': '/data/annotsv_control',
    'mutation_merged_dir': '/data/mutation/merged',
    'gnomad_file': '/data/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz',
    'sv_tolerance': 1000,
    'windows': '10,25,50,100',
    'single_window': 10,
    'mega_threshold': 50000000,
    'run_description': {},
}

FULL_CONFIG = {
    'vcf_dir': '/data/vcf_files',
    'annotations_dir': '/data/annotations',
    'mutation_types': ['SNV'],
    'doses': ['dA'],
    'genome_build': 'hg38',
    'sigprofiler_ref': 'GRCh38',
    'annotsv_dir': '/data/annotsv',
    'annotsv_control_dir': '/data/annotsv_control',
    'gnomad_file': '/data/gnomad.bgz',
    'sv_tolerance': 1000,
    'windows': '10,25,50,100',
    'single_window': 10,
    'mega_threshold': 50000000,
    'run_description': {},
}


# ── Mock RADICAL stack ───────────────────────────────────────────────

class _TaskRecorder:
    """Records every command string produced by executable_task calls.

    The ``executable_task`` decorator wraps a coroutine so that when
    called it returns an ``asyncio.Task`` (instead of a bare coroutine).
    This matches RADICAL's behaviour where task handles can be:
      - awaited multiple times (idempotent)
      - passed as positional dependency args to downstream tasks
    """

    def __init__(self):
        self.commands: list[str] = []

    def executable_task(self, fn):
        recorder = self

        async def _inner(*args, **kwargs):
            # Resolve any Task/Future dependency args first
            resolved = []
            for a in args:
                if isinstance(a, asyncio.Task) or asyncio.isfuture(a):
                    resolved.append(await a)
                elif asyncio.iscoroutine(a):
                    resolved.append(await a)
                else:
                    resolved.append(a)

            cmd = await fn(*resolved, **kwargs)
            recorder.commands.append(cmd)
            return cmd

        def wrapper(*args, **kwargs):
            # Return an asyncio.Task so the handle can be awaited
            # multiple times and passed as a dependency arg
            loop = asyncio.get_event_loop()
            return loop.create_task(_inner(*args, **kwargs))

        wrapper.__name__ = fn.__name__
        wrapper.__qualname__ = fn.__qualname__
        return wrapper


class _MockWorkflowEngine:
    """Minimal stand-in for ``radical.asyncflow.WorkflowEngine``."""

    def __init__(self):
        self.recorder = _TaskRecorder()
        self.executable_task = self.recorder.executable_task

    @classmethod
    async def create(cls, backend=None):
        return cls()

    async def shutdown(self):
        pass


def install_mock_radical():
    """Inject mock ``radical`` and ``rhapsody`` packages into sys.modules.

    Must be called **before** importing any workflow module.
    """
    # radical.asyncflow
    radical = types.ModuleType('radical')
    radical_asyncflow = types.ModuleType('radical.asyncflow')
    radical_asyncflow.WorkflowEngine = _MockWorkflowEngine
    radical_asyncflow_logging = types.ModuleType('radical.asyncflow.logging')
    radical_asyncflow_logging.init_default_logger = lambda *a, **k: None

    sys.modules['radical'] = radical
    sys.modules['radical.asyncflow'] = radical_asyncflow
    sys.modules['radical.asyncflow.logging'] = radical_asyncflow_logging

    # rhapsody
    rhapsody = types.ModuleType('rhapsody')
    rhapsody_backends = types.ModuleType('rhapsody.backends')

    async def _mock_backend(desc):
        return MagicMock()

    rhapsody_backends.DaskExecutionBackend = _mock_backend
    rhapsody_backends.RadicalExecutionBackend = _mock_backend
    sys.modules['rhapsody'] = rhapsody
    sys.modules['rhapsody.backends'] = rhapsody_backends
