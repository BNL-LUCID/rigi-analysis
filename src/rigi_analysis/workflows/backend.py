"""Execution backend factory for workflow orchestration.

Centralises backend selection so that all workflow modules share
a single configuration point.  Swap the backend class here to
switch the entire project between Dask, RADICAL-Pilot, Dragon, etc.
without touching individual pipeline files.

Usage inside a workflow::

    from .backend import get_backend

    backend = await get_backend(run_description)
    flow    = await WorkflowEngine.create(backend=backend)
"""

from rhapsody.backends import DaskExecutionBackend

# ── Default backend class ────────────────────────────────────────────
# Change this single assignment to switch every workflow at once.
# Available options (provided by rhapsody-py):
#   - DaskExecutionBackend       (Dask Distributed)
#   - RadicalExecutionBackend    (RADICAL-Pilot)
#   - DragonExecutionBackendV3   (Dragon runtime)
#   - Concurrent                 (thread / process pool)
_BACKEND_CLASS = DaskExecutionBackend


async def get_backend(run_description: dict):
    """Create and return an execution backend instance.

    Args:
        run_description: Dictionary with resource, runtime, sandbox,
            and any backend-specific keys consumed by the chosen
            backend class.

    Returns:
        An initialised backend ready to be passed to
        ``WorkflowEngine.create(backend=...)``.
    """
    return await _BACKEND_CLASS(run_description)
