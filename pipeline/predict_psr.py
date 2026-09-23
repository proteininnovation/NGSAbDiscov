"""Compatibility wrapper for the legacy PSR predictor module name.

The active implementation lives in :mod:`pipeline.predict_ml` because the
Delphi backend can run any configured prediction label, including PSR, SEC,
SPR, or future Delphi targets. This wrapper keeps older scripts/imports that
refer to ``pipeline.predict_psr`` working.
"""

from .predict_ml import *  # noqa: F401,F403
