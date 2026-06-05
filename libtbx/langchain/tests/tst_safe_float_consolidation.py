"""Regression tests for the _safe_float consolidation + sanity_checker fix.

Two things are guarded here:

1. The sanity_checker string-metric bug: when phenix.model_vs_data runs (model
   placement unnecessary), metric values arrive as strings, and
   `_check_metric_anomalies` did `curr - prev` on strings -> TypeError.  Fixed
   by coercing with `_safe_float`.

2. The _safe_float consolidation: there must be exactly ONE definition of
   `_safe_float` in the langchain tree (in utils/run_utils.py), and the modules
   that use it must import that one rather than re-defining their own.  This
   keeps the (previously 5+) copies from drifting or creeping back.

Source-scan + behavioral, no PHENIX import; 2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                      # .../libtbx/langchain
_AGENT = os.path.join(_ROOT, "agent")
_UTILS = os.path.join(_ROOT, "utils")

# Modules expected to USE _safe_float via import (not define it).
_CONSUMERS = [
    "sanity_checker.py",
    "validation_history.py",
    "metric_evaluator.py",
    "metrics_analyzer.py",
    "structure_model.py",
    "display_data_model.py",
]

_CANONICAL_IMPORT = "from libtbx.langchain.utils.run_utils import _safe_float"


def _agent_file(name):
    e = os.environ.get(name.upper().replace(".PY", "_PY"))
    if e and os.path.isfile(e):
        return e
    p = os.path.join(_AGENT, name)
    return p if os.path.isfile(p) else None


def _run_utils_path():
    e = os.environ.get("RUN_UTILS_PY")
    if e and os.path.isfile(e):
        return e
    p = os.path.join(_UTILS, "run_utils.py")
    return p if os.path.isfile(p) else None


# --- Consolidation structure -------------------------------------------------

def test_single_definition_in_run_utils():
    """run_utils.py defines _safe_float exactly once."""
    path = _run_utils_path()
    if path is None:
        print("  (skip) run_utils.py not found")
        return
    n = len(re.findall(r'^def _safe_float\(', open(path).read(), re.MULTILINE))
    assert n == 1, "run_utils.py must define _safe_float exactly once (found %d)" % n


def test_no_other_definitions_in_agent():
    """No agent module re-defines _safe_float (it must be imported)."""
    if not os.path.isdir(_AGENT):
        print("  (skip) agent dir not found")
        return
    offenders = []
    for fn in os.listdir(_AGENT):
        if not fn.endswith(".py"):
            continue
        src = open(os.path.join(_AGENT, fn)).read()
        if re.search(r'^def _safe_float\(', src, re.MULTILINE):
            offenders.append(fn)
    assert not offenders, \
        "these agent modules still DEFINE _safe_float instead of importing it: %s" \
        % ", ".join(offenders)


def test_consumers_import_canonical():
    """Each consumer that uses _safe_float imports the canonical one."""
    for fn in _CONSUMERS:
        path = _agent_file(fn)
        if path is None:
            print("  (skip) %s not found" % fn)
            continue
        src = open(path).read()
        if "_safe_float" not in src:
            continue  # module no longer uses it; fine
        assert _CANONICAL_IMPORT in src, \
            "%s uses _safe_float but does not import the canonical run_utils one" % fn
        assert not re.search(r'^def _safe_float\(', src, re.MULTILINE), \
            "%s must not define its own _safe_float" % fn


# --- Behavior of _safe_float -------------------------------------------------

def _load_run_utils():
    path = _run_utils_path()
    if path is None:
        return None
    # exec just the function to avoid importing run_utils' heavier deps.
    # Capture from the def line through the except-branch 'return None'.
    src = open(path).read()
    m = re.search(
        r'def _safe_float\(val\):.*?except \(ValueError, TypeError\):\s*\n\s*return None\n',
        src, re.DOTALL)
    if not m:
        return None
    ns = {}
    exec(m.group(0), ns)
    return ns.get("_safe_float")


def test_safe_float_behavior():
    sf = _load_run_utils()
    if sf is None:
        print("  (skip) could not load _safe_float")
        return
    assert sf("0.385") == 0.385      # the JSON-string case
    assert sf(0.385) == 0.385
    assert sf(0.0) == 0.0            # not None (falsy but valid)
    assert sf(None) is None
    assert sf("N/A") is None
    assert sf("") is None
    assert sf([1, 2]) is None        # non-convertible type -> None, no raise


# --- sanity_checker bug fix (source-scan; behavior covered above) ------------

def test_sanity_checker_coerces_metrics():
    """sanity_checker._check_metric_anomalies coerces metric values before
    arithmetic (guards the string-subtraction crash)."""
    path = _agent_file("sanity_checker.py")
    if path is None:
        print("  (skip) sanity_checker.py not found")
        return
    src = open(path).read()
    # both metric blocks must coerce via _safe_float, not use raw .get()
    assert '_safe_float(prev.get("r_free"))' in src
    assert '_safe_float(curr.get("r_free"))' in src
    assert '_safe_float(prev.get("map_cc"))' in src
    assert '_safe_float(curr.get("map_cc"))' in src
    # guards must use "is not None" (not the falsy "if prev and curr" which
    # would wrongly skip a legitimate 0.0)
    assert "prev_rfree is not None and curr_rfree is not None" in src
    assert "prev_cc is not None and curr_cc is not None" in src
    # the raw-subtraction form must be gone
    assert "change = curr_rfree - prev_rfree" in src  # arithmetic still present
    assert 'curr_rfree = curr.get("r_free")' not in src  # raw assignment gone


_TESTS = [
    test_single_definition_in_run_utils,
    test_no_other_definitions_in_agent,
    test_consumers_import_canonical,
    test_safe_float_behavior,
    test_sanity_checker_coerces_metrics,
]


def run_all_tests():
    for fn in _TESTS:
        fn()
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    p = f = 0
    for fn in _TESTS:
        try:
            fn(); print("  PASS: %s" % fn.__name__); p += 1
        except AssertionError as e:
            print("  FAIL: %s -- %s" % (fn.__name__, e)); f += 1
    print("\n%d passed, %d failed" % (p, f))
