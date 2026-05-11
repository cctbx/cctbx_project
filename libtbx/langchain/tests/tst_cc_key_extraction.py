"""Tests for v116.10 CC key typo fix in _generate_structure_report.

The bug: `_generate_structure_report` extracts the best CC value from
cycle metrics to decide whether the workflow achieved good results
(`_metrics_good`).  v115.05 introduced a typo in the lookup key:
it asked for "map_model_cc" but the cycle metrics dict uses
"model_map_cc" (the canonical key used everywhere else in
ai_agent.py: lines 3132, 8564, 9084, 9185, 9187).

The consequence: cryo-EM workflows that produced excellent
model-map CC (say 0.85) were extracted as `_best_cc = None`,
giving `_metrics_good = False`, and routed to
`generate_stopped_report` — which prints "SESSION STOPPED —
INCOMPLETE" / "No structure model available."

This test reproduces the metric-extraction logic in isolation
(without instantiating the full ai_agent class) and locks in:

  1. `model_map_cc` is recognized (the canonical key)
  2. `map_model_cc` is still recognized (legacy/typo data)
  3. `map_cc` is recognized (workflow_engine variant)
  4. `cc_mask` is recognized (real_space_refine output)
  5. The fallback structure_model path uses the same key set
  6. A successful cryo-EM cycle correctly sets `_metrics_good=True`

Fail-fast convention: prints PASS on success, raises on failure.
"""

from __future__ import absolute_import, division, print_function

import re
import os
import sys


# --- Locate ai_agent.py across PHENIX layouts ------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_HERE)


def _find_ai_agent_py():
    """Find phenix/programs/ai_agent.py across common PHENIX layouts."""
    tried = []
    cand = os.path.join(_PROJECT_ROOT, "phenix", "programs",
                        "ai_agent.py")
    tried.append(cand)
    if os.path.exists(cand):
        return cand
    cand = os.path.normpath(os.path.join(
        _PROJECT_ROOT, "..", "..", "..",
        "phenix", "phenix", "programs", "ai_agent.py"))
    tried.append(cand)
    if os.path.exists(cand):
        return cand
    try:
        from libtbx.env_config import unpickle
        env = unpickle()
        phenix_dist = env.find_dist_path('phenix', default=None)
        if phenix_dist:
            for sub in ("phenix/programs/ai_agent.py",
                        "programs/ai_agent.py"):
                cand = os.path.join(phenix_dist, sub)
                tried.append(cand)
                if os.path.exists(cand):
                    return cand
    except Exception as e:
        tried.append("libtbx.env_config: %s" % e)
    phenix_env = os.environ.get("PHENIX")
    if phenix_env:
        cand = os.path.join(
            phenix_env, "modules", "phenix", "phenix",
            "programs", "ai_agent.py")
        tried.append(cand)
        if os.path.exists(cand):
            return cand
    raise RuntimeError(
        "Cannot locate ai_agent.py. Tried:\n  %s"
        % "\n  ".join(tried))


# --- Source-level invariants -----------------------------------------------

def test_canonical_cc_key_is_recognized():
    """The cycle-metric extraction must look up `model_map_cc`.

    This is a source-level invariant: if a future refactor removes
    the `model_map_cc` lookup, the cryo-EM success-detection breaks
    silently.
    """
    print("Test: canonical_cc_key_is_recognized")
    with open(_find_ai_agent_py()) as f:
        source = f.read()

    # Find the `_cc = (` block that lives near the "model_map_cc"
    # extraction and assert it includes the canonical key.
    # Match a few lines around `_la.get("model_map_cc")`.
    m = re.search(
        r"_cc = \(\s*[^)]*?_la\.get\(\"model_map_cc\"\)[^)]*?\)",
        source, re.DOTALL)
    assert m is not None, (
        "Cycle metrics CC extraction missing `_la.get(\"model_map_cc\")`. "
        "This is the canonical key used by all writers in ai_agent.py; "
        "without it, cryo-EM workflows are misclassified as "
        "INCOMPLETE.")
    print("  PASS")


def test_fallback_structure_model_uses_canonical_cc_key():
    """The structure-model fallback path must also use `model_map_cc`."""
    print("Test: fallback_structure_model_uses_canonical_cc_key")
    with open(_find_ai_agent_py()) as f:
        source = f.read()

    # Find the for-loop that iterates CC fallback keys.
    m = re.search(
        r'for _key in \(([^)]*?(?:map_cc|cc_mask)[^)]*?)\):',
        source, re.DOTALL)
    assert m is not None, (
        "Could not find structure_model fallback CC key loop")
    keys_blob = m.group(1)
    assert '"model_map_cc"' in keys_blob, (
        "Structure-model fallback CC extraction missing "
        "\"model_map_cc\" in key list: %s" % keys_blob)
    print("  PASS")


# --- Behavioral invariants (extraction logic in isolation) -----------------

def _extract_best_cc(cycles_metrics):
    """Mimic the cycle-CC extraction logic from ai_agent.py.

    This must be kept in sync with the production code.  If the
    production code changes, update this helper and the test logic
    below together.

    Args:
      cycles_metrics: list of dicts representing the contents of
        session.data["cycles"][i]["metrics"] (or "log_analysis").

    Returns:
      float or None.  None means no CC was found.
    """
    best_cc = None
    for _la in cycles_metrics:
        if not isinstance(_la, dict):
            continue
        _cc = (
            _la.get("model_map_cc")
            or _la.get("map_model_cc")
            or _la.get("map_cc")
            or _la.get("cc_mask")
        )
        if _cc is not None:
            try:
                _cc_f = float(_cc)
                if best_cc is None or _cc_f > best_cc:
                    best_cc = _cc_f
            except (ValueError, TypeError):
                pass
    return best_cc


def test_canonical_key_extraction():
    """`model_map_cc` (canonical) is extracted."""
    print("Test: canonical_key_extraction")
    cycles = [{"model_map_cc": 0.85}]
    assert _extract_best_cc(cycles) == 0.85
    print("  PASS")


def test_legacy_typo_key_still_extracted():
    """`map_model_cc` (the pre-v116.10 typo spelling) is still
    extracted, so legacy session data isn't misclassified."""
    print("Test: legacy_typo_key_still_extracted")
    cycles = [{"map_model_cc": 0.72}]
    assert _extract_best_cc(cycles) == 0.72
    print("  PASS")


def test_workflow_engine_variant_key():
    """`map_cc` (workflow_engine spelling) is extracted."""
    print("Test: workflow_engine_variant_key")
    cycles = [{"map_cc": 0.68}]
    assert _extract_best_cc(cycles) == 0.68
    print("  PASS")


def test_cc_mask_key():
    """`cc_mask` (phenix.real_space_refine output) is extracted."""
    print("Test: cc_mask_key")
    cycles = [{"cc_mask": 0.91}]
    assert _extract_best_cc(cycles) == 0.91
    print("  PASS")


def test_best_value_wins():
    """If multiple cycles record CC, the BEST (highest) wins."""
    print("Test: best_value_wins")
    cycles = [
        {"model_map_cc": 0.45},
        {"model_map_cc": 0.78},  # best
        {"model_map_cc": 0.62},
    ]
    assert _extract_best_cc(cycles) == 0.78
    print("  PASS")


def test_no_cc_present_returns_none():
    """A cycle without any CC field returns None (not 0.0)."""
    print("Test: no_cc_present_returns_none")
    cycles = [{"r_free": 0.25}, {}, {"random_key": 5}]
    assert _extract_best_cc(cycles) is None
    print("  PASS")


def test_string_cc_value_coerced():
    """CC stored as a string is coerced to float (defense in
    depth — happens with some log parsers)."""
    print("Test: string_cc_value_coerced")
    cycles = [{"model_map_cc": "0.83"}]
    assert _extract_best_cc(cycles) == 0.83
    print("  PASS")


def test_garbage_cc_value_skipped():
    """A non-numeric CC string is skipped without crashing."""
    print("Test: garbage_cc_value_skipped")
    cycles = [
        {"model_map_cc": "not_a_number"},
        {"model_map_cc": 0.5},
    ]
    assert _extract_best_cc(cycles) == 0.5
    print("  PASS")


def test_motivating_case_cryoem_workflow():
    """The bug's motivating case: a successful cryo-EM workflow
    with metrics > 0.6 is classified as `_metrics_good = True`,
    NOT INCOMPLETE.

    This is the scenario from terwill's bug report: full cryo-EM
    workflow, successful, produced metrics — but the report
    falsely said SESSION STOPPED - INCOMPLETE because the CC
    extraction missed `model_map_cc` due to the v115.05 typo.
    """
    print("Test: motivating_case_cryoem_workflow")
    # Simulate a successful cryo-EM run's cycles
    cycles = [
        {"r_free": None, "model_map_cc": 0.55},  # early cycle
        {"r_free": None, "model_map_cc": 0.71},
        {"r_free": None, "model_map_cc": 0.83},  # converged
    ]
    best_cc = _extract_best_cc(cycles)
    # Metrics-good threshold from ai_agent.py
    metrics_good = (best_cc is not None and best_cc > 0.6)
    assert metrics_good is True, (
        "Cryo-EM workflow with model_map_cc=0.83 must be "
        "_metrics_good=True (the threshold is 0.6). "
        "Got best_cc=%s" % best_cc)
    print("  PASS — best_cc=%.2f, _metrics_good=True" % best_cc)


# --- Entry point -----------------------------------------------------------

def run_all_tests(verbose=False):
    tests = [
        test_canonical_cc_key_is_recognized,
        test_fallback_structure_model_uses_canonical_cc_key,
        test_canonical_key_extraction,
        test_legacy_typo_key_still_extracted,
        test_workflow_engine_variant_key,
        test_cc_mask_key,
        test_best_value_wins,
        test_no_cc_present_returns_none,
        test_string_cc_value_coerced,
        test_garbage_cc_value_skipped,
        test_motivating_case_cryoem_workflow,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            print("  FAIL: %s: %s" % (type(e).__name__, e))
            failed += 1
    print("\n%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
