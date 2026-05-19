"""Tests for v117.3 extension to stop-intent phrasing recognition.

The v117.3 plan adds new patterns and markers so the directive
extractor recognizes phrasings like:

    "stop the workflow immediately after phaser completes"
    "phaser is the last step"

These are real-world phrasings that the v117.2 regex helpers and
prompt schema missed, causing the `explicit_stop_after_phaser` LLM
test to fail 0/5 on openai.

This file holds the boundary tests (K1-K8b) that lock the v117.3
patches.  By design:

  Pre-v117.3 baseline (failing tests on v117.2):
    K1, K2  — _is_stop_after_requested returns False on new phrasings
    K6      — _imperative_marker_nearby returns False on
              explicit_stop_after_phaser fixture

  Pre-v117.3 baseline (passing — preservation tests):
    K4      — AF_7mjs's stripped advice continues to return False
              (the v117.1 K2 grounding test depends on this)
    K5      — descriptive "validation runs after refinement completes"
              continues to return False (false-positive guard)
    K7a/b   — existing phrasings ("density modify and stop",
              "refine and stop") continue to return True
    K8a/b   — existing negation cases continue to return False

Post-v117.3 (after applying patches), all 9 tests pass.

Test K3 was considered and dropped from the v117.3 plan because the
existing ',\\s*stop\\b' pattern already catches "after phaser
completes, stop" (verified during pre-implementation baseline check).
Tests K9 and K10 (cross-sentence safety) were dropped because the
v117.2 pre-existing handling of \\bstop\\s+if\\b on inputs like
"Stop if you see an error..." is a SEPARATE issue, not a v117.3
concern.  See next_steps.md for tracking.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOTS = (
    os.path.abspath(os.path.join(_HERE, "..")),
    os.path.abspath(_HERE),
)
for _root in _ROOTS:
    if _root not in sys.path:
        sys.path.insert(0, _root)

for mod in ("libtbx", "libtbx.langchain", "libtbx.langchain.agent"):
    if mod not in sys.modules:
        sys.modules[mod] = types.ModuleType(mod)

if "agent.program_registry" not in sys.modules:
    pr = types.ModuleType("agent.program_registry")
    class _Registry:
        def __init__(self, *a, **k): pass
        def get(self, *a, **k): return None
        def all_programs(self): return []
        def get_all_program_keys(self): return []
    pr.ProgramRegistry = _Registry
    sys.modules["agent.program_registry"] = pr
    sys.modules["libtbx.langchain.agent.program_registry"] = pr

try:
    from libtbx.langchain.agent import directive_extractor as de
except ImportError:
    from agent import directive_extractor as de


# =====================================================================
# Test inputs
# =====================================================================

# AF_7mjs preservation fixture (matches the C1 test infrastructure
# fixture; reproduced here for unit-test independence).  The stripped
# form (after _strip_preprocessor_stop_condition) is the input checked
# in K4 — this catches any regression where the new patterns would
# fire on AF_7mjs's mostly-descriptive preprocessed advice.
AF7MJS_FIXTURE = (
    "1. **Input Files Found**: 7mjs_23883_H.fa, 7mjs_23883_H_1.ccp4, "
    "7mjs_23883_H_2.ccp4, 7mjs_23883_H_full.ccp4, 7mjs_23883_H.pdb\n\n"
    "2. **Experiment Type**: cryo-EM (PredictAndBuild with half-maps)\n\n"
    "3. **Primary Goal**: Run PredictAndBuild for the Fab heavy chain "
    "using the sequence and cryo-EM half-maps, with PDB templates "
    "disabled, to trim/dock the predicted model, rebuild the "
    "mispredicted loop (residues ~100-120), and refine the model.\n\n"
    "4. **Key Parameters**:\n"
    "   - Wavelength: None\n"
    "   - Resolution limit: Not explicitly given (auto-detected from maps)\n"
    "   - Number of expected sites: N/A\n"
    "   - Heavy atom type: N/A\n"
    "   - Space group: N/A\n\n"
    "5. **Program Parameters**:\n"
    "alphafold.use_templates=False\n\n"
    "6. **Special Instructions**:\n"
    "- Input should be the Fab heavy chain sequence "
    "(7mjs_23883_H.fa) and the two half-maps.\n"
    "- In the GUI, uncheck \"Include templates from the PDB\" before "
    "running.\n"
    "- If running offline with supplied predicted models, enable the "
    "\"carry_on\" option.\n\n"
    "7. **Stop Condition**: None\n"
)

# explicit_stop_after_phaser fixture (used in K6; matches the C1 LLM
# test fixture).  Contains the "stop the workflow immediately after
# phaser completes" phrasing that triggered the openai 0/5 failure.
EXPLICIT_STOP_FIXTURE = (
    "1. **Input Files Found**: data.mtz, model.pdb, sequence.fa\n\n"
    "2. **Experiment Type**: X-ray crystallography\n\n"
    "3. **Primary Goal**: Run phenix.phaser for molecular replacement "
    "on the provided data. Stop after running phaser; I will look at "
    "the output before continuing.\n\n"
    "4. **Key Parameters**:\n"
    "   - Wavelength: 1.0\n"
    "   - Resolution limit: 2.5\n"
    "   - Number of expected sites: N/A\n"
    "   - Heavy atom type: N/A\n"
    "   - Space group: P 21 21 21\n\n"
    "5. **Program Parameters**:\n"
    "phaser.search_method=fast\n\n"
    "6. **Special Instructions**:\n"
    "- Use sequence.fa as the sequence file.\n"
    "- Stop the workflow immediately after phaser completes.\n\n"
    "7. **Stop Condition**: Stop after running phenix.phaser\n"
)


# =====================================================================
# K1, K2 — _is_stop_after_requested should recognize new phrasings
# =====================================================================

def test_K1_stop_the_workflow_immediately_after_completes():
    """v117.3 adds recognition for "stop the workflow ... after X completes"
    phrasing.

    On v117.2 baseline: returns False (the phrasing isn't caught).
    Post-v117.3: returns True (new patterns added)."""
    print("Test: K1_stop_the_workflow_immediately_after_completes")
    advice = "stop the workflow immediately after phaser completes"
    result = de._is_stop_after_requested(advice)
    assert result is True, (
        "Expected True (v117.3 should recognize 'stop the workflow ... "
        "after X completes'), got %r" % result)
    print("  PASS")


def test_K2_X_is_the_last_step():
    """v117.3 adds recognition for "X is the last step" phrasing.

    On v117.2 baseline: returns False.
    Post-v117.3: returns True."""
    print("Test: K2_X_is_the_last_step")
    advice = "phaser is the last step"
    result = de._is_stop_after_requested(advice)
    assert result is True, (
        "Expected True (v117.3 should recognize 'X is the last step'), "
        "got %r" % result)
    print("  PASS")


# =====================================================================
# K4 — AF_7mjs preservation (must NOT fire after v117.3)
# =====================================================================

def test_K4_AF_7mjs_preservation():
    """AF_7mjs's preprocessed advice (after _strip_preprocessor_stop_condition)
    must continue to return False from _is_stop_after_requested.  If
    this fires, the v117.1 grounding-exemption would incorrectly keep
    a fabricated after_program on AF_7mjs, breaking the v117.1 K2
    test and the AF_7mjs production behavior.

    Pre-v117.3 baseline: False (correct).
    Post-v117.3: still False (preservation)."""
    print("Test: K4_AF_7mjs_preservation")
    stripped = de._strip_preprocessor_stop_condition(AF7MJS_FIXTURE)
    result = de._is_stop_after_requested(stripped)
    assert result is False, (
        "AF_7mjs preservation broken: _is_stop_after_requested fired "
        "on AF_7mjs's stripped advice.  This would cascade into the "
        "v117.1 grounding-exemption keeping fabricated after_program "
        "for AF_7mjs.  Got %r" % result)
    print("  PASS")


# =====================================================================
# K5 — Descriptive phrasing without user-stop intent
# =====================================================================

def test_K5_descriptive_after_completes_no_intent():
    """The phrasing "X runs after Y completes" describes a sequence
    of operations but does NOT express user-stop intent.

    During v117.3 pre-implementation review (step 2), Option 3 was
    adopted: the over-permissive patterns r'\\bafter\\s+\\w+\\s+completes?\\b'
    and r'\\bafter\\s+\\w+\\s+(?:finishes|is\\s+done)\\b' were dropped
    from the v117.3 plan.  The "after X completes" / "after X
    finishes" phrasings are recognized only via the imperative-marker
    check (window-bounded near program name) and the LLM prompt
    schema documentation, not via global regex.

    Pre-v117.3 baseline: False.
    Post-v117.3: still False (no new regex matches this descriptive
    phrasing)."""
    print("Test: K5_descriptive_after_completes_no_intent")
    advice = "validation runs after refinement completes"
    result = de._is_stop_after_requested(advice)
    assert result is False, (
        "K5 false-positive guard broken: descriptive 'validation runs "
        "after refinement completes' should NOT be classified as user-"
        "stop intent.  Got %r.  This means a v117.3 regex pattern is "
        "too permissive — see Plan §3 Change 2 for the Option 3 "
        "constraints." % result)
    print("  PASS")


# =====================================================================
# K6 — _imperative_marker_nearby on explicit_stop_after_phaser fixture
# =====================================================================

def test_K6_imperative_marker_in_explicit_stop_fixture():
    """The explicit_stop_after_phaser fixture (the C1 LLM test that
    was failing 0/5 on openai) contains the phrase "stop the workflow
    immediately after phaser completes" in its Special Instructions
    section.  Pre-v117.3, _imperative_marker_nearby returned False
    because none of the 15 _IMPERATIVE_STOP_MARKERS matched within
    the 300-char window of "phaser".  Post-v117.3, "stop the workflow"
    is in the markers list and DOES match.

    This is the unit-test correlate of the LLM test fixing.

    Pre-v117.3 baseline: False (gap).
    Post-v117.3: True (new marker)."""
    print("Test: K6_imperative_marker_in_explicit_stop_fixture")
    stripped = de._strip_preprocessor_stop_condition(EXPLICIT_STOP_FIXTURE)
    variants = de._program_name_variants("phenix.phaser")
    _variant, position = de._find_variant_in_text(
        variants, stripped, lambda x: None)
    assert position is not None, (
        "Test fixture sanity: 'phaser' must appear in the stripped "
        "advice.  If this fails, the fixture has changed.")
    result = de._imperative_marker_nearby(stripped, position)
    assert result is True, (
        "Expected True (v117.3 marker 'stop the workflow' should "
        "match within 300 chars of 'phaser' position %d), got %r"
        % (position, result))
    print("  PASS")


# =====================================================================
# K7a, K7b — Existing phrasings still recognized
# =====================================================================

def test_K7a_density_modify_and_stop_preserved():
    """The C1 fixture "density modify and stop" must continue to
    return True post-v117.3.  Catches regressions where the new
    patterns/markers might shadow or interact with existing patterns."""
    print("Test: K7a_density_modify_and_stop_preserved")
    result = de._is_stop_after_requested("density modify and stop")
    assert result is True, (
        "Existing pattern '\\band\\s+stop\\b' broken: 'density modify "
        "and stop' should return True, got %r" % result)
    print("  PASS")


def test_K7b_refine_and_stop_preserved():
    """The C1 fixture "refine and stop" must continue to return True
    post-v117.3.  Same justification as K7a."""
    print("Test: K7b_refine_and_stop_preserved")
    result = de._is_stop_after_requested("refine and stop")
    assert result is True, (
        "Existing pattern '\\band\\s+stop\\b' broken: 'refine and stop' "
        "should return True, got %r" % result)
    print("  PASS")


# =====================================================================
# K8a, K8b — Existing negation cases still recognized
# =====================================================================

def test_K8a_negation_preserved():
    """The negation pattern _NEGATIVE_STOP_PATTERN must continue to
    suppress matches on inputs like "do not stop"."""
    print("Test: K8a_negation_preserved")
    result = de._is_stop_after_requested("do not stop after refinement")
    assert result is False, (
        "Negation guard broken: 'do not stop after refinement' should "
        "return False, got %r" % result)
    print("  PASS")


def test_K8b_stop_condition_none_preserved():
    """The _STOP_CONDITION_NONE pattern must continue to recognize
    "Stop Condition: None" as the absence of stop intent."""
    print("Test: K8b_stop_condition_none_preserved")
    result = de._is_stop_after_requested("**Stop Condition**: None")
    assert result is False, (
        "Stop Condition: None guard broken: should return False, "
        "got %r" % result)
    print("  PASS")


def run_all_tests():
    tests = [
        test_K1_stop_the_workflow_immediately_after_completes,
        test_K2_X_is_the_last_step,
        test_K4_AF_7mjs_preservation,
        test_K5_descriptive_after_completes_no_intent,
        test_K6_imperative_marker_in_explicit_stop_fixture,
        test_K7a_density_modify_and_stop_preserved,
        test_K7b_refine_and_stop_preserved,
        test_K8a_negation_preserved,
        test_K8b_stop_condition_none_preserved,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            print("  FAIL: %s" % e)
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
