"""Tests for v116.10 Phase 4a improvements to directive_validator.

Covers two changes:

  1. `_suggest_similar_programs` upgraded from substring+overlap
     heuristic to difflib edit-distance.  Realistic typos now produce
     correct suggestions; unrelated names produce no suggestions.

  2. The `after_program` typo check in `validate_directives` now
     includes "Did you mean: X?" hints in its issue messages,
     matching the user-advice-text branch.

Each test_* function is one test case, asserts inline, prints PASS
on success.  Follows the same fail-fast convention as
tst_skip_to_program.py and tst_general_resolver.py.
"""

from __future__ import absolute_import, division, print_function

import os
import sys


# --- Path setup ------------------------------------------------------------
# Allow the test to find directive_validator regardless of where it's run from.
_HERE = os.path.dirname(os.path.abspath(__file__))
_AGENT_DIR = os.path.normpath(os.path.join(_HERE, "..", "agent"))
if _AGENT_DIR not in sys.path:
    sys.path.insert(0, _AGENT_DIR)

try:
    from libtbx.langchain.agent.directive_validator import (
        _suggest_similar_programs, validate_directives,
    )
except ImportError:
    from directive_validator import (
        _suggest_similar_programs, validate_directives,
    )


# --- Test fixtures ---------------------------------------------------------
# Fixed available_programs dict — represents what the production registry
# returns.  Defined here so tests don't depend on the live YAML.
_AVAILABLE_PROGRAMS = {
    "phenix.refine": {},
    "phenix.real_space_refine": {},
    "phenix.phaser": {},
    "phenix.predict_and_build": {},
    "phenix.xtriage": {},
    "phenix.mtriage": {},
    "phenix.autosol": {},
    "phenix.autobuild": {},
    "phenix.molprobity": {},
    "phenix.model_vs_data": {},
    "phenix.map_correlations": {},
    "phenix.ligandfit": {},
    "phenix.polder": {},
    "phenix.process_predicted_model": {},
    "phenix.dock_in_map": {},
    "phenix.map_sharpening": {},
    "phenix.validation_cryoem": {},
    "phenix.pdbtools": {},
    "phenix.map_symmetry": {},
}
_AVAILABLE_SET = set(_AVAILABLE_PROGRAMS.keys())


# =====================================================================
# SECTION A: _suggest_similar_programs algorithm
# =====================================================================

def test_suggest_typos_with_one_char_change():
    """Common typos (missing/extra/transposed letter) suggest the
    correct program name."""
    print("Test: suggest_typos_with_one_char_change")
    cases = [
        ("phenix.predict_and_buld",  "phenix.predict_and_build"),
        ("phenix.predikt_and_build", "phenix.predict_and_build"),
        ("phenix.refien",            "phenix.refine"),
        ("phenix.xtraige",           "phenix.xtriage"),
        ("phenix.autosool",          "phenix.autosol"),
        ("phenix.molprobitiy",       "phenix.molprobity"),
        ("phenix.autobild",          "phenix.autobuild"),
    ]
    for typo, expected in cases:
        suggestions = _suggest_similar_programs(typo, _AVAILABLE_SET)
        assert expected in suggestions, \
            "%s should suggest %s, got %s" % (
                typo, expected, suggestions)
    print("  PASS")


def test_suggest_exact_match_included():
    """An exact match is returned in the suggestion list."""
    print("Test: suggest_exact_match_included")
    suggestions = _suggest_similar_programs(
        "phenix.refine", _AVAILABLE_SET)
    assert "phenix.refine" in suggestions, \
        "Exact match should be returned, got %s" % suggestions
    print("  PASS")


def test_suggest_gibberish_returns_empty():
    """Names with no meaningful similarity get no suggestions.

    The old algorithm over-suggested based on shared characters;
    difflib's cutoff correctly rejects these.
    """
    print("Test: suggest_gibberish_returns_empty")
    gibberish_cases = [
        "phenix.zxqwerty",
        "phenix.completely_made_up",
        "phenix.foo_bar_baz",
    ]
    for name in gibberish_cases:
        suggestions = _suggest_similar_programs(name, _AVAILABLE_SET)
        assert len(suggestions) == 0, \
            "%s should return [], got %s" % (name, suggestions)
    print("  PASS")


def test_suggest_caps_at_three_results():
    """Never returns more than 3 suggestions even when many are close."""
    print("Test: suggest_caps_at_three_results")
    # 'map' is a prefix of several programs
    suggestions = _suggest_similar_programs(
        "phenix.map", _AVAILABLE_SET)
    assert len(suggestions) <= 3, \
        "Should return at most 3 suggestions, got %d: %s" % (
            len(suggestions), suggestions)
    print("  PASS")


def test_suggest_works_without_phenix_prefix():
    """Bare program names (no 'phenix.' prefix) are accepted."""
    print("Test: suggest_works_without_phenix_prefix")
    suggestions = _suggest_similar_programs("refien", _AVAILABLE_SET)
    assert "phenix.refine" in suggestions, \
        "Bare 'refien' should suggest 'phenix.refine', got %s" % \
        suggestions
    print("  PASS")


def test_suggest_returns_canonical_form():
    """Suggestions come back as 'phenix.NAME' regardless of input form."""
    print("Test: suggest_returns_canonical_form")
    suggestions = _suggest_similar_programs(
        "phenix.refien", _AVAILABLE_SET)
    for s in suggestions:
        assert s.startswith("phenix."), \
            "Suggestion %r should start with 'phenix.'" % s
    print("  PASS")


# =====================================================================
# SECTION B: validate_directives with after_program suggestions
# =====================================================================

def test_validate_typo_after_program_returns_invalid():
    """A typo'd after_program produces valid=False."""
    print("Test: validate_typo_after_program_returns_invalid")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.predict_and_buld"},
    }
    result = validate_directives(
        user_advice="",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert result.valid is False, \
        "Typo should produce valid=False, got %r" % result.valid
    print("  PASS")


def test_validate_typo_issue_includes_did_you_mean():
    """The issue message for a typo'd after_program names the
    correct program."""
    print("Test: validate_typo_issue_includes_did_you_mean")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.predict_and_buld"},
    }
    result = validate_directives(
        user_advice="",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    typo_issue = next(
        (i for i in result.issues if "predict_and_buld" in i),
        None,
    )
    assert typo_issue is not None, \
        "Expected issue mentioning the typo, got %s" % result.issues
    assert "Did you mean" in typo_issue, \
        "Issue should include 'Did you mean': %r" % typo_issue
    assert "phenix.predict_and_build" in typo_issue, \
        "Issue should suggest the correct program: %r" % typo_issue
    print("  PASS")


def test_validate_gibberish_no_misleading_suggestion():
    """A gibberish after_program is rejected without a misleading hint."""
    print("Test: validate_gibberish_no_misleading_suggestion")
    directives = {
        "stop_conditions": {"after_program": "phenix.zxqwerty"},
    }
    result = validate_directives(
        user_advice="",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert result.valid is False
    gibberish_issue = next(
        (i for i in result.issues if "zxqwerty" in i), None)
    assert gibberish_issue is not None
    assert "Did you mean" not in gibberish_issue, \
        "Gibberish should not suggest anything: %r" % gibberish_issue
    print("  PASS")


def test_validate_correct_after_program_passes():
    """A correctly-spelled after_program produces valid=True with
    no issues."""
    print("Test: validate_correct_after_program_passes")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.predict_and_build"},
    }
    result = validate_directives(
        user_advice="",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert result.valid is True, \
        "Correct name should pass, got valid=%r issues=%s" % (
            result.valid, result.issues)
    assert not any(
        "predict_and_build" in i for i in result.issues), \
        "No issues should mention the program, got %s" % \
        result.issues
    print("  PASS")


def test_validate_message_contains_suggestion():
    """The user-facing message text includes the suggestion (this is
    what gets printed via vlog.quiet)."""
    print("Test: validate_message_contains_suggestion")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.predict_and_buld"},
    }
    result = validate_directives(
        user_advice="",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert "phenix.predict_and_build" in (result.message or ""), \
        "Message should contain the suggestion, got %r" % \
        result.message
    print("  PASS")


# =====================================================================
# SECTION C: End-to-end realistic scenarios
# =====================================================================

def test_predict_and_stop_typo_blocks_cleanly():
    """A 'predict and stop' request with a typo'd program name
    fails validation with a helpful hint."""
    print("Test: predict_and_stop_typo_blocks_cleanly")
    directives = {
        "intent": "task",
        "stop_conditions": {
            "after_program": "phenix.predict_and_buld",
            "skip_validation": True,
        },
        "program_settings": {
            "phenix.predict_and_buld": {
                "stop_after_predict": True},
        },
    }
    result = validate_directives(
        user_advice="predict and stop",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert result.valid is False
    assert "predict_and_buld" in (result.message or ""), \
        "Message should name the typo, got %r" % result.message
    assert "phenix.predict_and_build" in (result.message or ""), \
        "Message should suggest the correct name, got %r" % \
        result.message
    print("  PASS")


def test_predict_and_stop_correct_spelling_passes():
    """A correctly-spelled 'predict and stop' passes validation."""
    print("Test: predict_and_stop_correct_spelling_passes")
    directives = {
        "intent": "task",
        "stop_conditions": {
            "after_program": "phenix.predict_and_build",
            "skip_validation": True,
        },
        "program_settings": {
            "phenix.predict_and_build": {
                "stop_after_predict": True},
        },
    }
    result = validate_directives(
        user_advice="predict and stop",
        directives=directives,
        available_programs=_AVAILABLE_PROGRAMS,
    )
    assert result.valid is True, \
        "Correct spelling should pass, got issues=%s" % \
        result.issues
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """Run all tests using cctbx fail-fast convention."""
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        try:
            from tests.tst_utils import (
                run_tests_with_fail_fast)
        except ImportError:
            # Fall through to standalone runner
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    """Simple runner when tst_utils is not available."""
    test_fns = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
