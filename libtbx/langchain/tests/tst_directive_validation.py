"""K_H5_1_1: Directive-extractor validation cleanups (v119.H5.1.1).

Covers two latent bugs in agent/directive_extractor.py surfaced
during the §2.1 v118.9 leftovers review:

  §A   merge_directives _corrected_from sidecar protection (7 tests)
  §B   validate_directives boolean list-wrap defense (16 tests)

Total: 23 tests.

Sandbox-friendly: agent.directive_extractor imports cleanly
without libtbx (no group_args / Program dependency in the
two functions exercised here).  All tests should PASS both
in sandbox and under PHENIX.

Item 3 (§A) verifies:
  - merge_directives in agent/directive_extractor.py
  - Production callers: phenix_ai/run_ai_analysis.py:1276
    (engine wrapper) and programs/ai_agent.py:8376 (agent
    main path).  Both inherit the fix automatically since
    merge_directives is the shared helper.

Item 4 (§B) verifies:
  - validate_directives in agent/directive_extractor.py
  - Two sites: file_preferences booleans (prefer_anomalous,
    prefer_unmerged, prefer_merged) and workflow_preferences
    booleans (use_experimental_phasing,
    use_molecular_replacement, use_mr_sad, model_is_placed,
    wants_validation_only).
"""
from __future__ import absolute_import, division, print_function

import copy
import os
import sys


# =====================================================================
# Import helpers — sandbox-friendly
# =====================================================================

def _try_import_de_funcs():
    """Import merge_directives and validate_directives.

    agent/directive_extractor.py is sandbox-importable because
    merge_directives and validate_directives don't transitively
    require libtbx.  Returns ((merge_fn, validate_fn), None) on
    success, or ((None, None), error_msg) on failure.
    """
    # Prefer PHENIX-style import path first
    try:
        from libtbx.langchain.agent.directive_extractor import (
            merge_directives, validate_directives)
        return (merge_directives, validate_directives), None
    except ImportError:
        pass
    # Sandbox path
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from agent.directive_extractor import (
            merge_directives, validate_directives)
        return (merge_directives, validate_directives), None
    except ImportError as e:
        return (None, None), str(e)


# =====================================================================
# §A: merge_directives _corrected_from sidecar protection (7 tests)
# =====================================================================
#
# The _corrected_from sidecar in stop_conditions describes the
# CURRENT value of after_program after an
# experiment-type-mismatch correction.  Before H5.1.1, the
# default override-wins dict-merge had two pathologies:
#
#   1. If override.after_program == _corrected_from.from, the
#      merge silently REVERTED the correction (the LLM had been
#      corrected, the simple-extractor re-introduced the wrong
#      value, override won).
#   2. If override.after_program was a third (different) value,
#      the sidecar's `to` field became stale — pointing at a
#      transition that no longer described the merged
#      after_program.  ("Zombie metadata.")
#
# Fix: detect both cases at merge time.  Strip in case 1
# (preserve correction); blanket-wipe sidecar in case 2
# (clear stale tracking).


def _build_corrected_base(
        after_program="phenix.resolve_cryo_em",
        from_value="phenix.autobuild_denmod",
        extra=None):
    """Build a stop_conditions dict mimicking a corrected
    after_program with sidecar."""
    sc = {
        "after_program": after_program,
        "_corrected_from": {
            "from": from_value,
            "to": after_program,
            "reason": "experiment_type_mismatch",
            "experiment_type": "cryoem",
        },
    }
    if extra:
        sc.update(extra)
    return {"stop_conditions": sc}


def test_merge_strip_when_override_reintroduces_corrected_from():
    """Original bug case: override's after_program equals the
    corrected-away value.  Result must keep base's corrected
    after_program AND base's sidecar."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {
        "after_program": "phenix.autobuild_denmod"}}

    result = merge_fn(base, override)
    sc = result["stop_conditions"]

    assert sc["after_program"] == "phenix.resolve_cryo_em", (
        "Correction was reverted — expected "
        "phenix.resolve_cryo_em, got %r" % sc["after_program"])
    assert "_corrected_from" in sc, (
        "Sidecar was cleared in strip branch — should survive")
    assert sc["_corrected_from"]["from"] == "phenix.autobuild_denmod"
    assert sc["_corrected_from"]["to"] == "phenix.resolve_cryo_em"
    print("  PASS: test_merge_strip_when_override_reintroduces_corrected_from")


def test_merge_wipe_when_override_sets_third_value():
    """Zombie-metadata case: override sets after_program to a
    DIFFERENT value (not matching _corrected_from.from).
    Result must use override's value AND clear the now-stale
    sidecar."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {
        "after_program": "phenix.refine"}}  # Third value

    result = merge_fn(base, override)
    sc = result["stop_conditions"]

    assert sc["after_program"] == "phenix.refine", (
        "Override should win — expected phenix.refine, got %r"
        % sc["after_program"])
    assert "_corrected_from" not in sc, (
        "Stale sidecar should have been cleared (zombie metadata "
        "prevention); got %r" % sc.get("_corrected_from"))
    print("  PASS: test_merge_wipe_when_override_sets_third_value")


def test_merge_override_sidecar_takes_precedence():
    """If override brings its own _corrected_from, that's the
    NEW sidecar — base's is replaced via standard dict-merge."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {
        "after_program": "phenix.refine",
        "_corrected_from": {
            "from": "phenix.autobuild",
            "to": "phenix.refine",
            "reason": "test_override_sidecar",
            "experiment_type": "xray",
        }}}

    result = merge_fn(base, override)
    sc = result["stop_conditions"]

    assert sc["after_program"] == "phenix.refine"
    assert "_corrected_from" in sc
    # Override's sidecar should win — standard dict-merge
    assert sc["_corrected_from"]["reason"] == "test_override_sidecar", (
        "Override's sidecar should replace base's; got reason=%r"
        % sc["_corrected_from"].get("reason"))
    print("  PASS: test_merge_override_sidecar_takes_precedence")


def test_merge_no_sidecar_passthrough():
    """Edge: base has after_program but no sidecar (e.g., LLM
    never needed correction).  Override changes it.  Standard
    dict-merge behaviour expected — no new logic should fire."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = {"stop_conditions": {"after_program": "phenix.autosol"}}
    override = {"stop_conditions": {"after_program": "phenix.refine"}}

    result = merge_fn(base, override)
    sc = result["stop_conditions"]
    assert sc["after_program"] == "phenix.refine"
    assert "_corrected_from" not in sc
    print("  PASS: test_merge_no_sidecar_passthrough")


def test_merge_override_doesnt_touch_after_program():
    """If override only touches after_cycle (not after_program),
    base's corrected after_program AND sidecar survive
    unchanged, override's after_cycle merges in normally."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {"after_cycle": 5}}

    result = merge_fn(base, override)
    sc = result["stop_conditions"]
    assert sc["after_program"] == "phenix.resolve_cryo_em"
    assert "_corrected_from" in sc
    assert sc["after_cycle"] == 5
    print("  PASS: test_merge_override_doesnt_touch_after_program")


def test_merge_no_mutation_of_inputs():
    """Invariant check: the dict-comprehension rebuilds in the
    new code must not mutate the caller's input dicts.
    Without this, repeated calls or shared references would
    silently corrupt state."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {
        "after_program": "phenix.autobuild_denmod"}}

    # Take deep snapshots BEFORE the call so we can compare.
    base_snapshot = copy.deepcopy(base)
    override_snapshot = copy.deepcopy(override)

    _ = merge_fn(base, override)

    assert base == base_snapshot, (
        "merge_directives mutated 'base' input!\n"
        "Before: %r\n"
        "After:  %r" % (base_snapshot, base))
    assert override == override_snapshot, (
        "merge_directives mutated 'override' input!\n"
        "Before: %r\n"
        "After:  %r" % (override_snapshot, override))

    # Also exercise the wipe branch for the same invariant
    override2 = {"stop_conditions": {"after_program": "phenix.refine"}}
    override2_snapshot = copy.deepcopy(override2)
    base2 = _build_corrected_base()
    base2_snapshot = copy.deepcopy(base2)

    _ = merge_fn(base2, override2)
    assert base2 == base2_snapshot, (
        "wipe branch mutated base!  Before: %r  After: %r"
        % (base2_snapshot, base2))
    assert override2 == override2_snapshot, (
        "wipe branch mutated override!  Before: %r  After: %r"
        % (override2_snapshot, override2))
    print("  PASS: test_merge_no_mutation_of_inputs")


def test_merge_strip_preserves_other_override_keys():
    """When the strip branch fires (override's after_program
    matches _corrected_from.from), other keys in override's
    stop_conditions still merge in normally — only
    after_program is stripped."""
    (merge_fn, _), err = _try_import_de_funcs()
    if merge_fn is None:
        print("  SKIP: cannot import merge_directives (%s)" % err)
        return

    base = _build_corrected_base()
    override = {"stop_conditions": {
        "after_program": "phenix.autobuild_denmod",  # Will be stripped
        "after_cycle": 10,                            # Should survive
        "max_refine_cycles": 5,                       # Should survive
    }}

    result = merge_fn(base, override)
    sc = result["stop_conditions"]
    assert sc["after_program"] == "phenix.resolve_cryo_em"
    assert sc["after_cycle"] == 10
    assert sc["max_refine_cycles"] == 5
    assert "_corrected_from" in sc
    print("  PASS: test_merge_strip_preserves_other_override_keys")


# =====================================================================
# §B: validate_directives boolean list-wrap defense (16 tests)
# =====================================================================
#
# bool([False]) returns True in Python — a non-empty list is
# truthy regardless of contents.  If the LLM emits a
# list-wrapped boolean ([true] or [false]), the bare bool()
# coercion in validate_directives would silently flip [False]
# to True.
#
# Fix: explicit type-gated unwrap.  Only unwrap [bool] (not
# [int] or [str]).  Non-bool list elements get dropped + logged
# rather than silently coerced.
#
# Eight input shapes × 2 keys (one from each of the two boolean
# blocks: prefer_anomalous in file_preferences, use_mr_sad in
# workflow_preferences).


def _validate_pref_bool(validate_fn, key, value):
    """Helper: run a single-key file_preferences validation and
    return the validated value (or None if dropped)."""
    d = {"file_preferences": {key: value}}
    result = validate_fn(d)
    fp = result.get("file_preferences", {})
    return fp.get(key, None)


def _validate_wf_bool(validate_fn, key, value):
    """Helper: run a single-key workflow_preferences validation."""
    d = {"workflow_preferences": {key: value}}
    result = validate_fn(d)
    wf = result.get("workflow_preferences", {})
    return wf.get(key, None)


def test_validate_bool_bare_true_unchanged():
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", True) is True
    assert _validate_wf_bool(validate_fn, "use_mr_sad", True) is True
    print("  PASS: test_validate_bool_bare_true_unchanged")


def test_validate_bool_bare_false_unchanged():
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", False) is False
    assert _validate_wf_bool(validate_fn, "use_mr_sad", False) is False
    print("  PASS: test_validate_bool_bare_false_unchanged")


def test_validate_bool_list_wrapped_false():
    """The headline fix: [False] must become False, not True."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", [False]) is False, (
        "[False] should unwrap to False — pre-fix bug returned True")
    assert _validate_wf_bool(validate_fn, "use_mr_sad", [False]) is False
    print("  PASS: test_validate_bool_list_wrapped_false")


def test_validate_bool_list_wrapped_true():
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", [True]) is True
    assert _validate_wf_bool(validate_fn, "use_mr_sad", [True]) is True
    print("  PASS: test_validate_bool_list_wrapped_true")


def test_validate_bool_empty_list_dropped():
    """Empty list — invalid, key dropped from output."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", []) is None, (
        "Empty list should drop the key, not produce a False")
    assert _validate_wf_bool(validate_fn, "use_mr_sad", []) is None
    print("  PASS: test_validate_bool_empty_list_dropped")


def test_validate_bool_multi_element_list_dropped():
    """Multi-element list — invalid, key dropped from output."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", [True, False]) is None
    assert _validate_wf_bool(validate_fn, "use_mr_sad", [False, True]) is None
    print("  PASS: test_validate_bool_multi_element_list_dropped")


def test_validate_bool_list_with_non_bool_element_dropped():
    """[1] should NOT unwrap (type-gate rejects).  Even though
    int(1) would coerce to True, the type-gate requires
    isinstance(v[0], bool).  This is stricter than
    _coerce_setting_value would have been — deliberate per the
    Gemini critique."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", [1]) is None, (
        "[1] (int in list) should drop, not unwrap — type-gate "
        "must reject non-bool inner elements")
    assert _validate_wf_bool(validate_fn, "use_mr_sad", ["true"]) is None
    print("  PASS: test_validate_bool_list_with_non_bool_element_dropped")


def test_validate_bool_bare_int_preserved():
    """Bare int 0/1 — existing behaviour preserved (the fix is
    additive; non-list paths unchanged)."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", 0) is False
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", 1) is True
    assert _validate_wf_bool(validate_fn, "use_mr_sad", 0) is False
    assert _validate_wf_bool(validate_fn, "use_mr_sad", 1) is True
    print("  PASS: test_validate_bool_bare_int_preserved")


def test_validate_bool_bare_string_legacy_quirk():
    """Bare string preserved including the legacy quirk that
    bool('false') == True.  This is INTENTIONALLY not fixed
    by H5.1.1 (separate bug class — would balloon scope).
    The test pins the existing behaviour so a future fix
    catches it deliberately rather than accidentally."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    # Both 'true' and 'false' come out True under bool() — this
    # is the legacy quirk.
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", "true") is True
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", "false") is True
    # Empty string is the only falsy string for bool() — comes
    # out False (same as pre-fix; matches the str branch which
    # still applies bool() to non-list-wrapped strings).
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", "") is False, (
        "Empty string should produce False (existing bool('') "
        "behaviour); if H5.1.1 changed this, behaviour drift")
    print("  PASS: test_validate_bool_bare_string_legacy_quirk")


def test_validate_bool_dict_input_dropped():
    """Non-list, non-bool, non-int, non-str — drop with log."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous",
                                 {"k": "v"}) is None
    assert _validate_wf_bool(validate_fn, "use_mr_sad",
                              {"k": "v"}) is None
    print("  PASS: test_validate_bool_dict_input_dropped")


def test_validate_bool_none_input_dropped():
    """None — drop with log.  Pre-fix bool(None) == False
    would have silently produced a False value; the fix is
    stricter and drops."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    assert _validate_pref_bool(validate_fn, "prefer_anomalous", None) is None
    assert _validate_wf_bool(validate_fn, "use_mr_sad", None) is None
    print("  PASS: test_validate_bool_none_input_dropped")


def test_validate_bool_other_pref_keys_share_fix():
    """All three file_preferences boolean keys share the fix."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    for key in ("prefer_anomalous", "prefer_unmerged", "prefer_merged"):
        assert _validate_pref_bool(validate_fn, key, [False]) is False, (
            "key=%s did not get the list-wrap fix" % key)
        assert _validate_pref_bool(validate_fn, key, [True]) is True
    print("  PASS: test_validate_bool_other_pref_keys_share_fix")


def test_validate_bool_other_wf_keys_share_fix():
    """All five workflow_preferences boolean keys share the fix."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    for key in ("use_experimental_phasing", "use_molecular_replacement",
                "use_mr_sad", "model_is_placed", "wants_validation_only"):
        assert _validate_wf_bool(validate_fn, key, [False]) is False, (
            "key=%s did not get the list-wrap fix" % key)
        assert _validate_wf_bool(validate_fn, key, [True]) is True
    print("  PASS: test_validate_bool_other_wf_keys_share_fix")


def test_validate_non_bool_keys_unaffected():
    """Item 4 regression check: list-typed key (`exclude` in
    file_preferences) MUST keep its list — the fix only
    touches the boolean blocks."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    d = {"file_preferences": {"exclude": ["bad1.mtz", "bad2.mtz"]}}
    result = validate_fn(d)
    fp = result.get("file_preferences", {})
    assert fp.get("exclude") == ["bad1.mtz", "bad2.mtz"], (
        "exclude (list-typed) was flattened or modified — fix "
        "leaked outside boolean blocks.  Got: %r" % fp.get("exclude"))
    print("  PASS: test_validate_non_bool_keys_unaffected")


def test_validate_mixed_bool_keys_in_one_dict():
    """Realistic case: a directives dict with multiple
    boolean keys at once.  All get the fix; non-boolean keys
    are untouched."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    d = {
        "file_preferences": {
            "prefer_anomalous": [False],   # list-wrap: fix
            "prefer_unmerged": True,         # bare: unchanged
            "exclude": ["x.mtz"],            # list-typed: untouched
        },
        "workflow_preferences": {
            "use_mr_sad": [True],            # list-wrap: fix
            "use_experimental_phasing": False,  # bare: unchanged
        },
    }
    result = validate_fn(d)
    fp = result.get("file_preferences", {})
    wf = result.get("workflow_preferences", {})

    assert fp.get("prefer_anomalous") is False
    assert fp.get("prefer_unmerged") is True
    assert fp.get("exclude") == ["x.mtz"]
    assert wf.get("use_mr_sad") is True
    assert wf.get("use_experimental_phasing") is False
    print("  PASS: test_validate_mixed_bool_keys_in_one_dict")


def test_validate_no_mutation_of_input():
    """Invariant: validate_directives must not mutate its
    input dict."""
    (_, validate_fn), err = _try_import_de_funcs()
    if validate_fn is None:
        print("  SKIP: cannot import validate_directives (%s)" % err)
        return
    d = {
        "file_preferences": {"prefer_anomalous": [False]},
        "workflow_preferences": {"use_mr_sad": [True]},
    }
    snapshot = copy.deepcopy(d)
    _ = validate_fn(d)
    assert d == snapshot, (
        "validate_directives mutated input!  Before: %r  After: %r"
        % (snapshot, d))
    print("  PASS: test_validate_no_mutation_of_input")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: merge_directives sidecar protection (7)
    test_merge_strip_when_override_reintroduces_corrected_from()
    test_merge_wipe_when_override_sets_third_value()
    test_merge_override_sidecar_takes_precedence()
    test_merge_no_sidecar_passthrough()
    test_merge_override_doesnt_touch_after_program()
    test_merge_no_mutation_of_inputs()
    test_merge_strip_preserves_other_override_keys()

    # §B: validate_directives boolean list-wrap defense (16)
    test_validate_bool_bare_true_unchanged()
    test_validate_bool_bare_false_unchanged()
    test_validate_bool_list_wrapped_false()
    test_validate_bool_list_wrapped_true()
    test_validate_bool_empty_list_dropped()
    test_validate_bool_multi_element_list_dropped()
    test_validate_bool_list_with_non_bool_element_dropped()
    test_validate_bool_bare_int_preserved()
    test_validate_bool_bare_string_legacy_quirk()
    test_validate_bool_dict_input_dropped()
    test_validate_bool_none_input_dropped()
    test_validate_bool_other_pref_keys_share_fix()
    test_validate_bool_other_wf_keys_share_fix()
    test_validate_non_bool_keys_unaffected()
    test_validate_mixed_bool_keys_in_one_dict()
    test_validate_no_mutation_of_input()


if __name__ == "__main__":
    print("K_H5_1_1: Directive-extractor validation cleanups (v119.H5.1.1)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H5_1_1 complete.")
