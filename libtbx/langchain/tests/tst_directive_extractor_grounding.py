"""Tests for v116.19a (Option A) grounding guardrail.

Replaces the broader v116.19 test suite.  The Option A guardrail drops
LLM-extracted after_program / prefer_programs only in two situations:

  (1) Pure fabrication — program name absent from advice in any
      canonical spelling variant.  AF_7mjs's
      after_program=phenix.real_space_refine where
      "real_space_refine" never appears in the README is the
      canonical example.

  (2) Preprocessed advice + no global stop intent + no imperative
      marker near the program name.  AF_7mjs's
      after_program=phenix.predict_and_build (had it been extracted)
      would be caught here: "PredictAndBuild" IS in the advice, but
      the user wrote Stop Condition: None and the goal is multi-step.

Bare user input like "run xtriage on my data" passes the guardrail
(program name is present, advice is not preprocessed) so the
directive extractor's documented bare-imperative mappings continue
to work.  This is the deliberate trade-off vs the broader v116.19
prototype, which would have dropped these.
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOTS = (
    os.path.abspath(os.path.join(HERE, "..")),
    os.path.abspath(HERE),
)
for _root in _ROOTS:
    if _root not in sys.path:
        sys.path.insert(0, _root)

if "libtbx" not in sys.modules:
    import types as _types
    for _mod in (
        "libtbx", "libtbx.langchain",
        "libtbx.langchain.agent",
        "libtbx.langchain.knowledge",
        "libtbx.langchain.knowledge.yaml_loader",
        "libtbx.langchain.agent.intent_classifier",
        "libtbx.langchain.agent.program_registry",
    ):
        sys.modules[_mod] = _types.ModuleType(_mod)

    pr = sys.modules["libtbx.langchain.agent.program_registry"]
    class _PR:
        def __init__(self): pass
    pr.ProgramRegistry = _PR

    ic = sys.modules["libtbx.langchain.agent.intent_classifier"]
    ic.classify_intent = lambda a: None


from agent.directive_extractor import (
    _is_program_grounded,
    _program_name_variants,
    _find_variant_in_text,
    _imperative_marker_nearby,
    _advice_came_from_preprocessor,
    _advice_has_explicit_stop_intent,
    _validate_after_program_grounded,
)


def _quiet(msg):
    pass


# =====================================================================
# Helper unit tests
# =====================================================================

def test_variants_real_space_refine():
    print("[test] variants(phenix.real_space_refine)")
    v = _program_name_variants("phenix.real_space_refine")
    assert "phenix.real_space_refine" in v
    assert "real_space_refine" in v
    assert "real space refine" in v
    assert "real-space-refine" in v
    assert "RealSpaceRefine" in v
    print("  PASS — %r" % (v,))


def test_variants_predict_and_build_includes_camel():
    print("[test] variants(phenix.predict_and_build) includes CamelCase")
    v = _program_name_variants("phenix.predict_and_build")
    assert "PredictAndBuild" in v, \
        "CamelCase variant missing from %r" % (v,)
    print("  PASS")


def test_variants_bare_program():
    print("[test] variants('phenix.refine')")
    v = _program_name_variants("phenix.refine")
    assert "phenix.refine" in v
    assert "refine" in v
    print("  PASS")


def test_variants_none_safe():
    print("[test] variants(None) / variants('') are safe")
    assert _program_name_variants(None) == ()
    assert _program_name_variants("") == ()
    assert _program_name_variants(123) == ()
    print("  PASS")


def test_find_variant_word_boundary():
    print("[test] _find_variant_in_text respects word boundaries")
    text = "we will real_space_refine the model"
    variant, pos = _find_variant_in_text(
        _program_name_variants("phenix.real_space_refine"), text)
    assert variant is not None and pos >= 0, \
        "Should match 'real_space_refine' in %r, got %r at %d" % (text, variant, pos)
    print("  PASS — matched %r at pos %d" % (variant, pos))


def test_find_variant_does_not_match_substring():
    print("[test] 'refine' substring of 'real_space_refine' not matched as 'refine'")
    text = "we use real_space_refine to fix it"
    variants_refine = _program_name_variants("phenix.refine")
    variant, pos = _find_variant_in_text(variants_refine, text)
    assert variant is None, \
        "Should not match 'refine' as substring of 'real_space_refine'; " \
        "got variant=%r at %d" % (variant, pos)
    print("  PASS — no false-positive substring match")


def test_imperative_marker_finds_stop_after():
    print("[test] imperative marker 'stop after' detected within window")
    text = "Please stop after running phenix.phaser when done"
    found = _imperative_marker_nearby(text, position=33, window=300)
    assert found
    print("  PASS")


def test_imperative_marker_out_of_window():
    print("[test] imperative marker far away is not detected")
    text = ("stop after running this step.  " + "x" * 500 +
            "  Then do phenix.phaser later.")
    pos = text.find("phaser")
    found = _imperative_marker_nearby(text, position=pos, window=100)
    assert not found
    print("  PASS")


def test_preprocessor_signature_detection():
    print("[test] _advice_came_from_preprocessor detects signatures")
    af_7mjs = ("1. **Input Files Found**: foo.fa\n"
               "2. **Experiment Type**: cryo-EM\n")
    assert _advice_came_from_preprocessor(af_7mjs), \
        "Should detect preprocessor signature"
    raw_advice = "Just run xtriage on my data please."
    assert not _advice_came_from_preprocessor(raw_advice), \
        "Should NOT flag raw user advice as preprocessed"
    print("  PASS")


def test_preprocessor_signature_none_safe():
    print("[test] _advice_came_from_preprocessor handles None/empty")
    assert not _advice_came_from_preprocessor(None)
    assert not _advice_came_from_preprocessor("")
    assert not _advice_came_from_preprocessor(123)
    print("  PASS")


def test_stop_intent_detection():
    print("[test] _advice_has_explicit_stop_intent detects stop language")
    assert _advice_has_explicit_stop_intent("Stop after running phaser")
    assert _advice_has_explicit_stop_intent("Only run xtriage")
    assert _advice_has_explicit_stop_intent("Just run mtriage")
    assert _advice_has_explicit_stop_intent("Run phaser and stop")
    # Negative cases
    assert not _advice_has_explicit_stop_intent("Run xtriage on my data")
    assert not _advice_has_explicit_stop_intent("Rebuild and refine the model")
    print("  PASS")


def test_stop_intent_avoids_substring_false_positives():
    print("[test] stop intent uses word boundaries (no false positives)")
    # "stop after" inside compound: word boundary prevents match
    assert not _advice_has_explicit_stop_intent("nonstop afterburner")
    # "just run" inside "just running into trouble" — substring match
    # would falsely trigger; word boundary prevents
    assert not _advice_has_explicit_stop_intent(
        "We were just running into trouble.")
    print("  PASS")


# =====================================================================
# AF_7mjs regression (the original Bug A class)
# =====================================================================

AF_7MJS_PREPROCESSED = """1. **Input Files Found**: 7mjs_23883_H.fa, 7mjs_23883_H_1.ccp4

2. **Experiment Type**: cryo-EM (PredictAndBuild with half-maps)

3. **Primary Goal**: Run PredictAndBuild for the Fab heavy chain using the sequence and cryo-EM half-maps, with PDB templates disabled, to trim/dock the predicted model, rebuild the mispredicted loop (residues ~100-120), and refine the model.

4. **Key Parameters**:
   - Wavelength: None
   - Resolution limit: Not explicitly given (auto-detected from maps)

5. **Program Parameters**:
alphafold.use_templates=False

6. **Special Instructions**:
- Input should be the Fab heavy chain sequence and the two half-maps.

7. **Stop Condition**: None
"""


def test_af7mjs_drops_fabricated_real_space_refine():
    """Pure fabrication: program name absent from advice."""
    print("[test] AF_7mjs: drops after_program=phenix.real_space_refine "
          "(pure fabrication, name absent)")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.real_space_refine",
        },
        "workflow_preferences": {
            "prefer_programs": ["phenix.real_space_refine"],
        },
    }
    out = _validate_after_program_grounded(
        directives, AF_7MJS_PREPROCESSED, _quiet)
    assert "stop_conditions" not in out, \
        "Pure fabrication should drop. Got %r" % out
    assert "workflow_preferences" not in out, \
        "Pure fabrication should drop. Got %r" % out
    print("  PASS")


def test_af7mjs_drops_fabricated_predict_and_build():
    """Preprocessed + no-stop-intent + no nearby imperative:
    even though PredictAndBuild IS in the advice, the guardrail
    drops it because the rest of the gates fail.
    """
    print("[test] AF_7mjs: drops after_program=phenix.predict_and_build "
          "(name present, preprocessed, no stop intent, no imperative)")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.predict_and_build",
        },
    }
    out = _validate_after_program_grounded(
        directives, AF_7MJS_PREPROCESSED, _quiet)
    assert "stop_conditions" not in out, \
        "Should drop. Got %r" % out
    print("  PASS")


# =====================================================================
# Legitimate after_program is preserved (Option A wins)
# =====================================================================

def test_keeps_explicit_stop_after_phaser():
    """Raw advice with explicit 'Stop after running phaser': kept."""
    print("[test] keeps after_program when 'Stop after running phaser' "
          "is in raw advice")
    directives = {
        "stop_conditions": {"after_program": "phenix.phaser"},
    }
    advice = "Please solve the structure. Stop after running phaser."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert out.get("stop_conditions", {}).get("after_program") == "phenix.phaser"
    print("  PASS")


def test_keeps_only_run_xtriage():
    """Raw advice with 'Only run xtriage, nothing else': kept."""
    print("[test] keeps after_program with 'Only run xtriage, nothing else'")
    directives = {
        "stop_conditions": {"after_program": "phenix.xtriage"},
    }
    advice = "Only run xtriage, nothing else."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert out.get("stop_conditions", {}).get("after_program") == "phenix.xtriage"
    print("  PASS")


def test_keeps_just_do_then_stop():
    """Raw advice with 'just do ... then stop': kept."""
    print("[test] keeps after_program with 'Just do one refinement, "
          "then stop after phenix.refine'")
    directives = {
        "stop_conditions": {"after_program": "phenix.refine"},
    }
    advice = "Just do one refinement, then stop after phenix.refine."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert out.get("stop_conditions", {}).get("after_program") == "phenix.refine"
    print("  PASS")


# =====================================================================
# The Option A win: bare imperatives on RAW advice are kept
# =====================================================================

def test_keeps_bare_run_xtriage_raw():
    """Raw user input 'Run xtriage on my data' — keep (no preprocessor).

    This is the Option A win.  The original (broad) v116.19 would have
    dropped this because "xtriage" appears without "only"/"just"/"stop
    after" nearby.  Option A keeps it because the advice is not
    preprocessed.  The directive-extractor's documented mapping
    "run xtriage" → after_program=phenix.xtriage continues to apply.
    """
    print("[test] [Option A] keeps after_program=phenix.xtriage for "
          "bare 'Run xtriage' (raw, not preprocessed)")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.xtriage",
            "skip_validation": True,
        },
    }
    advice = "Run xtriage on my data."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert out.get("stop_conditions", {}).get(
        "after_program") == "phenix.xtriage", \
        "Bare 'Run xtriage' should be kept; got %r" % out
    print("  PASS")


def test_drops_indirect_reference_no_program_name():
    """Raw 'Check data quality' — the LLM might map to phenix.xtriage
    but the program name doesn't appear in advice.

    This is a known trade-off of Option A's Failure-1 rule (pure
    fabrication = name absent → drop).  The directive-extractor
    prompt documents "check data quality" → after_program=phenix.xtriage
    as a mapping, but Option A drops it because "xtriage" isn't in the
    advice.  The user gets the full workflow (xtriage + downstream)
    instead of stopping after xtriage.  This is a minor inconvenience
    rather than a hard failure; the agent does more than asked rather
    than less, and the user can rephrase as "Only check data quality"
    or "Just run xtriage" to get the explicit stop.

    If this trade-off proves too costly in production, the fix is to
    relax Failure 1 to allow well-known task→program mappings via a
    small lookup table.  For now we accept the trade-off as the price
    of catching AF_7mjs-class fabrications.
    """
    print("[test] [Option A trade-off] 'Check data quality' drops "
          "after_program=phenix.xtriage (name absent from advice)")
    directives = {
        "stop_conditions": {"after_program": "phenix.xtriage"},
    }
    advice = "Check data quality before doing anything."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    # Failure 1 fires — drop
    assert "stop_conditions" not in out, \
        "Indirect reference (program name absent) should drop under "  \
        "Failure 1.  Got %r" % out
    print("  PASS")


def test_keeps_bare_run_mtriage_raw():
    """Raw 'Run mtriage on my map' — keep."""
    print("[test] [Option A] keeps after_program=phenix.mtriage for "
          "bare 'Run mtriage on my map'")
    directives = {
        "stop_conditions": {"after_program": "phenix.mtriage"},
    }
    advice = "Run mtriage on my map."
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert out.get("stop_conditions", {}).get(
        "after_program") == "phenix.mtriage"
    print("  PASS")


# =====================================================================
# Pure fabrication is dropped regardless of preprocessing
# =====================================================================

def test_pure_fabrication_dropped_in_raw_advice():
    """Program name absent from RAW advice — drop.

    This case is the strict Failure-1 leg: if the program name was
    never mentioned, the LLM hallucinated it regardless of source.
    """
    print("[test] pure fabrication dropped even from raw advice")
    directives = {
        "stop_conditions": {"after_program": "phenix.real_space_refine"},
    }
    advice = "Rebuild the loop and refine the model."
    # "real_space_refine" not in advice (only verb "refine" is, which
    # is a different program suffix); even raw advice drops it.
    out = _validate_after_program_grounded(directives, advice, _quiet)
    assert "stop_conditions" not in out, \
        "Pure fabrication should drop even on raw advice. Got %r" % out
    print("  PASS")


# =====================================================================
# Preprocessed + has stop intent → KEEP (don't drop legitimate stops)
# =====================================================================

def test_preprocessed_with_stop_intent_keeps_directive():
    """If preprocessed advice has 'stop after X' anywhere, keep
    after_program=X.  This protects users who DID specify an explicit
    stop via the preprocessor template."""
    print("[test] preprocessed + explicit stop intent: keep after_program")
    advice = ("1. **Input Files Found**: foo.fa\n"
              "2. **Experiment Type**: cryo-EM\n"
              "3. **Primary Goal**: Stop after running phaser to validate.\n")
    directives = {
        "stop_conditions": {"after_program": "phenix.phaser"},
    }
    out = _validate_after_program_grounded(directives, advice, _quiet)
    # Imperative "stop after" is within 300 chars of "phaser" → grounded
    assert out.get("stop_conditions", {}).get(
        "after_program") == "phenix.phaser"
    print("  PASS")


# =====================================================================
# No-op cases (defensive)
# =====================================================================

def test_noop_when_no_after_program():
    print("[test] no-op when directives have no after_program / "
          "prefer_programs")
    directives = {
        "program_settings": {"phenix.refine": {"strategy": "tls"}},
    }
    out = _validate_after_program_grounded(directives, "some advice", _quiet)
    assert out == directives
    print("  PASS")


def test_noop_empty_advice():
    """Empty user_advice: directives unchanged (defensive)."""
    print("[test] empty user_advice: directives returned unchanged")
    directives = {"stop_conditions": {"after_program": "phenix.phaser"}}
    out = _validate_after_program_grounded(directives, "", _quiet)
    assert out == directives
    print("  PASS")


def test_noop_none_directives():
    print("[test] None directives is safe")
    out = _validate_after_program_grounded(None, "some advice", _quiet)
    assert out is None
    print("  PASS")


# =====================================================================
# Section cleanup
# =====================================================================

def test_empty_stop_conditions_section_removed():
    """When dropping the only key in stop_conditions, remove the section."""
    print("[test] empty stop_conditions section removed after drop")
    directives = {
        "stop_conditions": {"after_program": "phenix.real_space_refine"},
    }
    out = _validate_after_program_grounded(
        directives, AF_7MJS_PREPROCESSED, _quiet)
    assert "stop_conditions" not in out
    print("  PASS")


def test_other_stop_conditions_preserved():
    """Other stop_conditions fields preserved when after_program dropped."""
    print("[test] other stop_conditions fields preserved when "
          "after_program dropped")
    directives = {
        "stop_conditions": {
            "after_program": "phenix.real_space_refine",
            "max_refine_cycles": 3,
        },
    }
    out = _validate_after_program_grounded(
        directives, AF_7MJS_PREPROCESSED, _quiet)
    assert "stop_conditions" in out
    assert "after_program" not in out["stop_conditions"]
    assert out["stop_conditions"].get("max_refine_cycles") == 3
    print("  PASS")


def test_prefer_programs_selective_drop():
    """prefer_programs: one grounded, one not — keep grounded only."""
    print("[test] prefer_programs: selective drop (one in advice, one not)")
    directives = {
        "workflow_preferences": {
            "prefer_programs": ["phenix.phaser", "phenix.real_space_refine"],
        },
    }
    advice = ("1. **Input Files Found**: foo.fa\n"
              "2. **Experiment Type**: cryo-EM\n"
              "3. **Primary Goal**: Stop after running phaser to validate.\n")
    out = _validate_after_program_grounded(directives, advice, _quiet)
    kept = out.get("workflow_preferences", {}).get("prefer_programs", [])
    # phaser is grounded (has "stop after" nearby)
    # real_space_refine is pure fabrication (name absent from advice)
    assert "phenix.phaser" in kept
    assert "phenix.real_space_refine" not in kept
    print("  PASS — kept=%r" % kept)


# =====================================================================
# _is_program_grounded direct API
# =====================================================================

def test_is_program_grounded_af7mjs():
    """_is_program_grounded returns False for AF_7mjs's bogus directive."""
    print("[test] _is_program_grounded False for "
          "phenix.real_space_refine on AF_7mjs")
    assert not _is_program_grounded(
        "phenix.real_space_refine", AF_7MJS_PREPROCESSED, _quiet)
    print("  PASS")


def test_is_program_grounded_legitimate_raw():
    """_is_program_grounded returns True for bare 'run xtriage' raw."""
    print("[test] _is_program_grounded True for bare 'Run xtriage' "
          "(raw, not preprocessed)")
    assert _is_program_grounded(
        "phenix.xtriage", "Run xtriage on my data.", _quiet)
    print("  PASS")


def test_is_program_grounded_legitimate_stop_after():
    """_is_program_grounded True for 'Stop after phaser'."""
    print("[test] _is_program_grounded True for 'Stop after running phaser'")
    assert _is_program_grounded(
        "phenix.phaser", "Please stop after running phaser.", _quiet)
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    try:
        from libtbx.langchain.tests.tst_utils import run_tests_with_fail_fast
    except ImportError:
        try:
            from tests.tst_utils import run_tests_with_fail_fast
        except ImportError:
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    test_fns = [(k, v) for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for name, fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL [%s]: %s" % (name, e))
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)
    print("=" * 60)
    print("OK: tst_directive_extractor_grounding  (%d/%d)"
          % (passed, passed + failed))
    print("=" * 60)


if __name__ == "__main__":
    _standalone_runner()
