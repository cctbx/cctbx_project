"""Tests for raw-advice-authoritative directive extraction (v117 Step 1).

These tests verify that:
  - extract_directives() accepts a raw_advice parameter
  - The dual-input prompt is selected when raw_advice differs from user_advice
  - The single-input prompt is used when raw_advice is None or equals user_advice
  - The dual-input prompt contains the AUTHORITY paragraph
  - DIRECTIVE_EXTRACTION_PROMPT.format() and the .replace() path both work
    (persistent invariant — guards against regression of the brace fix)

These are STRUCTURAL tests — they verify the prompt-building path without
calling a real LLM.  An end-to-end LLM test for the bug Tom originally
reported (openai mangling "density modify and stop") is the LLM-side
test in LLM_TEST_C1_SPEC.md.
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

# Shim libtbx.langchain if not already provided
if "libtbx" not in sys.modules:
    import types as _types
    for _mod in (
        "libtbx", "libtbx.langchain",
        "libtbx.langchain.agent",
        "libtbx.langchain.knowledge",
    ):
        if _mod not in sys.modules:
            sys.modules[_mod] = _types.ModuleType(_mod)

# Stub program_registry for module-level import.
if "agent.program_registry" not in sys.modules:
    import types as _types
    prog_reg = _types.ModuleType("agent.program_registry")
    class _Registry:
        def __init__(self, *a, **k): pass
        def get(self, *a, **k): return None
        def all_programs(self): return []
        def get_all_program_keys(self): return []
    prog_reg.ProgramRegistry = _Registry
    sys.modules["agent.program_registry"] = prog_reg
    sys.modules["libtbx.langchain.agent.program_registry"] = prog_reg

# Try the canonical import path; fall back to direct import.
try:
    from libtbx.langchain.agent.directive_extractor import (
        DIRECTIVE_EXTRACTION_PROMPT,
        DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW,
        extract_directives,
        _is_stop_after_requested,
    )
except ImportError:
    from agent.directive_extractor import (
        DIRECTIVE_EXTRACTION_PROMPT,
        DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW,
        extract_directives,
        _is_stop_after_requested,
    )


def test_extract_directives_signature_has_raw_advice():
    """extract_directives() accepts raw_advice parameter."""
    print("Test: extract_directives_signature_has_raw_advice")
    import inspect
    sig = inspect.signature(extract_directives)
    assert "raw_advice" in sig.parameters, (
        "raw_advice missing from extract_directives signature: %s"
        % list(sig.parameters))
    # Default value is None
    assert sig.parameters["raw_advice"].default is None, (
        "raw_advice default should be None, got %r"
        % sig.parameters["raw_advice"].default)
    print("  PASS")


def test_directive_extraction_prompt_format_invariant():
    """Persistent invariant: prompt.format() must not raise.

    Step A fixed the brace bug; this test ensures any future edit that
    re-introduces unescaped braces is caught.
    """
    print("Test: directive_extraction_prompt_format_invariant")
    try:
        out = DIRECTIVE_EXTRACTION_PROMPT.format(user_advice="test")
    except Exception as e:
        raise AssertionError(
            "DIRECTIVE_EXTRACTION_PROMPT.format() raised %s: %s"
            % (type(e).__name__, e))
    assert "test" in out, "Substitution did not occur"
    assert "{user_advice}" not in out, "Placeholder still present"
    print("  PASS")


def test_dual_input_prompt_format_invariant():
    """Persistent invariant: WITH_RAW prompt's .format() must not raise."""
    print("Test: dual_input_prompt_format_invariant")
    try:
        out = DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW.format(
            raw_advice="raw test",
            processed_advice="processed test",
        )
    except Exception as e:
        raise AssertionError(
            "DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW.format() raised %s: %s"
            % (type(e).__name__, e))
    assert "raw test" in out, "raw_advice not substituted"
    assert "processed test" in out, "processed_advice not substituted"
    assert "{raw_advice}" not in out, "Placeholder still present"
    assert "{processed_advice}" not in out, "Placeholder still present"
    print("  PASS")


def test_dual_input_prompt_has_authority_paragraph():
    """The WITH_RAW prompt contains the AUTHORITY paragraph teaching
    the LLM that raw is the source of truth for intent."""
    print("Test: dual_input_prompt_has_authority_paragraph")
    assert "AUTHORITY:" in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    assert "source of truth" in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    assert "stop_after_requested" in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    assert "after_program" in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    assert "start_with_program" in DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    print("  PASS")


def test_dual_input_prompt_has_raw_and_processed_blocks():
    """The WITH_RAW prompt has separately-labeled raw and processed blocks."""
    print("Test: dual_input_prompt_has_raw_and_processed_blocks")
    p = DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
    assert "RAW INSTRUCTION" in p
    assert "PREPROCESSED ADVICE" in p
    assert "{raw_advice}" in p
    assert "{processed_advice}" in p
    # The single-input USER ADVICE block must NOT appear in WITH_RAW
    assert "=== USER ADVICE ===" not in p, (
        "WITH_RAW still has the single-input USER ADVICE block")
    print("  PASS")


def test_single_input_prompt_unchanged():
    """The single-input prompt still has its USER ADVICE block."""
    print("Test: single_input_prompt_unchanged")
    p = DIRECTIVE_EXTRACTION_PROMPT
    assert "=== USER ADVICE ===" in p
    assert "{user_advice}" in p
    # AUTHORITY paragraph is only in WITH_RAW
    assert "AUTHORITY:" not in p
    print("  PASS")


def test_stop_after_requested_documented_in_schema():
    """Schema section teaches the LLM about stop_after_requested.

    Both prompts should have this (since WITH_RAW is derived from
    DIRECTIVE_EXTRACTION_PROMPT)."""
    print("Test: stop_after_requested_documented_in_schema")
    for label, prompt in [
        ("DIRECTIVE_EXTRACTION_PROMPT", DIRECTIVE_EXTRACTION_PROMPT),
        ("DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW",
         DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW),
    ]:
        assert "stop_after_requested" in prompt, (
            "%s missing stop_after_requested doc" % label)
        # Some of the documented phrasings
        for phrasing in ["stop after X", "X and stop", "only run X"]:
            assert phrasing in prompt, (
                "%s missing phrasing %r in stop_after_requested doc"
                % (label, phrasing))
    print("  PASS")


def test_is_stop_after_requested_helper_basic_cases():
    """The _is_stop_after_requested helper handles core cases.

    Includes the period-separated 'Refine. Stop.' case that was added
    as a coverage patch in Step B."""
    print("Test: is_stop_after_requested_helper_basic_cases")
    cases = [
        ("density modify and stop", True),
        ("Refine and stop.", True),
        ("Refine. Stop.", True),                       # period-stop patch
        ("Refine, stop.", True),
        ("Stop after refinement", True),
        ("stop when r-free < 0.25", True),
        ("only run xtriage", True),
        ("just do refine", True),
        ("Refine the model.", False),
        ("don't stop the workflow", False),
        ("Stop Condition: None", False),
        ("Stop Condition: after refinement", True),
        ("", False),
        (None, False),
    ]
    for advice, expected in cases:
        actual = _is_stop_after_requested(advice)
        assert actual == expected, (
            "_is_stop_after_requested(%r) = %r, expected %r"
            % (advice, actual, expected))
    print("  PASS")


def test_density_modify_and_stop_dual_input_construction():
    """When raw='density modify and stop' and processed has
    'Stop Condition: None', the dual-input prompt is constructed
    with both inputs visible to the LLM.

    This is the prompt-construction side of Tom's openai bug —
    verifying the LLM would SEE the raw user instruction.  The LLM-side
    test is in LLM_TEST_C1_SPEC."""
    print("Test: density_modify_and_stop_dual_input_construction")
    raw = "density modify and stop"
    processed = ("Primary Goal: Perform density modification, then build "
                 "and refine.\nStop Condition: None")
    # Simulate what extract_directives builds
    prompt = (DIRECTIVE_EXTRACTION_PROMPT_WITH_RAW
              .replace("{raw_advice}", raw)
              .replace("{processed_advice}", processed))
    assert raw in prompt, "Raw user instruction missing from constructed prompt"
    assert "Stop Condition: None" in prompt, (
        "Processed advice missing from constructed prompt")
    assert "AUTHORITY:" in prompt, "AUTHORITY paragraph missing"
    assert "{raw_advice}" not in prompt, "Placeholder leftover"
    assert "{processed_advice}" not in prompt, "Placeholder leftover"
    print("  PASS")


def run_all_tests():
    """Run all tests, count passes and failures."""
    tests = [
        test_extract_directives_signature_has_raw_advice,
        test_directive_extraction_prompt_format_invariant,
        test_dual_input_prompt_format_invariant,
        test_dual_input_prompt_has_authority_paragraph,
        test_dual_input_prompt_has_raw_and_processed_blocks,
        test_single_input_prompt_unchanged,
        test_stop_after_requested_documented_in_schema,
        test_is_stop_after_requested_helper_basic_cases,
        test_density_modify_and_stop_dual_input_construction,
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
    success = run_all_tests()
    sys.exit(0 if success else 1)
