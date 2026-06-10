"""Tests for v116.15 sentinel-value rejection in directive extractor.

Bug: `_apply_crystal_symmetry_fallback` in `agent/directive_extractor.py`
ran a regex against the user advice to extract space_group as a defensive
fallback in case the directive-extraction LLM missed it.  The regex
matched any "space group <value>" pattern, including the standardized
preprocessor placeholders that the advice-preprocessor LLM emits when
the user did NOT specify a value.

A typical AF_7mjs / cryo-EM README produces preprocessed advice
containing:

    4. **Key Parameters**:
       - Wavelength: None
       - Resolution limit: Not specified
       - Number of expected sites: None
       - Heavy atom type: None
       - Space group: None

The regex above happily extracted `space_group="None"` from the last
line and stored it as a real directive.  `validate_directives()` already
rejects sentinel values but only on the LLM-extracted path — the regex
fallback ran AFTER validation and reintroduced the sentinel.

Fix: factor the sentinel check into a shared helper
`_is_symmetry_sentinel(value)` and call it from BOTH paths.
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
    ic.classify_intent = lambda advice: None


from agent.directive_extractor import (
    _apply_crystal_symmetry_fallback,
    _is_symmetry_sentinel,
    validate_directives,
)


def _quiet(msg):
    pass


def _sg_in(result):
    """Extract space_group from a directives dict, returning None if absent."""
    return (result.get("program_settings", {})
                  .get("default", {})
                  .get("space_group"))


# =============================================================================
# AF_7mjs regression: "Space group: None" must NOT be extracted
# =============================================================================

def test_af7mjs_space_group_none_rejected():
    print("[test] AF_7mjs regression: 'Space group: None' regex match")
    advice = """1. **Input Files Found**: 7mjs_23883_H.fa

4. **Key Parameters**:
   - Wavelength: None
   - Resolution limit: Not specified (auto from half-maps)
   - Space group: None
"""
    result = _apply_crystal_symmetry_fallback({}, advice, _quiet)
    assert _sg_in(result) is None, \
        "AF_7mjs regression: space_group=None should be REJECTED, got: %r" % result
    print("  PASS")


# =============================================================================
# Sentinel variants on the regex-fallback path
# =============================================================================

def test_sentinel_not_specified():
    print("[test] Sentinel: 'Not specified'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: Not specified", _quiet)
    assert _sg_in(result) is None, "'Not specified' should be rejected"
    print("  PASS")


def test_sentinel_na():
    print("[test] Sentinel: 'N/A'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: N/A", _quiet)
    assert _sg_in(result) is None, "'N/A' should be rejected"
    print("  PASS")


def test_sentinel_unknown_uppercase():
    print("[test] Sentinel: 'UNKNOWN' (case-insensitive)")
    result = _apply_crystal_symmetry_fallback({}, "Space group: UNKNOWN", _quiet)
    assert _sg_in(result) is None, "'UNKNOWN' should be rejected (case-insensitive)"
    print("  PASS")


def test_sentinel_tbd():
    print("[test] Sentinel: 'TBD'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: TBD", _quiet)
    assert _sg_in(result) is None, "'TBD' should be rejected"
    print("  PASS")


def test_sentinel_lowercase_none():
    print("[test] Sentinel: lowercase 'none'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: none", _quiet)
    assert _sg_in(result) is None, "lowercase 'none' should be rejected"
    print("  PASS")


# =============================================================================
# Real space groups must STILL be extracted (no regression on happy path)
# =============================================================================

def test_real_sg_p32_2_1():
    print("[test] Real space group: 'P 32 2 1'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: P 32 2 1", _quiet)
    got = _sg_in(result)
    assert got and got.startswith("P"), \
        "Real space group 'P 32 2 1' should be EXTRACTED, got: %r" % got
    print("  PASS — extracted as %r" % got)


def test_real_sg_p63():
    print("[test] Real space group: 'P63'")
    result = _apply_crystal_symmetry_fallback({}, "Space group: P63", _quiet)
    got = _sg_in(result)
    assert got == "P63", \
        "Real space group 'P63' should be EXTRACTED as P63, got: %r" % got
    print("  PASS")


def test_real_sg_c2_phrasing():
    print("[test] Real space group: 'space group is C 2 2 21'")
    result = _apply_crystal_symmetry_fallback({}, "space group is C 2 2 21", _quiet)
    got = _sg_in(result)
    # The existing regex captures only the first letter+digits chunk for
    # space-separated symbols.  This test documents that behavior is
    # unchanged by the sentinel fix — we just don't reject it.
    assert got and got[0] == "C", \
        "Real space group starting with 'C' should be EXTRACTED, got: %r" % got
    print("  PASS — extracted as %r (existing regex behavior)" % got)


# =============================================================================
# Helper function behavior
# =============================================================================

def test_helper_none_is_sentinel():
    print("[test] _is_symmetry_sentinel: 'None' → True")
    assert _is_symmetry_sentinel("None") is True
    print("  PASS")


def test_helper_strips_whitespace():
    print("[test] _is_symmetry_sentinel: '  None  ' → True (whitespace)")
    assert _is_symmetry_sentinel("  None  ") is True
    print("  PASS")


def test_helper_real_sg_not_sentinel():
    print("[test] _is_symmetry_sentinel: 'P 32 2 1' → False")
    assert _is_symmetry_sentinel("P 32 2 1") is False
    print("  PASS")


def test_helper_non_string_input():
    print("[test] _is_symmetry_sentinel: non-string input → False")
    assert _is_symmetry_sentinel(None) is False
    assert _is_symmetry_sentinel(123) is False
    assert _is_symmetry_sentinel([]) is False
    print("  PASS — None / int / list all return False")


# =============================================================================
# validate_directives also rejects sentinels (consistent paths)
# =============================================================================

def test_validate_directives_drops_sentinel():
    print("[test] validate_directives drops space_group='None' (LLM path)")
    bad = {"program_settings": {"default": {"space_group": "None"}}}
    out = validate_directives(bad, _quiet)
    assert not (out.get("program_settings", {})
                   .get("default", {})
                   .get("space_group")), \
        "validate_directives should drop space_group='None', got: %r" % out
    print("  PASS")


def test_validate_directives_keeps_real_sg():
    print("[test] validate_directives keeps space_group='P 32 2 1' (LLM path)")
    good = {"program_settings": {"default": {"space_group": "P 32 2 1"}}}
    out = validate_directives(good, _quiet)
    got = (out.get("program_settings", {})
              .get("default", {})
              .get("space_group"))
    assert got == "P 32 2 1", "Expected 'P 32 2 1', got: %r" % got
    print("  PASS")


# =============================================================================
# Test runner (matches the project's tst_*.py convention)
# =============================================================================

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
    print("=" * 60)
    print("OK: tst_directive_extractor_sentinel_space_group  (%d/%d)"
          % (passed, passed + failed))
    print("=" * 60)


if __name__ == "__main__":
    _standalone_runner()
