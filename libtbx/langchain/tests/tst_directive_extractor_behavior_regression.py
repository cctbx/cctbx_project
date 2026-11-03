"""v116.15 behavior regression — runs against ANY version.

Demonstrates the AF_7mjs bug by calling
`_apply_crystal_symmetry_fallback` directly with the preprocessed
advice from run 79.  Pre-v116.15 the bug causes
test_af7mjs_no_phantom_space_group to FAIL (extracts
space_group="None").  v116.15 makes all tests pass.

Kept separate from the main test file because it does NOT import the
new `_is_symmetry_sentinel` helper, so it can run against the
unfixed code to demonstrate the bug.
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


from agent.directive_extractor import _apply_crystal_symmetry_fallback


def _quiet(msg):
    pass


# Exact preprocessed advice excerpt from AF_7mjs run 79.
AF_7MJS_ADVICE = """1. **Input Files Found**: 7mjs_23883_H.fa, 7mjs_23883_H_1.ccp4

4. **Key Parameters**:
   - Wavelength: None
   - Resolution limit: Not specified (auto from half-maps)
   - Number of expected sites: None
   - Heavy atom type: None
   - Space group: None
"""


def test_af7mjs_no_phantom_space_group():
    print("[test] AF_7mjs preprocessed advice yields no space_group directive")
    result = _apply_crystal_symmetry_fallback({}, AF_7MJS_ADVICE, _quiet)
    got = (result.get("program_settings", {})
                  .get("default", {})
                  .get("space_group"))
    assert got is None, (
        "BUG (pre-v116.15): space_group should not be extracted from "
        "'Space group: None' placeholder, got: %r\nFull result: %r"
        % (got, result))
    print("  PASS")


def test_real_space_group_still_extracted():
    print("[test] Real space group is still extracted normally")
    result = _apply_crystal_symmetry_fallback(
        {}, "space group: P 21 21 21", _quiet)
    got = (result.get("program_settings", {})
                  .get("default", {})
                  .get("space_group"))
    assert got is not None and got.startswith("P"), \
        "Real space group should be extracted, got: %r" % got
    print("  PASS — extracted as %r" % got)


def test_existing_directive_preserved():
    print("[test] Existing space_group is not overwritten by sentinel")
    existing = {"program_settings": {"default": {"space_group": "P 1"}}}
    result = _apply_crystal_symmetry_fallback(
        existing, "Space group: None", _quiet)
    got = (result.get("program_settings", {})
                  .get("default", {})
                  .get("space_group"))
    assert got == "P 1", \
        "Existing space_group should be preserved, got: %r" % got
    print("  PASS — existing 'P 1' preserved")


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
    print("OK: tst_directive_extractor_behavior_regression  (%d/%d)"
          % (passed, passed + failed))
    print("=" * 60)


if __name__ == "__main__":
    _standalone_runner()
