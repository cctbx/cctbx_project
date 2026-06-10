"""
v119.H17.1: K-tests for the executor-side strip_flags application.

The executor in programs/ai_agent.py applies a robust regex to remove
flag=value pairs from retry commands.  This test verifies the regex
handles all the legitimate forms a flag value can take:
    flag=value          (basic)
    flag = value        (PHIL spacing — spaces around =)
    flag="path spaces"  (double-quoted with spaces)
    flag='path'         (single-quoted)
    flag at end of line (no trailing space)

A naive `\\S+` pattern would corrupt quoted-with-spaces and miss
PHIL-spacing (Gemini's critique on the H17 plan).

These tests pin the EXACT pattern shipped in ai_agent.py.  If the
pattern is changed there, this test must be updated in lockstep.

Run with: python tests/tst_h17_strip_executor.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import re

# Add parent to path so the test can be run standalone
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =============================================================================
# The strip pattern, copied verbatim from programs/ai_agent.py
# =============================================================================
# This is the contract: the executor uses this exact pattern.
# If you change it here, change it there, and vice versa.
def _strip_flag(command, flag_prefix):
    """Apply the executor's strip logic to a command string."""
    pattern = (r'(?:^|\s)' + re.escape(flag_prefix)
               + r'\s*=\s*(?:"[^"]*"|\'[^\']*\'|\S+)')
    new = re.sub(pattern, '', command).strip()
    new = re.sub(r'\s+', ' ', new)
    return new


# =============================================================================
# Tests
# =============================================================================

def test_h17_strip_basic_flag_equals_value():
    """flag=value (the common case)."""
    cmd = "phenix.autobuild data=foo.mtz map_file=raw.mtz nproc=4"
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result, "map_file= should be removed"
    assert "data=foo.mtz" in result, "other flags must survive"
    assert "nproc=4" in result, "trailing flags must survive"
    print("  PASS: test_h17_strip_basic_flag_equals_value")


def test_h17_strip_phil_spacing_around_equals():
    """flag = value (PHIL-style spaces around =)."""
    cmd = "phenix.autobuild data=foo.mtz map_file = raw.mtz nproc=4"
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result, (
        "PHIL-spaced map_file must be removed; got: %r" % result)
    assert "raw.mtz" not in result, (
        "the value associated with PHIL-spaced map_file must also "
        "go; got: %r" % result)
    assert "data=foo.mtz" in result
    assert "nproc=4" in result
    print("  PASS: test_h17_strip_phil_spacing_around_equals")


def test_h17_strip_double_quoted_value_with_spaces():
    """flag="path with spaces.mtz" (the Gemini-flagged failure mode)."""
    cmd = 'phenix.autobuild data=foo.mtz map_file="path with spaces.mtz" nproc=4'
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result, (
        "quoted map_file should be removed atomically; got: %r" % result)
    assert "spaces.mtz" not in result, (
        "the dangling fragment a naive regex would leave must NOT be "
        "present; got: %r" % result)
    assert "data=foo.mtz" in result
    assert "nproc=4" in result
    print("  PASS: test_h17_strip_double_quoted_value_with_spaces")


def test_h17_strip_single_quoted_value():
    """flag='path.mtz' (single-quoted form)."""
    cmd = "phenix.autobuild data=foo.mtz map_file='/some/path.mtz' nproc=4"
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result
    assert "path.mtz" not in result
    assert "data=foo.mtz" in result
    print("  PASS: test_h17_strip_single_quoted_value")


def test_h17_strip_flag_at_end_of_command():
    """flag=value at the very end (no trailing space)."""
    cmd = "phenix.autobuild data=foo.mtz map_file=raw.mtz"
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result, (
        "end-of-command flag must be removed; got: %r" % result)
    assert result.endswith("data=foo.mtz") or result.endswith("foo.mtz"), (
        "should not leave trailing whitespace; got: %r" % result)
    print("  PASS: test_h17_strip_flag_at_end_of_command")


def test_h17_strip_lysozyme_exact_command():
    """The exact failing command from Tom's lysozyme-MRSAD cycle 5."""
    cmd = ("phenix.autobuild "
           "data=/Users/terwill/Downloads/test4/AutoSol_run_1_/"
           "overall_best_refine_data.mtz "
           "seq_file=/Users/terwill/unix/PHENIX/build_39/modules/"
           "phenix_examples/lysozyme-MRSAD/hewl.seq "
           "model=/Users/terwill/unix/PHENIX/build_39/modules/"
           "phenix_examples/lysozyme-MRSAD/1fkq_prot.pdb "
           "map_file=/Users/terwill/unix/PHENIX/build_39/modules/"
           "phenix_examples/lysozyme-MRSAD/lyso2001_scala1.mtz "
           "resolution=2.0 rebuild_in_place=False nproc=4")
    result = _strip_flag(cmd, "map_file")
    assert "map_file" not in result, (
        "lysozyme map_file must be stripped; got: %r" % result)
    assert "lyso2001_scala1.mtz" not in result, (
        "the bad MTZ path must be gone; got: %r" % result)
    # The remaining four params must survive
    assert "data=" in result
    assert "seq_file=" in result
    assert "model=" in result
    assert "resolution=2.0" in result
    assert "rebuild_in_place=False" in result
    assert "nproc=4" in result
    print("  PASS: test_h17_strip_lysozyme_exact_command")


def test_h17_strip_multiple_variant_prefixes():
    """If the LLM ever uses 'input_map_file=' or
    'input_files.map_file=' instead of bare 'map_file=', the YAML
    strip_parameters covers all three.  Verify each variant strips."""
    for variant in ["map_file", "input_map_file", "input_files.map_file"]:
        cmd = "phenix.autobuild data=foo.mtz %s=raw.mtz nproc=4" % variant
        result = _strip_flag(cmd, variant)
        assert variant not in result, (
            "%s variant should strip; got: %r" % (variant, result))
        assert "data=foo.mtz" in result
        assert "nproc=4" in result
    print("  PASS: test_h17_strip_multiple_variant_prefixes")


def test_h17_strip_does_not_affect_unrelated_flags():
    """Stripping 'map_file' must not touch other flags whose names
    contain 'map' (e.g. 'map_coeffs_file', 'data_map')."""
    cmd = ("phenix.foo data=x.mtz map_coeffs_file=cc.mtz map_file=raw.mtz "
           "data_map=dm.ccp4 nproc=4")
    result = _strip_flag(cmd, "map_file")
    assert "map_file=raw.mtz" not in result, "target should be gone"
    # These superficially-similar flags must survive:
    assert "map_coeffs_file=cc.mtz" in result, (
        "different flag must survive; got: %r" % result)
    assert "data_map=dm.ccp4" in result, (
        "different flag must survive; got: %r" % result)
    print("  PASS: test_h17_strip_does_not_affect_unrelated_flags")


def test_h17_strip_idempotent_when_not_present():
    """Stripping a flag that isn't in the command is a no-op."""
    cmd = "phenix.autobuild data=foo.mtz nproc=4"
    result = _strip_flag(cmd, "map_file")
    # Whitespace is normalized; check semantic equivalence
    assert result == cmd or result == cmd.strip(), (
        "no-op strip should not modify the command; got: %r" % result)
    print("  PASS: test_h17_strip_idempotent_when_not_present")


# =============================================================================
# RUNNER
# =============================================================================

def run_all_tests():
    tests = [
        test_h17_strip_basic_flag_equals_value,
        test_h17_strip_phil_spacing_around_equals,
        test_h17_strip_double_quoted_value_with_spaces,
        test_h17_strip_single_quoted_value,
        test_h17_strip_flag_at_end_of_command,
        test_h17_strip_lysozyme_exact_command,
        test_h17_strip_multiple_variant_prefixes,
        test_h17_strip_does_not_affect_unrelated_flags,
        test_h17_strip_idempotent_when_not_present,
    ]
    passed = 0
    failed = []
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed.append((t.__name__, str(e)))
            print("  FAIL: %s — %s" % (t.__name__, e))
    print()
    print("Total: %d passed, %d failed" % (passed, len(failed)))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
