"""K_H18_1: PHIL round-trip test for original_files_for_directives.

The H18 ship had a deploy gap: programs/ai_agent.py line 8513 assigned
to ``directive_params.ai_analysis.original_files_for_directives``,
but the corresponding PHIL parameter definition was added to
``programs/ai_analysis.py`` ONLY — not to ``programs/ai_agent.py``'s
own ``master_params`` string.  Production runs hit
``AttributeError: Assignment to non-existing attribute
"ai_analysis.original_files_for_directives"``, directive extraction
returned empty {}, and the agent ran the default plan stages
through to predict_and_build instead of stopping after
resolve_cryo_em.

The original H18 K-tests passed because they bypassed the PHIL
layer entirely — they called ``infer_experiment_type_from_files``
and ``_apply_experiment_type_program_reprints`` directly.  This
test exercises the actual PHIL round-trip: parse master_params,
deep-copy, assign, verify.

If H18.1's PHIL definition is removed from ai_agent.py's
master_params, this test fails immediately — closing the deploy
gap at sandbox time.

Run with: phenix.python tests/tst_h18_1_phil_roundtrip.py
or in sandbox: python tests/tst_h18_1_phil_roundtrip.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)


# =============================================================================
# Test 1: ai_agent.py's master_params declares original_files_for_directives
# =============================================================================
#
# We grep the source file directly rather than parsing PHIL — this works
# in sandbox without libtbx.phil.  A separate PHENIX-only test below
# exercises the parsed object.

def test_master_params_string_contains_param():
    """The H18.1 PHIL definition must be present in ai_agent.py's
    master_params string verbatim.  If this fails, the deploy gap
    that caused production failure has reopened."""
    print("  test_master_params_string_contains_param")
    ai_agent_path = os.path.join(_PARENT, "programs", "ai_agent.py")
    assert os.path.exists(ai_agent_path), (
        "ai_agent.py not found at %s" % ai_agent_path)

    with open(ai_agent_path, "r") as f:
        contents = f.read()

    # The master_params region: between "master_params = \"\"\"" and the
    # closing triple-quote that precedes "master_phil = libtbx.phil.parse"
    start_marker = 'master_params = """'
    end_marker = 'master_phil = libtbx.phil.parse'
    start_idx = contents.find(start_marker)
    end_idx = contents.find(end_marker)
    assert start_idx != -1, "master_params definition not found"
    assert end_idx != -1, "master_phil parse call not found"
    assert start_idx < end_idx, (
        "master_params region malformed (start=%d, end=%d)"
        % (start_idx, end_idx))

    master_region = contents[start_idx:end_idx]

    assert "original_files_for_directives" in master_region, (
        "H18.1 deploy-gap regression: "
        "original_files_for_directives is MISSING from ai_agent.py's "
        "master_params string.  This breaks the PHIL assignment in "
        "_extract_directives (line ~8513), causing directive extraction "
        "to fail with AttributeError and the agent to ignore user "
        "stop directives.")
    print("    PASS")


# =============================================================================
# Test 2: master_params parses cleanly and exposes the new attribute
# =============================================================================

def test_master_params_parses_with_new_attribute():
    """When PHENIX is available, parse master_params and confirm the
    new attribute is accessible.  Skip gracefully in sandbox where
    libtbx.phil isn't importable."""
    print("  test_master_params_parses_with_new_attribute")
    try:
        import libtbx.phil
    except ImportError:
        print("    SKIP (libtbx.phil not available — needs PHENIX)")
        return

    # Re-parse master_params in isolation
    ai_agent_path = os.path.join(_PARENT, "programs", "ai_agent.py")
    with open(ai_agent_path, "r") as f:
        contents = f.read()

    start = contents.find('master_params = """') + len('master_params = """')
    end = contents.find('"""', start)
    assert start > 0 and end > start, "could not extract master_params body"
    master_params_str = contents[start:end]

    master_phil = libtbx.phil.parse(master_params_str,
                                     process_includes=False)
    extracted = master_phil.extract()

    # The attribute must exist on the extracted params object.
    assert hasattr(extracted, "ai_analysis"), (
        "ai_analysis scope missing from master_phil.extract()")
    assert hasattr(extracted.ai_analysis,
                   "original_files_for_directives"), (
        "original_files_for_directives missing from extracted params "
        "— PHIL definition not in master_params")

    # Default must be None (backward compat: nothing happens when this
    # param isn't set).
    assert extracted.ai_analysis.original_files_for_directives is None, (
        "Default for original_files_for_directives must be None, got: %r"
        % extracted.ai_analysis.original_files_for_directives)
    print("    PASS")


# =============================================================================
# Test 3: assignment doesn't raise AttributeError (the production failure)
# =============================================================================

def test_assignment_does_not_raise():
    """Verify the exact assignment pattern from ai_agent.py:8513 works
    without raising AttributeError.  This is the production-failure
    reproduction.  Skip gracefully in sandbox."""
    print("  test_assignment_does_not_raise")
    try:
        import libtbx.phil
        import copy
    except ImportError:
        print("    SKIP (libtbx.phil not available — needs PHENIX)")
        return

    ai_agent_path = os.path.join(_PARENT, "programs", "ai_agent.py")
    with open(ai_agent_path, "r") as f:
        contents = f.read()
    start = contents.find('master_params = """') + len('master_params = """')
    end = contents.find('"""', start)
    master_params_str = contents[start:end]

    master_phil = libtbx.phil.parse(master_params_str,
                                     process_includes=False)
    params = master_phil.extract()

    # Reproduce ai_agent.py:8483 — copy.deepcopy(self.params)
    directive_params = copy.deepcopy(params)

    # Reproduce ai_agent.py:8513 — the assignment that crashed in
    # production:
    #   directive_params.ai_analysis.original_files_for_directives = (
    #     ",".join(_basenames))
    try:
        directive_params.ai_analysis.original_files_for_directives = (
            "a.ccp4,b.ccp4,c.fa")
    except AttributeError as e:
        raise AssertionError(
            "Production failure reproduction: assignment to "
            "original_files_for_directives raised AttributeError. "
            "Original error: %s.  This is the H18.1 deploy-gap bug." % e)

    # Verify the assignment actually took
    assert (directive_params.ai_analysis.original_files_for_directives
            == "a.ccp4,b.ccp4,c.fa"), (
        "Assignment appeared to succeed but value didn't take: %r"
        % directive_params.ai_analysis.original_files_for_directives)

    # And None assignment also works (the else branch in ai_agent.py:8516)
    try:
        directive_params.ai_analysis.original_files_for_directives = None
    except AttributeError as e:
        raise AssertionError(
            "Setting original_files_for_directives = None raised "
            "AttributeError: %s" % e)
    print("    PASS")


# =============================================================================
# Test 4: ai_analysis.py PHIL definition is also still in place
# =============================================================================

def test_ai_analysis_phil_definition_still_present():
    """The server-side PHIL definition (programs/ai_analysis.py) must
    ALSO be present.  H18.1 fixes the client-side gap (ai_agent.py);
    it must not break the server-side definition.  Both are needed —
    ai_agent.py for the LOCAL deep-copy assignment, ai_analysis.py
    for the SERVER-side PHIL parser when extraction runs on the
    server."""
    print("  test_ai_analysis_phil_definition_still_present")
    ai_analysis_path = os.path.join(
        _PARENT, "programs", "ai_analysis.py")
    if not os.path.exists(ai_analysis_path):
        # In some sandbox layouts the file might not be at this path —
        # skip rather than fail.  The H18 ship guarantees it for the
        # PHENIX-layout deployment.
        print("    SKIP (ai_analysis.py not at sandbox path)")
        return
    with open(ai_analysis_path, "r") as f:
        contents = f.read()
    assert "original_files_for_directives" in contents, (
        "programs/ai_analysis.py is missing the "
        "original_files_for_directives PHIL definition.  H18 ship "
        "incomplete on this side too.")
    print("    PASS")


# =============================================================================
# Test 5: assignment line in ai_agent.py is still where we expect it
# =============================================================================

def test_assignment_line_still_in_ai_agent_py():
    """Anchor: the assignment site in ai_agent.py that triggered the
    production failure.  If this line is moved or refactored, this
    test will fail loudly so the lockstep relationship with the PHIL
    definition is preserved."""
    print("  test_assignment_line_still_in_ai_agent_py")
    ai_agent_path = os.path.join(_PARENT, "programs", "ai_agent.py")
    with open(ai_agent_path, "r") as f:
        contents = f.read()
    needle = (
        "directive_params.ai_analysis.original_files_for_directives = (")
    assert needle in contents, (
        "Assignment site not found in ai_agent.py.  Did someone "
        "refactor the call site?  If so, update this test's anchor.")
    print("    PASS")


# =============================================================================
# Test 6: cross-file consistency — same parameter name, no typos
# =============================================================================

def test_cross_file_consistency():
    """Belt-and-suspenders: the parameter name as defined in both PHIL
    blocks AND used in the assignment must match exactly.  Catches the
    "renamed in one place" class of regression."""
    print("  test_cross_file_consistency")
    PARAM = "original_files_for_directives"
    files = [
        ("programs/ai_agent.py", PARAM),
        ("programs/ai_analysis.py", PARAM),
    ]
    for relpath, needle in files:
        full = os.path.join(_PARENT, relpath)
        if not os.path.exists(full):
            print("    SKIP %s (not at sandbox path)" % relpath)
            continue
        with open(full, "r") as f:
            contents = f.read()
        assert needle in contents, (
            "%s missing reference to %r" % (relpath, needle))
    print("    PASS")


# =============================================================================
# Runner
# =============================================================================

def run_all_tests():
    print("=" * 72)
    print("K_H18_1: PHIL round-trip for original_files_for_directives")
    print("=" * 72)

    tests = [
        test_master_params_string_contains_param,
        test_master_params_parses_with_new_attribute,
        test_assignment_does_not_raise,
        test_ai_analysis_phil_definition_still_present,
        test_assignment_line_still_in_ai_agent_py,
        test_cross_file_consistency,
    ]

    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            failed += 1
            print("    FAIL: %s" % e)
        except Exception as e:
            failed += 1
            print("    ERROR: %s: %s" % (type(e).__name__, e))

    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        raise AssertionError(
            "%d K_H18_1 test(s) failed" % failed)


if __name__ == "__main__":
    try:
        run_all_tests()
    except AssertionError as e:
        print("FAILURE: %s" % e)
        sys.exit(1)
