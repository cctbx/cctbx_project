"""
Sandbox tests for v118 Section C-prime: directive-layer diagnostic
split.

C-prime adds two diagnostic concepts to _query_agent_for_command:

  - directives_user_intent: the canonical directives extracted from
    user advice (= session.data["directives"]).  Stable across
    cycles.  Captured BEFORE the per-cycle plan merge.

  - directives_effective_runtime: the per-cycle merge of user
    directives with the current plan stage.  This is what the LLM
    receives and what gets logged as "directives" in DIAG_PLAN.

These tests verify:

  K_C1  Source-code regression guard: _query_agent_for_command
        reads session.get_directives() but never writes back to
        session.data["directives"].  This guards against a future
        developer accidentally turning the local merge into a
        persistent side-effect.

  K_C2  Source-code regression guard: the user_intent snapshot is
        taken with deepcopy BEFORE the merge, so subsequent
        mutations of the runtime dict cannot retroactively change
        the snapshot.

  K_C3  Source-code regression guard: session_info includes BOTH
        "directives" (runtime) AND "directives_user_intent" so the
        server can log both layers.

  K_C4  Source-code regression guard: the diagnostic emits both
        labels distinctly (DIRECTIVE_LAYERS user_intent / DIRECTIVE_LAYERS
        effective_runtime).

  K_C5  Source-code regression guard: the diagnostic only fires
        when the two layers differ (no log noise when they're
        identical).

  K_C6  Source-code regression guard: the diagnostic is in a
        try/except so it cannot break the cycle.

These are static-source assertions because the runtime path runs
inside the agent loop and isn't unit-testable without launching
phenix.  The asserts read the source of programs/ai_agent.py and
check structural invariants.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

AI_AGENT_PATH = os.path.join(ROOT, "programs", "ai_agent.py")


def _read_ai_agent():
    with open(AI_AGENT_PATH) as f:
        return f.read()


def _extract_query_agent_for_command(src):
    """Return the source of _query_agent_for_command method."""
    import re
    m = re.search(
        r'  def _query_agent_for_command\(self.*?\n(?=  def [^\s]|\Z)',
        src,
        re.DOTALL,
    )
    if not m:
        raise AssertionError(
            "Could not locate _query_agent_for_command in ai_agent.py")
    return m.group(0)


def test_k_c1_no_persistent_writeback():
    """K_C1: _query_agent_for_command must not write session.data['directives']."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    # The function reads via session.get_directives() — that's fine
    assert "session.get_directives()" in fn_src, (
        "K_C1 failed: _query_agent_for_command should call "
        "session.get_directives()")

    # But it must NEVER do session.data["directives"] = ... — that
    # would be a persistent side-effect.
    import re
    bad_patterns = [
        r'session\.data\[.directives.\]\s*=',
        r'session\.set_directives\s*\(',
    ]
    for pat in bad_patterns:
        m = re.search(pat, fn_src)
        assert m is None, (
            "K_C1 failed: _query_agent_for_command contains "
            "persistent directive write pattern: %s" % pat)
    print("  PASS: K_C1 (no persistent writeback)")


def test_k_c2_deepcopy_before_merge():
    """K_C2: user_intent snapshot taken with deepcopy BEFORE merge."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    # Find the deepcopy line
    assert "_directives_user_intent" in fn_src, (
        "K_C2 failed: _directives_user_intent variable not defined")
    assert "deepcopy" in fn_src, (
        "K_C2 failed: should use deepcopy for the user_intent snapshot")

    # The deepcopy must appear BEFORE the merge_directives call
    deepcopy_pos = fn_src.find("_copy_module.deepcopy")
    merge_pos = fn_src.find("merge_directives(")
    assert deepcopy_pos > 0, "K_C2 failed: deepcopy call not found"
    assert merge_pos > 0, "K_C2 failed: merge_directives call not found"
    assert deepcopy_pos < merge_pos, (
        "K_C2 failed: deepcopy (%d) must precede merge_directives (%d) "
        "so the snapshot captures pre-merge state"
        % (deepcopy_pos, merge_pos))
    print("  PASS: K_C2 (deepcopy taken before merge)")


def test_k_c3_session_info_both_layers():
    """K_C3: session_info dict includes both 'directives' and 'directives_user_intent'."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    # Both keys must appear in the session_info construction
    assert '"directives": directives,' in fn_src, (
        "K_C3 failed: session_info must include 'directives' key "
        "with the runtime value")
    assert '"directives_user_intent": _directives_user_intent' in fn_src, (
        "K_C3 failed: session_info must include "
        "'directives_user_intent' key with the pre-merge snapshot")
    print("  PASS: K_C3 (session_info contains both layers)")


def test_k_c4_diagnostic_emits_both_labels():
    """K_C4: diagnostic log distinguishes user_intent vs effective_runtime."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    assert "[DIRECTIVE_LAYERS] user_intent" in fn_src, (
        "K_C4 failed: diagnostic should label user_intent layer")
    assert "[DIRECTIVE_LAYERS] effective_runtime" in fn_src, (
        "K_C4 failed: diagnostic should label effective_runtime layer")
    print("  PASS: K_C4 (diagnostic emits both labels)")


def test_k_c5_diagnostic_only_when_differ():
    """K_C5: diagnostic only fires when user_intent != effective_runtime."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    # The if-guard before the diagnostic must compare the two
    import re
    # Look for the comparison _directives_user_intent != directives
    assert re.search(
        r'_directives_user_intent\s*!=\s*directives',
        fn_src,
    ), (
        "K_C5 failed: diagnostic should be guarded by "
        "_directives_user_intent != directives")
    print("  PASS: K_C5 (diagnostic conditional on layer difference)")


def test_k_c6_diagnostic_in_try_except():
    """K_C6: diagnostic emission is in try/except so it can't break the cycle."""
    src = _read_ai_agent()
    fn_src = _extract_query_agent_for_command(src)

    # Locate the diagnostic block and confirm it sits in a try/except
    import re
    # Pattern: try: ... [DIRECTIVE_LAYERS] ... except
    block = re.search(
        r'try:\s*\n[^\n]*\n[^\n]*\n[^\n]*\n[^\n]*\[DIRECTIVE_LAYERS\]',
        fn_src,
        re.DOTALL,
    )
    if block is None:
        # Try a broader match
        diag_pos = fn_src.find("[DIRECTIVE_LAYERS]")
        assert diag_pos > 0, "K_C6 failed: diagnostic not found"
        # Look backward for try:
        before = fn_src[:diag_pos]
        last_try = before.rfind("try:")
        # And forward for except:
        after = fn_src[diag_pos:]
        first_except = after.find("except")
        assert last_try > 0, (
            "K_C6 failed: no try: before the diagnostic")
        assert first_except > 0, (
            "K_C6 failed: no except: after the diagnostic")
    print("  PASS: K_C6 (diagnostic guarded by try/except)")


def test_k_c7_no_persistent_writeback_global():
    """K_C7: global scan — only the documented sites write
    session.data['directives'].

    Per the investigation, the only legitimate writers are:
      - _extract_directives in programs/ai_agent.py
        (once-per-session, at extraction time)
      - _initialize_plan_inner in programs/ai_agent.py
        (once-per-session, three preprocessing/needs-plan branches)

    No NEW write sites should appear in _query_agent_for_command or
    anywhere else added by C-prime.
    """
    src = _read_ai_agent()
    import re
    # Find every site that writes session.data["directives"]
    writes = []
    for m in re.finditer(
            r'session\.data\[.directives.\]\s*=',
            src):
        line_no = src[:m.start()].count("\n") + 1
        writes.append(line_no)

    # Map each line to its enclosing function
    fn_starts = []
    for m in re.finditer(r'  def (_\w+)\(self', src):
        ln = src[:m.start()].count("\n") + 1
        fn_starts.append((ln, m.group(1)))

    def fn_for_line(line):
        owner = None
        for start, name in fn_starts:
            if start <= line:
                owner = name
            else:
                break
        return owner

    allowed_owners = {
        "_extract_directives",
        "_initialize_plan_inner",
    }
    for ln in writes:
        owner = fn_for_line(ln)
        assert owner in allowed_owners, (
            "K_C7 failed: session.data['directives'] is written at "
            "line %d, inside %s(), which is NOT in the allowed set "
            "%s.  If this is a legitimate new write site, update "
            "K_C7's allowed set and document the rationale."
            % (ln, owner, sorted(allowed_owners)))
    print("  PASS: K_C7 (writes confined to known startup sites: %s)"
          % sorted(allowed_owners))


def run_all_tests():
    tests = [
        ("K_C1_no_persistent_writeback", test_k_c1_no_persistent_writeback),
        ("K_C2_deepcopy_before_merge", test_k_c2_deepcopy_before_merge),
        ("K_C3_session_info_both_layers", test_k_c3_session_info_both_layers),
        ("K_C4_diagnostic_emits_both_labels", test_k_c4_diagnostic_emits_both_labels),
        ("K_C5_diagnostic_only_when_differ", test_k_c5_diagnostic_only_when_differ),
        ("K_C6_diagnostic_in_try_except", test_k_c6_diagnostic_in_try_except),
        ("K_C7_no_persistent_writeback_global", test_k_c7_no_persistent_writeback_global),
    ]
    passed = 0
    failed = 0
    for name, fn in tests:
        print("Test: %s" % name)
        try:
            fn()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            print("  ERROR: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
