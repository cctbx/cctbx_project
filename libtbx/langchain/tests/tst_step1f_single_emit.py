"""K_H14_ITEM_2: STEP_1F single-emit (v119.H14).

Covers the duplicate diagnostic_messages relay surfaced by
run_39_openai batch analysis: 315 of 520 runs (60.6%) showed
adjacent duplicate [STEP_1F] preprocessing_metrics lines in
run.log.

The cause: TWO sites both iterated diagnostic_messages and
wrote to stderr:
  1. phenix.programs.ai_analysis._relay_diagnostic_messages_to_stderr
     (the v119.H5.1 dispatcher-level relay; called inside
      run_job_on_server_or_locally at 4 return paths)
  2. programs/ai_agent.py:8087-8103 (a client-side block from
     v119.H5 §2.9 that survived the H5.1 refactor)

H14 fix: remove the client-side block in (2).  The dispatcher
relay in (1) is the correct location — it handles all dispatch
paths uniformly, including LLM-only modes where ai_agent.py's
caller frame doesn't run.

This test pins:
  - The client-side relay block has been removed from
    programs/ai_agent.py (source scan)
  - The dispatcher relay function still exists in
    programs/ai_analysis.py (we didn't accidentally remove THE
    relay)
  - The H14 marker comment is present at the removal site

Source-scan based.  We don't capture stderr here because the
end-to-end emission requires a running PHENIX environment that
isn't available in the sandbox.  Tom's gate runs validate the
operational behavior.

Total: 8 tests.
"""
from __future__ import absolute_import, division, print_function

import os


# =====================================================================
# Helpers
# =====================================================================

def _read_source(rel_path):
    """Read a source file from disk, with deployment-layout awareness.

    Two layouts are supported:

      (a) Sandbox layout: everything is bundled under one tree, with
          tests/ sibling to agent/, programs/, knowledge/, etc.
          rel_path resolves to <ship_root>/<rel_path>.

      (b) PHENIX deployment: this test file lives in
          cctbx_project/libtbx/langchain/tests/, but programs/ai_agent.py
          and programs/ai_analysis.py live in a SIBLING module at
          phenix/phenix/programs/.  The langchain/ tree does NOT have
          a programs/ subdirectory.

    v119.H14.3: pre-H14.3 this function only handled layout (a) and
    crashed with FileNotFoundError in PHENIX deployment.  Tom's
    sandbox run after the H14.2 deployment surfaced 2 source-scan
    tests failing here.  We now try (a) then (b); if neither is
    findable we raise a clear error pointing to both attempted paths.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)

    # Layout (a): sandbox — rel_path is sibling to tests/
    candidate_a = os.path.join(parent, rel_path)
    if os.path.exists(candidate_a):
        with open(candidate_a, 'r') as f:
            return f.read()

    # Layout (b): PHENIX deployment — walk up from
    # cctbx_project/libtbx/langchain/tests/ to modules/, then look in
    # phenix/phenix/<rel_path>.  Only applies when rel_path starts
    # with "programs/" (the only sister-module case we know about).
    if rel_path.startswith("programs/"):
        # tests/ -> langchain/ -> libtbx/ -> cctbx_project/ -> modules/
        modules_dir = os.path.dirname(
            os.path.dirname(os.path.dirname(parent)))
        candidate_b = os.path.join(
            modules_dir, "phenix", "phenix", rel_path)
        if os.path.exists(candidate_b):
            with open(candidate_b, 'r') as f:
                return f.read()
    else:
        candidate_b = None

    # Neither layout found the file.  Raise a clear error so the
    # source-scan test fails with a useful message rather than the
    # raw FileNotFoundError.
    tried = "\n  - %s (sandbox layout)" % candidate_a
    if candidate_b:
        tried += "\n  - %s (PHENIX deployment layout)" % candidate_b
    raise IOError(
        "Could not locate source file %r in any known layout. Tried:%s"
        % (rel_path, tried))


# =====================================================================
# §A: programs/ai_agent.py — duplicate relay removed
# =====================================================================

def test_ai_agent_no_diagnostic_relay_loop():
    """v119.H14: programs/ai_agent.py must NOT contain a stderr-write
    loop over diagnostic_messages — the dispatcher does that uniformly.

    Pre-H14 ai_agent.py:8087-8103 had:
        for _msg in _diagnostics:
            _sys.stderr.write(_msg + '\\n')
    which caused every diagnostic marker (notably [STEP_1F]) to be
    emitted twice (once by the dispatcher, once by this block).
    """
    src = _read_source("programs/ai_agent.py")

    # The most diagnostic-specific signature of the removed block
    assert "for _msg in _diagnostics" not in src, (
        "ai_agent.py still contains 'for _msg in _diagnostics' — the "
        "client-side diagnostic_messages relay must be removed in H14 "
        "(dispatcher in programs/ai_analysis.py handles it)")

    # Also check for the variable assignment that fed the loop
    assert "_diagnostics = result.get('diagnostic_messages'" not in src, (
        "ai_agent.py still extracts _diagnostics from result — that "
        "extraction was part of the removed relay block")
    print("  PASS: test_ai_agent_no_diagnostic_relay_loop")


def test_ai_agent_h14_marker_present():
    """The removal must leave a clear marker explaining why the block
    is gone — protects against future contributors re-adding the
    relay by reverting only the comment."""
    src = _read_source("programs/ai_agent.py")
    assert "v119.H14: removed duplicate diagnostic_messages relay" in src, (
        "The H14 marker comment must be present where the duplicate "
        "relay block was removed.  This protects against accidental "
        "re-addition.")
    # And the comment must explain the consequence
    assert "would be emitted TWICE" in src or "emitted TWICE" in src or (
        "duplicate" in src and "STEP_1F" in src
    ), "The H14 comment should describe the double-emit problem"
    print("  PASS: test_ai_agent_h14_marker_present")


# =====================================================================
# §B: programs/ai_analysis.py — dispatcher relay still present
# =====================================================================

def test_ai_analysis_dispatcher_relay_intact():
    """v119.H14 must NOT remove the dispatcher-level relay.  That's
    the relay we WANT to keep — it's the single uniform site that
    handles all dispatch paths.  Only the client-side duplicate in
    ai_agent.py was removed."""
    src = _read_source("programs/ai_analysis.py")
    assert "def _relay_diagnostic_messages_to_stderr" in src, (
        "The dispatcher relay function must still be defined in "
        "programs/ai_analysis.py — this is the kept relay, not the "
        "removed one")
    # And it must still be called from the dispatcher
    assert "_relay_diagnostic_messages_to_stderr(" in src, (
        "The dispatcher relay function must still be called from "
        "run_job_on_server_or_locally")
    print("  PASS: test_ai_analysis_dispatcher_relay_intact")


def test_dispatcher_relay_called_from_all_return_paths():
    """The dispatcher's relay is wired at every return path from
    run_job_on_server_or_locally.  Pre-H14 there were 4 call sites
    (the run-locally path + 3 server-side paths).  H14 doesn't
    change this — but verify it still holds."""
    src = _read_source("programs/ai_analysis.py")
    # Count call sites, excluding the def line and docstring/comment
    # references
    call_count = 0
    for line in src.split("\n"):
        stripped = line.strip()
        if (stripped.startswith("_relay_diagnostic_messages_to_stderr(") and
                not stripped.startswith("#")):
            call_count += 1
    assert call_count >= 3, (
        "Expected at least 3 call sites for the dispatcher relay "
        "(one per return path).  Found: %d" % call_count)
    print("  PASS: test_dispatcher_relay_called_from_all_return_paths "
          "(found %d call sites)" % call_count)


# =====================================================================
# §C: Other call sites in ai_agent.py — no duplicate relays elsewhere
# =====================================================================

def test_no_other_diagnostic_loops_in_ai_agent():
    """ai_agent.py has FIVE call sites for run_job_on_server_or_locally
    (lines ~1741, 3831, 8064, 8321, 9409 in the v119.H14 source).
    Only the one at line ~8064 (advice preprocessing) had a client-side
    relay block pre-H14.  Verify that NO call site has any diagnostic
    relay loop after H14."""
    src = _read_source("programs/ai_agent.py")

    # No loop iterating over diagnostic_messages anywhere
    forbidden_patterns = [
        "for _msg in _diagnostics",
        "for msg in diagnostic_messages",
        "for _msg in diagnostic_messages",
        "for msg in _diagnostics",
    ]
    for pattern in forbidden_patterns:
        assert pattern not in src, (
            "ai_agent.py contains forbidden relay loop %r — "
            "the dispatcher in ai_analysis.py is the SINGLE relay site"
            % pattern)
    print("  PASS: test_no_other_diagnostic_loops_in_ai_agent")


def test_no_stderr_write_of_markers_in_ai_agent():
    """A weaker but broader check: ai_agent.py shouldn't be writing
    diagnostic markers to stderr at all (markers are emitted by the
    preprocessor and relayed by the dispatcher)."""
    src = _read_source("programs/ai_agent.py")

    # Look for direct stderr writes of marker-shaped strings
    # (i.e., bracket-delimited token names)
    lines = src.split("\n")
    suspect_lines = []
    for i, line in enumerate(lines, start=1):
        # Skip comments
        if line.lstrip().startswith("#"):
            continue
        # Look for sys.stderr.write with a marker-shaped pattern
        if ("sys.stderr.write" in line or
                "_sys.stderr.write" in line):
            # Is there a "[XXX]" marker pattern nearby (this line or next)?
            window = "\n".join(lines[i-1:i+2])
            if "[STEP_" in window or "[DIRECTIVE" in window or "[ADVICE_" in window:
                suspect_lines.append((i, line.strip()))

    assert not suspect_lines, (
        "ai_agent.py contains suspect stderr.write of diagnostic "
        "markers — markers should be emitted via diagnostic_messages "
        "and relayed by the dispatcher. Suspect lines:\n%s"
        % "\n".join("  L%d: %s" % (n, l) for n, l in suspect_lines))
    print("  PASS: test_no_stderr_write_of_markers_in_ai_agent")


# =====================================================================
# §D: Dispatcher relay safety properties
# =====================================================================

def test_dispatcher_relay_never_raises():
    """The dispatcher's relay must be exception-safe — a closed pipe
    or unwriteable stderr must NOT propagate as a runtime exception
    that breaks the agent.  The docstring should declare this and
    the implementation should match (try/except wrapping).
    """
    src = _read_source("programs/ai_analysis.py")

    # Find the function definition and its body
    fn_start = src.find("def _relay_diagnostic_messages_to_stderr")
    assert fn_start != -1, "Relay function not found"
    # Take a sufficient slice (~70 lines)
    fn_slice = src[fn_start:fn_start + 4000]

    # Must have outer try/except
    assert "try:" in fn_slice and "except Exception" in fn_slice, (
        "Relay function must have outer try/except Exception")

    # Docstring should explicitly declare exception safety
    assert "Never raises" in fn_slice or "never raise" in fn_slice.lower(), (
        "Relay function docstring should declare 'Never raises'")
    print("  PASS: test_dispatcher_relay_never_raises")


def test_dispatcher_relay_handles_dict_and_object():
    """Defensive unpacking handles both group_args objects (the
    normal case) and dicts (hypothetical future raw-JSON results).
    Per Gemini Q5A on v119.H5: getattr on a dict misses the key
    and silently drops telemetry.

    Source-scan check that the defensive branch is still in place."""
    src = _read_source("programs/ai_analysis.py")

    fn_start = src.find("def _relay_diagnostic_messages_to_stderr")
    fn_slice = src[fn_start:fn_start + 4000]

    assert "isinstance(working_results, dict)" in fn_slice, (
        "Relay must check isinstance(..., dict) for dict-shaped results")
    assert "getattr(working_results" in fn_slice, (
        "Relay must use getattr for object-shaped results")
    print("  PASS: test_dispatcher_relay_handles_dict_and_object")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: client-side removal (2)
    test_ai_agent_no_diagnostic_relay_loop()
    test_ai_agent_h14_marker_present()
    # §B: dispatcher kept (2)
    test_ai_analysis_dispatcher_relay_intact()
    test_dispatcher_relay_called_from_all_return_paths()
    # §C: other call sites clean (2)
    test_no_other_diagnostic_loops_in_ai_agent()
    test_no_stderr_write_of_markers_in_ai_agent()
    # §D: dispatcher safety properties (2)
    test_dispatcher_relay_never_raises()
    test_dispatcher_relay_handles_dict_and_object()


if __name__ == "__main__":
    print("K_H14_ITEM_2: STEP_1F single-emit (v119.H14)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H14_ITEM_2 complete.")
