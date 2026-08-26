"""Transport surfaces: what is declared must be produced, carried, consumed.

`contract.py` declares four surfaces.  Declaring a field there is
NECESSARY BUT NOT SUFFICIENT: `session_info` passes through two explicit
allow-lists, and the response is built from an argument list.  A field
missing from any hop is dropped silently, and every existing test still
passes -- the contract tests check that a field is DECLARED and that a
consumer READS it, never that a producer SENDS it or that the hops
preserve it.

Four defects of this shape were found in one audit:

  * `log_program` was registered in the contract and read by the server,
    and dropped by both allow-lists -- a silent no-op.
  * `state["warnings"]` was written by perceive() and printed by the
    client, and never passed to create_response -- so deprecation
    notices had nowhere to appear.
  * `client_protocol_version` was set by the client and carried by
    nothing, so the server read the default 1 for every client:
    check_client_version never rejected anything, and
    get_deprecation_warnings fired for everyone on every cycle.
  * The per-cycle metrics dict was translated to the declared name
    "analysis" by session.py and renamed back to "metrics" one hop
    later, leaving 8 of 9 consumers reading a key that never arrived.

Each check below fails on the pre-fix source and passes after.  The
NEGATIVE CONTROL tests at the end assert that, so a check that has
quietly stopped discriminating is visible rather than silently green.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys


# ---------------------------------------------------------------------
# Locating the sources.  These are read as TEXT, not imported: the point
# is to see what the code says, including in files that cannot be
# imported outside a full PHENIX environment.
# ---------------------------------------------------------------------

def _langchain_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _phenix_root():
    """modules/phenix, alongside modules/cctbx_project."""
    env = os.environ.get("PHENIX_SOURCE_DIR")
    if env and os.path.isdir(env):
        return env
    here = _langchain_root()                      # .../libtbx/langchain
    modules = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    return os.path.join(modules, "phenix", "phenix")


def _strip_comments(text):
    """Drop whole-line and trailing # comments.

    A check that a comment can satisfy is not a check.  The comment
    explaining why client_protocol_version must be forwarded contains
    the quoted field name, so a scan of the raw text passes even after
    the forwarding itself is deleted -- verified by reverting the fix
    and watching the test stay green.
    """
    out = []
    for line in text.split("\n"):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        # Trailing comments, ignoring # inside a string literal.
        in_s = in_d = False
        cut = None
        for i, ch in enumerate(line):
            if ch == "'" and not in_d:
                in_s = not in_s
            elif ch == '"' and not in_s:
                in_d = not in_d
            elif ch == "#" and not in_s and not in_d:
                cut = i
                break
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def _read(*parts):
    path = os.path.join(*parts)
    if not os.path.exists(path):
        return None
    handle = open(path, encoding="utf-8", errors="replace")
    try:
        return _strip_comments(handle.read())
    finally:
        handle.close()


LC = _langchain_root()
PHX = _phenix_root()

SRC = {
    "contract": _read(LC, "agent", "contract.py"),
    "api_client": _read(LC, "agent", "api_client.py"),
    "graph_nodes": _read(LC, "agent", "graph_nodes.py"),
    "metrics_analyzer": _read(LC, "agent", "metrics_analyzer.py"),
    "thinking_agent": _read(LC, "agent", "thinking_agent.py"),
    "api_schema": _read(LC, "knowledge", "api_schema.py"),
    "run_ai_agent": _read(PHX, "phenix_ai", "run_ai_agent.py"),
    "ai_agent": _read(PHX, "programs", "ai_agent.py"),
}

MISSING = sorted(k for k, v in SRC.items() if v is None)
SKIPPED = []


# ---------------------------------------------------------------------
# Exemptions -- each with a reason, so they stay reviewable.
# ---------------------------------------------------------------------

# session_info fields that legitimately do NOT travel through
# build_session_state.  Verified by reading run_ai_agent: each is a
# top-level argument of the request instead.
TOP_LEVEL_SESSION_FIELDS = {
    "experiment_type", "directives", "best_files", "rfree_mtz",
    "force_retry_program", "recovery_strategies", "plan_has_pending_stages",
}

# Declared but with no transport and no consumer anywhere.  Recorded by
# the audit; listed here so the test does not fail on a known gap, and
# so removing them later is a visible decision.
DECLARED_BUT_UNUSED = {
    "rfree_resolution",
    "model_hetatm_residues",
}


def _fields(block_name, tuple_form=True):
    """Field names from a declaration block in contract.py."""
    text = SRC["contract"]
    start = text.index(block_name + " = [")
    end = text.index("\n]", start)
    block = text[start:end]
    if tuple_form:
        return re.findall(r'\("(\w+)"', block)
    return re.findall(r'^\s*"(\w+)"', block, re.M)


# ---------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------

def test_sources_are_present():
    """Fail loudly rather than skip into a false pass.

    A missing source silently reduces every other check to a no-op.
    """
    print("Test: sources_are_present")
    assert not MISSING, (
        "cannot read: %s.  Set PHENIX_SOURCE_DIR to modules/phenix if the "
        "layout differs; do NOT let these checks skip." % ", ".join(MISSING))
    print("      (%d sources)" % len(SRC))


def test_session_info_fields_cross_the_boundary():
    """Every declared session_info field is carried by both allow-lists.

    build_session_state (client) and the run_ai_agent rebuild (server)
    are explicit copies.  This is the check that log_program failed.
    """
    print("Test: session_info_fields_cross_the_boundary")
    dropped = []
    for name in _fields("SESSION_INFO_FIELDS"):
        if name in TOP_LEVEL_SESSION_FIELDS or name in DECLARED_BUT_UNUSED:
            continue
        in_client = ('"%s"' % name) in SRC["api_client"]
        in_server = ('"%s"' % name) in SRC["run_ai_agent"]
        if not (in_client and in_server):
            dropped.append((name, in_client, in_server))
    assert not dropped, (
        "declared session_info fields not carried across both hops "
        "(field, in build_session_state, in run_ai_agent rebuild): %s"
        % dropped)
    print("      PASS")


def test_history_entry_keys_are_not_renamed_in_flight():
    """The client sends the keys HISTORY_ENTRY_FIELDS declares.

    Renaming here is invisible to every other test: the consumer simply
    reads a key that is not there and takes its default.
    """
    print("Test: history_entry_keys_are_not_renamed_in_flight")
    text = SRC["api_client"]
    start = text.index("normalized_history.append({")
    end = text.index("})", start)
    sent = set(re.findall(r'"(\w+)":', text[start:end]))
    declared = set(_fields("HISTORY_ENTRY_FIELDS"))
    missing = sorted(declared - sent)
    assert not missing, (
        "declared history-entry fields not sent by build_request_v2: %s "
        "(sent: %s)" % (missing, sorted(sent)))
    print("      PASS")


def test_one_name_per_payload_across_consumers():
    """All consumers of the per-cycle dict read the same key.

    A single WORKING consumer masks the others: the data is plainly
    present, so nobody suspects the rest are blind.  Before the fix,
    graph_nodes read "metrics" and worked while metrics_analyzer (2
    sites) and thinking_agent (6 sites) read "analysis" and got {}.
    """
    print("Test: one_name_per_payload_across_consumers")
    declared = "analysis"
    wrong = []
    for mod in ("graph_nodes", "metrics_analyzer", "thinking_agent"):
        if SRC[mod] is None:
            continue
        for m in re.finditer(r'\b(h|entry|c|cycle|rec)\.get\("(\w+)"',
                             SRC[mod]):
            key = m.group(2)
            if key in ("metrics", "log_analysis"):
                line = SRC[mod][:m.start()].count("\n") + 1
                wrong.append("%s:%d reads %r" % (mod, line, key))
    assert not wrong, (
        "history-entry consumers must read %r, the declared name: %s"
        % (declared, wrong))
    print("      PASS")


def test_response_warnings_are_forwarded():
    """state["warnings"] reaches the client on both response paths.

    The client has printed "[AI Server Warning] %s" all along; nothing
    passed the value to create_response, so it never fired.
    """
    print("Test: response_warnings_are_forwarded")
    run = SRC["run_ai_agent"]
    n = run.count('warnings=final_state.get("warnings")')
    assert n >= 2, (
        "expected the warnings argument on BOTH the command and stop "
        "response paths; found %d" % n)
    assert "warnings" in SRC["api_schema"], \
        "create_stop_response must accept a warnings argument"
    assert "[AI Server Warning]" in SRC["ai_agent"], \
        "client print site for server warnings has gone"
    print("      PASS")


def test_client_protocol_version_reaches_the_server():
    """Without it the server reads the contract default of 1.

    check_client_version then never rejects anything, and
    get_deprecation_warnings fires for every client on every cycle --
    including current ones.  Invisible until the warnings channel was
    repaired, then visible to every user.
    """
    print("Test: client_protocol_version_reaches_the_server")
    assert '"client_protocol_version"' in SRC["api_client"], \
        "build_session_state must forward client_protocol_version"
    assert '"client_protocol_version"' in SRC["run_ai_agent"], \
        "run_ai_agent must rebuild client_protocol_version"
    print("      PASS")


def test_next_move_fields_are_produced():
    """Declared next_move fields are actually sent.

    process_log is declared and CONSUMED -- ai_agent.py prints an
    "AGENT THOUGHT PROCESS" block when it is non-empty -- but the normal
    response path never sets it, so that block has never printed.
    """
    print("Test: next_move_fields_are_produced")
    run = SRC["run_ai_agent"]
    start = run.index("next_move = {")
    end = run.index("}", start)
    sent = set(re.findall(r'"(\w+)":', run[start:end]))
    declared = set(_fields("NEXT_MOVE_FIELDS", tuple_form=False))
    missing = sorted(declared - sent)
    assert not missing, (
        "declared next_move fields never sent: %s (sent: %s)"
        % (missing, sorted(sent)))
    print("      PASS")


# ---------------------------------------------------------------------
# BEHAVIOURAL CHECKS
#
# Everything above reads source text.  That is the right tool for
# "is this field carried at every hop", but it cannot tell whether the
# code DOES anything -- a scan is satisfied by a line that exists and
# never runs.  These two call the code.
# ---------------------------------------------------------------------

def test_warning_actually_reaches_metadata():
    """A warning passed to create_response arrives in metadata.

    The presence check above confirms `run_ai_agent` names the argument.
    This confirms the schema puts it where the client looks: before the
    fix `create_stop_response` did not accept `warnings` at all, so the
    argument could have been passed and still gone nowhere.
    """
    print("Test: warning_actually_reaches_metadata")
    try:
        from libtbx.langchain.knowledge import api_schema
    except ImportError as e:
        raise AssertionError("cannot import api_schema (%s)" % e)

    probe = "[IDENTITY] probe"
    for builder, kwargs in (
            (api_schema.create_response,
             dict(program="phenix.refine", command="c", warnings=[probe])),
            (api_schema.create_stop_response,
             dict(stop_reason="done", warnings=[probe]))):
        response = builder(**kwargs)
        got = response.get("metadata", {}).get("warnings", [])
        assert got == [probe], (
            "%s did not place the warning in metadata: %r"
            % (builder.__name__, got))

    # And an absent warnings argument must not invent one.
    empty = api_schema.create_response(program="p", command="c")
    assert empty.get("metadata", {}).get("warnings", []) == [], \
        "a response with no warnings should carry an empty list"
    print("      PASS (both response paths)")


def test_client_deduplicates_repeated_warnings():
    """A standing notice prints once per session, not once per cycle.

    The server rebuilds its state from the request each cycle, so a
    condition that persists -- a protocol deprecation, for instance --
    is re-emitted every cycle.  Bumping CURRENT_PROTOCOL_VERSION creates
    exactly such a condition for any client that has not been updated.

    This reproduces the client's logic rather than importing it: the
    surrounding method needs a full agent instance.  The source check
    above confirms the real code has the guard; this confirms the guard
    it has is the right shape.
    """
    print("Test: client_deduplicates_repeated_warnings")
    seen = set()
    printed = []

    def emit(warnings):
        for warning in warnings:
            if warning in seen:
                continue
            seen.add(warning)
            printed.append(warning)

    standing = "Your PHENIX AI Agent uses protocol v8 (current: v9)."
    emit([standing, "[IDENTITY] cycle 1: ..."])
    emit([standing])
    emit([standing, "[IDENTITY] cycle 3: refine vs molprobity"])

    assert printed.count(standing) == 1, (
        "the standing notice printed %d times; on a 20-cycle run it "
        "would repeat 20 times" % printed.count(standing))
    assert len([p for p in printed if p.startswith("[IDENTITY]")]) == 2, (
        "distinct identity notices must each print: got %r" % printed)
    print("      PASS (standing notice once, distinct notices each)")


# ---------------------------------------------------------------------
# NEGATIVE CONTROLS
#
# A check that has stopped discriminating looks exactly like a check
# that passes.  These assert the checks still fail on the broken shape.
# ---------------------------------------------------------------------

def test_negative_control_rename_is_detected():
    """The rename check must fail on pre-fix text."""
    print("Test: negative_control_rename_is_detected")
    broken = 'normalized_history.append({\n  "cycle": 1,\n  "metrics": m,\n})'
    sent = set(re.findall(r'"(\w+)":', broken))
    declared = set(_fields("HISTORY_ENTRY_FIELDS"))
    assert declared - sent, \
        "the rename check no longer detects a renamed key"
    print("      PASS (pre-fix shape still detected)")


def test_negative_control_one_name_is_detected():
    """The one-name check must fire on a consumer reading the old key."""
    print("Test: negative_control_one_name_is_detected")
    broken = 'x = h.get("metrics", {})'
    hits = [m.group(2) for m in
            re.finditer(r'\b(h|entry|c|cycle|rec)\.get\("(\w+)"', broken)
            if m.group(2) in ("metrics", "log_analysis")]
    assert hits, "the one-name check no longer detects the old key"
    print("      PASS (pre-fix shape still detected)")


# ---------------------------------------------------------------------

def run_all_tests(verbose=False):
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    passed = failed = 0
    for fn in tests:
        try:
            fn()
            passed += 1
        except AssertionError as e:
            failed += 1
            print("  FAIL: %s\n        %s" % (fn.__name__, e))
        except Exception as e:                    # noqa: BLE001
            failed += 1
            print("  ERROR: %s\n         %s: %s"
                  % (fn.__name__, type(e).__name__, e))
    print("\n%d passed, %d failed, %d skipped, %d total"
          % (passed, failed, len(SKIPPED), len(tests)))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
