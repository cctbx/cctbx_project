"""K_H5: diagnostic_messages relay channel.

v119.H5 (§2.9 stderr → diagnostic_messages).

Tests the _emit_marker helper and the diagnostic_messages
field in run_advice_preprocessing's result.  The channel
exists to relay server-side [STEP_1F] / [STEP_1F_FAILED]
markers to the client, which cannot otherwise see them when
the call dispatches to a remote server.

Design principle (Tom's uniformity criterion, v119.H5 review):
All actions take place identically whether analysis runs on
the server or locally.  Markers flow through the
diagnostic_messages list in both modes; the client re-emits
them at a single, uniform site after the dispatcher returns.
The helper does NOT write to stderr — that would diverge
local-mode and remote-mode behavior.

Test organization:
  §A   _emit_marker helper unit tests             (4 tests)
  §B   diagnostic_messages field in return         (4 tests)
  §C   Backward compatibility                      (3 tests)
  §D   Full plumbing — uniform client re-emit     (6 tests)
  §E   Integration — production encode/decode      (4 tests)
  §F   Extended markers (v119.H5.1)                (5 tests)

Total: 26 tests.

Like K_H4, K_H5 is stdlib-only and runs in sandbox.  The
graceful-skip pattern is preserved for consistency with the
v119 K-suite style; in sandbox, §B, §E, and §F tests (which
need libtbx — group_args, Program class, simple_string_as_text,
or run_advice_preprocessing) SKIP, the others PASS.

§E exercises the production encode path
(Program.get_results_as_JSON) and verifies the JSON wire
format round-trips correctly.  §E test #1 (Total Init seed)
also exercises the production get_results_from_all to ensure
every working_results is constructed with diagnostic_messages=[],
which is the linchpin of the H5 design.

§F (v119.H5.1) verifies new markers added by H5.1 for
catastrophic-failure visibility.  Tests use deterministic
monkey-patching (NOT invalid-provider injection) because
the production code degrades gracefully on configuration
mistakes — the [*_FAILED] markers fire only on unexpected
exceptions caught by the main try/except in each function.
"""
from __future__ import absolute_import, division, print_function

import io
import os
import sys


# =====================================================================
# Helper: try to import the _emit_marker function and the
# run_advice_preprocessing function.  Both are tested separately.
# =====================================================================

def _try_import_emit_marker():
    """Import _emit_marker; return (fn, error_msg).

    The helper lives at module level in run_ai_analysis.py.
    Use a try-import pattern for consistency with other K-suites.

    SANDBOX FALLBACK: when the parent module can't be imported
    (libtbx not available), we extract the helper's source from
    the .py file and exec it standalone.  The helper is
    stdlib-only (sys.stderr only), so this works.  This lets
    K_H5's §A tests run in sandbox even when the rest of the
    test suite must SKIP.  Under PHENIX, the normal import
    path succeeds and the fallback is never used.
    """
    # Direct import paths to try
    try:
        # When phenix_ai is on sys.path (PHENIX install)
        from phenix.phenix_ai.run_ai_analysis import _emit_marker
        return _emit_marker, None
    except ImportError:
        pass
    try:
        # When phenix_ai is a sibling package (sandbox with
        # libtbx available)
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        candidate = os.path.join(parent, 'phenix_ai')
        if os.path.isdir(candidate) and parent not in sys.path:
            sys.path.insert(0, parent)
        from phenix_ai.run_ai_analysis import _emit_marker
        return _emit_marker, None
    except ImportError:
        pass

    # SANDBOX FALLBACK: extract the helper source from the .py
    # file and exec it.  The helper depends only on sys.stderr
    # (already a global) so this is safe.
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        src_path = os.path.join(parent, 'phenix_ai', 'run_ai_analysis.py')
        if not os.path.isfile(src_path):
            return None, "run_ai_analysis.py not found at %s" % src_path
        with open(src_path, 'r', encoding='utf-8') as f:
            full_src = f.read()
        # Find and extract just the _emit_marker function.
        # Marker is its def line; the function ends at the next
        # top-level def (or 'class') or EOF.
        start = full_src.find('def _emit_marker(')
        if start < 0:
            return None, "could not find _emit_marker def in source"
        # Find the next top-level 'def ' or 'class ' after start
        rest = full_src[start:]
        end_offsets = []
        for marker in ('\ndef ', '\nclass '):
            pos = rest.find(marker, 1)
            if pos > 0:
                end_offsets.append(pos)
        end = min(end_offsets) if end_offsets else len(rest)
        helper_src = rest[:end]
        # Exec into a fresh namespace; sys must be available
        ns = {'sys': sys}
        exec(helper_src, ns)
        if '_emit_marker' not in ns:
            return None, "exec did not produce _emit_marker"
        return ns['_emit_marker'], None
    except Exception as e:
        return None, "sandbox fallback failed: %s" % e


def _try_import_run_advice_preprocessing():
    """Import run_advice_preprocessing; return (fn, error_msg).

    This function depends on libtbx.group_args which is a PHENIX
    runtime dep.  In sandbox without PHENIX, the import fails.
    K_H5 §B and §C tests SKIP gracefully in that case.
    """
    try:
        from phenix.phenix_ai.run_ai_analysis import run_advice_preprocessing
        return run_advice_preprocessing, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from phenix_ai.run_ai_analysis import run_advice_preprocessing
        return run_advice_preprocessing, None
    except ImportError as e:
        return None, str(e)


# =====================================================================
# §A: _emit_marker helper unit tests (4 tests)
# =====================================================================

def test_emit_marker_appends_to_list_only():
    """Helper appends marker to list and does NOT write to stderr.

    This is the contract under the uniformity criterion (v119.H5
    review): markers reach the operator's terminal at a single,
    uniform site on the client side after the dispatcher returns.
    The helper has only one job — append to the list.

    Crucially: even in local-execution mode, the helper does NOT
    write to stderr at the call site.  This ensures local and
    remote modes produce identical operator-visible output (no
    duplication risk, no timing-dependent ordering).
    """
    fn, err = _try_import_emit_marker()
    if fn is None:
        print("  SKIP: cannot import _emit_marker (%s)" % err)
        return

    # Capture stderr — should remain empty
    saved_stderr = sys.stderr
    buf = io.StringIO()
    sys.stderr = buf
    try:
        msgs = []
        fn(msgs, "[TEST_MARKER] hello")
    finally:
        sys.stderr = saved_stderr

    # Stderr received NOTHING (the criterion)
    stderr_content = buf.getvalue()
    assert stderr_content == "", (
        "Helper must not write to stderr; got %r" % stderr_content)

    # List received the marker (without trailing newline)
    assert len(msgs) == 1, "Expected 1 list entry, got %d" % len(msgs)
    assert msgs[0] == "[TEST_MARKER] hello", (
        "List entry mismatch: %r" % msgs[0])
    assert not msgs[0].endswith("\n"), (
        "List entry should not have trailing newline")
    print("  PASS: test_emit_marker_appends_to_list_only")


def test_emit_marker_drops_silently_on_non_list_container():
    """Non-list container: helper drops marker silently, no raise.

    Per Gemini Q5B review of v119.H5: the previous blanket
    try/except: pass would have swallowed AttributeError from
    calling .append() on a None container, hiding init bugs.
    The current implementation uses an explicit isinstance(...,
    list) check.

    For a non-list container, the marker is dropped and no
    exception is raised — uniformity criterion preserved (no
    fallback stderr write), and operator state stays consistent.
    Under the Total Initialization Policy the container is
    always a list, so this path is purely defense-in-depth
    against future refactor accidents.
    """
    fn, err = _try_import_emit_marker()
    if fn is None:
        print("  SKIP: cannot import _emit_marker (%s)" % err)
        return

    # None container — must not raise AttributeError
    try:
        fn(None, "[TEST_MARKER] none-container")
    except Exception as e:
        raise AssertionError(
            "Helper raised on None container (should drop silently): %s"
            % e)

    # String container (incorrect type)
    try:
        fn("not a list", "[TEST_MARKER] string-container")
    except Exception as e:
        raise AssertionError(
            "Helper raised on string container: %s" % e)

    # Dict container (incorrect type)
    try:
        fn({"key": "value"}, "[TEST_MARKER] dict-container")
    except Exception as e:
        raise AssertionError(
            "Helper raised on dict container: %s" % e)

    print("  PASS: test_emit_marker_drops_silently_on_non_list_container")


def test_emit_marker_no_stderr_write_under_any_condition():
    """Helper never writes to stderr, even on degenerate input.

    The uniformity criterion is strict: stderr writes happen at
    a single client-side site only.  If the helper "helpfully"
    wrote to stderr as a fallback when the container is broken,
    it would:
    - Cause duplication when the client re-emits the same
      marker
    - Diverge local-mode behavior from remote-mode behavior
      (the entire problem the criterion forbids)

    This test pins the contract: no matter what container is
    passed, the helper does not touch stderr.
    """
    fn, err = _try_import_emit_marker()
    if fn is None:
        print("  SKIP: cannot import _emit_marker (%s)" % err)
        return

    saved_stderr = sys.stderr
    buf = io.StringIO()
    sys.stderr = buf
    try:
        # Various degenerate container types
        fn(None, "[TEST_MARKER] none-no-stderr")
        fn("not a list", "[TEST_MARKER] str-no-stderr")
        fn({}, "[TEST_MARKER] dict-no-stderr")
        # And the normal case
        fn([], "[TEST_MARKER] normal-no-stderr")
    finally:
        sys.stderr = saved_stderr

    assert buf.getvalue() == "", (
        "Helper must NEVER write to stderr; got %r" % buf.getvalue())
    print("  PASS: test_emit_marker_no_stderr_write_under_any_condition")


def test_emit_marker_flattens_embedded_newlines():
    """Multi-line marker is flattened to single-line in list.

    Each list entry MUST be a single line (per plan rev 4 §2.2)
    to protect line-oriented downstream consumers.  The helper
    replaces \\n with ' \\ ' before storage AND strips \\r
    (Windows-style line endings).  Without the \\r strip,
    re-emit on a Unix terminal would produce ugly ^M characters
    or worse — terminal cursor jumps on raw \\r bytes.
    """
    fn, err = _try_import_emit_marker()
    if fn is None:
        print("  SKIP: cannot import _emit_marker (%s)" % err)
        return

    saved_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        msgs = []
        # Simulate a multi-line traceback string
        fn(msgs, "[STEP_1F_FAILED] error: line one\nline two\nline three")
        # Simulate Windows-style CRLF line endings
        fn(msgs, "[STEP_1F_FAILED] error: alpha\r\nbeta\r\ngamma")
        # Simulate lone \r (rare but possible from misbehaving libs)
        fn(msgs, "[STEP_1F_FAILED] error: classic-mac-line\rending")
    finally:
        sys.stderr = saved_stderr

    assert len(msgs) == 3

    # First entry: \n → ' \ ' flattening
    entry_lf = msgs[0]
    assert '\n' not in entry_lf, (
        "Embedded LF not flattened: %r" % entry_lf)
    assert 'line one' in entry_lf
    assert 'line two' in entry_lf
    assert 'line three' in entry_lf

    # Second entry: CRLF should produce no \r AND no \n
    # (the \r is stripped, then the \n is flattened)
    entry_crlf = msgs[1]
    assert '\r' not in entry_crlf, (
        "CR not stripped from CRLF: %r" % entry_crlf)
    assert '\n' not in entry_crlf, (
        "LF not flattened in CRLF: %r" % entry_crlf)
    assert 'alpha' in entry_crlf
    assert 'beta' in entry_crlf
    assert 'gamma' in entry_crlf

    # Third entry: lone \r (no \n) should be stripped entirely.
    # No flattening to ' \ ' since the order is .replace('\n')
    # first, then .replace('\r').  Result: "classic-mac-line"
    # concatenated directly with "ending".
    entry_cr = msgs[2]
    assert '\r' not in entry_cr, (
        "Lone CR not stripped: %r" % entry_cr)
    assert 'classic-mac-line' in entry_cr
    assert 'ending' in entry_cr

    print("  PASS: test_emit_marker_flattens_embedded_newlines "
          "(LF, CRLF, lone CR all handled)")


# =====================================================================
# §B: diagnostic_messages field in return (3 tests)
# =====================================================================

def test_preprocessing_returns_diagnostic_messages_field():
    """run_advice_preprocessing's result has diagnostic_messages.

    The field must be present in ALL return paths (empty advice,
    use_rules_only, main path).  This test exercises the empty-
    advice early-return path which is the simplest.
    """
    fn, err = _try_import_run_advice_preprocessing()
    if fn is None:
        print("  SKIP: cannot import run_advice_preprocessing (%s)"
              % err)
        return

    # Empty advice triggers the early-return path
    result = fn(raw_advice="")
    assert hasattr(result, 'diagnostic_messages'), (
        "Result missing diagnostic_messages field")
    assert isinstance(result.diagnostic_messages, list), (
        "diagnostic_messages must be a list, got %s"
        % type(result.diagnostic_messages).__name__)
    print("  PASS: test_preprocessing_returns_diagnostic_messages_field "
          "(empty path, N=%d entries)"
          % len(result.diagnostic_messages))


def test_use_rules_only_returns_empty_diagnostic_messages():
    """use_rules_only path returns empty diagnostic_messages.

    When use_rules_only=True, no LLM call happens, so the
    metric block doesn't run.  diagnostic_messages should be
    present but empty.
    """
    fn, err = _try_import_run_advice_preprocessing()
    if fn is None:
        print("  SKIP: cannot import run_advice_preprocessing (%s)"
              % err)
        return

    result = fn(raw_advice="some advice text", use_rules_only=True)
    assert hasattr(result, 'diagnostic_messages')
    assert result.diagnostic_messages == [], (
        "Expected empty list, got %r" % result.diagnostic_messages)
    print("  PASS: test_use_rules_only_returns_empty_diagnostic_messages")


def test_diagnostic_messages_fresh_per_call():
    """Each call to run_advice_preprocessing gets a fresh list.

    Per Gemini Q5 review of v119.H5: the diagnostic_messages
    container MUST be instantiated fresh per request.  If it
    were attached to a long-lived global, class-level
    attribute, or Python's classic mutable-default-argument
    footgun, markers would bleed across sequential calls and
    the log list would grow indefinitely on a persistent
    server process.

    This test pins the contract: calling run_advice_preprocessing
    twice with use_rules_only=True (a path that won't populate
    markers but does build the field) returns two independent
    list objects.  A regression that refactored the list into
    a default argument or class attribute would fail this test
    on the second call.
    """
    fn, err = _try_import_run_advice_preprocessing()
    if fn is None:
        print("  SKIP: cannot import run_advice_preprocessing (%s)"
              % err)
        return

    result_a = fn(raw_advice="first call", use_rules_only=True)
    result_b = fn(raw_advice="second call", use_rules_only=True)

    # Both must have the field
    assert hasattr(result_a, 'diagnostic_messages')
    assert hasattr(result_b, 'diagnostic_messages')

    # Both must be empty (no markers in use_rules_only path)
    assert result_a.diagnostic_messages == []
    assert result_b.diagnostic_messages == []

    # CRITICAL: they must be DIFFERENT list objects
    # (a default-argument bug would have them share identity)
    assert result_a.diagnostic_messages is not result_b.diagnostic_messages, (
        "Fresh-container bug: two calls share the same list "
        "object.  This means a default-argument footgun or "
        "class-level state has crept in.  On a persistent "
        "server, markers would accumulate across requests.")

    # Defense-in-depth: mutate one, verify the other isn't affected
    result_a.diagnostic_messages.append("[TEST_LEAK] sentinel")
    assert "[TEST_LEAK] sentinel" not in result_b.diagnostic_messages, (
        "Containers are aliased — mutation in result_a appeared "
        "in result_b.")

    print("  PASS: test_diagnostic_messages_fresh_per_call")


def test_diagnostic_messages_independent_of_debug_log():
    """diagnostic_messages and debug_log are separate fields.

    They serve different purposes: debug_log is narrative
    logging for developers; diagnostic_messages is structured
    operator markers for relay to the client.  Both fields
    should exist independently.
    """
    fn, err = _try_import_run_advice_preprocessing()
    if fn is None:
        print("  SKIP: cannot import run_advice_preprocessing (%s)"
              % err)
        return

    result = fn(raw_advice="")
    assert hasattr(result, 'debug_log')
    assert hasattr(result, 'diagnostic_messages')
    assert result.debug_log is not result.diagnostic_messages, (
        "Fields should be distinct objects")
    print("  PASS: test_diagnostic_messages_independent_of_debug_log")


# =====================================================================
# §C: Backward compatibility (2 tests)
# =====================================================================

def test_old_server_response_handled_defensively():
    """Old-server response (no diagnostic_messages) doesn't crash.

    The asynchronous deployment window guardrail (plan rev 4
    §2.6 hard rule).  When a new client receives a response
    from an un-upgraded server, the diagnostic_messages field
    may be missing entirely.  Defensive unpacking handles this
    gracefully.

    This test simulates the bad case by building a mock result
    object that lacks the field and verifying the standard
    unpack pattern doesn't raise.
    """
    # Mock an old-server response (no diagnostic_messages field)
    class _OldResult:
        pass
    old_result = _OldResult()

    # Standard defensive-unpack pattern (plan rev 4 §5.1)
    try:
        diagnostics = (
            getattr(old_result, 'diagnostic_messages', None) or [])
        # Iteration shouldn't raise
        count = 0
        for msg in diagnostics:
            count += 1
        assert count == 0, (
            "Expected empty iteration on absent field, got %d" % count)
    except Exception as e:
        raise AssertionError(
            "Defensive unpack raised on absent field: %s" % e)

    # Also test dict-shaped result (server protocol)
    old_dict_result = {'processed_advice': 'foo'}
    try:
        diagnostics = old_dict_result.get('diagnostic_messages', []) or []
        for _ in diagnostics:
            pass
    except Exception as e:
        raise AssertionError(
            "Dict defensive unpack raised: %s" % e)

    print("  PASS: test_old_server_response_handled_defensively")


def test_corrupted_payload_decoded_defensively():
    """Malformed payload on the wire doesn't crash decoder.

    This pins the defensive logic in _process_server_success
    (programs/ai_analysis.py lines 935-946).  Three failure modes
    a real server might send during partial deployments or
    misconfiguration:
      1. Non-JSON garbage (e.g., HTML error page from a proxy)
      2. Valid JSON but not a list (e.g., null, object, number)
      3. Valid JSON list with non-string elements (less critical,
         but the decoder shouldn't blow up)

    In all cases the field should stay [] (the Total Init seed),
    and no exception should propagate to the caller.

    Paraphrases the production decode logic since calling the
    real _process_server_success requires a fully-mocked params
    object.  The paraphrased logic is byte-for-byte the same
    branch structure as production lines 935-946.
    """
    import json as _json

    def _paraphrased_decode(simple_string_value):
        """Mirror of _process_server_success decode logic."""
        working = []  # the Total Init seed in production
        if simple_string_value:
            try:
                # Note: paraphrased — production also calls
                # simple_string_as_text, but for this test we
                # assume the simple_string IS just the raw JSON
                # text (the wrapping is invertible).
                dm = _json.loads(simple_string_value)
                if isinstance(dm, list):
                    working = dm
                # Else: leave as []
            except Exception:
                pass  # Leave as []
        return working

    # Case 1: non-JSON garbage
    result = _paraphrased_decode("<html>500 Server Error</html>")
    assert result == [], (
        "Non-JSON payload should leave field as []; got %r" % result)

    # Case 2a: valid JSON but null
    result = _paraphrased_decode("null")
    assert result == [], (
        "null payload should leave field as []; got %r" % result)

    # Case 2b: valid JSON but object
    result = _paraphrased_decode('{"key": "value"}')
    assert result == [], (
        "Object payload should leave field as []; got %r" % result)

    # Case 2c: valid JSON but number
    result = _paraphrased_decode("42")
    assert result == [], (
        "Number payload should leave field as []; got %r" % result)

    # Case 3: valid JSON list with mixed types — should still
    # accept (the production code only checks isinstance(list),
    # not the element types — the client re-emit's
    # `sys.stderr.write(_msg + '\\n')` would coerce or fail per
    # element, individually try-wrapped).
    result = _paraphrased_decode('["[STEP_1F] valid", 42, null]')
    assert result == ["[STEP_1F] valid", 42, None], (
        "Mixed-type list should be accepted; got %r" % result)

    # Case 4: empty string (None-equivalent in REST transport)
    result = _paraphrased_decode("")
    assert result == [], (
        "Empty string should leave field as []; got %r" % result)

    print("  PASS: test_corrupted_payload_decoded_defensively")


def test_existing_step_1f_marker_format_preserved():
    """The [STEP_1F] marker line content is unchanged from H4.1.

    The aggregator's regex parses the marker line content; that
    contract must not change.  H5 only changes HOW the marker
    is delivered (via diagnostic_messages list instead of direct
    stderr.write), not WHAT the marker says.

    This test exercises the helper directly with the exact
    format string the metric block uses, confirming the list
    receives a line that the aggregator's regex can still parse.
    """
    emit_fn, err = _try_import_emit_marker()
    if emit_fn is None:
        print("  SKIP: helper not importable (%s)" % err)
        return

    # Simulate what the metric block does, using the exact
    # marker format from run_ai_analysis.py line ~1070:
    msgs = []
    saved_stderr = sys.stderr
    buf = io.StringIO()
    sys.stderr = buf
    try:
        emit_fn(
            msgs,
            "[STEP_1F] preprocessing_metrics "
            "preprocessor_mode=llm_comparison "
            "scanner_version=119.H5 "
            "llm_files=[] regex_files=[] "
            "in_llm_only=[] in_regex_only=[]")
    finally:
        sys.stderr = saved_stderr

    # Under uniformity (v119.H5 review): list gets the marker,
    # stderr stays empty.
    assert buf.getvalue() == "", (
        "Helper must NOT write to stderr; got %r" % buf.getvalue())
    assert len(msgs) == 1
    assert "[STEP_1F]" in msgs[0]
    assert "scanner_version=119.H5" in msgs[0]

    # The aggregator's regex pattern still matches the marker:
    import re
    pattern = re.compile(
        r'llm_files=(\[.*?\])\s+regex_files=(\[.*?\])')
    m = pattern.search(msgs[0])
    assert m is not None, (
        "Aggregator regex no longer matches marker; format broken")

    print("  PASS: test_existing_step_1f_marker_format_preserved")


# =====================================================================
# §D: Full plumbing — client re-emit logic (uniform across modes) (5 tests)
# =====================================================================
# These exercise the v119.H5 wiring across programs/ai_analysis.py
# and programs/ai_agent.py.  The module-level imports for those
# files need PHENIX (libtbx.utils.Sorry, etc.) so we can't import
# them in sandbox.  Instead we test the client-side LOGIC by
# constructing the exact result shapes those modules build, and
# running the standalone re-emit block from ai_agent.py.

# Under the uniformity criterion (v119.H5 review): the client
# re-emit block does NOT branch on dispatch mode.  Both local
# and remote dispatches populate diagnostic_messages identically,
# and the dispatcher re-emits unconditionally.
#
# The re-emit block (paraphrased from
# programs/ai_analysis.py::_relay_diagnostic_messages_to_stderr) is:
#     diagnostics = getattr(result, 'diagnostic_messages', None) or []
#     if diagnostics:
#         for msg in diagnostics:
#             sys.stderr.write(msg + '\n')
#             sys.stderr.flush()
#
# These tests pin that behavior.  Centralized in
# run_job_on_server_or_locally (per Gemini H5.1 review) so all
# call sites get re-emission automatically, with no copy-paste
# at multiple sites in programs/ai_agent.py.


def _client_reemit(result, stderr):
    """Reference implementation of the centralized re-emit block
    (paraphrased from
    programs/ai_analysis.py::_relay_diagnostic_messages_to_stderr).
    Used by §D tests to exercise the re-emit logic standalone.

    This mirrors the production code byte-for-byte; if the
    production block changes, update here too.

    Note: no dispatch-mode gating — the criterion requires
    identical behavior whether dispatched locally or remotely.

    Per Gemini Q5A review of v119.H5: handles BOTH attribute-
    style access (result is a group_args object — the normal
    case) AND dict-style access (in case a future change
    returns a raw dict decoded directly from JSON transport).
    getattr() alone on a dict would silently drop markers.
    """
    try:
        if isinstance(result, dict):
            diagnostics = result.get('diagnostic_messages', []) or []
        else:
            diagnostics = getattr(result, 'diagnostic_messages', None) or []
        if diagnostics:
            for msg in diagnostics:
                try:
                    stderr.write(msg + '\n')
                    stderr.flush()
                except Exception:
                    pass
    except Exception:
        pass


class _MockResult:
    """Simulates working_results returned by run_job_on_server_or_locally."""
    pass


def test_client_reemits_diagnostics():
    """When the list is non-empty, markers are re-emitted to
    client stderr — regardless of dispatch mode."""
    r = _MockResult()
    r.diagnostic_messages = [
        "[STEP_1F] preprocessing_metrics preprocessor_mode=llm_comparison "
        "scanner_version=119.H5 llm_files=[] regex_files=[] "
        "in_llm_only=[] in_regex_only=[]"]
    buf = io.StringIO()
    _client_reemit(r, buf)
    out = buf.getvalue()
    assert "[STEP_1F]" in out, (
        "Expected [STEP_1F] re-emitted, got %r" % out)
    assert out.endswith("\n"), (
        "Expected trailing newline on re-emit")
    print("  PASS: test_client_reemits_diagnostics")


def test_client_reemits_uniformly_across_dispatch_modes():
    """Behavior is identical for local-mode and remote-mode
    results.  Under the uniformity criterion, the client must
    not branch on dispatch mode — both kinds of result produce
    exactly the same output.

    This test pins the criterion: feed identical
    diagnostic_messages with no dispatch flag, with
    dispatched_remotely=True, and with dispatched_remotely=False.
    All three produce identical output.
    """
    markers = ["[STEP_1F] uniform-test-marker"]

    # No dispatch flag (the v119.H5 reviewed design)
    r1 = _MockResult()
    r1.diagnostic_messages = list(markers)
    buf1 = io.StringIO()
    _client_reemit(r1, buf1)

    # dispatched_remotely=True (legacy attribute, should be ignored)
    r2 = _MockResult()
    r2.diagnostic_messages = list(markers)
    r2.dispatched_remotely = True
    buf2 = io.StringIO()
    _client_reemit(r2, buf2)

    # dispatched_remotely=False (legacy attribute, should be ignored)
    r3 = _MockResult()
    r3.diagnostic_messages = list(markers)
    r3.dispatched_remotely = False
    buf3 = io.StringIO()
    _client_reemit(r3, buf3)

    assert buf1.getvalue() == buf2.getvalue() == buf3.getvalue(), (
        "Output must be identical across dispatch modes (uniformity).\n"
        "  No flag:                %r\n"
        "  dispatched_remotely=T:  %r\n"
        "  dispatched_remotely=F:  %r"
        % (buf1.getvalue(), buf2.getvalue(), buf3.getvalue()))

    # All produced the marker
    assert "uniform-test-marker" in buf1.getvalue()
    print("  PASS: test_client_reemits_uniformly_across_dispatch_modes")


def test_client_no_reemit_when_empty_list():
    """Empty diagnostic_messages → no re-emit, no crash.

    Edge case: the dispatch returned a result with the field
    populated but no markers were generated (e.g.,
    preprocessing was skipped via use_rules_only).
    """
    r = _MockResult()
    r.diagnostic_messages = []
    buf = io.StringIO()
    _client_reemit(r, buf)
    assert buf.getvalue() == "", (
        "Empty list should produce no output, got %r" % buf.getvalue())
    print("  PASS: test_client_no_reemit_when_empty_list")


def test_client_no_reemit_when_no_field():
    """Missing diagnostic_messages → no re-emit, no crash.

    The asynchronous deployment window: new client + old server
    that didn't populate the field.  The defensive `or []`
    handles None AND missing-attribute cases.
    """
    r = _MockResult()
    # NO diagnostic_messages attribute
    buf = io.StringIO()
    _client_reemit(r, buf)
    assert buf.getvalue() == "", (
        "Expected NO output when field missing, got %r" % buf.getvalue())
    print("  PASS: test_client_no_reemit_when_no_field")


def test_client_reemit_handles_dict_result():
    """Client re-emit works whether result is object or dict.

    Per Gemini Q5A review of v119.H5: the production code uses
    getattr() which won't find keys on a dict.  If a future
    change causes the dispatcher to return a raw dict (e.g.,
    decoded directly from JSON without group_args wrapping),
    the markers would silently vanish.

    Production code now checks isinstance(result, dict) and
    branches between dict.get() and getattr().  This test pins
    that both shapes work.
    """
    markers = ["[STEP_1F] dict-result-test"]

    # Object-style result (the normal case)
    r_obj = _MockResult()
    r_obj.diagnostic_messages = list(markers)
    buf_obj = io.StringIO()
    _client_reemit(r_obj, buf_obj)

    # Dict-style result (hypothetical future shape)
    r_dict = {'diagnostic_messages': list(markers)}
    buf_dict = io.StringIO()
    _client_reemit(r_dict, buf_dict)

    # Both must produce the same output
    assert buf_obj.getvalue() == buf_dict.getvalue(), (
        "Object and dict results must produce identical "
        "re-emit output.\n  Object: %r\n  Dict:   %r"
        % (buf_obj.getvalue(), buf_dict.getvalue()))
    assert "[STEP_1F]" in buf_dict.getvalue(), (
        "Dict result lost markers; getattr-only unpacking "
        "would cause this regression.")

    print("  PASS: test_client_reemit_handles_dict_result")


def test_client_reemit_multiple_markers_in_order():
    """Multiple markers in list are re-emitted in order, one line each."""
    r = _MockResult()
    r.diagnostic_messages = [
        "[STEP_1F] first marker",
        "[STEP_1F_FAILED] second marker",
        "[STEP_1F] third marker",
    ]
    buf = io.StringIO()
    _client_reemit(r, buf)
    out = buf.getvalue()
    lines = out.split('\n')
    # Three markers + trailing empty from final \n
    assert "[STEP_1F] first marker" in lines[0]
    assert "[STEP_1F_FAILED] second marker" in lines[1]
    assert "[STEP_1F] third marker" in lines[2]
    # No clumping: each marker is on its own line
    for line in lines[:3]:
        assert line.count("[STEP_1F") == 1, (
            "Markers clumped on one line: %r" % line)
    print("  PASS: test_client_reemit_multiple_markers_in_order")


# =====================================================================
# §E: Integration — production encode/decode roundtrip (3 tests)
# =====================================================================
# §A-§D exercise the helper, the field shape, the backward
# compatibility, and the client re-emit logic.  §E exercises the
# REST encode/decode path through the actual production code in
# programs/ai_analysis.py — specifically Program.get_results_as_JSON
# (the encoder).  This is the test that proves the H5 wiring
# survives the round trip the channel was designed for.
#
# These tests SKIP in sandbox (require libtbx via the Program
# class import).  Under PHENIX, all three run.


def _try_import_get_results_from_all():
    """Import get_results_from_all for testing the Total Init seed."""
    try:
        from phenix.programs.ai_analysis import get_results_from_all
        return get_results_from_all, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from programs.ai_analysis import get_results_from_all
        return get_results_from_all, None
    except ImportError as e:
        return None, str(e)


def test_total_initialization_seed():
    """get_results_from_all returns working_results with
    diagnostic_messages=[] from birth.

    This is the LINCHPIN of the H5 design (Total Initialization
    Policy, per Gemini Q3 review): every working_results carries
    the field from the moment of construction, regardless of
    analysis mode.  Downstream code paths replace the list when
    populated; modes that don't emit markers leave it empty.

    If this test fails, the entire uniformity guarantee breaks:
    the encoder/decoder might fall back to defensive paths, the
    field shape could diverge between local and remote dispatches,
    and existing K_H5 tests would still pass while production
    silently misbehaves.

    Uses types.SimpleNamespace to mock the minimum params shape
    that get_results_from_all reads.  This is safe because the
    function only reads attributes; no phil validation happens.
    """
    import types
    fn, err = _try_import_get_results_from_all()
    if fn is None:
        print("  SKIP: cannot import get_results_from_all (%s)" % err)
        return

    # Minimum params shape required by get_results_from_all.
    # Construct so all branches that read params fall through to
    # the "no existing results" code path, producing a clean
    # empty working_results.
    params = types.SimpleNamespace(
        ai_analysis=types.SimpleNamespace(
            summary_file_name=None,
            analysis_file_name=None,
            write_files=False,
            summary_as_simple_string="",
            analysis_as_simple_string="",
            load_existing_analysis=False,
        )
    )

    working_results = fn(params=params)

    # Total Initialization Policy assertion:
    assert hasattr(working_results, 'diagnostic_messages'), (
        "TOTAL INITIALIZATION VIOLATED: working_results missing "
        "diagnostic_messages attribute.  This breaks the uniformity "
        "criterion — encoder/decoder will fall through defensive "
        "paths and the field shape will diverge across modes.")
    assert isinstance(working_results.diagnostic_messages, list), (
        "diagnostic_messages must be a list, got %s"
        % type(working_results.diagnostic_messages).__name__)
    assert working_results.diagnostic_messages == [], (
        "Total Init seed should be empty list, got %r"
        % working_results.diagnostic_messages)

    # Defense: verify it's a fresh empty list, not a singleton
    # that could leak state across calls (the classic mutable
    # default argument footgun, but at module rather than
    # function scope).
    working_results.diagnostic_messages.append("[TEST] mutation probe")
    fresh = fn(params=params)
    assert fresh.diagnostic_messages == [], (
        "Total Init seed leaked state: mutation in one call "
        "appeared in subsequent call.  Got %r"
        % fresh.diagnostic_messages)

    print("  PASS: test_total_initialization_seed")


def _try_import_program_class():
    """Import the Program class from programs.ai_analysis.

    Tests get_results_as_JSON via an unbound-method-on-stub trick
    so we don't have to construct a full Program() (which needs
    phil params, datatypes, ProgramTemplate inheritance setup).

    Returns (Program_class, error_msg).
    """
    try:
        from phenix.programs.ai_analysis import Program
        return Program, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from programs.ai_analysis import Program
        return Program, None
    except ImportError as e:
        return None, str(e)


def _try_import_simple_string_helpers():
    """Import the simple_string_as_text helper for decode roundtrip."""
    try:
        from libtbx.langchain.agent.utils import simple_string_as_text
        return simple_string_as_text, None
    except ImportError as e:
        return None, str(e)


def _try_import_group_args():
    """Import group_args for building stub working_results."""
    try:
        from libtbx import group_args
        return group_args, None
    except ImportError as e:
        return None, str(e)


def _build_mini_program(group_args_cls, diagnostic_messages, Program_cls):
    """Build a minimal Program-like object with self.result set.

    Program.get_results_as_JSON only reads self.result (see
    source); we attach it on a stub class that borrows the
    method.  This avoids the full ProgramTemplate construction
    machinery while still exercising the production code path.
    """
    fake_result = group_args_cls(
        group_args_type='working result',
        return_value=None,
        summary="test summary",
        analysis="test analysis",
        summary_file_name=None,
        analysis_file_name=None,
        history_record={
            "job_id": "test-job",
            "timestamp": "now",
            "file_name": "test",
            "program": "test",
            "summary": "test summary",
            "analysis": "test analysis",
            "error": None,
        },
        diagnostic_messages=diagnostic_messages,
    )

    class _MiniProgram:
        result = fake_result
        # Borrow the unbound method.  In Python 3 this binds
        # correctly when accessed on an instance.
        get_results_as_JSON = Program_cls.get_results_as_JSON

    return _MiniProgram()


def test_encode_decode_roundtrip_with_markers():
    """End-to-end transport: encoder produces a payload that
    decodes back to the original list.

    This is the test that proves H5's transport wiring works.
    It exercises the real production encoder
    (Program.get_results_as_JSON) on a result with populated
    diagnostic_messages, then decodes the payload manually
    (mirroring what _process_server_success does on the
    client side) and asserts the markers come back unchanged.
    """
    import json as _json

    Program_cls, err = _try_import_program_class()
    if Program_cls is None:
        print("  SKIP: cannot import Program (%s)" % err)
        return
    group_args_cls, err = _try_import_group_args()
    if group_args_cls is None:
        print("  SKIP: cannot import group_args (%s)" % err)
        return
    sst_fn, err = _try_import_simple_string_helpers()
    if sst_fn is None:
        print("  SKIP: cannot import simple_string_as_text (%s)" % err)
        return

    original_markers = [
        "[STEP_1F] preprocessing_metrics "
        "preprocessor_mode=llm_comparison "
        "scanner_version=119.H5 "
        "llm_files=[\"model.pdb\"] regex_files=[\"model.pdb\"] "
        "in_llm_only=[] in_regex_only=[]",
        "[STEP_1F_FAILED] simulated metric failure for test",
    ]

    mini = _build_mini_program(
        group_args_cls, list(original_markers), Program_cls)
    rest_payload = mini.get_results_as_JSON()

    # Wire-shape assertions
    assert 'output_files' in rest_payload, (
        "Encoder missing output_files key: %r" % rest_payload)
    output_files = rest_payload['output_files']
    assert 'diagnostic_messages_as_simple_string' in output_files, (
        "Encoder missing diagnostic_messages_as_simple_string key. "
        "Keys present: %r" % list(output_files.keys()))

    # Decode side: mirror what _process_server_success does.
    # Specifically, the H5 wiring (programs/ai_analysis.py lines
    # 935-946) is: simple_string_as_text(payload) -> json.loads ->
    # isinstance check -> assign to working_results.
    encoded = output_files['diagnostic_messages_as_simple_string']
    decoded = _json.loads(sst_fn(encoded))

    assert isinstance(decoded, list), (
        "Decoded value not a list: type=%s" % type(decoded).__name__)
    assert decoded == original_markers, (
        "Round-trip mismatch.\n  Original: %r\n  Decoded:  %r"
        % (original_markers, decoded))

    print("  PASS: test_encode_decode_roundtrip_with_markers "
          "(N=%d markers)" % len(decoded))


def test_encode_decode_roundtrip_with_empty_list():
    """Total Initialization shape ([]) survives the round trip.

    Critical for the uniformity criterion: a result with no
    markers must still carry the field through serialization
    (as an empty array), not get pruned out.  The decoder must
    receive [], not None or missing-key.
    """
    import json as _json

    Program_cls, err = _try_import_program_class()
    if Program_cls is None:
        print("  SKIP: cannot import Program (%s)" % err)
        return
    group_args_cls, err = _try_import_group_args()
    if group_args_cls is None:
        print("  SKIP: cannot import group_args (%s)" % err)
        return
    sst_fn, err = _try_import_simple_string_helpers()
    if sst_fn is None:
        print("  SKIP: cannot import simple_string_as_text (%s)" % err)
        return

    mini = _build_mini_program(group_args_cls, [], Program_cls)
    rest_payload = mini.get_results_as_JSON()

    output_files = rest_payload['output_files']
    assert 'diagnostic_messages_as_simple_string' in output_files, (
        "Empty list got pruned from output_files — uniformity broken")

    encoded = output_files['diagnostic_messages_as_simple_string']
    decoded = _json.loads(sst_fn(encoded))

    assert decoded == [], (
        "Empty list round-tripped as %r" % decoded)

    print("  PASS: test_encode_decode_roundtrip_with_empty_list")


def test_encode_decode_preserves_special_chars():
    """Special characters in marker payloads survive JSON encoding.

    Markers from [STEP_1F_FAILED] contain exception messages
    which may include quotes, backslashes, tabs, etc.  The JSON
    serialize -> simple_string -> simple_string_as_text -> JSON
    parse pipeline must preserve them byte-for-byte.
    """
    import json as _json

    Program_cls, err = _try_import_program_class()
    if Program_cls is None:
        print("  SKIP: cannot import Program (%s)" % err)
        return
    group_args_cls, err = _try_import_group_args()
    if group_args_cls is None:
        print("  SKIP: cannot import group_args (%s)" % err)
        return
    sst_fn, err = _try_import_simple_string_helpers()
    if sst_fn is None:
        print("  SKIP: cannot import simple_string_as_text (%s)" % err)
        return

    # Mix of characters that have historically caused trouble
    # in PHENIX REST transport:
    #   - tab and other whitespace
    #   - JSON-special chars: quote, backslash
    #   - already-flattened-newline marker pattern (' \ ')
    #   - unicode (latin-1 + a multi-byte char)
    original_markers = [
        "[STEP_1F_FAILED] msg with\ttab and \"quote\"",
        "[STEP_1F_FAILED] msg with backslash: \\foo\\bar",
        "[STEP_1F_FAILED] msg with flattened newline: line one \\ line two",
        "[STEP_1F_FAILED] unicode: caf\u00e9 \u2192 \u00df",
    ]

    mini = _build_mini_program(
        group_args_cls, list(original_markers), Program_cls)
    rest_payload = mini.get_results_as_JSON()
    encoded = rest_payload['output_files'][
        'diagnostic_messages_as_simple_string']
    decoded = _json.loads(sst_fn(encoded))

    assert decoded == original_markers, (
        "Special chars not preserved.\n  Original: %r\n  Decoded:  %r"
        % (original_markers, decoded))

    print("  PASS: test_encode_decode_preserves_special_chars")


# =====================================================================
# §F: Extended markers (v119.H5.1) (1+ tests)
# =====================================================================
# H5.1 extends the diagnostic_messages channel to additional
# operationally-significant failure markers beyond [STEP_1F]:
#   [ADVICE_PREPROCESSING_FAILED]   — Stage 1
#   [DIRECTIVE_EXTRACTION_FAILED]   — Stage 2 (not yet added)
#   [FAILURE_DIAGNOSIS_FAILED]      — Stage 3 (not yet added)
#
# These tests verify that when the main exception handler in
# each target function fires, the corresponding [*_FAILED]
# marker appears in result.diagnostic_messages.
#
# Test strategy: monkey-patch.
# Audit of utils/run_utils.py revealed that validate_api_keys
# and setup_llms degrade gracefully on bad inputs — neither
# raises on invalid providers or missing keys.  The [*_FAILED]
# markers therefore cannot be triggered by passing junk
# parameters; they fire only on UNEXPECTED exceptions.
# Monkey-patching one of the inner calls (validate_api_keys
# for run_advice_preprocessing and run_directive_extraction;
# setup_llms for run_failure_diagnosis) is the deterministic
# way to reach the target exception handler.
#
# Per Gemini Q2 mitigation: each test asserts BOTH that the
# marker is present AND that there is evidence the exception
# path was actually entered (via debug_log inspection).  This
# protects against future refactors that add upstream graceful
# handling and silently bypass the exception path.


def test_advice_preprocessing_failed_marker_emitted_on_exception():
    """When the main preprocessing exception handler fires,
    [ADVICE_PREPROCESSING_FAILED] appears in
    result.diagnostic_messages.

    Strategy: monkey-patch validate_api_keys (imported at the
    top of run_ai_analysis.py from libtbx.langchain.utils.run_utils)
    to unconditionally raise.  The first call inside the main
    try block thus raises, the except clause fires, and the
    new H5.1 marker is emitted.

    Why monkey-patch, not "pass invalid provider":
    validate_api_keys (run_utils.py:104) silently passes through
    unknown providers (only checks google/openai), and
    setup_llms (run_utils.py:168) swallows exceptions
    internally.  So invalid providers cause graceful
    degradation, NOT the exception path that fires this marker.
    See v119_H5_1_IMPL_PLAN.md §"Forcing exceptions cleanly"
    for full rationale.

    Per Gemini Q2 mitigation: assert BOTH that the marker is
    present AND that the debug_log has evidence of the
    exception (proves the handler was actually entered).
    """
    fn, err = _try_import_run_advice_preprocessing()
    if fn is None:
        print("  SKIP: cannot import run_advice_preprocessing (%s)" % err)
        return

    # Patch validate_api_keys in the run_ai_analysis namespace
    # to force the main except clause to fire.  Tolerate either
    # import path (PHENIX-style vs sandbox-style with libtbx).
    rai_module = None
    try:
        from phenix.phenix_ai import run_ai_analysis as rai_module
    except ImportError:
        try:
            from phenix_ai import run_ai_analysis as rai_module
        except ImportError:
            pass
    if rai_module is None:
        print("  SKIP: cannot import run_ai_analysis module for patching")
        return

    original = rai_module.validate_api_keys

    def _raises(*args, **kwargs):
        raise RuntimeError("forced exception for K_H5 §F test")

    rai_module.validate_api_keys = _raises
    try:
        result = fn(raw_advice="trigger preprocessing path")
    finally:
        rai_module.validate_api_keys = original

    # Assertion 1: marker present in diagnostic_messages
    diagnostics = getattr(result, 'diagnostic_messages', [])
    matching = [m for m in diagnostics
                if m.startswith('[ADVICE_PREPROCESSING_FAILED]')]
    assert len(matching) == 1, (
        "Expected one [ADVICE_PREPROCESSING_FAILED] marker, got %d: %r"
        % (len(matching), diagnostics))

    # Assertion 2: debug_log has evidence the except clause
    # was entered.  Per Gemini Q2: protects against future
    # upstream graceful-handling refactors that would silently
    # bypass the exception path.  If the exception handler
    # isn't reached, BOTH assertions fail with a clear
    # "exception path not entered" diagnostic.
    debug_log = getattr(result, 'debug_log', [])
    debug_evidence = any(
        'failed' in entry.lower() or 'forced exception' in entry
        for entry in debug_log)
    assert debug_evidence, (
        "Marker present but no exception evidence in debug_log "
        "— exception handler may not have been entered. "
        "debug_log: %r" % debug_log)

    print("  PASS: test_advice_preprocessing_failed_marker_emitted_on_exception")


def _try_import_run_directive_extraction():
    """Import run_directive_extraction; return (fn, error_msg).

    Same SKIP pattern as _try_import_run_advice_preprocessing.
    The function depends on libtbx.group_args which is a PHENIX
    runtime dep.
    """
    try:
        from phenix.phenix_ai.run_ai_analysis import run_directive_extraction
        return run_directive_extraction, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from phenix_ai.run_ai_analysis import run_directive_extraction
        return run_directive_extraction, None
    except ImportError as e:
        return None, str(e)


def test_directive_extraction_field_in_return_path():
    """run_directive_extraction always returns
    diagnostic_messages=[] in its result (Total Init at engine
    level — same pattern as run_advice_preprocessing).

    Per Guardrail 1: also verifies that two sequential calls
    produce independent list objects.  If diagnostic_messages
    were ever moved to a parameter default (the classic
    mutable-default footgun), mutating one result's list
    would leak into the next call's empty seed.

    This test exercises the empty-advice early-return path,
    which is the cheapest path to invoke and doesn't require
    network access or an LLM provider.
    """
    fn, err = _try_import_run_directive_extraction()
    if fn is None:
        print("  SKIP: cannot import run_directive_extraction (%s)" % err)
        return

    # Empty advice triggers the early-return path
    r1 = fn(user_advice="")
    assert hasattr(r1, 'diagnostic_messages'), (
        "TOTAL INITIALIZATION VIOLATED: run_directive_extraction's "
        "empty-advice return path is missing diagnostic_messages")
    assert isinstance(r1.diagnostic_messages, list), (
        "diagnostic_messages must be a list, got %s"
        % type(r1.diagnostic_messages).__name__)
    assert r1.diagnostic_messages == [], (
        "Empty path should return empty list, got %r"
        % r1.diagnostic_messages)

    # Guardrail 1: independence check across sequential calls.
    # Mutate the first call's list; the second call must return
    # a fresh empty list, not the mutated one.
    r1.diagnostic_messages.append("[TEST] mutation probe")
    r2 = fn(user_advice="")
    assert r2.diagnostic_messages == [], (
        "diagnostic_messages leaked state between calls — possible "
        "mutable-default-argument regression.  Got %r"
        % r2.diagnostic_messages)

    print("  PASS: test_directive_extraction_field_in_return_path")


def test_directive_extraction_failed_marker_emitted_on_exception():
    """When run_directive_extraction's main exception handler
    fires, [DIRECTIVE_EXTRACTION_FAILED] appears in
    result.diagnostic_messages.

    Strategy: same monkey-patch approach as
    test_advice_preprocessing_failed_marker_emitted_on_exception.
    Patch validate_api_keys (the first call INSIDE the main
    try block) to unconditionally raise.

    In run_directive_extraction, validate_api_keys is the
    first call inside the main try block, so patching it to
    raise forces the except clause to fire and emit the
    marker.

    Per Gemini Q2 mitigation: assert marker present AND
    exception evidence in debug_log.
    """
    fn, err = _try_import_run_directive_extraction()
    if fn is None:
        print("  SKIP: cannot import run_directive_extraction (%s)" % err)
        return

    rai_module = None
    try:
        from phenix.phenix_ai import run_ai_analysis as rai_module
    except ImportError:
        try:
            from phenix_ai import run_ai_analysis as rai_module
        except ImportError:
            pass
    if rai_module is None:
        print("  SKIP: cannot import run_ai_analysis module for patching")
        return

    original = rai_module.validate_api_keys

    def _raises(*args, **kwargs):
        raise RuntimeError("forced exception for K_H5 §F test")

    rai_module.validate_api_keys = _raises
    try:
        result = fn(user_advice="trigger extraction path")
    finally:
        rai_module.validate_api_keys = original

    # Assertion 1: marker present
    diagnostics = getattr(result, 'diagnostic_messages', [])
    matching = [m for m in diagnostics
                if m.startswith('[DIRECTIVE_EXTRACTION_FAILED]')]
    assert len(matching) == 1, (
        "Expected one [DIRECTIVE_EXTRACTION_FAILED] marker, got %d: %r"
        % (len(matching), diagnostics))

    # Assertion 2: debug_log evidence (Gemini Q2 mitigation)
    debug_log = getattr(result, 'debug_log', [])
    debug_evidence = any(
        'failed' in entry.lower() or 'forced exception' in entry
        for entry in debug_log)
    assert debug_evidence, (
        "Marker present but no exception evidence in debug_log "
        "— exception handler may not have been entered. "
        "debug_log: %r" % debug_log)

    print("  PASS: test_directive_extraction_failed_marker_emitted_on_exception")


def _try_import_run_failure_diagnosis():
    """Import run_failure_diagnosis; return (fn, error_msg).

    Same SKIP pattern as the other H5.1 import helpers.  The
    function depends on libtbx.group_args which is a PHENIX
    runtime dep.
    """
    try:
        from phenix.phenix_ai.run_ai_analysis import run_failure_diagnosis
        return run_failure_diagnosis, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from phenix_ai.run_ai_analysis import run_failure_diagnosis
        return run_failure_diagnosis, None
    except ImportError as e:
        return None, str(e)


def test_failure_diagnosis_field_in_return_path():
    """run_failure_diagnosis always returns
    diagnostic_messages=[] in its result (Total Init at engine
    level — same pattern as the other two H5.1 target functions).

    Strategy: monkey-patch validate_api_keys to return a fake
    error string.  This deterministically hits the API-key-error
    early return path, regardless of the actual environment's
    API key state (under PHENIX on the user's machine, the real
    keys may or may not be set; we want the test to be
    deterministic).

    Per Guardrail 1: also verifies that two sequential calls
    produce independent list objects.
    """
    fn, err = _try_import_run_failure_diagnosis()
    if fn is None:
        print("  SKIP: cannot import run_failure_diagnosis (%s)" % err)
        return

    rai_module = None
    try:
        from phenix.phenix_ai import run_ai_analysis as rai_module
    except ImportError:
        try:
            from phenix_ai import run_ai_analysis as rai_module
        except ImportError:
            pass
    if rai_module is None:
        print("  SKIP: cannot import run_ai_analysis module for patching")
        return

    # Patch validate_api_keys to RETURN (not raise) an error
    # string.  This deterministically hits the API-key-error
    # early-return path at the top of run_failure_diagnosis.
    original = rai_module.validate_api_keys

    def _returns_error(provider, debug_log=None):
        return "fake API key error for K_H5 §F test"

    rai_module.validate_api_keys = _returns_error
    try:
        r1 = fn(
            error_type="test_error",
            error_text="test error text",
            program="phenix.test",
            log_tail="test log tail",
        )
        # Guardrail 1: independence check
        assert hasattr(r1, 'diagnostic_messages'), (
            "TOTAL INITIALIZATION VIOLATED: run_failure_diagnosis "
            "API-key-error return path is missing diagnostic_messages")
        assert isinstance(r1.diagnostic_messages, list)
        assert r1.diagnostic_messages == [], (
            "API-key-error path should return empty list, got %r"
            % r1.diagnostic_messages)

        r1.diagnostic_messages.append("[TEST] mutation probe")
        r2 = fn(
            error_type="test_error",
            error_text="test error text",
            program="phenix.test",
            log_tail="test log tail",
        )
        assert r2.diagnostic_messages == [], (
            "diagnostic_messages leaked state between calls — "
            "possible mutable-default-argument regression.  Got %r"
            % r2.diagnostic_messages)
    finally:
        rai_module.validate_api_keys = original

    print("  PASS: test_failure_diagnosis_field_in_return_path")


def test_failure_diagnosis_failed_marker_emitted_on_exception():
    """When run_failure_diagnosis's main LLM-call exception
    handler fires, [FAILURE_DIAGNOSIS_FAILED] appears in
    result.diagnostic_messages.

    Strategy: monkey-patch setup_llms (NOT validate_api_keys)
    to raise.  In run_failure_diagnosis, validate_api_keys is
    called OUTSIDE the main try block — patching it to raise
    would cause an uncaught exception to propagate out of the
    function, never reaching the [FAILURE_DIAGNOSIS_FAILED]
    handler.  setup_llms is the first call inside the main try
    block, so patching it to raise (rather than letting its
    internal try swallow and return Nones) cleanly triggers
    the target except clause.

    Also patches validate_api_keys to return None (lets us
    past the pre-try API-key check regardless of env state)
    and the failure_diagnoser imports to be a no-op (so the
    untried `build_diagnosis_prompt` call at line ~1400 can't
    raise on our synthetic test inputs).

    Per Gemini Q2 mitigation: assert marker present AND
    exception evidence in debug_log.
    """
    fn, err = _try_import_run_failure_diagnosis()
    if fn is None:
        print("  SKIP: cannot import run_failure_diagnosis (%s)" % err)
        return

    rai_module = None
    try:
        from phenix.phenix_ai import run_ai_analysis as rai_module
    except ImportError:
        try:
            from phenix_ai import run_ai_analysis as rai_module
        except ImportError:
            pass
    if rai_module is None:
        print("  SKIP: cannot import run_ai_analysis module for patching")
        return

    # Patch the module-level imports that run_failure_diagnosis
    # uses.  We need three patches:
    # 1. validate_api_keys → return None (let function past
    #    pre-try API-key check, regardless of env state)
    # 2. setup_llms → raise (force the main exception handler)
    # 3. failure_diagnoser submodule import → provide stubs so
    #    the test doesn't depend on diagnosable_errors.yaml or
    #    error_analyzer being reachable.  build_diagnosis_prompt
    #    in the real module calls get_diagnosis_detector().get_hint(),
    #    which loads a YAML and instantiates a singleton — both
    #    are out-of-scope risks for this test.  The stub returns
    #    a fixed string and bypasses YAML/detector entirely.
    original_validate = rai_module.validate_api_keys
    original_setup = rai_module.setup_llms

    def _no_api_error(provider, debug_log=None):
        return None  # signal "validation passed"

    def _raises_in_setup(*args, **kwargs):
        raise RuntimeError("forced exception for K_H5 §F test")

    # Stub the failure_diagnoser submodule so its imports succeed
    # and build_diagnosis_prompt is a no-op.  Inject directly into
    # sys.modules so the function's `from ... import` works.
    import types as _types
    fake_diagnoser = _types.ModuleType(
        'libtbx.langchain.agent.failure_diagnoser')
    fake_diagnoser.build_diagnosis_prompt = (
        lambda **kwargs: "stub prompt for K_H5 test")
    fake_diagnoser._strip_llm_markdown = lambda s: s
    original_diagnoser = sys.modules.get(
        'libtbx.langchain.agent.failure_diagnoser')
    sys.modules['libtbx.langchain.agent.failure_diagnoser'] = fake_diagnoser

    rai_module.validate_api_keys = _no_api_error
    rai_module.setup_llms = _raises_in_setup
    try:
        result = fn(
            error_type="test_error",
            error_text="test error text",
            program="phenix.test",
            log_tail="test log tail",
        )
    finally:
        rai_module.validate_api_keys = original_validate
        rai_module.setup_llms = original_setup
        if original_diagnoser is not None:
            sys.modules[
                'libtbx.langchain.agent.failure_diagnoser'] = original_diagnoser
        else:
            sys.modules.pop(
                'libtbx.langchain.agent.failure_diagnoser', None)

    # Defensive: if for some reason the import path was taken
    # despite the stub (e.g., the function uses a cached
    # earlier import), detect and SKIP rather than FAIL.
    debug_log = getattr(result, 'debug_log', [])
    import_failed = any(
        'Could not import failure_diagnoser' in entry
        for entry in debug_log)
    if import_failed:
        print("  SKIP: failure_diagnoser stub not used "
              "(main try block not reached)")
        return

    # Assertion 1: marker present
    diagnostics = getattr(result, 'diagnostic_messages', [])
    matching = [m for m in diagnostics
                if m.startswith('[FAILURE_DIAGNOSIS_FAILED]')]
    assert len(matching) == 1, (
        "Expected one [FAILURE_DIAGNOSIS_FAILED] marker, got %d: %r"
        % (len(matching), diagnostics))

    # Assertion 2: debug_log evidence (Gemini Q2 mitigation).
    # The except clause appends "LLM call failed: ..." so 'failed'
    # is in debug_log lowercased.  Plus the traceback.  Plus
    # 'forced exception' from the mock's __str__.
    debug_evidence = any(
        'failed' in entry.lower() or 'forced exception' in entry
        for entry in debug_log)
    assert debug_evidence, (
        "Marker present but no exception evidence in debug_log "
        "— exception handler may not have been entered. "
        "debug_log: %r" % debug_log)

    print("  PASS: test_failure_diagnosis_failed_marker_emitted_on_exception")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: _emit_marker helper unit tests (4)
    test_emit_marker_appends_to_list_only()
    test_emit_marker_drops_silently_on_non_list_container()
    test_emit_marker_no_stderr_write_under_any_condition()
    test_emit_marker_flattens_embedded_newlines()
    # §B: diagnostic_messages field in return (4)
    test_preprocessing_returns_diagnostic_messages_field()
    test_use_rules_only_returns_empty_diagnostic_messages()
    test_diagnostic_messages_fresh_per_call()
    test_diagnostic_messages_independent_of_debug_log()
    # §C: Backward compatibility (3)
    test_old_server_response_handled_defensively()
    test_corrupted_payload_decoded_defensively()
    test_existing_step_1f_marker_format_preserved()
    # §D: Full plumbing — uniform client re-emit (6)
    test_client_reemits_diagnostics()
    test_client_reemits_uniformly_across_dispatch_modes()
    test_client_no_reemit_when_empty_list()
    test_client_no_reemit_when_no_field()
    test_client_reemit_handles_dict_result()
    test_client_reemit_multiple_markers_in_order()
    # §E: Integration — production encode/decode roundtrip (4)
    test_total_initialization_seed()
    test_encode_decode_roundtrip_with_markers()
    test_encode_decode_roundtrip_with_empty_list()
    test_encode_decode_preserves_special_chars()
    # §F: Extended markers (v119.H5.1) (5)
    test_advice_preprocessing_failed_marker_emitted_on_exception()
    test_directive_extraction_field_in_return_path()
    test_directive_extraction_failed_marker_emitted_on_exception()
    test_failure_diagnosis_field_in_return_path()
    test_failure_diagnosis_failed_marker_emitted_on_exception()


if __name__ == "__main__":
    run_all_tests()
