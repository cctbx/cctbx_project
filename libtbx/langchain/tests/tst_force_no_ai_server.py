"""
FORCE_NO_AI_SERVER (v120 Task 1).

When the environment variable FORCE_NO_AI_SERVER == "1", phenix.ai_agent must
force communication.run_on_server = False, regardless of the PHIL default or an
explicit run_on_server=True.  Any other value (unset, "0", "true", "yes",
whitespace, ...) must leave run_on_server untouched.  This gives local-only
users the patch's behaviour without changing the shipped PHIL default.

Two layers:
  1. Behavioural — replicate the exact override block against lightweight stubs
     and assert the truth table (no Phenix imports; runs under plain python3).
  2. Source-scan — pin that programs/ai_agent.py actually contains the block,
     inside run(), after set_defaults() and before _handle_session_management(),
     so the behaviour can't silently drift away from the test.

Per the agent test guidelines: function-boundary search windows (not fixed
char offsets), file-exists guard, and AssertionError is never swallowed.
"""
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain
for _p in [_ROOT, os.path.join(_ROOT, "agent"), os.path.join(_ROOT, "programs")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

# programs/ai_agent.py and ai_analysis.py live in the phenix tree, not here.
# Locate them relative to plausible roots; source-scan tests skip (do not
# fail) if the file cannot be found, so behavioural tests still run anywhere.
# _ROOT is .../modules/cctbx_project/libtbx/langchain, so 3 ups reaches
# .../modules, then phenix/phenix/programs.  (Four ups overshoots modules.)
_AI_AGENT_CANDIDATES = [
    os.path.join(_ROOT, "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "programs", "ai_agent.py"),
]
_AI_ANALYSIS_CANDIDATES = [
    os.path.join(_ROOT, "programs", "ai_analysis.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_analysis.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "programs", "ai_analysis.py"),
]


def _find_ai_agent_source():
    for cand in _AI_AGENT_CANDIDATES:
        cand = os.path.abspath(cand)
        if os.path.isfile(cand):
            return cand
    env = os.environ.get("AI_AGENT_PY")
    if env and os.path.isfile(env):
        return env
    return None


def _find_ai_analysis_source():
    for cand in _AI_ANALYSIS_CANDIDATES:
        cand = os.path.abspath(cand)
        if os.path.isfile(cand):
            return cand
    env = os.environ.get("AI_ANALYSIS_PY")
    if env and os.path.isfile(env):
        return env
    return None


# ============================================================================
# Stubs
# ============================================================================

class _FakeVlog:
    def __init__(self):
        self.normal_lines = []
    def normal(self, msg):
        self.normal_lines.append(msg)


class _FakeCommParams:
    def __init__(self, run_on_server=True):
        self.run_on_server = run_on_server


class _FakeParams:
    def __init__(self, run_on_server=True):
        self.communication = _FakeCommParams(run_on_server)


def _apply_force_no_ai_server(params, vlog, environ):
    """Faithful replica of the override block in ai_agent.py::run().

    Kept byte-for-byte equivalent in *logic* to the shipped code; the
    source-scan tests below guarantee the real file still matches.
    """
    if environ.get("FORCE_NO_AI_SERVER", "").strip() == "1":
        if params.communication.run_on_server:
            vlog.normal(
                "NOTE: FORCE_NO_AI_SERVER=1 -- forcing run_on_server=False")
        params.communication.run_on_server = False


# ============================================================================
# Behavioural tests
# ============================================================================

def test_force_1_flips_true_to_false():
    params = _FakeParams(run_on_server=True)
    vlog = _FakeVlog()
    _apply_force_no_ai_server(params, vlog, {"FORCE_NO_AI_SERVER": "1"})
    assert params.communication.run_on_server is False, \
        "FORCE_NO_AI_SERVER=1 must force run_on_server=False"
    assert any("FORCE_NO_AI_SERVER=1" in m for m in vlog.normal_lines), \
        "expected a NOTE log line when overriding a True default"


def test_force_1_keeps_false_false_and_is_quiet():
    # Already local: still False, but no need to log a change.
    params = _FakeParams(run_on_server=False)
    vlog = _FakeVlog()
    _apply_force_no_ai_server(params, vlog, {"FORCE_NO_AI_SERVER": "1"})
    assert params.communication.run_on_server is False
    assert vlog.normal_lines == [], \
        "must not log a change when run_on_server was already False"


def test_unset_leaves_true_untouched():
    params = _FakeParams(run_on_server=True)
    vlog = _FakeVlog()
    _apply_force_no_ai_server(params, vlog, {})
    assert params.communication.run_on_server is True, \
        "unset FORCE_NO_AI_SERVER must not change run_on_server"
    assert vlog.normal_lines == []


def test_value_zero_does_not_override():
    params = _FakeParams(run_on_server=True)
    _apply_force_no_ai_server(params, _FakeVlog(), {"FORCE_NO_AI_SERVER": "0"})
    assert params.communication.run_on_server is True, \
        "FORCE_NO_AI_SERVER=0 must not override"


def test_value_true_word_does_not_override():
    # Only the literal "1" triggers the override; "true"/"yes" do not.
    for val in ("true", "True", "yes", "on"):
        params = _FakeParams(run_on_server=True)
        _apply_force_no_ai_server(
            params, _FakeVlog(), {"FORCE_NO_AI_SERVER": val})
        assert params.communication.run_on_server is True, \
            "FORCE_NO_AI_SERVER=%r must not override (only '1' does)" % val


def test_whitespace_padded_one_overrides():
    # .strip() means " 1 " counts as "1".
    params = _FakeParams(run_on_server=True)
    _apply_force_no_ai_server(
        params, _FakeVlog(), {"FORCE_NO_AI_SERVER": "  1  "})
    assert params.communication.run_on_server is False, \
        "whitespace-padded '1' must still override after strip()"


def test_empty_string_does_not_override():
    params = _FakeParams(run_on_server=True)
    _apply_force_no_ai_server(params, _FakeVlog(), {"FORCE_NO_AI_SERVER": ""})
    assert params.communication.run_on_server is True


# ============================================================================
# Source-scan tests (drift guards) — skip if the phenix-tree file is absent.
# ============================================================================

def _run_body():
    """Return the text of ai_agent.py::run(), or None if file not found.

    Uses function-boundary windows: from 'def run(self):' to the next
    top-level 'def ' at the same (2-space) method indent.
    """
    path = _find_ai_agent_source()
    if path is None:
        return None
    with open(path, "r") as fh:
        src = fh.read()
    m = re.search(r"\n  def run\(self\):\n", src)
    if not m:
        return None
    start = m.start()
    # Next method at the same 2-space indent ends run().
    nxt = re.search(r"\n  def [A-Za-z_]", src[m.end():])
    end = m.end() + nxt.start() if nxt else len(src)
    return src[start:end]


def test_source_block_present_in_run():
    body = _run_body()
    if body is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    assert "FORCE_NO_AI_SERVER" in body, \
        "run() must contain the FORCE_NO_AI_SERVER override"
    assert 'os.environ.get("FORCE_NO_AI_SERVER"' in body or \
           "os.environ.get('FORCE_NO_AI_SERVER'" in body, \
        "override must read os.environ.get('FORCE_NO_AI_SERVER')"
    assert "run_on_server = False" in body, \
        "override must set run_on_server = False"


def test_source_override_after_set_defaults_before_session_mgmt():
    body = _run_body()
    if body is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    i_defaults = body.find("set_defaults()")
    i_force = body.find("FORCE_NO_AI_SERVER")
    i_session = body.find("_handle_session_management")
    assert i_defaults != -1 and i_force != -1 and i_session != -1, \
        "expected set_defaults(), FORCE_NO_AI_SERVER, and " \
        "_handle_session_management all inside run()"
    assert i_defaults < i_force < i_session, \
        "override must sit AFTER set_defaults() and BEFORE " \
        "_handle_session_management() so all execution paths see it"


def test_source_strips_value():
    body = _run_body()
    if body is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    # The literal-"1" + .strip() contract the behavioural tests rely on.
    assert '.strip() == "1"' in body or ".strip() == '1'" in body, \
        "override must compare .strip() == '1' (literal one, whitespace-safe)"


def _ai_analysis_run_body():
    """Extract ai_analysis.py::run() body, or None if the file isn't found."""
    path = _find_ai_analysis_source()
    if path is None:
        return None
    with open(path, "r") as fh:
        src = fh.read()
    m = re.search(r"\n  def run\(self\):\n", src)
    if not m:
        return None
    start = m.start()
    nxt = re.search(r"\n  def [A-Za-z_]", src[m.end():])
    end = m.end() + nxt.start() if nxt else len(src)
    return src[start:end]


def test_source_block_present_in_ai_analysis_run():
    """ai_analysis.py::run() must carry the same FORCE_NO_AI_SERVER override
    so the env var applies to phenix.ai_analysis (all analysis_modes) as well
    as phenix.ai_agent."""
    body = _ai_analysis_run_body()
    if body is None:
        print("  (skip) ai_analysis.py not found in this checkout")
        return
    assert "FORCE_NO_AI_SERVER" in body, \
        "ai_analysis run() must contain the FORCE_NO_AI_SERVER override"
    assert 'os.environ.get("FORCE_NO_AI_SERVER"' in body or \
           "os.environ.get('FORCE_NO_AI_SERVER'" in body, \
        "override must read os.environ.get('FORCE_NO_AI_SERVER')"
    assert "run_on_server = False" in body, \
        "override must set run_on_server = False"
    assert '.strip() == "1"' in body or ".strip() == '1'" in body, \
        "override must compare .strip() == '1' (literal one, whitespace-safe)"
    # Must precede the server/local dispatch so every analysis_mode honors it.
    i_force = body.find("FORCE_NO_AI_SERVER")
    i_dispatch = body.find("run_job_on_server_or_locally")
    assert i_dispatch != -1, \
        "expected run_job_on_server_or_locally dispatch in ai_analysis run()"
    assert i_force < i_dispatch, \
        "override must sit BEFORE the run_job_on_server_or_locally dispatch"


# ============================================================================
# Runner
# ============================================================================

def _dispatcher_body(find_fn, func_name="run_job_on_server_or_locally"):
    """Extract a top-level function body from the resolved source file."""
    path = find_fn()
    if path is None:
        return None
    with open(path, "r") as fh:
        src = fh.read()
    m = re.search(r"\ndef %s\(" % re.escape(func_name), src)
    if not m:
        return None
    start = m.start()
    nxt = re.search(r"\ndef [A-Za-z_]", src[m.end():])
    end = m.end() + nxt.start() if nxt else len(src)
    return src[start:end]


def test_ai_analysis_dispatcher_has_absolute_guard():
    """run_job_on_server_or_locally in ai_analysis must refuse to contact the
    server when FORCE_NO_AI_SERVER=1 -- standard mode with no local RAG DB
    must error (Sorry), not silently fall back to the server."""
    body = _dispatcher_body(_find_ai_analysis_source)
    if body is None:
        print("  (skip) ai_analysis.py not found in this checkout")
        return
    assert "FORCE_NO_AI_SERVER" in body, \
        "ai_analysis dispatcher must guard on FORCE_NO_AI_SERVER"
    # It must raise (Sorry) rather than fall through to the server when local
    # execution is impossible, and must not leave run_on_server True.
    i_force = body.find("FORCE_NO_AI_SERVER")
    assert "Sorry(" in body[i_force:], \
        "dispatcher must raise Sorry when local execution is impossible"
    assert "run_on_server = False" in body[i_force:], \
        "dispatcher must force run_on_server=False under the flag"
    # Guard must precede the explicit-server branch.
    i_server_branch = body.find("if run_on_server:")
    assert i_force < i_server_branch, \
        "FORCE guard must come before the 'if run_on_server:' server branch"
    # CRITICAL: the guard must SHORT-CIRCUIT to local (run_job_locally +
    # return) inside the FORCE block, not merely set run_on_server=False and
    # fall through.  Falling through lets the `elif not has_database` branch
    # route a no-database LLM-only run to the server -- the hole this test
    # pins closed.  Check that run_job_locally and a return appear between the
    # FORCE check and the server-branch.
    guard_region = body[i_force:i_server_branch]
    assert "run_job_locally" in guard_region, \
        "FORCE guard must run locally directly (short-circuit), not fall through"
    assert "return working_results" in guard_region, \
        "FORCE guard must return after running locally (no fall-through to " \
        "the server/no-database branch)"


def test_ai_agent_dispatcher_has_absolute_guard():
    """The matching guard must exist in ai_agent's dispatcher too."""
    body = _dispatcher_body(_find_ai_agent_source)
    if body is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    assert "FORCE_NO_AI_SERVER" in body, \
        "ai_agent dispatcher must guard on FORCE_NO_AI_SERVER"
    i_force = body.find("FORCE_NO_AI_SERVER")
    assert "Sorry(" in body[i_force:], \
        "ai_agent dispatcher must raise Sorry when local execution impossible"
    assert "run_on_server = False" in body[i_force:], \
        "ai_agent dispatcher must force run_on_server=False under the flag"


_TESTS = [
    test_force_1_flips_true_to_false,
    test_force_1_keeps_false_false_and_is_quiet,
    test_unset_leaves_true_untouched,
    test_value_zero_does_not_override,
    test_value_true_word_does_not_override,
    test_whitespace_padded_one_overrides,
    test_empty_string_does_not_override,
    test_source_block_present_in_run,
    test_source_override_after_set_defaults_before_session_mgmt,
    test_source_strips_value,
    test_source_block_present_in_ai_analysis_run,
    test_ai_analysis_dispatcher_has_absolute_guard,
    test_ai_agent_dispatcher_has_absolute_guard,
]


def run_all_tests():
    for test_fn in _TESTS:
        test_fn()
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    passed = 0
    failed = 0
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        try:
            test_fn()
            print("  PASS: %s" % test_fn.__name__)
            passed += 1
        except Exception:
            import traceback
            print("  FAIL: %s" % test_fn.__name__)
            traceback.print_exc()
            failed += 1
    print()
    if failed:
        print("%d/%d tests FAILED." % (failed, passed + failed))
        sys.exit(1)
    else:
        print("All %d tests passed." % passed)
