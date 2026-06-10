"""Regression tests for the provider-lock ("provider-specific install") feature.

Prevent-mistakes level: phenix.install_ai_tools <provider> records a single
allowed provider in ai_config.json; the CLI programs and GUI steer to it so an
internal user does not ACCIDENTALLY ping a public LLM.

Source-extraction / source-scan style (no PHENIX import); 2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import json
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                       # .../libtbx/langchain


def _find(env, *rel_from_phenix):
    """Resolve a phenix-tree file via env override or common relative roots."""
    e = os.environ.get(env)
    if e and os.path.isfile(e):
        return e
    for base in (
        os.path.join(_ROOT, "..", "..", "..", "phenix", "phenix"),
        os.path.join(_ROOT, "..", "..", "..", "phenix"),
    ):
        p = os.path.abspath(os.path.join(base, *rel_from_phenix))
        if os.path.isfile(p):
            return p
    return None


# --- Config helpers (behavioral, against real source) ------------------------

def _load_config_helpers():
    """Extract the four ai_config helpers from utilities.py and exec them with a
    fake get_ai_dir pointing at a temp directory."""
    path = _find("UTILITIES_PY", "phenix_ai", "utilities.py")
    if path is None:
        return None, None
    src = open(path).read()
    block = re.search(
        r'(def get_ai_config_path.*?def write_ai_config.*?return None\n)',
        src, re.DOTALL)
    if not block:
        return None, None
    tmp = tempfile.mkdtemp()
    helpers = block.group(1).replace(
        "  try:\n    from libtbx.langchain.core.llm import SUPPORTED_PROVIDERS\n"
        "  except Exception:\n    SUPPORTED_PROVIDERS = "
        "[\"google\", \"openai\", \"ollama\", \"anthropic\", \"portkey\"]\n",
        "  SUPPORTED_PROVIDERS = "
        "['google','openai','ollama','anthropic','portkey']\n")
    preamble = ("import os, json, datetime\n"
                "def get_ai_dir():\n  return %r\n" % tmp)
    ns = {}
    exec(preamble + helpers, ns)
    return ns, tmp


def test_locked_provider_none_when_no_config():
    ns, _ = _load_config_helpers()
    if ns is None:
        print("  (skip) utilities.py helpers not found")
        return
    assert ns["get_locked_provider"]() is None


def test_write_and_read_lock():
    ns, _ = _load_config_helpers()
    if ns is None:
        print("  (skip)")
        return
    ns["write_ai_config"]("portkey")
    assert ns["get_locked_provider"]() == "portkey"


def test_all_clears_lock():
    ns, _ = _load_config_helpers()
    if ns is None:
        print("  (skip)")
        return
    ns["write_ai_config"]("portkey")
    assert ns["get_locked_provider"]() == "portkey"
    ns["write_ai_config"]("all")
    assert ns["get_locked_provider"]() is None
    assert not os.path.isfile(ns["get_ai_config_path"]())


def test_corrupt_and_unknown_are_none():
    ns, _ = _load_config_helpers()
    if ns is None:
        print("  (skip)")
        return
    open(ns["get_ai_config_path"](), "w").write("{not json")
    assert ns["get_locked_provider"]() is None
    open(ns["get_ai_config_path"](), "w").write(
        json.dumps({"allowed_provider": "bogus"}))
    assert ns["get_locked_provider"]() is None


def test_one_element_list_form():
    ns, _ = _load_config_helpers()
    if ns is None:
        print("  (skip)")
        return
    open(ns["get_ai_config_path"](), "w").write(
        json.dumps({"allowed_provider": ["google"]}))
    assert ns["get_locked_provider"]() == "google"


# --- install_ai_tools.py argument handling (source-scan) ---------------------

def test_install_requires_provider_arg():
    path = _find("INSTALL_AI_TOOLS_PY", "command_line", "install_ai_tools.py")
    if path is None:
        print("  (skip) install_ai_tools.py not found")
        return
    src = open(path).read()
    assert "if not args:" in src, "must fail when no provider arg is given"
    assert "_usage" in src, "must show usage on missing/invalid arg"
    assert "_VALID_CHOICES" in src and "portkey" in src
    # provider-specific key check replaced the hardcoded GOOGLE_API_KEY-only one
    assert "PORTKEY_AZURE_API_KEY" in src and "_PROVIDER_KEY_REQUIREMENTS" in src
    # forwards provider to the csh and writes the lock
    assert "build_file, provider" in src or "%s %s" in src
    assert "write_ai_config" in src


# --- CLI steering (source-scan in both programs) -----------------------------

def _assert_lock_steering(path, label):
    src = open(path).read()
    assert "get_locked_provider" in src, \
        "%s must consult get_locked_provider" % label
    assert "self.params.communication.provider = _locked" in src, \
        "%s must steer provider to the lock" % label
    # must be near the FORCE_NO_AI_SERVER override (same pre-dispatch point)
    i_force = src.find('FORCE_NO_AI_SERVER", "").strip() == "1"')
    i_lock = src.find("get_locked_provider")
    assert i_force != -1 and i_lock != -1 and i_lock > i_force, \
        "%s lock steering should follow the FORCE_NO_AI_SERVER block" % label


def test_ai_analysis_steers_to_lock():
    path = _find("AI_ANALYSIS_PY", "programs", "ai_analysis.py")
    if path is None:
        print("  (skip) ai_analysis.py not found")
        return
    _assert_lock_steering(path, "ai_analysis")


def test_ai_agent_steers_to_lock():
    path = _find("AI_AGENT_PY", "programs", "ai_agent.py")
    if path is None:
        print("  (skip) ai_agent.py not found")
        return
    _assert_lock_steering(path, "ai_agent")


# --- GUI branch on FORCE_NO_AI_SERVER (source-scan) --------------------------

def _assert_gui_omits_run_on_server(env_var, default_relpath, label):
    e = os.environ.get(env_var)
    path = e if (e and os.path.isfile(e)) else None
    if path is None:
        base = os.path.join(_ROOT, "..", "..", "..", "wxGUI2", "Programs")
        p = os.path.abspath(os.path.join(base, default_relpath))
        if os.path.isfile(p):
            path = p
    if path is None:
        print("  (skip) GUI %s not found" % label)
        return
    src = open(path).read()
    i_guard = src.find('FORCE_NO_AI_SERVER", "").strip() != "1"')
    i_draw = src.find("draw_phil_controls('communication.run_on_server')")
    assert i_guard != -1, "%s must branch on FORCE_NO_AI_SERVER" % label
    assert i_draw != -1 and i_draw > i_guard, \
        "%s run_on_server control must be inside the not-forced branch" % label


def test_gui_omits_run_on_server_when_forced():
    """AIAgent.py must omit the run_on_server control when FORCE_NO_AI_SERVER=1."""
    _assert_gui_omits_run_on_server("GUI_AIAGENT_PY", "AIAgent.py", "AIAgent")


def test_gui_analysis_omits_run_on_server_when_forced():
    """AIAnalysis.py must omit run_on_server under FORCE_NO_AI_SERVER=1 too
    (parity with AIAgent)."""
    _assert_gui_omits_run_on_server("GUI_AIANALYSIS_PY", "AIAnalysis.py",
                                    "AIAnalysis")


_TESTS = [
    test_locked_provider_none_when_no_config,
    test_write_and_read_lock,
    test_all_clears_lock,
    test_corrupt_and_unknown_are_none,
    test_one_element_list_form,
    test_install_requires_provider_arg,
    test_ai_analysis_steers_to_lock,
    test_ai_agent_steers_to_lock,
    test_gui_omits_run_on_server_when_forced,
    test_gui_analysis_omits_run_on_server_when_forced,
]


def run_all_tests():
    for fn in _TESTS:
        fn()
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    p = f = 0
    for fn in _TESTS:
        try:
            fn(); print("  PASS: %s" % fn.__name__); p += 1
        except AssertionError as e:
            print("  FAIL: %s -- %s" % (fn.__name__, e)); f += 1
    print("\n%d passed, %d failed" % (p, f))
