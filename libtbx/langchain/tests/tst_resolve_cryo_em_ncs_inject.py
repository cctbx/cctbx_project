"""Regression test: phenix.resolve_cryo_em must never receive ncs_file.

Bug (ai_agent_503): in a cryo-EM half-map density-modification workflow, a
stray .ncs_spec file produced by an earlier map-symmetry/segmentation step was
auto-mapped to `ncs_file=` and passed to phenix.resolve_cryo_em, which does not
accept a top-level ncs_file ("Some PHIL parameters are not recognized" /
"Unknown command line parameter definition: ncs_file = ...").  The reactive
bad_inject blacklist only learned to skip it AFTER one failed cycle.

Fix: a STATIC never-inject map (_STATIC_BAD_INJECT_PARAMS) in AgentSession,
unioned into get_bad_inject_params(), so ncs_file is excluded from
resolve_cryo_em on the FIRST command -- and both injection paths
(postprocess_command via get_bad_inject_params, and the client-side
_inject_missing_required_files) honor it.

Source-scan style (no PHENIX import needed); 2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain

_SESSION_CANDIDATES = [
    os.path.join(_ROOT, "agent", "session.py"),
]
_AI_AGENT_CANDIDATES = [
    os.path.join(_ROOT, "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "programs", "ai_agent.py"),
]
_POSTPROC_CANDIDATES = [
    os.path.join(_ROOT, "agent", "command_postprocessor.py"),
]


def _find(cands, env):
    for c in cands:
        c = os.path.abspath(c)
        if os.path.isfile(c):
            return c
    e = os.environ.get(env)
    if e and os.path.isfile(e):
        return e
    return None


def _build_session_stub():
    """Exec the real static map + the bad-inject methods into a stub class
    so the behavior is tested against the actual source, not a paraphrase."""
    path = _find(_SESSION_CANDIDATES, "SESSION_PY")
    if path is None:
        return None
    src = open(path).read()
    static_block = re.search(
        r'_STATIC_BAD_INJECT_PARAMS = \{.*?\}', src, re.DOTALL)
    get_method = re.search(
        r'    def get_bad_inject_params\(self.*?return learned \| static',
        src, re.DOTALL)
    get_all = re.search(
        r'    def get_all_bad_inject_params\(self.*?return merged',
        src, re.DOTALL)
    rec_method = re.search(
        r'    def record_bad_inject_param\(self.*?self\.save\(\)',
        src, re.DOTALL)
    if not (static_block and get_method and get_all and rec_method):
        return None
    cls = ("class S:\n"
           "    def __init__(self): self.data = {}\n"
           "    def save(self): pass\n"
           "    " + static_block.group(0).replace('\n', '\n    ') + "\n"
           + get_method.group(0) + "\n"
           + get_all.group(0) + "\n"
           + rec_method.group(0))
    ns = {}
    exec(cls, ns)
    return ns['S']


def test_static_blacklist_present_in_source():
    path = _find(_SESSION_CANDIDATES, "SESSION_PY")
    if path is None:
        print("  (skip) session.py not found")
        return
    src = open(path).read()
    assert "_STATIC_BAD_INJECT_PARAMS" in src, \
        "session.py must define a static never-inject map"
    assert '"phenix.resolve_cryo_em"' in src and "ncs_file" in src, \
        "static map must blacklist ncs_file for phenix.resolve_cryo_em"
    assert "return learned | static" in src, \
        "get_bad_inject_params must union learned with static"


def test_resolve_cryo_em_ncs_blacklisted_on_fresh_session():
    S = _build_session_stub()
    if S is None:
        print("  (skip) could not build session stub from source")
        return
    s = S()
    got = s.get_bad_inject_params("phenix.resolve_cryo_em")
    assert "ncs_file" in got, \
        "ncs_file must be blacklisted for resolve_cryo_em on a FRESH session"


def test_other_programs_unaffected():
    S = _build_session_stub()
    if S is None:
        print("  (skip)")
        return
    s = S()
    assert s.get_bad_inject_params("phenix.refine") == set(), \
        "programs with no static entry must start with an empty blacklist"


def test_learned_params_still_union():
    S = _build_session_stub()
    if S is None:
        print("  (skip)")
        return
    s = S()
    s.record_bad_inject_param("phenix.resolve_cryo_em", "another_key")
    got = s.get_bad_inject_params("phenix.resolve_cryo_em")
    assert got == {"ncs_file", "another_key"}, \
        "learned params must union with the static set"


def test_required_file_injector_honors_blacklist():
    """The client-side _inject_missing_required_files must skip blacklisted
    slots too (defense-in-depth across both injection paths)."""
    path = _find(_AI_AGENT_CANDIDATES, "AI_AGENT_PY")
    if path is None:
        print("  (skip) ai_agent.py not found")
        return
    src = open(path).read()
    assert "get_bad_inject_params(program_name)" in src, \
        "required-file injector must fetch the blacklist"
    assert "_bad_inject" in src and "continue" in src, \
        "required-file injector must skip blacklisted slots"


def test_merged_dict_includes_static_for_server():
    """get_all_bad_inject_params() must merge static + learned, so the static
    entries reach the SERVER (which only sees what the client transmits).
    This is the gap that made the original fix ineffective in server mode."""
    S = _build_session_stub()
    if S is None:
        print("  (skip)")
        return
    s = S()
    # Fresh session: static entry must already be in the merged dict.
    merged = s.get_all_bad_inject_params()
    assert merged.get("phenix.resolve_cryo_em") == ["ncs_file"], \
        "fresh merged dict must contain the static ncs_file entry"
    # Learned entries union in without dropping the static one or dupes.
    s.data["bad_inject_params"] = {"phenix.resolve_cryo_em": ["other"]}
    merged2 = s.get_all_bad_inject_params()
    assert set(merged2["phenix.resolve_cryo_em"]) == {"ncs_file", "other"}, \
        "merged dict must union static + learned with no duplicates"


def test_client_sends_merged_blacklist_to_server():
    """Source-scan: the client must send the MERGED blacklist in session_info,
    not the learned-only session.data['bad_inject_params'] (which would omit
    static entries and let the server inject them in server mode)."""
    path = _find(_AI_AGENT_CANDIDATES, "AI_AGENT_PY")
    if path is None:
        print("  (skip) ai_agent.py not found")
        return
    src = open(path).read()
    assert '"bad_inject_params": session.get_all_bad_inject_params()' in src, \
        "client must transmit get_all_bad_inject_params() (merged) to the server"


def test_json_roundtrip_then_sanitize_strips_ncs_file():
    """Integration gate (Gemini): the merged blacklist must survive a JSON
    serialization round-trip (client -> wire -> server) and still cause the
    server-side build extraction + sanitize_command to strip ncs_file.

    Exercises the actual serialization boundary plus the hardened
    {program:[keys]} -> per-program-set extraction used by the graph build
    node, then the real sanitize_command."""
    import json

    S = _build_session_stub()
    pp_path = _find(_POSTPROC_CANDIDATES, "POSTPROC_PY")
    if S is None or pp_path is None:
        print("  (skip) need session stub + command_postprocessor")
        return

    s = S()
    payload = s.get_all_bad_inject_params()           # client side
    wire = json.loads(json.dumps(payload))            # client -> server wire

    # Server build-node extraction (mirrors graph_nodes.py, hardened form):
    program = "phenix.resolve_cryo_em"
    all_bad = wire or {}
    if not isinstance(all_bad, dict):
        all_bad = {}
    prog_bad = all_bad.get(program) or []
    bad_set = set(prog_bad)
    assert "ncs_file" in bad_set, \
        "ncs_file must survive the JSON round-trip into the per-program set"

    # Real server-side sanitize must then strip the LLM-emitted ncs_file token.
    import importlib.util
    spec = importlib.util.spec_from_file_location("cmd_pp_rt", pp_path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:
        print("  (skip) could not import command_postprocessor: %s" % e)
        return
    cmd = ("phenix.resolve_cryo_em half1.ccp4 half2.ccp4 "
           "seq_file=seq.fa ncs_file=/path/x.ncs_spec")
    out = mod.sanitize_command(cmd, program_name=program,
                               bad_inject_params=bad_set)
    out_cmd = out[0] if isinstance(out, tuple) else out
    assert "ncs_file" not in out_cmd, \
        "after JSON round-trip, sanitize must still strip ncs_file; got: %s" % (
            out_cmd,)


def test_extraction_handles_null_payload():
    """Defensive (Gemini Risk 1): a bad_inject_params that deserializes to None,
    or a per-program value of None, must not raise in the build extraction."""
    # None whole-dict
    all_bad = None or {}
    if not isinstance(all_bad, dict):
        all_bad = {}
    assert set(all_bad.get("phenix.resolve_cryo_em") or []) == set()
    # None per-program value
    all_bad = {"phenix.resolve_cryo_em": None}
    prog_bad = all_bad.get("phenix.resolve_cryo_em") or []
    assert set(prog_bad) == set()  # must not raise on set(None)


def test_build_node_extraction_hardened_in_source():
    """Source-scan: the graph build node must guard the bad_inject extraction
    against None/non-dict (post-JSON) types."""
    cands = [
        os.path.join(_ROOT, "agent", "graph_nodes.py"),
    ]
    path = _find(cands, "GRAPH_NODES_PY")
    if path is None:
        print("  (skip) graph_nodes.py not found")
        return
    src = open(path).read()
    assert 'state.get("bad_inject_params") or {}' in src, \
        "build node must coerce a null bad_inject_params payload to {}"
    assert "isinstance(all_bad, dict)" in src, \
        "build node must guard against a non-dict payload"
    assert "all_bad.get(program) or []" in src, \
        "build node must coerce a null per-program value to []"


def test_sanitize_command_strips_blacklisted_ncs_file():
    """Behavioral: the REAL server-side sanitize_command must STRIP an existing
    ncs_file= token when 'ncs_file' is in bad_inject_params (not merely decline
    to inject it).  This is the server-mode guarantee -- the LLM emitted
    ncs_file= in its command, so it must be removed before execution."""
    path = _find(_POSTPROC_CANDIDATES, "POSTPROC_PY")
    if path is None:
        print("  (skip) command_postprocessor.py not found")
        return
    # Import the real module in isolation (stdlib + light deps only).
    import importlib.util
    spec = importlib.util.spec_from_file_location("cmd_pp", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:
        print("  (skip) could not import command_postprocessor: %s" % e)
        return
    if not hasattr(mod, "sanitize_command"):
        print("  (skip) sanitize_command not present")
        return
    cmd = ("phenix.resolve_cryo_em half1.ccp4 half2.ccp4 "
           "seq_file=seq.fa ncs_file=/path/shifted_used_ncs.ncs_spec")
    out = mod.sanitize_command(
        cmd, program_name="phenix.resolve_cryo_em",
        bad_inject_params={"ncs_file"})
    # sanitize_command may return a str or a tuple depending on signature;
    # normalize to the command string.
    out_cmd = out[0] if isinstance(out, tuple) else out
    assert "ncs_file" not in out_cmd, \
        "sanitize_command must STRIP the blacklisted ncs_file token; got: %s" % (
            out_cmd,)
    assert "seq_file=seq.fa" in out_cmd, \
        "sanitize_command must not strip legitimate params"


def test_sanitize_preserves_legitimate_symmetry_file():
    """The blacklist strips ncs_file but must PRESERVE the legitimate
    input_files.symmetry_file (the real way resolve_cryo_em takes symmetry).
    Guards against a short-key collision regression."""
    path = _find(_POSTPROC_CANDIDATES, "POSTPROC_PY")
    if path is None:
        print("  (skip) command_postprocessor.py not found")
        return
    import importlib.util
    spec = importlib.util.spec_from_file_location("cmd_pp_sym", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:
        print("  (skip) could not import command_postprocessor: %s" % e)
        return
    cmd = ("phenix.resolve_cryo_em half1.ccp4 half2.ccp4 seq_file=seq.fa "
           "input_files.symmetry_file=/path/sym.ncs_spec "
           "ncs_file=/path/stray.ncs_spec")
    out = mod.sanitize_command(cmd, program_name="phenix.resolve_cryo_em",
                               bad_inject_params={"ncs_file"})
    out_cmd = out[0] if isinstance(out, tuple) else out
    assert "ncs_file" not in out_cmd, "stray ncs_file must be stripped"
    assert "symmetry_file" in out_cmd, \
        "legitimate input_files.symmetry_file must be PRESERVED (no short-key " \
        "collision with the ncs_file blacklist)"


_TESTS = [
    test_static_blacklist_present_in_source,
    test_resolve_cryo_em_ncs_blacklisted_on_fresh_session,
    test_other_programs_unaffected,
    test_learned_params_still_union,
    test_required_file_injector_honors_blacklist,
    test_merged_dict_includes_static_for_server,
    test_client_sends_merged_blacklist_to_server,
    test_sanitize_command_strips_blacklisted_ncs_file,
    test_json_roundtrip_then_sanitize_strips_ncs_file,
    test_extraction_handles_null_payload,
    test_build_node_extraction_hardened_in_source,
    test_sanitize_preserves_legitimate_symmetry_file,
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
