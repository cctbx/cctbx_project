"""Regression tests for the dropped-resolution premature-stop bug.

At sub-2A resolution phenix.ai_agent could declare "R-free TARGET REACHED" after
a single refinement and stop, because the stop evaluator received resolution=None
(falling back to the resolution-INDEPENDENT R-free target 0.28) instead of the
known data resolution (e.g. 1.57 -> banded target 0.25).

Two fixes are guarded:
1. perceive() (graph_nodes.py) prefers state["session_resolution"] over the
   history-derived value.  [source-scan]
2. derive_metrics_from_history() (metrics_analyzer.py) recovers resolution from
   result text when the structured "analysis" dict omits it.  [behavioral]

The behavioral half drives the REAL metric_evaluator + the REAL metrics.yaml via
a faithful yaml_loader shim, so the 0.28-vs-0.25 band selection (the crux) is
exercised, not mocked away.  2-space indent; no PHENIX import.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import json
import types
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                       # .../libtbx/langchain
_AGENT = os.path.join(_ROOT, "agent")
_KNOWLEDGE = os.path.join(_ROOT, "knowledge")
_DATA = os.path.join(_HERE, "data")


def _find(env, *parts):
    e = os.environ.get(env)
    if e and os.path.isfile(e):
        return e
    p = os.path.join(*parts)
    return p if os.path.isfile(p) else None


# --- Fix 1: perceive() prefers session_resolution (source-scan) --------------

def test_perceive_prefers_session_resolution():
    path = _find("GRAPH_NODES_PY", _AGENT, "graph_nodes.py")
    if path is None:
        print("  (skip) graph_nodes.py not found")
        return
    src = open(path).read()
    # perceive() must resolve resolution through the shared resolver (which
    # prefers session_resolution), passing the metrics_history it has in scope.
    assert re.search(
        r'resolution,\s*_resolution_source\s*=\s*resolve_session_resolution\(\s*\n?\s*state,\s*metrics_history=metrics_history\)',
        src), \
        "perceive() must call resolve_session_resolution(state, metrics_history=metrics_history)"
    assert "from libtbx.langchain.utils.run_utils import _coerce_resolution" in src, \
        "graph_nodes must import _coerce_resolution (the shared coercion helper, " \
        "which applies _safe_float internally)"
    # The old bare history-only form and the old inline _safe_float-or form must be gone.
    assert not re.search(r'^\s*resolution = get_latest_resolution\(metrics_history\)\s*$',
                         src, re.MULTILINE), \
        "the bare 'resolution = get_latest_resolution(metrics_history)' must be replaced"
    # The nested find_resolution() closure must be gone (unified to module level).
    assert "def find_resolution():" not in src, \
        "the nested find_resolution() closure must be replaced by the module-level resolver"
    # build() routes through the shared resolver INDIRECTLY now: perceive()
    # resolves resolution (above) and stores it; the unified builder consumes
    # state["session_resolution"] rather than re-deriving it.  The legacy
    # build-side resolve_session_resolution(state, workflow_state=workflow_state)
    # call lived in the USE_NEW_COMMAND_BUILDER=False branch, which was removed
    # as dead code; build() now delegates to _build_with_new_builder().
    assert "_build_with_new_builder(state)" in src, \
        "build() must delegate to the unified _build_with_new_builder(state)"
    assert re.search(r'state\.get\(\s*["\']session_resolution["\']\s*\)', src), \
        "the unified build path must consume the resolver's output via " \
        "state['session_resolution'] (resolved upstream by perceive())"


def _load_resolver():
    """Extract resolve_session_resolution from graph_nodes and exec it with
    minimal deps (graph_nodes has heavy imports; the resolver itself needs only
    _safe_float, get_latest_resolution, and re)."""
    path = _find("GRAPH_NODES_PY", _AGENT, "graph_nodes.py")
    if path is None:
        return None
    src = open(path).read()
    m = re.search(r'def resolve_session_resolution\(.*?\n    return None, None\n',
                  src, re.DOTALL)
    if not m:
        return None

    def _safe_float(v):
        if v is None:
            return None
        try:
            return float(v)
        except (ValueError, TypeError):
            return None

    def get_latest_resolution(mh):
        for x in reversed(mh):
            r = _safe_float(x.get("resolution"))
            if r is not None:
                return r
        return None

    def _coerce_resolution(v):
        v = _safe_float(v)
        if v is not None and 0.5 < v < 20.0:
            return v
        return None

    ns = {"_safe_float": _safe_float, "_coerce_resolution": _coerce_resolution,
          "get_latest_resolution": get_latest_resolution, "re": re}
    exec(m.group(0), ns)
    return ns["resolve_session_resolution"]


def test_resolver_priority_tiers():
    """The unified resolver returns the right source in priority order."""
    R = _load_resolver()
    if R is None:
        print("  (skip) resolver not found")
        return
    # 1 session wins over everything below it
    assert R({"session_resolution": 1.57,
              "log_analysis": {"resolution": 3.0}},
             workflow_state={"resolution": 2.0},
             metrics_history=[{"resolution": 2.5}]) == (1.57, "session")
    # 2 workflow_state
    assert R({}, workflow_state={"resolution": 1.8}) == (1.8, "workflow_state")
    # 3 log_analysis
    assert R({"log_analysis": {"resolution": 1.9}}) == (1.9, "log_analysis")
    # 4 metrics_history (perceive's source)
    assert R({}, metrics_history=[{"resolution": 2.1}]) == (2.1, "metrics_history")
    # 5 previous command (last resort)
    assert R({"history": [{"command": "phenix.refine resolution=2.3"}]}) \
        == (2.3, "previous_command")
    # nothing known
    assert R({}) == (None, None)


def test_resolver_coercion_and_range_guard():
    """Non-numeric / out-of-range values are coerced or skipped, never returned."""
    R = _load_resolver()
    if R is None:
        print("  (skip)")
        return
    # string session is coerced
    assert R({"session_resolution": "1.57"}) == (1.57, "session")
    # garbage / 0.0 / negative session fall through to the next source
    for bad in ("abc", 0.0, -1.5, 25.0):
        assert R({"session_resolution": bad},
                 metrics_history=[{"resolution": 2.0}]) == (2.0, "metrics_history"), \
            "bad session %r must fall through" % (bad,)
    # string workflow_state coerced; out-of-range workflow_state skipped
    assert R({}, workflow_state={"resolution": "2.0"}) == (2.0, "workflow_state")
    assert R({}, workflow_state={"resolution": 25.0},
             metrics_history=[{"resolution": 2.0}]) == (2.0, "metrics_history")


def test_perceive_and_build_agree_on_shared_sources():
    """perceive() and build() (same resolver, different optional args) agree
    whenever resolution comes from a source they share (session / workflow_state
    / log_analysis / history)."""
    R = _load_resolver()
    if R is None:
        print("  (skip)")
        return
    as_perceive = lambda st: R(st, metrics_history=st.get("_mh", []))[0]
    as_build = lambda st: R(st, workflow_state=st.get("workflow_state", {}))[0]
    shared_cases = [
        {"session_resolution": 1.74, "workflow_state": {"resolution": 2.5}},
        {"workflow_state": {"resolution": 1.9}},
        {"log_analysis": {"resolution": 2.2}},
        {"history": [{"command": "phenix.refine resolution=2.8"}]},
        {},  # both None
    ]
    for st in shared_cases:
        assert as_perceive(st) == as_build(st), \
            "perceive and build disagree on %r: %s vs %s" % (
                st, as_perceive(st), as_build(st))


# --- Real-code harness for the metrics_analyzer / evaluator behavior ---------

def _load_metrics_stack():
    """Load the REAL metrics_analyzer + metric_evaluator.  Prefer the REAL
    yaml_loader (+ real metrics.yaml) so band selection is genuine end-to-end;
    fall back to a faithful yaml_loader shim if the loader isn't present.
    Returns (metrics_analyzer_module, metric_evaluator_module) or (None, None)."""
    ma_path = _find("METRICS_ANALYZER_PY", _AGENT, "metrics_analyzer.py")
    me_path = _find("METRIC_EVALUATOR_PY", _AGENT, "metric_evaluator.py")
    yaml_path = _find("METRICS_YAML", _KNOWLEDGE, "metrics.yaml")
    loader_path = _find("YAML_LOADER_PY", _KNOWLEDGE, "yaml_loader.py")
    if not (ma_path and me_path and yaml_path):
        return None, None
    try:
        import yaml
    except ImportError:
        return None, None

    # Register lightweight parent packages so the by-path module loads below
    # resolve.  GUARD with `if name not in sys.modules`: in a real PHENIX env
    # `libtbx`, `libtbx.langchain`, etc. are already real packages, and
    # UNCONDITIONALLY overwriting them replaces `libtbx.langchain` (path=None)
    # with a __path__-less module.  That persists in sys.modules for the rest of
    # the suite and makes a LATER test's `from libtbx.langchain.rag... import`
    # fail with "'libtbx.langchain' is not a package" (observed: tst_dependencies
    # reranker probe failing only inside run_all_tests).  Only inject a fake when
    # the real one is absent (bare sandbox).
    for name, path in [("libtbx", None), ("libtbx.langchain", None),
                       ("libtbx.langchain.agent", _AGENT),
                       ("libtbx.langchain.utils", None),
                       ("libtbx.langchain.knowledge", _KNOWLEDGE)]:
        if name not in sys.modules:
            m = types.ModuleType(name)
            if path:
                m.__path__ = [path]
            sys.modules[name] = m

    ru = types.ModuleType("libtbx.langchain.utils.run_utils")
    def _safe_float(v):
        if v is None:
            return None
        try:
            return float(v)
        except (ValueError, TypeError):
            return None
    ru._safe_float = _safe_float
    sys.modules["libtbx.langchain.utils.run_utils"] = ru

    def load(mod, path):
        spec = importlib.util.spec_from_file_location(mod, path)
        o = importlib.util.module_from_spec(spec)
        sys.modules[mod] = o
        spec.loader.exec_module(o)
        return o

    if loader_path is not None:
        # Use the REAL yaml_loader (resolves metrics.yaml relative to itself).
        try:
            load("libtbx.langchain.knowledge.yaml_loader", loader_path)
        except Exception:
            loader_path = None  # fall through to shim

    if loader_path is None:
        # Faithful shim mirroring get_metric_threshold band logic.
        data = yaml.safe_load(open(yaml_path))

        def get_metric_threshold(metric_name, level="acceptable", resolution=None):
            spec = data.get(metric_name, {})
            if resolution and "by_resolution" in spec:
                for band in spec["by_resolution"]:
                    lo, hi = band.get("range", [0, 999])
                    if lo <= resolution < hi:
                        return band.get(level)
            return spec.get("thresholds", {}).get(level)

        yl = types.ModuleType("libtbx.langchain.knowledge.yaml_loader")
        yl.load_metrics = lambda: data
        yl.get_metric_threshold = get_metric_threshold
        yl.get_metric_direction = lambda *a, **k: "minimize"
        yl.is_metric_good = lambda *a, **k: False
        yl.is_metric_acceptable = lambda *a, **k: False
        sys.modules["libtbx.langchain.knowledge.yaml_loader"] = yl

    me = load("libtbx.langchain.agent.metric_evaluator", me_path)
    ma = load("libtbx.langchain.agent.metrics_analyzer", ma_path)
    return ma, me


def test_band_target_selection():
    """Sanity-check the crux: None -> 0.28 (resolution-independent), 1.57 -> 0.25
    (the [1.5,2.0] band).  If this fails, the rest of the bug is moot."""
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip) metrics stack / pyyaml not available")
        return
    ev = me.MetricEvaluator()
    assert ev.get_target("r_free", None) == 0.28, "None must give the 0.28 default"
    assert ev.get_target("r_free", 1.57) == 0.25, "1.57 must give the banded 0.25"


def test_resolution_none_reports_target_reached_but_1p57_does_not():
    """The bug: with resolution=None the first refine reads 'TARGET REACHED';
    with the correct 1.57 it does not."""
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip)")
        return
    cycles = [{"cycle": 1, "program": "phenix.refine",
               "result": "R-work = 0.21  R-free = 0.252",
               "analysis": {"r_free": 0.252}}]
    mh = ma.derive_metrics_from_history(cycles)
    bug = ma.analyze_metrics_trend(mh, None, "xray", use_yaml_evaluator=True)
    fix = ma.analyze_metrics_trend(mh, 1.57, "xray", use_yaml_evaluator=True)
    assert "TARGET REACHED" in (bug.get("trend_summary") or ""), \
        "resolution=None should (wrongly) report TARGET REACHED"
    assert "TARGET REACHED" not in (fix.get("trend_summary") or ""), \
        "resolution=1.57 must NOT report TARGET REACHED after one refinement"


def test_resolution_recovered_from_result_text():
    """Fix 2: derive_metrics_from_history recovers resolution from result text
    when the structured analysis dict omits it.

    Uses the real xtriage form where 'Anomalous Resolution: 2.10' precedes the
    data 'Resolution: 1.74', so this also guards the WIRING: it passes only if
    derive_metrics_from_history routes through _extract_resolution.  A revert to
    a naive first-match _extract_float pattern would recover 2.10 here and fail.
    """
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip)")
        return
    # analysis dict has NO resolution; it appears only in result text, and the
    # anomalous line comes FIRST (the case a naive first-match gets wrong).
    cycles = [{"cycle": 1, "program": "phenix.xtriage",
               "result": "Anomalous Resolution: 2.10\nResolution: 1.74",
               "analysis": {}}]
    mh = ma.derive_metrics_from_history(cycles)
    assert mh[-1]["resolution"] == 1.74, \
        "resolution must be recovered as the DATA value 1.74, not anomalous 2.10 " \
        "(got %r) -- check derive_metrics_from_history uses _extract_resolution" \
        % mh[-1]["resolution"]
    assert ma.get_latest_resolution(mh) == 1.74


def test_extract_resolution_handles_real_hazards():
    """_extract_resolution must read the DATA resolution from real PHENIX text,
    not the anomalous resolution or the low-resolution end of a range."""
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip)")
        return
    er = ma._extract_resolution
    # The real p9 xtriage form: anomalous line precedes the data line.
    assert er("Anomalous Resolution: 2.10\nResolution: 1.74") == 1.74
    # Range form: take the high-resolution (second) limit.
    assert er("Resolution range: 20.0 - 1.74") == 1.74
    assert er("Resolution range: 50.00 \u2013 2.30") == 2.30  # en-dash
    # Simple forms.
    assert er("resolution = 1.57") == 1.57
    assert er("Resolution: 2.50") == 2.5
    # Multiple lines: the last (final summary) wins.
    assert er("Resolution: 3.00\nFinal Resolution: 1.85") == 1.85
    # Noise / absence -> None.
    assert er("Anomalous Resolution: 2.10") is None      # anomalous only
    assert er("R-free = 0.21") is None                   # no resolution
    assert er("Resolution: 25.0") is None                # out of sane range
    assert er("") is None and er(None) is None
    # Review-surfaced guards:
    # integer header must NOT be read as a 1.0 A resolution (decimal required).
    assert er("Completeness in resolution range: 1") is None
    # anomalous + data on the SAME line: strip anomalous, keep the data value.
    assert er("Anomalous Resolution: 2.10 Resolution: 1.74") == 1.74
    # boundary bounds (0.5 < res < 20, exclusive).
    assert er("Resolution: 0.50") is None
    assert er("Resolution: 20.0") is None
    assert er("Resolution: 0.51") == 0.51


def test_extract_resolution_beats_naive_first_match():
    """Negative control: a naive first-match `_extract_float`-style pattern would
    return the WRONG value on the real xtriage text; _extract_resolution fixes it.
    This guards against anyone reverting to the simple pattern."""
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip)")
        return
    import re
    text = "Anomalous Resolution: 2.10\nResolution: 1.74"
    # the old shipped pattern (first-match) -> grabs the anomalous 2.10
    naive = re.search(r"[Rr]esolution[:\s=]+([0-9.]+)", text, re.IGNORECASE)
    assert naive and float(naive.group(1)) == 2.10, \
        "sanity: the naive pattern really does grab the wrong (anomalous) value"
    # the corrected helper grabs the data resolution
    assert ma._extract_resolution(text) == 1.74, \
        "_extract_resolution must beat the naive first-match (got %r)" \
        % ma._extract_resolution(text)


def test_fixture_session_resolves_correctly():
    """End-to-end on the committed fixture: the session's authoritative
    resolution (1.57) plus the text-fallback both yield 1.57, and at 1.57 the
    first refine is not 'TARGET REACHED'."""
    ma, me = _load_metrics_stack()
    if ma is None:
        print("  (skip)")
        return
    fpath = _find("RFREE_FIXTURE", _DATA, "session_rfree_resolution_bug.json")
    if fpath is None:
        print("  (skip) fixture not found")
        return
    sess = json.load(open(fpath))
    mh = ma.derive_metrics_from_history(sess["cycles"])
    # the refine cycle's resolution is recovered (text-fallback) ...
    assert ma.get_latest_resolution(mh) == 1.57
    # ... and the authoritative session value agrees
    assert sess["session_resolution"] == 1.57
    trend = ma.analyze_metrics_trend(mh, sess["session_resolution"], "xray",
                                     use_yaml_evaluator=True)
    assert "TARGET REACHED" not in (trend.get("trend_summary") or "")


_TESTS = [
    test_perceive_prefers_session_resolution,
    test_resolver_priority_tiers,
    test_resolver_coercion_and_range_guard,
    test_perceive_and_build_agree_on_shared_sources,
    test_band_target_selection,
    test_resolution_none_reports_target_reached_but_1p57_does_not,
    test_resolution_recovered_from_result_text,
    test_extract_resolution_handles_real_hazards,
    test_extract_resolution_beats_naive_first_match,
    test_fixture_session_resolves_correctly,
]


def _snapshot_sysmodules():
    """Snapshot the sys.modules entries this test may add or overwrite, so
    run_all_tests() can restore them and never leak test-loaded modules into
    the rest of the suite (e.g. by-path metric_evaluator/metrics_analyzer or a
    fake yaml_loader shadowing the real package for a later test)."""
    return {k: v for k, v in sys.modules.items()
            if k == "libtbx" or k.startswith("libtbx.")
            or k in ("langchain_chroma",)}


def _restore_sysmodules(saved):
    """Restore sys.modules to the saved snapshot for the tracked keys: delete
    anything added since, put back anything overwritten."""
    for k in [k for k in list(sys.modules)
              if (k == "libtbx" or k.startswith("libtbx.")
                  or k == "langchain_chroma") and k not in saved]:
        sys.modules.pop(k, None)
    for k, v in saved.items():
        sys.modules[k] = v


def run_all_tests():
    saved = _snapshot_sysmodules()
    try:
        for fn in _TESTS:
            fn()
    finally:
        _restore_sysmodules(saved)
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    saved = _snapshot_sysmodules()
    p = f = 0
    try:
        for fn in _TESTS:
            try:
                fn(); print("  PASS: %s" % fn.__name__); p += 1
            except AssertionError as e:
                print("  FAIL: %s -- %s" % (fn.__name__, e)); f += 1
    finally:
        _restore_sysmodules(saved)
    print("\n%d passed, %d failed" % (p, f))
