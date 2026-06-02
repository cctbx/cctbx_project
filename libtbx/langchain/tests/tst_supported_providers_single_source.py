"""
SUPPORTED_PROVIDERS single source of truth (v120 Phase 2).

core/llm.py is the canonical home of SUPPORTED_PROVIDERS.  agent/graph_nodes.py
and phenix/programs/ai_agent.py both import it from there rather than keeping
their own copies, so a provider can never be wired into one layer but silently
rejected by another (the dual-literal drift hazard called out in the plan).

Two layers:
  1. Runtime identity — where the modules import cleanly, assert graph_nodes
     exposes the SAME list value as core.llm.  (graph_nodes pulls in `libtbx`,
     which is absent outside a Phenix build; that case is skipped, not failed.)
  2. Source-scan — always-on drift guard: both consumer files must import
     SUPPORTED_PROVIDERS from core.llm and must NOT define a local literal
     fallback list.

Per the agent test guidelines: file-exists guards, no fixed char offsets,
AssertionError never swallowed.
"""
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain
for _p in [_ROOT, os.path.join(_ROOT, "agent"), os.path.join(_ROOT, "core")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

_GRAPH_NODES = os.path.join(_ROOT, "agent", "graph_nodes.py")

_AI_AGENT_CANDIDATES = [
    # Sandbox layout: programs/ as a sibling of tests/ under langchain.
    os.path.join(_ROOT, "programs", "ai_agent.py"),
    # Real PHENIX build: ai_agent.py lives in the phenix tree.  _ROOT is
    # .../modules/cctbx_project/libtbx/langchain, so 3 ups reaches .../modules,
    # then into phenix/phenix/programs.
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "programs", "ai_agent.py"),
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


# ============================================================================
# Runtime identity
# ============================================================================

def test_core_llm_defines_supported_providers():
    from core.llm import SUPPORTED_PROVIDERS
    assert isinstance(SUPPORTED_PROVIDERS, list), \
        "SUPPORTED_PROVIDERS must be a list"
    assert SUPPORTED_PROVIDERS, "SUPPORTED_PROVIDERS must be non-empty"
    # v120 Phase 3: anthropic + portkey activated alongside the original
    # three.  This is the live provider set.
    assert set(SUPPORTED_PROVIDERS) == {
        "google", "openai", "ollama", "anthropic", "portkey"}, \
        "unexpected provider set: %r" % (SUPPORTED_PROVIDERS,)


def test_graph_nodes_shares_core_list():
    """graph_nodes must expose the identical value imported from core.llm."""
    from core.llm import SUPPORTED_PROVIDERS as core_list
    try:
        from agent.graph_nodes import SUPPORTED_PROVIDERS as gn_list
    except Exception as exc:  # libtbx (or other build dep) absent in sandbox
        if "libtbx" in str(exc) or isinstance(exc, ImportError):
            print("  (skip) graph_nodes import needs Phenix build: %s"
                  % type(exc).__name__)
            return
        raise
    assert gn_list == core_list, \
        "graph_nodes.SUPPORTED_PROVIDERS must equal core.llm's"
    # NOTE: we intentionally assert value-equality (==), not object identity
    # (is).  In a real Phenix build graph_nodes imports via
    # `libtbx.langchain.core.llm` while this test imports via `core.llm`;
    # Python treats those as two distinct module objects, each with its own
    # list instance.  Equal value is the contract that matters (same providers
    # everywhere); object identity is an implementation detail that does not
    # hold across the two import roots.


# ============================================================================
# Source-scan drift guards (always run)
# ============================================================================

def test_graph_nodes_imports_from_core_not_literal():
    assert os.path.isfile(_GRAPH_NODES), "graph_nodes.py not found"
    with open(_GRAPH_NODES) as fh:
        src = fh.read()
    assert "import SUPPORTED_PROVIDERS" in src, \
        "graph_nodes must import SUPPORTED_PROVIDERS"
    assert "core.llm import SUPPORTED_PROVIDERS" in src, \
        "graph_nodes must import SUPPORTED_PROVIDERS from core.llm"
    # No local literal assignment (the thing we removed).
    assert not re.search(r"SUPPORTED_PROVIDERS\s*=\s*\[", src), \
        "graph_nodes must NOT define a literal SUPPORTED_PROVIDERS list"


def test_ai_agent_imports_from_core_not_literal():
    path = _find_ai_agent_source()
    if path is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    with open(path) as fh:
        src = fh.read()
    assert "core.llm import SUPPORTED_PROVIDERS" in src, \
        "ai_agent must import SUPPORTED_PROVIDERS from core.llm"
    # The old hardcoded fallback literal must be gone.
    assert not re.search(r"SUPPORTED_PROVIDERS\s*=\s*\[", src), \
        "ai_agent must NOT define a literal SUPPORTED_PROVIDERS fallback"
    # And it must no longer import the list from graph_nodes (that route kept
    # graph_nodes as an unnecessary intermediary; core.llm is the source).
    assert "graph_nodes import SUPPORTED_PROVIDERS" not in src, \
        "ai_agent should import SUPPORTED_PROVIDERS from core.llm, " \
        "not via graph_nodes"


def test_core_llm_is_sole_literal_definition():
    """Exactly one literal `SUPPORTED_PROVIDERS = [...]` should exist in the
    shipped modules, and it must live in core/llm.py."""
    core_path = os.path.join(_ROOT, "core", "llm.py")
    assert os.path.isfile(core_path), "core/llm.py not found"
    with open(core_path) as fh:
        core_src = fh.read()
    assert re.search(r"SUPPORTED_PROVIDERS\s*=\s*\[", core_src), \
        "core/llm.py must hold the literal SUPPORTED_PROVIDERS definition"


# ============================================================================
# Runner
# ============================================================================

_TESTS = [
    test_core_llm_defines_supported_providers,
    test_graph_nodes_shares_core_list,
    test_graph_nodes_imports_from_core_not_literal,
    test_ai_agent_imports_from_core_not_literal,
    test_core_llm_is_sole_literal_definition,
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
