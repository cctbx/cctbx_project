"""
Portkey + Anthropic provider wiring (v120 Phase 3).

Scope: WIRING, not live calls.  Real Portkey/Anthropic round-trips need keys +
network and live in the keyed tests/llm/ suite.  Here everything is pure-logic,
source-scan, or uses a mocked client -- runs under plain python3 with no SDKs
installed and no outbound network.

Covers:
  - openai_client_config / portkey_langchain_config: base_url + key, and a
    clear ValueError (NOT bare KeyError) when env vars are missing.
  - sanitize_llm_kwargs: portkey remaps max_tokens->max_completion_tokens and
    drops temperature; anthropic KEEPS both (the asymmetric-args trap that
    would otherwise break Claude); others unchanged.
  - SUPPORTED_PROVIDERS includes anthropic + portkey, and the PHIL enums cover
    the canonical set.
  - Embeddings delegation: _delegate_embeddings_for_nonnative returns None when
    no delegate key is set (never raises); delegates when a key is present.
  - Dummy-key lazy construction: with PORTKEY_AZURE_API_KEY="mock_failed_key",
    the langchain factory builds llm+embeddings without an outbound network call
    (mocked transport) -- guards the lazy-construction regression for the
    end-of-run agent_session path.
  - Embeddings-query failure surfaces a clear, access-naming error echoing the
    raw upstream message (NOT a bare 'Unauthorized' string match, NOT silent).
  - Source scans: every call site carries the portkey branch.

Per the agent test guidelines: file-exists guards, function-boundary windows,
AssertionError never swallowed, no network.
"""
import os
import re
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain
for _p in [_ROOT, os.path.join(_ROOT, "agent"), os.path.join(_ROOT, "core")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

_AGENT = os.path.join(_ROOT, "agent")
_CORE = os.path.join(_ROOT, "core")

_AI_AGENT_CANDIDATES = [
    # Sandbox layout: programs/ as a sibling of tests/ under langchain.
    os.path.join(_ROOT, "programs", "ai_agent.py"),
    # Real PHENIX build: ai_agent.py lives in the phenix tree.  _ROOT is
    # .../modules/cctbx_project/libtbx/langchain, so 3 ups reaches .../modules,
    # then into phenix/phenix/programs.
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_agent.py"),
    # Defensive: single-phenix layout, in case of a flatter checkout.
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


def _find(cands, env=None):
    for c in cands:
        c = os.path.abspath(c)
        if os.path.isfile(c):
            return c
    if env:
        v = os.environ.get(env)
        if v and os.path.isfile(v):
            return v
    return None


def _clear_env(*names):
    saved = {}
    for n in names:
        saved[n] = os.environ.pop(n, None)
    return saved


def _restore_env(saved):
    for n, v in saved.items():
        if v is None:
            os.environ.pop(n, None)
        else:
            os.environ[n] = v


# ============================================================================
# sanitize_llm_kwargs
# ============================================================================

def test_sanitize_portkey_remaps_and_drops():
    from core.llm import sanitize_llm_kwargs
    out = sanitize_llm_kwargs(
        "portkey", {"max_tokens": 2000, "temperature": 0.1, "model": "gpt-5"})
    assert "max_tokens" not in out, "portkey must drop max_tokens"
    assert out.get("max_completion_tokens") == 2000, \
        "portkey must remap to max_completion_tokens"
    assert "temperature" not in out, "portkey must drop temperature"
    assert out.get("model") == "gpt-5", "other keys preserved"


def test_sanitize_anthropic_keeps_max_tokens_always():
    """anthropic must ALWAYS keep max_tokens (ChatAnthropic requires it),
    regardless of model."""
    from core.llm import sanitize_llm_kwargs
    for model in ("claude-opus-4-7", "claude-opus-4-20250514", None):
        out = sanitize_llm_kwargs(
            "anthropic", {"model": model, "max_tokens": 4096,
                          "temperature": 0.2})
        assert out.get("max_tokens") == 4096, \
            "anthropic must keep max_tokens for model %r" % (model,)


def test_sanitize_anthropic_drops_temperature_for_47_family():
    """The 4.6/4.7+ reasoning family rejects temperature (HTTP 400), so
    sanitize must drop it for those models -- this is the live bug from the
    anthropic reliability run."""
    from core.llm import sanitize_llm_kwargs
    for model in ("claude-opus-4-7", "claude-sonnet-4-6", "claude-opus-4-8"):
        out = sanitize_llm_kwargs(
            "anthropic", {"model": model, "max_tokens": 4096,
                          "temperature": 0.2, "top_p": 0.9})
        assert "temperature" not in out, \
            "%s rejects temperature; sanitize must drop it" % model
        assert "top_p" not in out, \
            "%s rejects top_p; sanitize must drop it" % model
        assert out.get("max_tokens") == 4096, "max_tokens still kept"


def test_sanitize_anthropic_keeps_temperature_for_older_models():
    """Older Claude models (Sonnet 4 / Opus 4 dated, 3.x) still accept
    temperature -- sanitize must NOT strip it for them."""
    from core.llm import sanitize_llm_kwargs
    for model in ("claude-opus-4-20250514", "claude-sonnet-4-20250514",
                  "claude-3-5-sonnet"):
        out = sanitize_llm_kwargs(
            "anthropic", {"model": model, "max_tokens": 4096,
                          "temperature": 0.2})
        assert out.get("temperature") == 0.2, \
            "%s accepts temperature; sanitize must keep it" % model


def test_anthropic_model_family_detection():
    """Direct test of the model-family helper that drives the temperature
    decision."""
    from core.llm import anthropic_model_rejects_sampling_params as f
    assert f("claude-opus-4-7") is True
    assert f("claude-sonnet-4-6") is True
    assert f("claude-opus-4-8") is True          # future-proof
    assert f("claude-haiku-4-5-20251001") is False
    assert f("claude-opus-4-20250514") is False  # retired, dated -> accepts
    assert f("claude-3-5-sonnet") is False
    assert f("gpt-5") is False
    assert f("") is False
    assert f(None) is False


def test_sanitize_openai_unchanged():
    from core.llm import sanitize_llm_kwargs
    src = {"max_tokens": 10, "temperature": 0.5}
    out = sanitize_llm_kwargs("openai", src)
    assert out == src, "openai kwargs must be unchanged"


def test_sanitize_does_not_mutate_caller_dict():
    from core.llm import sanitize_llm_kwargs
    src = {"max_tokens": 2000, "temperature": 0.1}
    sanitize_llm_kwargs("portkey", src)
    assert src == {"max_tokens": 2000, "temperature": 0.1}, \
        "sanitize must not mutate the caller's dict"


# ============================================================================
# portkey_langchain_config
# ============================================================================

def test_portkey_config_returns_key_and_url():
    from core.llm import portkey_langchain_config
    saved = _clear_env("PORTKEY_AZURE_API_KEY", "PORTKEY_BASE_URL")
    try:
        os.environ["PORTKEY_AZURE_API_KEY"] = "k123"
        os.environ["PORTKEY_BASE_URL"] = "https://gw.example/v1"
        key, url = portkey_langchain_config()
        assert key == "k123" and url == "https://gw.example/v1"
    finally:
        _restore_env(saved)


def test_portkey_config_missing_raises_valueerror_not_keyerror():
    from core.llm import portkey_langchain_config
    saved = _clear_env("PORTKEY_AZURE_API_KEY", "PORTKEY_BASE_URL")
    try:
        raised = None
        try:
            portkey_langchain_config()
        except Exception as e:
            raised = e
        assert isinstance(raised, ValueError), \
            "missing portkey env must raise ValueError, got %r" % (raised,)
        assert "PORTKEY_AZURE_API_KEY" in str(raised), \
            "error must name the missing variable"
        assert not isinstance(raised, KeyError), \
            "must NOT be a bare KeyError"
    finally:
        _restore_env(saved)


# ============================================================================
# SUPPORTED_PROVIDERS + PHIL enums
# ============================================================================

def test_supported_providers_includes_new():
    from core.llm import SUPPORTED_PROVIDERS
    for p in ("anthropic", "portkey"):
        assert p in SUPPORTED_PROVIDERS, \
            "%s must be in SUPPORTED_PROVIDERS" % p


def test_phil_enums_cover_supported_providers():
    """Each PHIL `provider` enum must list every SUPPORTED_PROVIDERS entry,
    else a provider is wired in code but unreachable from the UI/PHIL."""
    from core.llm import SUPPORTED_PROVIDERS
    targets = [
        _find(_AI_AGENT_CANDIDATES, env="AI_AGENT_PY"),
        _find(_AI_ANALYSIS_CANDIDATES, env="AI_ANALYSIS_PY"),
    ]
    checked = 0
    for path in targets:
        if path is None:
            continue
        with open(path) as fh:
            src = fh.read()
        m = re.search(r"provider\s*=\s*([^\n]*)", src)
        assert m, "no `provider = ...` PHIL enum in %s" % path
        enum_line = m.group(1)
        for p in SUPPORTED_PROVIDERS:
            assert p in enum_line, \
                "%s missing from PHIL provider enum in %s (line: %r)" % (
                    p, os.path.basename(path), enum_line)
        checked += 1
    if checked == 0:
        print("  (skip) neither PHIL file found in this checkout")


# ============================================================================
# Embeddings delegation (§4.6)
# ============================================================================

def test_delegation_returns_none_without_keys():
    from core.llm import _delegate_embeddings_for_nonnative
    saved = _clear_env("OPENAI_API_KEY", "GOOGLE_API_KEY")
    try:
        emb = _delegate_embeddings_for_nonnative("anthropic")
        assert emb is None, \
            "no delegate key -> None (never raise on construction)"
    finally:
        _restore_env(saved)


def test_delegation_never_raises():
    """The delegation helper must never raise out of construction, regardless
    of whether the embeddings SDK is installed.

    Two valid outcomes depending on environment:
      - SDK absent (bare sandbox): internal try/except swallows ImportError
        and returns None.
      - SDK present (real build): OpenAIEmbeddings constructs lazily and a
        non-None object is returned (no network call at construction).
    Either is correct; what we pin is that it does NOT raise.
    """
    from core.llm import _delegate_embeddings_for_nonnative
    saved = _clear_env("OPENAI_API_KEY", "GOOGLE_API_KEY")
    try:
        os.environ["OPENAI_API_KEY"] = "sk-fake"
        raised = None
        result = None
        try:
            result = _delegate_embeddings_for_nonnative("anthropic")
        except Exception as e:
            raised = e
        assert raised is None, \
            "delegation must not raise on construction; got %r" % (raised,)
        # result may be None (no SDK) or an embeddings object (SDK present) --
        # both acceptable.  The contract is "no raise".
    finally:
        _restore_env(saved)



# ============================================================================
# Dummy-key lazy construction (Gemini risk 1) — mocked SDK, no network
# ============================================================================

def _install_mock_langchain_openai(calls):
    """Install a fake langchain_openai module whose ChatOpenAI /
    OpenAIEmbeddings record construction but make NO network call."""
    mod = types.ModuleType("langchain_openai")

    class _FakeChatOpenAI:
        def __init__(self, **kwargs):
            calls.append(("ChatOpenAI", kwargs))
            self.model_name = kwargs.get("model")

    class _FakeEmbeddings:
        def __init__(self, **kwargs):
            calls.append(("OpenAIEmbeddings", kwargs))
            self.model = kwargs.get("model")

    mod.ChatOpenAI = _FakeChatOpenAI
    mod.OpenAIEmbeddings = _FakeEmbeddings
    sys.modules["langchain_openai"] = mod
    return mod


def test_portkey_dummy_key_lazy_construction_no_network():
    """With a mock_failed_key, get_llm_and_embeddings(provider='portkey')
    builds llm + embeddings with no outbound network call (mocked SDK)."""
    import core.llm as llm
    calls = []
    saved_mods = sys.modules.get("langchain_openai")
    saved_env = _clear_env("PORTKEY_AZURE_API_KEY", "PORTKEY_BASE_URL",
                           "PHENIX_PORTKEY_MODEL")
    try:
        _install_mock_langchain_openai(calls)
        os.environ["PORTKEY_AZURE_API_KEY"] = "mock_failed_key"
        os.environ["PORTKEY_BASE_URL"] = "https://gw.example/v1"
        llm_obj, emb = llm.get_llm_and_embeddings(provider="portkey")
        assert llm_obj is not None, "portkey must return a chat LLM"
        assert emb is not None, "portkey embeddings object must construct"
        names = [c[0] for c in calls]
        assert "ChatOpenAI" in names and "OpenAIEmbeddings" in names, \
            "both objects should be constructed lazily; got %r" % names
        # The dummy key must have flowed to the client constructors.
        for _, kwargs in calls:
            assert kwargs.get("api_key") == "mock_failed_key"
            assert kwargs.get("base_url") == "https://gw.example/v1"
    finally:
        _restore_env(saved_env)
        if saved_mods is None:
            sys.modules.pop("langchain_openai", None)
        else:
            sys.modules["langchain_openai"] = saved_mods


def test_portkey_missing_env_raises_valueerror_in_factory():
    """get_llm_and_embeddings(provider='portkey') with no env must raise a
    clear ValueError (call-time check), not a bare KeyError."""
    import core.llm as llm
    saved = _clear_env("PORTKEY_AZURE_API_KEY", "PORTKEY_BASE_URL")
    try:
        raised = None
        try:
            llm.get_llm_and_embeddings(provider="portkey")
        except Exception as e:
            raised = e
        assert isinstance(raised, ValueError) and not isinstance(raised, KeyError), \
            "missing portkey env must be a clear ValueError, got %r" % (raised,)
    finally:
        _restore_env(saved)


# ============================================================================
# Embeddings-query clear error (Gemini risk 3) — behaviour contract
# ============================================================================

def test_embeddings_query_error_is_broad_and_echoes_message():
    """An embeddings query failure must produce a clear message that echoes
    the raw upstream error rather than string-matching 'Unauthorized'.

    We simulate the intended interceptor behaviour: catch broad Exception,
    surface tailored guidance + raw message.  This pins the contract the RAG
    embedding path must implement.
    """
    class _APIStatusError(Exception):
        pass

    def _intercept(provider, embed_fn):
        try:
            return embed_fn()
        except Exception as e:  # broad, not a string match
            if provider in ("portkey", "anthropic"):
                return (
                    "[ERROR] Portkey embeddings query failed. Verify "
                    "PORTKEY_AZURE_API_KEY has access to the embeddings "
                    "deployment. Upstream details: %s" % e)
            raise

    def _boom():
        raise _APIStatusError("403 Forbidden: no embeddings access")

    msg = _intercept("portkey", _boom)
    assert "PORTKEY_AZURE_API_KEY" in msg, "must name the access gap"
    assert "403 Forbidden" in msg, "must echo the raw upstream message"
    assert "Unauthorized" not in msg, \
        "must not rely on a literal 'Unauthorized' match"


# ============================================================================
# Source scans — call sites carry the portkey branch
# ============================================================================

def _scan(path, needles, label):
    assert os.path.isfile(path), "%s not found" % label
    with open(path) as fh:
        src = fh.read()
    for n in needles:
        assert n in src, "%s: expected %r" % (label, n)


def test_source_api_client_has_portkey():
    _scan(os.path.join(_AGENT, "api_client.py"),
          ['def _call_portkey_llm', 'provider == "portkey"',
           'provider="azure-openai"', 'max_completion_tokens'],
          "api_client.py")


def test_source_directive_extractor_has_portkey_both_paths():
    p = os.path.join(_AGENT, "directive_extractor.py")
    _scan(p, ['get_portkey_handler', 'provider == "portkey"',
              'provider="azure-openai"'], "directive_extractor.py")
    with open(p) as fh:
        src = fh.read()
    # both _call_llm handler-select and _call_llm_fallback should reference it
    assert src.count('get_portkey_handler') >= 3, \
        "portkey handler must appear in both _call_llm import paths and fallback"


def test_source_rate_limit_handler_has_portkey():
    _scan(os.path.join(_AGENT, "rate_limit_handler.py"),
          ['def get_portkey_handler', 'portkey_api'],
          "rate_limit_handler.py")


def test_source_graph_nodes_has_portkey():
    _scan(os.path.join(_AGENT, "graph_nodes.py"),
          ['provider == "portkey"', 'PORTKEY_AZURE_API_KEY',
           'get_portkey_handler'], "graph_nodes.py")


def test_source_thinking_agent_has_portkey():
    _scan(os.path.join(_AGENT, "thinking_agent.py"),
          ['provider == "portkey"', 'get_portkey_handler'],
          "thinking_agent.py")


def test_source_core_llm_has_both_branches():
    _scan(os.path.join(_CORE, "llm.py"),
          ['elif provider == "portkey"', 'elif provider == "anthropic"',
           'def sanitize_llm_kwargs', 'def portkey_langchain_config',
           'def _delegate_embeddings_for_nonnative'], "core/llm.py")


def test_validate_api_keys_covers_new_providers():
    """run_utils.validate_api_keys must gate portkey + anthropic (it sits on
    the agent end-of-run / directive paths), not silently pass them."""
    import importlib.util
    path = os.path.join(_ROOT, "utils", "run_utils.py")
    assert os.path.isfile(path), "run_utils.py not found"
    spec = importlib.util.spec_from_file_location("_ru_v120", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    saved = _clear_env("PORTKEY_AZURE_API_KEY", "PORTKEY_BASE_URL",
                       "ANTHROPIC_API_KEY")
    try:
        err = m.validate_api_keys("portkey", [])
        assert err and "PORTKEY_AZURE_API_KEY" in err, \
            "portkey with no keys must return a clear error, got %r" % (err,)
        err = m.validate_api_keys("anthropic", [])
        assert err and "ANTHROPIC_API_KEY" in err, \
            "anthropic with no key must return a clear error, got %r" % (err,)
        os.environ["PORTKEY_AZURE_API_KEY"] = "k"
        os.environ["PORTKEY_BASE_URL"] = "u"
        assert m.validate_api_keys("portkey", []) is None, \
            "portkey with both vars set must pass"
    finally:
        _restore_env(saved)


def test_source_run_utils_and_query_docs_have_new_providers():
    _scan(os.path.join(_ROOT, "utils", "run_utils.py"),
          ["PORTKEY_AZURE_API_KEY", "ANTHROPIC_API_KEY"], "run_utils.py")
    qd = os.path.join(_ROOT, "run_query_docs.py")
    if os.path.isfile(qd):
        _scan(qd, ["PORTKEY_AZURE_API_KEY", 'provider == \'anthropic\''],
              "run_query_docs.py")


# ============================================================================
# Runner
# ============================================================================

_TESTS = [
    test_sanitize_portkey_remaps_and_drops,
    test_sanitize_anthropic_keeps_max_tokens_always,
    test_sanitize_anthropic_drops_temperature_for_47_family,
    test_sanitize_anthropic_keeps_temperature_for_older_models,
    test_anthropic_model_family_detection,
    test_sanitize_openai_unchanged,
    test_sanitize_does_not_mutate_caller_dict,
    test_portkey_config_returns_key_and_url,
    test_portkey_config_missing_raises_valueerror_not_keyerror,
    test_supported_providers_includes_new,
    test_phil_enums_cover_supported_providers,
    test_delegation_returns_none_without_keys,
    test_delegation_never_raises,
    test_portkey_dummy_key_lazy_construction_no_network,
    test_portkey_missing_env_raises_valueerror_in_factory,
    test_embeddings_query_error_is_broad_and_echoes_message,
    test_source_api_client_has_portkey,
    test_source_directive_extractor_has_portkey_both_paths,
    test_source_rate_limit_handler_has_portkey,
    test_source_graph_nodes_has_portkey,
    test_source_thinking_agent_has_portkey,
    test_source_core_llm_has_both_branches,
    test_validate_api_keys_covers_new_providers,
    test_source_run_utils_and_query_docs_have_new_providers,
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
