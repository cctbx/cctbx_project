"""K_H1: Default-model centralization regression tests.

v119.H1.  Source of truth for default model strings is core/llm.py.
These tests assert:

- All five role-tables are populated for expected providers
- No DEFAULTS value appears in RETIRED_MODELS
- default_model_for_provider() resolves correctly, normalizes
  provider keys, and raises on unknown inputs
- api_client.py, directive_extractor.py, and test_api_keys.py read
  from the helper (per-callsite source scans)
- AST-level: no runtime-code string-literal assignment in agent/
  modules contains a provider-model pattern (catches future
  contributors adding new hardcoded defaults)
- Retired-model detection helper recognises the right signatures
- [DIRECTIVE_EXTRACTION_MODEL_RETIRED] marker writes the expected
  fields to stderr
- Baseline fingerprint: the exact table values at H1 ship time
  (relax intentionally when a default changes)
"""
from __future__ import absolute_import, division, print_function

import os


# ---------- Import / file-discovery helpers ------------------------
#
# Two paths to find a file:
#   1. _module_file(name) tries `from libtbx.langchain.X import Y` and
#      falls back to `from X import Y`.  Returns the imported module
#      on success (callers needing actual function references use this).
#   2. _source_path(rel_path) finds the file on disk via known paths
#      under the langchain/ tree without importing.  Source-scan
#      tests use this so they work even when the module's transitive
#      imports aren't installed (e.g. running from a partial sandbox).


def _import_core_llm():
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm
    return llm


def _import_directive_extractor_for_helpers():
    """Import directive_extractor to call _looks_like_retired_model_error
    and _emit_retired_model_marker.  May fail if transitive deps are
    not installed; tests that use this should handle that gracefully.
    """
    try:
        from libtbx.langchain.agent import directive_extractor
    except ImportError:
        from agent import directive_extractor
    return directive_extractor


def _source_path(rel_path):
    """Return absolute path to a file under the langchain/ tree.

    rel_path is relative to langchain/, e.g. "agent/api_client.py".
    Tries the test's parent directory (sibling-of-tests layout)
    first, then walks up looking for the file.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "..", rel_path),                # tests/.. = langchain/
        os.path.join(here, "..", "..", rel_path),          # one more up if nested
    ]
    for c in candidates:
        c = os.path.normpath(c)
        if os.path.isfile(c):
            return c
    return None


def _source_text(rel_path):
    """Return contents of a file under langchain/, or None if missing."""
    path = _source_path(rel_path)
    if path is None:
        return None
    with open(path) as f:
        return f.read()


def _function_body(src, fn_name):
    """Slice src to the body of `def <fn_name>` (up to next def)."""
    idx = src.find("def %s" % fn_name)
    assert idx > 0, "%s not found in source" % fn_name
    next_def = src.find("\ndef ", idx + 1)
    if next_def < 0:
        return src[idx:]
    return src[idx:next_def]


# ---------- Table-shape tests --------------------------------------

def test_decision_models_table_complete():
    llm = _import_core_llm()
    expected = {"google", "openai", "anthropic", "ollama"}
    assert expected.issubset(set(llm.DECISION_MODEL_DEFAULTS)), (
      "DECISION_MODEL_DEFAULTS missing keys: %s" % (
        expected - set(llm.DECISION_MODEL_DEFAULTS)))
    print("  PASS: test_decision_models_table_complete")


def test_rag_models_table_minimum():
    """RAG path supports google/openai/ollama (no anthropic branch)."""
    llm = _import_core_llm()
    expected = {"google", "openai", "ollama"}
    assert expected.issubset(set(llm.RAG_MODEL_DEFAULTS))
    print("  PASS: test_rag_models_table_minimum")


def test_rag_embedding_table_minimum():
    llm = _import_core_llm()
    expected = {"google", "openai", "ollama"}
    assert expected.issubset(set(llm.RAG_EMBEDDING_DEFAULTS))
    print("  PASS: test_rag_embedding_table_minimum")


def test_expensive_models_table_minimum():
    llm = _import_core_llm()
    expected = {"google", "openai", "ollama"}
    assert expected.issubset(set(llm.EXPENSIVE_MODEL_DEFAULTS))
    print("  PASS: test_expensive_models_table_minimum")


def test_cheap_models_table_minimum():
    """CHEAP only requires ollama; remote providers fall through."""
    llm = _import_core_llm()
    assert "ollama" in llm.CHEAP_MODEL_DEFAULTS
    print("  PASS: test_cheap_models_table_minimum")


def test_no_default_in_retired_list():
    """No DEFAULTS value appears in RETIRED_MODELS."""
    llm = _import_core_llm()
    all_tables = [
        ("DECISION_MODEL_DEFAULTS", llm.DECISION_MODEL_DEFAULTS),
        ("RAG_MODEL_DEFAULTS", llm.RAG_MODEL_DEFAULTS),
        ("RAG_EMBEDDING_DEFAULTS", llm.RAG_EMBEDDING_DEFAULTS),
        ("EXPENSIVE_MODEL_DEFAULTS", llm.EXPENSIVE_MODEL_DEFAULTS),
        ("CHEAP_MODEL_DEFAULTS", llm.CHEAP_MODEL_DEFAULTS),
    ]
    for table_name, table in all_tables:
        for provider, model in table.items():
            assert model not in llm.RETIRED_MODELS, (
              "%s[%r]=%r is on RETIRED_MODELS list -- "
              "update to a current model." % (
                table_name, provider, model))
    print("  PASS: test_no_default_in_retired_list")


# ---------- Helper-function tests ----------------------------------

def test_default_model_for_provider_basic():
    """Default role is 'decision'."""
    llm = _import_core_llm()
    for provider, expected in llm.DECISION_MODEL_DEFAULTS.items():
        got = llm.default_model_for_provider(provider)
        assert got == expected, "%s != %s" % (got, expected)
    print("  PASS: test_default_model_for_provider_basic")


def test_default_model_for_provider_role():
    """All five roles resolve correctly."""
    llm = _import_core_llm()
    for role_name, table_attr in [
        ("decision", "DECISION_MODEL_DEFAULTS"),
        ("rag", "RAG_MODEL_DEFAULTS"),
        ("rag_embedding", "RAG_EMBEDDING_DEFAULTS"),
        ("expensive", "EXPENSIVE_MODEL_DEFAULTS"),
        ("cheap", "CHEAP_MODEL_DEFAULTS"),
    ]:
        table = getattr(llm, table_attr)
        for provider, expected in table.items():
            got = llm.default_model_for_provider(
                provider, role=role_name)
            assert got == expected, (
              "role=%s provider=%s: got %r expected %r" % (
                role_name, provider, got, expected))
    print("  PASS: test_default_model_for_provider_role")


def test_default_model_for_provider_normalizes_keys():
    """Provider key is case- and whitespace-insensitive."""
    llm = _import_core_llm()
    expected = llm.DECISION_MODEL_DEFAULTS["google"]
    for raw in ("Google", " GOOGLE ", "google\t", "  google"):
        assert llm.default_model_for_provider(raw) == expected, (
          "Provider normalization failed for %r" % raw)
    print("  PASS: test_default_model_for_provider_normalizes_keys")


def test_default_model_for_provider_unknown_role_raises():
    llm = _import_core_llm()
    raised = False
    try:
        llm.default_model_for_provider("google", role="bogus")
    except ValueError:
        raised = True
    assert raised
    print("  PASS: test_default_model_for_provider_unknown_role_raises")


def test_default_model_for_provider_unknown_provider_raises():
    llm = _import_core_llm()
    raised = False
    try:
        llm.default_model_for_provider("not-a-provider")
    except ValueError:
        raised = True
    assert raised
    print("  PASS: test_default_model_for_provider_unknown_provider_raises")


# ---------- Per-callsite source-scan tests -------------------------

def test_api_client_google_uses_central_default():
    src = _source_text("agent/api_client.py")
    assert src is not None, "agent/api_client.py not found"
    body = _function_body(src, "_call_google_llm")
    assert "default_model_for_provider" in body, (
      "_call_google_llm does not call default_model_for_provider")
    assert 'or "gemini-' not in body, (
      "_call_google_llm still has hardcoded gemini default")
    print("  PASS: test_api_client_google_uses_central_default")


def test_api_client_openai_uses_central_default():
    src = _source_text("agent/api_client.py")
    assert src is not None, "agent/api_client.py not found"
    body = _function_body(src, "_call_openai_llm")
    assert "default_model_for_provider" in body
    assert 'or "gpt-' not in body
    print("  PASS: test_api_client_openai_uses_central_default")


def test_api_client_anthropic_uses_central_default():
    src = _source_text("agent/api_client.py")
    assert src is not None, "agent/api_client.py not found"
    body = _function_body(src, "_call_anthropic_llm")
    assert "default_model_for_provider" in body
    assert 'or "claude-' not in body
    print("  PASS: test_api_client_anthropic_uses_central_default")


def test_api_client_ollama_uses_central_default():
    """v119.H13: ollama site now uses resolve_model_for_provider
    (which wraps default_model_for_provider with env-var precedence).
    The test verifies a centralized resolver is used — either function
    name — and that no hardcoded ollama model strings exist.
    """
    src = _source_text("agent/api_client.py")
    assert src is not None, "agent/api_client.py not found"
    body = _function_body(src, "_call_ollama_llm")
    assert ("default_model_for_provider" in body
            or "resolve_model_for_provider" in body), (
        "_call_ollama_llm should use a centralized model resolver")
    assert 'or "llama' not in body, (
        "_call_ollama_llm still has hardcoded llama default")
    print("  PASS: test_api_client_ollama_uses_central_default")


def test_directive_extractor_fallback_uses_central_defaults():
    """Single test for the fallback function -- covers all 4 providers.

    v119.H13: the ollama branch was switched to resolve_model_for_provider
    (env-var precedence wrapper).  Count BOTH function names — the
    centralization invariant is "every provider resolves via a central
    helper", not "every provider uses the exact same helper name".
    """
    src = _source_text("agent/directive_extractor.py")
    assert src is not None, "agent/directive_extractor.py not found"
    body = _function_body(src, "_call_llm_fallback")
    n = (body.count("default_model_for_provider")
         + body.count("resolve_model_for_provider"))
    assert n >= 4, (
      "_call_llm_fallback should resolve default for all 4 providers; "
      "found only %d call(s) (counting both default_model_for_provider "
      "and resolve_model_for_provider)" % n)
    for pattern in ('or "gemini-', 'or "gpt-',
                    'or "claude-', 'or "llama'):
        assert pattern not in body, (
          "_call_llm_fallback still has hardcoded default: %s"
          % pattern)
    print("  PASS: test_directive_extractor_fallback_uses_central_defaults")


def test_test_api_keys_uses_central_defaults():
    """test_api_keys.py reads from default_model_for_provider."""
    src = _source_text("tests/test_api_keys.py")
    if src is None:
        print("  SKIP: test_api_keys.py not found")
        return
    assert "default_model_for_provider" in src, (
      "test_api_keys.py doesn't import default_model_for_provider")
    # Should not have hardcoded model literals as string values in
    # API request bodies or comparison strings.
    import re
    hardcoded = re.findall(
      r'(?:model\s*=|"model":\s*)"(gemini-|gpt-|claude-|llama)[^"]*"',
      src)
    assert not hardcoded, (
      "test_api_keys.py still has hardcoded model strings: %s"
      % hardcoded)
    print("  PASS: test_test_api_keys_uses_central_defaults")


# ---------- AST-level structural regression test -------------------

def test_no_orphan_model_strings_in_agent_runtime_code():
    """AST scan: no runtime model-string assignments outside core/llm.py.

    Walks the AST of agent/api_client.py and agent/directive_extractor.py
    looking for any string-literal assignment whose value matches a
    provider-model pattern.  Module docstrings (which contain the
    illustrative example in directive_extractor.py) are exempt via
    ast.get_docstring(); other docstrings and comments are exempt
    because they aren't ast.Constant nodes in assignment contexts.

    If a future contributor adds `model_name = "gemini-3.0-pro"` to
    api_client.py, this test fires.
    """
    import ast
    import re

    model_pattern = re.compile(
      r"^(gemini-|gpt-|claude-|llama[0-9]|qwen|nomic-embed-|"
      r"text-embedding-)")

    files = [
        "agent/api_client.py",
        "agent/directive_extractor.py",
    ]

    for rel_path in files:
        path = _source_path(rel_path)
        assert path is not None, "%s not found" % rel_path
        with open(path) as f:
            source = f.read()
        tree = ast.parse(source)
        module_docstring = ast.get_docstring(tree)
        offenders = []
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Constant) and isinstance(
                    node.value, str)):
                continue
            if not model_pattern.match(node.value):
                continue
            # Exclude the module docstring (which contains the
            # illustrative example in directive_extractor.py).
            if (module_docstring is not None
                    and node.value in module_docstring):
                continue
            offenders.append((
              node.lineno, node.col_offset, node.value))
        assert not offenders, (
          "Orphan model strings in %s (line, col, value):\n%s\n"
          "All default model strings must live in core/llm.py."
          % (path,
             "\n".join("  %d:%d  %r" % o for o in offenders)))
    print("  PASS: test_no_orphan_model_strings_in_agent_runtime_code")


# ---------- Retired-model helper tests ----------------------------

def _try_import_directive_extractor_for_helpers():
    """Try to import directive_extractor; return (module, None) on
    success or (None, error_message) on failure.

    The helper-execution tests below skip with a clear message when
    the module's transitive dependencies aren't available (sandbox
    environments).  Under PHENIX these always succeed.
    """
    try:
        from libtbx.langchain.agent import directive_extractor
        return directive_extractor, None
    except ImportError:
        pass
    try:
        from agent import directive_extractor
        return directive_extractor, None
    except ImportError as e:
        return None, str(e)


def test_retired_error_detection_positive_cases():
    de, err = _try_import_directive_extractor_for_helpers()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    f = de._looks_like_retired_model_error
    cases = [
      ("404 Client Error: NOT_FOUND for model 'gemini-2.0-flash'",
       "gemini-2.0-flash"),
      ("Model gemini-2.0-flash is no longer available",
       "gemini-2.0-flash"),
      ("Error: not found for model fake-name",
       "fake-name"),
      ("404 NOT_FOUND: the model has been deprecated",
       "any-model"),
    ]
    for err_str, model in cases:
        assert f(err_str, model), "Should match: %r / %r" % (err_str, model)
    print("  PASS: test_retired_error_detection_positive_cases")


def test_retired_error_detection_negative_cases():
    de, err = _try_import_directive_extractor_for_helpers()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    f = de._looks_like_retired_model_error
    cases = [
      ("Rate limit exceeded", "gemini-2.5-flash-lite"),
      ("Connection timed out", "gemini-2.5-flash-lite"),
      ("", "gemini-2.5-flash-lite"),
      ("500 Internal Server Error", "gemini-2.5-flash-lite"),
      ("Invalid API key", "gemini-2.5-flash-lite"),
    ]
    for err_str, model in cases:
        assert not f(err_str, model), (
          "Should NOT match: %r / %r" % (err_str, model))
    print("  PASS: test_retired_error_detection_negative_cases")


def test_retired_marker_writes_expected_fields():
    import io
    import sys
    de, err = _try_import_directive_extractor_for_helpers()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    captured = io.StringIO()
    old_stderr = sys.stderr
    sys.stderr = captured
    try:
        de._emit_retired_model_marker(
            "google", "fake-retired-model", log=lambda m: None)
    finally:
        sys.stderr = old_stderr
    text = captured.getvalue()
    assert "[DIRECTIVE_EXTRACTION_MODEL_RETIRED]" in text
    assert "fake-retired-model" in text
    assert "core/llm.py" in text
    # Marker must mention at least one of the tables so the operator
    # can find the right one to edit
    assert any(t in text for t in (
        "DECISION_MODEL_DEFAULTS",
        "RAG_MODEL_DEFAULTS",
        "EXPENSIVE_MODEL_DEFAULTS"))
    print("  PASS: test_retired_marker_writes_expected_fields")


def test_fallback_property_documented():
    """Single-call property is documented in _call_llm_fallback docstring.

    Done via source scan (no module import) so this test works in
    sandbox environments without the full langchain dep tree.
    """
    src = _source_text("agent/directive_extractor.py")
    assert src is not None, "agent/directive_extractor.py not found"
    body = _function_body(src, "_call_llm_fallback")
    # Find the docstring (first triple-quoted string in the body).
    import re
    m = re.search(r'"""(.*?)"""', body, re.DOTALL)
    assert m is not None, "_call_llm_fallback has no docstring"
    doc = m.group(1).lower()
    has_property = (
      "no model-family variant cycling" in doc
      or ("single" in doc and "marker" in doc)
      or ("at most one" in doc and "marker" in doc))
    assert has_property, (
      "_call_llm_fallback docstring must document the single-call "
      "property to prevent future regression.  Current docstring:\n%s"
      % m.group(1))
    print("  PASS: test_fallback_property_documented")


# ---------- Baseline fingerprint test (v119.H1 ship snapshot) ------

def test_baseline_values_unchanged():
    """v119.H1 ship-time baseline -- catches accidental value drift.

    This is a canary against unintentional changes during the H1
    refactor itself.  Once H1 is stable, future legitimate default
    updates need to update both core/llm.py AND this baseline (and
    add the old value to RETIRED_MODELS).  That paired update is the
    intended developer workflow.
    """
    llm = _import_core_llm()
    # v120: anthropic + portkey entries added.  portkey fronts Azure OpenAI
    # (default gpt-5); anthropic chat defaults for decision/expensive roles.
    expected_decision = {
      "google":    "gemini-2.5-flash-lite",
      "openai":    "gpt-4o-mini",
      "anthropic": "claude-sonnet-4-6",
      "ollama":    "llama3.2",
      "portkey":   "gpt-5",
    }
    expected_rag = {
      "google":    "gemini-2.5-flash-lite",
      "openai":    "gpt-5-nano",
      "ollama":    "llama3.1:70b",
      "portkey":   "gpt-5",
    }
    expected_rag_embedding = {
      "google":    "gemini-embedding-001",
      "openai":    "text-embedding-3-small",
      "ollama":    "nomic-embed-text",
      # v120: portkey/Azure DEPLOYMENT name carries a trailing "-1"; the bare
      # name silently fell back to a chat model on the gateway.
      "portkey":   "text-embedding-3-small-1",
    }
    expected_expensive = {
      "google":    "gemini-2.5-pro",
      "openai":    "gpt-5",
      "ollama":    "qwen3:32b",
      "anthropic": "claude-opus-4-7",
      "portkey":   "gpt-5",
    }
    expected_cheap = {
      "ollama":    "qwen2.5:7b",
    }
    assert dict(llm.DECISION_MODEL_DEFAULTS) == expected_decision, (
      "DECISION_MODEL_DEFAULTS drifted from H1 baseline.\n"
      "  got:      %r\n  expected: %r" % (
        dict(llm.DECISION_MODEL_DEFAULTS), expected_decision))
    assert dict(llm.RAG_MODEL_DEFAULTS) == expected_rag, (
      "RAG_MODEL_DEFAULTS drifted from H1 baseline.\n"
      "  got:      %r\n  expected: %r" % (
        dict(llm.RAG_MODEL_DEFAULTS), expected_rag))
    assert dict(llm.RAG_EMBEDDING_DEFAULTS) == expected_rag_embedding, (
      "RAG_EMBEDDING_DEFAULTS drifted from H1 baseline.\n"
      "  got:      %r\n  expected: %r" % (
        dict(llm.RAG_EMBEDDING_DEFAULTS), expected_rag_embedding))
    assert dict(llm.EXPENSIVE_MODEL_DEFAULTS) == expected_expensive, (
      "EXPENSIVE_MODEL_DEFAULTS drifted from H1 baseline.\n"
      "  got:      %r\n  expected: %r" % (
        dict(llm.EXPENSIVE_MODEL_DEFAULTS), expected_expensive))
    assert dict(llm.CHEAP_MODEL_DEFAULTS) == expected_cheap, (
      "CHEAP_MODEL_DEFAULTS drifted from H1 baseline.\n"
      "  got:      %r\n  expected: %r" % (
        dict(llm.CHEAP_MODEL_DEFAULTS), expected_cheap))
    print("  PASS: test_baseline_values_unchanged")


# ---------- v119.H13 — resolve_model_for_provider env-var precedence -----

def test_resolve_model_for_provider_unset_env_falls_back_to_default():
    """When the provider's env-var override is unset, the resolver
    returns the central default (parity with default_model_for_provider).
    """
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    import os as _os
    saved = _os.environ.pop("OLLAMA_LLM_MODEL", None)
    try:
        got = llm.resolve_model_for_provider("ollama")
        expected = llm.default_model_for_provider("ollama")
        assert got == expected, (
            "resolve_model_for_provider with unset env-var should "
            "match default_model_for_provider.  Got %r, expected %r."
            % (got, expected))
    finally:
        if saved is not None:
            _os.environ["OLLAMA_LLM_MODEL"] = saved
    print("  PASS: test_resolve_model_for_provider_unset_env_falls_back_to_default")


def test_resolve_model_for_provider_set_env_overrides_default():
    """When OLLAMA_LLM_MODEL is set, resolver returns the env-var value,
    NOT the central default.  This is the contract H13 makes uniform
    across all callers (it was inconsistent pre-H13).
    """
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    import os as _os
    saved = _os.environ.pop("OLLAMA_LLM_MODEL", None)
    try:
        _os.environ["OLLAMA_LLM_MODEL"] = "qwen2.5:72b"
        got = llm.resolve_model_for_provider("ollama")
        assert got == "qwen2.5:72b", (
            "resolve_model_for_provider should honor OLLAMA_LLM_MODEL "
            "env-var.  Got %r." % got)
    finally:
        _os.environ.pop("OLLAMA_LLM_MODEL", None)
        if saved is not None:
            _os.environ["OLLAMA_LLM_MODEL"] = saved
    print("  PASS: test_resolve_model_for_provider_set_env_overrides_default")


def test_resolve_model_for_provider_empty_env_falls_through():
    """An EMPTY env-var (set to '') should NOT override the central
    default — empty string is treated as 'unset' to avoid the operator
    silently breaking their setup with `setenv OLLAMA_LLM_MODEL ""`.
    """
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    import os as _os
    saved = _os.environ.pop("OLLAMA_LLM_MODEL", None)
    try:
        _os.environ["OLLAMA_LLM_MODEL"] = ""
        got = llm.resolve_model_for_provider("ollama")
        expected = llm.default_model_for_provider("ollama")
        assert got == expected, (
            "Empty OLLAMA_LLM_MODEL should fall back to default. "
            "Got %r, expected %r." % (got, expected))
    finally:
        _os.environ.pop("OLLAMA_LLM_MODEL", None)
        if saved is not None:
            _os.environ["OLLAMA_LLM_MODEL"] = saved
    print("  PASS: test_resolve_model_for_provider_empty_env_falls_through")


def test_resolve_model_for_provider_unknown_provider_raises():
    """resolve_model_for_provider preserves default_model_for_provider's
    error contract: ValueError on unknown provider.
    """
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    raised = False
    try:
        llm.resolve_model_for_provider("xyz-not-a-provider")
    except ValueError:
        raised = True
    assert raised, "Unknown provider should raise ValueError"
    print("  PASS: test_resolve_model_for_provider_unknown_provider_raises")


# ---------- v119.H13 — normalize_ollama_openai_base_url ------------------

def test_normalize_ollama_openai_base_url_appends_v1():
    """Bare base URL gets /v1 appended (Tom's case)."""
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    cases = [
        ("http://localhost:11434", "http://localhost:11434/v1"),
        ("http://localhost:11434/", "http://localhost:11434/v1"),
        ("http://192.168.1.5:11434", "http://192.168.1.5:11434/v1"),
    ]
    for inp, want in cases:
        got = llm.normalize_ollama_openai_base_url(inp)
        assert got == want, (
            "normalize_ollama_openai_base_url(%r) = %r, expected %r"
            % (inp, got, want))
    print("  PASS: test_normalize_ollama_openai_base_url_appends_v1")


def test_normalize_ollama_openai_base_url_idempotent():
    """URL already ending in /v1 is unchanged (modulo trailing slash)."""
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    cases = [
        ("http://localhost:11434/v1", "http://localhost:11434/v1"),
        ("http://localhost:11434/v1/", "http://localhost:11434/v1"),
    ]
    for inp, want in cases:
        got = llm.normalize_ollama_openai_base_url(inp)
        assert got == want, (
            "normalize_ollama_openai_base_url(%r) = %r, expected %r"
            % (inp, got, want))
    print("  PASS: test_normalize_ollama_openai_base_url_idempotent")


def test_normalize_ollama_openai_base_url_passes_through_none():
    """None passes through unchanged (caller may pass None)."""
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    assert llm.normalize_ollama_openai_base_url(None) is None
    print("  PASS: test_normalize_ollama_openai_base_url_passes_through_none")


def test_normalize_ollama_openai_base_url_passes_through_empty():
    """Empty string passes through unchanged (does NOT become '/v1').

    Caller's default-handling decides what to do; we don't want to
    silently produce a useless '/v1' URL fragment from an empty
    OLLAMA_BASE_URL.
    """
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        from core import llm

    assert llm.normalize_ollama_openai_base_url("") == ""
    print("  PASS: test_normalize_ollama_openai_base_url_passes_through_empty")


# ---------- Runner ------------------------------------------------

def run_all_tests():
    # Table-shape (6)
    test_decision_models_table_complete()
    test_rag_models_table_minimum()
    test_rag_embedding_table_minimum()
    test_expensive_models_table_minimum()
    test_cheap_models_table_minimum()
    test_no_default_in_retired_list()
    # Helper-function (5)
    test_default_model_for_provider_basic()
    test_default_model_for_provider_role()
    test_default_model_for_provider_normalizes_keys()
    test_default_model_for_provider_unknown_role_raises()
    test_default_model_for_provider_unknown_provider_raises()
    # Per-callsite source-scan (6)
    test_api_client_google_uses_central_default()
    test_api_client_openai_uses_central_default()
    test_api_client_anthropic_uses_central_default()
    test_api_client_ollama_uses_central_default()
    test_directive_extractor_fallback_uses_central_defaults()
    test_test_api_keys_uses_central_defaults()
    # AST-level structural regression (1)
    test_no_orphan_model_strings_in_agent_runtime_code()
    # Retired-model helpers (4)
    test_retired_error_detection_positive_cases()
    test_retired_error_detection_negative_cases()
    test_retired_marker_writes_expected_fields()
    test_fallback_property_documented()
    # Baseline fingerprint (1)
    test_baseline_values_unchanged()
    # v119.H13 — resolve_model_for_provider + URL normalize (7)
    test_resolve_model_for_provider_unset_env_falls_back_to_default()
    test_resolve_model_for_provider_set_env_overrides_default()
    test_resolve_model_for_provider_empty_env_falls_through()
    test_resolve_model_for_provider_unknown_provider_raises()
    test_normalize_ollama_openai_base_url_appends_v1()
    test_normalize_ollama_openai_base_url_idempotent()
    test_normalize_ollama_openai_base_url_passes_through_none()
    test_normalize_ollama_openai_base_url_passes_through_empty()


if __name__ == "__main__":
    run_all_tests()
