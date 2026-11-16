"""
OpenAI / Google reasoning-model handling (temperature, output budget, reply
readers) and the matching provider defaults.

Background (live probe 2026-10-14): gpt-5.5 and gpt-6-luna reject a
non-default temperature (HTTP 400); gemini-3.5-flash-lite ignores temperature
with a warning; gemini-2.5-pro and gpt-6-luna spent a whole 2000-token cap on
reasoning and returned no text (finish MAX_TOKENS / length).

Scope: pure logic, source scans, and SDK modules faked through sys.modules.
No SDK imports, no keys, no network.  Replies are faked with
types.SimpleNamespace / small classes.

Covers:
  - P1 openai_model_is_reasoning: True for gpt-5*, gpt-6*, gpt-10*, o1,
    o3-mini, o4-mini; False for gpt-4o*, gpt-5-chat-latest, gpt-3.5-turbo,
    "" and None.
  - P2 google_model_ignores_sampling: True for Gemini 3+ (incl. "models/"
    prefix); False for Gemini 2.5 / 1.5, gemini-embedding-001, "" and None.
  - P3 openai_completion_kwargs: reasoning model -> only
    max_completion_tokens=REASONING_MAX_TOKENS; gpt-4o-mini -> temperature
    and max_tokens as given.
  - P4 google_generation_kwargs: max_output_tokens 16000 for both Gemini 3.5
    and 2.5; temperature only for 2.5.
  - P5 openai_reply_text: raises on finish_reason "length" (even with text)
    and on empty content; returns the text for "stop".
  - P6 google_reply_text: raises on finish_reason MAX_TOKENS (enum-like, even
    with text), on empty text, and when response.text raises ValueError;
    returns the text for STOP.
  - P7 get_llm_and_embeddings (langchain_google_genai / langchain_openai
    faked): temperature omitted for gemini-3.5-flash-lite and gpt-6-luna,
    present for gemini-2.5-flash-lite and gpt-4o-mini.
  - P8 api_client._call_openai_llm / _call_google_llm (openai / google.genai
    faked) send the reasoning-safe kwargs and return the reply text.
  - P9 DECISION / RAG / EXPENSIVE / embedding defaults for google and openai,
    and no default is in RETIRED_MODELS.
  - P10 Source scan: directive_extractor.py's google and openai fallback
    calls use google_generation_kwargs / openai_completion_kwargs and
    google_reply_text / openai_reply_text, with no temperature=0.1.
  - P11 openai_model_is_reasoning on gateway / fine-tune name forms
    ("openai/", "azure/", "ft:" prefixes) and non-reasoning look-alikes
    (chatgpt-4o-latest, gpt-oss:20b, gemma3, gpt-5-chat-latest).
  - P12 Cut-off / empty-reply RuntimeErrors from anthropic_reply_text,
    openai_reply_text and google_reply_text carry no 3+ digit number (a
    character count such as 1429 contains "429") and are not classified as
    rate limits by rate_limit_handler.is_rate_limit_error.
  - P13 google_reply_text raises RuntimeError naming SAFETY for a blocked
    reply, and raises (no crash) for an int finish_reason with empty text.
  - C1 response_text raises RuntimeError for a langchain reply whose
    response_metadata reports a cut-off (stop_reason max_tokens /
    model_context_window_exceeded, finish_reason length / MAX_TOKENS incl.
    enum-like), returns the text otherwise; errors carry no 3+ digit number
    and are not classified as rate limits.
  - C2 get_llm_and_embeddings raises a short timeout to
    HOSTED_LLM_MIN_TIMEOUT (180) for anthropic / openai / google chat
    clients (langchain faked); a longer timeout is kept.
  - C3 portkey_max_completion_tokens: reasoning models floored at 16000,
    others unchanged.
  - C4 api_client._call_portkey_llm (portkey_ai faked) sends
    max_completion_tokens 16000 for gpt-5, returns the text, and raises
    RuntimeError for finish_reason "length".
  - C5 Source scan: directive_extractor.py's portkey fallback uses
    portkey_max_completion_tokens( and returns openai_reply_text(response).
  - C6 analysis/summarizer.summarize_log_text (langchain_core /
    langchain_text_splitters faked, fake LLM): a map reply cut off with
    stop_reason max_tokens is skipped with a "Chunk N of batch M skipped"
    message; the run completes and the reduce step combines only the
    surviving chunk summaries.

Fail-fast: the first failing test stops the run with a non-zero exit;
AssertionError is never swallowed.
"""
import asyncio
import contextlib
import importlib.util
import io
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

_MISSING = object()


class _FakeModules(object):
    """Temporarily install fake modules in sys.modules and set/clear env
    vars (value None = clear); restore both on exit."""
    def __init__(self, modules, env=None):
        self.modules = modules
        self.env = env or {}

    def __enter__(self):
        self._saved_mods = {}
        for name, mod in self.modules.items():
            self._saved_mods[name] = sys.modules.get(name, _MISSING)
            sys.modules[name] = mod
        self._saved_env = {}
        for name, value in self.env.items():
            self._saved_env[name] = os.environ.get(name, _MISSING)
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
        return self

    def __exit__(self, *exc):
        for name, old in self._saved_mods.items():
            if old is _MISSING:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = old
        for name, old in self._saved_env.items():
            if old is _MISSING:
                os.environ.pop(name, None)
            else:
                os.environ[name] = old
        return False


def _ensure_libtbx_importable():
    """api_client imports libtbx.langchain.*; add the checkout root when
    libtbx is not already importable (plain python3 from a scratch dir)."""
    try:
        import libtbx  # noqa: F401
    except ImportError:
        sys.path.append(os.path.dirname(os.path.dirname(_ROOT)))


def _expect_runtime_error(fn, arg, what):
    try:
        out = fn(arg)
    except RuntimeError as e:
        return e
    raise AssertionError("expected RuntimeError for %s, got %r" % (what, out))


# ---------------------------------------------------------------------------
# P1. openai_model_is_reasoning
# ---------------------------------------------------------------------------

def test_openai_model_is_reasoning():
    from core.llm import openai_model_is_reasoning as f
    for name in ["gpt-5", "gpt-5.5", "gpt-5.4-nano", "gpt-6-luna",
                 "gpt-6.1-sol", "gpt-10-x", "o1", "o3-mini", "o4-mini"]:
        assert f(name) is True, "expected True for %r, got %r" % (name, f(name))
    for name in ["gpt-4o-mini", "gpt-4o", "gpt-5-chat-latest",
                 "gpt-3.5-turbo", "", None]:
        assert f(name) is False, \
            "expected False for %r, got %r" % (name, f(name))


# ---------------------------------------------------------------------------
# P2. google_model_ignores_sampling
# ---------------------------------------------------------------------------

def test_google_model_ignores_sampling():
    from core.llm import google_model_ignores_sampling as f
    for name in ["gemini-3.5-flash-lite", "gemini-3.8-flash",
                 "gemini-4-argon", "models/gemini-3.5-flash-lite"]:
        assert f(name) is True, "expected True for %r, got %r" % (name, f(name))
    for name in ["gemini-2.5-flash-lite", "gemini-2.5-pro", "gemini-1.5-pro",
                 "gemini-embedding-001", "", None]:
        assert f(name) is False, \
            "expected False for %r, got %r" % (name, f(name))


# ---------------------------------------------------------------------------
# P3. openai_completion_kwargs
# ---------------------------------------------------------------------------

def test_openai_completion_kwargs():
    from core.llm import openai_completion_kwargs
    from core.llm import REASONING_MAX_TOKENS
    assert REASONING_MAX_TOKENS == 16000, REASONING_MAX_TOKENS
    got = openai_completion_kwargs("gpt-6-luna", 0.1, 2000)
    assert got == {"max_completion_tokens": 16000}, got
    got = openai_completion_kwargs("gpt-4o-mini", 0.1, 2000)
    assert got == {"temperature": 0.1, "max_tokens": 2000}, got


# ---------------------------------------------------------------------------
# P4. google_generation_kwargs
# ---------------------------------------------------------------------------

def test_google_generation_kwargs():
    from core.llm import google_generation_kwargs
    got = google_generation_kwargs("gemini-3.5-flash-lite", 0.1, 2000)
    assert got.get("max_output_tokens") == 16000, got
    assert "temperature" not in got, got
    got = google_generation_kwargs("gemini-2.5-pro", 0.1, 2000)
    assert got.get("max_output_tokens") == 16000, got
    assert got.get("temperature") == 0.1, got


# ---------------------------------------------------------------------------
# P5. openai_reply_text
# ---------------------------------------------------------------------------

def _openai_reply(content, finish_reason):
    return types.SimpleNamespace(choices=[types.SimpleNamespace(
        finish_reason=finish_reason,
        message=types.SimpleNamespace(content=content))])


def test_openai_reply_text():
    from core.llm import openai_reply_text
    err = _expect_runtime_error(
        openai_reply_text, _openai_reply("partial answ", "length"),
        "finish_reason 'length' with text")
    assert "length" in str(err), str(err)
    _expect_runtime_error(
        openai_reply_text, _openai_reply("", "stop"), "empty content")
    _expect_runtime_error(
        openai_reply_text, _openai_reply(None, "stop"), "None content")
    out = openai_reply_text(_openai_reply("answer", "stop"))
    assert out == "answer", repr(out)


# ---------------------------------------------------------------------------
# P6. google_reply_text
# ---------------------------------------------------------------------------

class _FinishReason(object):
    """Enum-like finish_reason (google.genai FinishReason stand-in)."""
    def __init__(self, name):
        self.name = name
        self.value = name

    def __repr__(self):
        return "<FinishReason.%s>" % self.name

    __str__ = __repr__


def _google_reply(text, finish_name):
    return types.SimpleNamespace(
        text=text,
        candidates=[types.SimpleNamespace(
            finish_reason=_FinishReason(finish_name))])


class _TextRaises(object):
    """A reply whose .text accessor raises ValueError (no text parts)."""
    def __init__(self, finish_name):
        self.candidates = [types.SimpleNamespace(
            finish_reason=_FinishReason(finish_name),
            content=types.SimpleNamespace(parts=[]))]

    @property
    def text(self):
        raise ValueError("response has no text parts")


def test_google_reply_text():
    from core.llm import google_reply_text
    err = _expect_runtime_error(
        google_reply_text, _google_reply("partial answ", "MAX_TOKENS"),
        "finish_reason MAX_TOKENS with text")
    assert "MAX_TOKENS" in str(err), str(err)
    _expect_runtime_error(
        google_reply_text, _google_reply("", "STOP"), "empty text")
    _expect_runtime_error(
        google_reply_text, _TextRaises("STOP"), "response.text ValueError")
    out = google_reply_text(_google_reply("answer", "STOP"))
    assert out == "answer", repr(out)


# ---------------------------------------------------------------------------
# P7. get_llm_and_embeddings temperature handling (langchain faked)
# ---------------------------------------------------------------------------

def _recording_class(calls, tag):
    class _Rec(object):
        def __init__(self, **kw):
            calls.append((tag, kw))
    return _Rec


def _chat_kwargs(provider, model_name):
    from core.llm import get_llm_and_embeddings
    calls = []
    g = types.ModuleType("langchain_google_genai")
    g.ChatGoogleGenerativeAI = _recording_class(calls, "chat")
    g.GoogleGenerativeAIEmbeddings = _recording_class(calls, "emb")
    o = types.ModuleType("langchain_openai")
    o.ChatOpenAI = _recording_class(calls, "chat")
    o.OpenAIEmbeddings = _recording_class(calls, "emb")
    fakes = {"langchain_google_genai": g, "langchain_openai": o}
    env = {"GOOGLE_API_KEY": "dummy-not-a-key",
           "OPENAI_API_KEY": "dummy-not-a-key"}
    with _FakeModules(fakes, env=env):
        get_llm_and_embeddings(
            provider=provider, llm_model_name=model_name, temperature=0.3)
    chats = [kw for tag, kw in calls if tag == "chat"]
    assert len(chats) == 1, calls
    assert chats[0].get("model") == model_name, chats[0]
    return chats[0]


def test_get_llm_temperature_by_model():
    kw = _chat_kwargs("google", "gemini-3.5-flash-lite")
    assert "temperature" not in kw, kw
    kw = _chat_kwargs("google", "gemini-2.5-flash-lite")
    assert kw.get("temperature") == 0.3, kw
    kw = _chat_kwargs("openai", "gpt-6-luna")
    assert "temperature" not in kw, kw
    kw = _chat_kwargs("openai", "gpt-4o-mini")
    assert kw.get("temperature") == 0.3, kw


# ---------------------------------------------------------------------------
# P8. api_client direct calls (openai / google.genai faked)
# ---------------------------------------------------------------------------

def _openai_call_kwargs(model_name):
    _ensure_libtbx_importable()
    from agent.api_client import _call_openai_llm
    calls = []

    class _Completions(object):
        def create(self, **kw):
            calls.append(kw)
            return _openai_reply("pong", "stop")

    class OpenAI(object):
        def __init__(self, *a, **kw):
            self.chat = types.SimpleNamespace(completions=_Completions())

    fake = types.ModuleType("openai")
    fake.OpenAI = OpenAI
    with _FakeModules({"openai": fake}):
        out = _call_openai_llm("hi", model=model_name)
    assert out == "pong", repr(out)
    assert len(calls) == 1, calls
    assert calls[0].get("model") == model_name, calls[0]
    return calls[0]


def test_call_openai_llm_kwargs():
    kw = _openai_call_kwargs("gpt-6-luna")
    assert kw.get("max_completion_tokens") == 16000, kw
    assert "temperature" not in kw, kw
    assert "max_tokens" not in kw, kw
    kw = _openai_call_kwargs("gpt-4o-mini")
    assert "temperature" in kw, kw
    assert kw.get("max_tokens") == 2000, kw
    assert "max_completion_tokens" not in kw, kw


def _google_call_kwargs(model_name):
    _ensure_libtbx_importable()
    from agent.api_client import _call_google_llm
    gen_calls = []
    config_calls = []

    class GenerateContentConfig(object):
        def __init__(self, **kw):
            config_calls.append(kw)
            self.kw = kw

    class _Models(object):
        def generate_content(self, **kw):
            gen_calls.append(kw)
            return _google_reply("pong", "STOP")

    class Client(object):
        def __init__(self, *a, **kw):
            self.models = _Models()

    gtypes = types.ModuleType("google.genai.types")
    gtypes.GenerateContentConfig = GenerateContentConfig
    genai = types.ModuleType("google.genai")
    genai.Client = Client
    genai.types = gtypes
    google = types.ModuleType("google")
    google.__path__ = []
    google.genai = genai
    fakes = {"google": google, "google.genai": genai,
             "google.genai.types": gtypes}
    with _FakeModules(fakes, env={"GOOGLE_API_KEY": "dummy-not-a-key"}):
        out = _call_google_llm("hi", model=model_name)
    assert out == "pong", repr(out)
    assert len(gen_calls) == 1, gen_calls
    assert gen_calls[0].get("model") == model_name, gen_calls[0]
    assert len(config_calls) == 1, config_calls
    return config_calls[0]


def test_call_google_llm_kwargs():
    kw = _google_call_kwargs("gemini-3.5-flash-lite")
    assert kw.get("max_output_tokens") == 16000, kw
    assert "temperature" not in kw, kw
    kw = _google_call_kwargs("gemini-2.5-flash-lite")
    assert kw.get("max_output_tokens") == 16000, kw
    assert "temperature" in kw, kw


# ---------------------------------------------------------------------------
# P9. Defaults
# ---------------------------------------------------------------------------

def test_google_openai_defaults():
    from core.llm import DECISION_MODEL_DEFAULTS
    from core.llm import RAG_MODEL_DEFAULTS
    from core.llm import EXPENSIVE_MODEL_DEFAULTS
    from core.llm import RAG_EMBEDDING_DEFAULTS
    from core.llm import RETIRED_MODELS
    want = [
        (DECISION_MODEL_DEFAULTS, "google", "gemini-3.5-flash-lite"),
        (RAG_MODEL_DEFAULTS, "google", "gemini-3.5-flash-lite"),
        (EXPENSIVE_MODEL_DEFAULTS, "google", "gemini-3.8-flash"),
        (DECISION_MODEL_DEFAULTS, "openai", "gpt-6-luna"),
        (RAG_MODEL_DEFAULTS, "openai", "gpt-6-luna"),
        (EXPENSIVE_MODEL_DEFAULTS, "openai", "gpt-5.5"),
        (RAG_EMBEDDING_DEFAULTS, "google", "gemini-embedding-001"),
        (RAG_EMBEDDING_DEFAULTS, "openai", "text-embedding-3-small"),
    ]
    for table, provider, model in want:
        got = table.get(provider)
        assert got == model, "%s: expected %r, got %r" % (provider, model, got)
    for table in [DECISION_MODEL_DEFAULTS, RAG_MODEL_DEFAULTS,
                  EXPENSIVE_MODEL_DEFAULTS, RAG_EMBEDDING_DEFAULTS]:
        for provider, model in table.items():
            assert model not in RETIRED_MODELS, \
                "%s default %r is in RETIRED_MODELS" % (provider, model)


# ---------------------------------------------------------------------------
# P10. Source scan: directive_extractor google / openai fallback calls
# ---------------------------------------------------------------------------

def _fallback_provider_segments():
    """Return (google_code, openai_code) from _call_llm_fallback, comments
    stripped: google = 'if provider == "google":' up to the openai branch,
    openai = 'elif provider == "openai":' up to the next elif."""
    path = os.path.join(_AGENT, "directive_extractor.py")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("def _call_llm_fallback("):
            start = i
            break
    assert start is not None, "no _call_llm_fallback in directive_extractor"
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith("def "):
            end = i
            break
    body = [ln.split("#", 1)[0] for ln in lines[start:end]]
    g = o = a = None
    for i, line in enumerate(body):
        s = line.strip()
        if g is None and s == 'if provider == "google":':
            g = i
        elif g is not None and o is None and s == 'elif provider == "openai":':
            o = i
        elif o is not None and a is None and s.startswith("elif provider =="):
            a = i
    assert g is not None and o is not None and a is not None, (g, o, a)
    return "\n".join(body[g:o]), "\n".join(body[o:a])


def test_directive_extractor_google_openai_fallback():
    google_code, openai_code = _fallback_provider_segments()
    assert "google_generation_kwargs(" in google_code, \
        "google fallback does not use google_generation_kwargs"
    assert "return google_reply_text(response)" in google_code, \
        "google fallback does not return google_reply_text(response)"
    assert "openai_completion_kwargs(" in openai_code, \
        "openai fallback does not use openai_completion_kwargs"
    assert "return openai_reply_text(response)" in openai_code, \
        "openai fallback does not return openai_reply_text(response)"
    for label, code in [("google", google_code), ("openai", openai_code)]:
        assert not re.search(r"\btemperature\s*=\s*0\.1\b", code), \
            "%s fallback still passes temperature=0.1" % label
        assert "response.text" not in code.replace(
            "google_reply_text(response)", ""), \
            "%s fallback still returns response.text" % label
        assert "choices[0].message.content" not in code, \
            "%s fallback still reads choices[0].message.content" % label



# ---------------------------------------------------------------------------
# P11. openai_model_is_reasoning: gateway / fine-tune name forms
# ---------------------------------------------------------------------------

def test_openai_model_is_reasoning_prefixed_names():
    from core.llm import openai_model_is_reasoning as f
    for name in ["openai/gpt-6-luna", "azure/gpt-5", "ft:gpt-5-mini:org::abc",
                 "openai/o4-mini"]:
        assert f(name) is True, "expected True for %r, got %r" % (name, f(name))
    for name in ["openai/gpt-4o-mini", "chatgpt-4o-latest", "gpt-oss:20b",
                 "gemma3", "gpt-5-chat-latest"]:
        assert f(name) is False, \
            "expected False for %r, got %r" % (name, f(name))


# ---------------------------------------------------------------------------
# P12. Cut-off errors are not mistaken for rate limits
# ---------------------------------------------------------------------------

class _AnthropicBlock(object):
    def __init__(self, type, **kw):
        self.type = type
        for k, v in kw.items():
            setattr(self, k, v)


def _anthropic_reply(text, stop_reason):
    return types.SimpleNamespace(stop_reason=stop_reason, content=[
        _AnthropicBlock("text", text=text)])


def test_cutoff_errors_not_rate_limits():
    _ensure_libtbx_importable()
    from core.llm import anthropic_reply_text
    from core.llm import openai_reply_text
    from core.llm import google_reply_text
    from agent.rate_limit_handler import is_rate_limit_error
    long_text = "x" * 1429
    cases = [
        ("anthropic max_tokens", anthropic_reply_text,
         _anthropic_reply(long_text, "max_tokens")),
        ("anthropic model_context_window_exceeded", anthropic_reply_text,
         _anthropic_reply(long_text, "model_context_window_exceeded")),
        ("anthropic empty", anthropic_reply_text,
         _anthropic_reply("", "end_turn")),
        ("openai length", openai_reply_text,
         _openai_reply(long_text, "length")),
        ("openai empty", openai_reply_text, _openai_reply("", "stop")),
        ("google MAX_TOKENS", google_reply_text,
         _google_reply(long_text, "MAX_TOKENS")),
        ("google empty", google_reply_text, _google_reply("", "STOP")),
    ]
    for label, reader, reply in cases:
        err = _expect_runtime_error(reader, reply, label)
        msg = str(err)
        assert re.search(r"\d{3}", msg) is None, \
            "%s: message has a 3+ digit number: %r" % (label, msg)
        assert is_rate_limit_error(err) is False, \
            "%s: classified as a rate limit: %r" % (label, msg)


# ---------------------------------------------------------------------------
# P13. google_reply_text: SAFETY and int finish_reason
# ---------------------------------------------------------------------------

def test_google_reply_text_safety_and_int_finish():
    from core.llm import google_reply_text
    err = _expect_runtime_error(
        google_reply_text, _google_reply("", "SAFETY"), "SAFETY block")
    assert "SAFETY" in str(err), str(err)
    reply = types.SimpleNamespace(
        text="", candidates=[types.SimpleNamespace(finish_reason=2)])
    _expect_runtime_error(google_reply_text, reply, "int finish_reason")



# ---------------------------------------------------------------------------
# C1. response_text: langchain reply cut-off detection
# ---------------------------------------------------------------------------

def test_response_text_langchain_cutoff():
    _ensure_libtbx_importable()
    from core.llm import response_text
    from agent.rate_limit_handler import is_rate_limit_error
    cut = [
        {"stop_reason": "max_tokens"},
        {"stop_reason": "model_context_window_exceeded"},
        {"finish_reason": "length"},
        {"finish_reason": "MAX_TOKENS"},
        {"finish_reason": _FinishReason("MAX_TOKENS")},
    ]
    for meta in cut:
        reply = types.SimpleNamespace(content="partial", response_metadata=meta)
        err = _expect_runtime_error(response_text, reply, "metadata %r" % meta)
        msg = str(err)
        assert re.search(r"\d{3}", msg) is None, \
            "%r: message has a 3+ digit number: %r" % (meta, msg)
        assert is_rate_limit_error(err) is False, \
            "%r: classified as a rate limit: %r" % (meta, msg)
    ok = [
        {"stop_reason": "end_turn"},
        {"finish_reason": "stop"},
        {"finish_reason": "STOP"},
        {},
    ]
    for meta in ok:
        reply = types.SimpleNamespace(content="partial", response_metadata=meta)
        out = response_text(reply)
        assert out == "partial", "%r: got %r" % (meta, out)
    for reply in [types.SimpleNamespace(content="partial"),
                  types.SimpleNamespace(content="partial",
                                        response_metadata="max_tokens"),
                  types.SimpleNamespace(content="partial",
                                        response_metadata=None)]:
        out = response_text(reply)
        assert out == "partial", repr(out)


# ---------------------------------------------------------------------------
# C2. Hosted chat clients get a minimum timeout
# ---------------------------------------------------------------------------

def _chat_timeout(provider, timeout):
    from core.llm import get_llm_and_embeddings
    from core.llm import DECISION_MODEL_DEFAULTS
    calls = []
    a = types.ModuleType("langchain_anthropic")
    a.ChatAnthropic = _recording_class(calls, "chat")
    g = types.ModuleType("langchain_google_genai")
    g.ChatGoogleGenerativeAI = _recording_class(calls, "chat")
    g.GoogleGenerativeAIEmbeddings = _recording_class(calls, "emb")
    o = types.ModuleType("langchain_openai")
    o.ChatOpenAI = _recording_class(calls, "chat")
    o.OpenAIEmbeddings = _recording_class(calls, "emb")
    fakes = {"langchain_anthropic": a, "langchain_google_genai": g,
             "langchain_openai": o}
    env = {"GOOGLE_API_KEY": "dummy-not-a-key",
           "OPENAI_API_KEY": "dummy-not-a-key"}
    model = DECISION_MODEL_DEFAULTS[provider]
    with _FakeModules(fakes, env=env):
        get_llm_and_embeddings(
            provider=provider, llm_model_name=model, timeout=timeout)
    chats = [kw for tag, kw in calls if tag == "chat"]
    assert len(chats) == 1, calls
    return chats[0].get("timeout")


def test_hosted_llm_min_timeout():
    from core.llm import HOSTED_LLM_MIN_TIMEOUT
    assert HOSTED_LLM_MIN_TIMEOUT == 180, HOSTED_LLM_MIN_TIMEOUT
    for provider in ["anthropic", "openai", "google"]:
        got = _chat_timeout(provider, 60)
        assert got == 180, "%s timeout=60 -> %r" % (provider, got)
        got = _chat_timeout(provider, 300)
        assert got == 300, "%s timeout=300 -> %r" % (provider, got)


# ---------------------------------------------------------------------------
# C3. portkey_max_completion_tokens
# ---------------------------------------------------------------------------

def test_portkey_max_completion_tokens():
    from core.llm import portkey_max_completion_tokens as f
    for model, requested, want in [("gpt-5", 2000, 16000),
                                   ("gpt-6-luna", 20000, 20000),
                                   ("gpt-4o", 2000, 2000)]:
        got = f(model, requested)
        assert got == want, "(%r, %r): expected %r, got %r" % (
            model, requested, want, got)


# ---------------------------------------------------------------------------
# C4. api_client._call_portkey_llm (portkey_ai faked)
# ---------------------------------------------------------------------------

def _portkey_call(model_name, finish_reason):
    _ensure_libtbx_importable()
    from agent.api_client import _call_portkey_llm
    calls = []

    class _Completions(object):
        def create(self, **kw):
            calls.append(kw)
            return _openai_reply("pong", finish_reason)

    class Portkey(object):
        def __init__(self, *a, **kw):
            self.chat = types.SimpleNamespace(completions=_Completions())

    fake = types.ModuleType("portkey_ai")
    fake.Portkey = Portkey
    env = {"PORTKEY_AZURE_API_KEY": "dummy-not-a-key",
           "PORTKEY_BASE_URL": "http://127.0.0.1:9/not-used"}
    out = err = None
    with _FakeModules({"portkey_ai": fake}, env=env):
        try:
            out = _call_portkey_llm("hi", model=model_name)
        except RuntimeError as e:
            err = e
    assert len(calls) == 1, calls
    return calls[0], out, err


def test_call_portkey_llm():
    kw, out, err = _portkey_call("gpt-5", "stop")
    assert err is None, str(err)
    assert out == "pong", repr(out)
    assert kw.get("model") == "gpt-5", kw
    assert kw.get("max_completion_tokens") == 16000, kw
    kw, out, err = _portkey_call("gpt-5", "length")
    assert err is not None, \
        "expected RuntimeError for finish_reason 'length', got %r" % out


# ---------------------------------------------------------------------------
# C5. Source scan: directive_extractor portkey fallback
# ---------------------------------------------------------------------------

def _fallback_branch(provider):
    """Code (comments stripped) of the 'elif provider == "<provider>":'
    branch of _call_llm_fallback, up to the next provider branch."""
    path = os.path.join(_AGENT, "directive_extractor.py")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("def _call_llm_fallback("):
            start = i
            break
    assert start is not None, "no _call_llm_fallback in directive_extractor"
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith("def "):
            end = i
            break
    body = [ln.split("#", 1)[0] for ln in lines[start:end]]
    head = 'elif provider == "%s":' % provider
    b = None
    for i, line in enumerate(body):
        if line.strip() == head:
            b = i
            break
    assert b is not None, "no %r branch in _call_llm_fallback" % head
    indent = len(body[b]) - len(body[b].lstrip())
    e = len(body)
    for i in range(b + 1, len(body)):
        s = body[i]
        if s.strip() and (len(s) - len(s.lstrip())) <= indent \
                and s.strip().startswith(("elif ", "else")):
            e = i
            break
    return "\n".join(body[b:e])


def test_directive_extractor_portkey_fallback():
    code = _fallback_branch("portkey")
    assert "portkey_max_completion_tokens(" in code, \
        "portkey fallback does not use portkey_max_completion_tokens("
    assert "return openai_reply_text(response)" in code, \
        "portkey fallback does not return openai_reply_text(response)"
    assert "choices[0].message.content" not in code, \
        "portkey fallback still reads choices[0].message.content"



# ---------------------------------------------------------------------------
# C6. summarize_log_text skips a cut-off map chunk
# ---------------------------------------------------------------------------

_CHUNK_SPLIT = "\n=====FAKE-CHUNK-SPLIT=====\n"


def _fake_langchain_for_summarizer():
    """Minimal langchain_core / langchain_text_splitters stand-ins: enough for
    analysis/summarizer.py to import and run summarize_log_text."""
    class Document(object):
        def __init__(self, page_content, metadata=None):
            self.page_content = page_content
            self.metadata = metadata or {}

    class _Chain(object):
        def __init__(self, prompt, llm):
            self.prompt = prompt
            self.llm = llm

        async def ainvoke(self, inputs):
            return await self.llm.ainvoke(self.prompt.format(**inputs))

    class PromptTemplate(object):
        def __init__(self, template, input_variables=None, partials=None):
            self.template = template
            self.input_variables = input_variables or []
            self.partials = dict(partials or {})

        @classmethod
        def from_template(cls, template):
            return cls(template)

        def partial(self, **kw):
            merged = dict(self.partials)
            merged.update(kw)
            return PromptTemplate(self.template, self.input_variables, merged)

        def format(self, **kw):
            values = dict(self.partials)
            values.update(kw)
            out = self.template
            for k, v in values.items():
                out = out.replace("{%s}" % k, str(v))
            return out

        def __or__(self, llm):
            return _Chain(self, llm)

    class HumanMessage(object):
        def __init__(self, content):
            self.content = content

    class RecursiveCharacterTextSplitter(object):
        def __init__(self, chunk_size=None, chunk_overlap=None):
            pass

        def split_documents(self, docs):
            out = []
            for d in docs:
                for part in d.page_content.split(_CHUNK_SPLIT):
                    out.append(Document(page_content=part))
            return out

    core = types.ModuleType("langchain_core")
    core.__path__ = []
    documents = types.ModuleType("langchain_core.documents")
    documents.Document = Document
    prompts = types.ModuleType("langchain_core.prompts")
    prompts.PromptTemplate = PromptTemplate
    messages = types.ModuleType("langchain_core.messages")
    messages.HumanMessage = HumanMessage
    core.documents = documents
    core.prompts = prompts
    core.messages = messages
    splitters = types.ModuleType("langchain_text_splitters")
    splitters.RecursiveCharacterTextSplitter = RecursiveCharacterTextSplitter
    return {"langchain_core": core,
            "langchain_core.documents": documents,
            "langchain_core.prompts": prompts,
            "langchain_core.messages": messages,
            "langchain_text_splitters": splitters}


class _FakeSummaryLLM(object):
    """Map step: chunk containing CUT-CHUNK gets a reply cut off by
    max_tokens; other chunks get 'SUMMARY-OF-<tag>'.  Reduce step records
    the combine prompt and returns FINAL-SUMMARY."""
    model = "fake-model"

    def __init__(self):
        self.reduce_prompts = []

    async def ainvoke(self, prompt_text):
        for tag in ["CHUNK-A", "CHUNK-C"]:
            if tag in prompt_text:
                return types.SimpleNamespace(
                    content="SUMMARY-OF-%s" % tag,
                    response_metadata={"stop_reason": "end_turn"})
        assert "CUT-CHUNK" in prompt_text, prompt_text[-300:]
        return types.SimpleNamespace(
            content="TRUNCATED-PARTIAL-SUMMARY",
            response_metadata={"stop_reason": "max_tokens"})

    def invoke(self, messages):
        self.reduce_prompts.append(messages[0].content)
        return types.SimpleNamespace(
            content="FINAL-SUMMARY",
            response_metadata={"stop_reason": "end_turn"})


def test_summarizer_skips_cutoff_chunk():
    _ensure_libtbx_importable()
    path = os.path.join(_ROOT, "analysis", "summarizer.py")
    assert os.path.isfile(path), path
    filler = "phenix.refine R-work 0.21 R-free 0.25 Space group: P 21 21 21\n"
    chunks = []
    for tag in ["CHUNK-A", "CUT-CHUNK", "CHUNK-C"]:
        chunks.append(tag + "\n" + filler * 10)
    text = _CHUNK_SPLIT.join(chunks)
    llm = _FakeSummaryLLM()
    out = io.StringIO()
    with _FakeModules(_fake_langchain_for_summarizer()):
        # Load under a private name so sys.modules is not left holding a
        # summarizer bound to the fake langchain modules.
        spec = importlib.util.spec_from_file_location(
            "_tst_fake_summarizer", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with contextlib.redirect_stdout(out):
            result = asyncio.run(mod.summarize_log_text(
                text, llm, timeout=30, batch_size=3,
                use_throttling=False, provider="openai",
                program_name="phenix.refine"))
    printed = out.getvalue()
    assert result.error is None, result.error
    # Since the coverage note (2026-10-15) the summary says a part was skipped.
    assert result.log_summary.startswith("FINAL-SUMMARY\n\n**Note: partial coverage.**"), \
      repr(result.log_summary)
    assert "2 of 3 parts" in result.log_summary, repr(result.log_summary)
    assert (result.parts_total, result.parts_skipped) == (3, 1), (
      result.parts_total, result.parts_skipped)
    assert "Chunk 2 of batch 1 skipped" in printed, printed[-2000:]
    assert len(llm.reduce_prompts) == 1, llm.reduce_prompts
    combined = llm.reduce_prompts[0]
    assert "SUMMARY-OF-CHUNK-A" in combined, combined[-1000:]
    assert "SUMMARY-OF-CHUNK-C" in combined, combined[-1000:]
    assert "TRUNCATED-PARTIAL-SUMMARY" not in combined, combined[-1000:]


_TESTS = [
    test_openai_model_is_reasoning,
    test_google_model_ignores_sampling,
    test_openai_completion_kwargs,
    test_google_generation_kwargs,
    test_openai_reply_text,
    test_google_reply_text,
    test_get_llm_temperature_by_model,
    test_call_openai_llm_kwargs,
    test_call_google_llm_kwargs,
    test_google_openai_defaults,
    test_directive_extractor_google_openai_fallback,
    test_openai_model_is_reasoning_prefixed_names,
    test_cutoff_errors_not_rate_limits,
    test_google_reply_text_safety_and_int_finish,
    test_response_text_langchain_cutoff,
    test_hosted_llm_min_timeout,
    test_portkey_max_completion_tokens,
    test_call_portkey_llm,
    test_directive_extractor_portkey_fallback,
    test_summarizer_skips_cutoff_chunk,
]


def run_all_tests():
    """Run every test in _TESTS, fail-fast (used by run_all_tests.py)."""
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        test_fn()
        print("  PASS: %s" % test_fn.__name__)


if __name__ == "__main__":
    # Fail-fast: any exception (including AssertionError) propagates and
    # exits non-zero.
    run_all_tests()
    print()
    print("All %d tests passed." % len(_TESTS))
    print("OK")
