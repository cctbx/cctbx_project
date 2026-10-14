"""
Claude 5.x readiness for the Anthropic provider path.

Scope: pure logic and source scans.  No SDK imports (anthropic / langchain are
not needed), no keys, no network.  Reply objects are faked with
types.SimpleNamespace / small classes.

Covers:
  - anthropic_model_rejects_sampling_params: True for every Claude 5.x name
    (incl. provider-prefixed), still True for the 4.6+ family, False for older
    Claude, non-Claude, empty and None.
  - sanitize_llm_kwargs("anthropic", ...) drops temperature/top_p for a 5.x
    model and keeps max_tokens.
  - response_text: returns the text of a reply in every shape Claude 5.x can
    produce (thinking block first, SDK-like typed blocks, multiple text
    blocks, bare strings, None content) without raising.
  - The Anthropic decision/expensive defaults are the 5.x models, and both are
    rejected for sampling params.
  - Source scan: api_client.py and directive_extractor.py no longer read
    response.content[0].text.
  - anthropic_reply_text: returns the text after a thinking block; raises
    RuntimeError (naming stop_reason) when the reply has no text.
  - anthropic_model_rejects_sampling_params on Bedrock / Vertex / dated
    gateway forms.
  - tests/canary_expected.json defaults_fingerprint matches
    compute_defaults_fingerprint().
  - Source scan: the Anthropic SDK call paths in api_client.py and
    directive_extractor.py return anthropic_reply_text(response).
  - ANTHROPIC_MAX_TOKENS is at least 16000 (thinking shares the output
    budget on Claude 5.x; 2000 was cut off with stop_reason max_tokens).
  - get_llm_and_embeddings("anthropic") builds ChatAnthropic (faked via
    sys.modules) with max_tokens == ANTHROPIC_MAX_TOKENS; no temperature for
    a 5.x model, temperature for claude-sonnet-4-20250514.
  - api_client._call_anthropic_llm(max_tokens=2000) sends max_tokens >=
    16000 to messages.create (anthropic SDK faked via sys.modules).
  - anthropic_reply_text raises RuntimeError naming stop_reason='max_tokens'
    for a truncated reply even when its text is non-empty.
  - Source scan: directive_extractor.py's Anthropic messages.create call uses
    anthropic_request_max_tokens(...), not max_tokens=2000.
  - claude-sonnet-4-5 (bare and dated) is not in RETIRED_MODELS.
  - ANTHROPIC_MAX_TOKENS is at most 21333: the anthropic SDK refuses a
    non-streaming request whose max_tokens implies more than 600 s
    (3600 * max_tokens / 128000 > 600).
  - anthropic_request_max_tokens: floors at ANTHROPIC_MAX_TOKENS (None, 0,
    -5 and 2000 -> 16000), passes values in range through (18000.7 -> int
    18000), caps at 21333.
  - api_client._call_anthropic_llm sends 21333 for max_tokens=50000 and
    16000 for max_tokens=2000 (anthropic SDK faked via sys.modules).
  - anthropic_reply_text raises RuntimeError naming
    stop_reason='model_context_window_exceeded' for a reply with non-empty
    text; max_tokens still raises, end_turn still returns the text.
  - ANTHROPIC_MAX_TOKENS <= ANTHROPIC_NONSTREAMING_MAX_TOKENS <= 21333.

Fail-fast: the first failing test stops the run with a non-zero exit;
AssertionError is never swallowed.
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


# ---------------------------------------------------------------------------
# 1. anthropic_model_rejects_sampling_params
# ---------------------------------------------------------------------------

def test_rejects_claude5_names():
    from core.llm import anthropic_model_rejects_sampling_params as f
    for name in ["claude-sonnet-5-5", "claude-opus-5-5", "claude-opus-5",
                 "claude-sonnet-5", "claude-haiku-5-5", "claude-fable-5-1",
                 "anthropic.claude-opus-5-5"]:
        assert f(name) is True, "expected True for %r, got %r" % (name, f(name))


def test_rejects_claude46_plus_still():
    from core.llm import anthropic_model_rejects_sampling_params as f
    for name in ["claude-opus-4-7", "claude-sonnet-4-6", "claude-opus-4-8"]:
        assert f(name) is True, "expected True for %r, got %r" % (name, f(name))


def test_accepts_older_and_non_claude():
    from core.llm import anthropic_model_rejects_sampling_params as f
    for name in ["claude-sonnet-4-5", "claude-sonnet-4-5-20250929",
                 "claude-haiku-4-5-20251001", "claude-opus-4-20250514",
                 "claude-3-5-sonnet", "gpt-5", "", None]:
        assert f(name) is False, \
            "expected False for %r, got %r" % (name, f(name))


# ---------------------------------------------------------------------------
# 2. sanitize_llm_kwargs for a 5.x model
# ---------------------------------------------------------------------------

def test_sanitize_drops_sampling_for_claude5():
    from core.llm import sanitize_llm_kwargs
    src = {"model": "claude-sonnet-5-5", "max_tokens": 4096,
           "temperature": 0.2, "top_p": 0.9}
    out = sanitize_llm_kwargs("anthropic", src)
    assert "temperature" not in out, out
    assert "top_p" not in out, out
    assert out.get("max_tokens") == 4096, out
    assert out.get("model") == "claude-sonnet-5-5", out
    # Caller's dict is not mutated.
    assert src["temperature"] == 0.2 and src["top_p"] == 0.9, src


# ---------------------------------------------------------------------------
# 3. response_text
# ---------------------------------------------------------------------------

class _Block(object):
    """Anthropic-SDK-like typed content block."""
    def __init__(self, type, **kw):
        self.type = type
        for k, v in kw.items():
            setattr(self, k, v)


def test_response_text_plain_and_string_content():
    from core.llm import response_text
    assert response_text("hello") == "hello"
    msg = types.SimpleNamespace(content="answer")
    assert response_text(msg) == "answer"


def test_response_text_dict_blocks_thinking_first():
    from core.llm import response_text
    thinking = "SECRET-REASONING-should-not-appear"
    msg = types.SimpleNamespace(content=[
        {"type": "thinking", "thinking": thinking},
        {"type": "text", "text": "answer"},
    ])
    out = response_text(msg)
    assert out == "answer", repr(out)
    assert thinking not in out, repr(out)


def test_response_text_sdk_objects_thinking_first():
    from core.llm import response_text
    thinking = "SDK-THINKING-should-not-appear"
    msg = types.SimpleNamespace(content=[
        _Block("thinking", thinking=thinking, signature="sig"),
        _Block("text", text="answer"),
    ])
    out = response_text(msg)
    assert out == "answer", repr(out)
    assert thinking not in out, repr(out)


def test_response_text_two_text_blocks_concatenated():
    from core.llm import response_text
    msg = types.SimpleNamespace(content=[
        {"type": "text", "text": "first "},
        {"type": "text", "text": "second"},
    ])
    assert response_text(msg) == "first second", repr(response_text(msg))
    msg2 = types.SimpleNamespace(content=[
        _Block("text", text="A"), _Block("text", text="B")])
    assert response_text(msg2) == "AB", repr(response_text(msg2))


def test_response_text_bare_strings_and_none():
    from core.llm import response_text
    msg = types.SimpleNamespace(content=["foo", "bar"])
    assert response_text(msg) == "foobar", repr(response_text(msg))
    none_msg = types.SimpleNamespace(content=None)
    assert response_text(none_msg) == "", repr(response_text(none_msg))


# ---------------------------------------------------------------------------
# 4. Anthropic defaults
# ---------------------------------------------------------------------------

def test_anthropic_defaults_are_claude5_and_rejected():
    from core.llm import DECISION_MODEL_DEFAULTS
    from core.llm import EXPENSIVE_MODEL_DEFAULTS
    from core.llm import anthropic_model_rejects_sampling_params as f
    dec = DECISION_MODEL_DEFAULTS["anthropic"]
    exp = EXPENSIVE_MODEL_DEFAULTS["anthropic"]
    assert dec == "claude-sonnet-5-5", dec
    assert exp == "claude-opus-5-5", exp
    assert f(dec) is True, dec
    assert f(exp) is True, exp


# ---------------------------------------------------------------------------
# 5. Source scans
# ---------------------------------------------------------------------------

_CONTENT0_TEXT = re.compile(r"\bresponse\.content\s*\[\s*0\s*\]\s*\.text\b")


def _scan_no_content0_text(fname):
    path = os.path.join(_AGENT, fname)
    assert os.path.isfile(path), "missing source file: %s" % path
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    hits = []
    for i, line in enumerate(lines, 1):
        code = line.split("#", 1)[0]
        if _CONTENT0_TEXT.search(code):
            hits.append("%s:%d: %s" % (fname, i, line.strip()))
    assert not hits, "still reads response.content[0].text:\n  " + \
        "\n  ".join(hits)


def test_api_client_no_content0_text():
    _scan_no_content0_text("api_client.py")


def test_directive_extractor_no_content0_text():
    _scan_no_content0_text("directive_extractor.py")


# ---------------------------------------------------------------------------
# 6. anthropic_reply_text
# ---------------------------------------------------------------------------

def test_anthropic_reply_text_thinking_then_text():
    from core.llm import anthropic_reply_text
    msg = types.SimpleNamespace(stop_reason="end_turn", content=[
        _Block("thinking", thinking="private reasoning", signature="sig"),
        _Block("text", text="answer"),
    ])
    out = anthropic_reply_text(msg)
    assert out == "answer", repr(out)


def test_anthropic_reply_text_thinking_only_raises():
    from core.llm import anthropic_reply_text
    msg = types.SimpleNamespace(stop_reason="max_tokens", content=[
        _Block("thinking", thinking="ran out of budget", signature="sig"),
    ])
    raised = None
    try:
        anthropic_reply_text(msg)
    except RuntimeError as e:
        raised = e
    assert raised is not None, "expected RuntimeError for thinking-only reply"
    assert "stop_reason='max_tokens'" in str(raised), str(raised)
    # The stop_reason reported is the reply's own, not a fixed string.
    msg2 = types.SimpleNamespace(stop_reason="end_turn", content=[
        _Block("text", text="")])
    raised2 = None
    try:
        anthropic_reply_text(msg2)
    except RuntimeError as e:
        raised2 = e
    assert raised2 is not None, "expected RuntimeError for empty text block"
    assert "stop_reason='end_turn'" in str(raised2), str(raised2)


def test_anthropic_reply_text_blank_text_raises():
    from core.llm import anthropic_reply_text
    for blank in ["", "   \n\t "]:
        msg = types.SimpleNamespace(stop_reason="end_turn", content=[
            _Block("text", text=blank)])
        raised = None
        try:
            anthropic_reply_text(msg)
        except RuntimeError as e:
            raised = e
        assert raised is not None, \
            "expected RuntimeError for text block %r" % blank


# ---------------------------------------------------------------------------
# 7. Gateway / provider name forms
# ---------------------------------------------------------------------------

def test_rejects_gateway_name_forms():
    from core.llm import anthropic_model_rejects_sampling_params as f
    expected = [
        ("us.anthropic.claude-opus-4-7-v1:0", True),
        ("claude-opus-5-5@20260101", True),
        ("anthropic.claude-3-5-sonnet-20240620-v1:0", False),
        ("claude-3-7-sonnet-20250219", False),
        ("claude-opus-4-1@20250805", False),
        ("claude-opus-4-5-20251101", False),
    ]
    for name, want in expected:
        got = f(name)
        assert got is want, "%r: expected %r, got %r" % (name, want, got)


# ---------------------------------------------------------------------------
# 8. Canary fingerprint
# ---------------------------------------------------------------------------

def test_canary_fingerprint_matches_defaults():
    import json
    from core.llm import compute_defaults_fingerprint
    path = os.path.join(_HERE, "canary_expected.json")
    assert os.path.isfile(path), "missing %s" % path
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    want = data.get("defaults_fingerprint")
    got = compute_defaults_fingerprint()
    assert want == got, "canary %r != computed %r" % (want, got)


# ---------------------------------------------------------------------------
# 9. Source scan: SDK reply goes through anthropic_reply_text
# ---------------------------------------------------------------------------

def _anthropic_create_windows(fname):
    """Each 'client.messages.create(' call up to the next 'return' line."""
    path = os.path.join(_AGENT, fname)
    assert os.path.isfile(path), "missing source file: %s" % path
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    windows = []
    for i, line in enumerate(lines):
        if "client.messages.create(" in line.split("#", 1)[0]:
            for j in range(i, min(i + 30, len(lines))):
                if lines[j].strip().startswith("return"):
                    windows.append((i + 1, lines[j].strip()))
                    break
            else:
                windows.append((i + 1, None))
    return windows


def _check_reply_text_call(fname):
    windows = _anthropic_create_windows(fname)
    assert windows, "%s: no client.messages.create( call found" % fname
    for lineno, ret in windows:
        assert ret == "return anthropic_reply_text(response)", \
            "%s:%d: Anthropic SDK path returns %r" % (fname, lineno, ret)


def test_api_client_uses_anthropic_reply_text():
    _check_reply_text_call("api_client.py")


def test_directive_extractor_uses_anthropic_reply_text():
    _check_reply_text_call("directive_extractor.py")


# ---------------------------------------------------------------------------
# 10. Output budget: ANTHROPIC_MAX_TOKENS
# ---------------------------------------------------------------------------

def test_anthropic_max_tokens_at_least_16000():
    from core.llm import ANTHROPIC_MAX_TOKENS
    assert isinstance(ANTHROPIC_MAX_TOKENS, int), repr(ANTHROPIC_MAX_TOKENS)
    assert ANTHROPIC_MAX_TOKENS >= 16000, ANTHROPIC_MAX_TOKENS


class _FakeModules(object):
    """Temporarily install fake modules in sys.modules and clear env vars;
    restore both on exit."""
    def __init__(self, modules, clear_env=()):
        self.modules = modules
        self.clear_env = clear_env

    def __enter__(self):
        self._saved_mods = {}
        for name, mod in self.modules.items():
            self._saved_mods[name] = sys.modules.get(name, _MISSING)
            sys.modules[name] = mod
        self._saved_env = {}
        for name in self.clear_env:
            self._saved_env[name] = os.environ.pop(name, _MISSING)
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


_MISSING = object()

# Keys that would make the anthropic embeddings delegate try a real backend.
_EMBED_KEYS = ("OPENAI_API_KEY", "GOOGLE_API_KEY")


# ---------------------------------------------------------------------------
# 11. get_llm_and_embeddings("anthropic") ChatAnthropic kwargs
# ---------------------------------------------------------------------------

def _chat_anthropic_kwargs(model_name):
    from core.llm import get_llm_and_embeddings
    calls = []

    class ChatAnthropic(object):
        def __init__(self, **kw):
            calls.append(kw)

    fake = types.ModuleType("langchain_anthropic")
    fake.ChatAnthropic = ChatAnthropic
    with _FakeModules({"langchain_anthropic": fake}, clear_env=_EMBED_KEYS):
        llm, emb = get_llm_and_embeddings(
            provider="anthropic", llm_model_name=model_name, temperature=0.3)
    assert isinstance(llm, ChatAnthropic), type(llm)
    assert emb is None, "embeddings should be None with no keys: %r" % emb
    assert len(calls) == 1, calls
    return calls[0]


def test_chat_anthropic_max_tokens_and_temperature():
    from core.llm import ANTHROPIC_MAX_TOKENS
    kw = _chat_anthropic_kwargs("claude-sonnet-5-5")
    assert kw.get("model") == "claude-sonnet-5-5", kw
    assert kw.get("max_tokens") == ANTHROPIC_MAX_TOKENS, kw
    assert "temperature" not in kw, kw
    kw4 = _chat_anthropic_kwargs("claude-sonnet-4-20250514")
    assert kw4.get("max_tokens") == ANTHROPIC_MAX_TOKENS, kw4
    assert kw4.get("temperature") == 0.3, kw4


# ---------------------------------------------------------------------------
# 12. api_client._call_anthropic_llm sends a >= 16000 budget
# ---------------------------------------------------------------------------

def _ensure_libtbx_importable():
    """api_client imports libtbx.langchain.*; add the checkout root when
    libtbx is not already importable (plain python3 from a scratch dir)."""
    try:
        import libtbx  # noqa: F401
    except ImportError:
        sys.path.append(os.path.dirname(os.path.dirname(_ROOT)))


def test_call_anthropic_llm_max_tokens_floor():
    _ensure_libtbx_importable()
    from agent.api_client import _call_anthropic_llm
    calls = []

    class _Messages(object):
        def create(self, **kw):
            calls.append(kw)
            return types.SimpleNamespace(
                stop_reason="end_turn",
                content=[_Block("text", text="pong")])

    class Anthropic(object):
        def __init__(self, *a, **kw):
            self.messages = _Messages()

    fake = types.ModuleType("anthropic")
    fake.Anthropic = Anthropic
    with _FakeModules({"anthropic": fake}):
        out = _call_anthropic_llm(
            "hi", model="claude-sonnet-5-5", max_tokens=2000)
    assert out == "pong", repr(out)
    assert len(calls) == 1, calls
    assert calls[0].get("model") == "claude-sonnet-5-5", calls[0]
    assert calls[0].get("max_tokens", 0) >= 16000, calls[0]


# ---------------------------------------------------------------------------
# 13. anthropic_reply_text: truncated reply with text still raises
# ---------------------------------------------------------------------------

def test_anthropic_reply_text_max_tokens_with_text_raises():
    from core.llm import anthropic_reply_text
    msg = types.SimpleNamespace(stop_reason="max_tokens", content=[
        _Block("thinking", thinking="reasoning", signature="sig"),
        _Block("text", text="The answer is partially wri"),
    ])
    raised = None
    try:
        anthropic_reply_text(msg)
    except RuntimeError as e:
        raised = e
    assert raised is not None, \
        "expected RuntimeError for truncated reply with non-empty text"
    assert "stop_reason='max_tokens'" in str(raised), str(raised)
    ok = types.SimpleNamespace(stop_reason="end_turn", content=[
        _Block("text", text="complete answer")])
    assert anthropic_reply_text(ok) == "complete answer"


# ---------------------------------------------------------------------------
# 14. Source scan: directive_extractor Anthropic call budget
# ---------------------------------------------------------------------------

def test_directive_extractor_anthropic_max_tokens():
    path = os.path.join(_AGENT, "directive_extractor.py")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    windows = []
    for i, line in enumerate(lines):
        if "client.messages.create(" in line.split("#", 1)[0]:
            body = []
            for j in range(i, min(i + 30, len(lines))):
                body.append(lines[j].split("#", 1)[0])
                if lines[j].strip().startswith("return"):
                    break
            windows.append((i + 1, "\n".join(body)))
    assert windows, "no client.messages.create( call in directive_extractor"
    for lineno, body in windows:
        assert not re.search(r"max_tokens\s*=\s*2000\b", body), \
            "directive_extractor.py:%d still passes max_tokens=2000" % lineno
        assert re.search(
            r"max_tokens\s*=\s*anthropic_request_max_tokens\s*\(", body), \
            "directive_extractor.py:%d does not use " \
            "anthropic_request_max_tokens(...)" % lineno


# ---------------------------------------------------------------------------
# 15. Sonnet 4.5 is not retired
# ---------------------------------------------------------------------------

def test_sonnet_45_not_retired():
    from core.llm import RETIRED_MODELS
    for name in ["claude-sonnet-4-5", "claude-sonnet-4-5-20250929"]:
        assert name not in RETIRED_MODELS, \
            "%r is listed in RETIRED_MODELS" % name



# ---------------------------------------------------------------------------
# 16. Output budget upper bound (anthropic SDK non-streaming limit)
# ---------------------------------------------------------------------------

def test_anthropic_max_tokens_within_nonstreaming_limit():
    from core.llm import ANTHROPIC_MAX_TOKENS
    assert ANTHROPIC_MAX_TOKENS <= 21333, ANTHROPIC_MAX_TOKENS
    # The SDK's own rule: expected time 3600 * max_tokens / 128000 s must
    # not exceed 600 s for a non-streaming request.
    expected_s = 3600.0 * ANTHROPIC_MAX_TOKENS / 128000.0
    assert expected_s <= 600.0, (ANTHROPIC_MAX_TOKENS, expected_s)



# ---------------------------------------------------------------------------
# 17. anthropic_request_max_tokens: floor and non-streaming cap
# ---------------------------------------------------------------------------

def test_anthropic_request_max_tokens_floor_and_cap():
    from core.llm import anthropic_request_max_tokens as f
    for requested, want in [(None, 16000), (2000, 16000), (18000, 18000),
                            (21333, 21333), (50000, 21333), (0, 16000),
                            (-5, 16000), (18000.7, 18000)]:
        got = f(requested)
        assert got == want, "%r: expected %r, got %r" % (requested, want, got)
        assert type(got) is int, \
            "%r: expected int, got %r (%s)" % (requested, got, type(got))


# ---------------------------------------------------------------------------
# 18. api_client._call_anthropic_llm: capped and floored budget sent
# ---------------------------------------------------------------------------

def _anthropic_sent_max_tokens(requested):
    _ensure_libtbx_importable()
    from agent.api_client import _call_anthropic_llm
    calls = []

    class _Messages(object):
        def create(self, **kw):
            calls.append(kw)
            return types.SimpleNamespace(
                stop_reason="end_turn",
                content=[_Block("text", text="pong")])

    class Anthropic(object):
        def __init__(self, *a, **kw):
            self.messages = _Messages()

    fake = types.ModuleType("anthropic")
    fake.Anthropic = Anthropic
    with _FakeModules({"anthropic": fake}):
        out = _call_anthropic_llm(
            "hi", model="claude-sonnet-5-5", max_tokens=requested)
    assert out == "pong", repr(out)
    assert len(calls) == 1, calls
    return calls[0].get("max_tokens")


def test_call_anthropic_llm_max_tokens_cap():
    got = _anthropic_sent_max_tokens(50000)
    assert got == 21333, "max_tokens=50000 sent %r" % got
    got = _anthropic_sent_max_tokens(2000)
    assert got == 16000, "max_tokens=2000 sent %r" % got


# ---------------------------------------------------------------------------
# 19. anthropic_reply_text: context-window cut-off raises
# ---------------------------------------------------------------------------

def _reply_text_error(stop_reason, text):
    from core.llm import anthropic_reply_text
    msg = types.SimpleNamespace(stop_reason=stop_reason, content=[
        _Block("text", text=text)])
    try:
        out = anthropic_reply_text(msg)
    except RuntimeError as e:
        return None, e
    return out, None


def test_anthropic_reply_text_context_window_exceeded_raises():
    out, err = _reply_text_error("model_context_window_exceeded",
                                 "partial answer")
    assert err is not None, \
        "expected RuntimeError for model_context_window_exceeded, got %r" % out
    assert "stop_reason='model_context_window_exceeded'" in str(err), str(err)
    out, err = _reply_text_error("max_tokens", "partial answer")
    assert err is not None, "expected RuntimeError for max_tokens"
    assert "stop_reason='max_tokens'" in str(err), str(err)
    out, err = _reply_text_error("end_turn", "complete answer")
    assert err is None, str(err)
    assert out == "complete answer", repr(out)


# ---------------------------------------------------------------------------
# 20. Budget constants are ordered and within the SDK limit
# ---------------------------------------------------------------------------

def test_anthropic_budget_constants_ordered():
    from core.llm import ANTHROPIC_MAX_TOKENS
    from core.llm import ANTHROPIC_NONSTREAMING_MAX_TOKENS
    assert ANTHROPIC_MAX_TOKENS <= ANTHROPIC_NONSTREAMING_MAX_TOKENS, \
        (ANTHROPIC_MAX_TOKENS, ANTHROPIC_NONSTREAMING_MAX_TOKENS)
    assert ANTHROPIC_NONSTREAMING_MAX_TOKENS <= 21333, \
        ANTHROPIC_NONSTREAMING_MAX_TOKENS


_TESTS = [
    test_rejects_claude5_names,
    test_rejects_claude46_plus_still,
    test_accepts_older_and_non_claude,
    test_sanitize_drops_sampling_for_claude5,
    test_response_text_plain_and_string_content,
    test_response_text_dict_blocks_thinking_first,
    test_response_text_sdk_objects_thinking_first,
    test_response_text_two_text_blocks_concatenated,
    test_response_text_bare_strings_and_none,
    test_anthropic_defaults_are_claude5_and_rejected,
    test_api_client_no_content0_text,
    test_directive_extractor_no_content0_text,
    test_anthropic_reply_text_thinking_then_text,
    test_anthropic_reply_text_thinking_only_raises,
    test_anthropic_reply_text_blank_text_raises,
    test_rejects_gateway_name_forms,
    test_canary_fingerprint_matches_defaults,
    test_api_client_uses_anthropic_reply_text,
    test_directive_extractor_uses_anthropic_reply_text,
    test_anthropic_max_tokens_at_least_16000,
    test_chat_anthropic_max_tokens_and_temperature,
    test_call_anthropic_llm_max_tokens_floor,
    test_anthropic_reply_text_max_tokens_with_text_raises,
    test_directive_extractor_anthropic_max_tokens,
    test_sonnet_45_not_retired,
    test_anthropic_max_tokens_within_nonstreaming_limit,
    test_anthropic_request_max_tokens_floor_and_cap,
    test_call_anthropic_llm_max_tokens_cap,
    test_anthropic_reply_text_context_window_exceeded_raises,
    test_anthropic_budget_constants_ordered,
]


if __name__ == "__main__":
    # Fail-fast: any exception (including AssertionError) propagates and
    # exits non-zero.
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        test_fn()
        print("  PASS: %s" % test_fn.__name__)
    print()
    print("All %d tests passed." % len(_TESTS))
    print("OK")
