"""
Provider refusal and safety-block signals.

Scope: pure logic.  No keys, no network.  Replies are faked with
types.SimpleNamespace.

Covers:
  - response_text (LangChain path) raises RuntimeError for Anthropic
    response_metadata stop_reason "refusal"; the message says "declined"
    and carries no digits.
  - response_text raises a "safety filter" RuntimeError for each Gemini
    blocked finish reason (SAFETY, RECITATION, PROHIBITED_CONTENT,
    BLOCKLIST, SPII), given as a plain string (any case), as
    "FinishReason.X", or as an object whose .name is X.
  - anthropic_reply_text raises for stop_reason "refusal" even when the
    reply has non-empty text.
  - google_reply_text raises "blocked by a safety filter" for each blocked
    finish reason (enum-like) even when response.text is non-empty.
  - Controls (pass on base too): stop_reason "end_turn" and finish_reason
    "STOP" return the text; google_reply_text with STOP returns the text;
    MAX_TOKENS still raises "cut off".
  - No refusal / block message contains "429" or "503", and none is
    classified as a rate limit by rate_limit_handler.is_rate_limit_error.
  - Empty text: google_reply_text with a blocked reason and anthropic_reply_text
    with stop_reason "refusal" raise the block / "declined" message, not the
    generic "no text" one.
  - Negative controls (pass on base too): finish_reason "UNKNOWN_12",
    "IMAGE_SAFETY" and "LANGUAGE" with text return the text through
    response_text and google_reply_text; "content_filter" still gives the
    content-filter message, not the safety-filter one.
  - Real SDK objects, built offline (each skipped with a note if its package
    is missing): google.genai.types.FinishReason.SAFETY with
    google_reply_text, and anthropic.types.Message with stop_reason
    "refusal" with anthropic_reply_text.

Fail-fast: the first failing test stops the run with a non-zero exit;
AssertionError is never swallowed.
"""
import os
import re
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain

BLOCKED = ["SAFETY", "RECITATION", "PROHIBITED_CONTENT", "BLOCKLIST", "SPII"]


def _llm_module():
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from core import llm
    return llm


def _rate_limit_module():
    try:
        from libtbx.langchain.agent import rate_limit_handler
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from agent import rate_limit_handler
    return rate_limit_handler


def _reply(content="some text", **meta):
    return types.SimpleNamespace(content=content, response_metadata=dict(meta))


def _anthropic_reply(text, stop_reason):
    block = types.SimpleNamespace(type="text", text=text)
    return types.SimpleNamespace(stop_reason=stop_reason, content=[block])


def _google_reply(text, finish):
    cand = types.SimpleNamespace(finish_reason=finish)
    return types.SimpleNamespace(candidates=[cand], text=text)


def _expect_runtime_error(fn, *args):
    try:
        got = fn(*args)
    except RuntimeError as e:
        return e
    raise AssertionError("expected RuntimeError from %s, got %r"
                         % (fn.__name__, got))


def _finish_forms(name):
    """A finish reason as plain string (two cases), "FinishReason.X" and an
    enum-like object with .name."""
    return [name, name.lower(), "FinishReason." + name,
            types.SimpleNamespace(name=name)]


# ---------------------------------------------------------------------------
# 1. LangChain path: Anthropic refusal
# ---------------------------------------------------------------------------

def test_langchain_refusal_raises():
    llm = _llm_module()
    err = _expect_runtime_error(
        llm.response_text, _reply("I can't help with that.",
                                  stop_reason="refusal"))
    msg = str(err)
    assert "declined" in msg, msg
    assert re.search(r"\d", msg) is None, msg
    assert "cut off" not in msg, msg


# ---------------------------------------------------------------------------
# 2. LangChain path: Gemini blocked finish reasons
# ---------------------------------------------------------------------------

def test_langchain_gemini_blocked_reasons_raise():
    llm = _llm_module()
    for name in BLOCKED:
        for finish in _finish_forms(name):
            err = _expect_runtime_error(
                llm.response_text, _reply("partial", finish_reason=finish))
            msg = str(err)
            assert "safety filter" in msg, (name, finish, msg)
            assert name in msg, (name, finish, msg)
            assert "cut off" not in msg, (name, finish, msg)


# ---------------------------------------------------------------------------
# 3. anthropic_reply_text: refusal with text
# ---------------------------------------------------------------------------

def test_anthropic_reply_text_refusal_with_text_raises():
    llm = _llm_module()
    err = _expect_runtime_error(
        llm.anthropic_reply_text,
        _anthropic_reply("I won't do that, but here is something.",
                         "refusal"))
    assert "declined" in str(err), str(err)


# ---------------------------------------------------------------------------
# 4. google_reply_text: blocked reasons with text
# ---------------------------------------------------------------------------

def test_google_reply_text_blocked_with_text_raises():
    llm = _llm_module()
    for name in BLOCKED:
        err = _expect_runtime_error(
            llm.google_reply_text,
            _google_reply("non-empty partial text",
                          types.SimpleNamespace(name=name)))
        msg = str(err)
        assert "safety filter" in msg, (name, msg)
        assert name in msg, (name, msg)


# ---------------------------------------------------------------------------
# 5. Controls (pass on base too)
# ---------------------------------------------------------------------------

def test_langchain_normal_stop_returns_text():
    # Control: passes on base too.
    llm = _llm_module()
    assert llm.response_text(_reply("full", stop_reason="end_turn")) == "full"
    assert llm.response_text(_reply("full", finish_reason="STOP")) == "full"
    assert llm.response_text(
        _reply("full", finish_reason=types.SimpleNamespace(name="STOP"))) \
        == "full"
    assert llm.anthropic_reply_text(_anthropic_reply("full", "end_turn")) \
        == "full"


def test_google_reply_text_stop_returns_text():
    # Control: passes on base too.
    llm = _llm_module()
    assert llm.google_reply_text(
        _google_reply("hello", types.SimpleNamespace(name="STOP"))) == "hello"
    assert llm.google_reply_text(_google_reply("hello", "STOP")) == "hello"


def test_max_tokens_still_cut_off():
    # Control: passes on base too.
    llm = _llm_module()
    msg = str(_expect_runtime_error(
        llm.response_text, _reply("partial", finish_reason="MAX_TOKENS")))
    assert "cut off" in msg, msg
    msg = str(_expect_runtime_error(
        llm.google_reply_text,
        _google_reply("partial", types.SimpleNamespace(name="MAX_TOKENS"))))
    assert "cut off" in msg, msg
    msg = str(_expect_runtime_error(
        llm.anthropic_reply_text, _anthropic_reply("partial", "max_tokens")))
    assert "cut off" in msg, msg


# ---------------------------------------------------------------------------
# 6. Not rate limits
# ---------------------------------------------------------------------------

def test_refusal_errors_not_rate_limits():
    llm = _llm_module()
    is_rate_limit_error = _rate_limit_module().is_rate_limit_error
    errors = [
        _expect_runtime_error(llm.response_text,
                              _reply("x", stop_reason="refusal")),
        _expect_runtime_error(llm.anthropic_reply_text,
                              _anthropic_reply("x" * 1429, "refusal")),
    ]
    for name in BLOCKED:
        errors.append(_expect_runtime_error(
            llm.response_text, _reply("x", finish_reason=name)))
        errors.append(_expect_runtime_error(
            llm.google_reply_text,
            _google_reply("x" * 1503, types.SimpleNamespace(name=name))))
    for err in errors:
        msg = str(err)
        assert "429" not in msg and "503" not in msg, msg
        assert is_rate_limit_error(err) is False, msg


# ---------------------------------------------------------------------------
# 7. Empty text: the block / refusal message wins over "no text"
# ---------------------------------------------------------------------------

def test_google_reply_text_blocked_empty_text_raises_block():
    llm = _llm_module()
    for text in ("", None):
        err = _expect_runtime_error(
            llm.google_reply_text,
            _google_reply(text, types.SimpleNamespace(name="SAFETY")))
        msg = str(err)
        assert "blocked by a safety filter" in msg, (text, msg)
        assert "SAFETY" in msg, (text, msg)
        assert "no text" not in msg, (text, msg)


def test_anthropic_reply_text_refusal_empty_text_raises_declined():
    llm = _llm_module()
    for reply in (_anthropic_reply("", "refusal"),
                  types.SimpleNamespace(stop_reason="refusal", content=[])):
        err = _expect_runtime_error(llm.anthropic_reply_text, reply)
        msg = str(err)
        assert "declined" in msg, msg
        assert "no text" not in msg, msg


# ---------------------------------------------------------------------------
# 8. Negative controls (pass on base too)
# ---------------------------------------------------------------------------

NOT_BLOCKED = ["UNKNOWN_12", "IMAGE_SAFETY", "LANGUAGE"]


def test_unblocked_finish_reasons_return_text():
    # Control: passes on base too.  Names that merely contain "SAFETY" or are
    # unknown are not blocked.
    llm = _llm_module()
    for name in NOT_BLOCKED:
        for finish in _finish_forms(name):
            got = llm.response_text(_reply("full text", finish_reason=finish))
            assert got == "full text", (name, finish, got)
            got = llm.google_reply_text(_google_reply("full text", finish))
            assert got == "full text", (name, finish, got)


def test_content_filter_keeps_its_message():
    # Control: passes on base too.
    llm = _llm_module()
    for finish in ("content_filter", "CONTENT_FILTER",
                   types.SimpleNamespace(name="CONTENT_FILTER")):
        msg = str(_expect_runtime_error(
            llm.response_text, _reply("partial", finish_reason=finish)))
        assert "content filter" in msg, (finish, msg)
        assert "safety filter" not in msg, (finish, msg)
    choice = types.SimpleNamespace(
        message=types.SimpleNamespace(content="partial"),
        finish_reason="content_filter")
    msg = str(_expect_runtime_error(
        llm.openai_reply_text, types.SimpleNamespace(choices=[choice])))
    assert "content filter" in msg, msg
    assert "safety filter" not in msg, msg


# ---------------------------------------------------------------------------
# 9. Real SDK objects, built offline (skipped if the package is missing)
# ---------------------------------------------------------------------------

def test_google_genai_real_safety_enum():
    try:
        from google.genai import types as gtypes
    except ImportError as e:
        print("    SKIP: google.genai not importable (%s)" % e)
        return
    llm = _llm_module()
    safety = gtypes.FinishReason.SAFETY
    empty = gtypes.GenerateContentResponse(
        candidates=[gtypes.Candidate(finish_reason=safety)])
    with_text = gtypes.GenerateContentResponse(candidates=[gtypes.Candidate(
        finish_reason=safety,
        content=gtypes.Content(role="model",
                               parts=[gtypes.Part(text="partial")]))])
    assert with_text.text == "partial", with_text.text
    for resp in (empty, with_text):
        msg = str(_expect_runtime_error(llm.google_reply_text, resp))
        assert "blocked by a safety filter" in msg, msg
        assert "SAFETY" in msg, msg
        assert "no text" not in msg, msg


def test_anthropic_real_message_refusal():
    try:
        from anthropic.types import Message
    except ImportError as e:
        print("    SKIP: anthropic not importable (%s)" % e)
        return
    llm = _llm_module()
    for content in ([], [{"type": "text", "text": ""}],
                    [{"type": "text", "text": "I can't help with that."}]):
        msg_obj = Message.model_validate({
            "id": "msg_offline_test", "type": "message", "role": "assistant",
            "model": "claude-opus-5-5", "content": content,
            "stop_reason": "refusal", "stop_sequence": None,
            "usage": {"input_tokens": 1, "output_tokens": 0}})
        assert msg_obj.stop_reason == "refusal", msg_obj.stop_reason
        msg = str(_expect_runtime_error(llm.anthropic_reply_text, msg_obj))
        assert "declined" in msg, (content, msg)
        assert "no text" not in msg, (content, msg)


_TESTS = [
    test_langchain_refusal_raises,
    test_langchain_gemini_blocked_reasons_raise,
    test_anthropic_reply_text_refusal_with_text_raises,
    test_google_reply_text_blocked_with_text_raises,
    test_langchain_normal_stop_returns_text,
    test_google_reply_text_stop_returns_text,
    test_max_tokens_still_cut_off,
    test_refusal_errors_not_rate_limits,
    test_google_reply_text_blocked_empty_text_raises_block,
    test_anthropic_reply_text_refusal_empty_text_raises_declined,
    test_unblocked_finish_reasons_return_text,
    test_content_filter_keeps_its_message,
    test_google_genai_real_safety_enum,
    test_anthropic_real_message_refusal,
]


def run_all_tests():
    # Fail-fast: any exception (including AssertionError) propagates and
    # exits non-zero.
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        test_fn()
        print("  PASS: %s" % test_fn.__name__)
    print()
    print("All %d tests passed." % len(_TESTS))
    print("OK")


if __name__ == "__main__":
    run_all_tests()
