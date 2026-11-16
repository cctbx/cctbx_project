"""
Provider cut-off and content-filter signals, and log-summary coverage.

Scope: pure logic.  No keys, no network.  Replies are faked with
types.SimpleNamespace; the LLM given to summarize_log_text is a small
langchain_core Runnable that never calls a provider.

Covers:
  - response_text raises RuntimeError for an Ollama reply with
    response_metadata done_reason "length".
  - response_text raises a distinct RuntimeError ("blocked by the provider's
    content filter") for finish_reason "content_filter", given as a string
    (any case) or as an object whose .name is that.
  - Control: done_reason "stop" (and no metadata) does not raise.
  - openai_reply_text raises RuntimeError naming the content filter for
    finish_reason "content_filter", even when the reply text is non-empty.
    Control: finish_reason "stop" returns the text.
  - summary_coverage_note: "" when nothing was skipped; otherwise a note
    that names the summarized, total and skipped counts.
  - summarize_log_text end to end (fake LLM, provider "ollama", 2+ chunks):
    map replies with finish_reason "length" are counted as skipped parts;
    the summary ends with the coverage note; parts_total / parts_skipped are
    returned.  With no cut-offs there is no note and parts_skipped == 0.
  - summarize_log_text timeout path: a batch whose map call outlasts the
    timeout counts every chunk of that batch as skipped.
  - Agent direct Ollama calls (fake openai.OpenAI, no network):
    api_client._call_ollama_llm (and call_llm_simple(provider="ollama"))
    raises RuntimeError "cut off" for finish_reason "length"; control
    finish_reason "stop" returns the text.  directive_extractor's
    _call_llm_fallback Ollama branch returns None (its failure value) and
    logs "Ollama API call failed ... cut off" instead of returning the
    partial text.

Fail-fast: the first failing test stops the run with a non-zero exit;
AssertionError is never swallowed.
"""
import asyncio
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain


def _llm_module():
    try:
        from libtbx.langchain.core import llm
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from core import llm
    return llm


def _summarizer_module():
    try:
        from libtbx.langchain.analysis import summarizer
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from analysis import summarizer
    return summarizer


def _reply(content="some text", **meta):
    return types.SimpleNamespace(content=content, response_metadata=dict(meta))


def _expect_runtime_error(fn, *args):
    try:
        fn(*args)
    except RuntimeError as e:
        return str(e)
    raise AssertionError("expected RuntimeError from %s" % fn.__name__)


# ---------------------------------------------------------------------------
# 1. response_text: Ollama done_reason and content_filter
# ---------------------------------------------------------------------------

def test_ollama_done_reason_length_raises():
    llm = _llm_module()
    msg = _expect_runtime_error(llm.response_text,
                                _reply("partial answer", done_reason="length"))
    assert "cut off" in msg, msg
    assert "content filter" not in msg, msg


def test_ollama_done_reason_stop_does_not_raise():
    # Control: passes on base too.
    llm = _llm_module()
    assert llm.response_text(_reply("full answer", done_reason="stop")) \
        == "full answer"
    assert llm.response_text(_reply("plain")) == "plain"


def test_content_filter_string_raises_distinct_error():
    llm = _llm_module()
    for value in ["content_filter", "CONTENT_FILTER", "Content_Filter"]:
        msg = _expect_runtime_error(
            llm.response_text, _reply("withheld", finish_reason=value))
        assert "blocked by the provider's content filter" in msg, (value, msg)
        assert "cut off" not in msg, (value, msg)


def test_content_filter_enum_like_raises_distinct_error():
    llm = _llm_module()
    for name in ["CONTENT_FILTER", "content_filter"]:
        finish = types.SimpleNamespace(name=name)
        msg = _expect_runtime_error(
            llm.response_text, _reply("withheld", finish_reason=finish))
        assert "blocked by the provider's content filter" in msg, (name, msg)


# ---------------------------------------------------------------------------
# 2. openai_reply_text: content_filter
# ---------------------------------------------------------------------------

def _openai_completion(text, finish):
    choice = types.SimpleNamespace(
        message=types.SimpleNamespace(content=text), finish_reason=finish)
    return types.SimpleNamespace(choices=[choice])


def test_openai_reply_text_content_filter_raises_with_text():
    llm = _llm_module()
    msg = _expect_runtime_error(
        llm.openai_reply_text,
        _openai_completion("non-empty partial text", "content_filter"))
    assert "content filter" in msg, msg


def test_openai_reply_text_stop_returns_text():
    # Control: passes on base too.
    llm = _llm_module()
    assert llm.openai_reply_text(_openai_completion("hello", "stop")) == "hello"


# ---------------------------------------------------------------------------
# 3. summary_coverage_note
# ---------------------------------------------------------------------------

def test_summary_coverage_note():
    s = _summarizer_module()
    assert s.summary_coverage_note(5, 0) == ""
    note = s.summary_coverage_note(5, 2)
    assert note, note
    assert "3 of 5" in note, note
    assert "2 part" in note, note


# ---------------------------------------------------------------------------
# 4. summarize_log_text end to end with a fake LLM
# ---------------------------------------------------------------------------

CUT_MARK = "ZZCUTOFFMARKZZ"
SLOW_MARK = "ZZSLOWMARKZZ"


def _make_fake_llm(slow_seconds=2.0):
    from langchain_core.runnables import Runnable

    class FakeLLM(Runnable):
        """Map calls (prompt value in) reply per chunk; reduce (messages in)
        replies with a fixed final summary."""

        def _reply_for(self, input):
            if isinstance(input, (list, tuple)):
                return _reply("FINAL SUMMARY", finish_reason="stop")
            text = input.to_string() if hasattr(input, "to_string") \
                else str(input)
            if CUT_MARK in text:
                return _reply("chunk summary, cut", finish_reason="length")
            return _reply("chunk summary", finish_reason="stop")

        def invoke(self, input, config=None, **kwargs):
            return self._reply_for(input)

        async def ainvoke(self, input, config=None, **kwargs):
            text = input.to_string() if hasattr(input, "to_string") else ""
            if SLOW_MARK in text:
                await asyncio.sleep(slow_seconds)
            return self._reply_for(input)

    return FakeLLM()


def _section(tag, n_lines=1100):
    # ~70 chars per line -> ~77000 chars per section, so with the ollama
    # chunk size (80000, overlap 8000) each section is about one chunk.
    lines = ["phenix.refine cycle %s line %05d R-work 0.2 R-free 0.25 ....."
             % (tag, k) for k in range(n_lines)]
    mid = n_lines // 2
    lines[mid] = lines[mid] + " " + tag
    return "\n".join(lines)


def _run_summary(text, **kw):
    s = _summarizer_module()
    old = s.DEBUG_SUMMARIZER
    s.DEBUG_SUMMARIZER = False
    try:
        args = dict(llm=_make_fake_llm(), timeout=60, batch_size=3,
                    use_throttling=False, provider="ollama",
                    program_name="phenix.refine")
        args.update(kw)
        return asyncio.run(s.summarize_log_text(text, **args))
    finally:
        s.DEBUG_SUMMARIZER = old


def _chunks(text):
    s = _summarizer_module()
    old = s.DEBUG_SUMMARIZER
    s.DEBUG_SUMMARIZER = False
    try:
        return s._custom_log_chunker(text, provider="ollama")
    finally:
        s.DEBUG_SUMMARIZER = old


def test_summarize_log_text_counts_cut_off_chunks():
    text = "\n\n".join([_section("AAAA"), _section(CUT_MARK),
                        _section("BBBB"), _section(CUT_MARK + "2")])
    docs = _chunks(text)
    n_total = len(docs)
    n_cut = sum(1 for d in docs if CUT_MARK in d.page_content)
    assert n_total >= 2, n_total
    assert 1 <= n_cut < n_total, (n_cut, n_total)
    r = _run_summary(text)
    assert r.error is None, r.error
    assert r.parts_total == n_total, (r.parts_total, n_total)
    assert r.parts_skipped == n_cut, (r.parts_skipped, n_cut)
    s = _summarizer_module()
    note = s.summary_coverage_note(n_total, n_cut)
    assert r.log_summary.startswith("FINAL SUMMARY"), r.log_summary[:200]
    assert r.log_summary.endswith(note), r.log_summary[-400:]
    assert "%d of %d" % (n_total - n_cut, n_total) in r.log_summary, \
        r.log_summary[-400:]


def test_summarize_log_text_no_cut_off_no_note():
    text = "\n\n".join([_section("AAAA"), _section("BBBB"), _section("CCCC")])
    n_total = len(_chunks(text))
    assert n_total >= 2, n_total
    r = _run_summary(text)
    assert r.error is None, r.error
    assert r.parts_total == n_total, (r.parts_total, n_total)
    assert r.parts_skipped == 0, r.parts_skipped
    assert r.log_summary == "FINAL SUMMARY", r.log_summary[-400:]
    assert "partial coverage" not in r.log_summary.lower()


def test_summarize_log_text_timed_out_batch_counts_whole_batch():
    # batch_size=2: the first batch holds the slow chunk and times out, so
    # both of its chunks are skipped; the later batch is summarized.
    text = "\n\n".join([_section(SLOW_MARK), _section("AAAA"),
                        _section("BBBB")])
    docs = _chunks(text)
    n_total = len(docs)
    assert n_total >= 3, n_total
    assert SLOW_MARK in docs[0].page_content
    assert not any(SLOW_MARK in d.page_content for d in docs[2:])
    r = _run_summary(text, batch_size=2, timeout=0.3)
    assert r.error is None, r.error
    assert r.parts_total == n_total, (r.parts_total, n_total)
    assert r.parts_skipped == 2, r.parts_skipped
    assert "%d of %d" % (n_total - 2, n_total) in r.log_summary, \
        r.log_summary[-400:]


# ---------------------------------------------------------------------------
# 5. Agent direct Ollama calls (fake openai.OpenAI; no network)
# ---------------------------------------------------------------------------

def _api_client_module():
    try:
        from libtbx.langchain.agent import api_client
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from agent import api_client
    return api_client


def _directive_extractor_module():
    try:
        from libtbx.langchain.agent import directive_extractor
    except ImportError:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        from agent import directive_extractor
    return directive_extractor


class _FakeOpenAIClient(object):
    """Stands in for openai.OpenAI; create() returns a canned completion."""
    calls = []
    reply_text = "partial"
    reply_finish = "length"

    def __init__(self, base_url=None, api_key=None, **kwargs):
        _FakeOpenAIClient.calls.append(dict(base_url=base_url,
                                            api_key=api_key))
        self.chat = types.SimpleNamespace(
            completions=types.SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        return _openai_completion(_FakeOpenAIClient.reply_text,
                                  _FakeOpenAIClient.reply_finish)


def _with_fake_openai(text, finish, fn):
    """Run fn() with openai.OpenAI replaced by _FakeOpenAIClient and
    OLLAMA_BASE_URL pointed at an unused local port; restore both."""
    import openai
    old_cls = openai.OpenAI
    old_url = os.environ.get("OLLAMA_BASE_URL")
    _FakeOpenAIClient.calls = []
    _FakeOpenAIClient.reply_text = text
    _FakeOpenAIClient.reply_finish = finish
    openai.OpenAI = _FakeOpenAIClient
    os.environ["OLLAMA_BASE_URL"] = "http://127.0.0.1:9"
    try:
        result = fn()
    finally:
        openai.OpenAI = old_cls
        if old_url is None:
            os.environ.pop("OLLAMA_BASE_URL", None)
        else:
            os.environ["OLLAMA_BASE_URL"] = old_url
    assert len(_FakeOpenAIClient.calls) == 1, _FakeOpenAIClient.calls
    return result


PARTIAL = "partial directive text ZZPARTIALZZ"


def test_call_ollama_llm_length_raises():
    ac = _api_client_module()
    for fn in (lambda: ac._call_ollama_llm("prompt", model="fake-model"),
               lambda: ac.call_llm_simple("prompt", provider="ollama",
                                          model="fake-model")):
        try:
            got = _with_fake_openai(PARTIAL, "length", fn)
        except RuntimeError as e:
            assert "cut off" in str(e), str(e)
            assert PARTIAL not in str(e), str(e)
        else:
            raise AssertionError(
                "expected RuntimeError, got partial text %r" % (got,))


def test_call_ollama_llm_stop_returns_text():
    # Control: passes on base too.
    ac = _api_client_module()
    got = _with_fake_openai(
        "full reply", "stop",
        lambda: ac._call_ollama_llm("prompt", model="fake-model"))
    assert got == "full reply", got


def test_directive_fallback_ollama_length_returns_none():
    # Narrow entry point: _call_llm_fallback is the function holding the
    # Ollama branch.  The public extract path reaches it only when
    # api_client fails to import, so call it directly.  On any exception
    # the branch logs "DIRECTIVES: Ollama API call failed - ..." and
    # returns None; that is the observable failure result.
    de = _directive_extractor_module()
    logged = []
    got = _with_fake_openai(
        PARTIAL, "length",
        lambda: de._call_llm_fallback("prompt", "ollama", "fake-model",
                                      logged.append))
    assert got is None, "partial text returned as success: %r" % (got,)
    fail = [m for m in logged if "Ollama API call failed" in m]
    assert len(fail) == 1, logged
    assert "cut off" in fail[0], fail[0]


def test_directive_fallback_ollama_stop_returns_text():
    # Control: passes on base too.
    de = _directive_extractor_module()
    logged = []
    got = _with_fake_openai(
        "full reply", "stop",
        lambda: de._call_llm_fallback("prompt", "ollama", "fake-model",
                                      logged.append))
    assert got == "full reply", (got, logged)
    assert not any("failed" in m for m in logged), logged


_TESTS = [
    test_ollama_done_reason_length_raises,
    test_ollama_done_reason_stop_does_not_raise,
    test_content_filter_string_raises_distinct_error,
    test_content_filter_enum_like_raises_distinct_error,
    test_openai_reply_text_content_filter_raises_with_text,
    test_openai_reply_text_stop_returns_text,
    test_summary_coverage_note,
    test_summarize_log_text_counts_cut_off_chunks,
    test_summarize_log_text_no_cut_off_no_note,
    test_summarize_log_text_timed_out_batch_counts_whole_batch,
    test_call_ollama_llm_length_raises,
    test_call_ollama_llm_stop_returns_text,
    test_directive_fallback_ollama_length_returns_none,
    test_directive_fallback_ollama_stop_returns_text,
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
