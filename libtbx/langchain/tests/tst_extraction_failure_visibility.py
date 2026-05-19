"""
Sandbox tests for v118 Section E: LLM extraction-failure diagnostics.

Section E adds observability for silent LLM failures at two sites:

  E1+E2 — agent/directive_extractor.py extract_directives()
    On LLM exception:
      * log("CRITICAL_FALLBACK: ...") via existing log_func callback
        (developer visibility — Gemini Q3)
      * stderr "[DIRECTIVE_EXTRACTION_FAILED] ..." marker
        (operator visibility — Gemini Q1)

  E3+E4 — agent/advice_preprocessor.py preprocess_advice()
    On LLM exception:
      * print("CRITICAL_FALLBACK: ...", file=out) parallel message
      * stderr "[ADVICE_PREPROCESSING_FAILED] ..." marker
        (Gemini Q4 — same failure class as extractor)

Both diagnostics are observation-only.  Zero behavioral impact:
the underlying fallback paths (return {} for extractor, return
raw_advice for preprocessor) are preserved exactly as before.

Tests K_E1-K_E6 cover the extractor site.
Tests K_E7-K_E10 cover the preprocessor site.
K_E11 is a shared structural test.
"""

import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

DIRECTIVE_EXTRACTOR_PATH = os.path.join(
    ROOT, "agent", "directive_extractor.py")
ADVICE_PREPROCESSOR_PATH = os.path.join(
    ROOT, "agent", "advice_preprocessor.py")


# --------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------

def _import_module(path, module_name):
    """Dynamically import a module by path."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _capture_stderr():
    """Replace sys.stderr with a StringIO; return (stream, restore_fn)."""
    saved = sys.stderr
    captured = io.StringIO()
    sys.stderr = captured

    def restore():
        sys.stderr = saved
    return captured, restore


def _capture_log():
    """Return (log_func, captured_messages)."""
    msgs = []
    def log_func(m):
        msgs.append(m)
    return log_func, msgs


# --------------------------------------------------------------------
# Site 1: directive extractor (K_E1 - K_E6)
# --------------------------------------------------------------------

def test_k_e1_extractor_exception_emits_stderr_marker():
    """K_E1: extractor LLM exception → stderr contains marker."""
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de1")

    # Force _call_llm to raise — this triggers the exception handler
    # path at line 834.
    def fake_call_llm(prompt, provider, model, log):
        raise RuntimeError("simulated google api outage")
    mod._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    captured, restore = _capture_stderr()
    try:
        result = mod.extract_directives(
            "refine the model and stop",
            provider="google",
            log_func=log_func)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "[DIRECTIVE_EXTRACTION_FAILED]" in err_output, (
        "K_E1 failed: stderr should contain "
        "[DIRECTIVE_EXTRACTION_FAILED] marker. Got: %r" % err_output)
    assert result == {}, (
        "K_E1 failed: return value should still be {} on failure, "
        "got %r" % result)

    print("  PASS: K_E1 (extractor exception → stderr marker)")


def test_k_e2_extractor_success_no_stderr_marker():
    """K_E2: extractor LLM success → stderr does NOT contain marker."""
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de2")

    # Return valid JSON; bypass actual LLM
    def fake_call_llm(prompt, provider, model, log):
        return '{"stop_conditions": {"after_program": "phenix.refine"}}'
    mod._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    captured, restore = _capture_stderr()
    try:
        mod.extract_directives(
            "refine the model and stop",
            provider="google",
            log_func=log_func)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "[DIRECTIVE_EXTRACTION_FAILED]" not in err_output, (
        "K_E2 failed: stderr should NOT contain marker on success. "
        "Got: %r" % err_output)

    print("  PASS: K_E2 (extractor success → no stderr marker)")


def test_k_e3_extractor_empty_response_no_stderr_marker():
    """K_E3: extractor LLM returns empty → stderr does NOT contain marker.

    An empty LLM response is semantically valid — the model judged the
    advice non-actionable.  This is NOT a failure and should not fire
    the diagnostic (per Gemini Q2 — no alert fatigue).
    """
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de3")

    # Empty response — different code path (line ~648 "No response from
    # LLM" returns {} without raising).
    def fake_call_llm(prompt, provider, model, log):
        return ""
    mod._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    captured, restore = _capture_stderr()
    try:
        result = mod.extract_directives(
            "refine the model and stop",
            provider="google",
            log_func=log_func)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "[DIRECTIVE_EXTRACTION_FAILED]" not in err_output, (
        "K_E3 failed: empty response is not failure; stderr should "
        "NOT contain marker. Got: %r" % err_output)

    print("  PASS: K_E3 (extractor empty response → no stderr marker)")


def test_k_e4_extractor_marker_includes_provider_and_error():
    """K_E4: extractor marker line contains provider name + error msg."""
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de4")

    err_text = "specific test error message 12345"
    def fake_call_llm(prompt, provider, model, log):
        raise ValueError(err_text)
    mod._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    captured, restore = _capture_stderr()
    try:
        mod.extract_directives(
            "refine the model and stop",
            provider="openai",
            log_func=log_func)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "provider=openai" in err_output, (
        "K_E4 failed: marker should include provider name. "
        "Got: %r" % err_output)
    assert err_text in err_output, (
        "K_E4 failed: marker should include exception message. "
        "Got: %r" % err_output)
    assert "Falling back to rules-only resolver" in err_output, (
        "K_E4 failed: marker should explain fallback. "
        "Got: %r" % err_output)

    print("  PASS: K_E4 (extractor marker has provider + error)")


def test_k_e5_extractor_diagnostic_preserves_return_value():
    """K_E5: extractor diagnostic does NOT change the return value."""
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de5")

    def fake_call_llm(prompt, provider, model, log):
        raise RuntimeError("boom")
    mod._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    captured, restore = _capture_stderr()
    try:
        result = mod.extract_directives(
            "refine and stop",
            provider="google",
            log_func=log_func)
    finally:
        restore()

    assert result == {}, (
        "K_E5 failed: return value should be {} after failure, "
        "got %r (type %s)" % (result, type(result).__name__))

    print("  PASS: K_E5 (extractor diagnostic preserves return value)")


def test_k_e6_extractor_critical_fallback_via_log_func():
    """K_E6: CRITICAL_FALLBACK message goes through log_func callback."""
    mod = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de6")

    def fake_call_llm(prompt, provider, model, log):
        raise ConnectionError("network down")
    mod._call_llm = fake_call_llm

    log_func, log_msgs = _capture_log()
    _, restore = _capture_stderr()
    try:
        mod.extract_directives(
            "refine and stop",
            provider="google",
            log_func=log_func)
    finally:
        restore()

    critical_msgs = [m for m in log_msgs if "CRITICAL_FALLBACK" in m]
    assert len(critical_msgs) >= 1, (
        "K_E6 failed: CRITICAL_FALLBACK should be in log_func "
        "messages. Got: %r" % log_msgs)
    msg = critical_msgs[0]
    assert "LLM extraction failed" in msg, (
        "K_E6 failed: CRITICAL_FALLBACK should mention LLM "
        "extraction failure. Got: %r" % msg)
    assert "rule-based" in msg, (
        "K_E6 failed: CRITICAL_FALLBACK should mention fallback "
        "to rules. Got: %r" % msg)

    print("  PASS: K_E6 (extractor CRITICAL_FALLBACK via log_func)")


# --------------------------------------------------------------------
# Site 2: advice preprocessor (K_E7 - K_E10)
# --------------------------------------------------------------------

class _FakeLLM(object):
    """Minimal stand-in for an LLM client used by preprocess_advice."""
    def __init__(self, response=None, raises=None):
        self.response = response
        self.raises = raises

    def invoke(self, prompt):
        if self.raises is not None:
            raise self.raises
        class _R(object):
            def __init__(s, content):
                s.content = content
        return _R(self.response)


def test_k_e7_preprocessor_exception_emits_stderr_marker():
    """K_E7: preprocessor LLM exception → stderr contains marker."""
    mod = _import_module(ADVICE_PREPROCESSOR_PATH, "pp7")

    llm = _FakeLLM(raises=RuntimeError("simulated google credentials timeout"))

    out_buf = io.StringIO()
    captured, restore = _capture_stderr()
    try:
        result = mod.preprocess_advice(
            "refine the model and stop, generate r-free flags too please",
            experiment_type="xray",
            file_list=["nsf-d2.mtz", "nsf-d2_noligand.pdb"],
            llm=llm,
            out=out_buf)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "[ADVICE_PREPROCESSING_FAILED]" in err_output, (
        "K_E7 failed: stderr should contain "
        "[ADVICE_PREPROCESSING_FAILED] marker. Got: %r" % err_output)
    # Return value should be raw_advice unchanged
    assert "refine the model and stop" in result, (
        "K_E7 failed: should return raw_advice on failure. "
        "Got: %r" % result)

    print("  PASS: K_E7 (preprocessor exception → stderr marker)")


def test_k_e8_preprocessor_success_no_stderr_marker():
    """K_E8: preprocessor LLM success → stderr does NOT contain marker."""
    mod = _import_module(ADVICE_PREPROCESSOR_PATH, "pp8")

    # Return a valid-looking processed advice (length > 20)
    processed = (
        "1. **Input Files Found**: nsf-d2.mtz, nsf-d2_noligand.pdb\n"
        "2. **Experiment Type**: X-ray crystallography\n"
        "3. **Primary Goal**: Refine the model")
    llm = _FakeLLM(response=processed)

    out_buf = io.StringIO()
    captured, restore = _capture_stderr()
    try:
        mod.preprocess_advice(
            "refine the model and stop",
            experiment_type="xray",
            file_list=["nsf-d2.mtz"],
            llm=llm,
            out=out_buf)
    finally:
        restore()

    err_output = captured.getvalue()
    assert "[ADVICE_PREPROCESSING_FAILED]" not in err_output, (
        "K_E8 failed: stderr should NOT contain marker on success. "
        "Got: %r" % err_output)

    print("  PASS: K_E8 (preprocessor success → no stderr marker)")


def test_k_e9_preprocessor_critical_fallback_via_out():
    """K_E9: CRITICAL_FALLBACK message printed to out stream on failure."""
    mod = _import_module(ADVICE_PREPROCESSOR_PATH, "pp9")

    llm = _FakeLLM(raises=ConnectionError("network timeout"))

    out_buf = io.StringIO()
    _, restore = _capture_stderr()
    try:
        mod.preprocess_advice(
            "refine the model and stop",
            experiment_type="xray",
            file_list=["nsf-d2.mtz"],
            llm=llm,
            out=out_buf)
    finally:
        restore()

    out_text = out_buf.getvalue()
    assert "CRITICAL_FALLBACK" in out_text, (
        "K_E9 failed: CRITICAL_FALLBACK should be in out stream. "
        "Got: %r" % out_text)
    assert "LLM preprocessing failed" in out_text, (
        "K_E9 failed: CRITICAL_FALLBACK should mention LLM "
        "preprocessing failure. Got: %r" % out_text)
    assert "raw advice" in out_text, (
        "K_E9 failed: CRITICAL_FALLBACK should mention fallback to "
        "raw advice. Got: %r" % out_text)

    print("  PASS: K_E9 (preprocessor CRITICAL_FALLBACK via out)")


def test_k_e10_preprocessor_diagnostic_preserves_return_value():
    """K_E10: preprocessor diagnostic doesn't change the return value."""
    mod = _import_module(ADVICE_PREPROCESSOR_PATH, "pp10")

    llm = _FakeLLM(raises=ValueError("auth error"))
    raw_advice = "refine the model and stop, please use r-free flags"

    out_buf = io.StringIO()
    _, restore = _capture_stderr()
    try:
        result = mod.preprocess_advice(
            raw_advice,
            experiment_type="xray",
            file_list=["nsf-d2.mtz"],
            llm=llm,
            out=out_buf)
    finally:
        restore()

    # raw_advice should be returned unchanged on failure
    assert result == raw_advice, (
        "K_E10 failed: should return raw_advice unchanged on "
        "failure. Got: %r, expected: %r" % (result, raw_advice))

    print("  PASS: K_E10 (preprocessor diagnostic preserves return)")


# --------------------------------------------------------------------
# Shared structural test (K_E11)
# --------------------------------------------------------------------

class _BrokenStream(object):
    """A stream that raises on every write — for stderr robustness test."""
    def write(self, data):
        raise IOError("stderr unavailable")
    def flush(self):
        raise IOError("stderr unavailable")


def test_k_e11_both_sites_robust_to_broken_stderr():
    """K_E11: both sites complete normally if sys.stderr is broken.

    The diagnostic must NEVER break the path it's diagnosing.  If
    sys.stderr can't be written (closed, replaced with broken stream,
    etc.), the extraction/preprocessing fallback paths must still
    return their correct values.
    """
    # Site 1: extractor
    de = _import_module(DIRECTIVE_EXTRACTOR_PATH, "de11")
    def fake_call_llm(prompt, provider, model, log):
        raise RuntimeError("boom")
    de._call_llm = fake_call_llm

    log_func, _ = _capture_log()
    saved_stderr = sys.stderr
    sys.stderr = _BrokenStream()
    try:
        result1 = de.extract_directives(
            "refine and stop",
            provider="google",
            log_func=log_func)
    finally:
        sys.stderr = saved_stderr

    assert result1 == {}, (
        "K_E11 failed (extractor): return value should still be {} "
        "when stderr is broken, got %r" % result1)

    # Site 2: preprocessor
    pp = _import_module(ADVICE_PREPROCESSOR_PATH, "pp11")
    llm = _FakeLLM(raises=RuntimeError("boom"))
    raw = "refine and stop"
    out_buf = io.StringIO()

    sys.stderr = _BrokenStream()
    try:
        result2 = pp.preprocess_advice(
            raw,
            experiment_type="xray",
            file_list=["a.mtz"],
            llm=llm,
            out=out_buf)
    finally:
        sys.stderr = saved_stderr

    assert result2 == raw, (
        "K_E11 failed (preprocessor): return value should still be "
        "raw_advice when stderr is broken, got %r" % result2)

    print("  PASS: K_E11 (both sites robust to broken stderr)")


# --------------------------------------------------------------------
# Test runner
# --------------------------------------------------------------------

def run_all_tests():
    tests = [
        # Site 1 — directive extractor
        ("K_E1_extractor_exception_emits_stderr_marker",
         test_k_e1_extractor_exception_emits_stderr_marker),
        ("K_E2_extractor_success_no_stderr_marker",
         test_k_e2_extractor_success_no_stderr_marker),
        ("K_E3_extractor_empty_response_no_stderr_marker",
         test_k_e3_extractor_empty_response_no_stderr_marker),
        ("K_E4_extractor_marker_includes_provider_and_error",
         test_k_e4_extractor_marker_includes_provider_and_error),
        ("K_E5_extractor_diagnostic_preserves_return_value",
         test_k_e5_extractor_diagnostic_preserves_return_value),
        ("K_E6_extractor_critical_fallback_via_log_func",
         test_k_e6_extractor_critical_fallback_via_log_func),
        # Site 2 — preprocessor
        ("K_E7_preprocessor_exception_emits_stderr_marker",
         test_k_e7_preprocessor_exception_emits_stderr_marker),
        ("K_E8_preprocessor_success_no_stderr_marker",
         test_k_e8_preprocessor_success_no_stderr_marker),
        ("K_E9_preprocessor_critical_fallback_via_out",
         test_k_e9_preprocessor_critical_fallback_via_out),
        ("K_E10_preprocessor_diagnostic_preserves_return_value",
         test_k_e10_preprocessor_diagnostic_preserves_return_value),
        # Shared
        ("K_E11_both_sites_robust_to_broken_stderr",
         test_k_e11_both_sites_robust_to_broken_stderr),
    ]
    passed = 0
    failed = 0
    for name, fn in tests:
        print("Test: %s" % name)
        try:
            fn()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            import traceback
            print("  ERROR: %s" % e)
            print(traceback.format_exc())
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
