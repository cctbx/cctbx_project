"""
v119.H13 — Provider-error classification tests.

Pins behavior of core/llm._classify_provider_error: the three-way
404 sub-classification (RETIRED / UNAVAILABLE / FAILED) plus
AUTH_FAILED for 401-class errors.

Background: H13 was triggered by AIAgent_run_39a_ollama where an
ollama 404 "model not found" was opaque to operators.  The
classifier distinguishes:

  - MODEL_RETIRED:     provider has deprecated the model.
                       Action: update core/llm DEFAULT_MODELS.
  - MODEL_UNAVAILABLE: model unknown to this server (Ollama not
                       pulled, typo in OLLAMA_LLM_MODEL).  Action:
                       pull/rename.
  - AUTH_FAILED:       401-class.  Action: check API key.
  - FAILED:            generic / unrecognized.  Action: investigate.

Each test follows the sandbox-skip pattern used elsewhere in the
codebase (nested try/except for libtbx.langchain vs bare core).

Run with:
    PYTHONPATH=. python tests/tst_provider_error_classification.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal,
    assert_in,
    assert_true,
    run_tests_with_fail_fast,
)

# PHENIX/cctbx linter silencer
(assert_equal, assert_in, assert_true, run_tests_with_fail_fast)


def _import_classifier():
    """Sandbox-skip import of the classifier."""
    try:
        from libtbx.langchain.core.llm import _classify_provider_error
        return _classify_provider_error
    except ImportError:
        try:
            from core.llm import _classify_provider_error
            return _classify_provider_error
        except ImportError:
            return None


class _MockExc(Exception):
    """Plain exception class for synthetic test cases."""
    pass


class _MockResponseExc(Exception):
    """Exception with .status_code attribute (mocks library exceptions
    that surface HTTP status programmatically rather than only in __str__).
    """
    def __init__(self, message, status_code=None, body=None):
        super().__init__(message)
        if status_code is not None:
            self.status_code = status_code
        if body is not None:
            self.body = body


class _MockChainedResponseExc(Exception):
    """Exception with .response.status_code chain (mocks SDK clients
    that wrap an HTTP response object on the exception).
    """
    def __init__(self, message, status_code):
        super().__init__(message)

        class _Resp:
            pass
        self.response = _Resp()
        self.response.status_code = status_code


# =========================================================================
# RETIRED cases — provider has dropped the model
# =========================================================================

def test_classify_google_retired_model():
    """Google's verbose 404 'This model version is no longer available'
    classifies as MODEL_RETIRED.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        "404 Client Error: Not Found for url: ... "
        "models/gemini-1.0/...\nThis model version is no longer available.")
    tag, hint = fn(e, "gemini-1.0")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_RETIRED",
                 "Google retirement → MODEL_RETIRED")
    assert_true(hint and "gemini-1.0" in hint,
                "Hint mentions the model name")
    assert_in("DEFAULT_MODELS", hint,
              "Hint points operator at DEFAULT_MODELS")


def test_classify_openai_deprecated_model():
    """OpenAI's clean 404 'The model X has been deprecated' classifies
    as MODEL_RETIRED — 'deprecated' phrase co-occurs with 'model'.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        "Error code: 404 - {'error': {'message': "
        "'The model gpt-3.5-turbo-0301 has been deprecated.', "
        "'type': 'invalid_request_error'}}")
    tag, _ = fn(e, "gpt-3.5-turbo-0301")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_RETIRED",
                 "OpenAI deprecated → MODEL_RETIRED")


def test_classify_anthropic_retired_model():
    """Anthropic-style 'model has been retired' classifies as RETIRED."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        "model_not_found: The requested model 'claude-1' has been retired")
    tag, _ = fn(e, "claude-1")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_RETIRED",
                 "Anthropic retirement → MODEL_RETIRED")


# =========================================================================
# UNAVAILABLE cases — model name unknown to this server
# =========================================================================

def test_classify_ollama_model_not_pulled():
    """The Tom case: Ollama returns 404 with 'model X not found' but
    no retirement language.  Classifies as MODEL_UNAVAILABLE.

    THIS IS THE PRIMARY REGRESSION TEST for AIAgent_run_39a_ollama.
    Pre-H13 this was opaque ('DIRECTIVE_EXTRACTION_FAILED: 404 page
    not found').  Post-H13 the operator sees an actionable
    [DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE] marker with a 'pull X'
    hint.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        'Error code: 404 - {"error":{"message":"model \'llama3.2\' '
        'not found","type":"api_error","param":null,"code":null}}')
    tag, hint = fn(e, "llama3.2")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE",
                 "Ollama 'model not found' → MODEL_UNAVAILABLE")
    assert_true(hint and "llama3.2" in hint,
                "Hint mentions the model name")
    assert_in("ollama pull", hint.lower(),
              "Hint suggests 'ollama pull' as the operator action")


def test_classify_ollama_unavailable_with_qwen():
    """Variant of the Tom case with a different model name to ensure
    the regex isn't model-name-specific.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        '404 - {"error":{"message":"model \'qwen2.5:72b\' not found"}}')
    tag, hint = fn(e, "qwen2.5:72b")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE",
                 "Different model name same classification")
    assert_true("qwen2.5:72b" in hint,
                "Hint mentions the colon-containing model name")


def test_classify_unavailable_double_quoted_format():
    """Double-quoted format like 'model "X" not found' also matches."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc('404 - error: model "gemma:7b" not found on server')
    tag, _ = fn(e, "gemma:7b")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE",
                 "Double-quoted variant → MODEL_UNAVAILABLE")


# =========================================================================
# FAILED cases — load-bearing false-positive guards
# =========================================================================

def test_classify_edge_proxy_404_not_retired():
    """Edge-proxy 404 page mentioning 'endpoint is deprecated' must NOT
    be classified as MODEL_RETIRED.  The co-occurrence rule requires
    'model' to be near 'deprecated'; here 'endpoint' is near
    'deprecated' but 'model' is not.

    LOAD-BEARING: if this ever flips to MODEL_RETIRED, the
    co-occurrence rule has been weakened in a way that pollutes
    operator alerts with noise.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc(
        "404 Not Found: The requested resource at /api/v1/... has moved. "
        "This endpoint is deprecated; see new path.")
    tag, _ = fn(e, "some-model")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_FAILED",
                 "Edge-proxy 404 with 'endpoint deprecated' → FAILED "
                 "(NOT RETIRED — load-bearing co-occurrence sentinel)")


def test_classify_generic_404_unrecognized():
    """Plain 404 with no retirement language and no 'model not found'
    pattern falls through to FAILED.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc("404 Not Found")
    tag, _ = fn(e, "foo")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_FAILED",
                 "Bare 404 → FAILED")


def test_classify_empty_exception():
    """Empty / whitespace-only exception string handles gracefully."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc("")
    tag, _ = fn(e, "foo")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_FAILED",
                 "Empty exception → FAILED, no crash")


# =========================================================================
# AUTH cases
# =========================================================================

def test_classify_401_unauthorized():
    """401 status classifies as AUTH_FAILED."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc("401 Unauthorized")
    tag, hint = fn(e, "foo")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_AUTH_FAILED",
                 "401 → AUTH_FAILED")
    assert_true(hint and "API key" in hint,
                "Hint mentions API key")


def test_classify_invalid_api_key():
    """'invalid api key' phrase (no status code) classifies as AUTH_FAILED."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    e = _MockExc("Authentication error: invalid API key provided")
    tag, _ = fn(e, "foo")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_AUTH_FAILED",
                 "'invalid api key' → AUTH_FAILED")


# =========================================================================
# Class-attribute defense — exception with structured properties
# =========================================================================

def test_classify_via_status_code_attribute():
    """SDK exception with status_code=404 attribute classifies correctly
    even if __str__ doesn't contain '404'.

    Defense-in-depth (Gemini's H13 review): some library exceptions
    stringify to an internal error code, hiding the body.  The
    classifier interrogates .status_code programmatically.
    """
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    # __str__ is just the message — no '404' visible textually.
    # But .status_code=404 and the body contains retirement language.
    e = _MockResponseExc(
        message="<api_error>",
        status_code=404,
        body={"error": {"message": "The model has been deprecated."}})
    tag, _ = fn(e, "foo-model")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_RETIRED",
                 "Status-code attribute interrogation works "
                 "(stringification doesn't show 404)")


def test_classify_via_chained_response_status():
    """SDK exception with .response.status_code (chained) also works."""
    fn = _import_classifier()
    if fn is None:
        print("  SKIP (core.llm not importable)")
        return

    # __str__ doesn't show 404; chained attribute does.
    # Body would normally be in .response but mock keeps it simple:
    # message contains 'model not found' so the second pattern fires.
    e = _MockChainedResponseExc(
        message="api call: model 'foo-7b' not found",
        status_code=404)
    tag, _ = fn(e, "foo-7b")
    assert_equal(tag, "DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE",
                 "Chained response.status_code interrogation works")


def run_all_tests():
    """Run all provider-error classification tests."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    success = run_tests_with_fail_fast()
    if not success:
        sys.exit(1)
    print("\nOK")
