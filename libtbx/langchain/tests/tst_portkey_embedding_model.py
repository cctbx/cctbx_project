"""Regression tests for the portkey embedding model name + silent-fallback guard.

Two issues from the field:
1. The portkey/Azure embedding DEPLOYMENT name is "text-embedding-3-small-1"
   (trailing "-1"), NOT the bare OpenAI id "text-embedding-3-small".  With the
   bare name the Azure gateway silently fell back to a chat deployment
   (gpt-5-mini) instead of erroring, corrupting the database.
2. verify_embeddings() preflights the embeddings with one real call and prints
   a LOUD warning (does NOT raise / abort) when the result looks wrong --
   notably a dimension mismatch for the fixed-size OpenAI 1536 family.

Source-extraction style (no PHENIX import): the EMBEDDING_EXPECTED_DIM table
and verify_embeddings body are pulled from core/llm.py and exercised directly,
so the tests track the real source.  2-space indent.
"""

from __future__ import absolute_import, division, print_function

import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)                      # .../libtbx/langchain
_LLM_CANDIDATES = [
    os.path.join(_ROOT, "core", "llm.py"),
]


def _find_llm():
    for c in _LLM_CANDIDATES:
        c = os.path.abspath(c)
        if os.path.isfile(c):
            return c
    e = os.environ.get("LLM_PY")
    if e and os.path.isfile(e):
        return e
    return None


def _load_verify():
    """Extract EMBEDDING_EXPECTED_DIM + verify_embeddings from real source."""
    path = _find_llm()
    if path is None:
        return None, None
    src = open(path).read()
    dim = re.search(r'EMBEDDING_EXPECTED_DIM = \{.*?\}', src, re.DOTALL)
    fn = re.search(r'def verify_embeddings\(.*?\n  return \(True, observed\)',
                   src, re.DOTALL)
    if not (dim and fn):
        return None, None
    ns = {}
    exec(dim.group(0) + "\n" + fn.group(0), ns)
    return ns.get("verify_embeddings"), ns.get("EMBEDDING_EXPECTED_DIM")


# --- Issue 1: the model name -------------------------------------------------

def test_portkey_embedding_name_has_minus_one():
    path = _find_llm()
    if path is None:
        print("  (skip) core/llm.py not found")
        return
    src = open(path).read()
    block = re.search(r'RAG_EMBEDDING_DEFAULTS = \{.*?\}', src, re.DOTALL)
    assert block, "RAG_EMBEDDING_DEFAULTS not found"
    body = block.group(0)
    # portkey must use the Azure deployment alias WITH the trailing -1
    m = re.search(r'"portkey":\s*"([^"]+)"', body)
    assert m, "portkey embedding entry not found"
    assert m.group(1) == "text-embedding-3-small-1", \
        "portkey embedding must be 'text-embedding-3-small-1' (got %r)" % m.group(1)


def test_direct_openai_embedding_name_unchanged():
    """The bare OpenAI provider must keep the bare id (no -1)."""
    path = _find_llm()
    if path is None:
        print("  (skip)")
        return
    src = open(path).read()
    body = re.search(r'RAG_EMBEDDING_DEFAULTS = \{.*?\}', src, re.DOTALL).group(0)
    m = re.search(r'"openai":\s*"([^"]+)"', body)
    assert m and m.group(1) == "text-embedding-3-small", \
        "direct openai embedding must remain 'text-embedding-3-small'"


# --- Issue 2: the verify_embeddings guard ------------------------------------

class _Good:
    def embed_query(self, t):
        return [0.1] * 1536


class _WrongDim:
    def embed_query(self, t):
        return [0.0] * 8         # gpt-5-mini-style fallback


class _Raises:
    def embed_query(self, t):
        raise Exception("404 deployment not found")


class _Empty:
    def embed_query(self, t):
        return []


class _Gemini3072:
    def embed_query(self, t):
        return [0.1] * 3072      # configurable model, dim NOT asserted


def test_verify_ok_for_correct_1536():
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip) could not load verify_embeddings")
        return
    msgs = []
    ok, dim = verify(_Good(), embedding_model_name="text-embedding-3-small-1",
                     provider="portkey", log=msgs.append)
    assert ok is True and dim == 1536
    assert not any("WARNING" in m for m in msgs), "no warning on a correct model"


def test_verify_warns_on_dimension_mismatch():
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    msgs = []
    ok, dim = verify(_WrongDim(), embedding_model_name="text-embedding-3-small-1",
                     provider="portkey", log=msgs.append)
    assert ok is False, "must report not-ok on dimension mismatch"
    assert any("DIFFERENT" in m or "gpt-5-mini" in m for m in msgs), \
        "warning must explain the silent fallback"


def test_verify_does_not_raise_on_failure():
    """Warn-don't-fail: an embed call that raises must be caught, not propagated."""
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    msgs = []
    try:
        ok, dim = verify(_Raises(),
                         embedding_model_name="text-embedding-3-small-1",
                         provider="portkey", log=msgs.append)
    except Exception as e:
        raise AssertionError("verify_embeddings must NOT raise; got %r" % e)
    assert ok is False and dim is None
    assert any("trailing -1" in m or "deployment name" in m for m in msgs)


def test_verify_warns_on_empty_vector():
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    msgs = []
    ok, dim = verify(_Empty(), embedding_model_name="text-embedding-3-small-1",
                     provider="portkey", log=msgs.append)
    assert ok is False and dim is None
    assert any("no usable vector" in m for m in msgs)


def test_verify_no_false_warning_for_variable_dim_model():
    """gemini-embedding-001 is configurable (768/1536/3072); its dimension is
    NOT in EMBEDDING_EXPECTED_DIM, so verify must not warn on any size."""
    verify, table = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    assert "gemini-embedding-001" not in table, \
        "variable-dimension models must be absent from the expected-dim table"
    msgs = []
    ok, dim = verify(_Gemini3072(), embedding_model_name="gemini-embedding-001",
                     provider="google", log=msgs.append)
    assert ok is True and dim == 3072
    assert not any("WARNING" in m for m in msgs), \
        "must not false-warn on a configurable-dimension model"


class _StrReturn:
    def embed_query(self, t):
        return "error: deployment not found"   # pathological non-vector


def test_verify_accepts_numpy_vector():
    """A numpy float32 vector of the right size must be ACCEPTED (no false
    warning) -- some embedding backends return numpy arrays."""
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    try:
        import numpy as np
    except ImportError:
        print("  (skip) numpy not available")
        return

    class _NumpyGood:
        def embed_query(self, t):
            return np.zeros(1536, dtype=np.float32)

    msgs = []
    ok, dim = verify(_NumpyGood(), embedding_model_name="text-embedding-3-small-1",
                     provider="portkey", log=msgs.append)
    assert ok is True and dim == 1536, (ok, dim)
    assert not any("WARNING" in m for m in msgs), \
        "a valid numpy vector must not trigger a false warning"


def test_verify_warns_on_string_return():
    """A string return (not a vector) must be flagged, not mistaken for one."""
    verify, _ = _load_verify()
    if verify is None:
        print("  (skip)")
        return
    msgs = []
    ok, dim = verify(_StrReturn(), embedding_model_name="text-embedding-3-small-1",
                     provider="portkey", log=msgs.append)
    assert ok is False and dim is None
    assert any("no usable vector" in m for m in msgs)


def test_run_build_db_calls_verify_before_persist():
    """Source-scan: the build path must preflight embeddings AFTER constructing
    them and BEFORE persisting, so every caller (rebuild/update/...) is guarded
    at one canonical place."""
    cands = [
        os.path.join(_ROOT, "run_build_db.py"),
    ]
    path = None
    for c in cands:
        c = os.path.abspath(c)
        if os.path.isfile(c):
            path = c
            break
    if path is None and os.environ.get("RUN_BUILD_DB_PY"):
        path = os.environ["RUN_BUILD_DB_PY"]
    if not path or not os.path.isfile(path):
        print("  (skip) run_build_db.py not found")
        return
    src = open(path).read()
    assert "verify_embeddings" in src, \
        "run_build_db must call verify_embeddings"
    i_construct = src.find("get_llm_and_embeddings(")
    i_verify = src.find("verify_embeddings(")
    i_persist = src.find("create_and_persist_db(")
    assert -1 not in (i_construct, i_verify, i_persist), \
        "expected construct + verify + persist calls all present"
    assert i_construct < i_verify < i_persist, \
        "verify_embeddings must run AFTER construction and BEFORE persisting"


_TESTS = [
    test_portkey_embedding_name_has_minus_one,
    test_direct_openai_embedding_name_unchanged,
    test_verify_ok_for_correct_1536,
    test_verify_warns_on_dimension_mismatch,
    test_verify_does_not_raise_on_failure,
    test_verify_warns_on_empty_vector,
    test_verify_no_false_warning_for_variable_dim_model,
    test_verify_accepts_numpy_vector,
    test_verify_warns_on_string_return,
    test_run_build_db_calls_verify_before_persist,
]


def run_all_tests():
    for fn in _TESTS:
        fn()
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    p = f = 0
    for fn in _TESTS:
        try:
            fn(); print("  PASS: %s" % fn.__name__); p += 1
        except AssertionError as e:
            print("  FAIL: %s -- %s" % (fn.__name__, e)); f += 1
    print("\n%d passed, %d failed" % (p, f))
