#!/usr/bin/env python3
"""
Test whether OpenAI and Google API keys are configured and working.

Usage:
  python test_api_keys.py
  phenix.python test_api_keys.py
"""

from __future__ import print_function
import os
import sys
import json


# v119.H1: pull model names from the central tables in core/llm.py.
# This script is deliberately stdlib-only at the standard library
# level -- core/llm.py is also stdlib-only at module import (the
# langchain_* dependencies are deferred inside get_llm_and_embeddings),
# so importing the helper here doesn't break standalone usability.
# When invoked directly from tests/, sys.path[0] is the tests/
# directory, not langchain/.  Prepend langchain/ so 'core.llm'
# resolves whether the user runs 'python tests/test_api_keys.py'
# from langchain/ or anywhere else.
_LANGCHAIN_DIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
if _LANGCHAIN_DIR not in sys.path:
    sys.path.insert(0, _LANGCHAIN_DIR)

# Standard agent import-fallback pattern: try the libtbx-installed
# path first, then the direct path.
try:
  from libtbx.langchain.core.llm import default_model_for_provider
except ImportError:
  try:
    from core.llm import default_model_for_provider
  except ImportError:
    # Hard failure: if core/llm.py is unreachable, this diagnostic
    # script can't reliably check the right models.  Print a clear
    # error and exit non-zero so the user fixes their invocation.
    print(
      "test_api_keys: cannot import core.llm.default_model_for_provider.\n"
      "  Run via 'phenix.python tests/test_api_keys.py' or from the\n"
      "  langchain/ directory with 'python tests/test_api_keys.py'.",
      file=sys.stderr)
    sys.exit(2)

def test_openai(model=None):
  """Test OpenAI API key."""
  # v119.H1: resolve model from central defaults BEFORE the key
  # check so the diagnostic label is meaningful even with no key.
  if model is None:
    model = default_model_for_provider("openai", role="decision")

  key = os.environ.get("OPENAI_API_KEY", "")
  if not key:
    return "NOT SET", "OPENAI_API_KEY environment variable is not set", model
  masked = key[:8] + "..." + key[-4:]

  try:
    # Use urllib so there are no dependencies beyond stdlib
    if sys.version_info.major >= 3:
      from urllib.request import Request, urlopen
      from urllib.error import HTTPError, URLError
    else:
      from urllib2 import Request, urlopen, HTTPError, URLError

    body = json.dumps({
      "model": model,
      "messages": [{"role": "user", "content": "Reply with exactly: OK"}],
      "max_tokens": 5,
    }).encode("utf-8")

    req = Request(
      "https://api.openai.com/v1/chat/completions",
      data=body,
      headers={
        "Authorization": "Bearer %s" % key,
        "Content-Type": "application/json",
      },
    )
    resp = urlopen(req, timeout=30)
    data = json.loads(resp.read().decode("utf-8"))
    reply = data["choices"][0]["message"]["content"].strip()
    actual_model = data.get("model", "?")
    return "OK", "Key %s  model=%s  reply=%s" % (masked, actual_model, reply), model

  except HTTPError as e:
    body = ""
    try:
      body = e.read().decode("utf-8", errors="replace")[:200]
    except Exception:
      pass
    if e.code == 401:
      return "INVALID KEY", "Key %s  HTTP 401: %s" % (masked, body), model
    elif e.code == 429:
      return "RATE LIMITED", "Key %s  HTTP 429: %s" % (masked, body), model
    else:
      return "ERROR", "Key %s  HTTP %d: %s" % (masked, e.code, body), model
  except URLError as e:
    return "NETWORK ERROR", "Key %s  %s" % (masked, str(e.reason)), model
  except Exception as e:
    return "ERROR", "Key %s  %s" % (masked, str(e)), model


def test_google(model=None):
  """Test Google (Gemini) API key."""
  # v119.H1: resolve model from central defaults BEFORE the key
  # check so the diagnostic label is meaningful even with no key.
  # Tests the EXPENSIVE model (used by ai_analysis/ai_agent) -- the
  # cheaper decision model has a separate larger quota and would
  # show OK even when the expensive quota is exhausted.
  if model is None:
    model = default_model_for_provider("google", role="expensive")

  key = os.environ.get("GOOGLE_API_KEY", "")
  if not key:
    return "NOT SET", "GOOGLE_API_KEY environment variable is not set", model
  masked = key[:8] + "..." + key[-4:]

  try:
    if sys.version_info.major >= 3:
      from urllib.request import Request, urlopen
      from urllib.error import HTTPError, URLError
    else:
      from urllib2 import Request, urlopen, HTTPError, URLError

    url = (
      "https://generativelanguage.googleapis.com/v1beta"
      "/models/%s:generateContent?key=%s" % (model, key)
    )
    body = json.dumps({
      "contents": [
        {"parts": [{"text": "Reply with exactly: OK"}]}
      ],
      "generationConfig": {"maxOutputTokens": 5},
    }).encode("utf-8")

    req = Request(
      url, data=body,
      headers={"Content-Type": "application/json"},
    )
    resp = urlopen(req, timeout=30)
    data = json.loads(resp.read().decode("utf-8"))
    reply = (data.get("candidates", [{}])[0]
             .get("content", {})
             .get("parts", [{}])[0]
             .get("text", "").strip())
    return "OK", "Key %s  model=%s  reply=%s" % (masked, model, reply), model

  except HTTPError as e:
    body = ""
    try:
      body = e.read().decode("utf-8", errors="replace")[:200]
    except Exception:
      pass
    if e.code == 400:
      return "INVALID KEY", "Key %s  HTTP 400: %s" % (masked, body), model
    elif e.code == 403:
      return "FORBIDDEN", "Key %s  HTTP 403: %s" % (masked, body), model
    elif e.code == 429:
      return "QUOTA EXCEEDED", "Key %s  HTTP 429 (%s quota exhausted): %s" % (masked, model, body), model
    else:
      return "ERROR", "Key %s  HTTP %d: %s" % (masked, e.code, body), model
  except URLError as e:
    return "NETWORK ERROR", "Key %s  %s" % (masked, str(e.reason)), model
  except Exception as e:
    return "ERROR", "Key %s  %s" % (masked, str(e)), model


def main():
  print()
  print("=" * 60)
  print("  PHENIX AI Agent - API Key Test")
  print("=" * 60)

  print()
  status, detail, openai_model = test_openai()
  print("OpenAI (%s):" % openai_model)
  print("  Status: %s" % status)
  print("  %s" % detail)

  print()
  status_g, detail_g, google_model = test_google()
  print("Google (%s):" % google_model)
  print("  Status: %s" % status_g)
  print("  %s" % detail_g)

  print()
  print("-" * 60)
  ok = 0
  quota_exceeded = []
  for name, s in [("OpenAI", status), ("Google", status_g)]:
    if s == "OK":
      ok += 1
    elif s == "QUOTA EXCEEDED":
      quota_exceeded.append(name)
    elif s == "NOT SET":
      pass
    else:
      print("  WARNING: %s key check failed (%s)" % (name, s))
  if quota_exceeded:
    # google_model is the one whose quota is most likely to bite,
    # since the openai test uses the cheaper decision model.
    print("  NOTE: %s quota exhausted for %s." % (
      "/".join(quota_exceeded), google_model))
    print("  The agent will fail with this provider until quota resets.")
    if ok > 0:
      print("  Use the other working provider instead.")
  if ok == 0 and not quota_exceeded:
    print("  No working API keys found.")
    print("  Set OPENAI_API_KEY or GOOGLE_API_KEY and retry.")
  elif ok == 1:
    print("  One API is working. The agent needs at least one.")
  elif ok == 2:
    print("  Both APIs are working.")
  print()
  return 0 if ok > 0 else 1


if __name__ == "__main__":
  sys.exit(main())
