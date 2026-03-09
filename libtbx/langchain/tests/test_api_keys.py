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

def test_openai():
  """Test OpenAI API key."""
  key = os.environ.get("OPENAI_API_KEY", "")
  if not key:
    return "NOT SET", "OPENAI_API_KEY environment variable is not set"
  masked = key[:8] + "..." + key[-4:]

  try:
    # Use urllib so there are no dependencies beyond stdlib
    if sys.version_info.major >= 3:
      from urllib.request import Request, urlopen
      from urllib.error import HTTPError, URLError
    else:
      from urllib2 import Request, urlopen, HTTPError, URLError

    body = json.dumps({
      "model": "gpt-4o-mini",
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
    model = data.get("model", "?")
    return "OK", "Key %s  model=%s  reply=%s" % (masked, model, reply)

  except HTTPError as e:
    body = ""
    try:
      body = e.read().decode("utf-8", errors="replace")[:200]
    except Exception:
      pass
    if e.code == 401:
      return "INVALID KEY", "Key %s  HTTP 401: %s" % (masked, body)
    elif e.code == 429:
      return "RATE LIMITED", "Key %s  HTTP 429: %s" % (masked, body)
    else:
      return "ERROR", "Key %s  HTTP %d: %s" % (masked, e.code, body)
  except URLError as e:
    return "NETWORK ERROR", "Key %s  %s" % (masked, str(e.reason))
  except Exception as e:
    return "ERROR", "Key %s  %s" % (masked, str(e))


def test_google():
  """Test Google (Gemini) API key."""
  key = os.environ.get("GOOGLE_API_KEY", "")
  if not key:
    return "NOT SET", "GOOGLE_API_KEY environment variable is not set"
  masked = key[:8] + "..." + key[-4:]

  try:
    if sys.version_info.major >= 3:
      from urllib.request import Request, urlopen
      from urllib.error import HTTPError, URLError
    else:
      from urllib2 import Request, urlopen, HTTPError, URLError

    model = "gemini-2.5-flash"
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
    return "OK", "Key %s  model=%s  reply=%s" % (masked, model, reply)

  except HTTPError as e:
    body = ""
    try:
      body = e.read().decode("utf-8", errors="replace")[:200]
    except Exception:
      pass
    if e.code == 400:
      return "INVALID KEY", "Key %s  HTTP 400: %s" % (masked, body)
    elif e.code == 403:
      return "FORBIDDEN", "Key %s  HTTP 403: %s" % (masked, body)
    elif e.code == 429:
      return "RATE LIMITED", "Key %s  HTTP 429: %s" % (masked, body)
    else:
      return "ERROR", "Key %s  HTTP %d: %s" % (masked, e.code, body)
  except URLError as e:
    return "NETWORK ERROR", "Key %s  %s" % (masked, str(e.reason))
  except Exception as e:
    return "ERROR", "Key %s  %s" % (masked, str(e))


def main():
  print()
  print("=" * 60)
  print("  PHENIX AI Agent — API Key Test")
  print("=" * 60)

  print()
  print("OpenAI (GPT-4o-mini):")
  status, detail = test_openai()
  print("  Status: %s" % status)
  print("  %s" % detail)

  print()
  print("Google (Gemini 2.5 Flash):")
  status_g, detail_g = test_google()
  print("  Status: %s" % status_g)
  print("  %s" % detail_g)

  print()
  print("-" * 60)
  ok = 0
  for s in [status, status_g]:
    if s == "OK":
      ok += 1
    elif s == "NOT SET":
      pass
    else:
      print("  WARNING: at least one key failed")
  if ok == 0:
    print("  No working API keys found.")
    print("  Set OPENAI_API_KEY or GOOGLE_API_KEY and retry.")
  elif ok == 1:
    print("  One API is working. The agent needs at least one.")
  else:
    print("  Both APIs are working.")
  print()
  return 0 if ok > 0 else 1


if __name__ == "__main__":
  sys.exit(main())
