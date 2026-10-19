"""
LLM and Embeddings setup for the Phenix Crystallography Agent.

This module handles initialization of Language Models and Embedding
models from different providers (Ollama, Google, OpenAI).

Usage:
  from libtbx.langchain.core.llm import get_llm_and_embeddings

  # Use local Ollama (recommended)
  llm, embeddings = get_llm_and_embeddings(provider='ollama')

  # Or with Google
  llm, embeddings = get_llm_and_embeddings(provider='google')
"""
from __future__ import absolute_import, division, print_function

import hashlib
import json
import os


# Global verbosity setting for LLM module
# Can be set via set_llm_verbosity() or LLM_VERBOSITY environment variable
_llm_verbosity = None


# =====================================================================
# Default model strings -- single source of truth for the entire agent.
# =====================================================================
# Five role-labelled tables.  When a provider retires a model:
#   1. Update the relevant entry below.
#   2. Add the retired model name to RETIRED_MODELS.
#   3. That's the entire change -- every caller in agent/ and core/
#      reads from these tables via default_model_for_provider().
#
# Why five tables instead of one keyed by (role, provider)?  Flat
# constants are grep-friendlier ("grep DECISION_MODEL_DEFAULTS"
# returns every decision default in one shot) and let readers scan
# each role's lineup independently.
#
# Note on capability tier (counterintuitive): the DECISION path
# intentionally uses SMALLER/CHEAPER models than the RAG path.
# Decision-making is a narrow structured-JSON extraction over a
# short prompt and runs once per session -- speed and cost matter
# more than breadth.  The RAG path is broader semantic comprehension
# over PHENIX docs -- capability matters more.  Don't "fix" this
# inversion by promoting decision models to RAG models without
# revisiting the budget arithmetic.

DECISION_MODEL_DEFAULTS = {
  "google":    "gemini-3.5-flash-lite",   # was gemini-2.5-flash-lite (2.5: existing users only)
  "openai":    "gpt-6-luna",              # was gpt-4o-mini (snapshot retires 2027-01-14)
  "anthropic": "claude-sonnet-5-5",   # was claude-sonnet-4-6 (still active)
  "ollama":    "llama3.2",
  # portkey fronts Azure OpenAI; the upstream model is chosen by
  # PHENIX_PORTKEY_MODEL (default below).  gpt-5 is the deployment default.
  "portkey":   "gpt-5",
}

RAG_MODEL_DEFAULTS = {
  # anthropic has no entry: get_llm_and_embeddings delegates anthropic
  # embeddings (no native endpoint) and anthropic is not used for RAG.
  "google":    "gemini-3.5-flash-lite",   # was gemini-2.5-flash-lite
  "openai":    "gpt-6-luna",   # was gpt-5.4-nano (retires 2027-04-01); before that gpt-5-nano
  "ollama":    "llama3.1:70b",
  "portkey":   "gpt-5",
}

RAG_EMBEDDING_DEFAULTS = {
  "google":    "gemini-embedding-001",
  "openai":    "text-embedding-3-small",
  "ollama":    "nomic-embed-text",
  # portkey embeddings route through the same Azure gateway using an
  # OpenAI-compatible embedding deployment.  NOTE: the Azure DEPLOYMENT name
  # carries a trailing "-1" ("text-embedding-3-small-1"), which is NOT the same
  # as the bare OpenAI model id "text-embedding-3-small" used for the direct
  # openai provider above.  Using the bare name here causes the Azure gateway
  # to silently fall back to a default chat deployment (gpt-5-mini) instead of
  # erroring, which corrupts the embeddings.  Keep the "-1".
  "portkey":   "text-embedding-3-small-1",
}

# Expected embedding vector dimension, per embedding MODEL NAME, used to detect
# a silent provider/gateway fallback to the wrong model.  Only models whose
# dimension is FIXED and reliable are listed:
#   - The OpenAI text-embedding-3-small family is always 1536.  This is the
#     important one: when the portkey/Azure gateway is given an unknown
#     embedding deployment name it silently falls back to a chat deployment
#     (e.g. gpt-5-mini), which does NOT return a 1536-vector -- a dimension
#     check catches that.
# gemini-embedding-001 (default 3072 but truncatable to 768/1536 via MRL) and
# nomic-embed-text (varies by build) are deliberately OMITTED: their dimension
# is configurable/ambiguous, so a strict equality check would risk FALSE
# failures.  Absence here means "known vector, dimension not asserted".
EMBEDDING_EXPECTED_DIM = {
  "text-embedding-3-small":   1536,
  "text-embedding-3-small-1": 1536,   # portkey/Azure deployment alias
  "text-embedding-3-large":   3072,
  "text-embedding-ada-002":   1536,
}

EXPENSIVE_MODEL_DEFAULTS = {
  "google":    "gemini-3.8-flash",   # was gemini-2.5-pro (2.5: existing users only; no stable 3.x Pro)
  "openai":    "gpt-5.5",   # was gpt-5; gpt-5-2025-08-07 retires 2026-12-11
  "ollama":    "qwen3:32b",
  "anthropic": "claude-opus-5-5",   # was claude-opus-4-7 (still active)
  "portkey":   "gpt-5",
}

# CHEAP only contains explicit overrides -- remote providers fall
# through to RAG_MODEL_DEFAULTS via get_llm_and_embeddings.  Ollama
# overrides to a smaller local model.
CHEAP_MODEL_DEFAULTS = {
  "ollama":    "qwen2.5:7b",
}

# Known-retired model names across providers.  tst_default_models.py
# asserts no value in the five DEFAULTS tables above appears in this
# set.  Append as providers retire.
RETIRED_MODELS = frozenset([
  # LLM models -- Google
  "gemini-2.0-flash",         # retired pre-v118.8
  "gemini-1.5-flash",
  "gemini-1.5-pro",
  "gemini-1.0-pro",
  # LLM models -- OpenAI
  "gpt-3.5-turbo",
  "gpt-4-vision-preview",
  # OpenAI GPT-5 / o3 dated snapshots retired 2026-12-11 (notice 2026-06-11).
  # Bare aliases (gpt-5, gpt-5-nano) are intentionally NOT listed: bare "gpt-5"
  # is still the portkey/Azure deployment name in the DEFAULTS tables, and the
  # tst_default_models disjointness check forbids any DEFAULTS value here.
  "gpt-5-2025-08-07",
  "gpt-5-mini-2025-08-07",
  "gpt-5-nano-2025-08-07",
  "gpt-5-pro-2025-10-06",
  "o3-2025-04-16",
  "o3-pro-2025-06-10",
  # LLM models -- Anthropic
  "claude-2",
  "claude-instant-1",
  "claude-sonnet-4-20250514",  # Claude Sonnet 4, retired 2026-04-20
  "claude-opus-4-20250514",    # Claude Opus 4, retired 2026-04-20
  # LLM models -- Ollama (older defaults)
  "llama2",
  # Embedding models -- Google
  "text-embedding-004",       # superseded by gemini-embedding-001
  # Embedding models -- OpenAI
  "text-embedding-ada-002",
])


# =====================================================================
# Supported LLM providers -- single source of truth.
# =====================================================================
# Canonical list of provider names the agent accepts.  graph_nodes.py
# and phenix/programs/ai_agent.py both import THIS list rather than
# keeping their own copies, so a provider can never be wired into one
# layer but silently rejected by another.  Keep this in sync with the
# `provider` PHIL enums in ai_agent.py and ai_analysis.py (a test
# asserts the enum set covers this list).
#
# v120: anthropic (Claude) and portkey (Azure-OpenAI gateway) added.
SUPPORTED_PROVIDERS = ["google", "openai", "ollama", "anthropic", "portkey"]


def compute_defaults_fingerprint():
  """Return SHA-256 fingerprint of the centralized active defaults.

  v119.H2: exposed via the agent_build response field so operators
  can verify what defaults the server has loaded.  Changes
  whenever any default value in the five DEFAULTS tables changes.

  Excludes RETIRED_MODELS: that's tombstone data, not active
  capability.  Including it would force the H3 startup canary to
  update its expected value on every retirement, adding noise
  without operational signal.

  Format: 'sha256:<hex>' (RFC-style algorithm prefix).
  Determinism: sort_keys + fixed separators + ASCII-only.

  Returns:
    str: e.g. 'sha256:9f86d081884c7d659a2feaa0c55ad015a3...'
  """
  payload = {
    "decision":      dict(DECISION_MODEL_DEFAULTS),
    "rag":           dict(RAG_MODEL_DEFAULTS),
    "rag_embedding": dict(RAG_EMBEDDING_DEFAULTS),
    "expensive":     dict(EXPENSIVE_MODEL_DEFAULTS),
    "cheap":         dict(CHEAP_MODEL_DEFAULTS),
  }
  canonical = json.dumps(
    payload,
    sort_keys=True,
    separators=(",", ":"),
    ensure_ascii=True,
  )
  digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
  return "sha256:" + digest


_ROLE_TABLES = {
  "decision":       DECISION_MODEL_DEFAULTS,
  "rag":            RAG_MODEL_DEFAULTS,
  "rag_embedding":  RAG_EMBEDDING_DEFAULTS,
  "expensive":      EXPENSIVE_MODEL_DEFAULTS,
  "cheap":          CHEAP_MODEL_DEFAULTS,
}


def default_model_for_provider(provider, role="decision"):
  """Return the agent's default model name for the given provider and role.

  Args:
    provider: one of 'google', 'openai', 'anthropic', 'ollama'.
      Case- and whitespace-insensitive.
    role: which table to consult.  One of: 'decision', 'rag',
      'rag_embedding', 'expensive', 'cheap'.  Default: 'decision'
      (matches the original v119 H1 use case in api_client and
      directive_extractor).

  Returns:
    str: the default model name from the chosen table.

  Raises:
    ValueError: if `role` is not a known role, or if `provider`
      is not a key in the chosen table.  Callers must handle
      ValueError -- falling back silently to a hardcoded string
      would defeat the centralization (and is exactly the failure
      mode v119 H1 fixed).
  """
  if role not in _ROLE_TABLES:
    raise ValueError(
      "Unknown role %r; expected one of: %s" % (
        role, ", ".join(sorted(_ROLE_TABLES))))
  provider_norm = (provider or "").lower().strip()
  table = _ROLE_TABLES[role]
  if provider_norm not in table:
    raise ValueError(
      "Provider %r not in %s table.  Available: %s" % (
        provider, role, ", ".join(sorted(table))))
  return table[provider_norm]


# =====================================================================
# v119.H13 — Ollama provider robustness helpers
# =====================================================================
# Two related fixes for AIAgent_run_39a_ollama bugs:
#
# A. URL normalization: setting OLLAMA_BASE_URL=http://localhost:11434
#    (no /v1 suffix) is reasonable but breaks the OpenAI-SDK path
#    construction.  normalize_ollama_openai_base_url() idempotently
#    appends /v1 so both bare and /v1-suffixed URLs work.
#
# B. Env-var override precedence: callers had inconsistent behavior
#    around OLLAMA_LLM_MODEL — get_llm_and_embeddings honored it, but
#    directive_extractor and api_client did not.  resolve_model_for_provider()
#    centralizes the precedence rule: env-var wins, central default
#    fallback.

def normalize_ollama_openai_base_url(base_url):
  """Ensure base_url ends in /v1 for OpenAI-compat client use.

  Ollama's OpenAI-compatibility endpoint lives at
  /v1/chat/completions.  The OpenAI Python SDK appends
  /chat/completions to the configured base_url, so base_url
  must already include /v1.  Users naturally set
  OLLAMA_BASE_URL to the bare host (matching OLLAMA_HOST or
  the Ollama docs around the native /api/* endpoints); this
  normalization makes that input work.

  Args:
    base_url: a URL string, possibly with or without trailing
      slash and with or without /v1 suffix.  None or empty
      string passes through unchanged (caller's default handling
      decides what to do).

  Returns:
    str or None: normalized URL with trailing /v1, no trailing slash.
      Idempotent for inputs already ending in /v1.

  Examples:
    'http://localhost:11434'    -> 'http://localhost:11434/v1'
    'http://localhost:11434/'   -> 'http://localhost:11434/v1'
    'http://localhost:11434/v1' -> 'http://localhost:11434/v1'
    'http://localhost:11434/v1/' -> 'http://localhost:11434/v1'
    None                        -> None
    ''                          -> ''  (passthrough; never produces bare '/v1')
  """
  if not base_url:
    return base_url
  stripped = base_url.rstrip('/')
  if stripped.endswith('/v1'):
    return stripped
  return stripped + '/v1'


# Provider → optional env-var name mapping for model overrides.
# When the env-var is set, it takes precedence over the central
# DEFAULT_MODELS table for that provider.  Matches the existing
# precedence in get_llm_and_embeddings.
_PROVIDER_MODEL_ENV_OVERRIDES = {
  "ollama": "OLLAMA_LLM_MODEL",
  # portkey forwards this model string to the Azure-OpenAI upstream.
  # Matches the user's add_portkey_support.patch (PHENIX_PORTKEY_MODEL,
  # default gpt-5 in the DEFAULTS tables above).
  "portkey": "PHENIX_PORTKEY_MODEL",
  # anthropic model override for parity with OLLAMA_LLM_MODEL.
  "anthropic": "ANTHROPIC_LLM_MODEL",
  # google/openai models are not commonly overridden via env-var in
  # this codebase; if a future need arises, add the env-var name here.
}


def resolve_model_for_provider(provider, role="decision"):
  """Resolve model name for a provider with env-var precedence.

  Precedence order:
    1. Provider-specific env-var (e.g., OLLAMA_LLM_MODEL)
    2. Central DEFAULT_MODELS table (default_model_for_provider)

  This centralizes the rule that get_llm_and_embeddings has
  always followed ("Env-var overrides take precedence; fall
  back to the central default tables when the env vars are
  unset"), making it available to other call sites in
  directive_extractor and api_client.

  Args:
    provider: same as default_model_for_provider
    role: same as default_model_for_provider (default 'decision')

  Returns:
    str: the resolved model name.

  Raises:
    ValueError: from default_model_for_provider, if provider is
      unknown.  Env-var override path raises nothing — an unset
      or empty env-var falls through to the central default.
  """
  env_var = _PROVIDER_MODEL_ENV_OVERRIDES.get(
    (provider or "").lower().strip())
  if env_var:
    env_value = os.getenv(env_var)
    if env_value:
      return env_value
  return default_model_for_provider(provider, role=role)


# =====================================================================
# v119.H13 — Provider-error classification
# =====================================================================
# Tom's run_39a_ollama failure surfaced that 404 errors fall into at
# least three operationally-distinct sub-categories:
#
#   - MODEL_RETIRED:    provider has deprecated the model.  Action:
#                       update core/llm DEFAULT_MODELS.
#   - MODEL_UNAVAILABLE: model unknown to this server (e.g., Ollama
#                       has not pulled it, or wrong name).  Action:
#                       pull the model, or fix the name.
#   - FAILED:           unrecognized 404 cause.  Action: investigate.
#
# Plus AUTH_FAILED for 401-class errors.
#
# Classifier uses programmatic property interrogation (status_code,
# response, body) plus string-content matching to handle SDK
# stringification variations.  Retirement classification uses a
# co-occurrence rule (retirement phrase must appear near "model"
# word) to avoid false-positives on edge-proxy 404 pages that
# mention "deprecated endpoints".

import re as _re_h13  # local alias; module 're' may not be imported in this file

_RETIREMENT_PHRASES = ('no longer available', 'deprecated', 'retired')
_RETIREMENT_WINDOW = 80  # chars around a phrase to scan for 'model'


def _classify_provider_error(exc, model_name):
  """Classify an LLM-call exception into a diagnostic marker tag.

  Args:
    exc: the exception object caught from a provider call.
    model_name: the model name that was being requested when
      exc was raised (for inclusion in the hint message).

  Returns:
    (marker_tag, hint) tuple where:
      marker_tag is one of:
        'DIRECTIVE_EXTRACTION_MODEL_RETIRED'
        'DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE'
        'DIRECTIVE_EXTRACTION_AUTH_FAILED'
        'DIRECTIVE_EXTRACTION_FAILED'
      hint is None or a short operator-action string.
  """
  # 1. Programmatic status-code extraction (handles SDKs that
  #    hide retirement language behind opaque __str__).
  status_code = (
    getattr(exc, 'status_code', None)
    or getattr(getattr(exc, 'response', None), 'status_code', None))

  # 2. Combine str(exc) with body/message attributes so we don't
  #    miss content hidden behind brief __str__ representations.
  err_body = ""
  if hasattr(exc, 'body') and isinstance(getattr(exc, 'body'), dict):
    err_body = str(exc.body).lower()
  elif hasattr(exc, 'message'):
    err_body = str(exc.message).lower()
  err_all = (str(exc) + " " + err_body).lower()

  is_404 = (
    status_code == 404
    or '404' in err_all
    or 'not_found' in err_all)

  is_401 = (
    status_code == 401
    or '401' in err_all
    or 'unauthor' in err_all
    or 'invalid api key' in err_all)

  if is_401:
    hint = (
      "Authentication failed for model %s. "
      "Check the provider API key." % model_name)
    return ('DIRECTIVE_EXTRACTION_AUTH_FAILED', hint)

  if is_404:
    # 3a. Retirement check — phrase must co-occur with 'model'.
    #     Distinguishes "this MODEL is deprecated" from generic
    #     "this ENDPOINT is deprecated" (edge-proxy 404 pages).
    for phrase in _RETIREMENT_PHRASES:
      for m in _re_h13.finditer(_re_h13.escape(phrase), err_all):
        window = err_all[
          max(0, m.start() - _RETIREMENT_WINDOW)
          : m.end() + _RETIREMENT_WINDOW]
        if 'model' in window:
          hint = (
            "Model %s appears to be retired by the provider. "
            "Update core/llm DEFAULT_MODELS." % model_name)
          return ('DIRECTIVE_EXTRACTION_MODEL_RETIRED', hint)

    # 3b. Unavailable check — "model 'X' not found" pattern.
    #     Indicates the model name is unknown to this host
    #     (Ollama pull missing, or env-var typo).  No retirement
    #     language present.
    if _re_h13.search(r"model\s+'[^']+'\s+not\s+found", err_all) \
        or _re_h13.search(r'model\s+"[^"]+"\s+not\s+found', err_all):
      hint = (
        "Model %s is not available on this server. "
        "For ollama: run 'ollama pull %s' or set "
        "OLLAMA_LLM_MODEL to an installed model. "
        "For cloud providers: verify the model name in "
        "core/llm DEFAULT_MODELS." % (model_name, model_name))
      return ('DIRECTIVE_EXTRACTION_MODEL_UNAVAILABLE', hint)

  return ('DIRECTIVE_EXTRACTION_FAILED', None)


def set_llm_verbosity(level):
  """Set the verbosity level for LLM-related messages.

  Args:
    level: One of 'quiet', 'normal', 'verbose', 'debug'
  """
  global _llm_verbosity
  _llm_verbosity = level


def get_llm_verbosity():
  """Get the current verbosity level."""
  global _llm_verbosity
  if _llm_verbosity is not None:
    return _llm_verbosity
  return os.getenv("LLM_VERBOSITY", "normal")


def _llm_log(msg, level='normal'):
  """Print message if verbosity level permits."""
  levels = ['quiet', 'normal', 'verbose', 'debug']
  current = get_llm_verbosity()
  current_idx = (
    levels.index(current) if current in levels else 1)
  msg_idx = (
    levels.index(level) if level in levels else 1)
  if msg_idx <= current_idx:
    print(msg)


# =====================================================================
# v120 — Portkey (Azure-OpenAI gateway) + per-provider kwarg sanitizing
# =====================================================================
# Portkey is an OpenAI-SDK-compatible gateway.  The user's deployment
# (add_portkey_support.patch) fronts Azure OpenAI via the Portkey SDK
# with provider="azure-openai".  Two env vars are required at call time:
#   PORTKEY_AZURE_API_KEY  -- the Portkey gateway key
#   PORTKEY_BASE_URL       -- the Portkey gateway base URL
# and the upstream model is chosen by PHENIX_PORTKEY_MODEL (default gpt-5,
# resolved via resolve_model_for_provider / the DEFAULTS tables).

PORTKEY_DEFAULT_MODEL = "gpt-5"


def portkey_langchain_config():
  """Return (api_key, base_url) for the langchain ChatOpenAI/Embeddings
  portkey path, reading env vars at CALL time (never import time).

  Raises a clear ValueError naming the missing variable rather than a
  bare KeyError (the patch used os.environ[...] which raises KeyError).
  """
  api_key = os.getenv("PORTKEY_AZURE_API_KEY")
  base_url = os.getenv("PORTKEY_BASE_URL")
  missing = [name for name, val in (
    ("PORTKEY_AZURE_API_KEY", api_key),
    ("PORTKEY_BASE_URL", base_url)) if not val]
  if missing:
    raise ValueError(
      "Provider 'portkey' requires the environment variable(s): %s. "
      "Set them before running with provider=portkey." % ", ".join(missing))
  return api_key, base_url


def anthropic_model_rejects_sampling_params(model_name):
  """True if the given Anthropic model rejects temperature/top_p/top_k.

  The Claude 4.6 / 4.7 / 4.8 reasoning family and every Claude 5.x model
  (claude-sonnet-5-5, claude-opus-5-5, claude-opus-5, claude-haiku-5-5,
  claude-fable-5-1, ...) reject the sampling parameters with HTTP 400;
  adaptive thinking controls randomness instead.  Older Claude models
  (Sonnet 4 / Opus 4 / Sonnet 4.5 / Haiku 4.5 / 3.x) still accept
  temperature.  Detected by model-name family rather than a list, so later
  models in either generation are covered without another edit.

  Matches claude-<family>-<major>[-<minor>] with major >= 5, or major 4 with
  minor >= 6.  A dated snapshot (claude-opus-4-20250514) has an 8-digit
  "minor" and is the older, temperature-accepting form.
  """
  m = (model_name or "").lower()
  if "claude" not in m:
    return False
  import re as _re
  match = _re.search(r"claude-[a-z]+-(\d+)(?:-(\d+))?", m)
  if not match:
    return False
  major = int(match.group(1))
  minor = match.group(2)
  if major >= 5 and len(match.group(1)) < 5:
    return True
  if major != 4 or minor is None or len(minor) >= 5:
    return False
  return int(minor) >= 6


# Gemini finish reasons meaning the reply was withheld or stopped by a safety
# or policy filter (not cut off by the token budget).
GEMINI_BLOCKED_FINISH_REASONS = (
  "SAFETY", "RECITATION", "PROHIBITED_CONTENT", "BLOCKLIST", "SPII")


def _gemini_blocked_reason(finish):
  """The GEMINI_BLOCKED_FINISH_REASONS name for a finish_reason, or None.

  Accepts an enum (its .name), "FinishReason.SAFETY" or a plain string.
  """
  name = str(getattr(finish, "name", finish) or "").upper().split(".")[-1]
  return name if name in GEMINI_BLOCKED_FINISH_REASONS else None


# Gemini prompt_feedback.block_reason values meaning the prompt itself was
# blocked, so no candidate (reply) was produced.  Normal replies carry 0 or
# BLOCKED_REASON_UNSPECIFIED.
GEMINI_BLOCKED_PROMPT_REASONS = (
  "SAFETY", "OTHER", "BLOCKLIST", "PROHIBITED_CONTENT", "IMAGE_SAFETY",
  "MODEL_ARMOR", "JAILBREAK")


def _gemini_prompt_block_reason(feedback):
  """Name of the reason Gemini blocked the prompt, or None if it did not.

  feedback is a google.genai PromptFeedback (or None), or the dict
  langchain_google_genai puts in response_metadata["prompt_feedback"].
  block_reason may be an enum, a name ("SAFETY", "BlockedReason.SAFETY")
  or a number; 0, None and BLOCKED_REASON_UNSPECIFIED mean not blocked.
  A nonzero number or unknown name is reported as "UNRECOGNIZED" so that
  the message carries no digits; any other value type is ignored.
  """
  if feedback is None:
    return None
  if isinstance(feedback, dict):
    reason = feedback.get("block_reason")
  else:
    reason = getattr(feedback, "block_reason", None)
  if reason is None or isinstance(reason, bool):
    return None
  name_attr = getattr(reason, "name", None)
  if isinstance(name_attr, str):
    name = name_attr
  elif isinstance(reason, str):
    name = reason
  elif isinstance(reason, int):
    return "UNRECOGNIZED" if reason else None
  else:
    return None    # not a recognizable value (e.g. a test double)
  name = name.upper().split(".")[-1]
  if not name or name in ("0", "BLOCKED_REASON_UNSPECIFIED"):
    return None
  return name if name in GEMINI_BLOCKED_PROMPT_REASONS else "UNRECOGNIZED"


def _raise_if_langchain_reply_cut_off(response):
  """Raise RuntimeError if a LangChain message says it was cut off.

  LangChain puts the provider's own stop reason in response_metadata:
  Anthropic stop_reason "max_tokens" / "model_context_window_exceeded",
  OpenAI finish_reason "length", Gemini finish_reason "MAX_TOKENS", Ollama
  done_reason "length".  A cut-off reply used to pass as complete text.
  OpenAI finish_reason "content_filter" (reply withheld or truncated by the
  provider's filter) raises too, with its own message, as do Anthropic
  stop_reason "refusal" and the Gemini safety finish reasons
  (GEMINI_BLOCKED_FINISH_REASONS).  The messages carry no numbers
  (rate_limit_handler would read "429"/"503" in them as a rate limit).
  """
  meta = getattr(response, "response_metadata", None)
  if not isinstance(meta, dict):
    return
  stop_reason = meta.get("stop_reason")
  finish = meta.get("finish_reason")
  finish_name = str(getattr(finish, "name", finish) or "").upper()
  done_reason = meta.get("done_reason")
  if stop_reason == "refusal":
    raise RuntimeError(
      "LLM declined to answer (stop_reason='refusal')")
  prompt_blocked = _gemini_prompt_block_reason(meta.get("prompt_feedback"))
  if prompt_blocked:
    raise RuntimeError(
      "LLM prompt was blocked by the provider's safety filter "
      "(block_reason=%s)" % prompt_blocked)
  blocked = _gemini_blocked_reason(finish)
  if blocked:
    raise RuntimeError(
      "LLM reply was blocked by the provider's safety filter "
      "(finish_reason=%s)" % blocked)
  if finish_name == "CONTENT_FILTER":
    raise RuntimeError(
      "LLM reply was blocked by the provider's content filter "
      "(finish_reason='content_filter')")
  if (stop_reason in ANTHROPIC_CUTOFF_STOP_REASONS
      or finish_name in ("LENGTH", "MAX_TOKENS")
      or done_reason == "length"):
    raise RuntimeError(
      "LLM reply was cut off (stop_reason=%r, finish_reason=%r, "
      "done_reason=%r); raise the output budget or shorten the request"
      % (stop_reason, finish_name.lower() or None, done_reason))


def response_text(response):
  """Plain text of an LLM reply, whatever shape the provider returned.

  Accepts a LangChain message (``.content`` is a string, or a list of
  content blocks when the model also returned thinking blocks, as Claude
  5.x does), an Anthropic SDK Message (``.content`` is a list of typed
  blocks; the first may be a thinking block), or a plain string.  Returns
  the concatenated text blocks; thinking and tool blocks are skipped.
  Never raises on an unexpected shape: falls back to str().
  """
  if response is None:
    return ""
  if isinstance(response, str):
    return response
  _raise_if_langchain_reply_cut_off(response)
  content = getattr(response, "content", response)
  if content is None:
    return ""
  if isinstance(content, str):
    return content
  if isinstance(content, (list, tuple)):
    parts = []
    for block in content:
      if isinstance(block, str):
        parts.append(block)
      elif isinstance(block, dict):
        if block.get("type", "text") == "text":
          parts.append(block.get("text") or "")
      elif getattr(block, "type", None) == "text":
        parts.append(getattr(block, "text", "") or "")
    return "".join(parts)
  return str(content)


# Output-token cap for Claude requests.  On Claude 5.x models thinking is on
# by default and its tokens count toward max_tokens, so the cap must leave
# room for thinking plus the answer.  Measured 2026-10-14 (live probe): a
# short reasoning task used 2000-2900 output tokens (1000-1800 of them
# thinking) on claude-sonnet-5-5 / claude-opus-5-5, and a 2000 cap truncated
# both.  Only tokens actually generated are billed, so a higher cap costs
# nothing unless it is used.  Upper limit: these requests are not streamed,
# and the anthropic SDK (0.122) refuses a non-streaming request whose
# max_tokens implies more than 10 minutes (above 21333 tokens) when the
# client keeps its default timeout, so do not raise this past
# ANTHROPIC_NONSTREAMING_MAX_TOKENS without switching to streaming.  That
# ceiling applies only while the client uses the SDK's default timeout (both
# direct calls do); the SDK also caps a few retired Opus 4 / 4.1 models at
# 8192 non-streaming tokens, which this floor would exceed.
ANTHROPIC_MAX_TOKENS = 16000
ANTHROPIC_NONSTREAMING_MAX_TOKENS = 21333

# stop_reason values that mean the reply was cut off before it finished.
ANTHROPIC_CUTOFF_STOP_REASONS = ("max_tokens", "model_context_window_exceeded")


def anthropic_request_max_tokens(requested):
  """max_tokens to send on a non-streaming Claude request.

  At least ANTHROPIC_MAX_TOKENS (thinking shares the budget) and at most
  ANTHROPIC_NONSTREAMING_MAX_TOKENS (the SDK refuses larger non-streaming
  requests), whatever the caller asked for.
  """
  return min(max(int(requested or 0), ANTHROPIC_MAX_TOKENS),
             ANTHROPIC_NONSTREAMING_MAX_TOKENS)


def anthropic_reply_text(response):
  """Text of an Anthropic SDK reply; raises RuntimeError if it is unusable.

  On Claude 5.x models thinking is on by default and its tokens count
  toward max_tokens, so a reply can stop (stop_reason "max_tokens", or
  "model_context_window_exceeded" when the context window fills) with no
  text or with the text cut off.  Raising keeps that failure visible through
  the callers' existing error handling instead of returning "" or a
  truncated answer silently.
  """
  text = response_text(response)
  stop_reason = getattr(response, "stop_reason", None)
  if stop_reason == "refusal":
    raise RuntimeError(
      "Anthropic model declined to answer (stop_reason='refusal')")
  if stop_reason in ANTHROPIC_CUTOFF_STOP_REASONS:
    raise RuntimeError(
      "Anthropic reply was cut off (stop_reason=%r); on Claude 5.x models "
      "thinking counts toward max_tokens" % (stop_reason,))
  if not text.strip():
    raise RuntimeError(
      "Anthropic reply contained no text (stop_reason=%r); on Claude 5.x "
      "models thinking counts toward max_tokens"
      % (getattr(response, "stop_reason", None),))
  return text


# Output budget for the direct (non-LangChain) Google and OpenAI calls.
# Thinking / reasoning models spend output tokens before they answer: a live
# probe on 2026-10-14 saw gemini-2.5-pro and gpt-6-luna use all of a 2000
# token cap on a short reasoning task and return no text.  Only generated
# tokens are billed, so the larger cap costs nothing unless it is used.
REASONING_MAX_TOKENS = 16000

# Minimum request timeout (seconds) for hosted-provider LangChain clients.
HOSTED_LLM_MIN_TIMEOUT = 180


def openai_model_is_reasoning(model_name):
  """True for OpenAI reasoning models (gpt-5*, gpt-6* and later, o1/o3/o4).

  These reject a non-default temperature (HTTP 400 "Only the default (1)
  value is supported") and take max_completion_tokens.  The *-chat variants
  are not reasoning models.  langchain-openai 1.5 strips temperature only
  for gpt-5* names, so this check is needed for gpt-6* and later.
  """
  import re as _re
  m = (model_name or "").lower().strip()
  if "chat" in m:
    return False
  # Allow gateway / fine-tune prefixes: "openai/gpt-6-luna", "azure/gpt-5",
  # "ft:gpt-5-mini:org::id".
  match = _re.search(r"(?:^|[/:])gpt-(\d+)", m)
  if match:
    return int(match.group(1)) >= 5
  return bool(_re.search(r"(?:^|[/:])o\d", m))


def google_model_ignores_sampling(model_name):
  """True for Gemini 3 and later, which use fixed sampling defaults.

  They ignore temperature/top_p/top_k (langchain-google-genai warns), so
  PHENIX does not send them.
  """
  import re as _re
  match = _re.search(r"gemini-(\d+)", (model_name or "").lower())
  return bool(match) and int(match.group(1)) >= 3


def openai_reply_text(response):
  """Text of an OpenAI chat completion; raises RuntimeError if unusable.

  The error messages carry no numbers: rate_limit_handler treats any message
  containing "429" or "503" as a rate limit and would retry a cut-off call.

  finish_reason "length" means the reply was cut off (on reasoning models
  the reasoning tokens can use the whole budget and leave no text).
  finish_reason "content_filter" means the provider's filter withheld or
  truncated the reply.
  """
  choice = response.choices[0]
  text = getattr(choice.message, "content", None) or ""
  finish = getattr(choice, "finish_reason", None)
  if finish == "content_filter":
    raise RuntimeError(
      "OpenAI reply was blocked by the content filter "
      "(finish_reason='content_filter')")
  if finish == "length":
    raise RuntimeError(
      "OpenAI reply was cut off (finish_reason='length')")
  if not text.strip():
    raise RuntimeError(
      "OpenAI reply contained no text (finish_reason=%r)" % (finish,))
  return text


def google_reply_text(response):
  """Text of a Gemini generate_content reply; raises RuntimeError if unusable.

  finish_reason MAX_TOKENS means the reply was cut off (thinking can use the
  whole budget and leave no text); the GEMINI_BLOCKED_FINISH_REASONS mean a
  safety or policy filter withheld or stopped it.
  """
  candidates = getattr(response, "candidates", None) or []
  finish = getattr(candidates[0], "finish_reason", None) if candidates else None
  finish_name = getattr(finish, "name", None) or str(finish or "")
  prompt_blocked = _gemini_prompt_block_reason(
    getattr(response, "prompt_feedback", None))
  if prompt_blocked:
    raise RuntimeError(
      "Gemini prompt was blocked by a safety filter (block_reason=%s)"
      % prompt_blocked)
  try:
    text = response.text or ""
  except (ValueError, AttributeError):
    text = ""
  blocked = _gemini_blocked_reason(finish)
  if blocked:
    raise RuntimeError(
      "Gemini reply was blocked by a safety filter (finish_reason=%s)" % blocked)
  if "MAX_TOKENS" in finish_name.upper():
    raise RuntimeError(
      "Gemini reply was cut off (finish_reason=MAX_TOKENS)")
  if not text.strip():
    raise RuntimeError(
      "Gemini reply contained no text (finish_reason=%s)" % (finish_name or None,))
  return text


def google_generation_kwargs(model_name, temperature, max_tokens):
  """Generation settings for a direct Gemini call.

  Thinking models spend output tokens before answering, so the budget is
  at least REASONING_MAX_TOKENS; Gemini 3+ ignores temperature, so it is
  not sent there.
  """
  kwargs = {"max_output_tokens": max(int(max_tokens or 0), REASONING_MAX_TOKENS)}
  if not google_model_ignores_sampling(model_name):
    kwargs["temperature"] = temperature
  return kwargs


def openai_completion_kwargs(model_name, temperature, max_tokens):
  """Sampling and length settings for a direct OpenAI chat completion.

  Reasoning models reject a non-default temperature and need room for
  their reasoning tokens (max_completion_tokens); other models keep
  temperature and max_tokens as before.
  """
  if openai_model_is_reasoning(model_name):
    return {"max_completion_tokens": max(int(max_tokens or 0), REASONING_MAX_TOKENS)}
  return {"temperature": temperature, "max_tokens": max_tokens}


def portkey_max_completion_tokens(model_name, max_tokens):
  """max_completion_tokens for a Portkey (Azure OpenAI) chat call.

  Azure gpt-5 is a reasoning model: its reasoning tokens share the budget, so
  reasoning deployments get at least REASONING_MAX_TOKENS.
  """
  if openai_model_is_reasoning(model_name):
    return max(int(max_tokens or 0), REASONING_MAX_TOKENS)
  return max_tokens


def sanitize_llm_kwargs(provider, kwargs):
  """Reshape a shared kwargs dict to each provider's accepted schema.

  Keyed by provider so we never apply a blanket strip.  This exists so a
  *shared* kwargs dict (temperature / max_tokens assembled once) can't crash
  a strict client.  Returns a NEW dict; never mutates the caller's.

    - portkey  (Azure-OpenAI upstream, gpt-5): rename max_tokens ->
               max_completion_tokens, drop temperature (Azure gpt-5 rejects
               both forms).
    - anthropic: keep max_tokens (ChatAnthropic REQUIRES it).  For
               temperature, it depends on the model: the 4.6+/5.x reasoning
               family REJECTS temperature (HTTP 400), so drop it for those;
               older Claude models keep it.  Detected via
               anthropic_model_rejects_sampling_params() on the 'model' kwarg.
    - others (openai/google/ollama): unchanged.

  Known limitation: the portkey rule assumes the Azure-OpenAI upstream in
  use today.  If portkey is ever pointed at a non-Azure upstream that needs
  temperature, this rule must become upstream-aware (e.g. keyed on
  PHENIX_PORTKEY_MODEL).
  """
  clean = dict(kwargs)
  prov = (provider or "").lower().strip()
  if prov == "portkey":
    if "max_tokens" in clean:
      clean["max_completion_tokens"] = clean.pop("max_tokens")
    clean.pop("temperature", None)
  elif prov == "anthropic":
    # Keep max_tokens (required).  Drop temperature only for models that
    # reject it (4.6+/5.x reasoning family); older models keep it.
    model_name = clean.get("model") or clean.get("model_name")
    if anthropic_model_rejects_sampling_params(model_name):
      clean.pop("temperature", None)
      clean.pop("top_p", None)
      clean.pop("top_k", None)
  return clean


def _delegate_embeddings_for_nonnative(
  requesting_provider, embedding_model_name=None,
  batch_size=100, timeout=120):
  """Embeddings fallback for providers with no native embeddings endpoint
  (currently anthropic).

  Policy (v120 §4.6): never raise on construction.  Delegate to a default
  embeddings provider when a key is available; otherwise return None.  The
  agent's end-of-run assessment never queries embeddings, so None is safe
  there; only a real RAG query needs a working backend, and that path
  surfaces its own clear error if none exists.
  """
  # Prefer OpenAI, then Google, based on which key is present.
  if os.getenv("OPENAI_API_KEY"):
    try:
      from langchain_openai import OpenAIEmbeddings
      model = embedding_model_name or default_model_for_provider(
        "openai", role="rag_embedding")
      _llm_log(
        "[PROVIDER_DELEGATION] Fallback embedding triggered under "
        "non-native provider context (%s -> openai)" % requesting_provider,
        level='verbose')
      return OpenAIEmbeddings(model=model, chunk_size=batch_size)
    except Exception:
      pass  # fall through to next option / None
  if os.getenv("GOOGLE_API_KEY"):
    try:
      from langchain_google_genai import GoogleGenerativeAIEmbeddings
      model = embedding_model_name or default_model_for_provider(
        "google", role="rag_embedding")
      _llm_log(
        "[PROVIDER_DELEGATION] Fallback embedding triggered under "
        "non-native provider context (%s -> google)" % requesting_provider,
        level='verbose')
      return GoogleGenerativeAIEmbeddings(
        model=model, timeout=timeout, batch_size=batch_size,
        google_api_key=os.getenv("GOOGLE_API_KEY"))
    except Exception:
      pass
  # No delegate available: None is correct (never raise here).
  _llm_log(
    "[PROVIDER_DELEGATION] No embeddings backend available for "
    "non-native provider %s; returning None (RAG unavailable, chat "
    "unaffected)" % requesting_provider,
    level='verbose')
  return None


def get_llm_and_embeddings(
  provider: str = None,
  llm_model_name: str = None,
  embedding_model_name: str = None,
  temperature: float = 0.0,
  timeout: int = 120,
  batch_size: int = 100,
  ollama_base_url: str = None,
  json_mode: bool = False,
  num_ctx: int = 8192,
  seed: int = 42,
):
  """
  Initialize LLM and Embeddings from the specified provider.

  Args:
    provider: Which AI provider to use
      ('ollama', 'google', or 'openai'). If None, reads
      from LLM_PROVIDER env var, defaults to 'ollama'.
    llm_model_name: Specific model name, or None for default
    embedding_model_name: Specific embedding model, or None
      for default
    temperature: Sampling temperature
      (0.0 = deterministic, 1.0 = creative)
    timeout: Request timeout in seconds
    batch_size: Batch size for embedding requests
    ollama_base_url: Base URL for Ollama server
      (default from env or localhost)
    json_mode: If True, force JSON output (for structured
      responses). Only applies to Ollama provider.
    num_ctx: Context window size for Ollama models
      (default: 8192)

  Returns:
    tuple: (llm, embeddings) - initialized model objects

  Raises:
    ValueError: If provider is not supported
  """
  # Get provider from argument or environment variable
  if provider is None:
    provider = os.getenv("LLM_PROVIDER", "ollama")
  provider = provider.lower()
  # Hosted reasoning/thinking models can take well over a minute on long
  # prompts (live replies of 31-41 s on short ones, 2026-10-14), so never give
  # their clients less than HOSTED_LLM_MIN_TIMEOUT.  A timeout is a ceiling:
  # raising it only lengthens the wait for a call that would otherwise fail.
  if provider in ("anthropic", "google", "openai", "portkey"):
    timeout = max(timeout or 0, HOSTED_LLM_MIN_TIMEOUT)

  if provider == "ollama":
    from langchain_ollama import ChatOllama, OllamaEmbeddings

    # Env-var overrides take precedence; fall back to the central
    # default tables when the env vars are unset.
    # v119.H13: this inline pattern is the contract that
    # resolve_model_for_provider() centralizes for the OpenAI-compat
    # call sites in directive_extractor and api_client.  We leave
    # this inline because (a) it uses role='rag' not the helper's
    # default 'decision', and (b) the embed model has its own
    # env-var (OLLAMA_EMBED_MODEL) that the helper doesn't cover.
    # If the helper grows role + multi-env-var support, this
    # can be replaced with resolve_model_for_provider calls.
    if llm_model_name is None:
      llm_model_name = (
        os.getenv("OLLAMA_LLM_MODEL")
        or default_model_for_provider("ollama", role="rag"))
    if embedding_model_name is None:
      embedding_model_name = (
        os.getenv("OLLAMA_EMBED_MODEL")
        or default_model_for_provider("ollama", role="rag_embedding"))
    if ollama_base_url is None:
      ollama_base_url = os.getenv(
        "OLLAMA_BASE_URL", "http://localhost:11434")
      # NB: this path uses ChatOllama which hits Ollama's NATIVE
      # /api/chat endpoint — no /v1 needed.  Do NOT apply
      # normalize_ollama_openai_base_url() here (it would break
      # native ChatOllama).  /v1 normalization is ONLY for the
      # OpenAI-compat SDK paths in api_client._call_ollama_llm
      # and directive_extractor._call_llm_fallback.

    # Build kwargs conditionally
    llm_kwargs = {
      "model": llm_model_name,
      "base_url": ollama_base_url,
      "temperature": temperature,
      "num_ctx": num_ctx,
      "num_predict": 4096,
      "seed": seed,
      "think": False,
    }
    if json_mode:
      llm_kwargs["format"] = "json"

    llm = ChatOllama(**llm_kwargs)
    embeddings = OllamaEmbeddings(
      model=embedding_model_name,
      base_url=ollama_base_url,
    )
    mode_str = " (JSON mode)" if json_mode else ""
    _llm_log(
      f"Using Ollama at {ollama_base_url}{mode_str}",
      level='verbose')
    _llm_log(
      f"  LLM: {llm_model_name},"
      f" Embeddings: {embedding_model_name}",
      level='verbose')


  elif provider == "google":
    from langchain_google_genai import (
      ChatGoogleGenerativeAI,
      GoogleGenerativeAIEmbeddings,
    )

    if llm_model_name is None:
      llm_model_name = default_model_for_provider(
        "google", role="rag")
    if embedding_model_name is None:
      # gemini-embedding-001 is the current stable model
      # (released July 2025)
      # Note: text-embedding-004 is legacy and being
      # deprecated -- see RETIRED_MODELS for the policed list.
      embedding_model_name = default_model_for_provider(
        "google", role="rag_embedding")

    google_kwargs = {
      "model": llm_model_name,
      "timeout": timeout,
      "max_retries": 0,
    }
    # Gemini 3+ ignores sampling parameters; do not send them.
    if not google_model_ignores_sampling(llm_model_name):
      google_kwargs["temperature"] = temperature
    llm = ChatGoogleGenerativeAI(**google_kwargs)
    embeddings = GoogleGenerativeAIEmbeddings(
      model=embedding_model_name,
      timeout=timeout,
      batch_size=batch_size,
      google_api_key=os.getenv("GOOGLE_API_KEY")
    )
    _llm_log(
      f"Using Google Gemini: {llm_model_name}",
      level='verbose')

  elif provider == "openai":
    from langchain_openai import ChatOpenAI, OpenAIEmbeddings

    if llm_model_name is None:
      llm_model_name = default_model_for_provider(
        "openai", role="rag")
    if embedding_model_name is None:
      embedding_model_name = default_model_for_provider(
        "openai", role="rag_embedding")

    openai_kwargs = {
      "model": llm_model_name,
      "timeout": timeout,
      "max_retries": 2,
    }
    # Reasoning models reject a non-default temperature (HTTP 400).
    if not openai_model_is_reasoning(llm_model_name):
      openai_kwargs["temperature"] = temperature
    llm = ChatOpenAI(**openai_kwargs)
    embeddings = OpenAIEmbeddings(
      model=embedding_model_name,
      chunk_size=batch_size
    )
    _llm_log(
      f"Using OpenAI: {llm_model_name}",
      level='verbose')

  elif provider == "portkey":
    # Portkey is OpenAI-SDK-compatible; langchain's ChatOpenAI /
    # OpenAIEmbeddings talk to the Portkey gateway by pointing base_url at
    # it with the Portkey key.  (The raw-SDK call path in api_client /
    # directive_extractor uses the Portkey SDK directly with
    # provider="azure-openai"; this langchain path is the embeddings-capable
    # equivalent used by get_llm_and_embeddings / RAG.)
    # Check env FIRST (clear ValueError) before importing the SDK, so a
    # misconfigured deployment gets the friendly message rather than an
    # import error.
    api_key, base_url = portkey_langchain_config()  # call-time env check
    from langchain_openai import ChatOpenAI, OpenAIEmbeddings

    if llm_model_name is None:
      llm_model_name = resolve_model_for_provider("portkey", role="rag")
    if embedding_model_name is None:
      embedding_model_name = default_model_for_provider(
        "portkey", role="rag_embedding")

    # gpt-5 via Azure rejects temperature and uses max_completion_tokens;
    # ChatOpenAI is constructed without temperature accordingly.  We do not
    # pass max_tokens here (langchain default), so no remap is needed at
    # construction; sanitize_llm_kwargs covers any shared-kwargs call path.
    llm = ChatOpenAI(
      model=llm_model_name,
      api_key=api_key,
      base_url=base_url,
      timeout=timeout,
      max_retries=2,
    )
    # Embeddings object construction is lazy (no network call), so this is
    # safe even when the current key lacks embeddings access -- the agent's
    # end-of-run assessment constructs but never queries it.  An actual RAG
    # query against a no-access key surfaces a clear error at query time
    # (handled in the RAG path), never here.
    embeddings = OpenAIEmbeddings(
      model=embedding_model_name,
      api_key=api_key,
      base_url=base_url,
      chunk_size=batch_size,
    )
    _llm_log(
      f"Using Portkey (azure-openai): {llm_model_name}",
      level='verbose')

  elif provider == "anthropic":
    # Anthropic has no native embeddings endpoint.  Return a working chat
    # LLM and delegate embeddings to a default provider when one is
    # available, else None (the agent's end-of-run assessment never queries
    # embeddings, so None is harmless there).  Construction must never raise
    # on the embeddings side.
    from langchain_anthropic import ChatAnthropic

    if llm_model_name is None:
      llm_model_name = resolve_model_for_provider(
        "anthropic", role="expensive")

    # ChatAnthropic REQUIRES max_tokens.  temperature is model-dependent: the
    # 4.6+/5.x reasoning family rejects it (HTTP 400 "temperature is
    # deprecated for this model"), so only pass it for models that accept it.
    anthropic_kwargs = {
      "model": llm_model_name,
      "timeout": timeout,
      "max_tokens": ANTHROPIC_MAX_TOKENS,
      "max_retries": 2,
    }
    if not anthropic_model_rejects_sampling_params(llm_model_name):
      anthropic_kwargs["temperature"] = temperature
    llm = ChatAnthropic(**anthropic_kwargs)
    embeddings = _delegate_embeddings_for_nonnative(
      requesting_provider="anthropic",
      embedding_model_name=embedding_model_name,
      batch_size=batch_size,
      timeout=timeout,
    )
    _llm_log(
      f"Using Anthropic: {llm_model_name}",
      level='verbose')

  else:
    raise ValueError(
      f"Unsupported provider: '{provider}'."
      " Choose 'ollama', 'google', 'openai', 'anthropic', or 'portkey'.")

  return llm, embeddings


def verify_embeddings(embeddings, embedding_model_name=None, provider=None,
                      log=None):
  """Preflight a freshly-constructed embeddings object with ONE real call, to
  surface a silent fallback to the wrong model BEFORE building an entire vector
  database with corrupt vectors.

  Motivating failure: the portkey/Azure gateway, given an unknown embedding
  DEPLOYMENT name, does not error -- it silently routes to a default chat
  deployment (e.g. gpt-5-mini).  The embeddings object constructs fine (lazy,
  no network) and only the actual embed call reveals the problem.  Building the
  database would then "succeed" with meaningless vectors.

  This helper embeds a short probe string and prints a LOUD multi-line WARNING
  (it does NOT raise / does NOT abort the build) when something looks wrong:
    * the embed call failed, or returned nothing usable (gateway-fallback /
      no-access), or
    * the model has a known fixed dimension (EMBEDDING_EXPECTED_DIM -- the
      OpenAI 1536 family) and the returned dimension does not match (the
      gpt-5-mini-fallback case).
  Models with configurable/ambiguous dimensions (gemini, nomic) are checked
  only for "returned a non-empty numeric vector", never for an exact size, so
  this never warns spuriously for them.

  Warn-don't-fail is deliberate: a hard failure here could block a legitimate
  build if the expected-dimension table is ever stale or a provider changes
  dimensions.  The warning is loud enough to catch the eye in the build log
  without ever stopping a good run.

  Args:
    embeddings: a constructed embeddings object (LangChain Embeddings API).
    embedding_model_name: the requested model name (for the message and the
      dimension lookup).  If None, the dimension check is skipped.
    provider: provider name, for the message only.
    log: optional callable(str) for the success line (defaults to print).

  Returns:
    (ok, observed_dim): ok is False if a warning was emitted, True otherwise;
    observed_dim is the observed dimension or None if the call failed.
  """
  _emit = log if callable(log) else print

  def _warn(msg):
    bar = "!" * 70
    _emit("\n" + bar)
    _emit("WARNING: EMBEDDING SELF-CHECK (%s)" % who)
    for line in msg.splitlines():
      _emit("  " + line)
    _emit(bar + "\n")

  who = "provider=%s model=%s" % (provider or "?", embedding_model_name or "?")
  try:
    vec = embeddings.embed_query("phenix embedding self-check")
  except Exception as e:
    _warn(
      "The embedding call raised: %s: %s\n"
      "This usually means the embedding model/deployment name is wrong or not "
      "available on this endpoint.\n"
      "For portkey/Azure the DEPLOYMENT name must match EXACTLY -- e.g. "
      "'text-embedding-3-small-1' (with the trailing -1).\n"
      "If the build continues, the database may be empty or corrupt."
      % (type(e).__name__, e))
    return (False, None)

  # Validate we got a usable numeric vector.  Accept any non-string sequence
  # (list, tuple, or e.g. a numpy array) of numbers -- some embedding backends
  # return numpy arrays, and rejecting those would be a FALSE warning.  We probe
  # length and the first few elements defensively.
  bad_vector = False
  try:
    n = len(vec)
    if isinstance(vec, (str, bytes)) or n == 0:
      bad_vector = True
    else:
      head = list(vec[:8])
      try:
        # Coerce rather than isinstance: numpy scalars (np.float32) are not
        # Python int/float but are valid embedding components.
        for x in head:
          float(x)
      except (TypeError, ValueError):
        bad_vector = True
  except Exception:
    bad_vector = True   # not sequence-like at all
  if bad_vector:
    _warn(
      "The embedding call returned no usable vector (got: %r...).\n"
      "This is the signature of a SILENT gateway fallback to a non-embedding "
      "model.\nCheck the embedding model/deployment name before trusting the "
      "database." % (str(vec)[:80],))
    return (False, None)

  observed = len(vec)
  expected = EMBEDDING_EXPECTED_DIM.get(embedding_model_name)
  if expected is not None and observed != expected:
    _warn(
      "Expected a %d-dimensional vector for '%s' but got %d.\n"
      "This STRONGLY indicates the request silently fell back to a DIFFERENT "
      "model\n(e.g. the portkey/Azure gateway routing an unknown embedding "
      "deployment to a\nchat model such as gpt-5-mini).\n"
      "Verify the embedding model/deployment name before trusting the database."
      % (expected, embedding_model_name, observed))
    return (False, observed)

  _emit("Embedding self-check OK (%s, dim=%d)" % (who, observed))
  return (True, observed)

# Get expensive or cheap LLMs

def get_expensive_llm(
  provider = None, timeout = None,
  json_mode=False):
  try:
    if provider not in EXPENSIVE_MODEL_DEFAULTS:
      raise ValueError(
        "Sorry, unable to set up LLM."
        " No expensive default for provider (%s)" % provider)
    expensive_model = default_model_for_provider(
      provider, role="expensive")
    if provider == "google":
      expensive_llm, embeddings = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=expensive_model)
      _llm_log(
        "Using expensive model for analysis:"
        f" {expensive_llm.model}",
        level='verbose')
    elif provider == "openai":
      expensive_llm, embeddings = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=expensive_model)
      _llm_log(
        "Using expensive model for analysis:"
        f" {expensive_llm.model_name}",
        level='verbose')
    elif provider == "ollama":
      # Historical alternates explored during model evaluation
      # (kept as comments so the rationale is auditable):
      #   'llama3.1:70b', 'qwen2.5:72b', 'llama3.1:405b',
      #   'qwen2.5:72b'
      # Live setting comes from EXPENSIVE_MODEL_DEFAULTS["ollama"].
      expensive_llm, embeddings = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=expensive_model,
        json_mode=json_mode)
      _llm_log(
        "Using expensive model for analysis:"
        f" {expensive_llm.model}",
        level='verbose')
    elif provider == "portkey":
      expensive_llm, embeddings = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=expensive_model)
      _llm_log(
        "Using expensive model for analysis:"
        f" {expensive_llm.model_name}",
        level='verbose')
    elif provider == "anthropic":
      expensive_llm, embeddings = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=expensive_model)
      _llm_log(
        "Using expensive model for analysis:"
        f" {expensive_llm.model}",
        level='verbose')
    else:
      raise ValueError(
        "Sorry, unable to set up LLM."
        " Check llm provider (%s)" % provider)
  except Exception as e:
    raise ValueError(
      "Sorry, unable to set up LLM."
      " API keys for provider (%s)" % provider)
  return expensive_llm, embeddings

def get_cheap_llm(provider=None, timeout=None):
  if provider is None:
    provider = os.getenv("LLM_PROVIDER", "ollama")
  provider = (provider or "").lower().strip()

  try:
    # Remote providers fall through to RAG default (no entry in
    # CHEAP_MODEL_DEFAULTS); ollama overrides to a smaller local
    # model.  Historical alternates explored:
    #   'qwen2.5:14b', 'qwen2.5:32b', 'qwen3:32b', 'qwen2.5:72b',
    #   'llama3.1:8b'.  Live setting in CHEAP_MODEL_DEFAULTS["ollama"].
    if provider in CHEAP_MODEL_DEFAULTS:
      cheap_model = default_model_for_provider(
        provider, role="cheap")
      cheap_llm, _ = get_llm_and_embeddings(
        provider=provider, timeout=timeout,
        llm_model_name=cheap_model)
    else:
      cheap_llm, _ = get_llm_and_embeddings(
        provider=provider, timeout=timeout)

    # Handle different attribute names across providers
    model_name = (
      getattr(cheap_llm, 'model', None)
      or getattr(cheap_llm, 'model_name', 'unknown'))
    _llm_log(
      "Using cheap/fast model for summarization:"
      f" {model_name}",
      level='verbose')

    return cheap_llm

  except ValueError as e:
    _llm_log(str(e), level='quiet')
    raise ValueError(
      "Sorry, unable to set up LLM (%s)."
      " Check API keys." % provider)
