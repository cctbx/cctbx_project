"""Profile loader and Profile dataclass.

Section 13 of the design spec. JSON files; cascading override
(project > user > built-in); one level of based_on inheritance with cycle
detection; ${VAR} and ${env:NAME} expansion in mcp_servers strings.
"""

import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from libtbx.utils import Sorry


# Set of keys we know about at the top level. Unknown keys → warning only
# (forward compatibility).
_KNOWN_KEYS = {
  "name", "description", "based_on",
  "model", "max_tokens", "thinking", "vision_input",
  "system_prompt", "system_prompt_file",
  "skills", "mcp_servers", "tool_policy_default",
  "subagents",
  "backend", "claude",
}

# Backends supported by the agent factory. Validated at load time so
# typos surface before the chat window opens.
_KNOWN_BACKENDS = ("anthropic", "claude_code")


@dataclass
class McpServerConfig:
  """One MCP server entry from a profile's mcp_servers list.

  Built by _build_profile from the raw JSON dict so downstream consumers
  (ChatWindow, McpServerConnection) can use attribute access rather than
  hunting through nested dicts. Optional fields default to safe values so
  shorter JSON declarations still work."""
  name: str
  command: str = None
  args: list = field(default_factory=list)
  env: dict = field(default_factory=dict)
  tool_policy: dict = field(default_factory=dict)
  auto_start: bool = True


@dataclass
class ClaudeProfileOptions:
  """Backend-specific knobs for the claude_code backend.

  Only consulted when Profile.backend == 'claude_code'. Exposed with
  safe defaults on every Profile so call sites can read attributes
  unconditionally without `hasattr`/`getattr` guards."""
  cli_path: str = None              # Override PATH lookup of `claude`
  sdk_options: dict = field(default_factory=dict)  # ClaudeAgentOptions kwargs


@dataclass
class Profile:
  name: str
  model: str
  description: str = ""
  max_tokens: int = 8192
  thinking_enabled: bool = False
  thinking_budget: int = 0
  vision_input: bool = True
  system_prompt: str = ""                              # resolved (file content inlined)
  skills_additional: list = field(default_factory=list)
  skills_disabled: set = field(default_factory=set)
  mcp_servers: list = field(default_factory=list)      # list[McpServerConfig]; expanded
  tool_policy_default: str = "ask"
  subagents_enabled: bool = True
  subagents_max_depth: int = 1
  subagents_default_max_turns: int = 25
  subagents_default_model: str = "claude-opus-4-7"
  subagents_default_profile: str = None
  backend: str = "claude_code"                         # "anthropic" | "claude_code"
  claude: ClaudeProfileOptions = field(default_factory=ClaudeProfileOptions)
  source_path: Path = None                             # which file this came from


class ProfileLoader:
  """Resolve and load profiles with cascading override and inheritance."""

  def __init__(self, builtin_dir, user_dir=None, project_dir=None, log=None):
    # Search order: project > user > built-in (Section 13.4).
    self.search_dirs = [d for d in (project_dir, user_dir, builtin_dir) if d]
    self.log = log if log is not None else sys.stdout

  def load(self, name):
    """Resolve, parse, validate, expand. Raises Sorry on critical errors."""
    return self._load_with_chain(name, _seen=())

  def _load_with_chain(self, name, _seen):
    if name in _seen:
      raise Sorry("Profile inheritance cycle: %s -> %s"
                  % (" -> ".join(_seen), name))
    _seen = _seen + (name,)

    path = self._resolve_path(name)
    if path is None:
      raise Sorry("Profile not found: %s" % name)
    try:
      with open(path) as fh:
        data = json.load(fh)
    except Exception as e:
      raise Sorry("Profile %s parse error: %s" % (path, e))

    # Inheritance: parent first, child overrides.
    based_on = data.get("based_on")
    if based_on:
      parent = self._load_with_chain(based_on, _seen)
      merged = _profile_to_dict(parent)
      merged.update({k: v for k, v in data.items()
                     if k != "based_on" and v is not None})
      data = merged

    # Warn on unknown keys (forward compatibility).
    for k in data.keys():
      if k not in _KNOWN_KEYS:
        print("profile %s: unknown key '%s' (ignored)"
              % (path, k), file=self.log)

    # Required fields.
    if "name" not in data:
      raise Sorry("Profile %s: missing required field 'name'" % path)
    if "model" not in data:
      raise Sorry("Profile %s: missing required field 'model'" % path)

    # Validate backend at load time so typos fail before the chat window
    # opens, not after the user has been waiting for a streaming response.
    backend = data.get("backend", "claude_code")
    if backend not in _KNOWN_BACKENDS:
      raise Sorry(
        "Profile %s: unknown backend %r (expected one of: %s)"
        % (path, backend, ", ".join(repr(b) for b in _KNOWN_BACKENDS)))

    return _build_profile(data, source_path=path)

  def _resolve_path(self, name):
    for d in self.search_dirs:
      candidate = Path(d) / ("%s.json" % name)
      if candidate.exists():
        return candidate
    return None


# ---- helpers ---------------------------------------------------------------

def _build_profile(data, source_path):
  thinking = data.get("thinking", {}) or {}
  skills = data.get("skills", {}) or {}
  subagents = data.get("subagents", {}) or {}
  # The claude subtree is parsed on every load (cheap), but only the
  # claude_code backend will consult it; for anthropic profiles it just
  # carries default values so attribute access is always safe.
  claude_data = data.get("claude", {}) or {}
  claude_opts = ClaudeProfileOptions(
    cli_path=claude_data.get("cli_path"),
    sdk_options=dict(claude_data.get("sdk_options") or {}))

  # System prompt: inline xor file (file content inlined here if present).
  system_prompt = data.get("system_prompt", "") or ""
  system_prompt_file = data.get("system_prompt_file")
  if system_prompt and system_prompt_file:
    raise Sorry("Profile %s: system_prompt and system_prompt_file are mutually exclusive" % source_path)
  if system_prompt_file:
    # Per Section 13.3: system_prompt_file goes through the same ${VAR}
    # expansion as mcp_servers fields. E.g., "${PROFILE_DIR}/prompt.md".
    expanded = _expand_str(system_prompt_file, source_path)
    file_path = (Path(source_path).parent / expanded).resolve()
    if not file_path.exists():
      raise Sorry("Profile %s: system_prompt_file not found: %s"
                  % (source_path, file_path))
    system_prompt = file_path.read_text()

  return Profile(
    name=data["name"],
    model=data["model"],
    description=data.get("description", ""),
    max_tokens=int(data.get("max_tokens", 8192)),
    thinking_enabled=bool(thinking.get("enabled", False)),
    thinking_budget=int(thinking.get("budget_tokens", 0)),
    vision_input=bool(data.get("vision_input", True)),
    system_prompt=system_prompt,
    skills_additional=list(skills.get("additional", [])),
    skills_disabled=set(skills.get("disabled", [])),
    mcp_servers=_expand_mcp_servers(
      data.get("mcp_servers", []), source_path=source_path),
    tool_policy_default=data.get("tool_policy_default", "ask"),
    subagents_enabled=bool(subagents.get("enabled", True)),
    subagents_max_depth=int(subagents.get("max_depth", 1)),
    subagents_default_max_turns=int(subagents.get("default_max_turns", 25)),
    subagents_default_model=subagents.get(
      "default_model", "claude-opus-4-7"),
    subagents_default_profile=subagents.get("default_profile"),
    backend=data.get("backend", "claude_code"),
    claude=claude_opts,
    source_path=source_path,
  )


def _profile_to_dict(p):
  """Reverse of _build_profile for inheritance merging. Only emits the
  fields we expose; helpers / resolved file content stays inlined."""
  return {
    "name": p.name,
    "description": p.description,
    "model": p.model,
    "max_tokens": p.max_tokens,
    "thinking": {
      "enabled": p.thinking_enabled,
      "budget_tokens": p.thinking_budget,
    },
    "vision_input": p.vision_input,
    "system_prompt": p.system_prompt,
    "skills": {
      "additional": list(p.skills_additional),
      "disabled": sorted(p.skills_disabled),
    },
    "mcp_servers": list(p.mcp_servers),
    "tool_policy_default": p.tool_policy_default,
    "subagents": {
      "enabled": p.subagents_enabled,
      "max_depth": p.subagents_max_depth,
      "default_max_turns": p.subagents_default_max_turns,
      "default_profile": p.subagents_default_profile,
      "default_model": p.subagents_default_model,
    },
    "backend": p.backend,
    "claude": {
      "cli_path": p.claude.cli_path,
      "sdk_options": dict(p.claude.sdk_options),
    },
  }


_VAR_RE = re.compile(r"\$\{([^}]+)\}")


def _expand_str(value, source_path):
  def _sub(m):
    token = m.group(1)
    if token.startswith("env:"):
      return os.environ.get(token[4:], "")
    if token == "PROFILE_DIR":
      return str(Path(source_path).parent.resolve())
    if token == "PROJECT_DIR":
      return os.environ.get("PHENIX_PROJECT_DIR", "")
    if token == "CHAT_ROOT":
      return os.environ.get("PHENIX_CHAT_HOME", "")
    return m.group(0)                  # leave unknown tokens as-is
  return _VAR_RE.sub(_sub, value)


def _expand_in_obj(obj, source_path):
  if isinstance(obj, str):
    return _expand_str(obj, source_path)
  if isinstance(obj, list):
    return [_expand_in_obj(v, source_path) for v in obj]
  if isinstance(obj, dict):
    return {k: _expand_in_obj(v, source_path) for k, v in obj.items()}
  return obj


def _expand_mcp_servers(servers, source_path):
  out = []
  for raw in servers:
    expanded = _expand_in_obj(raw, source_path)
    if isinstance(expanded, McpServerConfig):
      out.append(expanded)
      continue
    if not isinstance(expanded, dict):
      raise Sorry(
        "Profile %s: mcp_servers entries must be objects, got: %r"
        % (source_path, expanded))
    name = expanded.get("name")
    if not name:
      raise Sorry(
        "Profile %s: every mcp_servers entry needs a 'name'" % source_path)
    out.append(McpServerConfig(
      name=name,
      command=expanded.get("command"),
      args=list(expanded.get("args") or []),
      env=dict(expanded.get("env") or {}),
      tool_policy=dict(expanded.get("tool_policy") or {}),
      auto_start=bool(expanded.get("auto_start", True)),
    ))
  return out
