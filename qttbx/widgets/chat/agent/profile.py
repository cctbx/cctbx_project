"""Profile loader and Profile dataclass.

JSON files; cascading override (project > user > built-in); one level
of ``based_on`` inheritance with cycle detection; ``${VAR}`` and
``${env:NAME}`` expansion in ``mcp_servers`` strings.
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
  "server_tools",
  "backend", "claude", "portkey",
}

# Backends supported by the agent factory. Validated at load time so
# typos surface before the chat window opens.
_KNOWN_BACKENDS = ("anthropic", "claude_code", "openai", "portkey", "google")

# User-facing assistant name per backend, for the chat UI (response label,
# input placeholder, question card). Portkey is a gateway whose underlying
# model varies, so it -- and anything unmapped -- gets the generic "Assistant".
_BACKEND_DISPLAY_NAMES = {
  "claude_code": "Claude",
  "anthropic": "Claude",
  "openai": "GPT",
  "google": "Gemini",
  "portkey": "Assistant",
}


def backend_display_name(backend):
  """Return the user-facing assistant name for a backend id.

  Parameters
  ----------
  backend : str or None
      A backend id (e.g. ``"openai"``). Unknown / empty / ``None`` map to
      ``"Assistant"``.

  Returns
  -------
  str
      The display name (``"Claude"``, ``"GPT"``, ``"Gemini"``, or
      ``"Assistant"``).
  """
  return _BACKEND_DISPLAY_NAMES.get(backend or "", "Assistant")

# Environment keys a profile-supplied MCP-server ``env`` block must not be
# able to set. A profile can be project-scoped (and therefore untrusted),
# so it must not control which binary or libraries load into the spawned
# subprocess (code execution) or override Phenix's own scoping. Two layers:
#   * the exact keys below -- executable lookup (PATH), shell startup-code
#     injection, and the per-interpreter module/loader hijacks for
#     python / node / perl / ruby; and
#   * the prefix families in _BLOCKED_SERVER_ENV_PREFIXES.
# Anything else (per-server API tokens, PYTHONUNBUFFERED, ...) passes
# through untouched -- this is a blocklist, not an allowlist, so that
# legitimate per-server configuration still works.
_BLOCKED_SERVER_ENV = frozenset([
  "PATH",                                        # which binary / helpers resolve
  "BASH_ENV", "ENV", "IFS",                      # shell startup-code injection
  "PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP",   # python module / loader hijack
  "NODE_OPTIONS", "NODE_PATH",                   # node --require / module hijack
  "PERL5LIB", "PERL5OPT",                        # perl module path / -M code
  "RUBYLIB", "RUBYOPT",                          # ruby load path / -r code
])

# Prefix families stripped wholesale -- no legitimate per-server env sets
# these: every dynamic-loader control (LD_PRELOAD, LD_AUDIT, LD_LIBRARY_PATH,
# DYLD_INSERT_LIBRARIES, DYLD_*_PATH, ...), exported shell functions
# (BASH_FUNC_*), and Phenix's own namespace (PHENIX_*, which includes
# PHENIX_PROJECT_DIR -- the project sandbox -- and PHENIX_TRUST_OTHER_ENV,
# which would otherwise switch off the dispatcher's inherited-env scrub).
_BLOCKED_SERVER_ENV_PREFIXES = ("LD_", "DYLD_", "BASH_FUNC_", "PHENIX_")


def sanitize_server_env(env):
  """Strip security-sensitive keys from an MCP-server ``env`` dict.

  Removes keys that would let an untrusted profile control which binary or
  libraries load into the spawned MCP-server subprocess (code execution)
  or override Phenix's own scoping. See ``_BLOCKED_SERVER_ENV`` and
  ``_BLOCKED_SERVER_ENV_PREFIXES`` for the exact set.

  Parameters
  ----------
  env : dict or None
      Raw per-server environment from a profile.

  Returns
  -------
  dict
      A new dict with the blocked keys removed; all other entries
      preserved.
  """
  if not env:
    return {}
  out = {}
  for k, v in env.items():
    if k in _BLOCKED_SERVER_ENV:
      continue
    if any(k.startswith(p) for p in _BLOCKED_SERVER_ENV_PREFIXES):
      continue
    out[k] = v
  return out


@dataclass
class McpServerConfig:
  """One MCP server entry from a profile's ``mcp_servers`` list.

  Built by ``_build_profile`` from the raw JSON dict so downstream
  consumers (``ChatWindow``, ``McpServerConnection``) can use attribute
  access rather than hunting through nested dicts. Optional fields default
  to safe values so shorter JSON declarations still work.
  """
  name: str
  command: str = None
  args: list = field(default_factory=list)
  env: dict = field(default_factory=dict)
  tool_policy: dict = field(default_factory=dict)
  auto_start: bool = True
  # Phenix-aware servers (default) receive PHENIX_PROJECT_DIR /
  # PHENIX_CHAT_HOME in their subprocess env; foreign servers (e.g. the
  # Coot bridge) set this False so no Phenix scoping leaks into them.
  inject_phenix_env: bool = True


@dataclass
class ClaudeProfileOptions:
  """Backend-specific knobs for the ``claude_code`` backend.

  Only consulted when ``Profile.backend == 'claude_code'``. Exposed with
  safe defaults on every ``Profile`` so call sites can read attributes
  unconditionally without ``hasattr``/``getattr`` guards.
  """
  cli_path: str = None              # Override PATH lookup of `claude`
  sdk_options: dict = field(default_factory=dict)  # ClaudeAgentOptions kwargs


@dataclass
class Profile:
  """A fully resolved agent profile.

  Produced by ``ProfileLoader.load`` after cascading override, one level
  of ``based_on`` inheritance, and ``${VAR}`` expansion. File-backed
  fields (e.g. ``system_prompt``) are inlined as resolved content.
  """
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
  # Provider-side tools the user opts into (canonical names, e.g.
  # "web_search", "web_fetch", "code_execution"). Empty by default so an
  # existing profile that omits the key gets pure-passthrough behavior --
  # the backend declares nothing extra. Only the anthropic backend
  # currently consults this; other backends carry the default empty list.
  server_tools: list = field(default_factory=list)
  subagents_enabled: bool = True
  subagents_max_depth: int = 1
  subagents_default_max_turns: int = 25
  subagents_default_model: str = "claude-opus-4-7"
  subagents_default_profile: str = None
  backend: str = "claude_code"  # anthropic|claude_code|openai|portkey|google
  claude: ClaudeProfileOptions = field(default_factory=ClaudeProfileOptions)
  # Backend-specific knobs for the "portkey" backend; only consulted when
  # backend == "portkey". Default None on every Profile so call sites can
  # read them unconditionally (the launcher falls back to these).
  portkey_virtual_key: str = None
  portkey_config: str = None
  source_path: Path = None                             # which file this came from


class ProfileLoader:
  """Resolve and load profiles with cascading override and inheritance.

  Directories are searched in the order project > user > built-in, so a
  same-named profile in an earlier directory overrides a later one.

  Parameters
  ----------
  builtin_dir : str or pathlib.Path
      Directory of profiles shipped with the application.
  user_dir : str or pathlib.Path, optional
      User-level profile directory; overrides ``builtin_dir``.
  project_dir : str or pathlib.Path, optional
      Project-level profile directory; overrides ``user_dir``.
  log : file-like, optional
      Destination for warnings. Defaults to ``sys.stdout``.
  """

  def __init__(self, builtin_dir, user_dir=None, project_dir=None, log=None):
    # Search order: project > user > built-in.
    self.search_dirs = [d for d in (project_dir, user_dir, builtin_dir) if d]
    self.log = log if log is not None else sys.stdout

  def load(self, name):
    """Resolve, parse, validate, and expand a profile by name.

    Parameters
    ----------
    name : str
        Profile name (without the ``.json`` extension).

    Returns
    -------
    Profile
        The fully resolved profile.

    Raises
    ------
    libtbx.utils.Sorry
        On critical errors: profile not found, parse error, inheritance
        cycle, missing required field, or unknown backend.
    """
    return self._load_with_chain(name, _seen=())

  def load_file(self, path):
    """Resolve, parse, validate, and expand a profile from an explicit file.

    Like :meth:`load`, but the root profile is read from ``path`` directly
    instead of by name. The file's own directory is searched first for its
    ``based_on`` parents (so it can reference sibling profiles, mirroring how
    ``${PROFILE_DIR}`` lets ``system_prompt_file`` reference siblings); names
    not found there fall through to the loader's configured search dirs, so an
    external file can still inherit from a shipped profile such as
    ``phenix_base``. The extra directory is scoped to this call -- the loader
    is not mutated, so a later :meth:`load` on the same instance is unaffected.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to a profile ``.json`` file.

    Returns
    -------
    Profile
        The fully resolved profile.

    Raises
    ------
    libtbx.utils.Sorry
        When ``path`` is not a file, or on any error :meth:`load` can raise
        (parse error, inheritance cycle, missing field, unknown backend).
    """
    p = Path(path)
    if not p.is_file():
      raise Sorry("Profile file not found: %s" % path)
    p = p.resolve()
    # Search the file's own directory first (so based_on can reach siblings),
    # threaded through resolution rather than mutating self.search_dirs.
    return self._load_path(p, _seen=(), extra_dirs=(p.parent,))

  def _load_with_chain(self, name, _seen, extra_dirs=()):
    """Resolve ``name`` to a path, then load it via :meth:`_load_path`.

    ``_seen`` carries the chain of already-visited names so an inheritance
    cycle raises rather than recursing forever. ``extra_dirs`` are searched
    ahead of ``self.search_dirs`` (threaded from :meth:`load_file`).
    """
    if name in _seen:
      raise Sorry("Profile inheritance cycle: %s -> %s"
                  % (" -> ".join(_seen), name))
    path = self._resolve_path(name, extra_dirs)
    if path is None:
      raise Sorry("Profile not found: %s" % name)
    return self._load_path(path, _seen + (name,), extra_dirs)

  def _load_path(self, path, _seen, extra_dirs=()):
    """Parse the profile at ``path``, resolve ``based_on``, validate, build.

    The shared body of :meth:`_load_with_chain` (name-resolved) and
    :meth:`load_file` (path-given). ``_seen`` is the inheritance chain so far;
    ``extra_dirs`` are forwarded to ``based_on`` resolution.
    """
    try:
      with open(path) as fh:
        data = json.load(fh)
    except Exception as e:
      raise Sorry("Profile %s parse error: %s" % (path, e))

    # Inheritance: parent first, child overrides.
    based_on = data.get("based_on")
    if based_on:
      parent = self._load_with_chain(based_on, _seen, extra_dirs)
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

  def _resolve_path(self, name, extra_dirs=()):
    """Return the first ``<dir>/<name>.json`` that exists, else ``None``.

    ``extra_dirs`` are searched before ``self.search_dirs`` (used by
    :meth:`load_file` so an external file's siblings resolve first), which
    keeps the loader itself stateless.
    """
    for d in tuple(extra_dirs) + tuple(self.search_dirs):
      candidate = Path(d) / ("%s.json" % name)
      if candidate.exists():
        return candidate
    return None


# ---- helpers ---------------------------------------------------------------

def _server_tools_list(value, source_path):
  """Validate a profile's ``server_tools`` field as a list of tool names.

  A bare string (e.g. ``"web_search"``) is a common mistake that plain
  ``list()`` would silently split into single characters -- yielding ~10 bogus
  tool names that each warn as unknown and enable nothing. Raise a clear Sorry
  instead.
  """
  if value is None:
    return []
  if not isinstance(value, list):
    raise Sorry("Profile %s: server_tools must be a list of tool names "
                "(e.g. [\"web_search\"]), not %r" % (source_path, value))
  return list(value)


def _build_profile(data, source_path):
  """Construct a ``Profile`` from a merged, validated profile dict.

  Flattens nested subtrees (``thinking``, ``skills``, ``subagents``,
  ``claude``) into ``Profile`` fields and inlines ``system_prompt_file``
  content with ``${VAR}`` expansion.

  Parameters
  ----------
  data : dict
      Merged profile data (parent overlaid by child for inheritance).
  source_path : str or pathlib.Path
      File the data came from; used for ``${PROFILE_DIR}`` expansion and
      error messages.

  Returns
  -------
  Profile
      The constructed profile.

  Raises
  ------
  libtbx.utils.Sorry
      When ``system_prompt`` and ``system_prompt_file`` are both set, or
      when ``system_prompt_file`` cannot be found.
  """
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

  # The portkey subtree is parsed on every load (cheap), but only the
  # portkey backend will consult it; for other backends it just yields
  # None defaults so attribute access on the Profile is always safe.
  portkey_data = data.get("portkey", {}) or {}

  # System prompt: inline xor file (file content inlined here if present).
  system_prompt = data.get("system_prompt", "") or ""
  system_prompt_file = data.get("system_prompt_file")
  if system_prompt and system_prompt_file:
    raise Sorry("Profile %s: system_prompt and system_prompt_file are mutually exclusive" % source_path)
  if system_prompt_file:
    # system_prompt_file goes through the same ${VAR} expansion as
    # mcp_servers fields. E.g., "${PROFILE_DIR}/prompt.md".
    expanded = _expand_str(system_prompt_file, source_path)
    base = Path(source_path).parent.resolve()
    file_path = (base / expanded).resolve()
    # Containment: a profile may be project-supplied (untrusted), so an
    # absolute path or ".." escape in system_prompt_file could inline an
    # arbitrary file (e.g. ~/.ssh/id_rsa) into the system prompt. Require
    # the resolved path to stay within the profile's own directory.
    if base not in file_path.parents:
      raise Sorry(
        "Profile %s: system_prompt_file escapes the profile directory: %s"
        % (source_path, system_prompt_file))
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
    server_tools=_server_tools_list(data.get("server_tools"), source_path),
    subagents_enabled=bool(subagents.get("enabled", True)),
    subagents_max_depth=int(subagents.get("max_depth", 1)),
    subagents_default_max_turns=int(subagents.get("default_max_turns", 25)),
    subagents_default_model=subagents.get(
      "default_model", "claude-opus-4-7"),
    subagents_default_profile=subagents.get("default_profile"),
    backend=data.get("backend", "claude_code"),
    claude=claude_opts,
    portkey_virtual_key=portkey_data.get("virtual_key"),
    portkey_config=portkey_data.get("config"),
    source_path=source_path,
  )


def _profile_to_dict(p):
  """Serialize a ``Profile`` back to a dict for inheritance merging.

  The inverse of ``_build_profile``. Only emits the fields we expose;
  helper state and resolved file content stay inlined.

  Parameters
  ----------
  p : Profile
      The parent profile being merged into a child.

  Returns
  -------
  dict
      Nested profile data suitable for overlaying with child values.
  """
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
    "server_tools": list(p.server_tools),
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
    "portkey": {
      "virtual_key": p.portkey_virtual_key,
      "config": p.portkey_config,
    },
  }


_VAR_RE = re.compile(r"\$\{([^}]+)\}")


def _expand_str(value, source_path):
  """Expand ``${VAR}`` tokens in a single string.

  Recognized tokens: ``${env:NAME}`` (environment variable, empty if
  unset), ``${PROFILE_DIR}`` (the profile file's directory),
  ``${PROJECT_DIR}`` (``PHENIX_PROJECT_DIR``), and ``${CHAT_ROOT}``
  (``PHENIX_CHAT_HOME``). Unknown tokens are left untouched.

  Parameters
  ----------
  value : str
      The string to expand.
  source_path : str or pathlib.Path
      Profile file path, used to resolve ``${PROFILE_DIR}``.

  Returns
  -------
  str
      The expanded string.
  """
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
  """Recursively expand ``${VAR}`` tokens in all strings within ``obj``.

  Walks strings, lists, and dicts; other types are returned unchanged.
  """
  if isinstance(obj, str):
    return _expand_str(obj, source_path)
  if isinstance(obj, list):
    return [_expand_in_obj(v, source_path) for v in obj]
  if isinstance(obj, dict):
    return {k: _expand_in_obj(v, source_path) for k, v in obj.items()}
  return obj


def _expand_mcp_servers(servers, source_path):
  """Expand ``${VAR}`` tokens and normalize entries to ``McpServerConfig``.

  Entries already of type ``McpServerConfig`` (e.g. carried through
  inheritance) pass through unchanged.

  Parameters
  ----------
  servers : list
      Raw ``mcp_servers`` entries (dicts or ``McpServerConfig``).
  source_path : str or pathlib.Path
      Profile file path, used for expansion and error messages.

  Returns
  -------
  list of McpServerConfig
      The normalized server configs.

  Raises
  ------
  libtbx.utils.Sorry
      When an entry is not an object or is missing a ``name``.
  """
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
      env=sanitize_server_env(expanded.get("env")),
      tool_policy=dict(expanded.get("tool_policy") or {}),
      auto_start=bool(expanded.get("auto_start", True)),
      inject_phenix_env=bool(expanded.get("inject_phenix_env", True)),
    ))
  return out
