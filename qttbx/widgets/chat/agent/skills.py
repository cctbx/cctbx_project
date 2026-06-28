"""Skill loader and ``Skill`` dataclass.

Skills are directories with a ``SKILL.md`` (YAML frontmatter + markdown
body) plus optional helper files. Built-in skills auto-load; user /
project skills are opt-in via ``profile.skills.additional``. Built-in
wins on name collision. Path safety in ``read_file`` uses realpath
containment so ``..``, absolute paths, and symlink escapes all reject
uniformly.
"""

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

from libtbx.utils import Sorry


@dataclass
class Skill:
  """A loaded skill: metadata plus its ``SKILL.md`` body.

  Parameters
  ----------
  name : str
      Skill name (from frontmatter); unique within a session.
  path : pathlib.Path
      Absolute path to the skill directory.
  description : str
      Short description surfaced in the system prompt.
  mode : str, optional
      ``"always"`` (body inlined in the prompt) or ``"on_demand"``
      (body fetched via the ``load_skill`` tool). Defaults to
      ``"always"``.
  body : str, optional
      The ``SKILL.md`` body with frontmatter stripped.
  requires : list, optional
      MCP server names the skill depends on.
  """
  name: str
  path: Path                                    # absolute path to skill dir
  description: str
  mode: str = "always"                          # "always" | "on_demand"
  body: str = ""                                # SKILL.md body (frontmatter stripped)
  requires: list = field(default_factory=list)  # MCP server names


class SkillLoader:
  """Discover, load, and surface skills for an agent session.

  Built-in skills auto-load; user / project skills are opt-in. Built-in
  wins on name collision. Also provides the skill-source tools that let
  the agent read and list files inside a loaded skill's directory.

  Parameters
  ----------
  builtin_path : str or pathlib.Path
      Directory of auto-loading built-in skills.
  project_path : str or pathlib.Path, optional
      Project-level skill search path.
  user_path : str or pathlib.Path, optional
      User-level skill search path.
  log : file-like, optional
      Stream for diagnostic messages. Defaults to ``sys.stdout``.
  """

  def __init__(self, builtin_path, project_path=None, user_path=None, log=None):
    self.builtin_path = Path(builtin_path)
    self.search_paths = [Path(p) for p in (project_path, user_path) if p]
    self.log = log if log is not None else sys.stdout

  # ---- discovery / loading -------------------------------------------------

  def load_default(self, additional=(), disabled=frozenset(),
                   mcp_servers=None):
    """Load built-in skills (minus ``disabled``), then ``additional``
    skills resolved against project/user paths or absolute paths.
    Built-in wins on name collision.

    Parameters
    ----------
    additional : iterable of str or dict, optional
        Extra skill references. A string is resolved against the
        loader's search paths; a dict with a ``path`` key is taken as
        an absolute / project-relative path.
    disabled : frozenset of str, optional
        Names of built-in skills to skip.
    mcp_servers : iterable of str, optional
        Names of MCP servers known to be available for this session
        (typically derived from ``profile.mcp_servers``, or ``[]`` when
        ``--no-mcp`` is set). When provided, skills whose ``requires``
        frontmatter lists a server outside this set are filtered out
        -- loading them would only mislead the agent into calling
        tools that don't exist. ``None`` (default) skips the check
        entirely, which suits callers that don't know which servers
        will be wired up.

    Returns
    -------
    list of Skill
        Loaded skills in the order they will appear in the system
        prompt.
    """
    if mcp_servers is None:
      available = None
    else:
      available = frozenset(mcp_servers)
    skills = self._load_builtins(disabled, available)
    seen = {s.name for s in skills}
    for ref in additional:
      skill = self._resolve_additional(ref, available)
      if skill is None:
        continue
      if skill.name in seen:
        print("skill '%s': built-in shadows additional"
              % skill.name, file=self.log)
        continue
      skills.append(skill)
      seen.add(skill.name)
    return skills

  def _load_builtins(self, disabled, available_servers):
    skills = []
    if not self.builtin_path.exists():
      return skills
    for entry in sorted(self.builtin_path.iterdir()):
      if not entry.is_dir():
        continue
      if entry.name in disabled:
        continue
      skill = self._load_checked(entry, available_servers)
      if skill is not None:
        skills.append(skill)
    return skills

  def _resolve_additional(self, ref, available_servers):
    # ref is a string (skill name) or dict with "path".
    path = None
    if hasattr(ref, "get"):
      raw = ref.get("path")
      if raw:
        path = Path(raw)
    else:
      for search in self.search_paths:
        candidate = search / ref
        if candidate.is_dir():
          path = candidate
          break
    if path is None:
      print("skill '%s' not found in any search path" % ref, file=self.log)
      return None
    return self._load_checked(path, available_servers)

  def _load_checked(self, path, available_servers):
    """Load the skill at ``path`` and apply the ``requires`` filter.

    The shared body of :meth:`_load_builtins` and :meth:`_resolve_additional`:
    load the ``SKILL.md`` (skipping with a log on a malformed / unreadable
    one), then drop it (with a log) when its ``requires`` lists an MCP server
    outside ``available_servers``.

    Parameters
    ----------
    path : pathlib.Path
        The skill directory to load.
    available_servers : frozenset of str or None
        Names of available MCP servers, or ``None`` to skip the check.

    Returns
    -------
    Skill or None
        The loaded skill, or ``None`` when it should be skipped.
    """
    try:
      skill = self._load_one(path)
    except (Sorry, OSError, ValueError) as e:
      # One malformed / unreadable skill (e.g. an undecodable SKILL.md raising
      # UnicodeDecodeError, a ValueError subclass) is skipped so it does not
      # abort the whole batch and drop every other skill.
      print("skill %s: skipping (%s)" % (path, e), file=self.log)
      return None
    if not self._requires_satisfied(skill, available_servers):
      missing = sorted(set(skill.requires) - available_servers)
      print("skill '%s': skipping (requires MCP servers not available: %s)"
            % (skill.name, missing), file=self.log)
      return None
    return skill

  def _requires_satisfied(self, skill, available_servers):
    """Report whether the skill's required MCP servers are available.

    Parameters
    ----------
    skill : Skill
        The skill whose ``requires`` declaration is checked.
    available_servers : frozenset of str or None
        Names of available MCP servers. ``None`` means the caller opted
        out of the check (typically because it does not know which
        servers will be available, e.g., test setup), in which case
        every skill loads regardless of its ``requires`` declaration.

    Returns
    -------
    bool
        True when the skill can usefully run.
    """
    if available_servers is None:
      return True
    if not skill.requires:
      return True
    return all(r in available_servers for r in skill.requires)

  def _load_one(self, path):
    skill_md = path / "SKILL.md"
    if not skill_md.exists():
      raise Sorry("SKILL.md missing in %s" % path)
    frontmatter, body = _parse_skill_md(skill_md)
    if "name" not in frontmatter:
      raise Sorry("SKILL.md missing required field 'name'")
    if "description" not in frontmatter:
      raise Sorry("SKILL.md missing required field 'description'")
    # `requires` may be a YAML list, a scalar, or absent. First normalize any
    # non-list (a bare scalar `requires: phenix`) to a one-element list -- list()
    # on a string would otherwise explode it into ['p','h','e','n','i','x'].
    # Then coerce EVERY element to str: requires entries are MCP server NAMES
    # (strings), and stringifying makes the list uniformly hashable. A YAML
    # MAPPING (`requires: {a: b}`) or a list element that is itself a dict/list
    # parses to an UNHASHABLE value; _requires_satisfied's `r in
    # available_servers` (frozenset membership) would then raise TypeError --
    # which is OUTSIDE _load_checked's (Sorry, OSError, ValueError) guard and
    # would abort the WHOLE skill batch. str() also folds non-string scalars
    # (`requires: 1`, `requires: true`) in uniformly (mirroring the str()
    # coercion of name/description/mode below); a stringified non-name simply
    # matches no server, so the skill is cleanly filtered rather than crashing.
    requires = frontmatter.get("requires") or []
    if not isinstance(requires, list):
      requires = [requires]
    requires = [str(r) for r in requires]
    return Skill(
      name=str(frontmatter["name"]),
      path=path.resolve(),
      description=str(frontmatter["description"]),
      mode=str(frontmatter.get("mode", "always")),
      body=body,
      requires=requires,
    )

  # ---- system prompt assembly ----------------------------------------------

  def assemble_system_prompt(self, base_prompt, skills):
    """Build the session's system prompt from the loaded skills.

    Re-callable so an agent can rebuild the prompt mid-session if the
    skill set changes. ``always``-mode skills inline their full body;
    ``on_demand`` skills emit name + description only (the body is
    fetched via the ``load_skill`` tool).

    Parameters
    ----------
    base_prompt : str
        Profile-level system prompt that precedes the skill sections.
    skills : list of Skill
        Skills to surface in the prompt.

    Returns
    -------
    str
        The assembled system prompt.
    """
    sections = [base_prompt, ""]
    if skills:
      sections.append("# Skills available\n")
      sections.append(
        "The following skills are loaded for this session. Use them when "
        "their description matches the user's request.\n")
      for s in skills:
        sections.append("## Skill: %s" % s.name)
        sections.append(s.description)
        if s.mode == "always":
          sections.append("")
          sections.append(s.body)
        sections.append("")
      sections.append("# Skill helper files\n")
      sections.append(
        "Use `read_skill_file(skill_name, relative_path)` to read any "
        "file inside a skill, `list_skill_files(skill_name)` to enumerate, "
        "and `load_skill(skill_name)` to fetch the body of an on-demand skill.\n")
    return "\n".join(sections)

  # ---- tools (registered with ToolRegistry) --------------------------------

  def tools(self, skills):
    """Return the skill-loader tools to register with a ``ToolRegistry``.

    ``read_skill_file`` and ``list_skill_files`` are always emitted
    when any skill is loaded; ``load_skill`` is only emitted when at
    least one skill is ``mode=on_demand`` (otherwise every skill's
    body is already inline in the system prompt and ``load_skill``
    has no purpose). All three are skill-source ``risk='read'``.

    Parameters
    ----------
    skills : list of Skill
        Skills the agent has access to.

    Returns
    -------
    list of (ToolSpec, callable)
        Pairs ready to feed into ``ToolRegistry.register_skill_tool``.
    """
    from qttbx.widgets.chat.agent.base import ToolSpec
    if not skills:
      return []
    by_name = {s.name: s for s in skills}
    skill_names = sorted(by_name)
    enum_clause = ", ".join("'%s'" % n for n in skill_names)

    def read_file_handler(name, input):
      skill = by_name.get(input.get("skill_name", ""))
      if skill is None:
        raise Sorry("Unknown skill: %s" % input.get("skill_name"))
      data = self.read_file(skill, input.get("relative_path", ""))
      try:
        return data.decode("utf-8")
      except UnicodeDecodeError:
        import base64
        return {"binary_base64": base64.b64encode(data).decode("ascii")}

    def list_files_handler(name, input):
      skill = by_name.get(input.get("skill_name", ""))
      if skill is None:
        raise Sorry("Unknown skill: %s" % input.get("skill_name"))
      return self.list_files(skill)

    def load_skill_handler(name, input):
      skill = by_name.get(input.get("skill_name", ""))
      if skill is None:
        raise Sorry("Unknown skill: %s" % input.get("skill_name"))
      return self.load_skill_body(skill)

    out = [
      (ToolSpec(
        name="read_skill_file",
        description=(
          "Read a file from inside a loaded skill's directory (text is "
          "returned as a string; binary is returned base64-encoded under "
          "binary_base64). Available skills: " + enum_clause + "."),
        input_schema={
          "type": "object",
          "properties": {
            "skill_name": {"type": "string", "enum": skill_names},
            "relative_path": {"type": "string"},
          },
          "required": ["skill_name", "relative_path"],
        }),
       read_file_handler),
      (ToolSpec(
        name="list_skill_files",
        description=(
          "List every file (relative path) inside a loaded skill's "
          "directory. Available skills: " + enum_clause + "."),
        input_schema={
          "type": "object",
          "properties": {
            "skill_name": {"type": "string", "enum": skill_names},
          },
          "required": ["skill_name"],
        }),
       list_files_handler),
    ]
    if any(s.mode == "on_demand" for s in skills):
      out.append((
        ToolSpec(
          name="load_skill",
          description=(
            "Fetch the full body of an on-demand skill. Always-mode skills "
            "are already inlined in the system prompt; on-demand skills "
            "need this call to surface their body. Available skills: "
            + enum_clause + "."),
          input_schema={
            "type": "object",
            "properties": {
              "skill_name": {"type": "string", "enum": skill_names},
            },
            "required": ["skill_name"],
          }),
        load_skill_handler))
    return out

  def read_file(self, skill, relative_path):
    safe = self._resolve_safe(skill, relative_path)
    if safe is None:
      raise Sorry("read_skill_file: path escapes skill directory")
    if not os.path.exists(safe):
      raise Sorry("read_skill_file: file not found: %s" % relative_path)
    with open(safe, "rb") as fh:
      return fh.read()

  def list_files(self, skill):
    out = []
    skill_root = str(skill.path)
    for dirpath, _dirnames, filenames in os.walk(skill.path):
      for fn in filenames:
        rel = os.path.relpath(os.path.join(dirpath, fn), skill_root)
        out.append(rel)
    return sorted(out)

  def load_skill_body(self, skill):
    return skill.body

  def _resolve_safe(self, skill, relative_path):
    """Resolve a skill-relative path, rejecting any escape from the dir.

    Rejects absolute paths, then resolves via realpath and verifies
    containment under ``realpath(skill.path)``. Catches ``..``,
    absolute, and symlink-escape uniformly.

    Parameters
    ----------
    skill : Skill
        The skill whose directory bounds the resolved path.
    relative_path : str
        Path supplied relative to the skill directory.

    Returns
    -------
    str or None
        The validated absolute path, or ``None`` if it escapes the
        skill directory.
    """
    if os.path.isabs(relative_path):
      return None
    skill_root = os.path.realpath(str(skill.path))
    candidate = os.path.realpath(os.path.join(skill_root, relative_path))
    # candidate must equal skill_root or live under it (with a separator).
    if candidate == skill_root:
      return candidate
    if candidate.startswith(skill_root + os.sep):
      return candidate
    return None


def _parse_skill_md(path):
  """Parse a ``SKILL.md`` into YAML frontmatter and markdown body.

  Parameters
  ----------
  path : pathlib.Path
      Path to the ``SKILL.md`` file.

  Returns
  -------
  tuple of (dict, str)
      The parsed frontmatter mapping and the markdown body.

  Raises
  ------
  libtbx.utils.Sorry
      If the frontmatter is missing, unterminated, not a mapping, or
      malformed, or if PyYAML is not installed.
  """
  text = path.read_text(encoding="utf-8")
  if not text.startswith("---"):
    raise Sorry("SKILL.md %s missing YAML frontmatter" % path)
  # Split after the first '---' on its own line, then find the next.
  lines = text.splitlines()
  end_idx = None
  for i, line in enumerate(lines[1:], start=1):
    if line.strip() == "---":
      end_idx = i
      break
  if end_idx is None:
    raise Sorry("SKILL.md %s frontmatter not terminated" % path)
  fm_text = "\n".join(lines[1:end_idx])
  body = "\n".join(lines[end_idx + 1:]).lstrip("\n")
  try:
    import yaml
  except ImportError:
    raise Sorry("PyYAML required to parse SKILL.md (install pyyaml)")
  try:
    fm = yaml.safe_load(fm_text) or {}
  except Exception as e:
    raise Sorry("SKILL.md %s frontmatter parse error: %s" % (path, e))
  if not isinstance(fm, dict):
    raise Sorry("SKILL.md %s frontmatter is not a mapping" % path)
  return fm, body
