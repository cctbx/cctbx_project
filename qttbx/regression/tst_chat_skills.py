import os
import shutil
import sys
import tempfile
from pathlib import Path

from libtbx.test_utils import Exception_expected
from libtbx.utils import format_cpu_times, Sorry, null_out

# SKILL.md frontmatter is parsed via PyYAML inside SkillLoader. Without
# it every _load_one() call raises Sorry and the skipping branch in
# _load_builtins() silently drops the skill, so the test would fail with
# an unhelpful 'len(skills) == 1' assertion. Match the PySide-skip
# pattern: cleanly exit before any test runs.
try:
  import yaml  # noqa: F401
except ImportError:
  print("PyYAML not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.chat.agent.skills import SkillLoader


def _make_skill(root, name, body="example body\n",
                description="example description",
                mode="always", requires=None):
  d = os.path.join(root, name)
  os.makedirs(d)
  requires_line = ""
  if requires is not None:
    requires_line = "requires: [%s]\n" % ", ".join(requires)
  skill_md = (
    "---\n"
    "name: %s\n"
    "description: %s\n"
    "mode: %s\n"
    "%s"
    "---\n\n%s" % (name, description, mode, requires_line, body))
  with open(os.path.join(d, "SKILL.md"), "w") as fh:
    fh.write(skill_md)
  return d


def exercise_load_builtin_only():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "refinement_workflows")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    assert len(skills) == 1
    assert skills[0].name == "refinement_workflows"
  finally:
    shutil.rmtree(tmp)


def exercise_builtin_wins_on_collision():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    user = os.path.join(tmp, "user")
    os.makedirs(builtin)
    os.makedirs(user)
    _make_skill(builtin, "same_name", body="from builtin\n")
    _make_skill(user, "same_name", body="from user\n")
    loader = SkillLoader(builtin_path=Path(builtin),
                         user_path=Path(user),
                         log=null_out())
    skills = loader.load_default(additional=["same_name"])
    assert len(skills) == 1
    assert "from builtin" in skills[0].body
  finally:
    shutil.rmtree(tmp)


def exercise_read_file_rejects_traversal():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "x")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skill = loader.load_default()[0]
    try:
      loader.read_file(skill, "../../etc/passwd")
    except Sorry:
      pass
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_read_file_rejects_absolute_path():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "x")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skill = loader.load_default()[0]
    try:
      loader.read_file(skill, "/etc/passwd")
    except Sorry:
      pass
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_read_file_rejects_symlink_escape():
  """Symlink inside skill dir pointing outside must be rejected. Different
  code path than '..' traversal (realpath check is what catches it)."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    outside = os.path.join(tmp, "outside_target.txt")
    os.makedirs(builtin)
    _make_skill(builtin, "x")
    with open(outside, "w") as fh:
      fh.write("secret\n")
    os.symlink(outside, os.path.join(builtin, "x", "leak"))
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skill = loader.load_default()[0]
    try:
      loader.read_file(skill, "leak")
    except Sorry:
      pass
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_read_file_succeeds_for_in_skill_file():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    skill_dir = _make_skill(builtin, "x")
    with open(os.path.join(skill_dir, "helper.md"), "w") as fh:
      fh.write("helper content")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skill = loader.load_default()[0]
    data = loader.read_file(skill, "helper.md")
    assert data == b"helper content"
  finally:
    shutil.rmtree(tmp)


def exercise_disabled_skill_excluded():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "a")
    _make_skill(builtin, "b")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default(disabled={"b"})
    names = {s.name for s in skills}
    assert names == {"a"}
  finally:
    shutil.rmtree(tmp)


def exercise_tools_returns_empty_for_no_skills():
  """Section 8.5: the skill tools are only useful when skills exist;
  with an empty list, tools() returns no entries so we don't pollute
  the registry with handlers that can't do anything."""
  loader = SkillLoader(builtin_path=Path("/nonexistent"), log=null_out())
  assert loader.tools([]) == []


def exercise_tools_returns_two_when_all_always_mode():
  """When every loaded skill is mode=always, the body is already in the
  system prompt so load_skill has nothing to fetch and is suppressed.
  Only read_skill_file and list_skill_files are returned."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "a")
    _make_skill(builtin, "b")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    tools = loader.tools(skills)
    names = sorted(spec.name for spec, _ in tools)
    assert names == ["list_skill_files", "read_skill_file"], names
  finally:
    shutil.rmtree(tmp)


def exercise_tools_returns_three_when_any_on_demand():
  """As soon as any skill is mode=on_demand, load_skill becomes the only
  way for the agent to see that body — must be registered."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "a", mode="always")
    _make_skill(builtin, "b", mode="on_demand")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    tools = loader.tools(skills)
    names = sorted(spec.name for spec, _ in tools)
    assert names == ["list_skill_files", "load_skill", "read_skill_file"], \
      names
  finally:
    shutil.rmtree(tmp)


def exercise_tools_input_schema_enumerates_loaded_skills():
  """Each tool's input_schema constrains skill_name to the loaded names
  so the agent gets typeahead-style guidance and can't ask for skills
  that don't exist."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "first", mode="on_demand")
    _make_skill(builtin, "second")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    tools = loader.tools(skills)
    for spec, _ in tools:
      enum = (spec.input_schema["properties"]["skill_name"]["enum"])
      assert sorted(enum) == ["first", "second"], (spec.name, enum)
  finally:
    shutil.rmtree(tmp)


def exercise_read_file_handler_returns_string():
  """The read_skill_file handler decodes UTF-8 text so the agent gets
  a string (not bytes) for a typical markdown helper."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    skill_dir = _make_skill(builtin, "x")
    with open(os.path.join(skill_dir, "note.md"), "wb") as fh:
      fh.write(b"citation list\n")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    tools = dict((spec.name, h) for spec, h in loader.tools(skills))
    out = tools["read_skill_file"](
      "read_skill_file", {"skill_name": "x", "relative_path": "note.md"})
    assert out == "citation list\n", repr(out)
  finally:
    shutil.rmtree(tmp)


def exercise_handlers_raise_sorry_for_unknown_skill_name():
  """Unknown skill names go through the loader's lookup, which raises
  Sorry rather than silently returning empty. Tested for read, list,
  and load_skill (since all three reach by_name)."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "x", mode="on_demand")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    tools = dict((spec.name, h) for spec, h in loader.tools(skills))
    for tool_name in ("read_skill_file", "list_skill_files", "load_skill"):
      try:
        if tool_name == "read_skill_file":
          tools[tool_name](tool_name, {
            "skill_name": "ghost", "relative_path": "x"})
        else:
          tools[tool_name](tool_name, {"skill_name": "ghost"})
      except Sorry:
        pass
      else:
        raise AssertionError(
          "expected Sorry from %s for unknown skill" % tool_name)
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_drops_skill_when_server_missing():
  """When the skill declares `requires: [phenix]` and the loader is told
  the phenix MCP server isn't available, the skill is filtered out and
  a message is logged. Without this, the agent would load guidance that
  references tools it can't actually call."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "needs_phenix", requires=["phenix"])
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default(mcp_servers=[])
    assert skills == [], [s.name for s in skills]
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_keeps_skill_when_server_present():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "needs_phenix", requires=["phenix"])
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default(mcp_servers=["phenix"])
    assert [s.name for s in skills] == ["needs_phenix"]
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_none_opts_out_of_check():
  """Backward compatibility: callers that don't pass mcp_servers (or
  pass None) bypass the check entirely and load every skill. Existing
  callers / tests that don't know which servers will run depend on this."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "needs_phenix", requires=["phenix"])
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()                    # no mcp_servers arg
    assert [s.name for s in skills] == ["needs_phenix"]
    skills2 = loader.load_default(mcp_servers=None)   # explicit None
    assert [s.name for s in skills2] == ["needs_phenix"]
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_passes_skill_with_no_requires():
  """Skills without a `requires` field are transport-agnostic and load
  even when no MCP servers are declared available."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "no_deps")                   # no requires kwarg
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default(mcp_servers=[])
    assert [s.name for s in skills] == ["no_deps"]
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_partial_match_drops_skill():
  """When a skill requires two servers and only one is available, it
  still gets dropped — the skill body assumes both."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "needs_both", requires=["phenix", "other"])
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default(mcp_servers=["phenix"])
    assert skills == []
  finally:
    shutil.rmtree(tmp)


def exercise_requires_filter_applies_to_additional_skills():
  """Same filter logic applies to non-builtin (project / user) skills
  resolved via `additional`."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    user = os.path.join(tmp, "user")
    os.makedirs(builtin)
    os.makedirs(user)
    _make_skill(user, "user_phenix", requires=["phenix"])
    loader = SkillLoader(
      builtin_path=Path(builtin), user_path=Path(user), log=null_out())
    skills = loader.load_default(
      additional=["user_phenix"], mcp_servers=[])
    assert skills == []
    skills_with = loader.load_default(
      additional=["user_phenix"], mcp_servers=["phenix"])
    assert [s.name for s in skills_with] == ["user_phenix"]
  finally:
    shutil.rmtree(tmp)


def exercise_assemble_system_prompt_includes_descriptions():
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "skill_a",
                description="does A things",
                body="UNIQUE_BODY_FOR_SKILL_A\n")
    _make_skill(builtin, "skill_b",
                description="does B things",
                mode="on_demand",
                body="UNIQUE_BODY_FOR_SKILL_B\n")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    prompt = loader.assemble_system_prompt("BASE PROMPT", skills)
    assert "BASE PROMPT" in prompt
    assert "does A things" in prompt
    assert "does B things" in prompt
    # always-mode skill body is inlined
    assert "UNIQUE_BODY_FOR_SKILL_A" in prompt
    # on_demand-mode skill body must NOT be inlined
    assert "UNIQUE_BODY_FOR_SKILL_B" not in prompt
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_load_builtin_only()
  exercise_builtin_wins_on_collision()
  exercise_read_file_rejects_traversal()
  exercise_read_file_rejects_absolute_path()
  exercise_read_file_rejects_symlink_escape()
  exercise_read_file_succeeds_for_in_skill_file()
  exercise_disabled_skill_excluded()
  exercise_assemble_system_prompt_includes_descriptions()
  exercise_tools_returns_empty_for_no_skills()
  exercise_tools_returns_two_when_all_always_mode()
  exercise_tools_returns_three_when_any_on_demand()
  exercise_tools_input_schema_enumerates_loaded_skills()
  exercise_read_file_handler_returns_string()
  exercise_handlers_raise_sorry_for_unknown_skill_name()
  exercise_requires_filter_drops_skill_when_server_missing()
  exercise_requires_filter_keeps_skill_when_server_present()
  exercise_requires_filter_none_opts_out_of_check()
  exercise_requires_filter_passes_skill_with_no_requires()
  exercise_requires_filter_partial_match_drops_skill()
  exercise_requires_filter_applies_to_additional_skills()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
