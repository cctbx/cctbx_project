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


def _make_undecodable_skill(root, name):
  """Write a SKILL.md whose bytes are not valid UTF-8 (a lone 0xFF), so
  reading it raises UnicodeDecodeError under any locale. Used to prove one
  unreadable skill is skipped rather than aborting the whole load."""
  d = os.path.join(root, name)
  os.makedirs(d)
  with open(os.path.join(d, "SKILL.md"), "wb") as fh:
    fh.write(b"---\nname: bad\ndescription: x\nmode: always\n---\n\n"
             b"\xff\xfe not valid utf-8\n")
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


def exercise_one_unreadable_skill_does_not_drop_the_rest():
  """A single malformed / undecodable SKILL.md must be SKIPPED, not abort
  the whole load. Before the fix the UnicodeDecodeError escaped the
  per-skill `except Sorry` guard and dropped every other skill too (F6).
  The bad skill sorts AFTER the good one, so a passing result also proves
  the good skill is not lost when a later skill blows up."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "good_skill")
    _make_undecodable_skill(builtin, "zz_bad_skill")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    assert [s.name for s in skills] == ["good_skill"], [s.name for s in skills]
  finally:
    shutil.rmtree(tmp)


def exercise_non_ascii_utf8_skill_loads():
  """A SKILL.md containing non-ASCII UTF-8 text loads cleanly and its body /
  description round-trip. Guards the encoding='utf-8' pin in _parse_skill_md
  (F6): without it, read_text() uses the locale default, so a non-ASCII
  SKILL.md raises UnicodeDecodeError under a C/POSIX locale."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    d = os.path.join(builtin, "unicode_skill")
    os.makedirs(d)
    skill_md = (
      "---\n"
      "name: unicode_skill\n"
      "description: résumé of ångström refinement\n"
      "mode: always\n"
      "---\n\n"
      "Body with non-ASCII: éåö — café.\n")
    with open(os.path.join(d, "SKILL.md"), "w", encoding="utf-8") as fh:
      fh.write(skill_md)
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()
    assert len(skills) == 1, [s.name for s in skills]
    assert skills[0].name == "unicode_skill"
    assert "café" in skills[0].body, skills[0].body
    assert "résumé" in skills[0].description, skills[0].description
  finally:
    shutil.rmtree(tmp)


def exercise_scalar_requires_coerced_to_single_element_list():
  """A SKILL.md with a SCALAR `requires: phenix` (not a YAML list) keeps that
  one requirement intact. PyYAML parses the scalar as the string "phenix", and
  the old `list(frontmatter.get("requires", []) or [])` would explode it into
  ['p','h','e','n','i','x'] -- none of which match a real MCP server, so the
  skill was silently dropped with a garbage 'missing servers' log (F3). The
  scalar must coerce to the one-element list ['phenix']: the skill loads when
  phenix is available and is filtered when it isn't."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    d = os.path.join(builtin, "scalar_req")
    os.makedirs(d)
    skill_md = (
      "---\n"
      "name: scalar_req\n"
      "description: needs phenix\n"
      "mode: always\n"
      "requires: phenix\n"               # scalar, NOT [phenix]
      "---\n\nbody\n")
    with open(os.path.join(d, "SKILL.md"), "w") as fh:
      fh.write(skill_md)
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    # requires must be the one-element list, not exploded into characters.
    skill = loader.load_default()[0]      # no mcp filter -> loads regardless
    assert skill.requires == ["phenix"], skill.requires
    # ... and the requires filter treats it as a single 'phenix' requirement.
    assert [s.name for s in loader.load_default(mcp_servers=["phenix"])] \
      == ["scalar_req"]
    assert loader.load_default(mcp_servers=[]) == []
  finally:
    shutil.rmtree(tmp)


def exercise_non_string_scalar_requires_does_not_drop_batch():
  """[F3] A NON-string scalar `requires` (e.g. `requires: 1` or `requires: true`)
  must coerce to a one-element list and NOT abort the whole skill load. An even
  earlier round's coercion handled only str; a non-str scalar still reached
  list(1) / list(True) -> TypeError, which is NOT in _load_checked's
  (Sorry, OSError, ValueError) guard, so it escaped _load_builtins ->
  load_default raised -> the session caught it as 'Failed to load skills
  (continuing without)' and lost EVERY skill.

  requires entries are server NAMES, so each element is now str()-coerced after
  the list-wrap: the scalar 1 becomes ['1'] and true becomes ['True'] (a single
  hashable, non-matching requirement) and the rest of the batch survives."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "good_skill")              # sorts before the scalar ones
    # Raw SKILL.md with a scalar int / bool `requires` (NOT a YAML list).
    for name, scalar in (("zz_int_req", "1"), ("zz_bool_req", "true")):
      d = os.path.join(builtin, name)
      os.makedirs(d)
      with open(os.path.join(d, "SKILL.md"), "w") as fh:
        fh.write("---\nname: %s\ndescription: d\nmode: always\n"
                 "requires: %s\n---\n\nbody\n" % (name, scalar))
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    skills = loader.load_default()                  # must NOT raise; keeps all
    by_name = {s.name: s for s in skills}
    assert "good_skill" in by_name, list(by_name)   # batch not dropped
    # Each scalar coerces to its single-element STRING list -- not exploded,
    # not dropped, and uniformly hashable for the requires membership test.
    assert by_name["zz_int_req"].requires == ["1"], by_name["zz_int_req"].requires
    assert by_name["zz_bool_req"].requires == ["True"], \
      by_name["zz_bool_req"].requires
  finally:
    shutil.rmtree(tmp)


def exercise_mapping_requires_coerced_and_does_not_abort_batch():
  """[F#1] A YAML MAPPING `requires: {phenix: null}` (or a list element that is
  itself a dict/list) parses to an UNHASHABLE requires value. This round's bare
  `[requires]` list-wrap left that dict in place, and _requires_satisfied's
  `r in available_servers` (frozenset membership) then raised TypeError --
  OUTSIDE _load_checked's (Sorry, OSError, ValueError) guard -- so it escaped
  _load_builtins and aborted the WHOLE skill batch (every skill lost), strictly
  worse than the pre-this-round clean single-skill drop. str()-coercing each
  element makes any value hashable: the skill loads (with a stringified,
  non-matching requirement) and the batch survives.

  Revert-proof: without the `[str(r) for r in requires]` coercion, the
  mapping skill's membership test raises TypeError out of load_default (it runs
  because mcp_servers is non-None), so the first load_default below raises and
  none of the asserts run -- good_skill is lost with it."""
  tmp = tempfile.mkdtemp()
  try:
    builtin = os.path.join(tmp, "builtin")
    os.makedirs(builtin)
    _make_skill(builtin, "good_skill")              # sorts before the zz_* ones
    # A MAPPING `requires` (unhashable dict element) ...
    d = os.path.join(builtin, "zz_map_req")
    os.makedirs(d)
    with open(os.path.join(d, "SKILL.md"), "w") as fh:
      fh.write("---\nname: zz_map_req\ndescription: d\nmode: always\n"
               "requires: {phenix: null}\n---\n\nbody\n")
    # ... and a LIST whose elements include a non-string scalar.
    d2 = os.path.join(builtin, "zz_list_nonstr")
    os.makedirs(d2)
    with open(os.path.join(d2, "SKILL.md"), "w") as fh:
      fh.write("---\nname: zz_list_nonstr\ndescription: d\nmode: always\n"
               "requires: [1, phenix]\n---\n\nbody\n")
    loader = SkillLoader(builtin_path=Path(builtin), log=null_out())
    # mcp_servers is non-None so the requires membership test actually RUNS --
    # that membership test is where the unhashable dict used to raise TypeError.
    skills = loader.load_default(mcp_servers=["phenix"])
    by_name = {s.name: s for s in skills}
    # The batch is NOT dropped: good_skill survives the malformed siblings.
    assert "good_skill" in by_name, list(by_name)
    # Each malformed skill is cleanly FILTERED (its stringified requirement
    # doesn't match a real server) rather than aborting the load.
    assert "zz_map_req" not in by_name, list(by_name)
    assert "zz_list_nonstr" not in by_name, list(by_name)
    # With the check opted out (mcp_servers=None) the skills load and EVERY
    # requires element is a string -- proving the load-time stringification.
    everything = {s.name: s for s in loader.load_default()}
    assert all(isinstance(r, str) for r in everything["zz_map_req"].requires), \
      everything["zz_map_req"].requires
    assert everything["zz_list_nonstr"].requires == ["1", "phenix"], \
      everything["zz_list_nonstr"].requires
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
  exercise_one_unreadable_skill_does_not_drop_the_rest()
  exercise_non_ascii_utf8_skill_loads()
  exercise_scalar_requires_coerced_to_single_element_list()
  exercise_non_string_scalar_requires_does_not_drop_batch()
  exercise_mapping_requires_coerced_and_does_not_abort_batch()
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
