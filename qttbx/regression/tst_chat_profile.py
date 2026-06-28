import json
import os
import shutil
import tempfile
from pathlib import Path

from libtbx.test_utils import Exception_expected
from libtbx.utils import format_cpu_times, Sorry, null_out
from qttbx.widgets.chat.agent.profile import ProfileLoader


def _write_profile(dirpath, name, body):
  Path(dirpath).mkdir(parents=True, exist_ok=True)
  with open(os.path.join(dirpath, name + ".json"), "w") as fh:
    json.dump(body, fh)


def exercise_minimal_profile_loads():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "minimal",
                   {"name": "minimal", "model": "claude-opus-4-7"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("minimal")
    assert p.name == "minimal"
    assert p.model == "claude-opus-4-7"
    # All other fields take defaults
    assert p.max_tokens == 8192
    assert p.tool_policy_default == "ask"
    assert p.vision_input is True
    assert p.skills_additional == []
    assert p.skills_disabled == set()
    assert p.subagents_enabled is True
    assert p.subagents_default_model == "claude-opus-4-7"
  finally:
    shutil.rmtree(tmp)


def exercise_missing_required_field_raises():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "bad", {"name": "bad"})  # missing model
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load("bad")
    except Sorry:
      pass
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_based_on_inheritance():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "parent", {
      "name": "parent",
      "model": "claude-opus-4-7",
      "max_tokens": 4096,
      "tool_policy_default": "allow",
    })
    _write_profile(tmp, "child", {
      "name": "child",
      "based_on": "parent",
      "model": "claude-haiku-4-5",   # override
      # max_tokens inherited from parent
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("child")
    assert p.name == "child"
    assert p.model == "claude-haiku-4-5"          # override
    assert p.max_tokens == 4096                   # inherited
    assert p.tool_policy_default == "allow"       # inherited
  finally:
    shutil.rmtree(tmp)


def exercise_based_on_child_system_prompt_file_overrides_parent_prompt():
  """A child may override an inherited inline system_prompt with its own
  system_prompt_file. The parent's resolved system_prompt must not survive
  the inheritance merge -- otherwise the merged profile carries BOTH a
  system_prompt and a system_prompt_file and trips the mutual-exclusivity
  gate in _build_profile, so the chat window never opens (F11)."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "parent", {
      "name": "parent",
      "model": "claude-opus-4-7",
      "system_prompt": "PARENT INLINE PROMPT",
    })
    with open(os.path.join(tmp, "child_prompt.md"), "w") as fh:
      fh.write("CHILD PROMPT FROM FILE")
    _write_profile(tmp, "child", {
      "name": "child",
      "based_on": "parent",
      "system_prompt_file": "${PROFILE_DIR}/child_prompt.md",
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("child")          # must not raise Sorry
    assert p.system_prompt == "CHILD PROMPT FROM FILE", p.system_prompt
  finally:
    shutil.rmtree(tmp)


def exercise_based_on_child_system_prompt_overrides_parent_prompt_file():
  """A parent's system_prompt_file is inlined into system_prompt when the
  PARENT loads, so by merge time the parent contributes only a resolved inline
  system_prompt (never a system_prompt_file key). A child supplying its own
  inline system_prompt therefore simply overrides that inherited prompt --
  ordinary child-over-parent. This is NOT the symmetric counterpart of the F11
  eviction: no system_prompt_file ever survives into the merged dict, so there
  is nothing to evict and the mutual-exclusivity gate is never even approached
  (this case passes with or without the eviction). The gate itself is covered
  by exercise_based_on_does_not_weaken_mutual_exclusivity_gate."""
  tmp = tempfile.mkdtemp()
  try:
    with open(os.path.join(tmp, "parent_prompt.md"), "w") as fh:
      fh.write("PARENT PROMPT FROM FILE")
    _write_profile(tmp, "parent", {
      "name": "parent",
      "model": "claude-opus-4-7",
      "system_prompt_file": "${PROFILE_DIR}/parent_prompt.md",
    })
    _write_profile(tmp, "child", {
      "name": "child",
      "based_on": "parent",
      "system_prompt": "CHILD INLINE PROMPT",
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("child")          # must not raise Sorry
    assert p.system_prompt == "CHILD INLINE PROMPT", p.system_prompt
  finally:
    shutil.rmtree(tmp)


def exercise_based_on_does_not_weaken_mutual_exclusivity_gate():
  """The based_on merge evicts the PARENT's inherited prompt so a child can
  switch prompt mechanisms (F11) -- but a profile that GENUINELY declares both
  system_prompt AND system_prompt_file must still raise. Covers a standalone
  profile and a based_on child that each declare both."""
  tmp = tempfile.mkdtemp()
  try:
    # (i) Standalone profile declaring both keys. The prompt file genuinely
    # exists and stays in-bounds, so the raised Sorry is unambiguously the
    # mutual-exclusivity gate -- not a missing-file or directory-escape error.
    with open(os.path.join(tmp, "both_prompt.md"), "w") as fh:
      fh.write("FILE PROMPT")
    _write_profile(tmp, "both", {
      "name": "both",
      "model": "claude-opus-4-7",
      "system_prompt": "INLINE PROMPT",
      "system_prompt_file": "${PROFILE_DIR}/both_prompt.md",
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load("both")
    except Sorry as e:
      assert "mutually exclusive" in str(e), str(e)
    else:
      raise Exception_expected

    # (ii) A based_on child that ITSELF declares both keys. The F11 eviction
    # pops only the PARENT's inherited system_prompt; the child's OWN
    # system_prompt then overlays back in next to its system_prompt_file, so
    # _build_profile still sees both and the gate must fire. The eviction
    # narrows what is inherited -- it does not suppress the gate.
    with open(os.path.join(tmp, "child_prompt.md"), "w") as fh:
      fh.write("CHILD FILE PROMPT")
    _write_profile(tmp, "parent", {
      "name": "parent",
      "model": "claude-opus-4-7",
      "system_prompt": "PARENT INLINE PROMPT",
    })
    _write_profile(tmp, "child", {
      "name": "child",
      "based_on": "parent",
      "system_prompt": "CHILD INLINE PROMPT",
      "system_prompt_file": "${PROFILE_DIR}/child_prompt.md",
    })
    try:
      loader.load("child")
    except Sorry as e:
      assert "mutually exclusive" in str(e), str(e)
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_based_on_cycle_detected():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "a",
                   {"name": "a", "model": "m", "based_on": "b"})
    _write_profile(tmp, "b",
                   {"name": "b", "model": "m", "based_on": "a"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load("a")
    except Sorry as e:
      assert "cycle" in str(e).lower() or "inherit" in str(e).lower()
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_load_file_from_explicit_path():
  """load_file reads a profile from an explicit path, not by name."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "ext", {"name": "ext", "model": "m"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load_file(os.path.join(tmp, "ext.json"))
    assert p.name == "ext"
    assert p.model == "m"
  finally:
    shutil.rmtree(tmp)


def exercise_load_file_based_on_sibling_wins_over_builtin():
  """A file's based_on resolves a sibling in the file's OWN directory, and that
  sibling takes PRECEDENCE over a same-named profile in a configured dir.

  Both dirs hold a profile named "base" with different values, so only search
  priority decides the outcome -- this locks in 'extra dir searched first',
  catching a regression that appended instead of prepended (or dropped it)."""
  bundle = tempfile.mkdtemp()
  builtin = tempfile.mkdtemp()
  try:
    _write_profile(builtin, "base", {"name": "base", "model": "m",
                                     "tool_policy_default": "ask"})
    _write_profile(bundle, "base", {"name": "base", "model": "m",
                                    "tool_policy_default": "allow"})
    _write_profile(bundle, "ext", {"name": "ext", "based_on": "base"})
    loader = ProfileLoader(builtin_dir=Path(builtin), log=null_out())
    p = loader.load_file(os.path.join(bundle, "ext.json"))
    assert p.name == "ext"
    # The SIBLING base (allow) must win over the builtin base (ask).
    assert p.tool_policy_default == "allow", p.tool_policy_default
  finally:
    shutil.rmtree(bundle)
    shutil.rmtree(builtin)


def exercise_load_file_does_not_mutate_loader_search_dirs():
  """load_file is stateless: the external file's directory is scoped to that
  call and does not leak into a later load() on the same loader instance."""
  bundle = tempfile.mkdtemp()
  builtin = tempfile.mkdtemp()
  try:
    _write_profile(bundle, "leak", {"name": "leak", "model": "m"})
    _write_profile(bundle, "ext", {"name": "ext", "model": "m"})
    _write_profile(builtin, "real", {"name": "real", "model": "m"})
    loader = ProfileLoader(builtin_dir=Path(builtin), log=null_out())
    loader.load_file(os.path.join(bundle, "ext.json"))   # bundle scoped to here
    # A later load() must NOT resolve names from the bundle dir.
    try:
      loader.load("leak")
    except Sorry as e:
      assert "not found" in str(e).lower(), str(e)
    else:
      raise Exception_expected
    # ... and a genuinely configured (builtin) profile still loads fine.
    assert loader.load("real").name == "real"
  finally:
    shutil.rmtree(bundle)
    shutil.rmtree(builtin)


def exercise_load_file_based_on_falls_through_to_builtin():
  """When the parent is not a sibling, based_on still resolves against the
  configured search dirs (so an external file can inherit from a builtin)."""
  bundle = tempfile.mkdtemp()
  builtin = tempfile.mkdtemp()
  try:
    _write_profile(builtin, "base", {"name": "base", "model": "m",
                                     "max_tokens": 4096})
    _write_profile(bundle, "ext", {"name": "ext", "based_on": "base"})
    loader = ProfileLoader(builtin_dir=Path(builtin), log=null_out())
    p = loader.load_file(os.path.join(bundle, "ext.json"))
    assert p.max_tokens == 4096   # inherited from the builtin parent
  finally:
    shutil.rmtree(bundle)
    shutil.rmtree(builtin)


def exercise_load_file_missing_raises():
  """load_file on a path that is not a file raises a clear Sorry."""
  tmp = tempfile.mkdtemp()
  try:
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load_file(os.path.join(tmp, "nope.json"))
    except Sorry as e:
      assert "not found" in str(e).lower(), str(e)
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_variable_expansion_env():
  tmp = tempfile.mkdtemp()
  saved = os.environ.get("PHENIX_TEST_VAR")
  os.environ["PHENIX_TEST_VAR"] = "hello"
  try:
    _write_profile(tmp, "p", {
      "name": "p", "model": "m",
      "mcp_servers": [{
        "name": "test",
        "command": "${env:PHENIX_TEST_VAR}",
        "args": ["--flag=${env:PHENIX_TEST_VAR}"],
      }],
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")
    assert len(p.mcp_servers) == 1
    srv = p.mcp_servers[0]
    assert srv.command == "hello"
    assert srv.args == ["--flag=hello"]
  finally:
    if saved is None:
      del os.environ["PHENIX_TEST_VAR"]
    else:
      os.environ["PHENIX_TEST_VAR"] = saved
    shutil.rmtree(tmp)


def exercise_project_overrides_user_overrides_builtin():
  builtin = tempfile.mkdtemp()
  user = tempfile.mkdtemp()
  project = tempfile.mkdtemp()
  try:
    _write_profile(builtin, "phenix_expert",
                   {"name": "phenix_expert", "model": "claude-opus-4-7",
                    "max_tokens": 1})
    _write_profile(user, "phenix_expert",
                   {"name": "phenix_expert", "model": "claude-opus-4-7",
                    "max_tokens": 2})
    _write_profile(project, "phenix_expert",
                   {"name": "phenix_expert", "model": "claude-opus-4-7",
                    "max_tokens": 3})
    loader = ProfileLoader(builtin_dir=Path(builtin),
                           user_dir=Path(user),
                           project_dir=Path(project),
                           log=null_out())
    p = loader.load("phenix_expert")
    assert p.max_tokens == 3  # project wins
  finally:
    shutil.rmtree(builtin)
    shutil.rmtree(user)
    shutil.rmtree(project)


def exercise_unknown_fields_warn_not_fail():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {
      "name": "p", "model": "m",
      "future_field_we_dont_know_about": True,
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")  # should not raise
    assert p.name == "p"
  finally:
    shutil.rmtree(tmp)


def exercise_system_prompt_file_with_expansion():
  """Section 13.3: system_prompt_file is run through ${VAR} expansion."""
  tmp = tempfile.mkdtemp()
  try:
    prompt_path = os.path.join(tmp, "my_prompt.md")
    with open(prompt_path, "w") as fh:
      fh.write("You are a helpful assistant.")
    _write_profile(tmp, "p", {
      "name": "p", "model": "m",
      "system_prompt_file": "${PROFILE_DIR}/my_prompt.md",
    })
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")
    assert p.system_prompt == "You are a helpful assistant."
  finally:
    shutil.rmtree(tmp)


def exercise_system_prompt_file_non_ascii_utf8_round_trips():
  """A system_prompt_file holding non-ASCII UTF-8 text loads and round-trips
  regardless of the platform's default text encoding. read_text() with no
  explicit encoding uses that default (ASCII under a C/POSIX locale), which
  raises UnicodeDecodeError on a non-ASCII byte -- so the prompt-file read is
  pinned to encoding='utf-8' (matching skills._parse_skill_md). The Path
  monkeypatch below simulates an ASCII default so this also fails on a UTF-8
  dev machine if the pin is dropped."""
  tmp = tempfile.mkdtemp()
  orig_read_text = Path.read_text

  def _ascii_default_read_text(self, encoding=None, *args, **kwargs):
    # Faithfully model the C/POSIX locale: a read with no explicit encoding
    # defaults to ASCII -- the bug condition the utf-8 pin guards against.
    if encoding is None:
      encoding = "ascii"
    return orig_read_text(self, encoding, *args, **kwargs)

  try:
    prompt_path = os.path.join(tmp, "prompt.md")
    non_ascii = "You are a résumé / ångström assistant — café.\n"
    with open(prompt_path, "w", encoding="utf-8") as fh:
      fh.write(non_ascii)
    # The bytes are genuinely non-ASCII, so an ASCII decode of the file fails;
    # that is exactly the failure the utf-8 pin prevents.
    try:
      Path(prompt_path).read_bytes().decode("ascii")
    except UnicodeDecodeError:
      pass
    else:
      raise Exception_expected
    _write_profile(tmp, "p", {
      "name": "p", "model": "m",
      "system_prompt_file": "${PROFILE_DIR}/prompt.md"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    Path.read_text = _ascii_default_read_text
    p = loader.load("p")
    assert p.system_prompt == non_ascii, repr(p.system_prompt)
  finally:
    Path.read_text = orig_read_text
    shutil.rmtree(tmp)


def exercise_profile_json_non_ascii_utf8_loads():
  """A profile JSON containing non-ASCII UTF-8 (e.g. an accented description)
  loads regardless of the platform's default text encoding. _load_path read
  the file with `open(path)` -- no encoding -- so under a C/POSIX locale
  (ASCII default) a non-ASCII profile raised UnicodeDecodeError, surfacing as
  a confusing 'Profile parse error' (F9). The read is pinned to
  encoding='utf-8' (matching profile.py's :442 read_text and
  skills._parse_skill_md). The builtins.open monkeypatch models an ASCII
  default so this also fails on a UTF-8 dev machine if the pin is dropped."""
  import builtins
  tmp = tempfile.mkdtemp()
  orig_open = builtins.open

  def _ascii_default_open(file, mode="r", *args, **kwargs):
    # Faithfully model the C/POSIX locale: a text read with no explicit
    # encoding defaults to ASCII -- the bug condition the utf-8 pin guards.
    if "b" not in mode and "encoding" not in kwargs and not args:
      kwargs["encoding"] = "ascii"
    return orig_open(file, mode, *args, **kwargs)

  try:
    desc = "résumé / ångström refinement — café"
    path = os.path.join(tmp, "uni.json")
    # ensure_ascii=False writes the description as LITERAL UTF-8 bytes (how a
    # human-authored profile naturally carries non-ASCII), not \uXXXX escapes
    # -- so the on-disk file genuinely needs a UTF-8 decode.
    with open(path, "w", encoding="utf-8") as fh:
      json.dump({"name": "uni", "model": "m", "description": desc}, fh,
                ensure_ascii=False)
    # The bytes are genuinely non-ASCII: an ASCII decode of the file fails;
    # that is exactly the failure the utf-8 pin prevents.
    try:
      Path(path).read_bytes().decode("ascii")
    except UnicodeDecodeError:
      pass
    else:
      raise Exception_expected
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    builtins.open = _ascii_default_open
    p = loader.load("uni")          # must not raise under an ASCII default
    assert p.description == desc, repr(p.description)
  finally:
    builtins.open = orig_open
    shutil.rmtree(tmp)


def exercise_profile_backend_defaults_to_claude_code():
  """Profiles with no `backend:` key default to claude_code so a user
  with `claude login` already done can chat with no extra setup. Users
  who want the API-key path must opt in with `backend: anthropic`."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "claude-opus-4-7"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")
    assert p.backend == "claude_code", p.backend
    # The claude.* subtree is exposed with safe defaults so callers can
    # read profile.claude.cli_path / sdk_options unconditionally.
    assert p.claude.cli_path is None, p.claude.cli_path
    assert p.claude.sdk_options == {}, p.claude.sdk_options
  finally:
    shutil.rmtree(tmp)


def exercise_profile_backend_claude_code_parsed():
  """backend: claude_code parses cleanly and the optional claude.*
  subtree is exposed on the loaded Profile."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {
      "name": "p", "model": "claude-opus-4-7",
      "backend": "claude_code",
      "claude": {"cli_path": "/opt/claude/bin/claude",
                 "sdk_options": {"max_turns": 20}}})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")
    assert p.backend == "claude_code", p.backend
    assert p.claude.cli_path == "/opt/claude/bin/claude", p.claude.cli_path
    assert p.claude.sdk_options == {"max_turns": 20}, p.claude.sdk_options
  finally:
    shutil.rmtree(tmp)


def exercise_openai_is_known_backend():
  """The openai backend must validate at load time (it is a real agent
  backend now), mirroring the claude_code/anthropic known-backend cases."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p",
                   {"name": "p", "model": "gpt-x", "backend": "openai"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    prof = loader.load("p")
    assert prof.backend == "openai", prof.backend
  finally:
    shutil.rmtree(tmp)


def exercise_portkey_is_known_backend():
  """portkey validates at load time, like the other real backends."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "gpt-x", "backend": "portkey"})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.backend == "portkey", prof.backend
  finally:
    shutil.rmtree(tmp)


def exercise_google_is_known_backend():
  """google validates at load time, like the other real backends."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "gemini-x", "backend": "google"})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.backend == "google", prof.backend
  finally:
    shutil.rmtree(tmp)


def exercise_portkey_block_parsed():
  """A `portkey` JSON sub-block populates portkey_virtual_key / portkey_config."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "gpt-x", "backend": "portkey",
                              "portkey": {"virtual_key": "vk", "config": "cfg"}})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.portkey_virtual_key == "vk", prof.portkey_virtual_key
    assert prof.portkey_config == "cfg", prof.portkey_config
  finally:
    shutil.rmtree(tmp)


def exercise_portkey_block_absent_defaults_none():
  """No `portkey` block -> the two attrs default to None (not missing/crash)."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "gpt-x", "backend": "portkey"})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert getattr(prof, "portkey_virtual_key", "MISSING") is None
    assert getattr(prof, "portkey_config", "MISSING") is None
  finally:
    shutil.rmtree(tmp)


def exercise_profile_backend_unknown_raises():
  """Unknown backend values must fail at load time, not at runtime, so
  the user gets the error before they wait for the chat window."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {
      "name": "p", "model": "claude-opus-4-7", "backend": "bogus"})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load("p")
    except Sorry as e:
      assert "backend" in str(e), str(e)
      assert "bogus" in str(e), str(e)
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_server_tools_string_raises_clear_sorry():
  """A bare string `server_tools` (a common mistake) raises a clear Sorry
  rather than silently splitting into single characters ('w', 'e', 'b', ...)
  that then warn as ~10 unknown tools and enable none."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "m",
                              "server_tools": "web_search"})
    try:
      ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
      raise Exception_expected
    except Sorry as e:
      assert "server_tools" in str(e) and "list" in str(e), str(e)
  finally:
    shutil.rmtree(tmp)


def exercise_server_tools_parsed_and_defaults_empty():
  """An opt-in `server_tools` list parses into Profile.server_tools; a
  profile that omits it gets an empty list (pure no-op default)."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "m",
                              "server_tools": ["web_search", "code_execution"]})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.server_tools == ["web_search", "code_execution"], prof.server_tools
    _write_profile(tmp, "q", {"name": "q", "model": "m"})
    assert ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("q").server_tools == []
  finally:
    shutil.rmtree(tmp)


def exercise_skills_scalar_additional_disabled_coerced_to_single_id():
  """Regression: a scalar `skills.additional: myskill` (or `disabled: bad`) is
  a single skill id, not a sequence of characters. Without a scalar guard the
  bare list()/set() explode it into ['m','y','s','k','i','l','l'] / a char-set,
  silently dropping the real skill (and disabling phantom one-char ids). Coerce
  a scalar to a one-element container so the id survives intact -- the same
  forgiving coercion already applied to a skill's `requires` frontmatter."""
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "m",
                              "skills": {"additional": "myskill",
                                         "disabled": "badskill"}})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.skills_additional == ["myskill"], prof.skills_additional
    assert prof.skills_disabled == {"badskill"}, prof.skills_disabled
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_server_env_cannot_override_path_or_phenix_vars():
  """A profile may be project-supplied (untrusted). Its per-server env
  must not override PATH / dynamic-loader / PHENIX_* vars and thereby
  redirect which binary runs or escape the project sandbox. Legitimate
  per-server vars are preserved."""
  from qttbx.widgets.chat.agent.profile import sanitize_server_env
  out = sanitize_server_env({
    "PATH": "/tmp/evil", "PHENIX_PROJECT_DIR": "/elsewhere",
    "LD_PRELOAD": "/tmp/x.so", "API_TOKEN": "keep"})
  assert out == {"API_TOKEN": "keep"}, out
  # Broader runtime-loader / interpreter families are stripped too,
  # including PHENIX_TRUST_OTHER_ENV (which would disable the dispatcher
  # env scrub). A harmless var (PYTHONUNBUFFERED) and a server's own
  # token survive -- we strip the dangerous keys, not all PYTHON*.
  out2 = sanitize_server_env({
    "PYTHONPATH": "/evil", "PYTHONHOME": "/evil",
    "NODE_OPTIONS": "--require /evil.js", "NODE_PATH": "/evil",
    "BASH_ENV": "/evil.sh", "ENV": "/evil.sh", "IFS": " ",
    "LD_AUDIT": "/evil.so", "DYLD_FRAMEWORK_PATH": "/evil",
    "PERL5OPT": "-M", "RUBYOPT": "-r/evil",
    "BASH_FUNC_x%%": "() { :; }", "PHENIX_TRUST_OTHER_ENV": "1",
    "PYTHONUNBUFFERED": "1", "MY_SERVER_TOKEN": "keep"})
  assert out2 == {"PYTHONUNBUFFERED": "1", "MY_SERVER_TOKEN": "keep"}, out2
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {
      "name": "p", "model": "m",
      "mcp_servers": [{
        "name": "x", "command": "phenix.mcp_server",
        "env": {"PATH": "/tmp/evil", "PHENIX_PROJECT_DIR": "/x",
                "API_TOKEN": "keep"}}]})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    p = loader.load("p")
    env = p.mcp_servers[0].env
    assert "PATH" not in env, env
    assert "PHENIX_PROJECT_DIR" not in env, env
    assert env.get("API_TOKEN") == "keep", env
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_server_inject_phenix_env_default_true():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "m",
      "mcp_servers": [{"name": "phenix", "command": "x"}]})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.mcp_servers[0].inject_phenix_env is True
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_server_inject_phenix_env_parsed_false():
  tmp = tempfile.mkdtemp()
  try:
    _write_profile(tmp, "p", {"name": "p", "model": "m",
      "mcp_servers": [{"name": "ext", "command": "x", "inject_phenix_env": False}]})
    prof = ProfileLoader(builtin_dir=Path(tmp), log=null_out()).load("p")
    assert prof.mcp_servers[0].inject_phenix_env is False
  finally:
    shutil.rmtree(tmp)


def exercise_system_prompt_file_outside_profile_dir_rejected():
  """system_prompt_file is inlined into the system prompt, so a project-
  supplied (untrusted) profile must not point it at an arbitrary file via
  an absolute path or a '..' escape (which would exfiltrate e.g.
  ~/.ssh/id_rsa). Both must raise Sorry."""
  tmp = tempfile.mkdtemp()
  secret_dir = tempfile.mkdtemp()
  try:
    secret = os.path.join(secret_dir, "secret.txt")
    with open(secret, "w") as fh:
      fh.write("PRIVATE KEY MATERIAL")
    _write_profile(tmp, "abs", {
      "name": "abs", "model": "m", "system_prompt_file": secret})
    loader = ProfileLoader(builtin_dir=Path(tmp), log=null_out())
    try:
      loader.load("abs")
    except Sorry as e:
      assert "system_prompt_file" in str(e), str(e)
    else:
      raise Exception_expected
    rel = os.path.relpath(secret, tmp)
    _write_profile(tmp, "rel", {
      "name": "rel", "model": "m", "system_prompt_file": rel})
    try:
      loader.load("rel")
    except Sorry as e:
      assert "system_prompt_file" in str(e), str(e)
    else:
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)
    shutil.rmtree(secret_dir)


def exercise():
  exercise_minimal_profile_loads()
  exercise_missing_required_field_raises()
  exercise_based_on_inheritance()
  exercise_based_on_child_system_prompt_file_overrides_parent_prompt()
  exercise_based_on_child_system_prompt_overrides_parent_prompt_file()
  exercise_based_on_does_not_weaken_mutual_exclusivity_gate()
  exercise_based_on_cycle_detected()
  exercise_load_file_from_explicit_path()
  exercise_load_file_based_on_sibling_wins_over_builtin()
  exercise_load_file_does_not_mutate_loader_search_dirs()
  exercise_load_file_based_on_falls_through_to_builtin()
  exercise_load_file_missing_raises()
  exercise_variable_expansion_env()
  exercise_project_overrides_user_overrides_builtin()
  exercise_unknown_fields_warn_not_fail()
  exercise_system_prompt_file_with_expansion()
  exercise_system_prompt_file_non_ascii_utf8_round_trips()
  exercise_profile_json_non_ascii_utf8_loads()
  exercise_profile_backend_defaults_to_claude_code()
  exercise_profile_backend_claude_code_parsed()
  exercise_openai_is_known_backend()
  exercise_portkey_is_known_backend()
  exercise_google_is_known_backend()
  exercise_portkey_block_parsed()
  exercise_portkey_block_absent_defaults_none()
  exercise_profile_backend_unknown_raises()
  exercise_server_tools_parsed_and_defaults_empty()
  exercise_server_tools_string_raises_clear_sorry()
  exercise_skills_scalar_additional_disabled_coerced_to_single_id()
  exercise_mcp_server_env_cannot_override_path_or_phenix_vars()
  exercise_mcp_server_inject_phenix_env_default_true()
  exercise_mcp_server_inject_phenix_env_parsed_false()
  exercise_system_prompt_file_outside_profile_dir_rejected()
  exercise_backend_display_name_maps_each_backend()


def exercise_backend_display_name_maps_each_backend():
  """backend_display_name turns a backend id into a user-facing assistant name
  (Claude / GPT / Gemini / …); unknown or empty falls back to 'Assistant'."""
  from qttbx.widgets.chat.agent.profile import backend_display_name
  assert backend_display_name("claude_code") == "Claude"
  assert backend_display_name("anthropic") == "Claude"
  assert backend_display_name("openai") == "GPT"
  assert backend_display_name("google") == "Gemini"
  assert backend_display_name("portkey") == "Assistant"
  assert backend_display_name("") == "Assistant"
  assert backend_display_name(None) == "Assistant"


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
