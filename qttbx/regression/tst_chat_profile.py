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


def exercise():
  exercise_minimal_profile_loads()
  exercise_missing_required_field_raises()
  exercise_based_on_inheritance()
  exercise_based_on_cycle_detected()
  exercise_variable_expansion_env()
  exercise_project_overrides_user_overrides_builtin()
  exercise_unknown_fields_warn_not_fail()
  exercise_system_prompt_file_with_expansion()
  exercise_profile_backend_defaults_to_claude_code()
  exercise_profile_backend_claude_code_parsed()
  exercise_profile_backend_unknown_raises()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
