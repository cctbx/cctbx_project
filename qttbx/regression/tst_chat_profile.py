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
  exercise_based_on_cycle_detected()
  exercise_variable_expansion_env()
  exercise_project_overrides_user_overrides_builtin()
  exercise_unknown_fields_warn_not_fail()
  exercise_system_prompt_file_with_expansion()
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
