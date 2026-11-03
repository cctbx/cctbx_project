"""Verify the environment guard the dispatchers emit before running Python.

The guard clears the variables through which a caller's environment can
inject code or libraries into the shipped Python: PYTHONPATH, the loader
search variables, GUI-toolkit module paths (GTK, GIO, GDK-pixbuf,
fontconfig, Qt). It also disables the per-user site-packages directory
(PYTHONNOUSERSITE). Setting either CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES or
PHENIX_TRUST_OTHER_ENV to a non-empty value opts out of the whole guard.
PYTHONHOME is cleared outside the guard, opt-out or not. LD_PRELOAD is
deliberately kept: VirtualGL's vglrun works by setting it. So are MPLBACKEND
and MATPLOTLIBRC (headless users set MPLBACKEND=Agg, and a clean sys.path
already stops a foreign backend loading) and GTK_IM_MODULE (it only picks
an input method among modules found through paths the guard clears).

The DivCon integration variables (QB_PYTHONPATH, QB_LD_LIBRARY_PATH,
QB_DYLD_LIBRARY_PATH) are appended after the guard so they survive it.

Every dispatcher gets the guard: the conda-package dispatchers of an
installation and the development-build dispatchers of every builder, with or
without conda. Only the development-build libtbx.scons is exempt, because it
builds against the caller's compiler environment; the installed one is
guarded like every other installed dispatcher.
"""

import os
import shutil
import subprocess
import sys
import tempfile

import libtbx.load_env
from libtbx import env_config
from libtbx import dispatcher_env
from libtbx.utils import format_cpu_times

IS_NT = (os.name == "nt")

PROBE = ["PYTHONPATH", "PYTHONHOME", "PYTHONNOUSERSITE", "LD_LIBRARY_PATH",
         "LD_PRELOAD", "GTK_MODULES", "QT_PLUGIN_PATH", "MPLBACKEND",
         "MATPLOTLIBRC", "CONDA_DLL_SEARCH_MODIFICATION_ENABLE",
         "GTK_EXE_PREFIX", "GTK_DATA_PREFIX", "GI_TYPELIB_PATH",
         "GTK_IM_MODULE", "SSL_CERT_FILE"]

POLLUTION = {
  "PYTHONPATH": "/pp",
  "PYTHONHOME": "/ph",
  "LD_LIBRARY_PATH": "/ll",
  "LD_PRELOAD": "/usr/lib/libvglfaker.so",
  "GTK_MODULES": "canberra-gtk-module",
  "GTK_IM_MODULE": "ibus",
  "QT_PLUGIN_PATH": "/usr/lib/qt/plugins",
  "MPLBACKEND": "Agg",
  "MATPLOTLIBRC": "/home/user/.config/matplotlib",
  "GTK_EXE_PREFIX": "/usr",
  "GTK_DATA_PREFIX": "/usr",
  "GI_TYPELIB_PATH": "/usr/lib/girepository-1.0",
}


def write_file(path, text):
  with open(path, "w") as f:
    f.write(text)


def parse_output(text):
  result = {}
  for line in text.splitlines():
    if "=" in line:
      key, _, value = line.partition("=")
      result[key.strip()] = value.strip()
  return result


def parse_env_dump(text):
  """The environment a child process received, from ``env`` output. Names
  it did not receive read "UNSET", so an exported empty value ("") stays
  distinguishable from an unset one, and a variable the script set without
  exporting it reads "UNSET" too."""
  result = {}
  for line in text.splitlines():
    key, sep, value = line.partition("=")
    if sep:
      result[key] = value
  for v in PROBE:
    result.setdefault(v, "UNSET")
  return result


def make_harness(tmp, shell):
  body = (env_config.environment_guard_lines(shell)
          + env_config.external_path_lines(shell))
  if shell == "bat":
    # cmd has no export step: set affects child processes directly, and
    # set "X=" removes X, so echoing in the same shell is faithful.
    header = ["@echo off", 'set "LIBTBX_PREFIX=%s"' % tmp]
    footer = ["if defined %s (echo %s=%%%s%%) else (echo %s=UNSET)"
              % (v, v, v, v) for v in PROBE]
    script = os.path.join(tmp, "harness.bat")
  else:
    # Read the environment from a child process, not the script's own
    # shell variables.
    header = ["#!/bin/sh", 'LIBTBX_PREFIX="%s"' % tmp, "export LIBTBX_PREFIX"]
    footer = ["env"]
    script = os.path.join(tmp, "harness.sh")
  write_file(script, "\n".join(header + body + footer) + "\n")
  if shell != "bat":
    os.chmod(script, 0o755)
  return script


def run_harness(script, extra_env):
  child_env = os.environ.copy()
  for v in PROBE + list(dispatcher_env.dispatcher_environment_opt_outs) \
      + ["QB_PYTHONPATH", "QB_LD_LIBRARY_PATH", "QB_DYLD_LIBRARY_PATH"]:
    child_env.pop(v, None)
  child_env.update(POLLUTION)
  child_env.update(extra_env)
  cmd = ["cmd", "/c", script] if IS_NT else [script]
  p = subprocess.run(cmd, env=child_env, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, universal_newlines=True)
  assert p.returncode == 0, (p.returncode, p.stdout, p.stderr)
  return parse_output(p.stdout) if IS_NT else parse_env_dump(p.stdout)


def exercise_guard_runtime():
  """Run the guard in a shell and check what the child process sees."""
  shell = "bat" if IS_NT else "sh"
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    script = make_harness(tmp, shell)

    # Default: the variables this platform's guard clears are gone (not
    # exported empty), the rest (e.g. LD_* and GTK_* on Windows) are
    # untouched, and user site-packages is off.
    cleared = (dispatcher_env.dispatcher_cleared_variables_win32 if IS_NT
               else dispatcher_env.dispatcher_cleared_variables)
    out = run_harness(script, {})
    for v, value in POLLUTION.items():
      if v in cleared:
        assert out[v] == "UNSET", (v, out)
      else:
        assert out[v] == value, (v, out)
    # LD_PRELOAD survives on every platform; vglrun depends on it.
    assert out["LD_PRELOAD"] == POLLUTION["LD_PRELOAD"], out
    # So do the caller's matplotlib settings, e.g. MPLBACKEND=Agg on a
    # headless node.
    for v in ("MPLBACKEND", "MATPLOTLIBRC"):
      assert out[v] == POLLUTION[v], (v, out)
    # And the input-method choice (ibus for CJK input): it only selects
    # among modules found through paths the guard already cleared.
    assert out["GTK_IM_MODULE"] == POLLUTION["GTK_IM_MODULE"], out
    if not IS_NT:
      # The bundled GTK 4 reads the GTK prefixes and GObject introspection
      # reads GI_TYPELIB_PATH; another GTK installation's values must go.
      for v in ("GTK_EXE_PREFIX", "GTK_DATA_PREFIX", "GI_TYPELIB_PATH"):
        assert out[v] == "UNSET", (v, out)
    assert out["PYTHONNOUSERSITE"] == "1", out
    # The guard sets nothing about DLL search on any platform.
    assert out["CONDA_DLL_SEARCH_MODIFICATION_ENABLE"] == "UNSET", out

    # Either opt-out name disables the guard entirely.
    for opt_out in dispatcher_env.dispatcher_environment_opt_outs:
      out = run_harness(script, {opt_out: "1"})
      for v, value in POLLUTION.items():
        assert out[v] == value, (opt_out, v, out)
      assert out["PYTHONNOUSERSITE"] == "UNSET", (opt_out, out)

    if not IS_NT:
      # DivCon paths are appended after the guard, without a leading
      # separator when the guard emptied the variable.
      out = run_harness(script, {"QB_PYTHONPATH": "/qb",
                                 "QB_LD_LIBRARY_PATH": "/qbl"})
      assert out["PYTHONPATH"] == "/qb", out
      assert out["LD_LIBRARY_PATH"] == "/qbl", out
      out = run_harness(script, {"QB_PYTHONPATH": "/qb",
                                 "PHENIX_TRUST_OTHER_ENV": "1"})
      assert out["PYTHONPATH"] == "/pp" + os.pathsep + "/qb", out
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_generated_dispatchers():
  """The conda-package dispatcher and the development-build include both
  carry the guard, ahead of the conda activation block."""
  env = libtbx.env
  env._dispatcher_include_at_start = []
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    bin_dir = os.path.join(tmp, "bin")
    os.makedirs(bin_dir)
    source_file = os.path.join(tmp, "tst_guard_src.py")
    write_file(source_file, "print('hello')\n")
    target_file = os.path.join(bin_dir, "tst_guard_disp")
    if IS_NT:
      target_file += ".bat"
    env.write_conda_dispatcher(
      source_file=env.as_relocatable_path(source_file),
      target_file=env.as_relocatable_path(target_file))
    with open(target_file) as f:
      text = f.read()
    for opt_out in dispatcher_env.dispatcher_environment_opt_outs:
      assert opt_out in text, (opt_out, text)
    assert "PYTHONNOUSERSITE=1" in text, text
    assert text.index("PYTHONNOUSERSITE=1") < text.index("activate.d"), text
    for v in ("MPLBACKEND", "MATPLOTLIBRC"):
      assert v not in text, (v, text)
    assert "CONDA_DLL_SEARCH_MODIFICATION_ENABLE" not in text, text
    if IS_NT:
      assert_cleared_before_guard(text, '@set "PYTHONHOME="',
                                  "@if not defined")
    if not IS_NT:
      assert "LD_PRELOAD" not in text, text
      assert "QB_PYTHONPATH" in text, text

    if not IS_NT:
      # development-build include (write_gui_dispatcher_include --use_conda)
      from libtbx.auto_build import write_gui_dispatcher_include
      build_dir = os.path.join(tmp, "build")
      base_dir = os.path.join(tmp, "base")
      os.makedirs(os.path.join(base_dir, "lib"))
      os.makedirs(build_dir)
      write_gui_dispatcher_include.run([
        "--build_dir=%s" % build_dir,
        "--base_dir=%s" % base_dir,
        "--suffix=tst",
        "--use_conda",
        "--ignore_missing_dirs",
        "--quiet",
      ])
      with open(os.path.join(build_dir, "dispatcher_include_tst.sh")) as f:
        text = f.read()
      # The guard comes from the dispatcher writers now, so the include
      # keeps only the DivCon paths; the BOOST_ADAPTBX_*_DEFAULT exports
      # went with it (nothing has read them since 2019).
      assert "PYTHONNOUSERSITE" not in text, text
      assert "libtbx.scons" not in text, text
      assert "# include before command" in text, text
      assert "QB_PYTHONPATH" in text, text
      assert "BOOST_ADAPTBX" not in text, text
      assert "include at start" not in text, text

      # The non-conda include (legacy base/ builds) must not carry its own
      # copy of the guard either: it would run again, after the writers'
      # guard and after the site includes, and it ignored one opt-out.
      if sys.platform == "darwin":
        os.makedirs(os.path.join(base_dir, "Python.framework", "Versions",
                                 "Current", "lib"))
      write_gui_dispatcher_include.run([
        "--build_dir=%s" % build_dir,
        "--base_dir=%s" % base_dir,
        "--suffix=tst_base",
        "--ignore_missing_dirs",
        "--quiet",
      ])
      with open(os.path.join(build_dir, "dispatcher_include_tst_base.sh")) as f:
        text = f.read()
      assert "include at start" not in text, text
      assert "PHENIX_TRUST_OTHER_ENV" not in text, text
      assert "BOOST_ADAPTBX" not in text, text
      assert "QB_PYTHONPATH" in text, text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_windows_guard_text():
  """The Windows guard is generated as plain text, so check it on every
  platform: it disables user site-packages and leaves DLL search alone."""
  text = "\n".join(dispatcher_env.environment_guard_lines("bat"))
  assert '@set "PYTHONNOUSERSITE=1"' in text, text
  assert "CONDA_DLL_SEARCH_MODIFICATION_ENABLE" not in text, text


def exercise_windows_dev_dispatcher_text():
  """The Windows development dispatcher is plain text too. Its own path
  variables must not end in ';%VAR%' once the guard has undefined VAR: cmd
  expands an undefined variable to nothing, leaving a trailing ';', and
  CPython turns that empty PYTHONPATH entry into the working directory,
  ahead of the standard library. The ';' is prepended only when the
  variable is defined."""
  import re
  env = libtbx.env
  env._dispatcher_include_at_start = []
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    source_file = os.path.join(tmp, "tst_win_dev_src.py")
    write_file(source_file, "print('hello')\n")
    target_file = os.path.join(tmp, "tst_win_dev.bat")
    env.write_win32_dispatcher(
      source_file=env.as_relocatable_path(source_file),
      target_file=env.as_relocatable_path(target_file))
    with open(target_file) as f:
      lines = f.read().splitlines()
    assert not [l for l in lines if re.match(r"^@set \w+=.*;%\w+%$", l)], lines
    own = [(i, m) for i, m in
           ((i, re.match(r'^@set "(\w+)=(.*)%(\w+)%"$', l))
            for i, l in enumerate(lines)) if m and m.group(1) == m.group(3)]
    assert {m.group(1) for i, m in own} >= {"PYTHONPATH", "PATH"}, lines
    for i, m in own:
      var, paths = m.group(1), m.group(2)
      assert not paths.endswith(";"), lines[i]
      assert lines[i - 1] == '@if defined %s @set "%s=;%%%s%%"' % (var, var, var), \
        (lines[i - 1], lines[i])
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_windows_conda_dispatcher_path():
  """The conda-package .bat lives in <prefix>\\Library\\bin, so LIBTBX_PREFIX
  is <prefix>\\Library and the PATH prepend must be relative to that: the
  prefix root (python.exe), Library\\mingw-w64\\bin, Library\\usr\\bin,
  Library\\bin (the dispatchers and DLLs), then Scripts and bin under the
  prefix. Generated as text on every platform."""
  env = libtbx.env
  env._dispatcher_include_at_start = []
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  tmp = os.path.realpath(tempfile.mkdtemp())
  saved_platform = sys.platform
  try:
    source_file = os.path.join(tmp, "tst_win_conda_src.py")
    write_file(source_file, "print('hello')\n")
    target_file = os.path.join(tmp, "tst_win_conda.bat")
    sys.platform = "win32"
    try:
      env.write_conda_dispatcher(
        source_file=env.as_relocatable_path(source_file),
        target_file=env.as_relocatable_path(target_file))
    finally:
      sys.platform = saved_platform
    with open(target_file) as f:
      lines = f.read().splitlines()
    expected = ("@set PATH=%LIBTBX_PREFIX%\\..;%LIBTBX_PREFIX%\\mingw-w64\\bin;"
                "%LIBTBX_PREFIX%\\usr\\bin;%LIBTBX_PREFIX%\\bin;"
                "%LIBTBX_PREFIX%\\..\\Scripts;%LIBTBX_PREFIX%\\..\\bin;%PATH%")
    assert expected in lines, [l for l in lines if l.startswith("@set PATH")]
    # the installer's phenix_bin copy is relocated from the writer's real
    # output, so the two cannot drift apart unnoticed
    from libtbx.auto_build.conda_build import create_custom_bin
    relocated = create_custom_bin.relocate_bat_dispatcher(lines)
    assert relocated != lines, "the LIBTBX_PREFIX block was not recognized"
    assert '@for %%F in ("%~dp0..\\Library") do @set LIBTBX_PREFIX=%%~fF' \
      in relocated, relocated[:8]
    assert not [l for l in relocated if l.startswith("@set LIBTBX_PREFIX=")], \
      relocated[:8]
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_bootstrap_configure_environment():
  """bootstrap removes PYTHONHOME, and unless the caller opted out with
  either variable the dispatchers honor also PYTHONPATH and the loader
  paths, from the environment handed to the steps that run the conda
  python (configure.py, the include generator, pip); the bootstrap
  process's own environment is left alone."""
  import importlib.util
  path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "auto_build", "bootstrap.py")
  spec = importlib.util.spec_from_file_location("tst_bootstrap", path)
  bootstrap = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(bootstrap)
  saved = os.environ.get("PYTHONPATH")
  os.environ["PYTHONPATH"] = os.path.join(os.sep, "site", "modules")
  try:
    cmd = bootstrap.ShellCommand(command=["true"], env={"PYTHONPATH": None})
    env = cmd.get_environment()
    assert "PYTHONPATH" not in env, env.get("PYTHONPATH")
    assert os.environ["PYTHONPATH"] == os.path.join(os.sep, "site", "modules")
  finally:
    if saved is None:
      del os.environ["PYTHONPATH"]
    else:
      os.environ["PYTHONPATH"] = saved
  # the clearing covers the module and library paths and PYTHONHOME;
  # either opt-out name, set to a non-empty value, keeps the paths but
  # PYTHONHOME goes regardless, as in the dispatchers
  env = bootstrap.conda_configure_environment({})
  for name in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH",
               "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH"):
    assert name in env and env[name] is None, (name, env)
  for name in ("CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES",
               "PHENIX_TRUST_OTHER_ENV"):
    env = bootstrap.conda_configure_environment({name: "1"})
    assert env == {"PYTHONHOME": None}, (name, env)
    assert "PYTHONPATH" in bootstrap.conda_configure_environment({name: ""})


def exercise_dev_dispatchers():
  """Every development-build dispatcher carries the guard, with or without
  conda, ahead of the site include files and of its own PYTHONPATH (so it
  neither undoes a site's settings nor erases the build's module paths).
  libtbx.scons is left without it."""
  env = libtbx.env
  env._dispatcher_include_at_start = ["CCTBX_TST_SITE_INCLUDE=1"]
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  saved_use_conda = env.build_options.use_conda
  writer = env.write_win32_dispatcher if IS_NT else env.write_bin_sh_dispatcher
  own_pythonpath = ('@set "PYTHONPATH=%LIBTBX_BUILD%' if IS_NT
                    else 'PYTHONPATH="$LIBTBX_BUILD')
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    bin_dir = os.path.join(tmp, "bin")
    os.makedirs(bin_dir)
    source_file = os.path.join(tmp, "tst_dev_guard_src.py")
    write_file(source_file, "print('hello')\n")

    def generate(name):
      target_file = os.path.join(bin_dir, name + (".bat" if IS_NT else ""))
      writer(source_file=env.as_relocatable_path(source_file),
             target_file=env.as_relocatable_path(target_file))
      with open(target_file) as f:
        return f.read()

    for use_conda in (False, True):
      env.build_options.use_conda = use_conda
      text = generate("tst_dev_guard")
      for opt_out in dispatcher_env.dispatcher_environment_opt_outs:
        assert opt_out in text, (use_conda, opt_out, text)
      guard = text.index("PYTHONNOUSERSITE=1")
      assert guard < text.index("CCTBX_TST_SITE_INCLUDE"), (use_conda, text)
      assert guard < text.index(own_pythonpath), (use_conda, text)
      if use_conda:
        assert guard < text.index("activate.d"), text
      text = generate("libtbx.scons")
      assert "PYTHONNOUSERSITE" not in text, (use_conda, text)
  finally:
    env.build_options.use_conda = saved_use_conda
    env._dispatcher_include_at_start = []
    shutil.rmtree(tmp, ignore_errors=True)


def assert_cleared_before_guard(text, line, guard_start):
  """``line`` appears exactly once, and before the guard, so an opt-out
  cannot skip it."""
  assert text.count(line) == 1, (line, text)
  assert text.index(line) < text.index(guard_start), (line, text)


def run_preamble(dispatcher_text, bin_dir, extra_env):
  """Run a generated POSIX dispatcher up to the point where it would start
  Python, then dump the environment a child process receives."""
  preamble = dispatcher_text[:dispatcher_text.index("\nLIBTBX_PYEXE=")]
  script = os.path.join(bin_dir, "tst_preamble")
  write_file(script, preamble + "\nenv\n")
  os.chmod(script, 0o755)
  child_env = os.environ.copy()
  for v in PROBE + list(dispatcher_env.dispatcher_environment_opt_outs):
    child_env.pop(v, None)
  child_env.update(extra_env)
  p = subprocess.run([script], env=child_env, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, universal_newlines=True)
  assert p.returncode == 0, (p.returncode, p.stdout, p.stderr)
  return parse_env_dump(p.stdout)


def exercise_interpreter_setup_ignores_opt_outs():
  """PYTHONHOME is always cleared and the dispatcher's own bin directory is
  always first on PATH, even when a caller opts out of the guard: a stray
  PYTHONHOME stops the interpreter from starting, and tools that ship only in
  that bin directory (xtb for qrefine) are found through PATH. Checked by
  running the installed and development POSIX dispatchers up to the exec,
  and on the Windows development dispatcher's text. The development
  dispatcher also points SSL_CERT_FILE at certifi's bundle, but only when
  the caller has not set one (a corporate CA bundle must win); the
  installed dispatcher leaves it alone."""
  env = libtbx.env
  env._dispatcher_include_at_start = []
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  saved_use_conda = env.build_options.use_conda
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    bin_dir = os.path.join(tmp, "bin")
    os.makedirs(bin_dir)
    source_file = os.path.join(tmp, "tst_setup_src.py")
    write_file(source_file, "print('hello')\n")

    def generate(writer, name):
      target_file = os.path.join(bin_dir, name)
      writer(source_file=env.as_relocatable_path(source_file),
             target_file=env.as_relocatable_path(target_file))
      with open(target_file) as f:
        return f.read()

    # The Windows development dispatcher is plain text on every platform.
    text = generate(env.write_win32_dispatcher, "tst_setup_dev.bat")
    assert_cleared_before_guard(text, '@set "PYTHONHOME="', "@if not defined")

    if not IS_NT:
      env.build_options.use_conda = True
      for writer, name in ((env.write_conda_dispatcher, "tst_setup_installed"),
                           (env.write_bin_sh_dispatcher, "tst_setup_dev")):
        text = generate(writer, name)
        for opt_out in [None] + list(dispatcher_env.dispatcher_environment_opt_outs):
          extra = {"PYTHONHOME": "/ph", "PATH": "/usr/bin:/bin"}
          if opt_out is not None:
            extra[opt_out] = "1"
          out = run_preamble(text, bin_dir, extra)
          assert out["PYTHONHOME"] == "UNSET", (name, opt_out, out)
          assert out["PATH"].split(":")[0] == bin_dir, (name, opt_out, out)
          if name == "tst_setup_dev":
            assert out["SSL_CERT_FILE"].endswith("cacert.pem"), (opt_out, out)
            extra["SSL_CERT_FILE"] = "/corp/ca-bundle.pem"
            out = run_preamble(text, bin_dir, extra)
            assert out["SSL_CERT_FILE"] == "/corp/ca-bundle.pem", (opt_out, out)
          else:
            assert out["SSL_CERT_FILE"] == "UNSET", (opt_out, out)
  finally:
    env.build_options.use_conda = saved_use_conda
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_opt_out_comment_names_every_opt_out():
  """The guard's comment is built from the opt-out names, so adding one
  neither breaks the guard nor leaves it undocumented."""
  saved = dispatcher_env.dispatcher_environment_opt_outs
  try:
    dispatcher_env.dispatcher_environment_opt_outs = saved + (
      "CCTBX_TST_THIRD_OPT_OUT",)
    for shell, prefix in (("sh", "# "), ("bat", "@rem ")):
      lines = dispatcher_env.environment_guard_lines(shell)
      comment = " ".join(l[len(prefix):] for l in lines if l.startswith(prefix))
      for name in dispatcher_env.dispatcher_environment_opt_outs:
        assert name in comment, (shell, name, lines)
  finally:
    dispatcher_env.dispatcher_environment_opt_outs = saved


def exercise_include_generator_ignores_foreign_libtbx():
  """bootstrap runs the include generator with whatever Python it has; a
  different libtbx importable there (e.g. a cctbx-base pulled in by another
  package) must neither break it nor be used in place of this checkout's."""
  if IS_NT:
    return
  script = os.path.join(os.path.dirname(os.path.abspath(dispatcher_env.__file__)),
                        "auto_build", "write_gui_dispatcher_include.py")
  for foreign_has_module in (False, True):
    tmp = os.path.realpath(tempfile.mkdtemp())
    try:
      foreign = os.path.join(tmp, "foreign", "libtbx")
      os.makedirs(foreign)
      write_file(os.path.join(foreign, "__init__.py"), "")
      if foreign_has_module:
        write_file(os.path.join(foreign, "dispatcher_env.py"),
                   "def external_path_lines(shell):\n"
                   "  return ['CCTBX_TST_FOREIGN_LIBTBX']\n")
      build_dir = os.path.join(tmp, "build")
      base_dir = os.path.join(tmp, "base")
      os.makedirs(os.path.join(base_dir, "lib"))
      os.makedirs(build_dir)
      child_env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"),
                   "HOME": os.environ.get("HOME", tmp),
                   "PYTHONPATH": os.path.join(tmp, "foreign")}
      p = subprocess.run(
        [sys.executable, script, "--build_dir=%s" % build_dir,
         "--base_dir=%s" % base_dir, "--suffix=tst", "--use_conda",
         "--ignore_missing_dirs", "--quiet"],
        env=child_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
      assert p.returncode == 0, (foreign_has_module, p.stdout, p.stderr)
      with open(os.path.join(build_dir, "dispatcher_include_tst.sh")) as f:
        text = f.read()
      assert "QB_PYTHONPATH" in text, text
      assert "CCTBX_TST_FOREIGN_LIBTBX" not in text, text
    finally:
      shutil.rmtree(tmp, ignore_errors=True)


def exercise():
  exercise_guard_runtime()
  exercise_generated_dispatchers()
  exercise_windows_guard_text()
  exercise_windows_dev_dispatcher_text()
  exercise_windows_conda_dispatcher_path()
  exercise_bootstrap_configure_environment()
  exercise_dev_dispatchers()
  exercise_interpreter_setup_ignores_opt_outs()
  exercise_opt_out_comment_names_every_opt_out()
  exercise_include_generator_ignores_foreign_libtbx()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
