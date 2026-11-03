"""Verify that the dispatchers run the conda activation scripts
(etc/conda/activate.d) without requiring "conda activate".

The activation scripts set environment variables (e.g. GSETTINGS_SCHEMA_DIR,
XML_CATALOG_FILES) that some packages need at runtime. The dispatcher sources
them, setting CONDA_PREFIX only while they run and then restoring it, so only
the variables the scripts export persist. Sourcing is skipped if this
environment is already active.

Both the conda-package dispatchers (installed distributions) and the
development-build dispatchers (write_bin_sh_dispatcher / write_win32_dispatcher)
emit this block; the setpaths scripts no longer do.
"""

import contextlib
import os
import shutil
import subprocess
import tempfile

import libtbx.load_env
from libtbx import env_config
from libtbx.utils import format_cpu_times


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


def exercise_runtime():
  """Wrap the activation block in a minimal dispatcher fragment and run it,
  checking the activation scripts are sourced and CONDA_PREFIX restored.

  Both forms the codebase emits are exercised: the dispatcher form that
  resolves the prefix from LIBTBX_PREFIX at runtime (conda-package
  dispatchers), and the baked-literal form that the
  development-build dispatchers emit. With LIBTBX_PREFIX pointed at the same
  prefix, both must behave identically."""
  is_nt = (os.name == "nt")
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    # On Windows the conda-package dispatchers live under <prefix>/Library/bin
    # so LIBTBX_PREFIX is <prefix>/Library and the prefix is LIBTBX_PREFIX/.. ;
    # elsewhere LIBTBX_PREFIX is the prefix itself.
    conda_prefix = tmp
    if is_nt:
      libtbx_prefix = os.path.join(tmp, "Library")
      os.makedirs(libtbx_prefix)
      activate_name = "zzz_sentinel.bat"
      activate_text = '@set "CCTBX_TST_ACTIVATE_SENTINEL=%CONDA_PREFIX%\\marker"\n'
      shell = "bat"
      script_name = "harness.bat"
      header = [
        "@echo off",
        'set "LIBTBX_PREFIX=%s"' % libtbx_prefix,
      ]
      footer = [
        "if defined CCTBX_TST_ACTIVATE_SENTINEL "
        "(echo SENTINEL=%CCTBX_TST_ACTIVATE_SENTINEL%) else (echo SENTINEL=MISSING)",
        "if defined CONDA_PREFIX "
        "(echo CONDA_PREFIX=%CONDA_PREFIX%) else (echo CONDA_PREFIX=UNSET)",
      ]
    else:
      libtbx_prefix = tmp
      activate_name = "zzz_sentinel.sh"
      activate_text = 'export CCTBX_TST_ACTIVATE_SENTINEL="${CONDA_PREFIX}/marker"\n'
      shell = "sh"
      script_name = "harness.sh"
      header = [
        "#!/bin/sh",
        'LIBTBX_PREFIX="%s"' % libtbx_prefix,
        "export LIBTBX_PREFIX",
      ]
      footer = [
        'echo "SENTINEL=${CCTBX_TST_ACTIVATE_SENTINEL:-MISSING}"',
        'echo "CONDA_PREFIX=${CONDA_PREFIX:-UNSET}"',
      ]

    activate_dir = os.path.join(conda_prefix, "etc", "conda", "activate.d")
    os.makedirs(activate_dir)
    write_file(os.path.join(activate_dir, activate_name), activate_text)

    expected_marker = os.path.join(conda_prefix, "marker")
    other = os.path.join(tmp, "other_env")
    script = os.path.join(tmp, script_name)

    def same_path(a, b):
      return os.path.normcase(a) == os.path.normcase(b)

    def run(conda_prefix_value):
      child_env = os.environ.copy()
      child_env.pop("CCTBX_TST_ACTIVATE_SENTINEL", None)
      if conda_prefix_value is None:
        child_env.pop("CONDA_PREFIX", None)
      else:
        child_env["CONDA_PREFIX"] = conda_prefix_value
      cmd = ["cmd", "/c", script] if is_nt else [script]
      p = subprocess.run(cmd, env=child_env, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
      assert p.returncode == 0, (p.returncode, p.stdout, p.stderr)
      return parse_output(p.stdout)

    def check(activation_lines):
      write_file(script, "\n".join(header + activation_lines + footer) + "\n")
      if not is_nt:
        os.chmod(script, 0o755)

      # 1) Not active: the activation scripts run (CONDA_PREFIX is the prefix
      #    while they run), then CONDA_PREFIX is restored to unset.
      out = run(conda_prefix_value=None)
      assert same_path(out["SENTINEL"], expected_marker), out
      assert out["CONDA_PREFIX"] == "UNSET", out

      # 2) Already active: the activation scripts are skipped and CONDA_PREFIX
      #    is left untouched.
      out = run(conda_prefix_value=conda_prefix)
      assert out["SENTINEL"] == "MISSING", out
      assert same_path(out["CONDA_PREFIX"], conda_prefix), out

      # 3) A different environment is active: the activation scripts run, then
      #    CONDA_PREFIX is restored to the original (different) value.
      out = run(conda_prefix_value=other)
      assert same_path(out["SENTINEL"], expected_marker), out
      assert same_path(out["CONDA_PREFIX"], other), out

    # conda-package dispatcher form: prefix resolved from LIBTBX_PREFIX.
    check(env_config.conda_activation_lines(shell))
    # development-build dispatcher form: baked-literal prefix.
    check(env_config.conda_activation_lines(
      shell, conda_prefix=conda_prefix))
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


@contextlib.contextmanager
def isolated_dispatcher_env():
  """Yield (env, bin_dir, source_file) for generating dispatchers into a
  temporary directory, isolated from any dispatcher_include*.sh in the build
  directory. The use_conda build option is restored on exit."""
  env = libtbx.env
  env._dispatcher_include_at_start = []
  env._dispatcher_include_before_command = []
  env._dispatcher_precall_commands = []
  saved_use_conda = env.build_options.use_conda
  tmp = os.path.realpath(tempfile.mkdtemp())
  try:
    bin_dir = os.path.join(tmp, "bin")
    os.makedirs(bin_dir)
    source_file = os.path.join(tmp, "tst_conda_src.py")
    write_file(source_file, "print('hello')\n")
    yield env, bin_dir, source_file
  finally:
    env.build_options.use_conda = saved_use_conda
    shutil.rmtree(tmp, ignore_errors=True)


def write_dispatcher(env, writer, bin_dir, source_file, name):
  """Generate a dispatcher with ``writer`` and return its text."""
  target_file = os.path.join(bin_dir, name)
  if os.name == "nt":
    target_file += ".bat"
  writer(
    source_file=env.as_relocatable_path(source_file),
    target_file=env.as_relocatable_path(target_file))
  with open(target_file) as f:
    return f.read()


def exercise_generated_dispatcher():
  """The conda-package dispatcher and the development-build dispatcher must
  both emit the activation block."""
  activate_d = os.path.join("etc", "conda", "activate.d")
  with isolated_dispatcher_env() as (env, bin_dir, source_file):
    # The conda-package dispatcher always emits the block.
    text = write_dispatcher(
      env, env.write_conda_dispatcher, bin_dir, source_file, "tst_conda_disp")
    assert activate_d in text, text

    # The development-build dispatcher emits it when use_conda is set.
    env.build_options.use_conda = True
    dev_writer = (env.write_win32_dispatcher if os.name == "nt"
                  else env.write_bin_sh_dispatcher)
    text = write_dispatcher(
      env, dev_writer, bin_dir, source_file, "tst_dev_disp")
    assert activate_d in text, text


def exercise_conda_bin_on_path():
  """The sh development-build dispatcher must put <conda_prefix>/bin on PATH
  directly after $LIBTBX_BUILD/bin when use_conda is set, matching the
  conda-package dispatchers (whose prefix bin is the dispatcher bin), and
  must not add it otherwise."""
  if os.name == "nt":
    return  # write_win32_dispatcher already adds the Library/bin dirs.
  with isolated_dispatcher_env() as (env, bin_dir, source_file):
    conda_bin = (env.as_relocatable_path(env_config.get_conda_prefix())
                 / "bin").sh_value()

    def path_assignment(name):
      text = write_dispatcher(
        env, env.write_bin_sh_dispatcher, bin_dir, source_file, name)
      lines = [l.strip() for l in text.splitlines()
               if l.strip().startswith('PATH="')]
      assert len(lines) == 2, lines  # if/else branches of the essential
      return lines

    env.build_options.use_conda = True
    lines = path_assignment("tst_dev_conda")
    expected = 'PATH="$LIBTBX_BUILD/bin:%s' % conda_bin
    for line in lines:
      assert line.startswith(expected), (line, expected)

    env.build_options.use_conda = False
    lines = path_assignment("tst_dev_noconda")
    for line in lines:
      assert conda_bin not in line, (line, conda_bin)
      assert line.startswith('PATH="$LIBTBX_BUILD/bin'), line


def exercise():
  exercise_runtime()
  exercise_generated_dispatcher()
  exercise_conda_bin_on_path()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
