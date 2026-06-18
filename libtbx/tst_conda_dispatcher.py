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


def exercise_generated_dispatcher():
  """The conda-package dispatcher and the development-build dispatcher must
  both emit the activation block."""
  env = libtbx.env
  # Isolate from any dispatcher_include*.sh present in the build directory.
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
    activate_d = os.path.join("etc", "conda", "activate.d")

    def emits_activation(writer, name):
      target_file = os.path.join(bin_dir, name)
      if os.name == "nt":
        target_file += ".bat"
      writer(
        source_file=env.as_relocatable_path(source_file),
        target_file=env.as_relocatable_path(target_file))
      with open(target_file) as f:
        text = f.read()
      return (activate_d in text, text)

    # The conda-package dispatcher always emits the block.
    ok, text = emits_activation(env.write_conda_dispatcher, "tst_conda_disp")
    assert ok, text

    # The development-build dispatcher emits it when use_conda is set.
    env.build_options.use_conda = True
    dev_writer = (env.write_win32_dispatcher if os.name == "nt"
                  else env.write_bin_sh_dispatcher)
    ok, text = emits_activation(dev_writer, "tst_dev_disp")
    assert ok, text
  finally:
    env.build_options.use_conda = saved_use_conda
    shutil.rmtree(tmp, ignore_errors=True)


def exercise():
  exercise_runtime()
  exercise_generated_dispatcher()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
