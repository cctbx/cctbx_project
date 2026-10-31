"""write_env_files.py writes the phenix_env.{sh,csh,bat} files of an
installer. The Windows file must put the prefix (without the trailing
backslash of %~dp0) in PHENIX and PHENIX_PREFIX, and the default bin
directory is Library\\bin on Windows, bin elsewhere."""
import os
import shutil
import tempfile

from libtbx.auto_build.conda_build import write_env_files

def read(path):
  with open(path) as f:
    return f.read()

def exercise_default_bin_dirs():
  prefix = os.path.join(os.sep, "p")
  assert write_env_files.default_bin_dirs(prefix, "win32") == \
    [os.path.join(prefix, "Library", "bin")]
  for platform in ("linux", "darwin"):
    assert write_env_files.default_bin_dirs(prefix, platform) == \
      [os.path.join(prefix, "bin")], platform

def exercise_windows_file():
  tmp = tempfile.mkdtemp()
  try:
    write_env_files.write_files(
      program="phenix", prefix=tmp, bin_dirs=[os.path.join(tmp, "phenix_bin")],
      version="2.3-1", destination=tmp, platform="win32")
    lines = read(os.path.join(tmp, "phenix_env.bat")).splitlines()
    assert "set PHENIX=%~dp0" in lines, lines
    assert "set PHENIX=%PHENIX:~0,-1%" in lines, lines
    assert "set PHENIX_PREFIX=%PHENIX%" in lines, lines
    assert "set PHENIX_VERSION=2.3-1" in lines, lines
    assert "set PATH=%s;%%PATH%%" % os.path.join(tmp, "phenix_bin") in lines, lines
    assert not os.path.exists(os.path.join(tmp, "phenix_env.sh"))
  finally:
    shutil.rmtree(tmp, ignore_errors=True)

def exercise_posix_files():
  tmp = tempfile.mkdtemp()
  try:
    write_env_files.write_files(
      program="phenix", prefix=tmp, bin_dirs=[os.path.join(tmp, "phenix_bin")],
      version="2.3-1", destination=tmp, platform="linux")
    sh = read(os.path.join(tmp, "phenix_env.sh")).splitlines()
    assert "export PHENIX=%s" % tmp in sh, sh
    assert "export PHENIX_PREFIX=%s" % tmp in sh, sh
    assert "export PATH=%s:$PATH" % os.path.join(tmp, "phenix_bin") in sh, sh
    csh = read(os.path.join(tmp, "phenix_env.csh")).splitlines()
    assert "setenv PHENIX_PREFIX %s" % tmp in csh, csh
    assert not os.path.exists(os.path.join(tmp, "phenix_env.bat"))
  finally:
    shutil.rmtree(tmp, ignore_errors=True)

def exercise_run_joins_relative_bin_dir():
  tmp = tempfile.mkdtemp()
  try:
    write_env_files.run(["--program", "phenix", "--prefix", tmp,
                         "--bin-dir", "phenix_bin", "--version", "2.3-1",
                         "--destination", tmp])
    name = "phenix_env.bat" if os.name == "nt" else "phenix_env.sh"
    text = read(os.path.join(tmp, name))
    assert os.path.join(tmp, "phenix_bin") in text, text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)

def exercise_run_with_a_relative_prefix():
  # the default bin directory is derived from the absolute prefix, so a
  # relative --prefix does not end up as <prefix>/<prefix>/bin
  tmp = tempfile.mkdtemp()
  cwd = os.getcwd()
  try:
    os.chdir(os.path.dirname(tmp))
    rel = os.path.basename(tmp)
    write_env_files.run(["--program", "phenix", "--prefix", rel,
                         "--version", "2.3-1", "--destination", tmp])
    name = "phenix_env.bat" if os.name == "nt" else "phenix_env.sh"
    text = read(os.path.join(tmp, name))
    expected = write_env_files.default_bin_dirs(os.path.abspath(rel))[0]
    assert expected in text, (expected, text)
    assert os.path.join(rel, rel) not in text, text
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
  exercise_default_bin_dirs()
  exercise_windows_file()
  exercise_posix_files()
  exercise_run_joins_relative_bin_dir()
  exercise_run_with_a_relative_prefix()
  print("OK")
