from __future__ import absolute_import, division, print_function
"""Shared-suite test of the central GuidedCoding source, libtbx/guided_coding.

Copies the files listed in SOURCE_MANIFEST.sha256 to a private temporary
directory and runs the package's own unit tests there (skill registration,
source verification, Claude Code version gate, evidence and bundle checks),
so that bytecode or other unlisted files in the installed package are never
imported by this test. It fails if the tests fail, if they create files in
the copy, or if any listed file in the package is missing or changed;
unlisted files already present in the package are reported. It also runs the
installer's precompile step, libtbx.py_compile_all -i, on another copy and
checks that the copy still verifies while an unlisted module is still
refused. The package and its tools are for macOS and Linux only: on Windows,
or under Python 2, this test prints a skip line and OK.
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE = os.path.join(HERE, "guided_coding")


def run(args, env, cwd):
  p = subprocess.Popen(args, cwd=cwd, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  out = p.communicate()[0].decode("utf-8", "replace")
  return p.returncode, out


def files_under(root):
  result = set()
  for dirpath, dirnames, filenames in os.walk(root):
    for name in filenames:
      path = os.path.join(dirpath, name)
      result.add(os.path.relpath(path, root).replace(os.sep, "/"))
  return result


def listed_names():
  names = set(
    line.split("  ", 1)[1].strip()[2:]
    for line in open(os.path.join(PACKAGE, "SOURCE_MANIFEST.sha256")))
  names.add("SOURCE_MANIFEST.sha256")
  return names


def copy_listed(dest):
  for rel in sorted(listed_names()):
    target = os.path.join(dest, rel)
    if not os.path.isdir(os.path.dirname(target)):
      os.makedirs(os.path.dirname(target))
    shutil.copyfile(os.path.join(PACKAGE, rel), target)
  return dest


def listed_problems():
  problems = []
  with open(os.path.join(PACKAGE, "SOURCE_MANIFEST.sha256")) as f:
    for line in f:
      digest, rel = line.rstrip("\n").split("  ", 1)
      path = os.path.join(PACKAGE, rel[2:] if rel.startswith("./") else rel)
      if not os.path.isfile(path):
        problems.append("missing " + rel)
        continue
      with open(path, "rb") as g:
        if hashlib.sha256(g.read()).hexdigest() != digest:
          problems.append("changed " + rel)
  return problems


def exercise_package_tests(tmp, env):
  copy = copy_listed(os.path.join(tmp, "tested", "guided_coding"))
  before = files_under(copy)
  rc, out = run([sys.executable, "-B", "-m", "unittest", "discover",
                 "-s", "tests", "-p", "tst_*.py", "-v"], env, copy)
  print(out)
  assert rc == 0, "package unit tests failed (exit %d)" % rc
  ran = re.search(r"^Ran (\d+) tests? in ", out, re.M)
  assert ran is not None and int(ran.group(1)) > 0, "no package tests ran"
  created = sorted(files_under(copy) - before)
  assert not created, "the package tests created files in the package: %s" % \
    ", ".join(created)


def exercise_real_precompile(tmp, env):
  """Run the installer's precompile step, libtbx.py_compile_all -i, on a copy
  of the listed files: the copy must still verify (reporting the bytecode it
  ignored), and an unlisted module must still be refused."""
  rc, out = run([sys.executable, "-B", "-c",
                 "import importlib.util; importlib.util.find_spec("
                 "'libtbx.command_line.py_compile_all').name"], env, tmp)
  if rc != 0:
    print("skip: libtbx.py_compile_all needs a libtbx environment")
    return
  if getattr(sys, "pycache_prefix", None):
    print("skip: this Python writes bytecode outside __pycache__ (pycache_prefix)")
    return
  copy = copy_listed(os.path.join(tmp, "precompiled", "guided_coding"))
  rc, out = run([sys.executable, "-m", "libtbx.command_line.py_compile_all",
                 "-i", copy], env, copy)
  assert rc == 0, out
  compiled = [name for dirpath, dirnames, names in os.walk(copy)
              for name in names if name.endswith(".pyc")]
  assert compiled, "libtbx.py_compile_all wrote no bytecode"
  checker = os.path.join(copy, "payload", "tools", "screen_check.py")
  rc, out = run([sys.executable, "-I", "-B", checker, "verify-source", copy],
                env, copy)
  print(out.strip())
  assert rc == 0 and "VERIFIED complete source" in out, out
  assert "ignored %d precompiled" % len(compiled) in out, out
  with open(os.path.join(copy, "payload", "tools", "extra.py"), "w") as f:
    f.write("print('unlisted')\n")
  rc, out = run([sys.executable, "-I", "-B", checker, "verify-source", copy],
                env, copy)
  assert rc == 2 and "extra source file" in out, out
  print("precompiled copy verifies; an unlisted module is still refused")


def exercise():
  extra = sorted(files_under(PACKAGE) - listed_names())
  if extra:
    print("note: unlisted files present in the package (not used by this test): %s"
          % ", ".join(extra))
  # The tools refuse paths through symbolic links; macOS /var -> /private/var.
  tmp = os.path.realpath(tempfile.mkdtemp(prefix="tst_guided_coding_"))
  try:
    env = dict(os.environ)
    env.pop("PYTHONPYCACHEPREFIX", None)
    env.update({"TMPDIR": tmp, "PYTHONDONTWRITEBYTECODE": "1"})
    exercise_package_tests(tmp, env)
    exercise_real_precompile(tmp, env)
  finally:
    shutil.rmtree(tmp, ignore_errors=True)
  problems = listed_problems()
  assert not problems, "listed source files: %s" % "; ".join(problems)
  print("listed source files unchanged")


def run_all():
  if sys.platform == "win32":
    print("skip: GuidedCoding tools are for macOS and Linux only")
  elif sys.version_info[0] < 3:
    print("skip: GuidedCoding tools need Python 3")
  else:
    exercise()
  print("OK")


if __name__ == "__main__":
  run_all()
