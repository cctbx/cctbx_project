from __future__ import absolute_import, division, print_function
"""Shared-suite test of the central GuidedCoding source, libtbx/guided_coding.

Runs the package's own unit tests (skill registration, source verification,
Claude Code version gate, evidence and bundle checks) in a child process with
a private TMPDIR. Afterwards it fails if the run created any file in the
package, or if any file listed in SOURCE_MANIFEST.sha256 is missing or
changed. Unlisted files that were already present before the run (for
example bytecode written by an installer's precompile step) are reported but
do not fail this test. The package and its tools are for macOS and Linux
only: on Windows, or under Python 2, this test prints a skip line and OK.
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


def run(args, env):
  p = subprocess.Popen(args, cwd=PACKAGE, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  out = p.communicate()[0].decode("utf-8", "replace")
  return p.returncode, out


def package_files():
  result = set()
  for dirpath, dirnames, filenames in os.walk(PACKAGE):
    for name in filenames:
      path = os.path.join(dirpath, name)
      result.add(os.path.relpath(path, PACKAGE).replace(os.sep, "/"))
  return result


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


def exercise():
  listed = set(
    line.split("  ", 1)[1].strip()[2:]
    for line in open(os.path.join(PACKAGE, "SOURCE_MANIFEST.sha256")))
  listed.add("SOURCE_MANIFEST.sha256")
  before = package_files()
  extra = sorted(before - listed)
  if extra:
    print("note: unlisted files already present in the package: %s"
          % ", ".join(extra))
  # The tools refuse paths through symbolic links; macOS /var -> /private/var.
  tmp = os.path.realpath(tempfile.mkdtemp(prefix="tst_guided_coding_"))
  try:
    env = dict(os.environ)
    env.update({"TMPDIR": tmp, "PYTHONDONTWRITEBYTECODE": "1"})
    rc, out = run([sys.executable, "-B", "-m", "unittest", "discover",
                   "-s", "tests", "-p", "tst_*.py", "-v"], env)
    print(out)
    assert rc == 0, "package unit tests failed (exit %d)" % rc
    ran = re.search(r"^Ran (\d+) tests? in ", out, re.M)
    assert ran is not None and int(ran.group(1)) > 0, "no package tests ran"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)
  created = sorted(package_files() - before)
  assert not created, "the test run created files in the package: %s" % \
    ", ".join(created)
  problems = listed_problems()
  assert not problems, "listed source files: %s" % "; ".join(problems)
  print("listed source files unchanged; no files created in the package")


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
