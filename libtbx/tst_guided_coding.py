from __future__ import absolute_import, division, print_function
"""Shared-suite test of the central GuidedCoding source, libtbx/guided_coding.

Runs the package's own unit tests (skill registration, source verification,
Claude Code version gate, evidence and bundle checks) in a child process with
a private TMPDIR, then checks that the source still verifies, so a test run
cannot leave stray files in the package. The package and its tools are for
macOS and Linux only: on Windows, or under Python 2, this test prints a skip
line and OK.
"""
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


def exercise():
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
    rc, out = run([sys.executable, "-I", "-B",
                   os.path.join("payload", "tools", "screen_check.py"),
                   "verify-source", "."], env)
    print(out.strip())
    assert rc == 0 and "VERIFIED complete source" in out, \
      "the package source no longer verifies after its tests"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


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
