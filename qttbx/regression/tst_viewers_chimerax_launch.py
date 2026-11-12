"""How ChimeraXViewer.start_viewer starts ChimeraX.

On macOS 26 and later ChimeraX runs as a launchd job, in a coalition of its
own (libtbx.launchd_process): as a child of a program run from a terminal it
would share the terminal's coalition, and quitting it could make macOS
SIGTERM the terminal's processes. Elsewhere it stays a Popen child."""

from __future__ import absolute_import, division, print_function

import os
import sys
import tempfile
from contextlib import redirect_stdout
from io import StringIO

from libtbx import launchd_process
from libtbx.test_utils.launchd import (
  kill_and_unload, launchd_tests_can_run, resource_coalition)
from qttbx.viewers import chimerax


# A stand-in for ChimeraX: answers every GET with 200 on the port named in
# its "--cmd 'remotecontrol rest start port N'" argument. It runs this
# Python, not /usr/bin/python3, which on a Mac without the Command Line Tools
# opens their install dialog.
FAKE_CHIMERAX = r'''#!/bin/sh
exec "%s" -c '
import http.server, sys
port = int(sys.argv[-1].split()[-1])
class H(http.server.BaseHTTPRequestHandler):
  def do_GET(self):
    self.send_response(200); self.end_headers(); self.wfile.write(b"ok")
  def log_message(self, *a):
    pass
http.server.HTTPServer(("127.0.0.1", port), H).serve_forever()
' "$@"
'''


def _fake_chimerax(tmp):
  path = os.path.join(tmp, "ChimeraX")
  with open(path, "w") as f:
    f.write(FAKE_CHIMERAX % sys.executable)
  os.chmod(path, 0o755)
  return path


def exercise_dispatch():
  """start_viewer uses launchd_process.launch exactly when
  launchd_process.use_launchd() says so, with ChimeraX's command line and
  its cleaned environment; otherwise Popen. When launchd refuses the job (no
  GUI login session), it says so and falls back to Popen. Both are faked,
  and the REST check reports ready at once, so nothing starts."""
  import subprocess
  import types
  def refuse(cmd, **kw):
    raise launchd_process.LaunchdRefused("launchctl bootstrap exited with "
                                         "status 125")
  with tempfile.TemporaryDirectory() as tmp:
    command = _fake_chimerax(tmp)
    for use, launch in ((False, None), (True, None), (True, refuse)):
      popened, launched, out = [], [], StringIO()
      saved = (subprocess.Popen, launchd_process.use_launchd,
               launchd_process.launch, chimerax.ChimeraXViewer._check_status)
      subprocess.Popen = lambda cmd, **kw: popened.append((cmd, kw)) or \
        types.SimpleNamespace(pid=1)
      launchd_process.use_launchd = lambda: use
      launchd_process.launch = launch or (lambda cmd, **kw: launched.append(
        (cmd, kw)) or types.SimpleNamespace(pid=2))
      def ready(viewer, *args, **kwargs):
        viewer._connected = True
      chimerax.ChimeraXViewer._check_status = ready
      saved_pythonpath = os.environ.get("PYTHONPATH")
      os.environ["PYTHONPATH"] = "/dev/null/not-for-chimerax"
      try:
        viewer = chimerax.ChimeraXViewer()
        viewer.command, viewer.port = command, 50123
        with redirect_stdout(out):
          viewer.start_viewer(timeout=1)
      finally:
        (subprocess.Popen, launchd_process.use_launchd, launchd_process.launch,
         chimerax.ChimeraXViewer._check_status) = saved
        if saved_pythonpath is None:
          del os.environ["PYTHONPATH"]
        else:
          os.environ["PYTHONPATH"] = saved_pythonpath
      want_cmd = [command, "--cmd", "remotecontrol rest start port 50123"]
      job = use and launch is None
      ((cmd, kw),) = launched if job else popened
      assert (popened if job else launched) == [], (use, launch)
      assert cmd == want_cmd, cmd
      assert "PYTHONPATH" not in kw["env"], "ChimeraX got PYTHONPATH"
      assert viewer.process.pid == (2 if job else 1), (use, launch)
      if launch is refuse:
        assert "status 125" in out.getvalue() and \
          "child process" in out.getvalue(), out.getvalue()


def exercise_launch_failure_is_sorry():
  """A ChimeraX that cannot be started is reported as a Sorry (a message for
  the user, raised from the launch error), not an internal error: a launchd
  job refused inside a GUI session, a ChimeraX that exits as soon as launchd
  starts it, or a missing program. A bug, such as a RecursionError, is not
  dressed up as one. launch is faked; nothing starts."""
  from libtbx.utils import Sorry
  def refused(cmd, **kw):
    raise launchd_process.LaunchError(
      "launchctl bootstrap exited with status 5: Bootstrap failed: 5: "
      "Input/output error")
  def bug(cmd, **kw):
    raise RecursionError("maximum recursion depth exceeded")
  def exited(cmd, **kw):
    raise launchd_process.ExitedAtOnce("ChimeraX", 1)
  def missing(cmd, **kw):
    raise FileNotFoundError(2, "No such file or directory", cmd[0])
  with tempfile.TemporaryDirectory() as tmp:
    command = _fake_chimerax(tmp)
    for launch, want in ((refused, "Input/output error"),
                         (exited, "exited as soon as it started"),
                         (missing, "No such file or directory"),
                         (bug, None)):
      saved = (launchd_process.use_launchd, launchd_process.launch)
      launchd_process.use_launchd = lambda: True
      launchd_process.launch = launch
      try:
        viewer = chimerax.ChimeraXViewer()
        viewer.command, viewer.port = command, 50123
        with redirect_stdout(StringIO()):
          viewer.start_viewer(timeout=1)
      except Exception as exc:
        error = exc
      else:
        error = None
      finally:
        launchd_process.use_launchd, launchd_process.launch = saved
      if want is None:
        assert isinstance(error, RecursionError), repr(error)
        continue
      assert isinstance(error, Sorry), repr(error)
      assert "ChimeraX" in str(error) and want in str(error), str(error)
      assert error.__cause__ is not None, "not raised from the launch error"


def exercise_own_coalition():
  """macOS 26+: the started ChimeraX (here the fake) answers on its REST
  port from a coalition of its own, and close_viewer's communicate() works
  on its handle."""
  if not launchd_tests_can_run():
    return
  with tempfile.TemporaryDirectory() as tmp:
    viewer = chimerax.ChimeraXViewer()
    viewer.command = _fake_chimerax(tmp)
    try:
      with redirect_stdout(StringIO()):
        viewer.start_viewer(timeout=20)
      assert viewer._connected
      assert resource_coalition(viewer.process.pid) != \
        resource_coalition(os.getpid()), "ChimeraX in the caller's coalition"
    finally:
      kill_and_unload(viewer.process)
    assert viewer.process.communicate(timeout=10) == (None, None)


def exercise():
  exercise_dispatch()
  exercise_launch_failure_is_sorry()
  exercise_own_coalition()


if __name__ == "__main__":
  exercise()
  print("OK")
