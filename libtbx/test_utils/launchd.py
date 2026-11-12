"""Test helpers for programs started through libtbx.launchd_process (macOS 26
and later), shared by the cctbx and Phenix tests of its callers."""

from __future__ import absolute_import, division, print_function

import os

from libtbx import launchd_process as lp

# Label prefix of the probe jobs session_runs_live_tests loads.
PROBE_PREFIX = "org.cctbx.test_utils.launchd_probe"

_can_run = []

def launchd_tests_can_run():
  """True on hosts that start programs through launchd (macOS 26 and later)
  where this session loads launchd jobs; see session_runs_live_tests, which
  is asked once."""
  if not _can_run:
    _can_run.append(lp.use_launchd() and session_runs_live_tests())
  return _can_run[0]


def session_runs_live_tests():
  """True if this login session loads a launchd job, checked by loading a
  probe job that never runs and unloading it again (probe jobs left by a
  killed run are unloaded first).

  False, saying that the live tests are skipped, where the session has no
  GUI domain (an SSH login with nobody logged in at the screen, headless
  CI) or launchctl does not answer. Raises AssertionError where the GUI
  domain exists but the probe is refused: Coot and ChimeraX could not start
  in that session either, so that is a failure, not a skip."""
  import plistlib
  import subprocess
  import tempfile
  domain = "gui/%d" % os.getuid()
  label = lp.job_label(PROBE_PREFIX, os.getpid())
  try:
    lp.unload_stale_jobs(PROBE_PREFIX)
    if lp.launchctl("print", domain)[0] != 0:
      print("launchd live tests skipped: this session has no GUI domain (%s)"
            % domain)
      return False
    fd, path = tempfile.mkstemp(suffix=".plist")
    try:
      with os.fdopen(fd, "wb") as f:
        plistlib.dump({"Label": label, "ProgramArguments": ["/usr/bin/true"],
                       "RunAtLoad": False, "KeepAlive": False}, f)
      status, _, err = lp.launchctl("bootstrap", domain, path)
    finally:
      os.remove(path)
  except subprocess.TimeoutExpired:
    print("launchd live tests skipped: launchctl timed out")
    return False
  if status != 0:
    lines = err.strip().splitlines()
    raise AssertionError(
      "this GUI session refuses launchd jobs (launchctl bootstrap exited with "
      "status %d: %s), so Coot and ChimeraX cannot start in it"
      % (status, lines[0] if lines else ""))
  try:
    lp.launchctl("bootout", lp.job_target(label))
  except subprocess.TimeoutExpired:
    pass                               # the next probe's sweep unloads it
  return True


def resource_coalition(pid):
  """macOS: the resource coalition id of ``pid``, read with
  ``proc_pidinfo(PROC_PIDCOALITIONINFO)`` (no root needed for our own uid)."""
  import ctypes
  lib = ctypes.CDLL("/usr/lib/libproc.dylib")
  lib.proc_pidinfo.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
                               ctypes.c_void_p, ctypes.c_int]
  info = (ctypes.c_uint64 * 5)()       # coalition_id[2] + 3 reserved
  n = lib.proc_pidinfo(pid, 20, 0, info, ctypes.sizeof(info))
  assert n == ctypes.sizeof(info), ("proc_pidinfo", pid, n)
  return info[0]


def job_loaded(label):
  """True if launchd still has a job named ``label``."""
  return lp.launchctl("list", label)[0] == 0


def waiter_blocked(proc):
  """True while a thread holds LaunchdProcess ``proc``'s lock: one blocked
  in ``wait()`` holds it until the program exits (a ``poll()`` recording
  the exit also holds it, for a moment)."""
  return proc._lock.locked()


def kill_and_unload(proc):
  """Test cleanup: SIGKILL a launched program's process group if it is still
  running (so the pid is still its own), then unload its job."""
  import signal
  if proc is None:
    return
  if proc.poll() is None:
    try:
      os.killpg(proc.pid, signal.SIGKILL)
    except OSError:
      pass
  label = getattr(proc, "label", None)
  if label:
    lp.launchctl("bootout", lp.job_target(label))
