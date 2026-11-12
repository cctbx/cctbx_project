"""Tests for libtbx.launchd_process.

The live tests start real launchd jobs that run /bin/sh and /bin/sleep (no
GUI program, so killing them can never trigger the macOS application-death
cleanup the module exists to avoid). They run only where launchd jobs are
used (macOS 26 and later) and this session loads launchd jobs. The tests
with a faked launchctl run on any macOS, and those that need neither
launchctl nor kqueue on every platform."""

from __future__ import absolute_import, division, print_function

import os
import sys
import tempfile
import time

from libtbx import launchd_process as lp
from libtbx.test_utils.launchd import (
  job_loaded, kill_and_unload, launchd_tests_can_run, resource_coalition)

TEST_PREFIX = "org.cctbx.tst_launchd_process"


def own_jobs(label_prefix=TEST_PREFIX):
  """Labels of the jobs with ``label_prefix`` this process has started."""
  mine = "%s.%d." % (label_prefix, os.getpid())
  out = lp.launchctl("list")[1]
  return [line.split("\t")[-1] for line in out.splitlines()
          if line.split("\t")[-1].startswith(mine)]


# A value as Python holds an undecodable environment byte: a surrogate
# escape. Built by hand, because os.fsdecode raises on Windows.
_UNDECODABLE = b"caf\xe9".decode("utf-8", "surrogateescape")


def _sh(body):
  return ["/bin/sh", "-c", body]


def _wait_exit(proc, seconds=10):
  deadline = time.monotonic() + seconds
  while proc.poll() is None and time.monotonic() < deadline:
    time.sleep(0.1)
  return proc.poll()


class _signals_only_to(object):
  """Context manager: ``os.kill`` sends a signal only to the pids in
  ``allowed``; any other raises AssertionError instead of being sent, so a
  test of stop()'s guards can never signal the terminal's processes. Signal
  0 (an existence check) goes through, except that every call for a pid in
  ``foreign`` raises PermissionError, as for another user's process."""
  def __init__(self, allowed, foreign=()):
    self.allowed, self.foreign = allowed, foreign
  def __enter__(self):
    self.saved = os.kill
    def kill(pid, sig):
      if pid in self.foreign:
        raise PermissionError(1, "Operation not permitted")
      if sig != 0 and pid not in self.allowed:
        raise AssertionError("signal %d to pid %d, not one of %s"
                             % (sig, pid, sorted(self.allowed)))
      return self.saved(pid, sig)
    os.kill = kill
    return self
  def __exit__(self, *exc):
    os.kill = self.saved


def _wait_for_text(path, text, seconds=10):
  deadline = time.monotonic() + seconds
  content = ""
  while text not in content and time.monotonic() < deadline:
    time.sleep(0.1)
    content = open(path).read() if os.path.exists(path) else ""
  return content


# --- tests -------------------------------------------------------------------

def exercise_use_launchd():
  """macOS 26 (Darwin 25) and later use launchd; older macOS and other
  platforms do not. The Darwin release decides, not platform.mac_ver(),
  which a Python built against an old SDK reads as 10.16; an unreadable
  release counts as new. Runs on every platform (os.uname is faked, and
  restored even where it does not exist)."""
  import platform
  import types
  saved = (sys.platform, platform.mac_ver, getattr(os, "uname", None))
  try:
    for plat, darwin, mac_ver, want in (
        ("darwin", "24.6.0", "15.7.1", False),
        ("darwin", "25.0.0", "26.0", True),
        ("darwin", "27.2.0", "10.16", True),
        ("darwin", "", "", True),
        ("linux", "6.8.0", "", False),
        ("win32", "", "", False)):
      sys.platform = plat
      platform.mac_ver = lambda *a, **k: (mac_ver, ("", "", ""), "arm64")
      os.uname = lambda: types.SimpleNamespace(release=darwin)
      assert lp.use_launchd() is want, (plat, darwin, mac_ver)
  finally:
    sys.platform, platform.mac_ver = saved[0], saved[1]
    if saved[2] is None:
      del os.uname
    else:
      os.uname = saved[2]


def exercise_own_coalition():
  """A launched program does not share the caller's resource coalition, so
  its exit cannot make macOS SIGTERM the caller's terminal processes."""
  if not launchd_tests_can_run():
    return
  proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
  try:
    # Already the leader of its own process group, so a killpg(pid) at once
    # reaches it (while xpcproxy is still setting it up, its group is 1).
    assert os.getpgid(proc.pid) == proc.pid, os.getpgid(proc.pid)
    assert resource_coalition(proc.pid) != resource_coalition(os.getpid())
  finally:
    kill_and_unload(proc)


def exercise_lifecycle():
  """poll() is None while the program runs. wait(), blocked in another
  thread, returns once it exits, even while poll() is called concurrently;
  poll() then reports the exit. The finished job is unloaded before the
  returncode is set, so a caller that sees the exit never finds the job
  still loaded (the unload is slowed down here to show that)."""
  if not launchd_tests_can_run():
    return
  import signal
  import threading
  bootout = lp._bootout_quietly
  def slow_bootout(label):
    time.sleep(0.5)
    bootout(label)
  proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
  lp._bootout_quietly = slow_bootout
  try:
    time.sleep(0.5)
    assert proc.poll() is None, "a running program must poll() as None"
    waited = []
    waiter = threading.Thread(target=lambda: waited.append(proc.wait()),
                              daemon=True)
    waiter.start()
    time.sleep(0.2)                     # wait() is now blocked
    os.killpg(proc.pid, signal.SIGTERM)
    assert _wait_exit(proc) is not None, "poll() never saw the exit"
    assert not job_loaded(proc.label), "exit seen with the job still loaded"
    waiter.join(5)
    assert not waiter.is_alive(), "wait() did not return after the exit"
    assert waited == [-15], ("SIGTERM must give returncode -15", waited)
  finally:
    lp._bootout_quietly = bootout
    kill_and_unload(proc)


def exercise_job_environment():
  """The program runs in the caller's working directory with the given
  environment, and its stdout and stderr are appended to log_path, which
  may be a pathlib.Path."""
  if not launchd_tests_can_run():
    return
  from pathlib import Path
  with tempfile.TemporaryDirectory() as tmp:
    workdir = Path(tmp) / "project"
    workdir.mkdir()
    log_path = Path(tmp) / "out.log"
    env = dict(os.environ, TST_LAUNCHD_VALUE="given")
    saved_cwd = os.getcwd()
    os.chdir(str(workdir))
    proc = None
    try:
      proc = lp.launch(
        _sh('echo "cwd=$(pwd -P) value=$TST_LAUNCHD_VALUE"\n'
            'echo "to-stderr" >&2\nexec /bin/sleep 60'),
        env=env, log_path=log_path, label_prefix=TEST_PREFIX)
      text = _wait_for_text(str(log_path), "to-stderr")
      assert "cwd=%s " % os.path.realpath(str(workdir)) in text, text
      assert "value=given\n" in text, text
      assert "to-stderr" in text, text
    finally:
      os.chdir(saved_cwd)
      kill_and_unload(proc)


def exercise_immediate_exit_leaves_no_job():
  """A program that exits as soon as it starts is reported, as a failed
  launch or as a handle whose poll() shows the exit with its exit status
  (as Popen's returncode would), and leaves no job."""
  if not launchd_tests_can_run():
    return
  proc = None
  try:
    try:
      proc = lp.launch(_sh("exit 3"), label_prefix=TEST_PREFIX)
    except lp.ExitedAtOnce:
      pass
    else:
      assert _wait_exit(proc) == 3, ("exit status lost", proc.returncode)
    assert own_jobs() == [], own_jobs()
  finally:
    kill_and_unload(proc)
    for label in own_jobs():
      lp.launchctl("bootout", lp.job_target(label))


def exercise_failures():
  """A launchd failure raises and leaves no job loaded.
  - no_gui_domain: with no GUI login session (no gui/<uid> domain: SSH,
    sudo), the refused bootstrap raises LaunchdRefused, so a caller can fall
    back to Popen.
  - refused_in_session: a bootstrap refused inside a GUI session raises a
    plain LaunchError with the first line of launchctl's error (not its
    advice to re-run as root): falling back to Popen there would bring back
    the coalition hazard.
  - domain_unknown: if launchctl cannot say whether the GUI domain exists,
    it is taken to exist: no fallback.
  - no_pid: a job that has not run (run count 0) and whose pid launchd never
    reports is unloaded on giving up.
  - exited_at_once, exited_zero_at_once: a program that ran and exited
    before its pid could be read raises ExitedAtOnce at once, naming it and
    its exit status, not after the full pid wait. That includes exit status
    0, which launchd's job table cannot tell from "not run yet" (0 before
    any run) but its run count can.
  - setup_fails: if watching the started program fails, it is stopped (its
    job unloaded) rather than left running untracked.
  - bootstrap_times_out, list_times_out: a hung launchctl raises a
    LaunchError saying it timed out; the job is booted out (it may be loaded
    and running), and a bootout that also fails does not replace the
    original error.
  - interrupted: a Ctrl-C after the bootstrap also boots the job out.
  Every failure of launchd's or launchctl's is a LaunchError, so a caller
  can tell a failed launch from a bug; running out of file descriptors
  (setup_fails) raises OSError, as from Popen. launchctl is faked; nothing
  is started. macOS only (kqueue, os.getuid)."""
  if sys.platform != "darwin":
    return
  import plistlib
  import select
  def boom_kqueue():
    raise OSError(24, "Too many open files")
  import subprocess
  domain = "gui/%d" % os.getuid()
  table = {"no_pid": "-\t0", "exited_at_once": "-\t3",
           "exited_zero_at_once": "-\t0",
           "setup_fails": "%d\t0" % os.getpid()}
  for case in ("no_gui_domain", "refused_in_session", "domain_unknown",
               "no_pid", "exited_at_once", "exited_zero_at_once",
               "setup_fails", "bootstrap_times_out", "list_times_out",
               "interrupted"):
    calls, labels = [], []
    def fake_launchctl(*args):
      calls.append(args)
      if args[0] == "bootstrap":
        with open(args[2], "rb") as f:
          labels.append(plistlib.load(f)["Label"])
        if case == "no_gui_domain":
          return 112, "", "Could not find domain for user gui: 501\n"
        if case in ("refused_in_session", "domain_unknown"):
          return 5, "", ("Bootstrap failed: 5: Input/output error\n"
                         "Try re-running the command as root for richer "
                         "errors.\n")
        if case == "bootstrap_times_out":
          raise subprocess.TimeoutExpired([lp.LAUNCHCTL, "bootstrap"], 10)
        return 0, "", ""
      if args[0] == "bootout" and case.endswith("times_out"):
        raise subprocess.TimeoutExpired([lp.LAUNCHCTL, "bootout"], 10)
      if args == ("print", domain):
        if case == "domain_unknown":
          raise subprocess.TimeoutExpired([lp.LAUNCHCTL, "print"], 10)
        return (112, "", "Could not find domain") if case == "no_gui_domain" \
          else (0, "%s = {\n}\n" % domain, "")
      if args == ("list",) and labels:   # the job table, once loaded
        if case == "list_times_out":
          raise subprocess.TimeoutExpired([lp.LAUNCHCTL, "list"], 10)
        if case == "interrupted":
          raise KeyboardInterrupt
        return 0, "PID\tStatus\tLabel\n%s\t%s\n" % (table[case], labels[0]), ""
      if args[0] == "print":
        ran = case in ("exited_at_once", "exited_zero_at_once")
        return 0, "%s = {\n\tstate = not running\n\truns = %d\n}\n" % (
          args[1], 1 if ran else 0), ""
      return 0, "", ""
    saved = (lp.launchctl, lp.PID_TIMEOUT, select.kqueue)
    lp.launchctl = fake_launchctl
    lp.PID_TIMEOUT = 0.2 if case == "no_pid" else 5.0
    if case == "setup_fails":
      select.kqueue = boom_kqueue
    error = None
    began = time.monotonic()
    try:
      lp.launch(["/usr/bin/true"], label_prefix=TEST_PREFIX, name="TstViewer")
    except BaseException as exc:
      error = exc
    finally:
      lp.launchctl, lp.PID_TIMEOUT, select.kqueue = saved
    took = time.monotonic() - began
    assert error is not None, case
    message = str(error)
    assert len(labels) == 1, (case, calls)
    booted_out = [a for a in calls if a[0] == "bootout"]
    if case not in ("interrupted", "setup_fails"):
      assert isinstance(error, lp.LaunchError), (case, repr(error))
    if case in ("no_gui_domain", "refused_in_session", "domain_unknown"):
      refused = isinstance(error, lp.LaunchdRefused)
      assert refused == (case == "no_gui_domain"), (case, repr(error))
      want = {"no_gui_domain": "Could not find domain",
              "refused_in_session": "Input/output error",
              "domain_unknown": "Input/output error"}[case]
      assert want in message, message
      # Only launchctl's first line: its advice to re-run as root would send
      # the user to sudo, where there is no GUI domain at all.
      assert "root" not in message and "\n" not in message, message
      assert booted_out == [], calls       # nothing was loaded
      continue
    assert not isinstance(error, lp.LaunchdRefused), (case, repr(error))
    if case == "interrupted":
      assert isinstance(error, KeyboardInterrupt), repr(error)
    if case.endswith("times_out"):
      # The job may be loaded and its program running: it is booted out,
      # and a failing bootout does not hide the original error.
      assert "timed out" in message and case.split("_")[0] in message, message
      assert len(booted_out) == 1, (case, calls)
      continue
    assert booted_out == [("bootout", lp.job_target(labels[0]))], (case, calls)
    if case == "no_pid":
      assert "did not report a running TstViewer" in message, message
    if case.startswith("exited"):
      assert isinstance(error, lp.ExitedAtOnce), repr(error)
      assert took < 2.0, ("waited for a pid after the exit", took)
      code = 3 if case == "exited_at_once" else 0
      assert error.returncode == code, (case, error.returncode)
      assert "TstViewer exited as soon as it started (exit status %d)" % code \
        in message, message
    if case == "setup_fails":
      assert "Too many open files" in message, message

def exercise_program_lookup():
  """The program is found as Popen finds it: a bare name on the PATH of the
  environment it is given (launchd would search its own; an empty PATH means
  the working directory), a path as given. A missing or non-executable
  program raises FileNotFoundError or PermissionError, as from Popen, before
  any job is loaded, also when the only match on PATH is not executable; a
  non-executable match earlier on PATH does not hide an executable one
  later. A found program is the job's Program, and argv is passed
  unchanged. launchctl is faked (it refuses the job after the plist is
  read). macOS only."""
  if sys.platform != "darwin":
    return
  import plistlib
  with tempfile.TemporaryDirectory() as tmp:
    bin_dir = os.path.join(tmp, "bin")
    os.mkdir(bin_dir)
    program = os.path.join(bin_dir, "tst-launchd-program")
    plain = os.path.join(bin_dir, "not-executable")
    other_dir = os.path.join(tmp, "other")
    os.mkdir(other_dir)
    shadowed = os.path.join(other_dir, "shadowed")
    for path, mode in ((program, 0o755), (plain, 0o644), (shadowed, 0o755),
                       (os.path.join(bin_dir, "shadowed"), 0o644)):
      with open(path, "w") as f:
        f.write("#!/bin/sh\n")
      os.chmod(path, mode)
    saved_cwd = os.getcwd()
    os.chdir(bin_dir)                  # where an empty PATH looks
    try:
      for argv, path, want in (
          (["tst-launchd-program", "-x"], bin_dir, program),
          ([program], bin_dir, program),
          (["tst-launchd-program"], "", os.path.join(os.getcwd(),
                                                     "tst-launchd-program")),
          (["tst-no-such-program"], bin_dir, FileNotFoundError),
          ([os.path.join(bin_dir, "missing")], bin_dir, FileNotFoundError),
          ([plain], bin_dir, PermissionError),
          (["not-executable"], bin_dir, PermissionError),
          (["shadowed"], os.pathsep.join([bin_dir, other_dir]), shadowed)):
        jobs = []
        def fake_launchctl(*args):
          if args[0] == "bootstrap":
            with open(args[2], "rb") as f:
              jobs.append(plistlib.load(f))
            return 5, "", "refused by the test"
          return 0, "", ""
        saved, lp.launchctl = lp.launchctl, fake_launchctl
        try:
          lp.launch(argv, env={"PATH": path}, label_prefix=TEST_PREFIX)
        except Exception as exc:
          error = exc
        finally:
          lp.launchctl = saved
        if isinstance(want, type):
          assert isinstance(error, want), (argv, repr(error))
          assert jobs == [], (argv, "a job was loaded")
          continue
        assert "refused by the test" in str(error), (argv, repr(error))
        (job,) = jobs
        assert job["Program"] == want, (argv, job["Program"])
        assert job["ProgramArguments"] == argv, (argv, job["ProgramArguments"])
    finally:
      os.chdir(saved_cwd)


def exercise_job_state_parsing():
  """job_state reads the job's row of launchd's job table: its pid, or None
  when it is not running, and its last exit as a Popen returncode (-N after
  signal N). Only the row labelled exactly with the label counts. job_runs
  reads the job's own run count, the last "runs" line of launchctl print,
  which comes after the job's arguments and environment: one of those can
  hold a line that reads the same (macOS only). launchctl is faked with its
  real formats."""
  table = ("PID\tStatus\tLabel\n-\t0\tcom.apple.x\n%s"
           "-\t9\ttst.longer\n")
  printed = ("gui/501/tst = {\n\tstate = not running\n\targuments = {\n"
             "\t\t/bin/sh\n\t\tx\n\truns = 7\n\t}\n\tenvironment = {\n"
             "\t\tSPOOF => v\n\truns = 8\n\t}\n\truns = 1\n}\n")
  saved = lp.launchctl
  try:
    for row, want in (("4242\t0\ttst\n", (4242, 0)),
                      ("-\t-15\ttst\n", (None, -15)),
                      ("-\t3\ttst\n", (None, 3)),
                      ("", (None, None))):
      lp.launchctl = lambda *args: (0, table % row, "")
      assert lp.job_state("tst") == want, (lp.job_state("tst"), want)
    if sys.platform == "darwin":       # job_runs names the job by uid
      lp.launchctl = lambda *args: (0, printed, "")
      assert lp.job_runs("tst") == 1, lp.job_runs("tst")
  finally:
    lp.launchctl = saved

def exercise_job_runs():
  """launchd's run count tells a job that has not run yet from one that has
  run and exited 0: the job table's last exit is 0 for both. launch()
  relies on it."""
  if not launchd_tests_can_run():
    return
  waiting = lp.job_label(TEST_PREFIX, os.getpid())
  ran = lp.job_label(TEST_PREFIX, os.getpid())
  try:
    _bootstrap_job(waiting, ["/usr/bin/true"], run_at_load=False)
    _bootstrap_job(ran, ["/usr/bin/true"])
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and (
        lp.job_runs(ran) == 0 or lp.job_state(ran)[0] is not None):
      time.sleep(0.1)
    assert lp.job_state(waiting) == lp.job_state(ran) == (None, 0), \
      (lp.job_state(waiting), lp.job_state(ran))
    assert (lp.job_runs(waiting), lp.job_runs(ran)) == (0, 1), \
      (lp.job_runs(waiting), lp.job_runs(ran))
    assert lp.job_runs(lp.job_label(TEST_PREFIX, os.getpid())) == 0
  finally:
    for label in (waiting, ran):
      lp.launchctl("bootout", lp.job_target(label))


def exercise_exit_recording():
  """Popen's behavior at the edges of recording an exit, on handles for
  plain child processes (launchctl is faked to report exit status 3):
  - wait(timeout=0) returns the returncode of a program that has already
    exited, instead of raising TimeoutExpired;
  - an interruption (Ctrl-C) while the exit is being recorded leaves a
    handle whose next poll() records it, not one whose poll() raises
    ValueError on its closed kqueue;
  - communicate(input) ignores the input, as Popen's does without a stdin
    pipe;
  - wait(timeout=inf) works, as Popen's does;
  - an exit whose kqueue event was taken but not recorded (a Ctrl-C at that
    moment) is still found: wait() rechecks that the process exists at least
    once a second, so it does not hang, even without a timeout and when the
    first check still finds the process (a zombie launchd has not reaped);
  - nor when its pid now belongs to another process, in another coalition
    than the job's; a coalition that cannot be read is not taken for that;
  - send_signal(), terminate() and kill() signal nothing once the exit is
    recorded (the pid may be another process's by then).
  macOS only (kqueue)."""
  if sys.platform != "darwin":
    return
  import subprocess
  import threading
  def fake_launchctl(*args):
    if args == ("list",):
      return 0, ("PID\tStatus\tLabel\n-\t3\t%s.fake\n-\t-15\t%s.term\n"
                 % (TEST_PREFIX, TEST_PREFIX)), ""
    return 0, "", ""
  def exited_handle(name="fake"):
    child = subprocess.Popen(["/bin/sleep", "0.2"])
    handle = lp.LaunchdProcess("%s.%s" % (TEST_PREFIX, name), child.pid)
    child.wait()
    return handle
  def interrupted():
    raise KeyboardInterrupt
  saved, lp.launchctl = lp.launchctl, fake_launchctl
  try:
    handle = exited_handle()
    assert handle.wait(timeout=0) == 3, handle.returncode
    assert exited_handle("term").wait() == -15    # signal 15, as from Popen
    handle = exited_handle()
    handle._exit_code = interrupted
    try:
      handle.poll()
    except KeyboardInterrupt:
      pass
    else:
      raise AssertionError("the interruption did not happen")
    del handle._exit_code
    assert handle.poll() == 3, handle.returncode
    # As with Popen and no stdin pipe, input is ignored.
    assert handle.communicate(input=b"ignored") == (None, None)
    assert exited_handle().wait(timeout=float("inf")) == 3
    handle = exited_handle()
    assert handle._kq.control(None, 1, 5), "no exit event"   # taken, lost
    began = time.monotonic()
    assert handle.wait(timeout=5) == 3, handle.returncode
    assert time.monotonic() - began < 2.5, "the lost exit took the timeout"
    # Lost, and the first check still finds the (zombie) process.
    handle = exited_handle()
    assert handle._kq.control(None, 1, 5), "no exit event"
    real_pid_exists, checks = lp._pid_exists, []
    def zombie_once(pid):
      checks.append(pid)
      return len(checks) == 1 or real_pid_exists(pid)
    lp._pid_exists = zombie_once
    try:
      got = []
      waiter = threading.Thread(target=lambda: got.append(handle.wait()),
                                daemon=True)
      waiter.start()
      waiter.join(5)
    finally:
      lp._pid_exists = real_pid_exists
    assert got == [3] and len(checks) >= 2, (got, checks)
    # Lost, and the pid now belongs to a process in another coalition.
    handle = exited_handle()
    assert handle._kq.control(None, 1, 5), "no exit event"
    handle.coalition = 999999999       # as for a launchd job
    saved_helpers = lp._pid_exists, lp._resource_coalition
    lp._pid_exists = lambda pid: True
    lp._resource_coalition = lambda pid: 7
    try:
      began = time.monotonic()
      assert handle.wait(timeout=5) == 3, handle.returncode
      assert time.monotonic() - began < 2.5, "the reused pid took the timeout"
    finally:
      lp._pid_exists, lp._resource_coalition = saved_helpers
    child = subprocess.Popen(["/bin/sleep", "30"])
    try:
      handle = lp.LaunchdProcess(TEST_PREFIX + ".fake", child.pid)
      handle.coalition = 999999999
      lp._resource_coalition = lambda pid: None
      assert handle.poll() is None, "an unreadable coalition read as an exit"
    finally:
      lp._resource_coalition = saved_helpers[1]
      child.kill()
      child.wait()
    # Recorded: no signal goes to the pid.
    handle = exited_handle()
    assert handle.poll() == 3
    import signal
    with _signals_only_to(set()):
      handle.send_signal(signal.SIGUSR1)
      handle.terminate()
      handle.kill()
  finally:
    lp.launchctl = saved


def exercise_concurrent_wait_timeouts():
  """Two threads waiting with a timeout on a program that keeps running
  each get TimeoutExpired close to their own deadline. The second, started
  0.3 s later, gets the handle's lock when the first gives up, and must not
  then wait a full timeout again. test_utils.launchd.waiter_blocked sees a
  thread blocked in wait(), and no thread before. A handle for a plain child
  process; nothing is loaded. macOS only (kqueue)."""
  if sys.platform != "darwin":
    return
  import subprocess
  import threading
  child = subprocess.Popen(["/bin/sleep", "30"])
  handle = lp.LaunchdProcess(TEST_PREFIX + ".fake", child.pid)
  took = []
  def waiter():
    began = time.monotonic()
    try:
      handle.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
      took.append(time.monotonic() - began)
  from libtbx.test_utils.launchd import waiter_blocked
  threads = [threading.Thread(target=waiter) for _ in range(2)]
  try:
    assert not waiter_blocked(handle), "blocked before any wait()"
    threads[0].start()
    time.sleep(0.3)
    assert waiter_blocked(handle), "a thread in wait() not seen"
    threads[1].start()
    for thread in threads:
      thread.join(10)
  finally:
    # A failed assertion must not leave a waiter to record the exit later,
    # through whatever launchctl a later test has faked.
    saved, lp.launchctl = lp.launchctl, lambda *args: (1, "", "")
    try:
      child.kill()
      child.wait()
      for thread in threads:
        if thread.is_alive():
          thread.join(10)
    finally:
      lp.launchctl = saved
  assert len(took) == 2 and max(took) < 1.25, took


def exercise_live_test_gate():
  """The live tests run only where this session loads a launchd job, which
  their gate, test_utils.launchd.launchd_tests_can_run, checks by loading a
  probe job (never run) and unloading it again; probe jobs left by an
  earlier, killed run are unloaded first. Where launchd is not used, it is
  False without asking launchctl.
  - loads: run.
  - no GUI domain (SSH with nobody at the screen, headless CI) or launchctl
    timing out: skip, saying so.
  - refused although the GUI domain exists: fail. Coot and ChimeraX could
    not start in that session either, and a skip there would be reported
    as OK.
  launchctl is faked. macOS only (os.getuid)."""
  if sys.platform != "darwin":
    return
  import contextlib
  import subprocess
  from io import StringIO
  from libtbx.test_utils import launchd as tu
  gone = subprocess.Popen(["/usr/bin/true"])
  gone.wait()                          # the pid of a run that has exited
  stale = lp.job_label(tu.PROBE_PREFIX, gone.pid)
  domain = "gui/%d" % os.getuid()
  for case, want in (("loads", True), ("no_domain", False),
                     ("times_out", False), ("refused", AssertionError)):
    calls = []
    def fake_launchctl(*args):
      calls.append(args)
      if args == ("print", domain):
        return (112, "", "") if case == "no_domain" else (0, "", "")
      if args == ("list",):
        return 0, "PID\tStatus\tLabel\n-\t0\t%s\n" % stale, ""
      if args[0] == "bootstrap":
        if case == "times_out":
          raise subprocess.TimeoutExpired([lp.LAUNCHCTL, "bootstrap"], 10)
        if case != "refused":
          return 0, "", ""
        return 5, "", ("Bootstrap failed: 5\nTry re-running the command "
                       "as root for richer errors.\n")
      return 0, "", ""
    saved = (lp.launchctl, lp.use_launchd, list(tu._can_run))
    lp.launchctl, lp.use_launchd = fake_launchctl, lambda: True
    del tu._can_run[:]                 # ask afresh, through the gate
    out = StringIO()
    try:
      with contextlib.redirect_stdout(out):
        got = tu.launchd_tests_can_run()
    except AssertionError as exc:
      got = exc
    finally:
      lp.launchctl, lp.use_launchd = saved[0], saved[1]
      tu._can_run[:] = saved[2]
    if want is AssertionError:
      assert isinstance(got, AssertionError), (case, got)
      assert "refuses" in str(got) and "root" not in str(got), str(got)
    else:
      assert got is want, (case, got, calls)
    probes = [a[2] for a in calls if a[0] == "bootstrap"]
    booted_out = [a[1] for a in calls if a[0] == "bootout"]
    assert lp.job_target(stale) in booted_out, (case, "stale probe kept")
    if case == "loads":
      assert len(probes) == 1 and len(booted_out) == 2, calls
    if case == "no_domain":
      assert probes == [], calls
    if want is False:
      assert "skipped" in out.getvalue(), (case, out.getvalue())
  calls = []
  saved = (lp.launchctl, lp.use_launchd, list(tu._can_run))
  lp.launchctl = lambda *args: calls.append(args) or (0, "", "")
  lp.use_launchd = lambda: False
  del tu._can_run[:]
  try:
    assert tu.launchd_tests_can_run() is False
  finally:
    lp.launchctl, lp.use_launchd = saved[0], saved[1]
    tu._can_run[:] = saved[2]
  assert calls == [], calls

def exercise_describe_returncode():
  """Exit codes read as Popen gives them: -N is signal N, and a Windows crash
  code reads in hex, as Windows documents it (0xC0000005 is an access
  violation). Runs on every platform."""
  for code, want in ((0, "exit status 0"), (3, "exit status 3"),
                     (-15, "signal 15"),
                     (3221225477, "exit status 0xC0000005")):
    got = lp.describe_returncode(code)
    assert got == want, (code, got)


def exercise_exited_at_once_pickles():
  """ExitedAtOnce survives pickling (as an exception crossing a process
  boundary must), keeping its message, name, returncode and any notes added
  to it. Runs on every platform."""
  import pickle
  original = lp.ExitedAtOnce("Coot", 3)
  if hasattr(original, "add_note"):    # Python 3.11 and later
    original.add_note("while restarting")
  exc = pickle.loads(pickle.dumps(original))
  assert isinstance(exc, lp.ExitedAtOnce), repr(exc)
  assert (str(exc), exc.name, exc.returncode) == (
    "Coot exited as soon as it started (exit status 3)", "Coot", 3), repr(exc)
  assert getattr(exc, "__notes__", None) == \
    getattr(original, "__notes__", None), getattr(exc, "__notes__", None)


def exercise_spawn():
  """spawn() starts the program as a launchd job where use_launchd() says
  so, and otherwise with the caller's own popen(). When launchd refuses the
  job (no GUI login session: SSH, sudo), it warns and falls back to popen(),
  as before launchd was used; any other launch failure is raised. launch,
  use_launchd and popen are faked; nothing starts. Runs on every platform."""
  saved = (lp.use_launchd, lp.launch)
  def refused(argv, **kw):
    raise lp.LaunchdRefused("launchctl bootstrap exited with status 125")
  def broken(argv, **kw):
    raise RuntimeError("TstViewer exited as soon as it started")
  try:
    for use, launch, want in ((False, None, "popen"), (True, None, "job"),
                              (True, refused, "popen"),
                              (True, broken, RuntimeError)):
      launched, popened, warnings = [], [], []
      lp.use_launchd = lambda: use
      lp.launch = launch or (lambda argv, **kw: launched.append(
        (argv, kw)) or "job")
      try:
        got = lp.spawn(["viewer", "-x"], lambda: popened.append(1) or "popen",
                       env={"A": "1"}, log_path="out.log",
                       stderr_path="err.log", label_prefix="org.tst",
                       warn=warnings.append, name="TstViewer")
      except Exception as exc:
        got = exc
      if isinstance(want, type):
        assert isinstance(got, want) and popened == [], (repr(got), popened)
        continue
      assert got == want, (use, got)
      assert len(popened) == (1 if want == "popen" else 0), (use, popened)
      if launch is refused:
        assert len(warnings) == 1 and "status 125" in warnings[0] and \
          "TstViewer" in warnings[0], warnings
      else:
        assert warnings == [], warnings
      if want == "job":
        ((argv, kw),) = launched
        assert argv == ["viewer", "-x"], argv
        assert kw == dict(env={"A": "1"}, log_path="out.log",
                          stderr_path="err.log", label_prefix="org.tst",
                          warn=warnings.append, name="TstViewer"), sorted(kw)
  finally:
    lp.use_launchd, lp.launch = saved


def exercise_signals():
  """send_signal(), terminate() and kill() signal the program as Popen's do,
  and do nothing once it has exited."""
  if not launchd_tests_can_run():
    return
  import signal
  for stop, want in ((lambda p: p.terminate(), -signal.SIGTERM),
                     (lambda p: p.kill(), -signal.SIGKILL),
                     (lambda p: p.send_signal(signal.SIGUSR1),
                      -signal.SIGUSR1)):
    proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
    try:
      stop(proc)
      assert proc.wait(timeout=10) == want, (want, proc.returncode)
      stop(proc)                       # exited: nothing to signal, no error
    finally:
      kill_and_unload(proc)


# A program that starts a copy of itself in a new session (out of its
# process group, as a self-relaunch can be) and exits, as Firefox does when
# it restarts to apply an update. The copy's pid goes to argv[1]; it ignores
# SIGTERM, so only stop()'s SIGKILL ends it.
_RELAUNCHER = """import subprocess, sys, time
copy = subprocess.Popen(["/bin/sh", "-c", "trap '' TERM; exec /bin/sleep 60"],
                        start_new_session=True)
open(sys.argv[1], "w").write(str(copy.pid))
time.sleep(float(sys.argv[2]))
"""


def exercise_stop_reaches_relaunched_copy():
  """stop() ends the program and everything it started that is still in its
  job's coalition, even a copy in a new session that outlives it under a
  new pid (as Firefox restarts itself to apply an update), which neither the
  handle's pid, its process group nor unloading the job reaches; SIGKILL
  follows for what ignores SIGTERM. job_pids() lists them, and never the
  caller."""
  if not launchd_tests_can_run():
    return
  with tempfile.TemporaryDirectory() as tmp:
    pid_file = os.path.join(tmp, "copy.pid")
    for main_lives in (False, True):
      proc, copy = None, None
      try:
        proc = lp.launch([sys.executable, "-c", _RELAUNCHER, pid_file,
                          "60" if main_lives else "0.5"],
                         label_prefix=TEST_PREFIX)
        deadline = time.monotonic() + 10
        while not os.path.exists(pid_file) and time.monotonic() < deadline:
          time.sleep(0.05)
        time.sleep(0.1)
        copy = int(open(pid_file).read())
        if not main_lives:
          assert proc.wait(timeout=10) == 0, proc.returncode
        assert os.getpgid(copy) != proc.pid, "the copy is in the job's group"
        with _signals_only_to({proc.pid, copy}):
          pids = proc.job_pids()
          assert copy in pids and os.getpid() not in pids, pids
          assert (proc.pid in pids) == main_lives, (main_lives, pids)
          proc.stop(timeout=1)
        assert not lp._pid_exists(copy), "the relaunched copy survived"
        assert proc.poll() is not None and proc.job_pids() == [], \
          proc.job_pids()
      finally:
        for pid in (copy,):
          if pid is not None and lp._pid_exists(pid):
            os.kill(pid, 9)
        kill_and_unload(proc)
        if os.path.exists(pid_file):
          os.remove(pid_file)


def exercise_stop_own_coalition_guard():
  """A handle whose program shares the caller's coalition (not a launchd
  job: here a plain child) never treats that coalition as its job's, which
  would make stop() signal the terminal's processes: job_pids() is the
  program alone, and stop() ends just it. macOS only (kqueue)."""
  if sys.platform != "darwin":
    return
  import subprocess
  saved = lp.launchctl
  label = TEST_PREFIX + ".fake"
  child = subprocess.Popen(["/bin/sleep", "60"])
  # launchd even lists it as the job's program: the coalition still is not
  # taken, being the caller's.
  lp.launchctl = lambda *args: (
    (0, "PID\tStatus\tLabel\n%d\t0\t%s\n" % (child.pid, label), "")
    if args == ("list",) else (1, "", ""))
  try:
    handle = lp.LaunchdProcess(label, child.pid)
    assert handle.coalition is None, handle.coalition
    with _signals_only_to({child.pid}):
      assert handle.job_pids() == [child.pid], handle.job_pids()
      handle.stop(timeout=5)
    assert child.wait(timeout=5) == -15, child.returncode
  finally:
    lp.launchctl = saved
    if child.poll() is None:
      child.kill()
      child.wait()


def exercise_stop_needs_the_jobs_own_program():
  """A handle takes a coalition as its job's only when launchd reports its
  pid as the running program of job ``label``: one built around a pid in
  some third coalition (here another job's program, under a label that is
  not that job's) reaches that program alone. And a coalition listing that
  holds this process or its parent is never the job's: job_pids() falls
  back to the program alone. No signal can go anywhere else."""
  if not launchd_tests_can_run():
    return
  import signal
  proc = None
  try:
    proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
    assert proc.coalition == resource_coalition(proc.pid), proc.coalition
    stranger = lp.LaunchdProcess(TEST_PREFIX + ".not.its.job", proc.pid)
    assert stranger.coalition is None, stranger.coalition
    saved = lp._coalition_pids
    with _signals_only_to({proc.pid}):
      assert stranger.job_pids() == [proc.pid], stranger.job_pids()
      try:
        for caller in (os.getpid(), os.getppid()):   # a listing gone wrong
          lp._coalition_pids = lambda coalition: sorted([proc.pid, caller])
          assert proc.job_pids() == [proc.pid], (caller, proc.job_pids())
        assert proc.stop(timeout=5) == -signal.SIGTERM, proc.returncode
      finally:
        lp._coalition_pids = saved
  finally:
    kill_and_unload(proc)


def exercise_stop_skips_foreign_processes():
  """stop() neither stops at nor waits for a process in the job's coalition
  that it may not signal (another user's, such as a setuid helper's or
  sudo's): job_pids() leaves it out, and a PermissionError while signalling
  (one that became so meanwhile) is passed over like an exit."""
  if not launchd_tests_can_run():
    return
  import signal
  import subprocess
  with tempfile.TemporaryDirectory() as tmp:
    pid_file = os.path.join(tmp, "copy.pid")
    proc, copy = None, None
    try:
      proc = lp.launch([sys.executable, "-c", _RELAUNCHER, pid_file, "60"],
                       label_prefix=TEST_PREFIX)
      deadline = time.monotonic() + 10
      while not os.path.exists(pid_file) and time.monotonic() < deadline:
        time.sleep(0.05)
      time.sleep(0.1)
      copy = int(open(pid_file).read())
      with _signals_only_to({proc.pid}, foreign={copy}):
        assert copy not in proc.job_pids(), proc.job_pids()
        began = time.monotonic()
        assert proc.stop(timeout=5) == -signal.SIGTERM, proc.returncode
        took = time.monotonic() - began
      assert took < 2.5, "stop() waited for what it may not signal: %.1f s" \
        % took
    finally:
      if copy is not None and lp._pid_exists(copy):
        os.kill(copy, 9)
      kill_and_unload(proc)
  saved = lp.launchctl
  lp.launchctl = lambda *args: (1, "", "")
  child = subprocess.Popen(["/bin/sleep", "60"])
  other = subprocess.Popen(["/bin/sleep", "60"])
  try:
    handle = lp.LaunchdProcess(TEST_PREFIX + ".fake", child.pid)
    listings = [[child.pid, other.pid]]
    handle.job_pids = lambda: listings.pop() if listings else []
    with _signals_only_to({child.pid}, foreign={other.pid}):
      handle.stop(timeout=5)
    assert child.wait(timeout=5) == -signal.SIGTERM, child.returncode
  finally:
    lp.launchctl = saved
    for p in (child, other):
      if p.poll() is None:
        p.kill()
        p.wait()


def exercise_stop_bounded():
  """stop() never hangs or reads a failure as success:
  - a coalition listing that fails (proc_listallpids returns -1) raises
    OSError, before anything is signalled, instead of reading as 'nothing
    left' and then waiting for a program that still runs;
  - a program that has not exited STOP_WAIT seconds after the signals
    raises subprocess.TimeoutExpired;
  - stop(timeout=None) waits for the SIGTERM without a limit, rather than
    raising TypeError once the SIGTERM has gone out.
  Handles for plain child processes. macOS only (kqueue)."""
  if sys.platform != "darwin":
    return
  import signal
  import subprocess
  import threading
  import types
  saved = lp.launchctl, lp._libproc, lp.STOP_WAIT
  lp.launchctl = lambda *args: (1, "", "")
  children = []
  def child():
    children.append(subprocess.Popen(["/bin/sleep", "60"]))
    return children[-1], lp.LaunchdProcess(TEST_PREFIX + ".fake",
                                           children[-1].pid)
  try:
    p, handle = child()
    handle.coalition = 999999999       # as for a launchd job
    lp._libproc = lambda: types.SimpleNamespace(
      proc_listallpids=lambda buf, size: -1)
    with _signals_only_to({p.pid}):
      for call in (handle.job_pids, handle.stop):
        try:
          call()
        except OSError:
          pass
        else:
          raise AssertionError("%s: a failed listing read as empty" % call)
    lp._libproc = saved[1]
    assert p.poll() is None, "signalled despite the failed listing"
    p, handle = child()
    handle.job_pids = lambda: []       # misses the program
    lp.STOP_WAIT = 0.3
    raised = []
    def stop():
      try:
        handle.stop(timeout=0.1)
      except subprocess.TimeoutExpired as exc:
        raised.append(exc)
    stopper = threading.Thread(target=stop, daemon=True)
    stopper.start()
    stopper.join(10)
    assert len(raised) == 1, "stop() hung or returned for a running program"
    p, handle = child()
    with _signals_only_to({p.pid}):
      handle.stop(timeout=None)
    assert p.wait(timeout=5) == -signal.SIGTERM, p.returncode
  finally:
    lp.launchctl, lp._libproc, lp.STOP_WAIT = saved
    for p in children:
      if p.poll() is None:
        p.kill()
        p.wait()


def exercise_undecodable_environment():
  """A launchd job's environment is a plist of Unicode strings, so a
  variable whose value is not valid UTF-8 (held by Python as surrogate
  escapes) cannot be passed. It is left out with a warning naming it; the
  program still starts and gets the other variables."""
  if not launchd_tests_can_run():
    return
  with tempfile.TemporaryDirectory() as tmp:
    log_path = os.path.join(tmp, "out.log")
    env = dict(os.environ, TST_BAD=_UNDECODABLE, TST_GOOD="kept")
    warnings, proc = [], None
    try:
      proc = lp.launch(_sh('echo "good=$TST_GOOD"\nexec /bin/sleep 60'),
                       env=env, log_path=log_path, warn=warnings.append,
                       label_prefix=TEST_PREFIX)
      assert len(warnings) == 1 and "TST_BAD" in warnings[0], warnings
      assert "good=kept\n" in _wait_for_text(log_path, "good=")
    finally:
      kill_and_unload(proc)


def exercise_children_outlive():
  """As with Popen, a process the program started keeps running after the
  program exits. launchd's default is to kill whatever is left in a job's
  process group when the job exits."""
  if not launchd_tests_can_run():
    return
  import signal
  with tempfile.TemporaryDirectory() as tmp:
    child_file = os.path.join(tmp, "child.pid")
    proc, child = None, None
    try:
      proc = lp.launch(
        _sh('/bin/sleep 60 &\necho $! > "%s"\n/bin/sleep 0.5\nexit 0'
            % child_file),
        label_prefix=TEST_PREFIX)
      assert _wait_exit(proc) is not None, "the program did not exit"
      child = int(open(child_file).read())
      time.sleep(0.5)                  # give launchd time to reap the group
      os.kill(child, 0)                # raises if launchd killed it
    finally:
      if child is not None:
        try:
          os.kill(child, signal.SIGKILL)
        except OSError:
          pass
      kill_and_unload(proc)


def _bootstrap_job(label, argv, run_at_load=True):
  """Test setup: load a launchd job ``label`` that runs ``argv`` once (or,
  without ``run_at_load``, is loaded and never run)."""
  import plistlib
  with tempfile.NamedTemporaryFile(suffix=".plist", delete=False) as f:
    plistlib.dump({"Label": label, "ProgramArguments": argv,
                   "RunAtLoad": run_at_load, "KeepAlive": False}, f)
  try:
    status, _, err = lp.launchctl("bootstrap", "gui/%d" % os.getuid(), f.name)
    assert status == 0, err
  finally:
    os.remove(f.name)


def exercise_cleans_up_stale_jobs():
  """A program is usually left running when its starter exits, so its job
  can outlive the starter that would have unloaded it. Launching unloads
  the finished jobs (same label prefix) of starters that have exited. It
  leaves alone a program still running after its starter exited, and the
  jobs of a live starter."""
  if not launchd_tests_can_run():
    return
  import subprocess
  gone = subprocess.Popen(["/usr/bin/true"])
  gone.wait()                                  # a starter pid that has exited
  stale = lp.job_label(TEST_PREFIX, gone.pid)
  orphan_running = lp.job_label(TEST_PREFIX, gone.pid)
  live_starter = lp.job_label(TEST_PREFIX, os.getpid())
  proc = None
  try:
    _bootstrap_job(stale, ["/usr/bin/true"])
    _bootstrap_job(live_starter, ["/usr/bin/true"])
    _bootstrap_job(orphan_running, ["/bin/sleep", "60"])
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and (
        lp.job_state(stale)[0] is not None or
        lp.job_state(live_starter)[0] is not None):
      time.sleep(0.1)
    proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
    assert not job_loaded(stale), "stale finished job kept"
    assert job_loaded(live_starter), "a live starter's job was unloaded"
    assert lp.job_state(orphan_running)[0] is not None, \
      "a running program whose starter exited was stopped"
  finally:
    kill_and_unload(proc)
    for label in (stale, live_starter, orphan_running):
      lp.launchctl("bootout", lp.job_target(label))


def exercise_umask():
  """The program gets the caller's umask, as a Popen child would; launchd's
  own is 022, which would give files the program saves the wrong
  permissions (readable by others under a private 077, not group-writable
  under a shared 002)."""
  if not launchd_tests_can_run():
    return
  with tempfile.TemporaryDirectory() as tmp:
    log_path = os.path.join(tmp, "out.log")
    saved, proc = os.umask(0o077), None
    try:
      proc = lp.launch(_sh("umask\nexec /bin/sleep 60"), log_path=log_path,
                       label_prefix=TEST_PREFIX)
    finally:
      os.umask(saved)
    try:
      text = _wait_for_text(log_path, "\n")
      assert text.strip() == "0077", text
    finally:
      kill_and_unload(proc)


def exercise_separate_stderr():
  """stderr can go to a file of its own (stderr_path), as a Popen caller
  can give stdout and stderr different destinations; without one it goes
  to log_path with stdout."""
  if not launchd_tests_can_run():
    return
  with tempfile.TemporaryDirectory() as tmp:
    out_path = os.path.join(tmp, "out.log")
    err_path = os.path.join(tmp, "err.log")
    proc = None
    try:
      proc = lp.launch(_sh('echo to-stdout\necho to-stderr >&2\n'
                           'exec /bin/sleep 60'),
                       log_path=out_path, stderr_path=err_path,
                       label_prefix=TEST_PREFIX)
      assert _wait_for_text(err_path, "to-stderr\n") == "to-stderr\n"
      assert _wait_for_text(out_path, "to-stdout\n") == "to-stdout\n"
    finally:
      kill_and_unload(proc)


def exercise_umask_read():
  """Reading the umask means setting it for a moment. The value set
  meanwhile is 077, so a file another thread creates in that moment is
  private, never readable by others. The umask is then restored. Not on
  Windows, which keeps only a umask's read/write bits (no launchd there)."""
  if sys.platform == "win32":
    return
  calls, real = [], os.umask
  def recording(mask):
    calls.append(mask)
    return real(mask)
  before = real(0o027)
  try:
    os.umask = recording
    assert lp._current_umask() == 0o027
  finally:
    os.umask = real
    after = real(before)
  assert after == 0o027 and calls == [0o077, 0o027], (oct(after), calls)


def exercise_utf8_environment_names():
  """A variable dropped for not being UTF-8 is named in the warning even
  when its name has non-ASCII characters (the warning itself must not
  fail). Runs on every platform."""
  warnings = []
  kept = lp.utf8_environment(
    {u"CAF\u00c9_DIR": _UNDECODABLE, "PLAIN": "x"},
    warn=warnings.append)
  assert kept == {"PLAIN": "x"}, kept
  assert len(warnings) == 1 and "CAF\\xc9_DIR" in warnings[0], warnings


def exercise_deleted_working_directory():
  """A caller whose working directory has been deleted can still launch, as
  with Popen; the program then starts in the root directory."""
  if not launchd_tests_can_run():
    return
  saved_cwd, proc = os.getcwd(), None
  tmp = tempfile.mkdtemp()
  try:
    os.chdir(tmp)
    os.rmdir(tmp)
    proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
    assert proc.poll() is None, "not running"
  finally:
    os.chdir(saved_cwd)
    kill_and_unload(proc)


def exercise_popen_timeouts():
  """wait(timeout) and communicate(timeout) behave as Popen's: they raise
  subprocess.TimeoutExpired while the program runs, and once it has exited
  return its returncode and (None, None) (there are no pipes)."""
  if not launchd_tests_can_run():
    return
  import signal
  import subprocess
  proc = lp.launch(["/bin/sleep", "60"], label_prefix=TEST_PREFIX)
  try:
    for call in (lambda: proc.wait(timeout=0.3),
                 lambda: proc.communicate(timeout=0.3)):
      began = time.monotonic()
      try:
        call()
      except subprocess.TimeoutExpired:
        pass
      else:
        raise AssertionError("no TimeoutExpired while running")
      assert time.monotonic() - began < 2, "timeout not honored"
    os.killpg(proc.pid, signal.SIGTERM)
    assert proc.communicate(timeout=10) == (None, None)
    assert proc.wait(timeout=10) == -15, proc.returncode
  finally:
    kill_and_unload(proc)


def exercise():
  exercise_use_launchd()
  exercise_own_coalition()
  exercise_lifecycle()
  exercise_job_environment()
  exercise_immediate_exit_leaves_no_job()
  exercise_failures()
  exercise_program_lookup()
  exercise_job_state_parsing()
  exercise_job_runs()
  exercise_exit_recording()
  exercise_concurrent_wait_timeouts()
  exercise_live_test_gate()
  exercise_describe_returncode()
  exercise_exited_at_once_pickles()
  exercise_spawn()
  exercise_separate_stderr()
  exercise_umask_read()
  exercise_signals()
  exercise_stop_reaches_relaunched_copy()
  exercise_stop_own_coalition_guard()
  exercise_stop_needs_the_jobs_own_program()
  exercise_stop_skips_foreign_processes()
  exercise_stop_bounded()
  exercise_undecodable_environment()
  exercise_children_outlive()
  exercise_cleans_up_stale_jobs()
  exercise_umask()
  exercise_utf8_environment_names()
  exercise_deleted_working_directory()
  exercise_popen_timeouts()


if __name__ == "__main__":
  exercise()
  print("OK")
