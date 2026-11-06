"""LocalBackend tests with real subprocesses.

The ``exercise_windows_*`` tests run the Windows code path on a POSIX
host, with ``_FakeJobObjects`` in place of the Win32 job-object calls.
On Windows the other tests exercise the real job objects.
"""

import os
import signal
import sys
import tempfile
import time

from libtbx.utils import Sorry

WINDOWS = sys.platform == "win32"


def _manager(**kw):
  from libtbx.jobs.backends.local import LocalBackend
  from libtbx.jobs.manager import JobManager
  kw.setdefault("reaper_interval", 3600)
  backend = LocalBackend(max_concurrent=kw.pop("max_concurrent", 4))
  return JobManager(backend, **kw), backend


def _py_spec(tmp, code, **kw):
  from libtbx.jobs.manager import JobSpec
  kw.setdefault("cwd", tmp)
  kw.setdefault("log_path", os.path.join(tmp, "job.log"))
  return JobSpec(argv=[sys.executable, "-c", code], **kw)


def _read(path):
  with open(path) as fh:
    return fh.read()


def _await(manager, job, states, timeout=30):
  deadline = time.time() + timeout
  while time.time() < deadline:
    manager.reaper_iter()
    if job.state in states:
      return
    time.sleep(0.05)
  raise AssertionError("job stuck in %s" % job.state)


def _await_text(path, text, timeout=10):
  deadline = time.time() + timeout
  while not (os.path.exists(path) and text in _read(path)):
    assert time.time() < deadline, "%r never appeared in %s" % (text, path)
    time.sleep(0.05)


def _process_gone(pid, timeout=10):
  """True once ``pid`` no longer exists (POSIX only)."""
  deadline = time.time() + timeout
  while time.time() < deadline:
    try:
      os.kill(pid, 0)
    except ProcessLookupError:
      return True
    time.sleep(0.05)
  return False


class _FakeJobObjects(object):
  """The Win32 job-object layer, faked on a POSIX host.

  ``terminate`` sends SIGTERM, which ``_FAKE_WINDOWS_CHILD`` turns into
  exit code 1, the code a real ``TerminateJobObject`` call leaves.
  """

  creation_flags = 0             # POSIX Popen refuses anything else

  def __init__(self, fail=None):
    self.fail = fail
    self.calls = []
    self._pids = {}

  def adopt(self, pid):
    self.calls.append(("adopt", pid))
    if self.fail is not None:
      raise self.fail
    job_object = "job-object-%d" % len(self._pids)
    self._pids[job_object] = pid
    return job_object

  def terminate(self, job_object):
    self.calls.append(("terminate", job_object))
    os.kill(self._pids[job_object], signal.SIGTERM)

  def close(self, job_object):
    self.calls.append(("close", job_object))


# setsid: without a session of its own, a stray os.killpg in the code
# under test would signal this test's process group
_FAKE_WINDOWS_CHILD = """
import os, signal, time
try:
  os.setsid()
except OSError:
  pass
signal.signal(signal.SIGTERM, lambda *a: os._exit(1))
print('ready', flush=True)
time.sleep(60)
"""


def _fake_windows_manager(fail=None, **kw):
  manager, backend = _manager(**kw)
  fake = _FakeJobObjects(fail=fail)
  backend._job_objects = fake
  return manager, backend, fake


def exercise_exit_codes_and_log_capture():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    ok = manager.submit(_py_spec(
      tmp, "import sys; print('out'); print('err', file=sys.stderr)"))
    assert ok.state == "running" and ok.handle.pid > 0
    assert manager.get_status(ok.job_id)["pid"] == ok.handle.pid
    _await(manager, ok, ("finished",))
    assert ok.exit_code == 0 and ok.reason is None
    log = _read(os.path.join(tmp, "job.log"))
    assert "out\n" in log and "err\n" in log, log
    bad = manager.submit(_py_spec(tmp, "import sys; sys.exit(3)",
                                  log_path=os.path.join(tmp, "bad.log")))
    _await(manager, bad, ("failed",))
    assert bad.exit_code == 3, bad.exit_code
  finally:
    manager.stop()


def exercise_env_override_reaches_the_child():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  os.environ["LIBTBX_JOBS_ATTEMPT"] = "3"     # as if inside a requeued batch job
  try:
    job = manager.submit(_py_spec(
      tmp, "import os; print(os.environ['LIBTBX_JOBS_TEST']); "
           "print('attempt', os.environ.get('LIBTBX_JOBS_ATTEMPT', 'unset'))",
      env={"LIBTBX_JOBS_TEST": "value-1"}))
    _await(manager, job, ("finished",))
    log = _read(os.path.join(tmp, "job.log"))
    assert "value-1" in log
    assert "attempt unset" in log, "the batch attempt count leaked to a local child"
  finally:
    os.environ.pop("LIBTBX_JOBS_ATTEMPT", None)
    manager.stop()


def exercise_children_do_not_inherit_the_managers_stdin():
  """The manager may be an MCP server whose stdin is the JSON-RPC pipe. A
  child that reads stdin would swallow protocol messages, so children get
  the null device."""
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  r, w = os.pipe()
  os.write(w, b"a JSON-RPC message\n")
  os.close(w)
  saved = os.dup(0)
  os.dup2(r, 0)
  os.close(r)
  try:
    job = manager.submit(_py_spec(
      tmp, "import sys; print('STDIN=%r' % sys.stdin.read())"))
    _await(manager, job, ("finished",))
    log = _read(os.path.join(tmp, "job.log"))
    assert "STDIN=''" in log, log
  finally:
    os.dup2(saved, 0)
    os.close(saved)
    manager.stop()


def exercise_sigterm_cancel_is_reported_cancelled():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(tmp, "import time; time.sleep(30)"))
    d = manager.cancel(job.job_id)
    assert d["state"] == "cancelling", d
    assert job.handle.signalled is True
    _await(manager, job, ("cancelled",))
    if WINDOWS:
      assert job.exit_code == 1, job.exit_code
      assert job.reason == "terminated", job.reason
    else:
      assert job.exit_code == -signal.SIGTERM, job.exit_code
      assert job.reason == "killed by signal %d" % signal.SIGTERM
  finally:
    manager.stop()


_GRANDCHILD = """
import sys, time
while True:
  with open(sys.argv[1], "a") as fh:
    fh.write("x")
  time.sleep(0.05)
"""


def exercise_cancel_ends_the_whole_process_tree():
  """A Phenix job is a runner that starts the program as its own child;
  cancelling must stop both, or the program keeps running unseen."""
  tmp = tempfile.mkdtemp()
  beat = os.path.join(tmp, "beat")
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(
      tmp, "import subprocess, sys, time; "
           "subprocess.Popen([sys.executable, '-c', %r, %r]); "
           "time.sleep(60)" % (_GRANDCHILD, beat)))
    _await_text(beat, "x")
    manager.cancel(job.job_id)
    _await(manager, job, ("cancelled",))
    time.sleep(0.5)
    size = os.path.getsize(beat)
    time.sleep(0.5)
    assert os.path.getsize(beat) == size, "the grandchild outlived the cancel"
  finally:
    manager.stop()


def exercise_sigkill_escalation_for_a_child_ignoring_sigterm():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager(sigkill_after=0.3)
  try:
    job = manager.submit(_py_spec(
      tmp, "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
           "print('ready', flush=True); time.sleep(60)"))
    deadline = time.time() + 10
    while "ready" not in _read(os.path.join(tmp, "job.log")):
      assert time.time() < deadline, "child never started"
      time.sleep(0.05)
    manager.cancel(job.job_id)
    _await(manager, job, ("cancelled",))
    assert job.exit_code == -signal.SIGKILL, job.exit_code
    assert job._kill_sent
  finally:
    manager.stop()


def exercise_natural_end_before_cancel_is_finished():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(tmp, "pass"))
    job.handle.proc.wait(timeout=30)
    r = backend.poll(job)
    assert r == PollResult("exited", exit_code=0, cancelled=False), r
    d = manager.cancel(job.job_id)
    assert d["state"] == "finished", d
  finally:
    manager.stop()


def exercise_missing_executable_is_a_failed_job():
  from libtbx.jobs.manager import JobSpec
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(JobSpec(argv=["/nonexistent/binary"], cwd=tmp,
                                 log_path=os.path.join(tmp, "job.log")))
    assert job.state == "failed", job.state
    assert "/nonexistent/binary" in job.reason, job.reason
    assert "SORRY: failed to submit job" in _read(os.path.join(tmp, "job.log"))
  finally:
    manager.stop()


def exercise_spawn_error_names_the_program():
  """Windows spawn errors name no file ("[WinError 267] The directory
  name is invalid"), and a POSIX one for a bad cwd names only the cwd."""
  from libtbx.jobs.manager import JobSpec
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(JobSpec(
      argv=[sys.executable, "-c", "pass"], cwd=os.path.join(tmp, "gone"),
      log_path=os.path.join(tmp, "job.log")))
    assert job.state == "failed", job.state
    assert sys.executable in job.reason, job.reason
  finally:
    manager.stop()


def exercise_unwritable_log_still_runs_the_job():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(
      tmp, "print('hello')", log_path=os.path.join(tmp, "missing", "job.log")))
    assert job.state == "running", job.state
    _await(manager, job, ("finished",))
    assert not os.path.exists(os.path.join(tmp, "missing"))
  finally:
    manager.stop()


def exercise_on_shutdown_terminates_live_groups():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(tmp, "import time; time.sleep(30)"))
    manager.shutdown()
    job.handle.proc.wait(timeout=10)
    expected = 1 if WINDOWS else -signal.SIGTERM
    assert job.handle.proc.returncode == expected, job.handle.proc.returncode
  finally:
    manager.stop()


def exercise_reattach_raises_sorry():
  from libtbx.jobs.backends.local import LocalBackend
  try:
    LocalBackend().reattach({"pid": 1})
  except Sorry as exc:
    assert "local" in str(exc), str(exc)
  else:
    raise AssertionError("reattach did not raise")


def exercise_windows_cancel_terminates_the_job_object():
  tmp = tempfile.mkdtemp()
  manager, backend, fake = _fake_windows_manager()
  try:
    job = manager.submit(_py_spec(tmp, _FAKE_WINDOWS_CHILD))
    pid = job.handle.pid
    assert fake.calls == [("adopt", pid)], fake.calls
    _await_text(os.path.join(tmp, "job.log"), "ready")
    manager.cancel(job.job_id)
    _await(manager, job, ("cancelled",))
    assert job.exit_code == 1, job.exit_code
    assert job.reason == "terminated", job.reason
    assert fake.calls == [("adopt", pid), ("terminate", "job-object-0"),
                          ("close", "job-object-0")], fake.calls
  finally:
    manager.stop()


def exercise_windows_exit_code_1_without_a_cancel_is_a_failure():
  tmp = tempfile.mkdtemp()
  manager, backend, fake = _fake_windows_manager()
  try:
    job = manager.submit(_py_spec(tmp, "import sys; sys.exit(1)"))
    _await(manager, job, ("finished", "failed", "cancelled"))
    assert job.state == "failed" and job.exit_code == 1, (job.state, job.exit_code)
    assert job.reason is None, job.reason
    assert fake.calls == [("adopt", job.handle.pid),
                          ("close", "job-object-0")], fake.calls
  finally:
    manager.stop()


def exercise_windows_kill_terminates_the_job_object():
  """Windows has no SIGKILL; kill ends the job object like terminate."""
  tmp = tempfile.mkdtemp()
  manager, backend, fake = _fake_windows_manager()
  try:
    job = manager.submit(_py_spec(tmp, _FAKE_WINDOWS_CHILD))
    _await_text(os.path.join(tmp, "job.log"), "ready")
    backend.kill(job)
    assert fake.calls[1:] == [("terminate", "job-object-0")], fake.calls
    assert job.handle.signalled is True
    assert job.handle.proc.wait(timeout=10) == 1
  finally:
    manager.stop()


def exercise_windows_on_shutdown_terminates_the_job_object():
  tmp = tempfile.mkdtemp()
  manager, backend, fake = _fake_windows_manager()
  try:
    job = manager.submit(_py_spec(tmp, _FAKE_WINDOWS_CHILD))
    _await_text(os.path.join(tmp, "job.log"), "ready")
    manager.shutdown()
    assert fake.calls[1:] == [("terminate", "job-object-0")], fake.calls
    assert job.handle.proc.wait(timeout=10) == 1
  finally:
    manager.stop()


def exercise_windows_adopt_failure_fails_the_job_and_ends_the_child():
  """On Windows the child is still suspended when adopt fails; it must
  not be left behind."""
  tmp = tempfile.mkdtemp()
  import ctypes
  for error in (OSError("Access is denied"),
                ctypes.ArgumentError("argument 1: wrong type")):
    manager, backend, fake = _fake_windows_manager(fail=error)
    try:
      job = manager.submit(_py_spec(tmp, _FAKE_WINDOWS_CHILD))
      assert job.state == "failed", job.state
      assert str(error) in job.reason, job.reason
      assert job.handle is None
      pid = fake.calls[0][1]
      assert _process_gone(pid, timeout=0.5), \
        "the child %d was left running after %r" % (pid, error)
    finally:
      manager.stop()


def exercise_job_object_limits_have_the_win64_layout():
  """SetInformationJobObject takes the struct by size; one wrong field
  and every Windows submit fails. Sizes and offsets are the Win64 SDK's."""
  import ctypes
  from libtbx.jobs.backends import local
  if ctypes.sizeof(ctypes.c_void_p) != 8:
    return
  assert ctypes.sizeof(local._BasicLimitInformation) == 64
  assert local._BasicLimitInformation.LimitFlags.offset == 16
  assert local._BasicLimitInformation.ActiveProcessLimit.offset == 40
  assert local._BasicLimitInformation.SchedulingClass.offset == 60
  assert local._ExtendedLimitInformation.IoInfo.offset == 64
  assert local._ExtendedLimitInformation.ProcessMemoryLimit.offset == 112
  assert ctypes.sizeof(local._ExtendedLimitInformation) == 144


def run():
  exercise_job_object_limits_have_the_win64_layout()
  exercise_exit_codes_and_log_capture()
  exercise_env_override_reaches_the_child()
  exercise_children_do_not_inherit_the_managers_stdin()
  exercise_sigterm_cancel_is_reported_cancelled()
  exercise_cancel_ends_the_whole_process_tree()
  if not WINDOWS:
    # a child that ignores SIGTERM; Windows has neither signal
    exercise_sigkill_escalation_for_a_child_ignoring_sigterm()
  exercise_natural_end_before_cancel_is_finished()
  exercise_missing_executable_is_a_failed_job()
  exercise_spawn_error_names_the_program()
  exercise_unwritable_log_still_runs_the_job()
  exercise_on_shutdown_terminates_live_groups()
  exercise_reattach_raises_sorry()
  if not WINDOWS:
    exercise_windows_cancel_terminates_the_job_object()
    exercise_windows_exit_code_1_without_a_cancel_is_a_failure()
    exercise_windows_kill_terminates_the_job_object()
    exercise_windows_on_shutdown_terminates_the_job_object()
    exercise_windows_adopt_failure_fails_the_job_and_ends_the_child()
  print("OK")


if __name__ == "__main__":
  run()
