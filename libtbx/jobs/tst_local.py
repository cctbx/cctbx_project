"""LocalBackend tests with real subprocesses."""

import os
import signal
import sys
import tempfile
import time

from libtbx.utils import Sorry


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


def exercise_sigterm_cancel_is_reported_cancelled():
  tmp = tempfile.mkdtemp()
  manager, backend = _manager()
  try:
    job = manager.submit(_py_spec(tmp, "import time; time.sleep(30)"))
    d = manager.cancel(job.job_id)
    assert d["state"] == "cancelling", d
    assert job.handle.signalled is True
    _await(manager, job, ("cancelled",))
    assert job.exit_code == -signal.SIGTERM, job.exit_code
    assert job.reason == "killed by signal %d" % signal.SIGTERM
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
    assert "nonexistent" in job.reason, job.reason
    assert "SORRY: failed to submit job" in _read(os.path.join(tmp, "job.log"))
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
    assert job.handle.proc.returncode == -signal.SIGTERM
  finally:
    manager.stop()


def exercise_reattach_raises_sorry_and_windows_refuses():
  from libtbx.jobs.backends.local import LocalBackend
  try:
    LocalBackend().reattach({"pid": 1})
  except Sorry as exc:
    assert "local" in str(exc), str(exc)
  else:
    raise AssertionError("reattach did not raise")
  real = sys.platform
  sys.platform = "win32"
  try:
    try:
      LocalBackend()
    except Sorry as exc:
      assert "POSIX" in str(exc), str(exc)
    else:
      raise AssertionError("LocalBackend() did not raise on win32")
  finally:
    sys.platform = real


def run():
  if sys.platform == "win32":
    print("Skipping tst_local on Windows")
    print("OK")
    return
  exercise_exit_codes_and_log_capture()
  exercise_env_override_reaches_the_child()
  exercise_sigterm_cancel_is_reported_cancelled()
  exercise_sigkill_escalation_for_a_child_ignoring_sigterm()
  exercise_natural_end_before_cancel_is_finished()
  exercise_missing_executable_is_a_failed_job()
  exercise_unwritable_log_still_runs_the_job()
  exercise_on_shutdown_terminates_live_groups()
  exercise_reattach_raises_sorry_and_windows_refuses()
  print("OK")


if __name__ == "__main__":
  run()
