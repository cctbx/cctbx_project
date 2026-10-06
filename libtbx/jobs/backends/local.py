"""Run jobs as local subprocesses (spec 7)."""

import os
import subprocess
import sys
import threading

from libtbx.utils import Sorry

from libtbx.jobs.backends.base import Backend, PollResult, SubmitError


class _LocalHandle(object):
  """The subprocess behind a job.

  Attributes
  ----------
  proc : subprocess.Popen
  pid : int
  signalled : bool
      True once this backend has sent the process group a signal.
  """

  __slots__ = ("proc", "pid", "signalled")

  def __init__(self, proc):
    self.proc = proc
    self.pid = proc.pid
    self.signalled = False


class LocalBackend(Backend):
  """Spawn each job in its own session on this machine.

  Parameters
  ----------
  max_concurrent : int, optional
      Slots the manager may fill (default 8, minimum 1).
  """

  name = "local"
  handle_keys = ("pid",)

  def __init__(self, max_concurrent=8):
    if sys.platform == "win32":
      raise Sorry("LocalBackend requires a POSIX platform")
    self.max_concurrent = max(1, int(max_concurrent))
    self._lock = threading.Lock()

  def capacity(self):
    return self.max_concurrent

  def submit(self, job):
    spec = job.spec
    env = dict(os.environ)
    # Do not leak a surrounding batch job's attempt count to children.
    env.pop("LIBTBX_JOBS_ATTEMPT", None)
    env.update(spec.env)
    try:
      out = open(spec.log_path, "a", encoding="utf-8")
    except OSError:
      out = None
    try:
      proc = subprocess.Popen(
        list(spec.argv),
        cwd=spec.cwd,
        env=env,
        stdout=out if out is not None else subprocess.DEVNULL,
        stderr=subprocess.STDOUT if out is not None else subprocess.DEVNULL,
        start_new_session=True,
      )
    except OSError as exc:
      raise SubmitError(str(exc))
    finally:
      if out is not None:
        out.close()
    return "running", _LocalHandle(proc)

  def poll(self, job):
    handle = job.handle
    with self._lock:
      rc = handle.proc.poll()
    if rc is None:
      return PollResult("running")
    reason = "killed by signal %d" % -rc if rc < 0 else None
    if not handle.signalled:
      cancelled = False
    elif rc < 0:
      cancelled = True
    else:
      cancelled = None
    return PollResult("exited", exit_code=rc, reason=reason,
                      cancelled=cancelled)

  def _signal(self, job, signum):
    handle = job.handle
    with self._lock:
      if handle.proc.poll() is not None:
        return
      try:
        os.killpg(os.getpgid(handle.pid), signum)
        handle.signalled = True
      except ProcessLookupError:
        pass

  def terminate(self, job):
    import signal
    self._signal(job, signal.SIGTERM)

  def kill(self, job):
    import signal
    self._signal(job, signal.SIGKILL)

  def on_shutdown(self, jobs):
    import signal
    for job in jobs:
      handle = job.handle
      if handle is None:
        continue
      try:
        os.killpg(os.getpgid(handle.pid), signal.SIGTERM)
      except OSError:
        pass

  def handle_info(self, job):
    return {"pid": job.handle.pid if job.handle is not None else None}
