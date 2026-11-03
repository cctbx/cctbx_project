"""Tests for libtbx.jobs.manager driven by a scripted fake backend.

Each exercise drives ``reaper_iter`` by hand; the reaper interval is long
enough that the thread never ticks on its own.
"""

import sys
import tempfile
import threading
import time
import os

from libtbx.utils import Sorry


def exercise_poll_result_validates_phase():
  from libtbx.jobs.backends.base import PollResult
  for phase in ("pending", "running", "exited"):
    assert PollResult(phase).phase == phase
  try:
    PollResult("done")
  except ValueError as exc:
    assert "done" in str(exc), str(exc)
  else:
    raise AssertionError("PollResult accepted an unknown phase")
  r = PollResult("exited", exit_code=3, reason="FAILED", cancelled=False)
  assert (r.exit_code, r.reason, r.cancelled) == (3, "FAILED", False)


def exercise_backend_defaults():
  from libtbx.jobs.backends.base import Backend
  b = Backend()
  assert b.capacity() is None
  assert b.handle_info(None) == {}
  assert b.kill(None) is None
  assert b.on_shutdown([]) is None
  assert b.refresh([]) is None
  assert b.validate(None) is None
  try:
    b.reattach("x")
  except Sorry as exc:
    assert "abstract" in str(exc), str(exc)
  else:
    raise AssertionError("base reattach did not raise Sorry")
  for name in ("submit", "poll", "terminate"):
    try:
      getattr(b, name)(None)
    except NotImplementedError:
      pass
    else:
      raise AssertionError("%s did not raise NotImplementedError" % name)


def _spec(**kw):
  from libtbx.jobs.manager import JobSpec
  base = dict(argv=["true"], cwd="/tmp", log_path="/tmp/job.log")
  base.update(kw)
  return JobSpec(**base)


def exercise_jobspec_validation():
  from libtbx.jobs.manager import JobSpec
  s = _spec(name="", env={"A": "1"}, resources={"cpus": 2, "partition": "gpu"},
            metadata={"program": "x"})
  assert s.name == "job", s.name
  assert s.argv == ["true"]
  bad = [
    dict(argv=[]),
    dict(argv="true"),
    dict(cwd="relative/dir"),
    dict(log_path="job.log"),
    dict(env={"A": 1}),
    dict(env={"A;x": "1"}),                     # shell metacharacter
    dict(env={"X=$(id)": "2"}),
    dict(env={"1ABC": "x"}),                    # leading digit
    dict(resources={"cpus": 0}),
    dict(resources={"gpus": -1}),
    dict(resources={"time_minutes": "10"}),
    dict(resources={"partition": True}),           # bool is an int subclass
    dict(resources={"partition": "a b"}),          # whitespace
    dict(resources={"partition": "x" * 257}),      # too long
    dict(resources={"partition": "a\nb"}),         # newline
    dict(resources={"partition": "gpu\n"}),        # trailing newline
    dict(resources={"partition": None}),
    dict(log_path="/tmp/job\0.log"),                # NUL: open() would raise
    dict(cwd="/tmp/a\0b"),
    dict(argv=["true", "a\0b"]),
    dict(env={"A": "x\0y"}),
    dict(metadata={1: "x"}),
    dict(metadata={"obj": object()}),
  ]
  for kw in bad:
    try:
      _spec(**kw)
    except ValueError:
      pass
    else:
      raise AssertionError("JobSpec accepted %r" % (kw,))
  try:
    s.name = "other"
  except AttributeError:
    pass
  else:
    raise AssertionError("JobSpec is not frozen")
  assert JobSpec(argv=("a", "b"), cwd="/x", log_path="/x/l").argv == ["a", "b"]


def exercise_job_to_dict_layers_core_handle_and_metadata():
  from libtbx.jobs.manager import Job, CORE_STATUS_KEYS, _new_job_id
  from libtbx.jobs.backends.base import Backend

  class B(Backend):
    handle_keys = ("pid",)
    def handle_info(self, job):
      return {"pid": 42}

  spec = _spec(metadata={"program": "refine", "args": ["a"]})
  job = Job(job_id=_new_job_id(), spec=spec)
  assert job.job_id.startswith("j_") and len(job.job_id) == 14, job.job_id
  d = job.to_dict(B())
  for key in CORE_STATUS_KEYS:
    assert key in d, key
  assert d["state"] == "queued" and d["queue_position"] == 0
  assert d["pid"] == 42 and d["program"] == "refine" and d["args"] == ["a"]
  assert d["started_at"] is None and d["exit_code"] is None
  assert d["reason"] is None and d["log_path"] == "/tmp/job.log"
  assert d["args"] is not spec.metadata["args"], "metadata lists must be copied"


class FakeBackend(object):
  """Scripted backend: records calls and asserts the lock discipline.

  ``results[handle]`` is a PollResult or a list of PollResults consumed one
  per poll (the last one repeats). ``manager`` is set by ``_make``.
  """

  name = "fake"
  handle_keys = ("fake_handle",)

  def __init__(self, capacity=None, submit_phase="running",
               submit_delay=0.0):
    self.cap = capacity
    self.submit_phase = submit_phase
    self.submit_delay = submit_delay
    self.results = {}
    self.calls = []
    self.fail_submit_with = None
    self.counter = 0
    self.manager = None
    self.shutdown_seen = None
    self.reattach_raises = None
    self.reattach_bad = {}
    self.crash_submit = False
    self.crash_submit_with = None
    self.violations = []
    self.submit_gate = None
    self.submit_entered = threading.Event()
    self.capacity_gate = None
    self.capacity_entered = threading.Event()
    self.refresh_hook = None
    self.polled = []

  def _lock(self, name, owned):
    # Record as well as assert: the manager guards several backend calls
    # with ``except Exception``, which would otherwise swallow the assert.
    if self.manager is not None:
      held = self.manager.holds_lock()
      if held != owned:
        self.violations.append("%s called with lock owned=%s" % (name, held))
      assert held == owned, (
        "%s called with lock owned=%s" % (name, held))

  def validate(self, spec):
    self._lock("validate", False)

  def check(self, spec):
    self.validate(spec)

  def submit(self, job):
    from libtbx.jobs.backends.base import SubmitError, PollResult
    self._lock("submit", False)
    self.calls.append(("submit", job.job_id))
    if self.submit_delay:
      time.sleep(self.submit_delay)
    if self.submit_gate is not None:
      self.submit_entered.set()
      self.submit_gate.wait(10)
    if self.fail_submit_with:
      raise SubmitError(self.fail_submit_with)
    if self.crash_submit:
      raise RuntimeError("boom")
    if self.crash_submit_with is not None:
      raise self.crash_submit_with
    self.counter += 1
    handle = "h%d" % self.counter
    self.results.setdefault(handle, PollResult(
      "running" if self.submit_phase == "running" else "pending"))
    return self.submit_phase, handle

  def refresh(self, jobs, force=False):
    self._lock("refresh", False)
    self.calls.append(("refresh", tuple(j.job_id for j in jobs), force))
    if self.refresh_hook is not None:
      self.refresh_hook(jobs)

  def poll(self, job):
    from libtbx.jobs.backends.base import PollResult
    self._lock("poll", True)
    self.polled.append(job.job_id)
    r = self.results.get(job.handle, PollResult("running"))
    if isinstance(r, list):
      r = r.pop(0) if len(r) > 1 else r[0]
    if isinstance(r, Exception):
      raise r
    return r

  def terminate(self, job):
    self._lock("terminate", False)
    self.calls.append(("terminate", job.job_id))

  def kill(self, job):
    self._lock("kill", False)
    self.calls.append(("kill", job.job_id))

  def capacity(self):
    self._lock("capacity", False)
    if self.capacity_gate is not None:
      self.capacity_entered.set()
      self.capacity_gate.wait(10)
    return self.cap

  def on_shutdown(self, jobs):
    self.shutdown_seen = [j.job_id for j in jobs]

  def handle_info(self, job):
    self._lock("handle_info", True)
    return {"fake_handle": job.handle}

  def reattach(self, handle):
    self._lock("reattach", False)
    if self.reattach_raises:
      raise self.reattach_raises
    if handle in self.reattach_bad:
      raise self.reattach_bad[handle]
    return "re-" + str(handle)


def _make(capacity=None, **kw):
  from libtbx.jobs.manager import JobManager
  backend = FakeBackend(capacity=capacity)
  kw.setdefault("reaper_interval", 3600)
  manager = JobManager(backend, **kw)
  backend.manager = manager
  return manager, backend


def _finish(manager, backend):
  """Stop the manager and fail on any recorded lock-discipline violation."""
  manager.stop()
  assert not backend.violations, backend.violations


def _tmp_spec(tmp, **kw):
  kw.setdefault("argv", ["true"])
  kw.setdefault("cwd", tmp)
  kw.setdefault("log_path", os.path.join(tmp, "job.log"))
  return _spec(**kw)


def exercise_submit_runs_immediately_and_reports_running():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=2)
  try:
    job = manager.submit(_tmp_spec(tmp, metadata={"program": "p"}))
    assert job.state == "running", job.state
    assert job.handle == "h1" and not job._submitting
    assert job.started_at is not None and job.queue_position == 0
    d = manager.get_status(job.job_id)
    assert d["state"] == "running" and d["fake_handle"] == "h1"
    assert d["program"] == "p" and d["reason"] is None
    assert manager.list_jobs() == [d]
    backend.results["h1"] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    d = manager.get_status(job.job_id)
    assert d["state"] == "finished" and d["exit_code"] == 0, d
    assert d["finished_at"] is not None
    assert job._terminal_event.is_set()
    backend.results["h1"] = PollResult("running")
    manager.reaper_iter()
    assert manager.get_status(job.job_id)["state"] == "finished", \
      "a poll result must not move a terminal job"
  finally:
    _finish(manager, backend)


def exercise_capacity_queues_and_drains_in_order():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    c = manager.submit(_tmp_spec(tmp))
    assert (a.state, b.state, c.state) == ("running", "queued", "queued")
    assert (b.queue_position, c.queue_position) == (1, 2)
    assert b.handle is None
    backend.results["h1"] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    assert b.state == "running" and b.handle == "h2", (b.state, b.handle)
    assert c.queue_position == 1
    # A slot free but the FIFO non-empty: a new job queues behind.
    backend.cap = 2
    d = manager.submit(_tmp_spec(tmp))
    assert d.state == "queued" and d.queue_position == 2, (d.state, d.queue_position)
    manager.reaper_iter()
    assert c.state == "running" and d.state == "queued" and d.queue_position == 1
  finally:
    _finish(manager, backend)


def exercise_capacity_none_mirrors_scheduler_pending():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=None)
  backend.submit_phase = "pending"
  try:
    job = manager.submit(_tmp_spec(tmp))
    assert job.state == "queued" and job.queue_position == 0 and job.handle == "h1"
    assert job.started_at is None
    backend.results["h1"] = PollResult("pending", reason="PENDING: Priority")
    manager.reaper_iter()
    d = manager.get_status(job.job_id)
    assert d["state"] == "queued" and d["reason"] == "PENDING: Priority", d
    backend.results["h1"] = PollResult("running", reason=None)
    manager.reaper_iter()
    assert job.state == "running" and job.started_at is not None
    started = job.started_at
    backend.results["h1"] = PollResult("pending", reason="REQUEUED")
    manager.reaper_iter()
    assert job.state == "queued" and job.started_at is None, "requeue clears started_at"
    backend.results["h1"] = PollResult("running")
    manager.reaper_iter()
    assert job.started_at is not None and job.started_at >= started
    backend.results["h1"] = PollResult("exited", exit_code=None, reason="NODE_FAIL")
    manager.reaper_iter()
    d = manager.get_status(job.job_id)
    assert d["state"] == "failed" and d["reason"] == "NODE_FAIL", d
    with open(os.path.join(tmp, "job.log")) as fh:
      log = fh.read()
    assert "SORRY: job ended without an exit code: NODE_FAIL" in log, log
  finally:
    _finish(manager, backend)


def exercise_external_cancellation_and_signal_exit():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    backend.results[a.handle] = PollResult("exited", exit_code=None,
                                           reason="CANCELLED", cancelled=True)
    backend.results[b.handle] = PollResult("exited", exit_code=-15,
                                           reason="killed by signal 15")
    manager.reaper_iter()
    assert a.state == "cancelled" and a.reason == "CANCELLED"
    assert b.state == "failed" and b.exit_code == -15
  finally:
    _finish(manager, backend)


def exercise_submit_failure_is_a_failed_job_not_an_exception():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  backend.fail_submit_with = "invalid partition"
  try:
    job = manager.submit(_tmp_spec(tmp))
    d = manager.get_status(job.job_id)
    assert d["state"] == "failed" and d["exit_code"] is None, d
    assert "invalid partition" in d["reason"], d
    with open(os.path.join(tmp, "job.log")) as fh:
      log = fh.read()
    assert "SORRY: failed to submit job: invalid partition" in log, log
    assert not job._submitting and job._terminal_event.is_set()
  finally:
    _finish(manager, backend)


def exercise_backend_submit_crash_is_a_failed_job():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  backend.crash_submit = True
  try:
    job = manager.submit(_tmp_spec(tmp))
    d = manager.get_status(job.job_id)
    assert d["state"] == "failed" and "RuntimeError" in d["reason"], d
    assert not job._submitting and job._terminal_event.is_set()
    backend.crash_submit = False
    second = manager.submit(_tmp_spec(tmp))
    assert second.state == "running" and second.handle == "h1", \
      "a crashed submit leaked its capacity slot"
    backend.results["h1"] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    assert second.state == "finished"
  finally:
    _finish(manager, backend)


def exercise_an_interrupt_inside_submit_frees_the_slot():
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  backend.crash_submit_with = KeyboardInterrupt()
  try:
    try:
      manager.submit(_tmp_spec(tmp))
    except KeyboardInterrupt:
      pass
    else:
      raise AssertionError("the interrupt was swallowed")
    job = list(manager._jobs.values())[0]
    assert job.state == "failed" and not job._submitting, (job.state, job._submitting)
    assert "interrupted by KeyboardInterrupt" in job.reason, job.reason
    assert job._terminal_event.is_set()
    with open(os.path.join(tmp, "job.log")) as fh:
      assert "SORRY: submission interrupted" in fh.read()
    backend.crash_submit_with = None
    later = manager.submit(_tmp_spec(tmp))
    assert later.state == "running", "the interrupted submission kept its slot"
  finally:
    _finish(manager, backend)


def exercise_unknown_submit_phase_fails_and_terminates_the_job():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  backend.submit_phase = "started"
  try:
    job = manager.submit(_tmp_spec(tmp))
    assert job.state == "failed" and "started" in job.reason, (job.state, job.reason)
    assert job.handle == "h1", "the handle must be kept so terminate can use it"
    assert backend.calls[-1] == ("terminate", job.job_id), backend.calls
    with open(os.path.join(tmp, "job.log")) as fh:
      assert "SORRY: backend reported unknown phase" in fh.read()
  finally:
    _finish(manager, backend)


def _returns_within(fn, seconds):
  """Run ``fn`` on a thread; True when it returns within ``seconds``."""
  done = threading.Event()
  threading.Thread(target=lambda: (fn(), done.set()), daemon=True).start()
  return done.wait(seconds)


def exercise_slow_submit_does_not_block_get_status():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    first = manager.submit(_tmp_spec(tmp))
    backend.submit_gate = threading.Event()
    t = threading.Thread(target=manager.submit, args=(_tmp_spec(tmp),))
    t.start()
    assert backend.submit_entered.wait(5), "submit never reached the backend"
    try:
      assert _returns_within(lambda: manager.get_status(first.job_id), 2.0), \
        "get_status waited behind a slow submit"
      assert _returns_within(manager.list_jobs, 2.0), \
        "list_jobs waited behind a slow submit"
    finally:
      backend.submit_gate.set()
      t.join()
  finally:
    _finish(manager, backend)


def exercise_metadata_collisions_are_rejected():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    for key in ("state", "fake_handle"):
      try:
        manager.submit(_tmp_spec(tmp, metadata={key: 1}))
      except ValueError as exc:
        assert key in str(exc), str(exc)
      else:
        raise AssertionError("metadata key %r accepted" % key)
    assert manager.list_jobs() == []
  finally:
    _finish(manager, backend)


def exercise_unknown_job_raises_sorry():
  manager, backend = _make()
  try:
    for fn in (manager.get_status, manager.cancel,
               lambda jid: manager.wait_for_terminal(jid, 0.01)):
      try:
        fn("j_deadbeefdead")
      except Sorry as exc:
        assert "j_deadbeefdead" in str(exc), str(exc)
      else:
        raise AssertionError("unknown job id accepted")
  finally:
    _finish(manager, backend)


def exercise_describe_reports_backend_and_capacity():
  manager, backend = _make(capacity=3)
  try:
    assert manager.describe() == {"backend": "fake", "capacity": 3}
    assert manager.backend is backend
  finally:
    _finish(manager, backend)


def exercise_cancel_queued_job_without_handle():
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    c = manager.submit(_tmp_spec(tmp))
    d = manager.cancel(b.job_id)
    assert d["state"] == "cancelled" and d["fake_handle"] is None, d
    assert b.finished_at is not None and b._terminal_event.is_set()
    assert c.queue_position == 1, c.queue_position
    assert ("terminate", b.job_id) not in backend.calls
    # Cancel of a cancelled job: unchanged snapshot.
    assert manager.cancel(b.job_id)["state"] == "cancelled"
    manager.reaper_iter()
    assert ("submit", b.job_id) not in backend.calls, "cancelled job was submitted"
    assert a.state == "running"
  finally:
    _finish(manager, backend)


def exercise_cancel_running_job_terminates_then_kills_once():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(sigkill_after=10.0)
  try:
    job = manager.submit(_tmp_spec(tmp))
    t0 = time.monotonic()
    d = manager.cancel(job.job_id)
    assert d["state"] == "cancelling", d
    assert isinstance(job._cancel_sent_at, float)
    assert t0 <= job._cancel_sent_at <= time.monotonic()
    assert backend.calls[-1] == ("terminate", job.job_id), backend.calls
    # cancel of a cancelling job: no second terminate
    manager.cancel(job.job_id)
    assert backend.calls.count(("terminate", job.job_id)) == 1
    manager.reaper_iter()
    assert ("kill", job.job_id) not in backend.calls, "killed before the deadline"
    job._cancel_sent_at -= 20.0
    manager.reaper_iter()
    manager.reaper_iter()
    assert backend.calls.count(("kill", job.job_id)) == 1, backend.calls
    backend.results[job.handle] = PollResult("exited", exit_code=-9)
    manager.reaper_iter()
    d = manager.get_status(job.job_id)
    assert d["state"] == "cancelled" and d["exit_code"] == -9, d
  finally:
    _finish(manager, backend)


def exercise_cancel_of_a_job_that_already_exited_reports_its_exit():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    job = manager.submit(_tmp_spec(tmp))
    backend.results[job.handle] = PollResult("exited", exit_code=0)
    d = manager.cancel(job.job_id)
    assert d["state"] == "finished" and d["exit_code"] == 0, d
    assert ("terminate", job.job_id) not in backend.calls
    # A queued job whose poll says running gets started_at on cancel.
    backend.submit_phase = "pending"
    queued = manager.submit(_tmp_spec(tmp))
    assert queued.state == "queued" and queued.started_at is None
    backend.results[queued.handle] = PollResult("running")
    d = manager.cancel(queued.job_id)
    assert d["state"] == "cancelling" and d["started_at"] is not None, d
    assert backend.calls[-1] == ("terminate", queued.job_id)
  finally:
    _finish(manager, backend)


def exercise_cancelling_job_that_ended_on_its_own_uses_exit_code():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    job = manager.submit(_tmp_spec(tmp))
    manager.cancel(job.job_id)
    backend.results[job.handle] = PollResult("exited", exit_code=0, cancelled=False)
    manager.reaper_iter()
    assert job.state == "finished", job.state
    other = manager.submit(_tmp_spec(tmp))
    manager.cancel(other.job_id)
    backend.results[other.handle] = PollResult("exited", exit_code=0, cancelled=None)
    manager.reaper_iter()
    assert other.state == "cancelled", other.state
  finally:
    _finish(manager, backend)


def exercise_cancel_during_submission_is_honoured_afterwards():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  backend.submit_gate = threading.Event()
  try:
    holder = {}
    def go():
      holder["job"] = manager.submit(_tmp_spec(tmp))
    t = threading.Thread(target=go)
    t.start()
    assert backend.submit_entered.wait(5), "submit never reached the backend"
    job = manager.list_jobs()[0]
    assert job["state"] == "queued" and job["queue_position"] == 0, job
    d = manager.cancel(job["job_id"])
    assert d["state"] == "queued", d
    backend.submit_gate.set()
    t.join()
    job = holder["job"]
    assert job.state == "cancelling", job.state
    assert ("terminate", job.job_id) in backend.calls
    assert not job._cancel_requested
    # The same, but the submission fails: failed wins.
    backend.fail_submit_with = "refused"
    backend.submit_gate = threading.Event()
    backend.submit_entered.clear()
    t = threading.Thread(target=go)
    t.start()
    assert backend.submit_entered.wait(5), "submit never reached the backend"
    pending = [j for j in manager.list_jobs() if j["state"] == "queued"][0]
    manager.cancel(pending["job_id"])
    backend.submit_gate.set()
    t.join()
    d = manager.get_status(pending["job_id"])
    assert d["state"] == "failed" and "refused" in d["reason"], d
  finally:
    _finish(manager, backend)


def exercise_log_lines_are_written_with_the_lock_released():
  # SORRY lines are written outside the lock, so a stalled log mount never
  # stalls get_status; the terminal event waits for the line.
  import libtbx.jobs.manager as manager_module
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  backend.fail_submit_with = "refused"
  gate = threading.Event()
  entered = threading.Event()
  real_open = open
  def blocking_open(path, *args, **kwargs):
    if str(path).endswith("job.log"):
      entered.set()
      gate.wait(10)
    return real_open(path, *args, **kwargs)
  manager_module.open = blocking_open
  try:
    holder = {}
    t = threading.Thread(
      target=lambda: holder.update(job=manager.submit(_tmp_spec(tmp))))
    t.start()
    assert entered.wait(5), "the log line was never written"
    assert _returns_within(manager.list_jobs, 2.0), \
      "list_jobs waited behind a blocked log write"
    job = manager.list_jobs()[0]
    assert job["state"] == "failed", job
    live = [j for j in manager._jobs.values()][0]
    assert not live._terminal_event.is_set(), \
      "the terminal event fired before the log line was written"
    assert not _returns_within(
      lambda: manager.wait_for_terminal(job["job_id"], 5), 0.5), \
      "wait_for_terminal returned before the log line was written"
    gate.set()
    t.join()
    assert live._terminal_event.is_set()
    assert manager.wait_for_terminal(job["job_id"], 0)["state"] == "failed"
    with open(os.path.join(tmp, "job.log")) as fh:
      assert "SORRY: failed to submit job: refused" in fh.read()
  finally:
    del manager_module.open
    _finish(manager, backend)


def exercise_eviction_waits_for_the_log_line():
  # Eviction skips a job whose SORRY line is still being written, so a
  # tombstone never precedes the line.
  from libtbx.jobs.backends.base import PollResult
  import libtbx.jobs.manager as manager_module
  tmp = tempfile.mkdtemp()
  manager, backend = _make(max_retained_jobs=1)
  gate = threading.Event()
  entered = threading.Event()
  real_open = open
  def blocking_open(path, *args, **kwargs):
    if str(path).endswith("a.log"):
      entered.set()
      gate.wait(10)
    return real_open(path, *args, **kwargs)
  manager_module.open = blocking_open
  try:
    a = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "a.log")))
    b = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "b.log")))
    backend.results[a.handle] = PollResult("exited", exit_code=None,
                                           reason="NODE_FAIL")
    backend.results[b.handle] = PollResult("exited", exit_code=0)
    t = threading.Thread(target=manager.reaper_iter)
    t.start()
    assert entered.wait(5), "the reaper never wrote a's line"
    # Eviction from another thread mid-write skips a; b is within the cap.
    with manager._lock:
      manager._evict_terminal_locked()
    assert a.job_id in manager._jobs and a.job_id not in manager._evicted, \
      "a job was evicted before its log line was written"
    assert b.job_id in manager._jobs and not manager._evicted
    gate.set()
    t.join()
    assert a._terminal_event.is_set()
    manager.reaper_iter()
    assert a.job_id in manager._evicted and a.job_id not in manager._jobs
    assert b.job_id in manager._jobs, "the newer job was evicted instead"
    with open(os.path.join(tmp, "a.log")) as fh:
      assert "NODE_FAIL" in fh.read()
  finally:
    del manager_module.open
    _finish(manager, backend)


def exercise_poll_failures_are_reported_once_and_forgotten_on_eviction():
  from libtbx.jobs.backends.base import PollResult
  import io
  tmp = tempfile.mkdtemp()
  manager, backend = _make(max_retained_jobs=1)
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    a = manager.submit(_tmp_spec(tmp))
    backend.results[a.handle] = RuntimeError("poll broke")
    manager.reaper_iter()
    manager.reaper_iter()
    assert a.state == "running", a.state
    assert err.getvalue().count("poll failed for %s" % a.job_id) == 1
    assert a.job_id in manager._poll_errors_reported
    # A result that is not a PollResult is a poll failure too, and the
    # pass still submits the FIFO jobs it marked.
    backend.cap = 1
    c = manager.submit(_tmp_spec(tmp))
    assert c.state == "queued", c.state
    backend.results[a.handle] = "done"
    manager.reaper_iter()
    assert a.state == "running", a.state
    assert err.getvalue().count("poll failed for %s" % a.job_id) == 1, \
      "a second kind of poll failure on the same job was reported again"
    fresh = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "f.log")))
    assert fresh.state == "queued", "capacity is 1"
    backend.results[a.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    assert a.state == "finished" and c.state == "running", (a.state, c.state)
    backend.results[c.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    assert fresh.state == "running", fresh.state
    backend.results[fresh.handle] = "done"
    manager.reaper_iter()
    assert err.getvalue().count("poll failed for %s" % fresh.job_id) == 1
    assert "not a PollResult" in err.getvalue(), err.getvalue()
    backend.results[fresh.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    manager.reaper_iter()
    assert a.job_id in manager._evicted, "a was not evicted"
    assert a.job_id not in manager._poll_errors_reported
  finally:
    sys.stderr = real
    _finish(manager, backend)


def exercise_the_cap_line_holds_the_terminal_event():
  # A cap line still in flight when the job exits holds back both the
  # terminal event and eviction.
  from libtbx.jobs.backends.base import PollResult
  import libtbx.jobs.manager as manager_module
  tmp = tempfile.mkdtemp()
  manager, backend = _make(max_job_seconds=1, max_retained_jobs=1)
  gate = threading.Event()
  entered = threading.Event()
  real_open = open
  def blocking_open(path, *args, **kwargs):
    if str(path).endswith("capped.log"):
      entered.set()
      gate.wait(10)
    return real_open(path, *args, **kwargs)
  manager_module.open = blocking_open
  try:
    job = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "capped.log")))
    other = manager.submit(_tmp_spec(tmp))
    job._started_monotonic -= 5.0
    t = threading.Thread(target=manager.reaper_iter)   # queues the cap line
    t.start()
    assert entered.wait(5), "the cap line was never written"
    assert job.state == "cancelling" and job._pending_lines == 1
    # The job exits on another pass while the cap line is in flight.
    backend.results[job.handle] = PollResult("exited", exit_code=-15,
                                             cancelled=True)
    backend.results[other.handle] = PollResult("exited", exit_code=0)
    with manager._lock:
      manager._apply_poll_locked(job, backend.results[job.handle],
                                 manager_module._utcnow(), time.monotonic())
      manager._apply_poll_locked(other, backend.results[other.handle],
                                 manager_module._utcnow(), time.monotonic())
      manager._evict_terminal_locked()
    assert job.state == "cancelled" and not job._terminal_event.is_set(), \
      "the event fired with the cap line still in flight"
    assert job.job_id in manager._jobs, "evicted with a line in flight"
    gate.set()
    t.join()
    assert job._terminal_event.is_set() and job._pending_lines == 0
    with open(os.path.join(tmp, "capped.log")) as fh:
      assert "wall-clock cap" in fh.read()
  finally:
    del manager_module.open
    _finish(manager, backend)


def exercise_a_job_wakes_when_its_own_line_is_written():
  # Two jobs' lines in one flush: the first job's waiters wake as soon
  # as its line is written, not after the second job's slow write.
  from libtbx.jobs.backends.base import PollResult
  import libtbx.jobs.manager as manager_module
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  gate = threading.Event()
  entered = threading.Event()
  real_open = open
  def blocking_open(path, *args, **kwargs):
    if str(path).endswith("slow.log"):
      entered.set()
      gate.wait(10)
    return real_open(path, *args, **kwargs)
  manager_module.open = blocking_open
  try:
    fast = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "fast.log")))
    slow = manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "slow.log")))
    for job in (fast, slow):
      backend.results[job.handle] = PollResult("exited", exit_code=None,
                                               reason="NODE_FAIL")
    t = threading.Thread(target=manager.reaper_iter)
    t.start()
    assert entered.wait(5), "the slow write was never reached"
    assert fast._terminal_event.is_set(), \
      "the first job waited for the second job's write"
    assert _returns_within(lambda: manager.wait_for_terminal(fast.job_id, 5), 2.0)
    assert not slow._terminal_event.is_set()
    gate.set()
    t.join()
    assert slow._terminal_event.is_set()
  finally:
    del manager_module.open
    _finish(manager, backend)


def exercise_an_interrupted_flush_keeps_the_lines_owed():
  import libtbx.jobs.manager as manager_module
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  backend.fail_submit_with = "refused"
  strikes = []
  real_open = open
  def interrupting_open(path, *args, **kwargs):
    if str(path).endswith("job.log") and not strikes:
      strikes.append(path)
      raise KeyboardInterrupt()
    return real_open(path, *args, **kwargs)
  manager_module.open = interrupting_open
  try:
    try:
      job = manager.submit(_tmp_spec(tmp))
    except KeyboardInterrupt:
      job = manager.list_jobs()[0]
      job = manager._jobs[job["job_id"]]
    assert job.state == "failed" and job._pending_lines == 1, job._pending_lines
    assert not job._terminal_event.is_set()
    assert manager._log_queue and manager._log_queue[0][0] is job, \
      "the interrupted line was not put back on the queue"
    manager._flush_logs()
    assert job._terminal_event.is_set() and job._pending_lines == 0
    with open(os.path.join(tmp, "job.log")) as fh:
      assert "SORRY: failed to submit job: refused" in fh.read()
  finally:
    del manager_module.open
    _finish(manager, backend)


def exercise_a_broken_stderr_does_not_stop_the_reaper():
  # Every stderr note is guarded: with stderr broken, a poll failure, a
  # backend failure and a shutdown still complete.
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  class Broken(object):
    attempts = 0
    def write(self, text):
      Broken.attempts += 1
      raise OSError("Broken pipe")
    def flush(self):
      raise OSError("Broken pipe")
  real = sys.stderr
  sys.stderr = Broken()
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    backend.results[a.handle] = RuntimeError("poll broke")
    manager.reaper_iter()                       # the poll note
    backend.results[a.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    assert a.state == "finished" and b.state == "running", (a.state, b.state)
    real_terminate = backend.terminate
    backend.terminate = lambda job: (_ for _ in ()).throw(RuntimeError("no"))
    manager.cancel(b.job_id)                    # the backend-failed note
    backend.terminate = real_terminate
    assert b.state == "cancelling"
    backend.on_shutdown = lambda jobs: (_ for _ in ()).throw(RuntimeError("no"))
    manager.shutdown()                          # never raises
    # A refresh failure in a reaper pass and in adopt_many, a submit
    # that raises, and a malformed poll result: every note guarded.
    backend.refresh_hook = lambda jobs: (_ for _ in ()).throw(RuntimeError("no"))
    manager.reaper_iter()
    adopted = manager.adopt(_tmp_spec(tmp), "abc")   # its refresh note is guarded
    assert adopted.state == "running", adopted.state
    backend.refresh_hook = None
    backend.results[adopted.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()                       # frees the slot it took
    assert adopted.state == "finished", adopted.state
    backend.results[b.handle] = PollResult("exited", exit_code=-15, cancelled=True)
    manager.reaper_iter()                       # frees the single slot
    assert b.state == "cancelled", b.state
    backend.crash_submit = True
    crashed = manager.submit(_tmp_spec(tmp))
    backend.crash_submit = False
    assert crashed.state == "failed", crashed.state
    d = manager.submit(_tmp_spec(tmp))
    assert d.state == "running", d.state
    backend.results[d.handle] = "not a result"
    manager.reaper_iter()
    # With stderr missing altogether nothing may fall back to stdout,
    # which for a stdio server is the protocol channel.
    import io
    sys.stderr = None
    out = io.StringIO()
    real_out = sys.stdout
    sys.stdout = out
    try:
      backend.results[d.handle] = RuntimeError("poll broke again")
      manager._poll_errors_reported.discard(d.job_id)
      manager.reaper_iter()
    finally:
      sys.stdout = real_out
    assert out.getvalue() == "", out.getvalue()
    # The reaper thread itself survives a failed note.
    sys.stderr = Broken()
    live_manager, live_backend = _make(reaper_interval=0.05)
    try:
      c = live_manager.submit(_tmp_spec(tmp, log_path=os.path.join(tmp, "c.log")))
      before = Broken.attempts
      live_backend.results[c.handle] = RuntimeError("poll broke")
      deadline = time.time() + 5
      while Broken.attempts == before and time.time() < deadline:
        time.sleep(0.02)
      assert Broken.attempts > before, "the reaper thread never tried the note"
      time.sleep(0.2)
      assert live_manager._reaper_thread.is_alive(), "the reaper thread died"
      live_backend.results[c.handle] = PollResult("exited", exit_code=0)
      assert live_manager.wait_for_terminal(c.job_id, 5)["state"] == "finished"
    finally:
      live_manager.stop(timeout=5)
      assert not live_backend.violations, live_backend.violations
  finally:
    sys.stderr = real
    _finish(manager, backend)


def exercise_the_reaper_submits_marked_jobs_whatever_else_raises():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    backend.results[a.handle] = PollResult("exited", exit_code=0)
    real_flush = manager._flush_logs
    def broken_flush():
      manager._flush_logs = real_flush
      raise RuntimeError("flush broke")
    manager._flush_logs = broken_flush
    try:
      manager.reaper_iter()
    except RuntimeError:
      pass
    else:
      raise AssertionError("the injected failure was swallowed")
    assert b.state == "running" and b.handle == "h2", (b.state, b.handle)
    # A raise inside eviction leaves the pass's terminate and kill intact.
    manager.cancel(b.job_id)
    b._cancel_sent_at -= 20.0
    real_evict = manager._evict_terminal_locked
    def broken_evict():
      manager._evict_terminal_locked = real_evict
      raise RuntimeError("evict broke")
    manager._evict_terminal_locked = broken_evict
    import io
    err = io.StringIO()
    real = sys.stderr
    sys.stderr = err
    try:
      manager.reaper_iter()
    finally:
      sys.stderr = real
    assert ("kill", b.job_id) in backend.calls, backend.calls
    assert "evict broke" in err.getvalue()
    # Even a raise that escapes the pass (the flush again) lets the
    # terminate and kill calls out before the FIFO is drained.
    c = manager.submit(_tmp_spec(tmp))
    d = manager.submit(_tmp_spec(tmp))
    assert (c.state, d.state) == ("queued", "queued"), (c.state, d.state)
    backend.results[b.handle] = PollResult("exited", exit_code=-9)
    manager.reaper_iter()
    assert b.state == "cancelled" and c.state == "running", (b.state, c.state)
    manager.cancel(c.job_id)
    c._cancel_sent_at -= 20.0
    manager._flush_logs = broken_flush
    try:
      manager.reaper_iter()
    except RuntimeError:
      pass
    assert ("kill", c.job_id) in backend.calls, backend.calls
  finally:
    _finish(manager, backend)


def exercise_a_failing_log_write_still_wakes_waiters_and_frees_slots():
  from libtbx.jobs.backends.base import PollResult
  import libtbx.jobs.manager as manager_module
  import io
  tmp = tempfile.mkdtemp()
  manager, backend = _make(capacity=1)
  real_open = open
  def broken_open(path, *args, **kwargs):
    if str(path).endswith("job.log"):
      raise ValueError("embedded null byte")
    return real_open(path, *args, **kwargs)
  manager_module.open = broken_open
  err = io.StringIO()
  real_err = sys.stderr
  sys.stderr = err
  try:
    backend.fail_submit_with = "refused"
    job = manager.submit(_tmp_spec(tmp))        # must not raise
    assert job.state == "failed" and job._terminal_event.is_set()
    assert manager.wait_for_terminal(job.job_id, 5)["state"] == "failed"
    backend.fail_submit_with = None
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    assert (a.state, b.state) == ("running", "queued")
    backend.results[a.handle] = PollResult("exited", exit_code=None,
                                           reason="NODE_FAIL")
    manager.reaper_iter()                        # queues a line for a
    assert a.state == "failed" and a._terminal_event.is_set()
    assert b.state == "running", "a failed log write stranded the FIFO"
  finally:
    sys.stderr = real_err
    del manager_module.open
    _finish(manager, backend)
  assert "embedded null byte" in err.getvalue(), err.getvalue()


def exercise_stop_is_checked_under_the_lock():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  # A submit that passed the early check and is reading the capacity when
  # stop() lands must not register a job nothing will ever drain.
  backend.capacity_gate = threading.Event()
  outcome = {}
  def go():
    try:
      outcome["job"] = manager.submit(_tmp_spec(tmp))
    except RuntimeError as exc:
      outcome["error"] = str(exc)
  t = threading.Thread(target=go)
  t.start()
  assert backend.capacity_entered.wait(5), "submit never read the capacity"
  manager.stop()
  backend.capacity_gate.set()
  t.join()
  assert "stop" in outcome.get("error", ""), outcome
  assert manager.list_jobs() == [], manager.list_jobs()
  try:
    manager.adopt(_tmp_spec(tmp), "late")
  except RuntimeError:
    pass
  else:
    raise AssertionError("adopt after stop accepted")
  assert not backend.violations, backend.violations


def exercise_a_job_made_terminal_during_refresh_is_not_polled():
  # The reaper releases the lock for refresh; a cancel that ends the job
  # meanwhile must keep it from being polled afterwards.
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    job = manager.submit(_tmp_spec(tmp))
    backend.results[job.handle] = PollResult("exited", exit_code=0)
    def cancel_during_refresh(jobs):
      backend.refresh_hook = None
      d = manager.cancel(job.job_id)
      assert d["state"] == "finished", d
      backend.results[job.handle] = PollResult("running")
      backend.polled = []
    backend.refresh_hook = cancel_during_refresh
    manager.reaper_iter()
    assert backend.polled == [], "the reaper polled a job that had ended"
    assert job.state == "finished" and job.exit_code == 0, (job.state, job.exit_code)
  finally:
    _finish(manager, backend)


def exercise_wall_clock_cap_cancels_and_leaves_short_jobs_alone():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(max_job_seconds=1)
  try:
    job = manager.submit(_tmp_spec(tmp))
    manager.reaper_iter()
    assert job.state == "running"
    job._started_monotonic -= 2.0
    manager.reaper_iter()
    assert job.state == "cancelling" and job.reason == "wall-clock cap", (job.state, job.reason)
    assert backend.calls[-1] == ("terminate", job.job_id)
    with open(os.path.join(tmp, "job.log")) as fh:
      log = fh.read()
    assert "exceeded the wall-clock cap (1 s)" in log, log
    backend.results[job.handle] = PollResult(
      "exited", exit_code=-15, reason="killed by signal 15")
    manager.reaper_iter()
    assert job.state == "cancelled", job.state
    assert job.reason == "wall-clock cap (killed by signal 15)", job.reason
  finally:
    _finish(manager, backend)


def exercise_wait_for_terminal_and_eviction_race():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make(max_retained_jobs=1, max_evicted_tombstones=5)
  try:
    a = manager.submit(_tmp_spec(tmp, metadata={"args": ["a"]}))
    assert manager.wait_for_terminal(a.job_id, 0.05)["state"] == "running"
    backend.results[a.handle] = PollResult("exited", exit_code=0)
    def finish_later():
      time.sleep(0.1)
      manager.reaper_iter()
    worker = threading.Thread(target=finish_later)
    worker.start()
    d = manager.wait_for_terminal(a.job_id, 5)
    worker.join()
    assert d["state"] == "finished", d
    # A waiter whose job is evicted in the very pass that ends it still
    # gets the job's final status, not "Unknown job".
    b = manager.submit(_tmp_spec(tmp))
    c = manager.submit(_tmp_spec(tmp))
    for j in (b, c):
      backend.results[j.handle] = PollResult("exited", exit_code=1)
    result = {}
    class GatedEvent(threading.Event):
      entered = threading.Event()
      def wait(self, timeout=None):
        self.entered.set()
        return threading.Event.wait(self, timeout)
    b._terminal_event = GatedEvent()
    waiter = threading.Thread(
      target=lambda: result.update(manager.wait_for_terminal(b.job_id, 5)))
    waiter.start()
    assert GatedEvent.entered.wait(5), "the waiter never reached the event"
    manager.reaper_iter()
    waiter.join(5)
    assert not waiter.is_alive(), "wait_for_terminal did not wake"
    assert a.job_id not in manager._jobs and b.job_id not in manager._jobs
    assert result["state"] == "failed" and result["exit_code"] == 1, result
    assert manager.get_status(a.job_id)["state"] == "finished"
    assert manager.wait_for_terminal(a.job_id, 0.01)["state"] == "finished"
    manager.get_status(a.job_id)["args"].append("mutated")
    assert manager.get_status(a.job_id)["args"] == ["a"], "tombstone shared"
    assert manager.cancel(c.job_id)["state"] == "failed"
    assert len(manager.list_jobs()) == 1
  finally:
    _finish(manager, backend)


def exercise_adopt_and_adopt_many():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    backend.results["re-abc"] = PollResult("running")
    backend.results["re-def"] = PollResult("exited", exit_code=0)
    backend.results["re-ghi"] = PollResult("pending", reason="PENDING")
    jobs = manager.adopt_many([
      (_tmp_spec(tmp), "abc"),
      (_tmp_spec(tmp), "def", "j_000000000001"),
      (_tmp_spec(tmp), "ghi"),
    ])
    assert [j.state for j in jobs] == ["running", "finished", "queued"], \
      [j.state for j in jobs]
    assert jobs[1].job_id == "j_000000000001"
    assert jobs[0].started_at is not None and jobs[0].handle == "re-abc"
    refreshes = [c for c in backend.calls if c[0] == "refresh"]
    assert len(refreshes) == 1 and refreshes[0][2] is True, refreshes
    assert set(refreshes[0][1]) == set(j.job_id for j in jobs)
    one = manager.adopt(_tmp_spec(tmp), "abc")
    assert one.state == "running"
    # A bad handle adopts nothing, even when the good ones come first.
    backend.reattach_bad = {"bad": Sorry("bad handle")}
    before = len(manager.list_jobs())
    try:
      manager.adopt_many([(_tmp_spec(tmp), "ok"), (_tmp_spec(tmp), "ok2"),
                          (_tmp_spec(tmp), "bad")])
    except Sorry as exc:
      assert "bad handle" in str(exc), str(exc)
    else:
      raise AssertionError("adopt_many did not raise")
    assert len(manager.list_jobs()) == before
    assert not [c for c in backend.calls if c[0] == "refresh"
                and set(c[1]) - {j.job_id for j in jobs} - {one.job_id}], \
      "a refresh ran for a batch that adopted nothing"
    backend.reattach_bad = {}
    try:
      manager.adopt(_tmp_spec(tmp), "x", job_id=one.job_id)
    except Sorry as exc:
      assert one.job_id in str(exc)
    else:
      raise AssertionError("duplicate job id accepted")
    before = len(manager.list_jobs())
    try:
      manager.adopt_many([(_tmp_spec(tmp), "p", "j_00000000dup1"),
                          (_tmp_spec(tmp), "q", "j_00000000dup1")])
    except Sorry as exc:
      assert "j_00000000dup1" in str(exc), str(exc)
    else:
      raise AssertionError("the same job id twice in one batch accepted")
    assert len(manager.list_jobs()) == before
  finally:
    _finish(manager, backend)


def exercise_shutdown_passes_live_jobs_and_never_blocks():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    a = manager.submit(_tmp_spec(tmp))
    b = manager.submit(_tmp_spec(tmp))
    backend.results[b.handle] = PollResult("exited", exit_code=0)
    manager.reaper_iter()
    manager._lock.acquire()
    try:
      manager.shutdown()
    finally:
      manager._lock.release()
    assert backend.shutdown_seen == [a.job_id], backend.shutdown_seen
  finally:
    _finish(manager, backend)


def exercise_reaper_thread_starts_once_and_recovers_from_start_failure():
  tmp = tempfile.mkdtemp()
  manager, backend = _make()
  try:
    manager.submit(_tmp_spec(tmp))
    thread = manager._reaper_thread
    manager.submit(_tmp_spec(tmp))
    assert manager._reaper_thread is thread and thread.is_alive()
  finally:
    manager.stop(timeout=5)
    assert not backend.violations, backend.violations
  manager2, backend2 = _make()
  real_start = threading.Thread.start
  def boom(self):
    raise OSError("no threads")
  threading.Thread.start = boom
  try:
    try:
      manager2.submit(_tmp_spec(tmp))
    except OSError:
      pass
    else:
      raise AssertionError("thread start failure was swallowed")
    assert manager2._reaper_started is False
  finally:
    threading.Thread.start = real_start
  job = manager2.submit(_tmp_spec(tmp))
  assert manager2._reaper_started is True and job.state == "running"
  manager2.stop(timeout=5)
  assert not backend2.violations, backend2.violations
  # stop is final: new work is refused, existing jobs stay readable.
  for fn in (lambda: manager2.submit(_tmp_spec(tmp)),
             lambda: manager2.adopt(_tmp_spec(tmp), "x")):
    try:
      fn()
    except RuntimeError as exc:
      assert "stop" in str(exc), str(exc)
    else:
      raise AssertionError("a stopped manager accepted new work")
  assert manager2.get_status(job.job_id)["state"] == "running"
  manager2.stop()


def run():
  if sys.platform == "win32":
    print("Skipping tst_manager on Windows")
    print("OK")
    return
  exercise_poll_result_validates_phase()
  exercise_backend_defaults()
  exercise_jobspec_validation()
  exercise_job_to_dict_layers_core_handle_and_metadata()
  exercise_submit_runs_immediately_and_reports_running()
  exercise_capacity_queues_and_drains_in_order()
  exercise_capacity_none_mirrors_scheduler_pending()
  exercise_external_cancellation_and_signal_exit()
  exercise_submit_failure_is_a_failed_job_not_an_exception()
  exercise_backend_submit_crash_is_a_failed_job()
  exercise_unknown_submit_phase_fails_and_terminates_the_job()
  exercise_an_interrupt_inside_submit_frees_the_slot()
  exercise_slow_submit_does_not_block_get_status()
  exercise_metadata_collisions_are_rejected()
  exercise_unknown_job_raises_sorry()
  exercise_describe_reports_backend_and_capacity()
  exercise_cancel_queued_job_without_handle()
  exercise_cancel_running_job_terminates_then_kills_once()
  exercise_cancel_of_a_job_that_already_exited_reports_its_exit()
  exercise_cancelling_job_that_ended_on_its_own_uses_exit_code()
  exercise_cancel_during_submission_is_honoured_afterwards()
  exercise_a_job_made_terminal_during_refresh_is_not_polled()
  exercise_log_lines_are_written_with_the_lock_released()
  exercise_a_failing_log_write_still_wakes_waiters_and_frees_slots()
  exercise_eviction_waits_for_the_log_line()
  exercise_poll_failures_are_reported_once_and_forgotten_on_eviction()
  exercise_the_cap_line_holds_the_terminal_event()
  exercise_a_job_wakes_when_its_own_line_is_written()
  exercise_an_interrupted_flush_keeps_the_lines_owed()
  exercise_the_reaper_submits_marked_jobs_whatever_else_raises()
  exercise_a_broken_stderr_does_not_stop_the_reaper()
  exercise_stop_is_checked_under_the_lock()
  exercise_wall_clock_cap_cancels_and_leaves_short_jobs_alone()
  exercise_wait_for_terminal_and_eviction_race()
  exercise_adopt_and_adopt_many()
  exercise_shutdown_passes_live_jobs_and_never_blocks()
  exercise_reaper_thread_starts_once_and_recovers_from_start_failure()
  print("OK")


if __name__ == "__main__":
  run()
