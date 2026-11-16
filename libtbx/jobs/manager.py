"""Job data model and the ``JobManager`` state machine.

The contract is spec sections 4 and 5
(``docs/superpowers/specs/2026-09-24-libtbx-job-manager-design.md`` in
the phenix repository).

Do not interrupt a thread inside the manager with an exception-raising
signal: the signal can land while a lock is held but outside its
protected region. A SIGINT handler should set a flag; the host then
calls ``stop`` and ``shutdown`` from normal control flow.
"""

import contextlib
import copy
import json
import os
import re
import sys
import threading
import time
import traceback
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone

from libtbx.jobs.backends.base import PollResult, SubmitError
from libtbx.utils import Sorry


TERMINAL_STATES = frozenset(("finished", "failed", "cancelled"))
STATES = ("queued", "running", "cancelling", "finished", "failed", "cancelled")
CORE_STATUS_KEYS = ("job_id", "state", "queue_position", "exit_code",
                    "reason", "log_path", "started_at", "finished_at")
RESOURCE_KEYS = ("cpus", "mem_mb", "time_minutes", "gpus")
_RESOURCE_MINIMUM = {"cpus": 1, "mem_mb": 1, "time_minutes": 1, "gpus": 0}
PLACEHOLDER_VALUE_RE = re.compile(r"^[A-Za-z0-9_.,:/=+@-]{1,256}$")
ENV_KEY_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _utcnow():
  """Current UTC time as a tz-naive ``datetime``.

  Tz-naive so status-dict timestamps carry no ``+00:00`` suffix.
  """
  return datetime.now(timezone.utc).replace(tzinfo=None)


def _new_job_id():
  """Return a fresh id of the form ``j_<12 hex>``."""
  return "j_" + uuid.uuid4().hex[:12]


def validate_placeholder_value(key, value):
  """Check one free-form placeholder value (spec 4.1).

  Parameters
  ----------
  key : str
      The placeholder name, for the error message.
  value : object
      An int, a float, or a str of 1-256 characters drawn from
      ``A-Z a-z 0-9 _ . , : / = + @ -``.

  Raises
  ------
  ValueError
      For any other value, including ``bool``.
  """
  if isinstance(value, bool):
    raise ValueError("placeholder %r must not be a bool" % key)
  if isinstance(value, (int, float)):
    return
  if isinstance(value, str) and PLACEHOLDER_VALUE_RE.fullmatch(value):
    return
  raise ValueError(
    "placeholder %r must be an int, a float, or a string of 1-256 "
    "characters from A-Z a-z 0-9 _ . , : / = + @ - (got %r)" % (key, value))


@dataclass(frozen=True)
class JobSpec:
  """What to run.

  Parameters
  ----------
  argv : list of str
      The command; non-empty.
  cwd : str
      Absolute working directory; the caller ensures it exists.
  log_path : str
      Absolute path; the job's stdout and stderr are appended here.
  env : dict, optional
      Environment overrides applied on top of the inherited environment.
  name : str, optional
      Short label for scheduler job names; empty becomes ``"job"``.
  resources : dict, optional
      ``cpus``, ``mem_mb``, ``time_minutes``, ``gpus`` as ints, plus
      free-form placeholder values for batch templates.
  metadata : dict, optional
      Host-owned, JSON-serialisable; merged into the status dict.
  """

  argv: list
  cwd: str
  log_path: str
  env: dict = field(default_factory=dict)
  name: str = "job"
  resources: dict = field(default_factory=dict)
  metadata: dict = field(default_factory=dict)

  def __post_init__(self):
    if isinstance(self.argv, str) or not isinstance(self.argv, (list, tuple)) \
       or not self.argv or not all(isinstance(a, str) for a in self.argv):
      raise ValueError("argv must be a non-empty list of str")
    object.__setattr__(self, "argv", list(self.argv))
    for label in ("cwd", "log_path"):
      value = getattr(self, label)
      if not isinstance(value, str) or not os.path.isabs(value):
        raise ValueError("%s must be an absolute path (got %r)" % (label, value))
    if not isinstance(self.env, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in self.env.items()):
      raise ValueError("env must map str to str")
    for label, values in (("argv", self.argv), ("cwd", [self.cwd]),
                          ("log_path", [self.log_path]),
                          ("env", self.env.values())):
      if any("\0" in v for v in values):
        raise ValueError("%s must not contain a NUL character" % label)
    for key in self.env:
      if not ENV_KEY_RE.fullmatch(key):
        raise ValueError(
          "env key %r is not a valid environment variable name" % key)
    object.__setattr__(self, "env", dict(self.env))
    if not isinstance(self.name, str):
      raise ValueError("name must be a str")
    if not self.name:
      object.__setattr__(self, "name", "job")
    if not isinstance(self.resources, dict):
      raise ValueError("resources must be a dict")
    for key, value in self.resources.items():
      if not isinstance(key, str):
        raise ValueError("resources keys must be str (got %r)" % (key,))
      if key in RESOURCE_KEYS:
        if isinstance(value, bool) or not isinstance(value, int) \
           or value < _RESOURCE_MINIMUM[key]:
          raise ValueError(
            "resources[%r] must be an int >= %d (got %r)"
            % (key, _RESOURCE_MINIMUM[key], value))
      else:
        validate_placeholder_value(key, value)
    object.__setattr__(self, "resources", dict(self.resources))
    if not isinstance(self.metadata, dict) or not all(
        isinstance(k, str) for k in self.metadata):
      raise ValueError("metadata must be a dict with str keys")
    try:
      json.dumps(self.metadata)
    except (TypeError, ValueError) as exc:
      raise ValueError("metadata must be JSON-serialisable: %s" % exc)
    object.__setattr__(self, "metadata", copy.deepcopy(self.metadata))


@dataclass
class Job:
  """A job the manager owns; mutated only under the manager lock.

  Public fields follow spec 4.2; underscore fields are bookkeeping.
  ``_cancel_sent_at`` (a ``time.monotonic()`` float) is read by an MCP
  test, so keep its name and type.
  """

  job_id: str
  spec: JobSpec
  state: str = "queued"
  queue_position: int = 0
  submitted_at: datetime = None
  started_at: datetime = None
  finished_at: datetime = None
  exit_code: int = None
  reason: str = None
  handle: object = None
  _submitting: bool = False
  _cancel_requested: bool = False
  _cancel_sent_at: float = None
  _started_monotonic: float = None
  _kill_sent: bool = False
  _cancel_reason: str = None
  _pending_lines: int = 0
  _terminal_event: threading.Event = field(
    init=False, default_factory=threading.Event, repr=False, compare=False)

  def to_dict(self, backend):
    """Serialise to the status dict of spec 4.4."""
    d = {
      "job_id": self.job_id,
      "state": self.state,
      "queue_position": self.queue_position,
      "exit_code": self.exit_code,
      "reason": self.reason,
      "log_path": self.spec.log_path,
      "started_at": self.started_at.isoformat() if self.started_at else None,
      "finished_at": self.finished_at.isoformat() if self.finished_at else None,
    }
    d.update(backend.handle_info(self))
    d.update(copy.deepcopy(self.spec.metadata))
    return d


class JobManager(object):
  """Owns jobs, the FIFO and the reaper; drives every state transition.

  Parameters
  ----------
  backend : Backend
      The execution backend.
  max_retained_jobs : int, optional
      Terminal jobs kept in the live table before eviction (default 200).
  max_evicted_tombstones : int, optional
      Final status dicts kept for evicted jobs so ``get_status`` still
      answers for them (default 2000; 0 disables).
  max_job_seconds : int, optional
      Wall-clock cap after which a running job is cancelled; 0 = off.
  reaper_interval : float, optional
      Seconds between reaper passes (default 0.5).
  sigkill_after : float, optional
      Seconds after a cancel before the backend's ``kill`` (default 10).
  """

  def __init__(self, backend, max_retained_jobs=200,
               max_evicted_tombstones=2000, max_job_seconds=0,
               reaper_interval=0.5, sigkill_after=10.0):
    self.backend = backend
    self.max_retained_jobs = max(1, int(max_retained_jobs))
    self.max_evicted_tombstones = max(0, int(max_evicted_tombstones))
    self.max_job_seconds = max(0, int(max_job_seconds))
    self.reaper_interval = float(reaper_interval)
    self.sigkill_after = float(sigkill_after)
    self._lock = threading.Lock()
    self._lock_owner = None
    self._jobs = {}
    self._evicted = {}
    self._queue = []
    self._reaper_started = False
    self._reaper_thread = None
    self._stop_event = threading.Event()
    self._poll_errors_reported = set()
    self._log_queue = []

  # ---- lock helpers --------------------------------------------------------

  @contextlib.contextmanager
  def _owned(self):
    """Record the lock owner for ``holds_lock`` (diagnostics, tests).

    Used as ``with self._lock, self._owned():`` so it runs inside the
    lock's own context manager and an exception here still releases it.
    """
    self._lock_owner = threading.get_ident()
    try:
      yield
    finally:
      self._lock_owner = None

  def holds_lock(self):
    """True when the calling thread holds the manager lock."""
    return self._lock_owner == threading.get_ident()

  # ---- public API ------------------------------------------------------------

  def describe(self):
    """Return ``{"backend": name, "capacity": capacity}``."""
    return {"backend": self.backend.name, "capacity": self.backend.capacity()}

  def submit(self, spec):
    """Create a job and submit it now or queue it (spec 5.3).

    Returns
    -------
    Job
        The new job, ``running``, ``queued`` (scheduler-side or in the
        FIFO) or ``failed`` when the backend refused it.
    """
    self._prepare_spec(spec)
    job = Job(job_id=_new_job_id(), spec=spec, submitted_at=_utcnow())
    self._ensure_reaper_started()
    capacity = self.backend.capacity()
    with self._lock, self._owned():
      self._check_not_stopped()
      self._jobs[job.job_id] = job
      if not self._queue and self._slots_free_locked(capacity):
        # _submitting holds the slot while backend.submit runs unlocked,
        # so a concurrent submit or reaper pass cannot overfill capacity.
        job._submitting = True
        immediate = True
      else:
        self._queue.append(job.job_id)
        job.queue_position = len(self._queue)
        immediate = False
    if immediate:
      self._submit_job(job)
    return job

  def get_status(self, job_id):
    """Status dict for a live or evicted job; ``Sorry`` when unknown."""
    with self._lock, self._owned():
      job = self._jobs.get(job_id)
      if job is not None:
        return job.to_dict(self.backend)
      tomb = self._evicted.get(job_id)
      if tomb is not None:
        return copy.deepcopy(tomb)
      raise Sorry("Unknown job: %s" % job_id)

  def list_jobs(self):
    """Status dicts for every live job."""
    with self._lock, self._owned():
      return [j.to_dict(self.backend) for j in self._jobs.values()]

  def cancel(self, job_id):
    """Request cancellation (spec 4.3, 5.3); returns the status dict."""
    need_terminate = False
    with self._lock, self._owned():
      job = self._jobs.get(job_id)
      if job is None:
        raise Sorry("Unknown job: %s" % job_id)
      if job.state in TERMINAL_STATES or job.state == "cancelling":
        return job.to_dict(self.backend)
      if job._submitting:
        job._cancel_requested = True
        return job.to_dict(self.backend)
      if job.handle is None:
        if job_id in self._queue:
          self._queue.remove(job_id)
        job.state = "cancelled"
        job.queue_position = 0
        job.finished_at = _utcnow()
        job._terminal_event.set()
        self._renumber_queue_locked()
        return job.to_dict(self.backend)
      result = self._safe_poll_locked(job)
      if result is not None:
        self._apply_poll_locked(job, result, _utcnow(), time.monotonic())
      if job.state in TERMINAL_STATES:
        snapshot = job.to_dict(self.backend)
      else:
        self._begin_cancel_locked(job)
        need_terminate = True
    if not need_terminate:
      self._flush_logs()
      return snapshot
    self._call_backend("terminate", job)
    with self._lock, self._owned():
      return job.to_dict(self.backend)

  def wait_for_terminal(self, job_id, timeout):
    """Block until the job is done or ``timeout`` seconds pass.

    Done means terminal with every queued log line written, so the
    caller can read the log on return. After a timeout the state may be
    terminal while a line is still being written.
    """
    with self._lock, self._owned():
      job = self._jobs.get(job_id)
      if job is None:
        tomb = self._evicted.get(job_id)
        if tomb is not None:
          return copy.deepcopy(tomb)
        raise Sorry("Unknown job: %s" % job_id)
      event = job._terminal_event
      # Check the event, not the state: it is set only after the job's
      # SORRY line is written (spec 5.2).
      if event.is_set():
        return job.to_dict(self.backend)
    event.wait(timeout=timeout)
    with self._lock, self._owned():
      return job.to_dict(self.backend)

  def adopt(self, spec, handle, job_id=None):
    """Recover one job the backend already knows; see ``adopt_many``."""
    return self.adopt_many([(spec, handle, job_id)])[0]

  def adopt_many(self, entries):
    """Recover jobs the backend already knows, with one scheduler query.

    Parameters
    ----------
    entries : list of tuple
        ``(spec, handle)`` or ``(spec, handle, job_id)``. All handles are
        reattached first, so a bad handle adopts nothing.

    Returns
    -------
    list of Job
        The adopted jobs with accurate states.
    """
    prepared = []
    for entry in entries:
      spec, handle = entry[0], entry[1]
      job_id = entry[2] if len(entry) > 2 else None
      self._prepare_spec(spec)
      prepared.append((spec, self.backend.reattach(handle), job_id))
    self._ensure_reaper_started()
    jobs = []
    with self._lock, self._owned():
      self._check_not_stopped()
      seen = set()
      for spec, handle, job_id in prepared:
        if job_id is not None and job_id in self._jobs:
          raise Sorry("Job id already in use: %s" % job_id)
        if job_id is not None and job_id in seen:
          raise Sorry("Job id given twice in one adopt_many call: %s"
                      % job_id)
        seen.add(job_id)
      for spec, handle, job_id in prepared:
        job = Job(job_id=job_id or _new_job_id(), spec=spec,
                  submitted_at=_utcnow(), handle=handle)
        self._jobs[job.job_id] = job
        jobs.append(job)
    if jobs:
      try:
        self.backend.refresh(jobs, force=True)
      except Exception:
        self._print_exc()
    now = _utcnow()
    mono = time.monotonic()
    with self._lock, self._owned():
      for job in jobs:
        result = self._safe_poll_locked(job)
        if result is not None:
          self._apply_poll_locked(job, result, now, mono)
    self._flush_logs()
    return jobs

  def shutdown(self):
    """Hand live jobs to ``backend.on_shutdown``; signal-safe, never raises."""
    snapshot = []
    got = self._lock.acquire(blocking=False)
    try:
      try:
        snapshot = list(self._jobs.values())
      except RuntimeError:
        try:
          snapshot = list(self._jobs.values())
        except RuntimeError:
          snapshot = []
    finally:
      if got:
        self._lock.release()
    live = [j for j in snapshot
            if j.handle is not None and j.state in ("running", "cancelling", "queued")]
    try:
      self.backend.on_shutdown(live)
    except Exception:
      self._print_exc()

  def reaper_iter(self):
    """One reaper pass (spec 5.4)."""
    with self._lock, self._owned():
      pollable = [j for j in self._jobs.values()
                  if j.handle is not None and not j._submitting
                  and j.state not in TERMINAL_STATES]
    if pollable:
      try:
        self.backend.refresh(pollable)
      except Exception:
        self._print_exc()
    now = _utcnow()
    mono = time.monotonic()
    capacity = self.backend.capacity()
    to_cancel, to_kill, to_submit = [], [], []
    with self._lock, self._owned():
      for job in pollable:
        if job.state in TERMINAL_STATES:
          continue
        result = self._safe_poll_locked(job)
        if result is not None:
          self._apply_poll_locked(job, result, now, mono)
      cap = self.max_job_seconds
      for job in self._jobs.values():
        if job.state == "running" and cap > 0 \
           and job._started_monotonic is not None \
           and mono - job._started_monotonic > cap:
          self._begin_cancel_locked(job, reason="wall-clock cap")
          self._queue_line_locked(
            job, "SORRY: job exceeded the wall-clock cap (%d s); cancelling."
            % cap)
          to_cancel.append(job)
        elif job.state == "cancelling" and not job._kill_sent \
             and job._cancel_sent_at is not None \
             and mono - job._cancel_sent_at > self.sigkill_after:
          job._kill_sent = True
          to_kill.append(job)
      while self._queue and self._slots_free_locked(capacity):
        job = self._jobs.get(self._queue.pop(0))
        if job is None or job.state != "queued" or job.handle is not None:
          continue
        job._submitting = True
        job.queue_position = 0
        to_submit.append(job)
      self._renumber_queue_locked()
    # Whatever the flush or the eviction raises, the terminate and kill
    # calls go out (a kill is sent once and never retried) and the jobs
    # marked _submitting are submitted, so no slot is leaked.
    try:
      self._flush_logs()
      try:
        with self._lock, self._owned():
          # After the flush, so no job is evicted before its lines are
          # written and its event set.
          self._evict_terminal_locked()
      except Exception:
        # A backend's handle_info can raise while a tombstone is built.
        self._print_exc()
    finally:
      try:
        for job in to_cancel:
          self._call_backend("terminate", job)
        for job in to_kill:
          self._call_backend("kill", job)
      finally:
        for job in to_submit:
          self._submit_job(job)

  # ---- submission ------------------------------------------------------------

  def _prepare_spec(self, spec):
    if not isinstance(spec, JobSpec):
      raise TypeError("submit expects a JobSpec (got %r)" % type(spec))
    self.backend.validate(spec)
    taken = set(CORE_STATUS_KEYS) | set(self.backend.handle_keys)
    clash = sorted(set(spec.metadata) & taken)
    if clash:
      raise ValueError(
        "metadata keys collide with status keys: %s" % ", ".join(clash))

  def _submit_job(self, job):
    """Call ``backend.submit`` with the lock released, then apply."""
    try:
      phase, handle = self.backend.submit(job)
    except SubmitError as exc:
      with self._lock, self._owned():
        job._submitting = False
        self._fail_locked(job, "failed to submit job: %s" % exc)
      self._flush_logs()
      return
    except Exception as exc:
      # Includes a malformed return value; failing the job frees its slot.
      self._print_exc()
      with self._lock, self._owned():
        job._submitting = False
        self._fail_locked(
          job, "backend.submit raised %s: %s" % (type(exc).__name__, exc))
      self._flush_logs()
      return
    except BaseException as exc:
      # Interrupt or SystemExit: fail the job to free its slot and
      # re-raise. A child the backend already started is not tracked.
      with self._lock, self._owned():
        job._submitting = False
        self._fail_locked(
          job, "submission interrupted by %s" % type(exc).__name__)
      self._flush_logs()
      raise
    need_terminate = False
    with self._lock, self._owned():
      job._submitting = False
      if job.state in TERMINAL_STATES:
        return
      job.handle = handle
      if phase == "running":
        job.state = "running"
        job.queue_position = 0
        job.started_at = _utcnow()
        job._started_monotonic = time.monotonic()
      elif phase == "pending":
        job.state = "queued"
        job.queue_position = 0
      else:
        # Fail the job and terminate whatever the backend started.
        self._fail_locked(job, "backend reported unknown phase %r" % (phase,))
        need_terminate = True
      if job._cancel_requested and not need_terminate:
        job._cancel_requested = False
        self._begin_cancel_locked(job)
        need_terminate = True
    self._flush_logs()
    if need_terminate:
      self._call_backend("terminate", job)

  def _slots_free_locked(self, capacity):
    """True when ``capacity`` has a free slot.

    Callers read ``capacity`` before taking the lock because only
    ``poll`` and ``handle_info`` may be called on the backend under it.
    """
    if capacity is None:
      return True
    used = sum(1 for j in self._jobs.values()
               if j._submitting or j.state in ("running", "cancelling"))
    return used < capacity

  # ---- transitions (lock held) ---------------------------------------------

  def _safe_poll_locked(self, job):
    try:
      result = self.backend.poll(job)
      if not isinstance(result, PollResult):
        raise TypeError("backend.poll returned %r, not a PollResult"
                        % (result,))
      return result
    except Exception:
      if job.job_id not in self._poll_errors_reported:
        self._poll_errors_reported.add(job.job_id)
        self._warn("poll failed for %s" % job.job_id)
        self._print_exc()
      return None

  def _apply_poll_locked(self, job, result, now, mono):
    if job.state in TERMINAL_STATES:
      return
    if result.phase == "exited":
      self._apply_exited_locked(job, result, now)
    elif result.phase == "running":
      if job.state == "queued":
        job.state = "running"
        job.queue_position = 0
        job.started_at = now
        job._started_monotonic = mono
      if job.state == "running":
        job.reason = result.reason
    elif result.phase == "pending":
      if job.state == "running":
        job.state = "queued"
        job.queue_position = 0
        job.started_at = None
        job._started_monotonic = None
      if job.state == "queued":
        job.reason = result.reason

  def _apply_exited_locked(self, job, result, now):
    job.exit_code = result.exit_code
    if job.state == "cancelling" and job._cancel_reason is not None:
      # Keep the manager's reason (wall-clock cap), append the backend's.
      job.reason = job._cancel_reason
      if result.reason is not None:
        job.reason = "%s (%s)" % (job._cancel_reason, result.reason)
    elif result.reason is not None or job.state != "cancelling":
      job.reason = result.reason
    job.finished_at = now
    if job.state == "cancelling":
      if result.cancelled is False:
        state = "finished" if result.exit_code == 0 else "failed"
      else:
        state = "cancelled"
    elif result.cancelled is True:
      state = "cancelled"
    else:
      state = "finished" if result.exit_code == 0 else "failed"
    job.state = state
    if state == "failed" and result.exit_code is None and result.reason:
      self._queue_line_locked(
        job, "SORRY: job ended without an exit code: %s" % result.reason)
    # Otherwise _flush_logs sets the event once the queued lines are
    # written, so a waiter finds them in the log.
    if job._pending_lines <= 0:
      job._terminal_event.set()

  def _begin_cancel_locked(self, job, reason=None):
    job.state = "cancelling"
    job.queue_position = 0
    job._cancel_sent_at = time.monotonic()
    if reason is not None:
      job.reason = reason
      job._cancel_reason = reason

  def _fail_locked(self, job, reason):
    if job.job_id in self._queue:
      self._queue.remove(job.job_id)
      self._renumber_queue_locked()
    job.state = "failed"
    job.queue_position = 0
    job.exit_code = None
    job.reason = reason
    job.finished_at = _utcnow()
    self._queue_line_locked(job, "SORRY: %s" % reason)

  def _queue_line_locked(self, job, line):
    """Queue a line for ``_flush_logs``; the terminal event waits for it."""
    job._pending_lines += 1
    self._log_queue.append((job, line))

  def _renumber_queue_locked(self):
    for pos, jid in enumerate(self._queue, start=1):
      job = self._jobs.get(jid)
      if job is not None:
        job.queue_position = pos

  def _evict_terminal_locked(self):
    # A terminal job whose event is unset still has a line being written;
    # keep it so no tombstone precedes its line.
    terminal = [j for j in self._jobs.values()
                if j.state in TERMINAL_STATES and j._terminal_event.is_set()]
    if len(terminal) <= self.max_retained_jobs:
      return
    terminal.sort(key=lambda j: j.finished_at or datetime.min)
    for job in terminal[:len(terminal) - self.max_retained_jobs]:
      if self.max_evicted_tombstones > 0:
        self._evicted.pop(job.job_id, None)
        self._evicted[job.job_id] = job.to_dict(self.backend)
        while len(self._evicted) > self.max_evicted_tombstones:
          self._evicted.pop(next(iter(self._evicted)))
      self._jobs.pop(job.job_id, None)
      self._poll_errors_reported.discard(job.job_id)

  # ---- helpers (lock released) ---------------------------------------------

  @staticmethod
  def _print_exc():
    """``traceback.print_exc`` to stderr only; never raises."""
    stream = sys.stderr                # read once: it may be reassigned
    if stream is None:
      return
    try:
      traceback.print_exc(file=stream)
    except Exception:
      pass

  @staticmethod
  def _warn(text):
    """Print a note to stderr only; never raises."""
    stream = sys.stderr
    if stream is None:
      return
    try:
      print("libtbx.jobs: %s" % text, file=stream)
    except Exception:
      pass

  def _call_backend(self, name, job):
    try:
      getattr(self.backend, name)(job)
    except Exception:
      self._warn("backend.%s failed for %s" % (name, job.job_id))
      self._print_exc()

  def _flush_logs(self):
    """Write queued log lines with the lock released.

    Lines are queued under the lock and written here so a slow or hung
    file system never stalls the lock (spec 5.2).
    """
    remaining = []
    try:
      with self._lock, self._owned():
        remaining, self._log_queue = self._log_queue, []
      while remaining:
        job, line = remaining[0]
        try:
          with open(job.spec.log_path, "a", encoding="utf-8") as fh:
            fh.write("\n%s\n" % line)
        except OSError:
          pass
        except Exception:
          self._print_exc()
        # Written or not, the line is done; wake the job's waiters now,
        # not after the rest of the batch.
        with self._lock, self._owned():
          job._pending_lines -= 1
          remaining.pop(0)
          if job.state in TERMINAL_STATES and job._pending_lines <= 0:
            job._terminal_event.set()
    finally:
      if remaining:
        # Interrupted: requeue the unwritten lines to keep the queue
        # consistent.
        with self._lock, self._owned():
          self._log_queue = remaining + self._log_queue

  # ---- reaper thread ---------------------------------------------------------

  def _check_not_stopped(self):
    if self._stop_event.is_set():
      raise RuntimeError(
        "JobManager.stop() was called; create a new manager for new jobs")

  def _ensure_reaper_started(self):
    self._check_not_stopped()
    with self._lock, self._owned():
      if self._reaper_started:
        return
      self._reaper_started = True
    thread = threading.Thread(
      target=self._reaper_loop, name="libtbx-jobs-reaper", daemon=True)
    try:
      thread.start()
    except Exception:
      with self._lock, self._owned():
        self._reaper_started = False
      raise
    self._reaper_thread = thread

  def _reaper_loop(self):
    while not self._stop_event.wait(self.reaper_interval):
      try:
        self.reaper_iter()
      except Exception:
        self._print_exc()

  def stop(self, timeout=None):
    """Ask the reaper thread to exit; join it when ``timeout`` is given.

    Final: later ``submit`` or ``adopt`` calls raise ``RuntimeError``,
    since nothing would poll or drain the FIFO. A pass under way
    finishes; after that job states freeze (no polling, no FIFO
    submissions), though ``cancel`` still reaches the backend. Meant
    for a host that is exiting.
    """
    with self._lock, self._owned():
      self._stop_event.set()
    thread = self._reaper_thread
    if thread is not None and timeout is not None and thread.is_alive():
      thread.join(timeout)
