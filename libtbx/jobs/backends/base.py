"""The contract between ``JobManager`` and an execution backend."""

from dataclasses import dataclass

from libtbx.utils import Sorry


_PHASES = ("pending", "running", "exited")


@dataclass(frozen=True)
class PollResult:
  """What a backend knows about one job at one moment.

  Parameters
  ----------
  phase : str
      ``"pending"`` (accepted, not yet running), ``"running"`` or
      ``"exited"``.
  exit_code : int, optional
      Meaningful only when ``phase == "exited"``; ``None`` means unknown.
  reason : str, optional
      Scheduler state, scheduler reason, or a short explanation.
  cancelled : bool, optional
      For ``exited`` only. ``True`` when the backend knows the job was
      cancelled, ``False`` when it knows the job ended on its own,
      ``None`` when it cannot tell.
  """

  phase: str
  exit_code: int = None
  reason: str = None
  cancelled: bool = None

  def __post_init__(self):
    if self.phase not in _PHASES:
      raise ValueError(
        "PollResult.phase must be one of %s (got %r)"
        % (", ".join(_PHASES), self.phase))


class SubmitError(Exception):
  """The backend could not accept a job.

  ``JobManager`` turns this into a failed job whose ``reason`` is the
  message; it never propagates to the caller of ``submit``.
  """


class Backend:
  """Abstract execution backend.

  ``poll`` and ``handle_info`` are called with the manager lock held and
  must be quick and non-blocking. Every other method is called with the
  lock released and may run subprocesses. ``on_shutdown`` may run inside a
  signal handler and must not block on any lock.
  """

  name = "abstract"
  handle_keys = ()

  def validate(self, spec):
    """Raise ``ValueError`` or a subclass when ``spec`` cannot run here."""

  def check(self, spec):
    """Like ``validate`` but allowed to contact the scheduler."""
    self.validate(spec)

  def submit(self, job):
    """Start ``job``.

    Returns
    -------
    tuple
        ``(phase, handle)`` where ``phase`` is ``"running"`` or
        ``"pending"`` and ``handle`` is whatever ``poll`` needs later.
        Must not mutate ``job``. Raises ``SubmitError`` on refusal.
    """
    raise NotImplementedError

  def refresh(self, jobs, force=False):
    """Fetch fresh state for ``jobs`` in one batch (optional)."""

  def poll(self, job):
    """Return a ``PollResult`` for ``job`` without blocking."""
    raise NotImplementedError

  def terminate(self, job):
    """Ask the job to stop politely."""
    raise NotImplementedError

  def kill(self, job):
    """Stop the job forcefully; default no-op."""

  def capacity(self):
    """Number of jobs that may run at once, or ``None`` when unbounded."""
    return None

  def on_shutdown(self, jobs):
    """The owning process is exiting with these jobs still live."""

  def handle_info(self, job):
    """Extra status-dict keys for ``job`` (must be in ``handle_keys``)."""
    return {}

  def reattach(self, handle):
    """Turn a persisted handle back into one ``poll`` accepts."""
    raise Sorry(
      "The %s backend does not support reattaching jobs" % self.name)
