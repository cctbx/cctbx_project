"""SLURM backend (spec 9)."""

import os
import re
import shutil
import tempfile

from libtbx.jobs.backends.base import SubmitError
from libtbx.jobs.backends.batch import BatchBackend, SchedulerTimeout


_PENDING = frozenset(("PENDING", "CONFIGURING", "REQUEUED", "REQUEUE_HOLD",
                      "REQUEUE_FED", "RESV_DEL_HOLD", "SPECIAL_EXIT"))
_RUNNING = frozenset(("RUNNING", "COMPLETING", "SUSPENDED", "STAGE_OUT",
                      "SIGNALING", "STOPPED", "RESIZING"))
_EXITED = frozenset(("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
                     "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL",
                     "DEADLINE", "REVOKED", "LAUNCH_FAILED"))
_REQUEUED = frozenset(("REQUEUED", "REQUEUE_HOLD", "REQUEUE_FED",
                       "SPECIAL_EXIT"))
# Scheduler-log markers per raw state; a line explains only its own state.
_HINTS = {
  "TIMEOUT": ("DUE TO TIME LIMIT",),
  "OUT_OF_MEMORY": ("oom-kill", "Out Of Memory"),
  "PREEMPTED": ("DUE TO PREEMPTION",),
  "CANCELLED": ("CANCELLED AT",),
}
_ID_RE = re.compile(r"\s*(\d+)(?:;([A-Za-z0-9_.-]+))?\s*$")


def _word(raw_state):
  """First word of a raw state, upper-cased (``CANCELLED by 1234``)."""
  words = (raw_state or "").split()
  return words[0].upper() if words else ""


class SlurmBackend(BatchBackend):
  """Submit through ``sbatch``, poll with ``squeue`` and ``sacct``.

  Accepts the ``BatchBackend`` keyword arguments. The generated
  ``--gpus`` needs SLURM 19.05+; older sites use ``--gres=gpu:{{gpus}}``
  in their template. With ``--clusters`` in the template, the cluster
  ``sbatch --parsable`` reports is kept in the handle and passed to every
  later ``squeue``, ``sacct`` and ``scancel``.
  """

  name = "slurm"
  directive_prefix = "#SBATCH"
  script_extension = "sbatch"

  @classmethod
  def default_template_path(cls):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "slurm_default.sbatch")

  def submit_command(self, job, script_path):
    job_dir = os.path.dirname(script_path)
    scheduler_log = os.path.join(job_dir, "scheduler.log")
    cmd = ["sbatch", "--parsable",
           "--chdir=%s" % job.spec.cwd,
           "--output=%s" % scheduler_log,
           "--error=%s" % scheduler_log,
           # Keep earlier attempts' output across a requeue.
           "--open-mode=append",
           "--job-name=%s" % self.job_name(job.spec)]
    if self.requeue is not None:
      cmd.append("--requeue" if self.requeue else "--no-requeue")
    return cmd + [script_path]

  def parse_submit_output(self, stdout, stderr):
    # Only the first non-empty line: a site wrapper's trailing banner must
    # not orphan an accepted job.
    lines = [l for l in (stdout or "").splitlines() if l.strip()]
    match = _ID_RE.match(lines[0]) if lines else None
    if not match:
      raise SubmitError("unexpected sbatch output: %r %r" % (stdout, stderr))
    return match.group(1), match.group(2)

  @staticmethod
  def _cluster_args(cluster):
    return ["--clusters=%s" % cluster] if cluster else []

  def state_phase(self, raw_state):
    word = _word(raw_state)
    if word in _PENDING:
      return "pending"
    if word in _EXITED:
      return "exited"
    return "running"

  def is_cancelled_state(self, raw_state):
    return _word(raw_state) == "CANCELLED"

  def is_requeue_state(self, raw_state):
    return _word(raw_state) in _REQUEUED

  def is_completed_state(self, raw_state):
    return _word(raw_state) == "COMPLETED"

  def is_failed_state(self, raw_state):
    return _word(raw_state) == "FAILED"

  @staticmethod
  def _submitting_user():
    """The user name of the process's uid (not LOGNAME).

    ``None`` for a uid without a passwd entry (a container's ``-u``): the
    query then covers all users and the client-side id filter suffices.
    """
    try:
      uid = os.getuid()
    except AttributeError:
      import getpass
      return getpass.getuser()
    try:
      import pwd
      return pwd.getpwuid(uid).pw_name
    except (ImportError, KeyError):
      return None

  def _squeue_command(self, cluster=None):
    cmd = ["squeue", "--noheader", "--states=all"]
    user = self._submitting_user()
    if user:
      cmd.append("--user=%s" % user)
    cmd.append("--format=%i|%T|%r")
    return cmd + self._cluster_args(cluster)

  def attempt_variable(self):
    return "SLURM_RESTART_COUNT"

  def query_live(self, ids, cluster=None):
    cmd = self._squeue_command(cluster)
    try:
      proc = self.run_command(cmd)
    except (SchedulerTimeout, OSError) as exc:
      raise RuntimeError("squeue: %s" % exc)
    if proc.returncode != 0:
      raise RuntimeError("squeue exited %d: %s"
                         % (proc.returncode, proc.stderr.strip()))
    wanted = set(ids)
    out = {}
    for line in proc.stdout.splitlines():
      parts = line.split("|", 2)
      if len(parts) < 2:
        continue
      job_id = parts[0].strip()
      if job_id not in wanted:
        continue
      raw = parts[1].strip()
      reason = parts[2].strip() if len(parts) > 2 else ""
      if reason in ("", "None", "(null)"):
        reason = None
      phase = self.state_phase(raw)
      out[job_id] = (phase, raw, reason if phase == "pending" else None)
    return out

  def query_terminal(self, ids, cluster=None):
    cmd = ["sacct", "--noheader", "--parsable2", "--jobs=%s" % ",".join(ids),
           "--format=JobID,State,ExitCode"] + self._cluster_args(cluster)
    try:
      proc = self.run_command(cmd)
    except SchedulerTimeout:
      raise
    except OSError as exc:
      raise RuntimeError("sacct: %s" % exc)
    if proc.returncode != 0:
      raise RuntimeError("sacct exited %d: %s"
                         % (proc.returncode, proc.stderr.strip()))
    wanted = set(ids)
    out = {}
    for line in proc.stdout.splitlines():
      parts = line.split("|")
      if len(parts) < 3:
        continue
      job_id, state, exitcode = (p.strip() for p in parts[:3])
      if job_id not in wanted:
        continue
      code_text, _, signal_text = exitcode.partition(":")
      try:
        code = int(code_text)
        signum = int(signal_text or "0")
      except ValueError:
        code, signum = None, 0
      out[job_id] = (code, _word(state), "signal %d" % signum if signum else None)
    return out

  def cancel_command(self, job):
    handle = job.handle
    return ["scancel"] + self._cluster_args(handle.cluster) \
      + [handle.scheduler_job_id]

  def resource_directives(self, resources):
    lines = []
    if "cpus" in resources:
      lines.append("#SBATCH --cpus-per-task=%d" % resources["cpus"])
    if "mem_mb" in resources:
      lines.append("#SBATCH --mem=%dM" % resources["mem_mb"])
    if "time_minutes" in resources:
      lines.append("#SBATCH --time=%s" % self.format_time(resources["time_minutes"]))
    if "gpus" in resources:
      lines.append("#SBATCH --gpus=%d" % resources["gpus"])
    return lines

  def failure_hint(self, text, raw_state):
    markers = _HINTS.get(_word(raw_state), ())
    for line in reversed(text.splitlines()):
      # A requeue notice explains an abandoned attempt, never the end.
      if "DUE TO JOB REQUEUE" in line:
        continue
      if any(marker in line for marker in markers):
        return line.strip()
    return None

  def check(self, spec):
    """``validate``, then ``sbatch --test-only`` with the real submit
    command line, so its options (such as the requeue flag) are checked."""
    from libtbx.jobs.manager import Job
    self.validate(spec)
    tmp = tempfile.mkdtemp(prefix="libtbx_jobs_check_")
    try:
      path = os.path.join(tmp, "job.sbatch")
      with open(path, "w", encoding="utf-8") as fh:
        fh.write(self.render(spec))
      cmd = self.submit_command(Job(job_id="j_check", spec=spec), path)
      cmd.insert(len(cmd) - 1, "--test-only")   # before the script
      try:
        proc = self.run_command(cmd)
      except (SchedulerTimeout, OSError) as exc:
        raise SubmitError("sbatch: %s" % exc)
    finally:
      shutil.rmtree(tmp, ignore_errors=True)
    if proc.returncode != 0:
      raise SubmitError(proc.stderr.strip() or proc.stdout.strip()
                        or "sbatch --test-only exited %d" % proc.returncode)
