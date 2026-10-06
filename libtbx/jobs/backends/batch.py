"""Shared machinery for batch-scheduler backends (spec 8).

A ``BatchBackend`` renders a site template into a batch script, writes a
wrapper (working directory, environment, live log, exit-code sentinel),
submits through adapter hooks a subclass fills in, and polls in batches.

Template placeholders
---------------------
A site template is shell text with ``{{name}}`` placeholders; everything
else, including ``$SLURM_JOB_ID`` and ``$(hostname)``, passes through.
``{{command}}`` is required. The others:

``{{resource_directives}}``
    One directive line per requested resource (``cpus``, ``mem_mb``,
    ``time_minutes``, ``gpus``). Put it below any literal directive for
    the same resource so the job's request wins.
``{{environment}}``
    ``export`` lines for the job's environment overrides.
``{{job_name}}``, ``{{job_id}}``, ``{{work_dir}}``, ``{{log_path}}``,
``{{scheduler_log}}``
    Per-job values.
``{{cpus}}``, ``{{mem_mb}}``, ``{{time_minutes}}``, ``{{gpus}}``,
``{{time}}``
    The job's resource values; ``{{time}}`` is ``time_minutes`` as
    ``H:MM:SS``. Referencing one anywhere in the template takes that
    resource over: the backend emits no directive for it and every job
    must request it, or rendering fails. Use them only in a directive
    you write yourself, such as ``--gres=gpu:{{gpus}}``.
Any other name
    A site default or a free-form job resource value (spec 4.1).

End the script with ``{{command}}`` if the scheduler's accounting should
report the program's exit code; the wrapper's sentinel is authoritative
either way. Comment lines must not hold placeholders.
"""

import math
import os
import re
import shlex
import stat
import string
import subprocess
import sys
import threading
import time
from dataclasses import dataclass

from libtbx.jobs.backends.base import Backend, PollResult, SubmitError
from libtbx.jobs.manager import RESOURCE_KEYS, validate_placeholder_value
from libtbx.utils import Sorry


def _warn(text):
  """Note on stderr; never raises, never falls back to stdout (which may
  be a protocol channel)."""
  stream = sys.stderr                  # read once: it may be reassigned
  if stream is None:
    return
  try:
    print("libtbx.jobs: %s" % text, file=stream)
  except Exception:
    pass


class TemplateError(ValueError):
  """A site template cannot be rendered for this job."""


class SchedulerTimeout(RuntimeError):
  """A scheduler command did not finish within ``command_timeout``."""


class BraceTemplate(string.Template):
  """``string.Template`` that recognises only ``{{name}}``.

  Shell syntax (``$VAR``, ``${VAR}``, ``$(cmd)``, ``$?``) passes through.
  Whitespace inside the braces is allowed; there is no escape.
  """

  delimiter = "{{"
  pattern = r"""
    \{\{\s*(?P<braced>[_a-zA-Z][_a-zA-Z0-9]*)\s*\}\}
    | (?P<escaped>(?!))
    | (?P<named>(?!))
    | (?P<invalid>(?!))
  """


def referenced_placeholders(text):
  """Return the set of placeholder names ``text`` references."""
  return {m.group("braced") for m in BraceTemplate.pattern.finditer(text)}


RESERVED_PLACEHOLDERS = frozenset((
  "command", "resource_directives", "environment", "job_name", "job_id",
  "work_dir", "log_path", "scheduler_log",
  "cpus", "mem_mb", "time_minutes", "gpus", "time"))
_TAKEOVER = {"cpus": "cpus", "mem_mb": "mem_mb", "gpus": "gpus",
             "time": "time_minutes", "time_minutes": "time_minutes"}
PATH_RE = re.compile(r"^[A-Za-z0-9_./+@:,=-]+$")
SCHEDULER_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]*")
CLUSTER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
JOB_DIR_NAME = ".libtbx_jobs"
ATTEMPT_MARKER = "libtbx.jobs: attempt"
_NAME_BAD = re.compile(r"[^A-Za-z0-9_.-]")


@dataclass(frozen=True)
class BatchHandle:
  """What identifies a submitted batch job.

  Parameters
  ----------
  scheduler_job_id : str
      The scheduler's id.
  job_dir : str
      ``<cwd>/.libtbx_jobs/<job_id>``, holding the sentinel and the
      scheduler log.
  cluster : str, optional
      The cluster the scheduler routed the job to, if named; later
      commands for the job are addressed to it.
  """

  scheduler_job_id: str
  job_dir: str
  cluster: str = None


class BatchBackend(Backend):
  """Base class for scheduler backends; subclasses fill in the hooks.

  Parameters
  ----------
  template : str, optional
      Path of the site template. ``None`` uses the shipped default.
  template_text : str, optional
      Inline template text; exclusive with ``template``.
  defaults : dict, optional
      Site placeholder values (spec 4.1 value rules); reserved names and
      resource names are rejected.
  poll_interval : float, optional
      Seconds between scheduler queries (minimum 1).
  sentinel_grace : float, optional
      Seconds to wait for a late exit-code sentinel before trusting the
      scheduler's raw state or accounting's exit code.
  command_timeout : float, optional
      Seconds before a scheduler command is killed and reported as a
      failure (default 60, at most a day).
  requeue : bool, optional
      Whether the scheduler may rerun a job from the start after a node
      failure or preemption. ``None`` (default) keeps the scheduler's
      policy; ``False`` forbids it (pass it if the host cannot rerun a
      job in place); ``True`` asks for it. A rerun repeats the command in
      the same directory against the same log.
  """

  name = "batch"
  handle_keys = ("scheduler", "scheduler_job_id", "job_dir",
                 "scheduler_cluster")
  directive_prefix = "#BATCH"
  script_extension = "batch"
  path_re = PATH_RE

  def __init__(self, template=None, template_text=None, defaults=None,
               poll_interval=10.0, sentinel_grace=60.0,
               command_timeout=60.0, requeue=None):
    if template is not None and template_text is not None:
      raise ValueError("give template or template_text, not both")
    self.template_path = None
    if template is not None:
      self.template_path = str(template)
    elif template_text is None:
      self.template_path = self.default_template_path()
    if self.template_path is not None:
      with open(self.template_path, encoding="utf-8") as fh:
        template_text = fh.read()
    self.template_text = template_text
    self._referenced = referenced_placeholders(template_text)
    leftover = BraceTemplate.pattern.sub("", template_text)
    if "{{" in leftover:
      raise TemplateError(
        "template contains '{{' that is not a placeholder: a placeholder is "
        "mistyped or the template uses '{{' literally")
    for line in template_text.splitlines():
      stripped = line.lstrip()
      if stripped.startswith("#") and not stripped.startswith("#!") \
         and not stripped.startswith(self.directive_prefix) \
         and BraceTemplate.pattern.search(stripped):
        raise TemplateError(
          "placeholders are not allowed in comment lines, because they are "
          "substituted there too; describe them in words (line %r)" % line)
    self.defaults = {}
    for key, value in (defaults or {}).items():
      if not isinstance(key, str):
        raise ValueError("defaults keys must be str (got %r)" % (key,))
      if key in RESERVED_PLACEHOLDERS:
        raise ValueError(
          "site default %r is a reserved placeholder name; write a literal "
          "directive in the template instead" % key)
      validate_placeholder_value(key, value)
      self.defaults[key] = value
    for label, value in (("poll_interval", poll_interval),
                         ("sentinel_grace", sentinel_grace),
                         ("command_timeout", command_timeout)):
      if not math.isfinite(float(value)):
        raise ValueError("%s must be a finite number (got %r)" % (label, value))
    self.poll_interval = max(1.0, float(poll_interval))
    self.sentinel_grace = max(0.0, float(sentinel_grace))
    # A huge timeout overflows subprocess's selector; cap it at a day.
    self.command_timeout = min(max(1.0, float(command_timeout)), 86400.0)
    if requeue is not None and not isinstance(requeue, bool):
      # "false" or 0 from a configuration file must not mean True.
      raise ValueError("requeue must be True, False or None (got %r)"
                       % (requeue,))
    self.requeue = requeue
    self.accounting_available = True
    self._error_streaks = set()
    self._cache_lock = threading.Lock()
    self._cache = {}
    self._stamps = {}
    self._submitted = {}
    self._first_not_live = {}
    self._reattached = set()
    self._last_refresh = None

  # ---- adapter hooks (spec 8.5) --------------------------------------------

  @classmethod
  def default_template_path(cls):
    """Path of the shipped default template."""
    raise NotImplementedError

  def submit_command(self, job, script_path):
    """argv that submits ``script_path``, whose directory holds
    ``scheduler.log``."""
    raise NotImplementedError

  def parse_submit_output(self, stdout, stderr):
    """The scheduler job id, or ``(job_id, cluster)`` for a routed job;
    raise ``SubmitError`` if there is none."""
    raise NotImplementedError

  def query_live(self, ids, cluster=None):
    """Map each reported id to ``(phase, raw_state, reason)``; terminal
    states have phase ``"exited"``. ``cluster=None`` is the default."""
    raise NotImplementedError

  def query_terminal(self, ids, cluster=None):
    """Map scheduler id to ``(exit_code, raw_state, reason)`` from
    accounting; raise when accounting is unavailable."""
    return {}

  def state_phase(self, raw_state):
    """``"pending"``, ``"running"`` or ``"exited"`` for a raw state."""
    raise NotImplementedError

  def is_cancelled_state(self, raw_state):
    raise NotImplementedError

  def is_completed_state(self, raw_state):
    raise NotImplementedError

  def is_failed_state(self, raw_state):
    raise NotImplementedError

  def cancel_command(self, job):
    """argv that cancels ``job``."""
    raise NotImplementedError

  def resource_directives(self, resources):
    """Directive lines for the given subset of ``RESOURCE_KEYS``."""
    raise NotImplementedError

  def format_time(self, minutes):
    """Render a time limit; default ``H:MM:SS``."""
    return "%d:%02d:00" % (minutes // 60, minutes % 60)

  def failure_hint(self, text, raw_state):
    """A scheduler-log line in ``text`` explaining ``raw_state``, or
    ``None``; only a line matching that state counts."""
    return None

  def is_requeue_state(self, raw_state):
    """Whether a pending raw state says the scheduler requeued the job."""
    return False

  def attempt_variable(self):
    """The scheduler's variable counting a job's restarts, or ``None``."""
    return None

  def job_name(self, spec):
    """Scheduler-safe job name: ``[A-Za-z0-9_.-]``, at most 64 chars."""
    return _NAME_BAD.sub("_", spec.name)[:64] or "job"

  def run_command(self, cmd, **kwargs):
    """Run a scheduler command with ``command_timeout``.

    Returns the ``CompletedProcess`` with text output. Raises
    ``SchedulerTimeout`` on overrun and ``OSError`` if it cannot start.
    """
    kwargs.setdefault("capture_output", True)
    kwargs.setdefault("text", True)
    kwargs.setdefault("timeout", self.command_timeout)
    try:
      return subprocess.run(cmd, **kwargs)
    except subprocess.TimeoutExpired:
      raise SchedulerTimeout("%s did not finish within %g s"
                             % (cmd[0], kwargs["timeout"]))

  # ---- rendering -------------------------------------------------------------

  def job_dir(self, spec, job_id):
    """``<cwd>/.libtbx_jobs/<job_id>``."""
    return os.path.join(spec.cwd, JOB_DIR_NAME, job_id)

  def _context(self, spec, job_id, job_dir):
    resources = spec.resources
    mapping = dict(self.defaults)
    mapping.update(resources)
    taken = {_TAKEOVER[n] for n in self._referenced if n in _TAKEOVER}
    for_directives = {k: resources[k] for k in RESOURCE_KEYS
                      if k in resources and k not in taken}
    lines = self.resource_directives(for_directives) if for_directives else []
    env_lines = ["export %s=%s" % (k, shlex.quote(v))
                 for k, v in spec.env.items()]
    mapping.update({
      "command": "bash %s" % shlex.quote(os.path.join(job_dir, "run.sh")),
      "resource_directives": "\n".join(lines),
      "environment": "\n".join(env_lines),
      "job_name": self.job_name(spec),
      "job_id": job_id,
      "work_dir": spec.cwd,
      "log_path": spec.log_path,
      "scheduler_log": os.path.join(job_dir, "scheduler.log"),
    })
    for key in RESOURCE_KEYS:
      if key in resources:
        mapping[key] = resources[key]
    if "time_minutes" in resources:
      mapping["time"] = self.format_time(resources["time_minutes"])
    return mapping

  def _render(self, spec, job_id, job_dir):
    mapping = self._context(spec, job_id, job_dir)
    try:
      text = BraceTemplate(self.template_text).substitute(mapping)
    except KeyError as exc:
      name = exc.args[0]
      if name in RESERVED_PLACEHOLDERS:
        raise TemplateError(
          "template references {{%s}} but the job did not request that "
          "resource" % name)
      raise TemplateError("template references unknown placeholder {{%s}}"
                          % name)
    except ValueError as exc:
      raise TemplateError("template syntax error: %s" % exc)
    return text

  def render(self, spec_or_job):
    """Return the rendered script text without submitting anything."""
    if hasattr(spec_or_job, "spec"):
      spec, job_id = spec_or_job.spec, spec_or_job.job_id
    else:
      spec, job_id = spec_or_job, "j_render"
    return self._render(spec, job_id, self.job_dir(spec, job_id))

  def validate(self, spec):
    for key in spec.resources:
      if key in RESERVED_PLACEHOLDERS and key not in RESOURCE_KEYS:
        raise ValueError(
          "resources key %r is a reserved placeholder name" % key)
    reserved_env = {"LIBTBX_JOBS_ATTEMPT", "libtbx_jobs_attempt",
                    self.attempt_variable()}
    for key in spec.env:
      if key in reserved_env:
        raise ValueError(
          "env key %r is reserved for the wrapper's attempt count" % key)
    for label in ("cwd", "log_path"):
      value = getattr(spec, label)
      if not PATH_RE.fullmatch(value):
        raise ValueError(
          "%s must contain only A-Z a-z 0-9 _ . / + @ : , = - for a batch "
          "backend, because it appears unquoted in directives (got %r)"
          % (label, value))
    if "command" not in self._referenced:
      raise TemplateError("template does not reference {{command}}")
    self._render(spec, "j_validate", self.job_dir(spec, "j_validate"))

  # ---- files -----------------------------------------------------------------

  def _write_wrapper(self, spec, job_id, job_dir):
    """Write ``run.sh`` (spec 8.3) and return its path."""
    q = shlex.quote
    sentinel = os.path.join(job_dir, "exit_code")
    lines = [
      "#!/bin/bash",
      "# Generated by libtbx.jobs for job %s. Do not edit." % job_id,
      "cd %s || { echo 'libtbx.jobs: cannot chdir to' %s >&2; "
      "echo 127 > %s; exit 127; }" % (q(spec.cwd), q(spec.cwd), q(sentinel)),
      # A requeued job must not inherit the previous attempt's exit code.
      "rm -f %s" % q(sentinel),
    ]
    variable = self.attempt_variable()
    if variable:
      # Read before spec.env is applied, so a host cannot steer the count.
      lines.append("libtbx_jobs_attempt=\"${%s:-0}\"" % variable)
    for key, value in spec.env.items():
      lines.append("export %s=%s" % (key, q(value)))
    if variable:
      # LIBTBX_JOBS_ATTEMPT counts earlier attempts (0 on the first); set
      # after spec.env so a host cannot override it. A rerun marks both
      # logs on a fresh line, as the killed attempt may have left a partial.
      lines += [
        "export LIBTBX_JOBS_ATTEMPT=\"$libtbx_jobs_attempt\"",
        "if [ \"$LIBTBX_JOBS_ATTEMPT\" != 0 ]; then "
        "marker=\"%s $((LIBTBX_JOBS_ATTEMPT + 1)) of this job after a "
        "requeue\"; printf '\\n%%s\\n' \"$marker\" >> %s; "
        "printf '\\n%%s\\n' \"$marker\"; fi"
        % (ATTEMPT_MARKER, q(spec.log_path)),
      ]
    lines.append("%s >> %s 2>&1" % (" ".join(q(a) for a in spec.argv),
                                   q(spec.log_path)))
    lines += ["rc=$?", "echo \"$rc\" > %s" % q(sentinel), "exit \"$rc\"", ""]
    path = os.path.join(job_dir, "run.sh")
    with open(path, "w", encoding="utf-8") as fh:
      fh.write("\n".join(lines))
    os.chmod(path, stat.S_IRWXU)
    return path

  def _remove_sentinel(self, handle):
    try:
      os.remove(os.path.join(handle.job_dir, "exit_code"))
    except OSError:
      pass

  def _read_sentinel(self, handle):
    """The exit code the wrapper wrote, or ``None`` if not there yet."""
    try:
      with open(os.path.join(handle.job_dir, "exit_code"),
                encoding="utf-8") as fh:
        return int(fh.read().strip())
    except (OSError, ValueError):
      return None

  # ---- Backend contract ------------------------------------------------------

  def submit(self, job):
    spec = job.spec
    job_dir = self.job_dir(spec, job.job_id)
    os.makedirs(job_dir, exist_ok=True)
    self._write_wrapper(spec, job.job_id, job_dir)
    script_path = os.path.join(job_dir, "job.%s" % self.script_extension)
    with open(script_path, "w", encoding="utf-8") as fh:
      fh.write(self._render(spec, job.job_id, job_dir))
    cmd = self.submit_command(job, script_path)
    try:
      proc = self.run_command(cmd, cwd=spec.cwd)
    except SchedulerTimeout as exc:
      raise SubmitError("%s; the scheduler may still have accepted the job"
                        % exc)
    except OSError as exc:
      raise SubmitError("%s: %s" % (cmd[0], exc))
    if proc.returncode != 0:
      raise SubmitError(
        (proc.stderr or proc.stdout).strip()
        or "%s exited with status %d" % (cmd[0], proc.returncode))
    parsed = self.parse_submit_output(proc.stdout, proc.stderr)
    if isinstance(parsed, tuple):
      scheduler_job_id, cluster = parsed
    else:
      scheduler_job_id, cluster = parsed, None
    handle = self._make_handle(scheduler_job_id, job_dir, cluster)
    now = time.monotonic()
    with self._cache_lock:
      self._cache[handle] = PollResult("pending")
      self._stamps[handle] = now
      self._submitted[handle] = now
    return "pending", handle

  def poll(self, job):
    with self._cache_lock:
      return self._cache.get(job.handle, PollResult("pending"))

  def terminate(self, job):
    if not isinstance(job.handle, BatchHandle):
      return
    cmd = self.cancel_command(job)
    try:
      self.run_command(cmd)
    except (SchedulerTimeout, OSError) as exc:
      _warn("%s failed: %s" % (cmd[0], exc))

  def kill(self, job):
    """No-op: schedulers escalate from SIGTERM to SIGKILL themselves."""

  def on_shutdown(self, jobs):
    for job in jobs:
      handle = job.handle
      if isinstance(handle, BatchHandle):
        where = handle.scheduler_job_id
        if handle.cluster:
          where = "%s on cluster %s" % (where, handle.cluster)
        _warn("%s job %s is still live (job dir %s)"
              % (self.name, where, handle.job_dir))

  def handle_info(self, job):
    handle = job.handle
    if not isinstance(handle, BatchHandle):
      return {"scheduler": self.name, "scheduler_job_id": None,
              "job_dir": None, "scheduler_cluster": None}
    return {"scheduler": self.name,
            "scheduler_job_id": handle.scheduler_job_id,
            "job_dir": handle.job_dir,
            "scheduler_cluster": handle.cluster}

  @staticmethod
  def _make_handle(scheduler_job_id, job_dir, cluster=None):
    """Build a ``BatchHandle`` from checked parts (``ValueError`` if bad)."""
    scheduler_job_id = str(scheduler_job_id)
    if not SCHEDULER_ID_RE.fullmatch(scheduler_job_id):
      raise ValueError("scheduler job id %r must match %s"
                       % (scheduler_job_id, SCHEDULER_ID_RE.pattern))
    if cluster is not None and str(cluster) == "":
      cluster = None                 # a persisted empty field
    if cluster is not None:
      cluster = str(cluster)
      if not CLUSTER_RE.fullmatch(cluster):
        raise ValueError("cluster name %r must match %s"
                         % (cluster, CLUSTER_RE.pattern))
    return BatchHandle(scheduler_job_id, str(job_dir), cluster)

  def reattach(self, handle):
    if isinstance(handle, BatchHandle):
      parts = (handle.scheduler_job_id, handle.job_dir, handle.cluster)
    elif isinstance(handle, dict):
      try:
        parts = (handle["scheduler_job_id"], handle["job_dir"],
                 handle.get("scheduler_cluster"))
      except KeyError as exc:
        raise Sorry("cannot reattach: handle lacks %s" % exc)
    else:
      raise Sorry("cannot reattach: expected a BatchHandle or a dict with "
                  "scheduler_job_id and job_dir (got %r)" % (handle,))
    if parts[0] is None or parts[0] == "":
      raise Sorry("cannot reattach: empty scheduler job id")
    try:
      result = self._make_handle(*parts)
    except ValueError as exc:
      raise Sorry("cannot reattach: %s" % exc)
    if not os.path.isabs(result.job_dir):
      raise Sorry("cannot reattach: job_dir must be absolute (got %r)"
                  % result.job_dir)
    with self._cache_lock:
      self._submitted.setdefault(result, time.monotonic())
      self._reattached.add(result)
    return result

  # ---- state refresh (spec 8.4) ----------------------------------------------

  def refresh(self, jobs, force=False):
    handles = [j.handle for j in jobs if isinstance(j.handle, BatchHandle)]
    if not handles:
      return
    started = time.monotonic()
    with self._cache_lock:
      if not force and self._last_refresh is not None \
         and started - self._last_refresh < self.poll_interval:
        return
      self._last_refresh = started
      submitted = dict(self._submitted)
      first_not_live = dict(self._first_not_live)
      phases = {h: r.phase for h, r in self._cache.items()}
      reattached = set(self._reattached)
    consumed = set()
    # One query per cluster; a failed one leaves only its handles stale.
    live = {}
    for cluster, ids in self._by_cluster(handles):
      try:
        live[cluster] = self.query_live(ids, cluster=cluster)
      except Exception as exc:
        self._note_error("query_live", exc, cluster)
        continue
      self._clear_error("query_live", cluster)
      # A sentinel seen while the job is pending is from an abandoned
      # attempt. Remove it on a requeue state, or if the handle was seen
      # running or exited, or was just reattached; not on a first run's
      # plain pending, where the unlink would plant an NFS negative entry
      # and a lagging query could take a fresh sentinel. Done before the
      # other clusters' queries, so a rerun finishing meanwhile keeps its own.
      for handle in handles:
        if handle.cluster != cluster:
          continue
        consumed.add(handle)              # its cluster has answered
        entry = live[cluster].get(handle.scheduler_job_id)
        if entry is None or entry[0] != "pending":
          continue
        if self.is_requeue_state(entry[1]) or handle in reattached \
           or phases.get(handle) in ("running", "exited"):
          self._remove_sentinel(handle)
    handle_set = set(handles)
    handles = [h for h in handles if h.cluster in live]
    with self._cache_lock:
      # Spend only the reattach marks this refresh saw whose cluster
      # answered; one set by an adopt in flight since stays.
      self._reattached -= consumed & reattached
    if not handles:
      return
    results = {}
    candidates = []
    came_back_live = set()
    reset_first = set()
    for handle in handles:
      entry = live[handle.cluster].get(handle.scheduler_job_id)
      if entry is not None and entry[0] in ("pending", "running"):
        results[handle] = PollResult(
          entry[0], reason=self._reason_text(entry[1], entry[2]))
        came_back_live.add(handle)
        continue
      raw = entry[1] if entry is not None else None
      # The sentinel is authoritative but never read while the job is
      # live: an early miss would leave an NFS negative-cache entry.
      code = self._read_sentinel(handle)
      if code is not None:
        results[handle] = PollResult(
          "exited", exit_code=code, reason=raw,
          cancelled=self._cancelled_flag(raw))
        continue
      candidates.append((handle, entry))
    rows, accounting_timed_out = self._accounting_rows(candidates)
    new_first = {}
    absent_grace = max(2.0 * self.poll_interval, 30.0)
    finishing = PollResult("running",
                           reason="finishing: waiting for the exit code")
    for handle, entry in candidates:
      first = first_not_live.get(handle)
      if first is None:
        first = new_first.setdefault(handle, started)
      grace_over = started - first >= self.sentinel_grace
      row = rows.get(handle)
      if row is not None and self.state_phase(row[1]) == "exited":
        # Accounting reports the batch script's code, the program's only
        # if {{command}} ends the script; give the authoritative sentinel
        # sentinel_grace before trusting a completed or failed row.
        if grace_over or not (self.is_completed_state(row[1])
                              or self.is_failed_state(row[1])):
          results[handle] = self._from_accounting(handle, row)
        else:
          results[handle] = finishing
        continue
      if entry is None and started - submitted.get(handle, 0.0) < absent_grace:
        new_first.pop(handle, None)
        reset_first.add(handle)
        results[handle] = PollResult("pending")
        continue
      if not grace_over:
        results[handle] = finishing
        continue
      if handle.cluster in accounting_timed_out \
         and (entry is None or self.is_completed_state(entry[1])
              or self.is_failed_state(entry[1])) \
         and started - first < 2.0 * self.sentinel_grace:
        # Accounting timed out but could still add a purged job's state or
        # an exit code: wait at most another sentinel_grace. Cancelled and
        # scheduler-killed jobs resolve the same way without it.
        results[handle] = PollResult(
          "running", reason="finishing: waiting for the accounting query")
        continue
      results[handle] = self._fallback(handle, entry)
    with self._cache_lock:
      for handle, result in results.items():
        self._apply_locked(handle, result, started,
                           from_live=handle in came_back_live)
      for handle, when in new_first.items():
        # Not when a newer query already saw the handle live again.
        if self._stamps.get(handle, 0.0) <= started:
          self._first_not_live.setdefault(handle, when)
      for handle in came_back_live | reset_first:
        self._first_not_live.pop(handle, None)
      # Prune only here: a forced refresh (adopt) may see a partial set.
      # Entries written by a newer query are left alone.
      if force:
        return
      for handle in list(self._cache):
        if handle not in handle_set and self._cache[handle].phase == "exited" \
           and self._stamps.get(handle, 0.0) < started:
          for table in (self._cache, self._stamps, self._submitted,
                        self._first_not_live):
            table.pop(handle, None)

  @staticmethod
  def _by_cluster(handles):
    """``[(cluster, sorted ids), ...]`` for one scheduler query per cluster."""
    groups = {}
    for handle in handles:
      groups.setdefault(handle.cluster, set()).add(handle.scheduler_job_id)
    return [(cluster, sorted(ids))
            for cluster, ids in sorted(groups.items(),
                                       key=lambda kv: kv[0] or "")]

  def _accounting_rows(self, candidates):
    """``(rows by handle, clusters that timed out)``. A timeout is
    transient; any other failure turns accounting off for the process."""
    rows = {}
    timed_out = set()
    if not candidates or not self.accounting_available:
      return rows, timed_out
    handles = [h for h, _ in candidates]
    for cluster, ids in self._by_cluster(handles):
      try:
        found = self.query_terminal(ids, cluster=cluster)
      except SchedulerTimeout as exc:
        self._note_error("query_terminal", exc, cluster)
        timed_out.add(cluster)
        continue
      except Exception as exc:
        self.accounting_available = False
        _warn("accounting query failed (%s); relying on the exit-code "
              "sentinel from now on" % exc)
        break
      self._clear_error("query_terminal", cluster)
      for handle in handles:
        if handle.cluster == cluster and handle.scheduler_job_id in found:
          rows[handle] = found[handle.scheduler_job_id]
    return rows, timed_out

  def _apply_locked(self, handle, result, started, from_live=False):
    """Cache ``result`` unless a newer query wrote the entry. ``exited``
    is sticky except against a scheduler-reported ``pending`` (a requeue):
    once the manager stops passing a finished job, an out-of-order
    ``running`` would leave an entry the pruning rule never removes."""
    stamp = self._stamps.get(handle)
    if stamp is not None and stamp > started:
      return
    previous = self._cache.get(handle)
    if previous is not None and previous.phase == "exited" \
       and not (result.phase == "pending" and from_live):
      return
    self._cache[handle] = result
    self._stamps[handle] = started

  def _from_accounting(self, handle, row):
    """Resolve a terminal accounting row.

    Only a completed or failed row carries the program's exit code; a
    scheduler-killed job shows ``0:0`` and resolves as in ``_fallback``.
    """
    code, raw, reason = row
    if self.is_completed_state(raw) or self.is_failed_state(raw):
      if self.is_failed_state(raw) and not code:
        code = 1
      return PollResult("exited", exit_code=code,
                        reason=self._reason_text(raw, reason),
                        cancelled=False)
    return self._fallback(handle, ("exited", raw, reason))

  def _fallback(self, handle, entry):
    if entry is None:
      return PollResult("exited", reason="unknown to scheduler")
    raw = entry[1]
    if self.is_completed_state(raw):
      return PollResult("exited", exit_code=0, reason=raw, cancelled=False)
    if self.is_failed_state(raw):
      return PollResult("exited", exit_code=1,
                        reason="%s, exit code unavailable" % raw,
                        cancelled=False)
    if self.is_cancelled_state(raw):
      return PollResult("exited", reason=raw, cancelled=True)
    hint = self._hint(handle, raw)
    reason = "%s: %s" % (raw, hint) if hint else raw
    return PollResult("exited", reason=reason, cancelled=False)

  def _hint(self, handle, raw_state):
    try:
      with open(os.path.join(handle.job_dir, "scheduler.log"),
                encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    except OSError:
      return None
    # Only the last attempt (after the last rerun marker) explains the end.
    for index in range(len(lines) - 1, -1, -1):
      if lines[index].startswith(ATTEMPT_MARKER):
        lines = lines[index + 1:]
        break
    return self.failure_hint("\n".join(lines[-20:]), raw_state)

  @staticmethod
  def _reason_text(raw, reason):
    if not raw:
      return reason or None
    return "%s: %s" % (raw, reason) if reason else raw

  def _cancelled_flag(self, raw):
    if raw is None:
      return None
    return bool(self.is_cancelled_state(raw))

  def _note_error(self, what, exc, cluster=None):
    """Warn once per streak of failures of one query kind and cluster."""
    if cluster is not None:
      what = "%s (cluster %s)" % (what, cluster)
    if (what, cluster) not in self._error_streaks:
      _warn("%s: %s" % (what, exc))
      self._error_streaks.add((what, cluster))

  def _clear_error(self, what, cluster=None):
    if cluster is not None:
      what = "%s (cluster %s)" % (what, cluster)
    self._error_streaks.discard((what, cluster))
