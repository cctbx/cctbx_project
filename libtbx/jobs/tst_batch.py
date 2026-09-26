"""BatchBackend tests with an in-memory fake scheduler adapter."""

import os
import subprocess
import sys
import tempfile
import time


def _spec(tmp, **kw):
  from libtbx.jobs.manager import JobSpec
  kw.setdefault("argv", ["true"])
  kw.setdefault("cwd", tmp)
  kw.setdefault("log_path", os.path.join(tmp, "job.log"))
  return JobSpec(**kw)


class FakeAdapter(object):
  """Mixed into BatchBackend by ``_backend``; an in-memory scheduler.

  ``sched`` maps scheduler id to ``(phase, raw_state, reason)`` as
  query_live returns it; ``rows`` maps id to ``(exit_code, raw_state,
  reason)`` as query_terminal returns it.
  """

  name = "fake"
  directive_prefix = "#FAKE"
  script_extension = "fake"
  _EXITED = ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT")

  @classmethod
  def default_template_path(cls):
    return os.path.join(cls._template_dir, "default.fake")

  def submit_command(self, job, script_path):
    return list(self.submit_argv) + [script_path]

  def parse_submit_output(self, stdout, stderr):
    self.counter += 1
    if self.submit_cluster is not None:
      return str(self.counter), self.submit_cluster
    return str(self.counter)

  def query_live(self, ids, cluster=None):
    self.live_queries.append(list(ids) if cluster is None
                             else (cluster, list(ids)))
    if self.live_raises or cluster in self.live_raises_for:
      raise RuntimeError("squeue down%s" % (" on %s" % cluster if cluster else ""))
    if self.live_gate is not None:
      self.live_entered.set()
      self.live_gate.wait(10)
    table = self.sched if cluster is None else self.sched_by_cluster[cluster]
    return {i: table[i] for i in ids if i in table}

  def query_terminal(self, ids, cluster=None):
    self.terminal_queries.append(list(ids) if cluster is None
                                 else (cluster, list(ids)))
    if self.terminal_raises:
      raise RuntimeError("no sacct")
    if self.terminal_times_out:
      from libtbx.jobs.backends.batch import SchedulerTimeout
      raise SchedulerTimeout("sacct did not finish within 1 s")
    table = self.rows if cluster is None else self.rows_by_cluster[cluster]
    return {i: table[i] for i in ids if i in table}

  def state_phase(self, raw_state):
    if raw_state in ("PENDING",):
      return "pending"
    if raw_state in self._EXITED:
      return "exited"
    return "running"

  def is_cancelled_state(self, raw_state):
    return raw_state == "CANCELLED"

  def is_completed_state(self, raw_state):
    return raw_state == "COMPLETED"

  def is_failed_state(self, raw_state):
    return raw_state == "FAILED"

  def cancel_command(self, job):
    handle = job.handle
    self.cancelled.append(handle.scheduler_job_id if handle.cluster is None
                          else (handle.cluster, handle.scheduler_job_id))
    return list(self.cancel_argv)

  def resource_directives(self, resources):
    return ["#FAKE --%s=%s" % (k, v) for k, v in sorted(resources.items())]

  def failure_hint(self, text, raw_state):
    if raw_state != "TIMEOUT":
      return None
    for line in text.splitlines():
      if "TIME LIMIT" in line:
        return line.strip()
    return None

  def is_requeue_state(self, raw_state):
    return raw_state == "REQUEUED"

  def attempt_variable(self):
    return "FAKE_RESTART_COUNT"


def _backend(**kw):
  """Build a FakeAdapter-backed BatchBackend."""
  from libtbx.jobs.backends.batch import BatchBackend

  class Fake(FakeAdapter, BatchBackend):
    pass

  Fake._template_dir = tempfile.mkdtemp()
  with open(os.path.join(Fake._template_dir, "default.fake"), "w") as fh:
    fh.write("#!/bin/bash\n{{resource_directives}}\n{{command}}\n")
  b = Fake(**kw)
  b.counter = 0
  b.sched = {}
  b.rows = {}
  b.live_queries = []
  b.terminal_queries = []
  b.live_raises = False
  b.live_raises_for = set()
  b.live_gate = None
  b.live_entered = None
  b.terminal_raises = False
  b.terminal_times_out = False
  b.sched_by_cluster = {}
  b.rows_by_cluster = {}
  b.submit_cluster = None
  b.cancelled = []
  b.submit_argv = ["true"]
  b.cancel_argv = ["true"]
  return b


def _expect_template_error(fn, needle):
  from libtbx.jobs.backends.batch import TemplateError
  try:
    fn()
  except TemplateError as exc:
    assert needle in str(exc), (needle, str(exc))
  else:
    raise AssertionError("expected TemplateError containing %r" % needle)


def exercise_brace_template_leaves_shell_syntax_alone():
  from libtbx.jobs.backends.batch import BraceTemplate, referenced_placeholders
  text = ('echo "$SLURM_JOB_ID ${HOME} $(hostname) $? {{ partition }}"\n'
          "{{command}}\n")
  assert referenced_placeholders(text) == {"partition", "command"}
  out = BraceTemplate(text).substitute({"partition": "gpu", "command": "run"})
  assert out == 'echo "$SLURM_JOB_ID ${HOME} $(hostname) $? gpu"\nrun\n', out


def exercise_render_defaults_resources_and_takeover():
  tmp = tempfile.mkdtemp()
  b = _backend(template_text=(
    "#!/bin/bash\n#FAKE --partition={{partition}}\n#FAKE --account=x\n"
    "#FAKE --gres=gpu:{{gpus}}\n{{resource_directives}}\n{{environment}}\n"
    "{{command}}\n"), defaults={"partition": "normal"})
  spec = _spec(tmp, resources={"cpus": 4, "gpus": 2, "time_minutes": 90,
                               "partition": "gpu"},
               env={"A": "it's"})
  text = b.render(spec)
  assert "#FAKE --partition=gpu" in text, text
  assert "#FAKE --gres=gpu:2" in text, text
  assert "#FAKE --gpus=" not in text, "gpus taken over by the template"
  assert "#FAKE --cpus=4\n#FAKE --time_minutes=90" in text, text
  assert "export A='it'\"'\"'s'" in text, text
  assert text.rstrip().endswith("bash %s" % os.path.join(tmp, ".libtbx_jobs",
                                                          "j_render", "run.sh")), text
  # site default alone
  plain = b.render(_spec(tmp, resources={"gpus": 1}))
  assert "#FAKE --partition=normal" in plain
  # time placeholder formats minutes
  b2 = _backend(template_text="#FAKE --time={{time}}\n{{command}}\n")
  assert "#FAKE --time=1:30:00" in b2.render(_spec(tmp, resources={"time_minutes": 90}))
  # Every resource placeholder takes its directive over, and only its own.
  every = {"cpus": 3, "mem_mb": 512, "time_minutes": 45, "gpus": 1}
  for name, directive in (("cpus", "--cpus=3"), ("mem_mb", "--mem_mb=512"),
                          ("time_minutes", "--time_minutes=45"),
                          ("time", "--time_minutes=45"), ("gpus", "--gpus=1")):
    b3 = _backend(template_text="#FAKE --site-%s={{%s}}\n{{resource_directives}}\n"
                                "{{command}}\n" % (name, name))
    text = b3.render(_spec(tmp, resources=dict(every)))
    assert "#FAKE %s" % directive not in text, (name, text)
    others = {"--cpus=3", "--mem_mb=512", "--time_minutes=45", "--gpus=1"} - {directive}
    for other in others:
      assert "#FAKE %s" % other in text, (name, other, text)
    _expect_template_error(lambda: b3.validate(_spec(tmp)), name)


def exercise_render_errors():
  tmp = tempfile.mkdtemp()
  b = _backend(template_text="{{nope}}\n{{command}}\n")
  _expect_template_error(lambda: b.validate(_spec(tmp)), "nope")
  b = _backend(template_text="#FAKE --cpus={{cpus}}\n{{command}}\n")
  _expect_template_error(lambda: b.validate(_spec(tmp)), "cpus")
  b.validate(_spec(tmp, resources={"cpus": 2}))
  _expect_template_error(
    lambda: _backend(template_text="echo {{a b}}\n{{command}}\n"), "{{")
  for text in ("#!/bin/bash\n# runs {{command}}\n{{command}}\n",
               "  # {{resource_directives}} goes below\n{{command}}\n",
               "##FAKE --cpus={{cpus}}\n{{command}}\n"):
    _expect_template_error(lambda: _backend(template_text=text),
                           "not allowed in comment lines")
  _backend(template_text="#!/bin/bash\n#FAKE --p={{partition}}\n"
                         "# a comment\n{{command}}\n")
  b = _backend(template_text="echo hi\n")
  _expect_template_error(lambda: b.validate(_spec(tmp)), "command")
  assert _backend(command_timeout=1e300).command_timeout == 86400.0
  assert _backend().requeue is None and _backend(requeue=False).requeue is False
  for bad in (dict(requeue="false"), dict(requeue=0), dict(requeue=1),
              dict(poll_interval=float("inf")), dict(sentinel_grace=float("nan")),
              dict(command_timeout=float("inf")),
              dict(defaults={"command": "x"}), dict(defaults={"time_minutes": 5}),
              dict(defaults={"partition": "a b"}), dict(defaults={"gpus": True}),
              dict(template="/no/such/file"), dict(template="/x", template_text="y")):
    try:
      _backend(**bad)
    except (ValueError, OSError):
      pass
    else:
      raise AssertionError("constructor accepted %r" % (bad,))
  b = _backend()
  for bad in ({"resources": {"time": "1:00:00"}}, {"resources": {"job_id": "x"}}):
    try:
      b.validate(_spec(tmp, **bad))
    except ValueError as exc:
      assert "reserved" in str(exc), str(exc)
    else:
      raise AssertionError("validate accepted %r" % (bad,))
  for label, value in (("cwd", "/tmp/with space"), ("log_path", "/tmp/a'b"),
                       ("cwd", "/tmp/dir\n")):
    try:
      b.validate(_spec(tmp, **{label: value}))
    except ValueError as exc:
      assert label in str(exc), str(exc)
    else:
      raise AssertionError("validate accepted %s=%r" % (label, value))
  try:
    _spec(tmp, env={"A;cmd": "1"})
  except ValueError:
    pass
  else:
    raise AssertionError("env key with shell metacharacters accepted")
  for key in ("LIBTBX_JOBS_ATTEMPT", "libtbx_jobs_attempt", "FAKE_RESTART_COUNT"):
    try:
      b.validate(_spec(tmp, env={key: "9"}))
    except ValueError as exc:
      assert "attempt count" in str(exc), str(exc)
    else:
      raise AssertionError("env key %s accepted; it steers the attempt count" % key)
  # a value containing {{ is data, not a stray placeholder
  b_env = _backend(template_text="{{environment}}\n{{command}}\n")
  text = b_env.render(_spec(tmp, env={"X": "{{"}))
  assert "export X='{{'" in text, text
  try:
    b.validate(_spec(tmp, resources={"partition": "a\nb"}))
  except ValueError:
    pass
  else:
    raise AssertionError("newline in a placeholder value accepted")


def exercise_job_name_is_sanitised():
  tmp = tempfile.mkdtemp()
  b = _backend()
  assert b.job_name(_spec(tmp, name="phenix.refine run 1")) == "phenix.refine_run_1"
  assert b.job_name(_spec(tmp, name="x" * 80)) == "x" * 64
  assert b.job_name(_spec(tmp, name="")) == "job"


def _job(b, spec, job_id="j_000000000001"):
  from libtbx.jobs.manager import Job
  return Job(job_id=job_id, spec=spec)


def exercise_wrapper_runs_the_command_and_writes_the_sentinel():
  tmp = tempfile.mkdtemp()
  b = _backend()
  weird = "it's $HOME and spaces"
  spec = _spec(tmp, argv=[sys.executable, "-c",
                          "import sys; print(sys.argv[1]); sys.exit(7)", weird],
               env={"LIBTBX_JOBS_X": "1 2"})
  job = _job(b, spec)
  phase, handle = b.submit(job)
  job.handle = handle
  assert phase == "pending" and handle.scheduler_job_id == "1"
  assert handle.job_dir == os.path.join(tmp, ".libtbx_jobs", job.job_id)
  run_sh = os.path.join(handle.job_dir, "run.sh")
  script = os.path.join(handle.job_dir, "job.fake")
  assert os.path.exists(run_sh) and os.path.exists(script)
  assert os.stat(run_sh).st_mode & 0o777 == 0o700
  with open(script) as fh:
    assert fh.read().rstrip().endswith("bash %s" % run_sh)
  _write_sentinel(handle, 3)                   # a previous attempt's code
  rc = subprocess.call(["bash", run_sh])
  assert rc == 7, rc
  # The wrapper clears a previous attempt's sentinel before the command
  # runs, so a requeued job never reports the abandoned attempt's code.
  sentinel = os.path.join(tmp, ".libtbx_jobs", "j_000000000003", "exit_code")
  spec3 = _spec(tmp, argv=["sh", "-c", "test ! -e %s" % sentinel],
                log_path=os.path.join(tmp, "third.log"))
  phase, handle3 = b.submit(_job(b, spec3, job_id="j_000000000003"))
  _write_sentinel(handle3, 3)
  assert subprocess.call(["bash", os.path.join(handle3.job_dir, "run.sh")]) == 0, \
    "the previous attempt's sentinel was still there when the command ran"
  assert b._read_sentinel(handle3) == 0
  # A requeued attempt: the wrapper exports the attempt count and marks
  # the log; a first attempt does neither.
  spec4 = _spec(tmp, argv=["sh", "-c", "echo attempt=$LIBTBX_JOBS_ATTEMPT"],
                log_path=os.path.join(tmp, "fourth.log"))
  phase, handle4 = b.submit(_job(b, spec4, job_id="j_000000000004"))
  env = dict(os.environ)
  env.pop("FAKE_RESTART_COUNT", None)
  subprocess.call(["bash", os.path.join(handle4.job_dir, "run.sh")], env=env)
  with open(os.path.join(tmp, "fourth.log")) as fh:
    first = fh.read()
  assert first == "attempt=0\n", first
  env["FAKE_RESTART_COUNT"] = "1"
  proc = subprocess.run(["bash", os.path.join(handle4.job_dir, "run.sh")],
                        env=env, capture_output=True, text=True)
  with open(os.path.join(tmp, "fourth.log")) as fh:
    text = fh.read()
  marker = "libtbx.jobs: attempt 2 of this job after a requeue"
  assert text == "attempt=0\n\n%s\nattempt=1\n" % marker, text
  # The scheduler log (the wrapper's own stdout) gets the marker too.
  assert proc.stdout == "\n%s\n" % marker, proc.stdout
  # A host's env can neither override the attempt count nor steer it
  # through the scheduler's own variable.
  spec5 = _spec(tmp, argv=["sh", "-c", "echo attempt=$LIBTBX_JOBS_ATTEMPT"],
                log_path=os.path.join(tmp, "fifth.log"),
                env={"LIBTBX_JOBS_ATTEMPT": "42", "FAKE_RESTART_COUNT": "7"})
  phase, handle5 = b.submit(_job(b, spec5, job_id="j_000000000005"))
  subprocess.call(["bash", os.path.join(handle5.job_dir, "run.sh")], env=env)
  with open(os.path.join(tmp, "fifth.log")) as fh:
    assert fh.read().endswith("attempt=1\n")
  with open(os.path.join(tmp, "job.log")) as fh:
    assert fh.read() == weird + "\n"
  assert b._read_sentinel(handle) == 7
  assert b.poll(job).phase == "pending"
  # A cwd that vanished makes the wrapper exit 127, and its message
  # survives shell metacharacters in the path unharmed.
  import shutil
  gone = os.path.join(tmp, "a$b(c)")
  os.makedirs(gone)
  spec2 = _spec(tmp, cwd=gone, log_path=os.path.join(tmp, "gone.log"))
  job2 = _job(b, spec2, job_id="j_000000000002")
  phase, handle2 = b.submit(job2)
  copy = os.path.join(tmp, "run_copy.sh")
  shutil.copy(os.path.join(handle2.job_dir, "run.sh"), copy)
  shutil.rmtree(gone)
  proc = subprocess.run(["bash", copy], capture_output=True, text=True)
  assert proc.returncode == 127, proc
  assert "cannot chdir to" in proc.stderr and "a$b(c)" in proc.stderr, proc.stderr


def exercise_submit_failures_become_submit_errors():
  from libtbx.jobs.backends.base import SubmitError
  tmp = tempfile.mkdtemp()
  b = _backend()
  b.submit_argv = ["false"]
  try:
    b.submit(_job(b, _spec(tmp)))
  except SubmitError:
    pass
  else:
    raise AssertionError("failing submit command accepted")
  b.submit_argv = ["/nonexistent/sbatch"]
  try:
    b.submit(_job(b, _spec(tmp)))
  except SubmitError as exc:
    assert "nonexistent" in str(exc), str(exc)
  else:
    raise AssertionError("missing submit command accepted")
  slow = _backend(command_timeout=1.0)
  slow.submit_argv = ["sh", "-c", "sleep 5", "_"]
  t0 = time.monotonic()
  try:
    slow.submit(_job(slow, _spec(tmp)))
  except SubmitError as exc:
    assert "did not finish within 1 s" in str(exc), str(exc)
    assert "may still have accepted" in str(exc), str(exc)
  else:
    raise AssertionError("hung submit command accepted")
  assert time.monotonic() - t0 < 4, "timeout did not kill the submit command"


def exercise_terminate_shutdown_handle_info_and_reattach():
  from libtbx.jobs.backends.batch import BatchHandle
  from libtbx.utils import Sorry
  import io
  tmp = tempfile.mkdtemp()
  b = _backend()
  job = _job(b, _spec(tmp))
  assert b.handle_info(job) == {"scheduler": "fake", "scheduler_job_id": None,
                                "job_dir": None, "scheduler_cluster": None}
  assert set(b.handle_info(job)) == set(b.handle_keys)
  phase, handle = b.submit(job)
  job.handle = handle
  info = b.handle_info(job)
  assert info == {"scheduler": "fake", "scheduler_job_id": "1",
                  "job_dir": handle.job_dir, "scheduler_cluster": None}, info
  b.terminate(job)
  assert b.cancelled == ["1"]
  b.cancel_argv = ["sh", "-c", "sleep 5"]
  b.command_timeout = 1.0
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  t0 = time.monotonic()
  try:
    b.terminate(job)
  finally:
    sys.stderr = real
  assert time.monotonic() - t0 < 4, "timeout did not kill the cancel command"
  assert "did not finish within 1 s" in err.getvalue(), err.getvalue()
  b.cancel_argv = ["true"]
  assert b.kill(job) is None
  # A missing stderr: the note is dropped, never sent to stdout.
  import io as _io
  b.cancel_argv = ["/nonexistent/scancel"]
  sys.stderr = None
  out = _io.StringIO()
  real_out = sys.stdout
  sys.stdout = out
  try:
    b.terminate(job)
  finally:
    sys.stdout = real_out
    sys.stderr = real
  assert out.getvalue() == "", out.getvalue()
  b.cancel_argv = ["true"]
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.on_shutdown([job])
  finally:
    sys.stderr = real
  assert "1" in err.getvalue() and handle.job_dir in err.getvalue(), err.getvalue()
  again = b.reattach(info)
  assert again == handle
  assert b.reattach(handle) == handle
  assert b._submitted[handle] > 0
  for bad in ({"scheduler_job_id": "1"}, {"scheduler_job_id": "", "job_dir": "/x"},
              {"scheduler_job_id": "1", "job_dir": "relative"}, 42,
              {"scheduler_job_id": "--user=alice", "job_dir": "/x"},
              {"scheduler_job_id": "1 2", "job_dir": "/x"},
              {"scheduler_job_id": None, "job_dir": "/x"},
              {"scheduler_job_id": "1", "job_dir": "/x", "scheduler_cluster": "-M"},
              BatchHandle("-9", "/x")):
    try:
      b.reattach(bad)
    except Sorry as exc:
      assert "cannot reattach" in str(exc), str(exc)
    else:
      raise AssertionError("reattach accepted %r" % (bad,))
  assert b.reattach(BatchHandle("9", "/x")).scheduler_job_id == "9"
  assert b.reattach(BatchHandle("123_4", "/x")).scheduler_job_id == "123_4"
  with_cluster = b.reattach({"scheduler_job_id": "7", "job_dir": "/x",
                             "scheduler_cluster": "gpu"})
  assert with_cluster == BatchHandle("7", "/x", "gpu"), with_cluster
  # An empty persisted field (a Django CharField) means no cluster.
  assert b.reattach({"scheduler_job_id": "7", "job_dir": "/x",
                     "scheduler_cluster": ""}) == BatchHandle("7", "/x")


def _submitted(b, tmp, n):
  jobs = []
  for i in range(n):
    job = _job(b, _spec(tmp, log_path=os.path.join(tmp, "job%d.log" % i)),
               job_id="j_%012d" % (i + 1))
    phase, job.handle = b.submit(job)
    jobs.append(job)
  return jobs


def _write_sentinel(handle, code):
  with open(os.path.join(handle.job_dir, "exit_code"), "w") as fh:
    fh.write("%d\n" % code)


def exercise_refresh_rate_limit_empty_list_and_live_states():
  tmp = tempfile.mkdtemp()
  b = _backend(poll_interval=1.0)
  jobs = _submitted(b, tmp, 2)
  a, c = jobs
  b.refresh([])
  assert b.live_queries == [], "an empty refresh ran a scheduler command"
  b.sched = {"1": ("pending", "PENDING", "Priority"), "2": ("running", "RUNNING", None)}
  b.refresh(jobs)
  assert b.live_queries == [["1", "2"]]
  assert b.poll(a).phase == "pending" and b.poll(a).reason == "PENDING: Priority"
  assert b.poll(c).phase == "running" and b.poll(c).reason == "RUNNING"
  b.sched["1"] = ("running", "RUNNING", None)
  b.refresh(jobs)
  assert len(b.live_queries) == 1, "refresh ignored the poll interval"
  assert b.poll(a).phase == "pending"
  b.refresh(jobs, force=True)
  assert b.poll(a).phase == "running"


def exercise_sentinel_is_read_only_after_the_scheduler_reports_terminal():
  tmp = tempfile.mkdtemp()
  b = _backend()
  (job,) = _submitted(b, tmp, 1)
  _write_sentinel(job.handle, 0)
  b.sched = {"1": ("running", "RUNNING", None)}
  b.refresh([job], force=True)
  assert b.poll(job).phase == "running", "sentinel consulted while live"
  b.sched = {"1": ("exited", "COMPLETED", None)}
  b.refresh([job], force=True)
  r = b.poll(job)
  assert (r.phase, r.exit_code, r.reason, r.cancelled) == ("exited", 0, "COMPLETED", False), r
  assert b.terminal_queries == [], "accounting asked although the sentinel answered"


def exercise_accounting_resolves_and_stale_live_rows_do_not():
  tmp = tempfile.mkdtemp()
  b = _backend()
  a, c = _submitted(b, tmp, 2)
  b.sched = {}
  b.rows = {"1": (3, "FAILED", None), "2": (None, "RUNNING", None)}
  b._submitted[a.handle] -= 100
  b._submitted[c.handle] -= 100
  b.refresh([a, c], force=True)
  ra, rc = b.poll(a), b.poll(c)
  # A failed row carries an exit code the sentinel may still contradict:
  # it is trusted only once sentinel_grace has passed.
  assert ra.phase == "running" and ra.reason.startswith("finishing"), ra
  assert rc.phase == "running" and rc.reason.startswith("finishing"), rc
  with b._cache_lock:
    b._first_not_live[a.handle] -= 100
  b.refresh([a, c], force=True)
  ra = b.poll(a)
  assert (ra.phase, ra.exit_code, ra.reason, ra.cancelled) == ("exited", 3, "FAILED", False), ra
  b.rows["2"] = (None, "CANCELLED", "signal 15")
  b.refresh([a, c], force=True)
  rc = b.poll(c)
  assert rc.phase == "exited" and rc.cancelled is True and rc.reason == "CANCELLED", rc


def exercise_accounting_code_yields_to_a_sentinel_within_the_grace():
  # The batch script ends with `echo done`, so accounting reports 0 for a
  # program that exited 2; the sentinel arrives late (NFS) and must win.
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=0.3)
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("exited", "COMPLETED", None)}
  b.rows = {"1": (0, "COMPLETED", None)}
  b.refresh([job], force=True)
  r = b.poll(job)
  assert r.phase == "running" and r.reason.startswith("finishing"), r
  assert b.terminal_queries == [["1"]], b.terminal_queries
  _write_sentinel(job.handle, 2)
  b.refresh([job], force=True)
  r = b.poll(job)
  assert (r.phase, r.exit_code, r.reason, r.cancelled) == ("exited", 2, "COMPLETED", False), r
  # Cancelled and scheduler-killed rows carry no program exit code and
  # resolve at once, without waiting for a sentinel that will never come.
  k, x = _submitted(b, tempfile.mkdtemp(), 2)
  b.sched = {}
  b.rows = {"2": (0, "CANCELLED", "signal 15"), "3": (0, "TIMEOUT", None)}
  for j in (k, x):
    b._submitted[j.handle] -= 100
  b.refresh([k, x], force=True)
  rk, rx = b.poll(k), b.poll(x)
  assert (rk.phase, rk.exit_code, rk.cancelled, rk.reason) == ("exited", None, True, "CANCELLED"), rk
  assert (rx.phase, rx.exit_code, rx.reason) == ("exited", None, "TIMEOUT"), rx


def exercise_a_requeue_removes_the_abandoned_attempts_sentinel():
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=0.0)
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("running", "RUNNING", None)}
  b.refresh([job], force=True)
  _write_sentinel(job.handle, 3)              # the attempt that died
  b.sched = {"1": ("pending", "PENDING", None)}   # a lagging query, or a requeue
  b.refresh([job], force=True)
  assert b.poll(job).phase == "pending"
  assert b._read_sentinel(job.handle) is None, "the stale sentinel survived"
  # The rerun dies before its wrapper starts: no sentinel, and the raw
  # state decides rather than the abandoned attempt's code.
  b.sched = {"1": ("exited", "NODE_FAIL", None)}
  b.refresh([job], force=True)
  r = b.poll(job)
  assert (r.phase, r.exit_code, r.reason) == ("exited", None, "NODE_FAIL"), r
  # exited to pending (a requeue after the scheduler reported the end):
  # the sentinel goes, and the entry is allowed back to pending.
  (again,) = _submitted(b, tempfile.mkdtemp(), 1)
  _write_sentinel(again.handle, 0)
  b.sched = {"2": ("exited", "COMPLETED", None)}
  b.refresh([again], force=True)
  assert b.poll(again).phase == "exited"
  b.sched = {"2": ("pending", "PENDING", None)}   # plain pending after exited
  b.refresh([again], force=True)
  assert b.poll(again).phase == "pending"
  assert b._read_sentinel(again.handle) is None, "the stale sentinel survived"
  # An adopted handle first seen pending: the backend never saw it run,
  # and any sentinel it finds is still stale.
  adopted = b.reattach({"scheduler_job_id": "77", "job_dir": tempfile.mkdtemp()})
  os.makedirs(adopted.job_dir, exist_ok=True)
  _write_sentinel(adopted, 5)
  from libtbx.jobs.manager import Job
  holder = Job(job_id="j_00000000adop", spec=_spec(tmp), handle=adopted)
  # The mark is not spent by a refresh in which its cluster did not
  # answer, nor by one that did not carry the handle.
  b.live_raises = True
  b.refresh([holder], force=True)
  b.live_raises = False
  b.refresh([job], force=True)
  assert adopted in b._reattached, "the reattached mark was spent early"
  b.sched = {"77": ("pending", "PENDING", "Priority")}
  b.refresh([holder], force=True)
  assert b.poll(holder).phase == "pending"
  assert b._read_sentinel(adopted) is None, "an adopted handle's stale sentinel survived"
  assert adopted not in b._reattached
  # A first run's plain pending never unlinks: a lagging query must not
  # take a fresh sentinel, and the unlink would plant a negative entry.
  (fresh,) = _submitted(b, tempfile.mkdtemp(), 1)
  _write_sentinel(fresh.handle, 4)
  b.sched = {"3": ("pending", "PENDING", "Priority")}
  b.refresh([fresh], force=True)
  assert b.poll(fresh).phase == "pending"
  assert b._read_sentinel(fresh.handle) == 4, "a first run's sentinel was removed"
  # A requeue state removes it even on a handle never seen live.
  b.sched = {"3": ("pending", "REQUEUED", None)}
  b.refresh([fresh], force=True)
  assert b._read_sentinel(fresh.handle) is None, "a requeued job's sentinel survived"
  # The removal happens once the handle's own cluster answers, so a rerun
  # finishing during another cluster's query keeps its fresh sentinel.
  b.submit_cluster = "gpu"
  routed = _job(b, _spec(tempfile.mkdtemp()), job_id="j_000000000009")
  phase, routed.handle = b.submit(routed)
  b.submit_cluster = None
  b.sched = {"3": ("pending", "REQUEUED", None)}
  b.sched_by_cluster = {"gpu": {"4": ("running", "RUNNING", None)}}
  _write_sentinel(fresh.handle, 1)
  seen = []
  real_live = b.query_live
  def write_during_gpu_query(ids, cluster=None):
    if cluster == "gpu":
      seen.append(b._read_sentinel(fresh.handle))     # removed by now
      _write_sentinel(fresh.handle, 9)                 # the rerun writes
    return real_live(ids, cluster=cluster)
  b.query_live = write_during_gpu_query
  b.refresh([fresh, routed], force=True)
  del b.query_live
  assert seen == [None], seen
  assert b._read_sentinel(fresh.handle) == 9, "a fresh sentinel was removed"


def exercise_late_sentinel_wins_over_the_fallback():
  # With accounting off, a FAILED job without a sentinel stays finishing
  # until the sentinel supplies the real code, not the fallback's 1.
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=60.0)
  b.terminal_raises = True
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("exited", "FAILED", None)}
  import io
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.refresh([job], force=True)
    b.refresh([job], force=True)
  finally:
    sys.stderr = real
  r = b.poll(job)
  assert r.phase == "running" and r.reason.startswith("finishing"), r
  _write_sentinel(job.handle, 5)
  b.refresh([job], force=True)
  r = b.poll(job)
  assert (r.phase, r.exit_code, r.reason, r.cancelled) == ("exited", 5, "FAILED", False), r


def exercise_accounting_timeout_is_transient():
  import io
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=60.0)
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("exited", "CANCELLED", None)}
  b.terminal_times_out = True
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.refresh([job], force=True)
    b.refresh([job], force=True)
  finally:
    sys.stderr = real
  assert b.accounting_available is True, "a timeout switched accounting off"
  assert err.getvalue().count("did not finish") == 1, err.getvalue()
  assert "finishing" in b.poll(job).reason
  # Past the grace, a timed-out accounting query does not hold a job the
  # scheduler killed: accounting could not change that answer.
  with b._cache_lock:
    b._first_not_live[job.handle] -= 100
  b.refresh([job], force=True)
  r = b.poll(job)
  assert (r.phase, r.cancelled, r.reason) == ("exited", True, "CANCELLED"), r
  # A purged handle is held, because accounting could still name its
  # state, but only for another sentinel_grace; then the fallback runs.
  (purged,) = _submitted(b, tempfile.mkdtemp(), 1)
  b.sched = {}
  with b._cache_lock:
    b._submitted[purged.handle] -= 100
    b._first_not_live[purged.handle] = time.monotonic() - 70
  b.refresh([purged], force=True)
  r = b.poll(purged)
  assert r.phase == "running" and "accounting" in r.reason, r
  with b._cache_lock:
    b._first_not_live[purged.handle] = time.monotonic() - 130
  b.refresh([purged], force=True)
  r = b.poll(purged)
  assert (r.phase, r.reason) == ("exited", "unknown to scheduler"), r
  # One cluster's timeout leaves another cluster's rows usable.
  b.terminal_times_out = False
  b.submit_cluster = "gpu"
  routed = _job(b, _spec(tempfile.mkdtemp()), job_id="j_000000000009")
  phase, routed.handle = b.submit(routed)
  b.sched_by_cluster = {"gpu": {}}
  b.rows_by_cluster = {"gpu": {"3": (0, "CANCELLED", "signal 15")}}
  real_query = b.query_terminal
  def timeout_default_only(ids, cluster=None):
    if cluster is None:
      from libtbx.jobs.backends.batch import SchedulerTimeout
      raise SchedulerTimeout("sacct did not finish within 1 s")
    return real_query(ids, cluster=cluster)
  b.query_terminal = timeout_default_only
  b.submit_cluster = None
  with b._cache_lock:
    b._submitted[routed.handle] -= 100
  (again,) = _submitted(b, tempfile.mkdtemp(), 1)
  with b._cache_lock:
    b._submitted[again.handle] -= 100
  b.refresh([routed, again], force=True)
  assert b.poll(routed).cancelled is True, b.poll(routed)
  assert "finishing" in b.poll(again).reason, b.poll(again)
  del b.query_terminal
  b.rows = {"1": (0, "CANCELLED", "signal 15")}
  b.refresh([job], force=True)
  r = b.poll(job)
  assert r.phase == "exited" and r.cancelled is True, r
  assert b.terminal_queries[-1] == ["1"], b.terminal_queries


def exercise_handles_on_other_clusters_are_queried_there():
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=0.2)
  (plain,) = _submitted(b, tmp, 1)
  b.submit_cluster = "gpu"
  job = _job(b, _spec(tmp, log_path=os.path.join(tmp, "gpu.log")),
             job_id="j_000000000002")
  phase, job.handle = b.submit(job)
  assert job.handle.cluster == "gpu" and job.handle.scheduler_job_id == "2"
  info = b.handle_info(job)
  assert info["scheduler_cluster"] == "gpu", info
  assert b.reattach(info) == job.handle
  b.sched = {"1": ("running", "RUNNING", None)}
  b.sched_by_cluster = {"gpu": {"2": ("pending", "PENDING", "Priority")}}
  b.refresh([plain, job], force=True)
  assert b.live_queries[-2:] == [["1"], ("gpu", ["2"])], b.live_queries
  assert b.poll(plain).phase == "running"
  assert b.poll(job).phase == "pending" and b.poll(job).reason == "PENDING: Priority"
  b.terminate(job)
  assert b.cancelled[-1] == ("gpu", "2"), b.cancelled
  import io
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.on_shutdown([plain, job])
  finally:
    sys.stderr = real
  assert "job 2 on cluster gpu" in err.getvalue(), err.getvalue()
  # Accounting is asked per cluster too, and the same id on two clusters
  # stays apart.
  b.sched_by_cluster = {"gpu": {}}
  b.sched = {}
  b.rows = {"2": (0, "COMPLETED", None)}
  b.rows_by_cluster = {"gpu": {"2": (0, "CANCELLED", "signal 15")}}
  b._submitted[job.handle] -= 100
  b.refresh([job], force=True)
  assert b.terminal_queries[-1] == ("gpu", ["2"]), b.terminal_queries
  r = b.poll(job)
  assert r.phase == "exited" and r.cancelled is True, r


def exercise_one_failing_cluster_does_not_freeze_the_others():
  import io
  tmp = tempfile.mkdtemp()
  b = _backend()
  (plain,) = _submitted(b, tmp, 1)
  b.submit_cluster = "oldgpu"
  routed = _job(b, _spec(tmp, log_path=os.path.join(tmp, "r.log")),
                job_id="j_000000000002")
  phase, routed.handle = b.submit(routed)
  b.sched = {"1": ("running", "RUNNING", None)}
  b.sched_by_cluster = {"oldgpu": {"2": ("running", "RUNNING", None)}}
  b.refresh([plain, routed], force=True)
  assert b.poll(routed).phase == "running"
  b.live_raises_for = {"oldgpu"}
  _write_sentinel(plain.handle, 0)
  b.sched = {"1": ("exited", "COMPLETED", None)}
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.refresh([plain, routed], force=True)
    b.refresh([plain, routed], force=True)
  finally:
    sys.stderr = real
  assert b.poll(plain).phase == "exited", "a failing cluster froze another"
  assert b.poll(routed).phase == "running", "a failed query changed the cache"
  assert err.getvalue().count("squeue down on oldgpu") == 1, err.getvalue()
  assert "(cluster oldgpu)" in err.getvalue(), err.getvalue()
  # The failing cluster's handle is never pruned, even when exited under
  # a full refresh, and resolves once its controller answers.
  assert routed.handle in b._cache
  b.live_raises_for = set()
  _write_sentinel(routed.handle, 0)
  b.sched_by_cluster = {"oldgpu": {"2": ("exited", "COMPLETED", None)}}
  b.refresh([plain, routed], force=True)
  assert b.poll(routed).phase == "exited"
  b.live_raises_for = {"oldgpu"}
  b._last_refresh = None
  b.refresh([plain, routed])
  assert routed.handle in b._cache and routed.handle in b._submitted, \
    "a failing cluster's exited entry was pruned"
  assert b.poll(routed).phase == "exited"
  b.live_raises_for = set()
  b._last_refresh = None
  b.refresh([plain])
  assert routed.handle not in b._cache, \
    "an exited entry absent from a full refresh survived"


def exercise_accounting_codes_only_count_for_completed_and_failed():
  # sacct reports 0:0 for a job the scheduler killed (TIMEOUT, NODE_FAIL,
  # ...); that zero is not the program's exit code.
  tmp = tempfile.mkdtemp()
  b = _backend(sentinel_grace=0.0)
  t, f, c = _submitted(b, tmp, 3)
  with open(os.path.join(t.handle.job_dir, "scheduler.log"), "w") as fh:
    fh.write("slurmstepd: error: *** JOB 1 CANCELLED DUE TO TIME LIMIT ***\n")
  b.sched = {}
  b.rows = {"1": (0, "TIMEOUT", None), "2": (0, "FAILED", None),
            "3": (0, "COMPLETED", None)}
  for job in (t, f, c):
    b._submitted[job.handle] -= 100
  b.refresh([t, f, c], force=True)
  rt, rf, rc = b.poll(t), b.poll(f), b.poll(c)
  assert rt.phase == "exited" and rt.exit_code is None, rt
  assert rt.reason.startswith("TIMEOUT") and "TIME LIMIT" in rt.reason, rt
  assert rt.cancelled is False, rt
  assert (rf.phase, rf.exit_code, rf.cancelled) == ("exited", 1, False), rf
  assert rf.reason.startswith("FAILED"), rf
  assert (rc.phase, rc.exit_code, rc.reason) == ("exited", 0, "COMPLETED"), rc
  b2 = _backend(sentinel_grace=0.0)
  (bare,) = _submitted(b2, tempfile.mkdtemp(), 1)
  b2.sched = {}
  b2.rows = {"1": (0, "TIMEOUT", None)}
  b2._submitted[bare.handle] -= 100
  b2.refresh([bare], force=True)
  r = b2.poll(bare)
  assert (r.phase, r.exit_code, r.reason) == ("exited", None, "TIMEOUT"), r


def exercise_grace_periods_and_fallbacks():
  import io
  tmp = tempfile.mkdtemp()
  b = _backend(poll_interval=1.0, sentinel_grace=0.2)
  a, c, d = _submitted(b, tmp, 3)
  b.terminal_raises = True
  b.sched = {}
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.refresh([a], force=True)
    assert b.poll(a).phase == "pending", "absent handle right after submission"
    assert b.accounting_available is False
    b._submitted[a.handle] -= 60
    b.refresh([a], force=True)
  finally:
    sys.stderr = real
  assert err.getvalue().count("no sacct") == 1, err.getvalue()
  assert b.terminal_queries == [["1"]], b.terminal_queries
  assert b.poll(a).phase == "running" and "finishing" in b.poll(a).reason
  time.sleep(0.25)
  b.refresh([a], force=True)
  r = b.poll(a)
  assert (r.phase, r.exit_code, r.reason) == ("exited", None, "unknown to scheduler"), r
  b.sched = {"2": ("exited", "COMPLETED", None), "3": ("exited", "TIMEOUT", None)}
  with open(os.path.join(d.handle.job_dir, "scheduler.log"), "w") as fh:
    # An earlier attempt's hint sits above the marker and must not be
    # read as the final attempt's.
    fh.write("slurmstepd: error: *** JOB 3 ON n1 CANCELLED DUE TO TIME LIMIT ***\n"
             "\nlibtbx.jobs: attempt 2 of this job after a requeue\n"
             "slurmstepd: error: *** JOB 3 ON n2 CANCELLED DUE TO TIME LIMIT ***\n")
  b.refresh([c, d], force=True)
  assert "finishing" in b.poll(c).reason
  time.sleep(0.25)
  b.refresh([c, d], force=True)
  rc, rd = b.poll(c), b.poll(d)
  assert (rc.phase, rc.exit_code, rc.reason) == ("exited", 0, "COMPLETED"), rc
  assert rd.exit_code is None and rd.reason.startswith("TIMEOUT: ") and "TIME LIMIT" in rd.reason, rd
  assert "ON n2" in rd.reason and "ON n1" not in rd.reason, rd
  (e,) = _submitted(b, tmp, 1)
  b.sched = {"4": ("exited", "FAILED", None)}
  b.refresh([e], force=True)
  time.sleep(0.25)
  b.refresh([e], force=True)
  re_ = b.poll(e)
  assert (re_.exit_code, re_.reason) == (1, "FAILED, exit code unavailable"), re_


def exercise_query_failure_keeps_the_cache():
  import io
  tmp = tempfile.mkdtemp()
  b = _backend()
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("running", "RUNNING", None)}
  b.refresh([job], force=True)
  b.live_raises = True
  err = io.StringIO()
  real = sys.stderr
  sys.stderr = err
  try:
    b.refresh([job], force=True)
    b.refresh([job], force=True)
  finally:
    sys.stderr = real
  assert b.poll(job).phase == "running"
  assert err.getvalue().count("squeue down") == 1, err.getvalue()
  # A streak ends with a success and a new failure is reported again;
  # streaks of the two query kinds are independent.
  b.live_raises = False
  b.refresh([job], force=True)
  b.terminal_times_out = True
  b.sched = {"1": ("exited", "COMPLETED", None)}
  sys.stderr = err
  try:
    b.refresh([job], force=True)
    b.live_raises = True
    b.refresh([job], force=True)
    b.refresh([job], force=True)
    b.live_raises = False
    b.refresh([job], force=True)
    b.refresh([job], force=True)
  finally:
    sys.stderr = real
  text = err.getvalue()
  assert text.count("squeue down") == 2, text
  assert text.count("did not finish") == 1, text
  b.terminal_times_out = False


def exercise_stamps_sticky_exited_and_pruning():
  from libtbx.jobs.backends.base import PollResult
  tmp = tempfile.mkdtemp()
  b = _backend()
  a, c = _submitted(b, tmp, 2)
  b.sched = {"1": ("running", "RUNNING", None), "2": ("running", "RUNNING", None)}
  b.refresh([a, c], force=True)
  # An older query finishing after a newer one must not win.
  later = time.monotonic() + 5
  with b._cache_lock:
    b._cache[a.handle] = PollResult("running", reason="newer")
    b._stamps[a.handle] = later
  b.sched["1"] = ("pending", "PENDING", None)
  b.refresh([a, c], force=True)
  assert b.poll(a).reason == "newer"
  # exited is sticky against running but yields to pending.
  _write_sentinel(c.handle, 0)
  b.sched["2"] = ("exited", "COMPLETED", None)
  b.refresh([a, c], force=True)
  assert b.poll(c).phase == "exited"
  b.sched["2"] = ("running", "RUNNING", None)
  b.refresh([a, c], force=True)
  assert b.poll(c).phase == "exited", "a running result overwrote exited"
  b.sched["2"] = ("pending", "PENDING", "Requeued")
  b.refresh([a, c], force=True)
  assert b.poll(c).phase == "pending"
  # Pruning: exited entries absent from a full refresh go, live ones stay.
  _write_sentinel(c.handle, 0)
  b.sched["2"] = ("exited", "COMPLETED", None)
  b.refresh([a, c], force=True)
  with b._cache_lock:
    b._stamps[a.handle] = 0.0
  b._last_refresh = None
  b.refresh([a])
  assert c.handle not in b._cache and c.handle not in b._submitted
  b._last_refresh = None
  b.refresh([c])                      # a refresh without `a`
  assert a.handle in b._cache, "a live entry was pruned"


def exercise_pruning_survives_interleaved_refreshes():
  tmp = tempfile.mkdtemp()
  b = _backend()
  a, c = _submitted(b, tmp, 2)
  _write_sentinel(a.handle, 0)
  b.sched = {"1": ("exited", "COMPLETED", None), "2": ("running", "RUNNING", None)}
  b.refresh([a])
  assert b.poll(a).phase == "exited"
  # A forced partial refresh (adopt) must not prune an unconsumed entry.
  b.refresh([c], force=True)
  assert a.handle in b._cache and a.handle in b._submitted
  assert b.poll(a).phase == "exited"
  # Nor may a full refresh that started before the entry was written.
  with b._cache_lock:
    b._stamps[a.handle] = time.monotonic() + 5
  b._last_refresh = None
  b.refresh([c])
  assert a.handle in b._cache, "an entry newer than the refresh was pruned"
  with b._cache_lock:
    b._stamps[a.handle] = 0.0
  b._last_refresh = None
  b.refresh([c])
  assert a.handle not in b._cache and a.handle not in b._submitted
  # The absent-handle grace does not undo a sticky exited entry.
  (e,) = _submitted(b, tmp, 1)
  _write_sentinel(e.handle, 0)
  b.sched = {"3": ("exited", "COMPLETED", None)}
  b.refresh([e], force=True)
  assert b.poll(e).phase == "exited"
  os.remove(os.path.join(e.handle.job_dir, "exit_code"))
  b.sched = {}
  b.refresh([e], force=True)
  r = b.poll(e)
  assert (r.phase, r.exit_code) == ("exited", 0), r


def _returns_within(fn, seconds):
  """Run ``fn`` on a thread; True when it returns within ``seconds``."""
  import threading
  done = threading.Event()
  threading.Thread(target=lambda: (fn(), done.set()), daemon=True).start()
  return done.wait(seconds)


def exercise_poll_does_not_wait_for_a_slow_query():
  import threading
  tmp = tempfile.mkdtemp()
  b = _backend()
  (job,) = _submitted(b, tmp, 1)
  b.sched = {"1": ("running", "RUNNING", None)}
  b.live_gate = threading.Event()
  b.live_entered = threading.Event()
  t = threading.Thread(target=b.refresh, args=([job],), kwargs={"force": True})
  t.start()
  assert b.live_entered.wait(5), "refresh never reached the scheduler query"
  try:
    assert _returns_within(lambda: b.poll(job), 2.0), \
      "poll waited behind a scheduler query"
  finally:
    b.live_gate.set()
    t.join()
  assert b.poll(job).phase == "running"


def exercise_adopt_finds_the_sentinel_in_the_original_job_dir():
  from libtbx.jobs.manager import JobManager
  tmp = tempfile.mkdtemp()
  b = _backend()
  (job,) = _submitted(b, tmp, 1)
  _write_sentinel(job.handle, 0)
  b.sched = {"1": ("exited", "COMPLETED", None)}
  manager = JobManager(b, reaper_interval=3600)
  try:
    adopted = manager.adopt(job.spec, {"scheduler_job_id": "1",
                                       "job_dir": job.handle.job_dir})
    assert adopted.state == "finished" and adopted.exit_code == 0, adopted.state
    assert adopted.job_id != job.job_id
  finally:
    manager.stop()


def run():
  if sys.platform == "win32":
    print("Skipping tst_batch on Windows")
    print("OK")
    return
  exercise_brace_template_leaves_shell_syntax_alone()
  exercise_render_defaults_resources_and_takeover()
  exercise_render_errors()
  exercise_job_name_is_sanitised()
  exercise_wrapper_runs_the_command_and_writes_the_sentinel()
  exercise_submit_failures_become_submit_errors()
  exercise_terminate_shutdown_handle_info_and_reattach()
  exercise_refresh_rate_limit_empty_list_and_live_states()
  exercise_sentinel_is_read_only_after_the_scheduler_reports_terminal()
  exercise_accounting_resolves_and_stale_live_rows_do_not()
  exercise_accounting_code_yields_to_a_sentinel_within_the_grace()
  exercise_a_requeue_removes_the_abandoned_attempts_sentinel()
  exercise_late_sentinel_wins_over_the_fallback()
  exercise_accounting_timeout_is_transient()
  exercise_handles_on_other_clusters_are_queried_there()
  exercise_one_failing_cluster_does_not_freeze_the_others()
  exercise_accounting_codes_only_count_for_completed_and_failed()
  exercise_grace_periods_and_fallbacks()
  exercise_query_failure_keeps_the_cache()
  exercise_stamps_sticky_exited_and_pruning()
  exercise_pruning_survives_interleaved_refreshes()
  exercise_poll_does_not_wait_for_a_slow_query()
  exercise_adopt_finds_the_sentinel_in_the_original_job_dir()
  print("OK")


if __name__ == "__main__":
  run()
