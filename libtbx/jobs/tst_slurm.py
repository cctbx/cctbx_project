"""SlurmBackend tests against the fake command-line shims."""

import os
import sys
import tempfile
import time


_KNOBS = ("FAKE_SLURM_RUN", "FAKE_SLURM_SBATCH_FAIL", "FAKE_SLURM_TEST_ONLY_FAIL",
          "FAKE_SLURM_SACCT_MISSING", "FAKE_SLURM_EXTRA_ROWS",
          "FAKE_SLURM_CLUSTER_SUFFIX")


def _fresh(**knobs):
  """Install shims with a fresh state file; set the given knobs only."""
  from libtbx.jobs.fake_slurm import install_shims
  tmp = tempfile.mkdtemp()
  state = os.path.join(tmp, "state.json")
  install_shims(os.path.join(tmp, "bin"), state)
  for key in _KNOBS:
    os.environ.pop(key, None)
  os.environ.update(knobs)
  return tmp, state


def _spec(tmp, **kw):
  from libtbx.jobs.manager import JobSpec
  kw.setdefault("argv", [sys.executable, "-c", "print('hi')"])
  kw.setdefault("cwd", tmp)
  kw.setdefault("log_path", os.path.join(tmp, "job.log"))
  return JobSpec(**kw)


def _job(spec, job_id="j_000000000001"):
  from libtbx.jobs.manager import Job
  return Job(job_id=job_id, spec=spec)


def _backend(**kw):
  from libtbx.jobs.backends.slurm import SlurmBackend
  kw.setdefault("poll_interval", 1.0)
  return SlurmBackend(**kw)


def _await(manager, job, states, timeout=60):
  deadline = time.time() + timeout
  while time.time() < deadline:
    manager.reaper_iter()
    if job.state in states:
      return
    time.sleep(0.2)
  raise AssertionError("job stuck in %s (%s)" % (job.state, job.reason))


def exercise_submit_command_and_id_parsing():
  import json
  tmp, state = _fresh()
  b = _backend()
  job = _job(_spec(tmp, name="phenix refine", resources={"cpus": 2, "mem_mb": 4000,
                                                       "time_minutes": 61, "gpus": 1}))
  phase, handle = b.submit(job)
  assert phase == "pending" and handle.scheduler_job_id == "100", handle
  assert handle.cluster is None, handle
  with open(state, encoding="utf-8") as fh:
    recorded = json.load(fh)["jobs"]["100"]
  assert recorded["name"] == "phenix_refine", recorded
  assert recorded["output"] == os.path.join(handle.job_dir, "scheduler.log")
  assert recorded["error"] == recorded["output"], recorded
  assert recorded["open_mode"] == "append", recorded
  assert recorded["requeue"] is None, "the scheduler's own policy by default"
  assert recorded["chdir"] == tmp
  for wanted, flag in ((False, "--no-requeue"), (True, "--requeue")):
    cmd = _backend(requeue=wanted).submit_command(job, "/x/job.sbatch")
    assert flag in cmd and cmd[-1] == "/x/job.sbatch", cmd
  with open(os.path.join(handle.job_dir, "run.sh"), encoding="utf-8") as fh:
    wrapper = fh.read()
  assert 'libtbx_jobs_attempt="${SLURM_RESTART_COUNT:-0}"' in wrapper, wrapper
  assert 'export LIBTBX_JOBS_ATTEMPT="$libtbx_jobs_attempt"' in wrapper, wrapper
  with open(os.path.join(handle.job_dir, "job.sbatch"), encoding="utf-8") as fh:
    script = fh.read()
  for line in ("#SBATCH --cpus-per-task=2", "#SBATCH --mem=4000M",
               "#SBATCH --time=1:01:00", "#SBATCH --gpus=1"):
    assert line in script, (line, script)
  assert script.startswith("#!/bin/bash"), script
  assert b.parse_submit_output("4242\n", "") == ("4242", None)
  assert b.parse_submit_output("4242;gpu\n", "") == ("4242", "gpu")
  assert b.parse_submit_output("\n4242\nsite banner\n", "") == ("4242", None)
  from libtbx.jobs.backends.base import SubmitError
  for garbage in ("garbage", "4242;", "4242 extra", "4242;a b"):
    try:
      b.parse_submit_output(garbage, "")
    except SubmitError:
      pass
    else:
      raise AssertionError("sbatch output %r accepted" % garbage)


def exercise_a_routed_job_is_addressed_to_its_cluster():
  from libtbx.jobs.manager import JobManager
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh(FAKE_SLURM_CLUSTER_SUFFIX="gpu")
  b = _backend(poll_interval=1.0)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp))
    assert job.handle.cluster == "gpu", job.handle
    d = manager.get_status(job.job_id)
    assert d["scheduler_cluster"] == "gpu", d
    # The shims show a routed job only to commands carrying --clusters.
    assert b.query_live([job.handle.scheduler_job_id]) == {}
    set_job(state, job.handle.scheduler_job_id, state="RUNNING", reason=None)
    manager.reaper_iter()
    assert job.state == "running", (job.state, job.reason)
    manager.cancel(job.job_id)
    time.sleep(1.1)
    manager.reaper_iter()
    assert job.state == "cancelled", (job.state, job.reason)
    assert b.query_terminal(["100"], cluster="gpu")["100"][1] == "CANCELLED"
    assert b.query_terminal(["100"]) == {}
  finally:
    manager.stop()


def exercise_a_job_purged_from_squeue_resolves_through_sacct():
  from libtbx.jobs.manager import JobManager
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh()
  b = _backend(poll_interval=1.0, sentinel_grace=0.2)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp))
    set_job(state, job.handle.scheduler_job_id, state="RUNNING", reason=None)
    manager.reaper_iter()
    assert job.state == "running", (job.state, job.reason)
    # Older than MinJobAge: gone from squeue, still in accounting.
    set_job(state, job.handle.scheduler_job_id, state="FAILED", exit_code=3,
            purged=True)
    assert b.query_live([job.handle.scheduler_job_id]) == {}
    with b._cache_lock:
      b._submitted[job.handle] -= 100
    time.sleep(1.1)
    manager.reaper_iter()
    assert job.state == "running" and job.reason.startswith("finishing"), \
      (job.state, job.reason)
    time.sleep(1.1)
    manager.reaper_iter()
    assert job.state == "failed" and job.exit_code == 3, (job.state, job.exit_code, job.reason)
    assert job.reason == "FAILED", job.reason
  finally:
    manager.stop()


def exercise_default_template_renders_each_directive_once():
  tmp = tempfile.mkdtemp()
  b = _backend()
  text = b.render(_spec(tmp, resources={"cpus": 2, "mem_mb": 4000,
                                        "time_minutes": 61, "gpus": 1}))
  lines = text.splitlines()
  for directive in ("#SBATCH --cpus-per-task=2", "#SBATCH --mem=4000M",
                    "#SBATCH --time=1:01:00", "#SBATCH --gpus=1"):
    assert lines.count(directive) == 1, (directive, text)
  assert not [l for l in lines if l.startswith("# #SBATCH")], text
  assert "{{" not in text, text
  assert sum(1 for l in lines if l.startswith("#SBATCH")) == 4, text


def exercise_accounting_zero_for_a_killed_job_is_not_success():
  from libtbx.jobs.manager import JobManager
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh()
  b = _backend(poll_interval=1.0)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp))
    set_job(state, job.handle.scheduler_job_id, state="RUNNING", reason=None)
    manager.reaper_iter()
    assert job.state == "running", (job.state, job.reason)
    set_job(state, job.handle.scheduler_job_id, state="TIMEOUT",
            exit_code=None)
    assert b.query_terminal([job.handle.scheduler_job_id])[
      job.handle.scheduler_job_id] == (0, "TIMEOUT", None)
    time.sleep(1.1)
    manager.reaper_iter()
    assert job.state == "failed", (job.state, job.exit_code, job.reason)
    assert job.exit_code is None, job.exit_code
    assert job.reason.startswith("TIMEOUT"), job.reason
  finally:
    manager.stop()


def exercise_live_states_reasons_and_filtering():
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh(FAKE_SLURM_EXTRA_ROWS="999|RUNNING|None;1000|SUSPENDED|None")
  b = _backend()
  assert b.query_live(["1", "2"]) == {}, "empty queue must give an empty map"
  set_job(state, "1", state="PENDING", reason="Resources")
  set_job(state, "2", state="RUNNING", reason=None)
  set_job(state, "3", state="COMPLETING", reason="(null)")
  set_job(state, "4", state="WEIRD_NEW_STATE", reason=None)
  set_job(state, "5", state="TIMEOUT", reason=None)
  live = b.query_live(["1", "2", "3", "4", "5"])
  assert set(live) == {"1", "2", "3", "4", "5"}, live
  assert live["1"] == ("pending", "PENDING", "Resources"), live["1"]
  assert live["2"] == ("running", "RUNNING", None)
  assert live["3"] == ("running", "COMPLETING", None), live["3"]
  assert live["4"] == ("running", "WEIRD_NEW_STATE", None), live["4"]
  assert live["5"] == ("exited", "TIMEOUT", None)
  set_job(state, "6", state="PENDING",
          reason="ReqNodeNotAvail, UnavailableNodes:n[1-3]")
  set_job(state, "7", state="PENDING", reason="(null)")
  live = b.query_live(["6", "7"])
  assert live["6"] == ("pending", "PENDING",
                       "ReqNodeNotAvail, UnavailableNodes:n[1-3]"), live["6"]
  assert live["7"] == ("pending", "PENDING", None), live["7"]
  assert "999" not in live and "1000" not in live
  assert b.state_phase("CANCELLED") == "exited" and b.is_cancelled_state("CANCELLED")
  assert b.is_completed_state("COMPLETED") and b.is_failed_state("FAILED")
  assert b.state_phase("REQUEUE_HOLD") == "pending"
  assert b.state_phase("STOPPED") == "running"
  # squeue asks about the process's own user, not LOGNAME; a uid with no
  # passwd entry drops --user and relies on the id filter.
  import pwd
  try:
    expected = pwd.getpwuid(os.getuid()).pw_name
  except KeyError:
    expected = None
  saved = os.environ.get("LOGNAME")
  os.environ["LOGNAME"] = "somebody_else"
  real_getpwuid = pwd.getpwuid
  try:
    assert b._submitting_user() == expected
    if expected is not None:
      assert "--user=%s" % expected in b._squeue_command(), b._squeue_command()
      assert "--user=somebody_else" not in b._squeue_command()
    def no_entry(uid):
      raise KeyError(uid)
    pwd.getpwuid = no_entry
    assert b._submitting_user() is None
    cmd = b._squeue_command("gpu")
    assert not [a for a in cmd if a.startswith("--user")], cmd
    assert cmd[-2:] == ["--format=%i|%T|%r", "--clusters=gpu"], cmd
    # The no-user command works against the shims end to end.
    assert b.query_live(["1", "2"])["2"] == ("running", "RUNNING", None)
  finally:
    pwd.getpwuid = real_getpwuid
    if saved is None:
      os.environ.pop("LOGNAME", None)
    else:
      os.environ["LOGNAME"] = saved


def exercise_terminal_rows_from_sacct():
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh()
  b = _backend()
  set_job(state, "7", state="CANCELLED", exit_code=None)
  set_job(state, "8", state="FAILED", exit_code=3)
  set_job(state, "9", state="RUNNING")
  rows = b.query_terminal(["7", "8", "9"])
  assert rows["7"] == (0, "CANCELLED", "signal 15"), rows["7"]
  assert rows["8"] == (3, "FAILED", None), rows["8"]
  assert rows["9"] == (0, "RUNNING", None), rows["9"]
  assert "7.batch" not in rows
  os.environ["FAKE_SLURM_SACCT_MISSING"] = "1"
  try:
    b.query_terminal(["7"])
  except RuntimeError as exc:
    assert "sacct" in str(exc), str(exc)
  else:
    raise AssertionError("missing sacct did not raise")


def exercise_refusal_check_and_failure_hint():
  from libtbx.jobs.manager import JobManager
  from libtbx.jobs.backends.base import SubmitError
  tmp, state = _fresh(FAKE_SLURM_SBATCH_FAIL="Invalid partition name specified")
  b = _backend()
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp))
    assert job.state == "failed" and "Invalid partition" in job.reason, (job.state, job.reason)
  finally:
    manager.stop()
  os.environ.pop("FAKE_SLURM_SBATCH_FAIL")
  b.check(_spec(tmp))
  # The check carries the submission's own options.
  seen = []
  def record(cmd, **kw):
    seen.append(list(cmd))
    return type(b).run_command(checker, cmd, **kw)
  checker = _backend(requeue=False)
  checker.run_command = record
  checker.check(_spec(tmp))
  assert seen and seen[-1][0] == "sbatch" and seen[-1][-2] == "--test-only", seen
  assert "--no-requeue" in seen[-1] and "--open-mode=append" in seen[-1], seen[-1]
  os.environ["FAKE_SLURM_TEST_ONLY_FAIL"] = "Requested time limit is invalid"
  try:
    b.check(_spec(tmp))
  except SubmitError as exc:
    assert "time limit" in str(exc), str(exc)
  else:
    raise AssertionError("--test-only failure did not raise")
  limit = ("x\nslurmstepd: error: *** JOB 5 ON n1 CANCELLED AT "
           "2026-09-24T10:00:00 DUE TO TIME LIMIT ***\n")
  hint = b.failure_hint(limit, "TIMEOUT")
  assert hint and "DUE TO TIME LIMIT" in hint, hint
  assert b.failure_hint("nothing here", "TIMEOUT") is None
  # A line is a hint only for the state it explains: an abandoned
  # attempt's time-limit line says nothing about a final NODE_FAIL.
  assert b.failure_hint(limit, "NODE_FAIL") is None
  assert b.failure_hint("slurmstepd: error: Detected 1 oom-kill event(s)\n",
                        "OUT_OF_MEMORY")
  assert b.failure_hint("slurmstepd: error: *** JOB 5 CANCELLED AT x DUE TO "
                        "PREEMPTION ***\n", "PREEMPTED")
  # A requeue notice explains an abandoned attempt, never the final state.
  assert b.failure_hint("slurmstepd: error: *** JOB 5 ON n1 CANCELLED AT "
                        "2026-09-24T10:00:00 DUE TO JOB REQUEUE ***\n",
                        "CANCELLED") is None
  for state in ("REQUEUED", "REQUEUE_HOLD", "REQUEUE_FED", "SPECIAL_EXIT"):
    assert b.is_requeue_state(state) and b.state_phase(state) == "pending", state
  assert not b.is_requeue_state("PENDING") and not b.is_requeue_state("RUNNING")


def exercise_end_to_end_run_cancel_and_racing_completion():
  from libtbx.jobs.manager import JobManager
  from libtbx.jobs.fake_slurm import set_job
  tmp, state = _fresh(FAKE_SLURM_RUN="1")
  b = _backend(poll_interval=1.0)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp, argv=[sys.executable, "-c",
                                          "print('from slurm'); import sys; sys.exit(0)"]))
    assert job.state == "queued", job.state
    d = manager.get_status(job.job_id)
    assert d["scheduler"] == "slurm" and d["scheduler_job_id"] == "100"
    assert d["job_dir"] == job.handle.job_dir
    _await(manager, job, ("finished", "failed", "cancelled"))
    assert job.state == "finished" and job.exit_code == 0, (job.state, job.reason)
    with open(os.path.join(tmp, "job.log"), encoding="utf-8") as fh:
      assert "from slurm" in fh.read()
    assert b._read_sentinel(job.handle) == 0
    long = manager.submit(_spec(tmp, argv=[sys.executable, "-c", "import time; time.sleep(60)"],
                                log_path=os.path.join(tmp, "long.log")))
    _await(manager, long, ("running",))
    manager.cancel(long.job_id)
    assert long.state == "cancelling"
    _await(manager, long, ("finished", "failed", "cancelled"))
    assert long.state == "cancelled", (long.state, long.reason)
  finally:
    manager.stop()
  # scancel racing a job that completed: the exit code wins.
  os.environ.pop("FAKE_SLURM_RUN")
  b = _backend(poll_interval=1.0)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp, log_path=os.path.join(tmp, "race.log")))
    set_job(state, job.handle.scheduler_job_id, state="RUNNING", reason=None)
    manager.reaper_iter()
    assert job.state == "running"
    set_job(state, job.handle.scheduler_job_id, state="COMPLETED", exit_code=0)
    with open(os.path.join(job.handle.job_dir, "exit_code"), "w",
              encoding="utf-8") as fh:
      fh.write("0\n")
    d = manager.cancel(job.job_id)
    assert d["state"] == "cancelling", d
    time.sleep(1.1)
    manager.reaper_iter()
    assert job.state == "finished" and job.exit_code == 0, (job.state, job.reason)
  finally:
    manager.stop()


def exercise_real_cluster_if_requested(original_env):
  """Run a trivial job on a real cluster when LIBTBX_JOBS_SLURM_TEST is set.

  The shims the other exercises installed are removed from PATH first,
  so this talks to the site's own sbatch, squeue, sacct and scancel.
  """
  if not os.environ.get("LIBTBX_JOBS_SLURM_TEST"):
    return
  for key in _KNOBS + ("FAKE_SLURM_STATE",):
    os.environ.pop(key, None)
  os.environ["PATH"] = original_env["PATH"]
  import shutil
  for tool in ("sbatch", "squeue", "sacct", "scancel"):
    found = shutil.which(tool)
    assert found, "%s is not on PATH" % tool
    with open(found, "rb") as fh:
      head = fh.read(512)
    assert b"fake_slurm" not in head, \
      "%s resolves to the fake shim %r" % (tool, found)
  from libtbx.jobs.manager import JobManager
  tmp = tempfile.mkdtemp(dir=os.getcwd())
  b = _backend(poll_interval=5.0)
  manager = JobManager(b, reaper_interval=3600)
  try:
    job = manager.submit(_spec(tmp))
    _await(manager, job, ("finished", "failed", "cancelled"), timeout=600)
    assert job.state == "finished", (job.state, job.reason)
  finally:
    manager.stop()


def run():
  if sys.platform == "win32":
    print("Skipping tst_slurm on Windows")
    print("OK")
    return
  original_env = {"PATH": os.environ.get("PATH", "")}
  exercise_submit_command_and_id_parsing()
  exercise_a_routed_job_is_addressed_to_its_cluster()
  exercise_a_job_purged_from_squeue_resolves_through_sacct()
  exercise_default_template_renders_each_directive_once()
  exercise_live_states_reasons_and_filtering()
  exercise_accounting_zero_for_a_killed_job_is_not_success()
  exercise_terminal_rows_from_sacct()
  exercise_refusal_check_and_failure_hint()
  exercise_end_to_end_run_cancel_and_racing_completion()
  exercise_real_cluster_if_requested(original_env)
  print("OK")


if __name__ == "__main__":
  run()
