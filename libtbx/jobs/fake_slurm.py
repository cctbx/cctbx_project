"""Fake ``sbatch``, ``squeue``, ``sacct`` and ``scancel`` for tests.

Run as ``python fake_slurm.py <tool> <args...>``; ``install_shims`` writes
one shell wrapper per tool onto ``PATH``. State lives in the JSON file
named by ``FAKE_SLURM_STATE``. Knobs, all environment variables:

- ``FAKE_SLURM_RUN=1``: ``sbatch`` runs the script in the background;
  ``squeue`` reports RUNNING, then the terminal state; ``scancel`` kills.
- ``FAKE_SLURM_SBATCH_FAIL=<msg>``: ``sbatch`` prints ``<msg>``, exits 1.
- ``FAKE_SLURM_TEST_ONLY_FAIL=<msg>``: ``sbatch --test-only`` fails.
- ``FAKE_SLURM_SACCT_MISSING=1``: ``sacct`` behaves as if not installed.
- ``FAKE_SLURM_EXTRA_ROWS=<rows>``: extra ``squeue`` lines, ``;``-separated.
- ``FAKE_SLURM_CLUSTER_SUFFIX=<name>``: ``sbatch`` routes to cluster
  ``<name>`` and prints ``id;<name>``; later tools see the job only with
  ``--clusters=<name>``.

The shims exit 2 when a flag the backend relies on is missing
(``--parsable``, ``--noheader``, ``--states=all``, ``--parsable2``,
``--open-mode=append``, the exact ``--format``), so a command-line
regression fails the tests instead of passing by accident;
``squeue --user`` is optional, as in the backend. A ``purged`` job is
left out of ``squeue`` but answered by ``sacct``, like one older than
``MinJobAge``.

Standard library only: the wrappers run it with the test's interpreter.
"""

import json
import os
import subprocess
import sys
import time

TOOLS = ("sbatch", "squeue", "sacct", "scancel")
TERMINAL = ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL",
            "OUT_OF_MEMORY", "PREEMPTED")


class _Lock(object):

  def __init__(self, path):
    self.path = path + ".lock"

  def __enter__(self):
    deadline = time.time() + 10
    while True:
      try:
        os.close(os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return self
      except FileExistsError:
        if time.time() > deadline:
          raise RuntimeError("fake_slurm: stale lock %s" % self.path)
        time.sleep(0.01)

  def __exit__(self, *exc):
    try:
      os.remove(self.path)
    except OSError:
      pass


def _load(path):
  if not os.path.exists(path):
    return {"next_id": 100, "jobs": {}}
  with open(path, encoding="utf-8") as fh:
    return json.load(fh)


def _save(path, state):
  tmp = path + ".tmp"
  with open(tmp, "w", encoding="utf-8") as fh:
    json.dump(state, fh)
  os.replace(tmp, path)


def _require(args, tool, *flags):
  """Exit 2 unless every flag (``--name`` or ``--name=value``) is present."""
  for flag in flags:
    name, _, value = flag.partition("=")
    if value:
      present = flag in args
    else:
      present = any(a == name or a.startswith(name + "=") for a in args)
    if not present:
      print("fake_slurm: %s called without %s" % (tool, flag), file=sys.stderr)
      sys.exit(2)


def _cluster(args):
  """The ``--clusters=<name>`` / ``-M <name>`` value, or ``None``."""
  for i, arg in enumerate(args):
    if arg.startswith("--clusters="):
      return arg[len("--clusters="):]
    if arg in ("-M", "--clusters") and i + 1 < len(args):
      return args[i + 1]
  return None


def _visible(job, cluster):
  """Whether a query addressed to ``cluster`` sees ``job``."""
  return job.get("cluster") == cluster


def _alive(pid):
  try:
    os.kill(pid, 0)
  except ProcessLookupError:
    return False
  except PermissionError:
    return True
  return True


def _settle(state):
  """RUN mode: fold finished background scripts into terminal states.

  An rc file means completed or failed; a process gone without one was
  killed (NODE_FAIL). ``scancel`` marks its own jobs CANCELLED directly.
  """
  for job in state["jobs"].values():
    if job["state"] != "RUNNING" or not job.get("rc_file"):
      continue
    rc_file = job["rc_file"]
    if os.path.exists(rc_file):
      with open(rc_file, encoding="utf-8") as fh:
        text = fh.read().strip()
      rc = int(text) if text.isdigit() else 1
      job["exit_code"] = rc
      if job.get("cancel_requested"):
        job["state"] = "CANCELLED"
      else:
        job["state"] = "COMPLETED" if rc == 0 else "FAILED"
      job["reason"] = None
    elif job.get("pid") and not _alive(job["pid"]):
      job["exit_code"] = None
      job["state"] = "CANCELLED" if job.get("cancel_requested") else "NODE_FAIL"
      job["reason"] = None


def set_job(state_path, job_id, **fields):
  """Test helper: overwrite fields of one job in the state file."""
  with _Lock(state_path):
    state = _load(state_path)
    state["jobs"].setdefault(str(job_id), {
      "state": "PENDING", "reason": "Priority", "exit_code": None,
      "pid": None, "rc_file": None, "cancel_requested": False,
      "cluster": None, "purged": False})
    state["jobs"][str(job_id)].update(fields)
    _save(state_path, state)


def install_shims(bindir, state_path):
  """Write the four wrappers into ``bindir`` and put it first on PATH."""
  os.makedirs(bindir, exist_ok=True)
  here = os.path.abspath(__file__)
  for tool in TOOLS:
    path = os.path.join(bindir, tool)
    with open(path, "w", encoding="utf-8") as fh:
      fh.write('#!/bin/sh\nexec "%s" "%s" %s "$@"\n'
               % (sys.executable, here, tool))
    os.chmod(path, 0o700)
  os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")
  os.environ["FAKE_SLURM_STATE"] = state_path


def _sbatch(args, state_path):
  _require(args, "sbatch", "--parsable", "--chdir", "--output", "--error",
           "--open-mode=append", "--job-name")
  if "--test-only" in args:
    msg = os.environ.get("FAKE_SLURM_TEST_ONLY_FAIL")
    if msg:
      print("sbatch: error: %s" % msg, file=sys.stderr)
      return 1
    print("sbatch: Job 1 to start at 2026-09-24T00:00:00", file=sys.stderr)
    return 0
  msg = os.environ.get("FAKE_SLURM_SBATCH_FAIL")
  if msg:
    print("sbatch: error: %s" % msg, file=sys.stderr)
    return 1
  opts, script = {}, None
  for arg in args:
    if arg.startswith("--"):
      key, _, value = arg[2:].partition("=")
      opts[key] = value or True
    else:
      script = arg
  with _Lock(state_path):
    state = _load(state_path)
    job_id = str(state["next_id"])
    state["next_id"] += 1
    cluster = os.environ.get("FAKE_SLURM_CLUSTER_SUFFIX") or None
    job = {"state": "PENDING", "reason": "Priority", "exit_code": None,
           "pid": None, "rc_file": None, "cancel_requested": False,
           "script": script, "name": opts.get("job-name"),
           "output": opts.get("output"), "error": opts.get("error"),
           "open_mode": opts.get("open-mode"), "chdir": opts.get("chdir"),
           "requeue": (True if "requeue" in opts else
                       False if "no-requeue" in opts else None),
           "cluster": cluster, "purged": False}
    if os.environ.get("FAKE_SLURM_RUN") and script:
      rc_file = script + ".rc"
      output = opts.get("output") or os.devnull
      proc = subprocess.Popen(
        ["bash", "-c", 'bash "$1" > "$2" 2>&1; echo $? > "$3"', "_",
         script, output, rc_file],
        cwd=opts.get("chdir") or None,
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, start_new_session=True)
      job.update(state="RUNNING", reason=None, pid=proc.pid, rc_file=rc_file)
    state["jobs"][job_id] = job
    _save(state_path, state)
  print(job_id + (";%s" % cluster if cluster else ""))
  return 0


def _squeue(args, state_path):
  # --user is optional: the backend omits it for a uid with no passwd entry.
  _require(args, "squeue", "--noheader", "--states=all", "--format=%i|%T|%r")
  cluster = _cluster(args)
  with _Lock(state_path):
    state = _load(state_path)
    _settle(state)
    _save(state_path, state)
  if cluster:
    print("CLUSTER: %s" % cluster)
  for job_id, job in state["jobs"].items():
    if job.get("purged") or not _visible(job, cluster):
      continue
    print("%s|%s|%s" % (job_id, job["state"], job.get("reason") or "None"))
  extra = os.environ.get("FAKE_SLURM_EXTRA_ROWS")
  if extra:
    for row in extra.split(";"):
      print(row)
  return 0


def _sacct(args, state_path):
  if os.environ.get("FAKE_SLURM_SACCT_MISSING"):
    print("sacct: command not found", file=sys.stderr)
    return 127
  _require(args, "sacct", "--noheader", "--parsable2", "--jobs",
           "--format=JobID,State,ExitCode")
  cluster = _cluster(args)
  ids = []
  for arg in args:
    if arg.startswith("--jobs="):
      ids = arg[len("--jobs="):].split(",")
  with _Lock(state_path):
    state = _load(state_path)
    _settle(state)
    _save(state_path, state)
  for job_id in ids:
    job = state["jobs"].get(job_id)
    if job is None or not _visible(job, cluster):
      continue
    if job["state"] not in TERMINAL:
      print("%s|%s|0:0" % (job_id, job["state"]))
      continue
    if job["state"] == "CANCELLED":
      print("%s|CANCELLED by 1234|0:15" % job_id)
      print("%s.batch|CANCELLED|0:15" % job_id)
    else:
      # As in sacct, a scheduler-killed job without its own code shows 0:0.
      code = job.get("exit_code") or 0
      print("%s|%s|%d:0" % (job_id, job["state"], code))
      print("%s.batch|%s|%d:0" % (job_id, job["state"], code))
  return 0


def _scancel(args, state_path):
  job_id = args[-1]
  cluster = _cluster(args)
  for arg in args[:-1]:
    if arg.startswith("-") and not arg.startswith("--clusters=") \
       and arg not in ("-M", "--clusters"):
      print("scancel: error: unknown option %s" % arg, file=sys.stderr)
      return 2
  if not job_id.isdigit():
    print("scancel: error: Invalid job id %s" % job_id, file=sys.stderr)
    return 1
  with _Lock(state_path):
    state = _load(state_path)
    job = state["jobs"].get(job_id)
    if job is None or not _visible(job, cluster):
      print("scancel: error: Invalid job id %s" % job_id, file=sys.stderr)
      return 1
    if job["state"] in TERMINAL:
      print("scancel: error: Kill job error on job id %s: Job/step already "
            "completing or completed" % job_id, file=sys.stderr)
      return 1
    if job.get("pid"):
      import signal
      job["cancel_requested"] = True
      try:
        os.killpg(os.getpgid(job["pid"]), signal.SIGTERM)
      except OSError:
        pass
      # Terminal at once, like a real controller; an unreaped zombie would
      # keep _settle from seeing the process gone.
      job["state"] = "CANCELLED"
      job["reason"] = None
      job["exit_code"] = None
    else:
      job["state"] = "CANCELLED"
      job["reason"] = None
    _save(state_path, state)
  return 0


def main(argv):
  tool, args = argv[0], argv[1:]
  state_path = os.environ["FAKE_SLURM_STATE"]
  return {"sbatch": _sbatch, "squeue": _squeue, "sacct": _sacct,
          "scancel": _scancel}[tool](args, state_path)


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
