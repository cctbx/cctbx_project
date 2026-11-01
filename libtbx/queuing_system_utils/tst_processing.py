"""
Regression test for queue submission through libtbx.queuing_system_utils.

Covers the Phenix GUI path (processing.qsub) and the older generic module
with fake queue commands, so no queuing system is needed:

- the job id and polling output are text, not bytes, and bytes that are
  not valid UTF-8 (a colleague's accented job name in qstat output) do
  not abort polling
- the job script sources the submitting build's setpaths.sh, or
  $PHENIX/phenix_env.sh where there is none (Phenix 2.x installations
  ship no setpaths.sh), and a missing script is reported instead of
  being sourced on the cluster
- the generated job script runs under /bin/sh and dash with a minimal
  environment, so the environment line has to work
"""
import gc
import glob
import os
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import warnings

from libtbx.utils import Sorry

# Fake queue commands. Submission commands save the job script they receive
# on stdin to submitted_script.txt and print what the real command prints.
# Some output carries a Latin-1 byte (\351) that is not valid UTF-8.
FAKE_COMMANDS = {
  "qsub": """\
#!/bin/sh
cat > submitted_script.txt
terse=0
cwd=0
for a in "$@"; do
  [ "$a" = "-terse" ] && terse=1
  [ "$a" = "-cwd" ] && cwd=1
done
if [ $cwd = 1 ] && [ $terse = 0 ]; then
  printf 'Your job 12345 ("Jos\\351") has been submitted\\n'
else
  echo 12345
fi
""",
  "qsub_slow": """\
#!/bin/sh
cat > submitted_script.txt
echo 'Your job 12345 ("test") has been submitted'
sleep 2
touch qsub_finished
""",
  "sbatch": """\
#!/bin/sh
cat > submitted_script.txt
echo "Submitted batch job 67890"
""",
  "bsub": """\
#!/bin/sh
cat > submitted_script.txt
echo "Job <4242> is submitted to default queue <normal>."
""",
  "condor_submit": """\
#!/bin/sh
cat > /dev/null
echo "1 job(s) submitted to cluster 777."
""",
  "qstat": """\
#!/bin/sh
case "$1" in
  -j) printf 'job_number: %s\\njob_name: Jos\\351\\n' "$2" ;;
  -f) echo "job_state = R" ;;
  *)
    echo "Job id            Name             User     Time Use S Queue"
    echo "----------------  ---------------- -------- -------- - -----"
    printf '12345             Jos\\351             user     00:00:00 R workq\\n'
    ;;
esac
""",
  "squeue": """\
#!/bin/sh
echo R
""",
  "bjobs": """\
#!/bin/sh
echo "4242 user RUN normal host host test"
""",
  "condor_q": """\
#!/bin/sh
echo '<?xml version="1.0"?>'
echo '<classads><c><a n="ClusterId"><i>777</i></a></c></classads>'
""",
}

# platform, submit command, polling command with arguments, expected job id
PLATFORMS = [
  ("sge", "qsub", ("qstat", "-j"), "12345"),
  ("pbs", "qsub", ("qstat", "-f"), "12345"),
  ("pbspro", "qsub", ("qstat",), "12345"),
  ("lsf", "bsub", ("bjobs",), "4242"),
  ("slurm", "sbatch", ("squeue", "-o", "%t", "-h", "-j"), "67890"),
  ("condor", "condor_submit", ("condor_q", "-xml"), "777"),
]

class write_marker(object):
  """Picklable job target that writes a file so the test can see it ran."""

  def __init__(self, path):
    self.path = path

  def __call__(self):
    with open(self.path, "w") as f:
      f.write("ran\n")
    print("marker written")

def get_target(path):
  """
  Return a write_marker whose class pickles by its importable module name,
  so a separate libtbx.python process can unpickle it.
  """
  from libtbx.queuing_system_utils import tst_processing
  return tst_processing.write_marker(path)

def dispatcher_dir():
  """Directory holding the libtbx.python dispatcher, or None."""
  import libtbx.load_env
  candidates = [libtbx.env.under_build("bin")]
  candidates.extend(os.environ.get("PATH", "").split(os.pathsep))
  for d in candidates:
    if os.path.isfile(os.path.join(d, "libtbx.python")):
      return d
  return None

class sandbox(object):
  """Temporary working directory containing the fake queue commands."""

  def __init__(self):
    # a space and a dollar in the path: unquoted or double-quoted paths in
    # the job script would break
    self.dir = tempfile.mkdtemp(prefix="tst_processing $x ")
    self.old_cwd = os.getcwd()
    self.old_env = dict(
      (k, os.environ.get(k)) for k in ("PHENIX", "PHENIX_CUSTOM_ENV"))
    for name, text in FAKE_COMMANDS.items():
      path = os.path.join(self.dir, name)
      with open(path, "w") as f:
        f.write(text)
      os.chmod(path, stat.S_IRWXU)
    os.chdir(self.dir)

  def command(self, name, *args):
    """Command line for a fake command, as accepted by processing.qsub."""
    return " ".join(['"%s"' % os.path.join(self.dir, name)] + list(args))

  def phenix_env(self, bin_dir=None):
    """Write a fake phenix_env.sh here and point PHENIX at this directory."""
    path = os.path.join(self.dir, "phenix_env.sh")
    with open(path, "w") as f:
      f.write("# fake phenix_env.sh\n")
      if bin_dir is not None:
        f.write('export PATH="%s:$PATH"\n' % bin_dir)
    os.environ["PHENIX"] = self.dir
    os.environ.pop("PHENIX_CUSTOM_ENV", None)
    return path

  def submitted_script(self):
    with open(os.path.join(self.dir, "submitted_script.txt")) as f:
      return f.read()

  def run_submitted_script(self, shell):
    """
    Run the saved job script under `shell` with a minimal environment, the
    way a queue without environment export would, returning its output.
    """
    env = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", self.dir)}
    # Popen rather than easy_run: an argv list, no shell, and this exact
    # environment are the point of the check
    p = subprocess.Popen(
      [shell, os.path.join(self.dir, "submitted_script.txt")],
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env,
      encoding="utf-8", errors="replace")
    out, err = p.communicate()
    return p.returncode, out, err

  def cleanup(self):
    os.chdir(self.old_cwd)
    for k, v in self.old_env.items():
      if v is None:
        os.environ.pop(k, None)
      else:
        os.environ[k] = v
    shutil.rmtree(self.dir, ignore_errors=True)

class without_setpaths(object):
  """Hide the build's setpaths.sh, as in a Phenix 2.x installation."""

  def __init__(self, directory):
    self.missing = os.path.join(directory, "no_such_dir", "setpaths.sh")

  def __enter__(self):
    import libtbx.load_env
    original = libtbx.env.under_build
    libtbx.env.under_build = lambda path: (
      self.missing if path == "setpaths.sh" else original(path))
    return self.missing

  def __exit__(self, *args):
    import libtbx.load_env
    del libtbx.env.under_build

def exercise_text_job_ids():
  """Job ids and polling output are text on every platform."""
  box = sandbox()
  try:
    box.phenix_env() # an environment script whatever the build layout
    from libtbx.queuing_system_utils import processing
    for platform, submit, poll, expected in PLATFORMS:
      job = processing.qsub(
        target=get_target("marker"),
        name="tst_" + platform,
        platform=platform,
        command=box.command(submit),
        polling_command=box.command(*poll))
      job.start()
      assert isinstance(job.jobid, str), (platform, repr(job.jobid))
      assert job.jobid == expected, (platform, repr(job.jobid))
      # the GUI formats the id into messages and into a qdel command line
      assert "%s" % job.jobid == expected, (platform, repr(job.jobid))
      if platform == "condor":
        with open(job.name + ".condor.script") as f:
          assert f.read().startswith("#!/bin/sh\n"), "shebang must come first"
      # polling output must be decoded, and tolerate non-UTF-8 bytes
      assert job.is_alive() is True, platform
      job.status.cleanup()
  finally:
    box.cleanup()

def expect(exception, call, *args, **kwargs):
  """Call `call`, which must raise `exception`, and return the message."""
  try:
    call(*args, **kwargs)
  except exception as e:
    return str(e)
  raise AssertionError("%s did not raise %s" % (
    getattr(call, "__name__", call), exception.__name__))

def sourced_with_source(script):
  """Lines using the non-POSIX `source` command, which dash rejects."""
  return [line for line in script.splitlines() if line.startswith("source ")]

def exercise_environment_script():
  """
  Jobs source the submitting build's setpaths.sh, or $PHENIX/phenix_env.sh
  where there is none, and a script that does not exist is an error at
  submission time.
  """
  box = sandbox()
  try:
    import libtbx.load_env
    from libtbx import runtime_utils
    from libtbx.queuing_system_utils import processing
    os.environ.pop("PHENIX", None)
    os.environ.pop("PHENIX_CUSTOM_ENV", None)
    setpaths = libtbx.env.under_build("setpaths.sh")
    env_script = box.phenix_env()
    if os.path.isfile(setpaths):
      # the build's own script wins even with PHENIX exported for an
      # installation, so a developer's queue jobs run the developer's build
      assert runtime_utils.get_environment_script() == setpaths
      assert processing.get_libtbx_env_setpaths() == setpaths
      # PHENIX_CUSTOM_ENV skips phenix_env.sh.  It is an undocumented
      # developer override from 2010 that nothing has honoured since 2012,
      # and it will be removed in the future; these checks go with it.
      os.environ["PHENIX_CUSTOM_ENV"] = "1"
      assert runtime_utils.get_environment_script() == setpaths
      del os.environ["PHENIX_CUSTOM_ENV"]
    with without_setpaths(box.dir):
      # a Phenix 2.x installation: phenix_env.sh is all there is
      assert runtime_utils.get_environment_script() == env_script
      assert processing.get_libtbx_env_setpaths() == env_script
      job = processing.qsub(
        target=get_target("marker"),
        platform="sge",
        command=box.command("qsub"),
        polling_command=box.command("qstat", "-j"))
      job.start()
      script = box.submitted_script()
      assert (". %s\n" % shlex.quote(env_script)) in script, script
      assert not sourced_with_source(script), script
      job.status.cleanup()
      # detached (non-queue) GUI jobs use the same script
      run_script = os.path.join(box.dir, "run.sh")
      runtime_utils.write_run_script(run_script, ["echo", "hi"])
      with open(run_script) as f:
        text = f.read()
      assert (". %s\n" % shlex.quote(env_script)) in text, text
      assert not sourced_with_source(text), text
      os.environ["PHENIX_CUSTOM_ENV"] = "1" # to be removed, see above
      message = expect(Sorry, runtime_utils.get_environment_script)
      assert "PHENIX_CUSTOM_ENV" in message, message
      del os.environ["PHENIX_CUSTOM_ENV"]
      # a script that does not exist is reported here, instead of being
      # sourced on the cluster where the job dies with a shell error
      os.environ["PHENIX"] = os.path.join(box.dir, "no_such_dir")
      message = expect(Sorry, runtime_utils.get_environment_script)
      assert "phenix_env.sh" in message and "setpaths.sh" in message, message
      # a plain cctbx tree has no phenix_env.sh to point the user at
      del os.environ["PHENIX"]
      message = expect(Sorry, runtime_utils.get_environment_script)
      assert "setpaths.sh" in message, message
      assert "phenix_env.sh" not in message, message
      # the error leaves nothing behind: no half-written run script and no
      # pickled target that join would never remove
      half = os.path.join(box.dir, "half.sh")
      expect(Sorry, runtime_utils.write_run_script, half, ["echo"])
      assert not os.path.exists(half)
      from libtbx.queuing_system_utils import generic
      job = generic.qsub(
        target=get_target("marker"),
        name="tst_missing",
        platform="sge",
        command=os.path.join(box.dir, "qsub"))
      expect(Sorry, job.start)
      assert not glob.glob(os.path.join(box.dir, "tst_missing*"))
      # the scheduling layer reports it like any other queue setup problem
      from libtbx.scheduling import SetupError
      from libtbx.scheduling import cluster_handler
      message = expect(SetupError, cluster_handler.JobFactory,
        platform="sge",
        command=box.command("qsub"),
        poller=processing.sge_evaluate(command=box.command("qstat", "-j")))
      assert "setpaths.sh" in message, message
  finally:
    box.cleanup()

def exercise_job_script_runs():
  """
  The generated job scripts run under /bin/sh and dash with a minimal
  environment and execute the target, so the environment line must work.
  """
  box = sandbox()
  try:
    bin_dir = dispatcher_dir()
    if bin_dir is None:
      print("Skipping job script check: libtbx.python not found")
      return
    box.phenix_env(bin_dir=bin_dir)
    from libtbx.queuing_system_utils import processing
    from libtbx.queuing_system_utils import generic
    shells = [s for s in ("/bin/sh", "/bin/dash") if os.path.exists(s)]
    # hide the build's setpaths.sh so the fake phenix_env.sh is what runs
    with without_setpaths(box.dir):
      # processing: the Phenix GUI path
      marker = os.path.join(box.dir, "processing.marker")
      job = processing.qsub(
        target=get_target(marker),
        platform="sge",
        command=box.command("qsub"),
        polling_command=box.command("qstat", "-j"))
      job.start()
      for shell in shells:
        rc, out, err = box.run_submitted_script(shell)
        assert os.path.isfile(marker), (shell, rc, out, err)
        assert "marker written" in out, (shell, rc, out, err)
        assert "not found" not in err, (shell, rc, out, err)
        assert "exit_status 0" in err, (shell, rc, out, err)
        os.remove(marker)
      job.status.cleanup()
      # generic: the older module used by phenix.autosol
      marker = os.path.join(box.dir, "generic.marker")
      job = generic.qsub(
        target=get_target(marker),
        name="tst_generic",
        platform="sge",
        command=os.path.join(box.dir, "qsub"))
      job.start()
      assert job.jobid == 12345, repr(job.jobid)
      for shell in shells:
        rc, out, err = box.run_submitted_script(shell)
        assert rc == 0 and os.path.isfile(marker), (shell, rc, out, err)
        assert "marker written" in out, (shell, rc, out, err)
        assert "not found" not in err, (shell, rc, out, err)
        os.remove(marker)
      # joining releases the submission pipes, so nothing leaks
      with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        job.join()
        del job
        gc.collect()
      leaked = [w for w in caught if issubclass(w.category, ResourceWarning)]
      assert not leaked, [str(w.message) for w in leaked]
  finally:
    box.cleanup()

def exercise_synchronous_submission():
  """
  With qsub -sync y, as phenix.autosol's run_parallel submits, qsub stays
  alive until the cluster job ends.  start must still return at once so
  that the jobs run in parallel, and join releases the pipe afterwards.
  """
  box = sandbox()
  try:
    box.phenix_env()
    from libtbx.queuing_system_utils import generic
    qsub_slow = os.path.join(box.dir, "qsub_slow")
    finished = os.path.join(box.dir, "qsub_finished") # touched as qsub exits
    makers = [
      # the base Job, as run_parallel builds it; it never parses a job id
      (lambda: generic.Job(
        name="tst_sync",
        target=get_target("marker"),
        qinterface=generic.sge_interface(command=qsub_slow)), None),
      # SGEJob reads the job id from the same pipe
      (lambda: generic.qsub(
        target=get_target("marker"),
        name="tst_sync_sge",
        platform="sge",
        command=qsub_slow,
        asynchronous=False), 12345),
      ]
    for make, jobid in makers:
      job = make()
      job.start()
      assert not os.path.exists(finished), "start waited for qsub to finish"
      assert job.is_alive()
      assert job.jobid == jobid, repr(job.jobid)
      with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        job.join()
        job.join() # a second join is harmless
        del job
        gc.collect()
      leaked = [w for w in caught if issubclass(w.category, ResourceWarning)]
      assert not leaked, [str(w.message) for w in leaked]
      assert os.path.exists(finished), "join returned before qsub exited"
      os.remove(finished)
  finally:
    box.cleanup()

def run():
  if sys.platform == "win32":
    print("Skipping tst_processing: needs /bin/sh")
    print("OK")
    return
  exercise_text_job_ids()
  exercise_environment_script()
  exercise_job_script_runs()
  exercise_synchronous_submission()
  print("OK")

if __name__ == "__main__":
  run()
