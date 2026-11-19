"""How libtbx.xmlrpc_utils.external_program_thread starts its program (the
Phenix GUI's Coot and PyMOL).

On macOS 26 and later the program runs as a launchd job, in a coalition of
its own (libtbx.launchd_process); elsewhere it stays a Popen child with
shell=True."""

from __future__ import absolute_import, division, print_function

import io
import os
import time
from io import StringIO

from libtbx import launchd_process, xmlrpc_utils
from libtbx.test_utils.launchd import (
  launchd_tests_can_run, resource_coalition)


def exercise_dispatch():
  """The thread uses launchd_process.launch exactly when use_launchd() says
  so, running its command through /bin/sh -c as shell=True would (a command
  given as one string included), with the caller's environment under
  use_env, and naming the program by its program_id; otherwise Popen with
  shell=True. When launchd refuses the job (no GUI login session), it warns
  and falls back to Popen. Both are faked and report the program already
  exited."""
  import subprocess
  import types
  def refuse(argv, **kw):
    raise launchd_process.LaunchdRefused("launchctl bootstrap exited with "
                                         "status 125")
  command = "phenix.start_coot --no-state-script"
  for use, launch, command_args in ((False, None, [command]),
                                    (True, None, [command]),
                                    (True, None, command),
                                    (True, refuse, [command])):
    popened, launched, log = [], [], StringIO()
    done = types.SimpleNamespace(pid=1, poll=lambda: 0, stdout=None)
    saved = (subprocess.Popen, launchd_process.use_launchd,
             launchd_process.launch)
    subprocess.Popen = lambda **kw: popened.append(kw) or done
    launchd_process.use_launchd = lambda: use
    launchd_process.launch = launch or (lambda argv, **kw: launched.append(
      (argv, kw)) or done)
    try:
      thread = xmlrpc_utils.external_program_thread(
        command_args=command_args, program_id="Coot", log=log,
        intercept_output=False, use_env=True)
      thread.run()
    finally:
      (subprocess.Popen, launchd_process.use_launchd,
       launchd_process.launch) = saved
    assert not thread.is_alive()
    if not use or launch is refuse:
      assert launched == [] and len(popened) == 1, len(launched)
      assert popened[0]["shell"] is True, sorted(popened[0])
      if launch is refuse:
        assert "status 125" in log.getvalue() and \
          "child process" in log.getvalue(), log.getvalue()
      continue
    assert popened == [], [sorted(kw) for kw in popened]
    ((argv, kw),) = launched
    assert argv == ["/bin/sh", "-c", command], argv
    assert kw["name"] == "Coot", kw["name"]
    assert kw["env"] is not None and kw["env"].get("PATH") == \
      os.environ.get("PATH"), "use_env must pass the caller's environment"


def exercise_own_coalition_and_output():
  """macOS 26+: the program runs from a coalition of its own, and with
  intercept_output its output still reaches the thread's log (through a
  file, since a launchd job cannot hand back a pipe)."""
  if not launchd_tests_can_run():
    return
  log = StringIO()
  thread = xmlrpc_utils.external_program_thread(
    command_args=["echo pid=$$; exec /bin/sleep 3"], program_id="TST",
    log=log, intercept_output=True)
  thread.start()
  deadline = time.monotonic() + 10
  while "pid=" not in log.getvalue() and time.monotonic() < deadline:
    time.sleep(0.1)
  text = log.getvalue()
  assert "pid=" in text, text
  pid = int(text.split("pid=")[1].split()[0])
  assert resource_coalition(pid) != resource_coalition(os.getpid()), \
    "the program shares the caller's coalition"
  thread.join(20)
  assert not thread.is_alive(), "the thread did not see the program exit"


class _redirected_fds(object):
  """Context manager: file descriptors 1 and 2 write to files ``out`` and
  ``err`` (a launched program's output goes there, as a Popen child's
  would); restored afterwards."""
  def __init__(self, out, err):
    self.paths = (out, err)
  def __enter__(self):
    self.saved = [os.dup(1), os.dup(2)]
    for fd, path in zip((1, 2), self.paths):
      target = os.open(path, os.O_WRONLY | os.O_CREAT, 0o600)
      os.dup2(target, fd)
      os.close(target)
  def __exit__(self, *exc):
    for fd, copy in zip((1, 2), self.saved):
      os.dup2(copy, fd)
      os.close(copy)


def _read(path):
  with open(path) as f:
    return f.read()


def exercise_output_reaches_stdout_and_stderr():
  """macOS 26+: without intercept_output, the program's stdout and stderr
  reach this process's stdout and stderr (file descriptors 1 and 2), as a
  Popen child's would, and while it still runs: the Phenix GUI never
  intercepts, and messages such as "Sorry: the coot command could not be
  located" or a Coot traceback must reach the terminal it was started
  from. Each descriptor is pointed at a file for the test."""
  if not launchd_tests_can_run():
    return
  import tempfile
  with tempfile.TemporaryDirectory() as tmp:
    paths = [os.path.join(tmp, name) for name in ("out", "err")]
    with _redirected_fds(*paths):
      thread = xmlrpc_utils.external_program_thread(
        command_args=["echo to-stdout; echo to-stderr >&2; exec /bin/sleep 3"],
        program_id="TST", log=StringIO(), intercept_output=False)
      thread.start()
      deadline = time.monotonic() + 2.5
      while time.monotonic() < deadline and \
          [_read(path) for path in paths] != ["to-stdout\n", "to-stderr\n"]:
        time.sleep(0.1)
      running = thread.is_alive()
      text = [_read(path) for path in paths]
      thread.join(20)
    assert text == ["to-stdout\n", "to-stderr\n"], text
    assert running, "the output arrived only after the program exited"
    assert not thread.is_alive(), "the thread did not see the program exit"


def exercise_output_after_exit():
  """Output written just before the program exits still arrives, read once
  the exit is seen (the polling loop never runs):
  - when launch() returns a handle whose program has already exited;
  - when launch() raises ExitedAtOnce, because the program exited before its
    pid was read. The thread takes that as the program having run, like a
    Popen child that exits at once: no Popen fallback and no traceback.
  Intercepted output that ends in a cut-off UTF-8 sequence ends in U+FFFD,
  also when that sequence was read before the exit and nothing came after
  it; if writing to the log fails, the output files are still removed. Any
  other launch failure (a refusal inside a GUI session, or no room in the
  temporary directory for the output files) is written to the log, after
  whatever the program wrote, not raised out of the thread, where nobody
  would see it. launch and use_launchd are faked; nothing starts. Runs on
  every platform."""
  import tempfile
  import types
  import subprocess
  def fake_launch(how):
    def launch(argv, **kw):
      seen.append(kw)
      with open(kw["log_path"], "ab") as f:
        f.write(b"caf\xc3" if how.startswith("cut") else
                b"partial" if how == "refused_partial" else b"to-stdout\n")
      with open(kw["stderr_path"], "ab") as f:
        f.write(b"" if how.startswith("cut") or how == "refused_partial"
                else b"to-stderr\n")
      if how == "raise":
        raise launchd_process.ExitedAtOnce("TST", 1)
      if how in ("refused", "refused_partial"):
        raise RuntimeError("launchctl bootstrap exited with status 5: "
                           "Bootstrap failed: 5: Input/output error")
      polls = [None] if how == "cut_early" else []
      return types.SimpleNamespace(
        pid=1, poll=lambda: polls.pop() if polls else 0)
    return launch
  class FailingLog(object):
    def write(self, text):
      raise OSError(28, "No space left on device")
    def flush(self):
      pass
  saved = (subprocess.Popen, launchd_process.use_launchd,
           launchd_process.launch)
  popened = []
  subprocess.Popen = lambda **kw: popened.append(kw)
  launchd_process.use_launchd = lambda: True
  try:
    with tempfile.TemporaryDirectory() as tmp:
      paths = [os.path.join(tmp, name) for name in ("out", "err")]
      for how in ("exited", "raise"):
        seen = []
        launchd_process.launch = fake_launch(how)
        thread = xmlrpc_utils.external_program_thread(
          command_args=["tst"], program_id="TST", log=StringIO(),
          intercept_output=False)
        with _redirected_fds(*paths):
          thread.run()                 # must not raise
        assert [_read(path) for path in paths] == \
          ["to-stdout\n", "to-stderr\n"], (how, [_read(p) for p in paths])
        for path in paths:
          os.remove(path)
      for how in ("cut", "cut_early"):
        seen, log = [], StringIO()
        launchd_process.launch = fake_launch(how)
        xmlrpc_utils.external_program_thread(
          command_args=["tst"], program_id="TST", log=log,
          intercept_output=True).run()
        assert log.getvalue() == u"caf\ufffd", (how, repr(log.getvalue()))
      seen, log = [], StringIO()
      launchd_process.launch = fake_launch("refused")
      xmlrpc_utils.external_program_thread(
        command_args=["tst"], program_id="TST", log=log,
        intercept_output=True).run()       # must not raise
      assert log.getvalue() == "to-stdout\nto-stderr\nWARNING: could not " \
        "start TST: launchctl bootstrap exited with status 5: Bootstrap " \
        "failed: 5: Input/output error\n", log.getvalue()
      seen, log = [], StringIO()
      launchd_process.launch = fake_launch("refused_partial")
      xmlrpc_utils.external_program_thread(
        command_args=["tst"], program_id="TST", log=log,
        intercept_output=True).run()
      assert log.getvalue().startswith("partial\nWARNING: could not start"), \
        repr(log.getvalue())
      # Not intercepted, the output goes to file descriptor 1, which is
      # where the log (sys.stdout, in the GUI) writes too.
      with _redirected_fds(*paths):
        log = io.open(1, "w", closefd=False)
        launchd_process.launch = fake_launch("refused_partial")
        xmlrpc_utils.external_program_thread(
          command_args=["tst"], program_id="TST", log=log,
          intercept_output=False).run()
        log.flush()
      assert _read(paths[0]).startswith("partial\nWARNING: could not start"), \
        repr(_read(paths[0]))
      for path in paths:
        os.remove(path)
      seen, log = [], StringIO()
      saved_mkstemp = tempfile.mkstemp
      def full(*args, **kwargs):
        raise OSError(28, "No space left on device")
      tempfile.mkstemp = full
      try:
        xmlrpc_utils.external_program_thread(
          command_args=["tst"], program_id="TST", log=log,
          intercept_output=False).run()    # must not raise
      finally:
        tempfile.mkstemp = saved_mkstemp
      assert "could not start TST" in log.getvalue() and \
        "No space left on device" in log.getvalue(), log.getvalue()
      assert seen == [], "launched without its output files"
      seen = []
      launchd_process.launch = fake_launch("exited")
      try:
        xmlrpc_utils.external_program_thread(
          command_args=["tst"], program_id="TST", log=FailingLog(),
          intercept_output=True).run()
      except OSError:
        pass
      else:
        raise AssertionError("the log error was lost")
      left = [kw[k] for kw in seen for k in ("log_path", "stderr_path")
              if os.path.exists(kw[k])]
      assert left == [], left
  finally:
    (subprocess.Popen, launchd_process.use_launchd,
     launchd_process.launch) = saved
  assert popened == [], "fell back to Popen"


def exercise_intercepted_output_not_utf8():
  """macOS 26+: intercepted output that is not valid UTF-8 still reaches
  the log, with replacement characters, instead of stopping the thread."""
  if not launchd_tests_can_run():
    return
  log = StringIO()
  thread = xmlrpc_utils.external_program_thread(
    command_args=["printf 'caf\\351\\n'; echo done"], program_id="TST",
    log=log, intercept_output=True)
  thread.start()
  thread.join(20)
  assert not thread.is_alive(), "the thread did not see the program exit"
  assert log.getvalue() == u"caf\ufffd\ndone\n", repr(log.getvalue())


def exercise():
  exercise_dispatch()
  exercise_own_coalition_and_output()
  exercise_output_reaches_stdout_and_stderr()
  exercise_output_after_exit()
  exercise_intercepted_output_not_utf8()


if __name__ == "__main__":
  exercise()
  print("OK")
