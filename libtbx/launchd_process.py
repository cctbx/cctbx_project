"""Start a GUI program as a macOS launchd job, in a coalition of its own.

A process inherits its parent's resource coalition. From macOS 26, when an
application quits, the system acts on the other processes left in its
coalition: macOS 26 prompts the user, and macOS 27 terminates them when the
application's background activity is off (System Settings > General > Login
Items & Extensions). A viewer such as Coot or ChimeraX started as a child of
a program run from a terminal shares the terminal's coalition, so quitting
the viewer can SIGTERM tmux and every shell in it, taking the program that
started the viewer down with them. launchd starts each job in a coalition of
its own.

:func:`spawn` starts a program that way where it matters (:func:`use_launchd`)
and otherwise with the caller's own ``Popen``, as it also does where there is
no GUI login session to start a job in (an SSH login with nobody logged in
at the screen, or sudo; nothing there can quit as an application).
:func:`launch` starts a launchd job and returns a :class:`LaunchdProcess`,
which has the part of the ``subprocess.Popen`` interface a launcher needs:
``args``, ``pid``, ``returncode``, ``poll()``, ``wait()``, ``communicate()``,
``send_signal()``, ``terminate()`` and ``kill()``. Its ``stop()`` also ends
whatever the program left running in the job's coalition, such as a copy of
itself relaunched under a new pid (a browser applying an update).

Accepted consequences: a launchd job is its own responsible process, so the
program asks for privacy (TCC) access, e.g. to ~/Documents, itself instead of
using the terminal's grants; and it gets launchd's open-file limit (256), as
a program started from the Dock does.
"""

from __future__ import absolute_import, division, print_function

import os
import sys

LAUNCHCTL = "/bin/launchctl"

# Seconds one launchctl call may take before it counts as failed.
LAUNCHCTL_TIMEOUT = 10.0

# Seconds to wait for launchd to report the pid of a new job.
PID_TIMEOUT = 5.0

# Seconds to wait, once launchd has reported the pid, for the program to lead
# its own process group (xpcproxy sets it up for a few milliseconds first).
PGID_TIMEOUT = 1.0

# Seconds LaunchdProcess.stop() waits for the job to end after its SIGKILLs.
STOP_WAIT = 5.0


class LaunchError(RuntimeError):
  """launch() could not start the program: launchd refused or lost the job,
  or launchctl did not answer. (A program that cannot be found raises
  FileNotFoundError or PermissionError, as from ``Popen``.)"""


class LaunchdRefused(LaunchError):
  """launchd refused the job because this login session has no GUI domain
  (gui/<uid>) to start it in: an SSH login with nobody logged in at the
  screen, or sudo. Nothing was loaded or started."""


class ExitedAtOnce(LaunchError):
  """The program exited before launchd reported its pid, so there is no
  process to hand back; ``returncode`` is its exit, as ``Popen`` gives it."""
  def __init__(self, name, returncode):
    super(ExitedAtOnce, self).__init__(
      "%s exited as soon as it started (%s)"
      % (name, describe_returncode(returncode)))
    self.name = name
    self.returncode = returncode

  def __reduce__(self):
    # The state carries what was added since, such as add_note()'s notes.
    return ExitedAtOnce, (self.name, self.returncode), self.__dict__


def use_launchd():
  """True on macOS 26 (Darwin 25) and later.

  Reads the Darwin release rather than ``platform.mac_ver()``, which a Python
  built against an old SDK sees as 10.16 on any macOS 11 or later (SDK
  compatibility mode). An unreadable release counts as new: a launchd job
  works on every macOS."""
  if sys.platform != "darwin":
    return False
  try:
    return int(os.uname().release.split(".")[0]) >= 25
  except ValueError:
    return True


def launchctl(*args):
  """Run ``launchctl args``; return ``(returncode, stdout, stderr)``."""
  import subprocess
  result = subprocess.run([LAUNCHCTL] + [str(a) for a in args],
                          stdin=subprocess.DEVNULL, capture_output=True,
                          text=True, timeout=LAUNCHCTL_TIMEOUT)
  return result.returncode, result.stdout, result.stderr


def job_target(label):
  """The ``gui/<uid>/<label>`` service target of a job in this login session."""
  return "gui/%d/%s" % (os.getuid(), label)


def job_label(label_prefix, owner_pid):
  """A unique job label that records the pid of the process starting it."""
  import uuid
  return "%s.%d.%s" % (label_prefix, owner_pid, uuid.uuid4().hex[:12])


def job_state(label):
  """``(pid, last)`` of job ``label``, from its row of launchd's job table:
  the pid of its running process or None, and its last exit as ``Popen``
  gives a returncode (-N after signal N; 0 before any run, so only a
  non-zero value proves the job ran and exited). ``(None, None)`` if it is
  not loaded. The table, unlike a job's own listing, has no field a job's
  arguments could imitate."""
  status, out, _ = launchctl("list")
  if status != 0:
    return None, None
  for line in out.splitlines():
    fields = line.split("\t")             # PID, last exit status, label
    if len(fields) == 3 and fields[2] == label:
      try:
        last = int(fields[1])
      except ValueError:
        last = None
      return (int(fields[0]) if fields[0].isdigit() else None), last
  return None, None


def job_runs(label):
  """How many times launchd has started job ``label`` (0 before its first
  run, or if it is not loaded): what tells a job that has not run yet from
  one that ran and exited 0. Read from ``launchctl print``, whose own "runs"
  line follows the job's arguments and environment (either could hold a
  line that reads the same), and which shows that environment, so its
  output is parsed here and never shown."""
  import re
  status, out, _ = launchctl("print", job_target(label))
  runs = re.findall(r"^\truns = (\d+)$", out, re.M) if status == 0 else []
  return int(runs[-1]) if runs else 0


def _gui_domain_exists():
  """True if this login session has a GUI domain (gui/<uid>) to start jobs
  in; there is none over SSH without a login at the screen, or for root."""
  try:
    return launchctl("print", "gui/%d" % os.getuid())[0] == 0
  except Exception:
    return True                        # unknown: do not assume it is absent


def describe_returncode(code):
  """'exit status N' or 'signal N' for a Popen-style returncode; a Windows
  crash code (above 0xFFFF) in hex, as Windows documents it."""
  if code < 0:
    return "signal %d" % -code
  if code > 0xFFFF:
    return "exit status 0x%08X" % code
  return "exit status %d" % code


def _pid_exists(pid):
  """True if a process ``pid`` exists (owned by anyone)."""
  try:
    os.kill(pid, 0)
  except ProcessLookupError:
    return False
  except OSError:                      # EPERM: another user's process
    return True
  return True


def unload_stale_jobs(label_prefix):
  """Unload finished jobs labelled ``label_prefix`` whose starter has exited.

  A viewer is usually left running when the program that started it exits,
  so when the viewer exits there is no one left to unload its job. A job
  whose program is still running, or whose starter is still alive, is left
  alone."""
  status, out, _ = launchctl("list")
  if status != 0:
    return
  head = label_prefix + "."
  for line in out.splitlines():
    fields = line.split("\t")             # PID, last exit status, label
    if len(fields) != 3 or fields[0] != "-" or not fields[2].startswith(head):
      continue
    owner_pid = fields[2][len(head):].split(".")[0]
    if owner_pid.isdigit() and not _pid_exists(int(owner_pid)):
      launchctl("bootout", job_target(fields[2]))


_libproc_handle = []

def _libproc():
  """macOS: libproc, loaded once. Raises OSError where there is none."""
  if not _libproc_handle:
    import ctypes
    lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
    lib.proc_pidinfo.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
                                 ctypes.c_void_p, ctypes.c_int]
    lib.proc_listallpids.argtypes = [ctypes.c_void_p, ctypes.c_int]
    _libproc_handle.append(lib)
  return _libproc_handle[0]


def _resource_coalition(pid):
  """The resource coalition id of process ``pid`` (of any user), or None if
  it cannot be read."""
  import ctypes
  try:
    info = (ctypes.c_uint64 * 5)()     # coalition_id[2] + 3 reserved
    size = ctypes.sizeof(info)
    if _libproc().proc_pidinfo(pid, 20, 0, info, size) != size:  # COALITION
      return None
  except (OSError, AttributeError):
    return None
  return int(info[0])


def _may_signal(pid):
  """True if this process may signal process ``pid`` (an existing process
  of its own user's)."""
  try:
    os.kill(pid, 0)
  except OSError:
    return False
  return True


def _coalition_pids(coalition):
  """The pids of the live processes in resource coalition ``coalition``
  that this process may signal. Raises OSError if the processes cannot be
  listed."""
  import ctypes
  lib = _libproc()
  count = lib.proc_listallpids(None, 0)
  if count >= 0:
    pids = (ctypes.c_int * (count + 256))()
    count = lib.proc_listallpids(pids, ctypes.sizeof(pids))
  if count < 0:
    code = ctypes.get_errno()
    raise OSError(code, "proc_listallpids failed: %s" % os.strerror(code))
  return sorted(pid for pid in pids[:count]
                if pid > 0 and _resource_coalition(pid) == coalition and
                _may_signal(pid))


def _bootout_quietly(label):
  """Unload job ``label``, stopping its program if it runs. Best effort: a
  failure here must not replace the error that made the caller clean up."""
  try:
    launchctl("bootout", job_target(label))
  except Exception:
    pass


def _current_umask():
  """The process umask. Reading it means setting it, so for the moment
  between the two calls it is 077: a file another thread creates then is
  private, never readable by others."""
  mask = os.umask(0o077)
  os.umask(mask)
  return mask


def utf8_environment(env, warn=None):
  """``env`` without the variables a launchd plist cannot carry: plist strings
  are Unicode, so a name or value that is not valid UTF-8 (held by Python as
  surrogate escapes) is left out, and ``warn`` gets one message naming them."""
  kept, dropped = {}, []
  for key, value in env.items():
    try:
      key.encode("utf-8")
      value.encode("utf-8")
    except UnicodeEncodeError:
      dropped.append(key.encode("ascii", "backslashreplace").decode("ascii"))
      continue
    kept[key] = value
  if dropped and warn is not None:
    warn("WARNING: not passing these environment variables, whose names or "
         "values are not valid UTF-8: %s" % ", ".join(dropped))
  return kept


class LaunchdProcess(object):
  """The program a launchd job runs, with the part of the ``Popen`` interface
  a launcher uses: ``args``, ``pid``, ``returncode``, ``poll()``,
  ``wait(timeout)`` and ``communicate(timeout)`` (no pipes, so it returns
  ``(None, None)``); the timeouts raise ``subprocess.TimeoutExpired``.

  launchd, not the caller, is the program's parent, so its exit is seen
  through a kqueue NOTE_EXIT filter on the pid instead of ``waitpid``. As with
  ``Popen``, a ``wait()`` blocked in another thread holds the lock, and
  ``poll()`` returns None until that waiter has recorded the exit.
  ``returncode`` follows ``Popen``: the exit code, or -N after signal N, read
  from launchd's record of the exit (0 if launchd has none, as ``Popen``
  gives when the status is lost). Recording the exit also unloads the
  finished job from launchd, before ``returncode`` is set.

  ``send_signal()``, ``terminate()`` and ``kill()`` signal the program as
  ``Popen``'s do. ``job_pids()`` and ``stop()`` reach the job's resource
  coalition (``coalition``), which everything the program starts inherits:
  a copy of itself it relaunches under a new pid, even in a new session,
  stays there, while its pid, process group and job do not cover it (a copy
  started through LaunchServices gets a coalition of its own and is not
  reached). So does a daemon a member starts on demand (tmux, gpg-agent),
  which ``stop()`` then ends too. The coalition is taken as the job's only
  if launchd reports ``pid`` as job ``label``'s running program and it is
  not the caller's; otherwise ``coalition`` is None and only the program
  itself is reached, never the caller's other processes."""

  def __init__(self, label, pid, args=None):
    import errno
    import select
    import threading
    self.label = label
    self.pid = pid
    self.args = args
    self.returncode = None
    self.coalition = _resource_coalition(pid)
    if self.coalition is not None and (
        self.coalition == _resource_coalition(os.getpid()) or
        not self._is_jobs_program()):
      self.coalition = None            # not a launchd job: never ours to stop
    self._lock = threading.Lock()
    self._kq = select.kqueue()
    try:
      self._kq.control([select.kevent(pid, filter=select.KQ_FILTER_PROC,
                                      flags=select.KQ_EV_ADD,
                                      fflags=select.KQ_NOTE_EXIT)], 0, 0)
    except OSError as exc:
      if exc.errno != errno.ESRCH:
        self._kq.close()
        raise
      self._exited()                   # the program has already exited

  def _is_jobs_program(self):
    try:
      return job_state(self.label)[0] == self.pid
    except Exception:                  # unknown: not ours to stop
      return False

  def _gone(self):
    # The program has exited: no process has its pid, or the one that has it
    # now is in another coalition than the job's (the pid was reused).
    if not _pid_exists(self.pid):
      return True
    if self.coalition is None:
      return False
    now = _resource_coalition(self.pid)
    return now is not None and now != self.coalition

  def _exited(self):
    # returncode is set last, so whoever sees it finds the job unloaded. An
    # interruption part way leaves the kqueue closed and returncode None;
    # the next _check() starts again.
    self._kq.close()
    code = self._exit_code()
    _bootout_quietly(self.label)
    self.returncode = code

  def _exit_code(self):
    """The returncode from launchd's record of the exit, once launchd has
    recorded it (it may lag the kqueue event by a moment)."""
    import time
    deadline = time.monotonic() + 1.0
    while True:
      try:
        pid, last = job_state(self.label)
      except Exception:
        return 0
      if pid != self.pid or time.monotonic() >= deadline:
        break
      time.sleep(0.02)
    return last or 0

  def _check(self, timeout):
    # Caller holds self._lock. The exit event is waited for in slices of at
    # most a second, each after checking that the process still exists: an
    # event taken but not recorded (a Ctrl-C at that moment) is not sent
    # again.
    if timeout is None or timeout > 1.0:
      timeout = 1.0
    if self.returncode is None and (
        self._kq.closed or self._gone() or
        self._kq.control(None, 1, timeout)):
      self._exited()
    return self.returncode

  def poll(self):
    if self.returncode is not None or not self._lock.acquire(False):
      return self.returncode
    try:
      return self._check(0)
    finally:
      self._lock.release()

  def wait(self, timeout=None):
    import subprocess
    import threading
    import time
    deadline = None if timeout is None else time.monotonic() + timeout
    def left():
      return None if deadline is None else max(0.0, deadline - time.monotonic())
    while self.returncode is None:
      # Another thread's wait() may hold the lock until the exit, or until
      # its own deadline; whatever time is left counts from now.
      if self._lock.acquire(timeout=-1 if deadline is None else
                            min(left(), threading.TIMEOUT_MAX)):
        try:
          self._check(left())
        finally:
          self._lock.release()
      if self.returncode is None and left() == 0.0:
        raise subprocess.TimeoutExpired(self.args, timeout)
    return self.returncode

  def communicate(self, input=None, timeout=None):
    # As Popen's without pipes: no output, and no stdin for the input.
    self.wait(timeout)
    return None, None

  def send_signal(self, sig):
    """Send signal ``sig`` to the program, unless it has exited."""
    if self.poll() is None:
      try:
        os.kill(self.pid, sig)
      except ProcessLookupError:
        pass                           # it exited meanwhile

  def terminate(self):
    import signal
    self.send_signal(signal.SIGTERM)

  def kill(self):
    import signal
    self.send_signal(signal.SIGKILL)

  def job_pids(self):
    """The pids of the job's live processes that this process may signal:
    the program, while it runs, and whatever it started that is still in
    the job's coalition. Without a coalition of its own, or if the
    coalition holds this process or its parent, the program alone. Raises
    OSError if the processes cannot be listed."""
    if self.coalition is not None:
      pids = _coalition_pids(self.coalition)
      if os.getpid() not in pids and os.getppid() not in pids:
        return pids
    return [self.pid] if self.poll() is None and _pid_exists(self.pid) \
      else []

  def stop(self, timeout=10.0):
    """End the job: SIGTERM to every process in ``job_pids()``, SIGKILL to
    those left after ``timeout`` seconds (None: wait for the SIGTERM
    without a limit), then record the program's exit. Returns its
    returncode. Raises subprocess.TimeoutExpired if the program has not
    exited STOP_WAIT seconds after that, and OSError if the job's processes
    cannot be listed."""
    import signal
    import time
    for sig, wait_s in ((signal.SIGTERM, timeout), (signal.SIGKILL,
                                                    STOP_WAIT)):
      deadline = None if wait_s is None else time.monotonic() + wait_s
      pids = self.job_pids()
      for pid in pids:
        try:
          os.kill(pid, sig)
        except (ProcessLookupError, PermissionError):
          pass                         # exited, or no longer ours
      while pids and (deadline is None or time.monotonic() < deadline):
        time.sleep(0.05)
        pids = self.job_pids()
      if not pids:
        break
    return self.wait(STOP_WAIT)


def _find_program(program, env):
  """The file ``Popen`` would run for ``program`` with environment ``env``:
  a name with a directory part as given; a bare name looked up on ``env``'s
  PATH as ``Popen`` looks it up (``os.get_exec_path``; an empty entry is the
  working directory). launchd would search its own PATH. Raises
  FileNotFoundError, or PermissionError if a match is not executable, as
  ``Popen`` does."""
  import errno
  program = str(program)
  if os.path.dirname(program):
    candidates = [program]
  else:
    candidates = [os.path.join(d, program) for d in os.get_exec_path(env)]
  denied = False
  for path in candidates:
    if os.path.isfile(path) and os.access(path, os.X_OK):
      try:
        return os.path.abspath(path)
      except OSError:                  # relative, and the cwd is deleted
        return path
    denied = denied or os.path.exists(path)
  if denied:
    raise PermissionError(errno.EACCES, os.strerror(errno.EACCES), program)
  raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), program)


def _await_own_group(pid):
  """Wait (at most PGID_TIMEOUT) until ``pid`` leads its own process group.
  launchd reports the pid while xpcproxy is still setting the program up in
  group 1, and a ``killpg(pid)`` then fails."""
  import time
  deadline = time.monotonic() + PGID_TIMEOUT
  while time.monotonic() < deadline:
    try:
      if os.getpgid(pid) == pid:
        return
    except OSError:                    # it has exited
      return
    time.sleep(0.001)


def launch(argv, env=None, log_path=None, label_prefix="org.cctbx.launchd",
           warn=None, stderr_path=None, name=None):
  """Start ``argv`` as a launchd job and return a :class:`LaunchdProcess`.

  As with ``Popen``, ``argv[0]`` is found on the PATH of ``env`` (default:
  the caller's environment), and the program runs in the caller's working
  directory (the root directory if that has been deleted) with the caller's
  umask and ``env``; what it starts outlives it. stdin is the null device;
  stdout is appended to ``log_path`` and stderr to ``stderr_path`` (default:
  ``log_path``; both default to the null device). ``warn`` gets warnings;
  ``name`` (default: the program's file name) names the program in errors.
  It returns once the program leads its own process group, waiting at most
  PGID_TIMEOUT. Raises on failure, leaving no job loaded: an OSError where
  ``Popen`` would raise one (FileNotFoundError or PermissionError for the
  program, or running out of file descriptors), and otherwise a
  LaunchError: LaunchdRefused if there is no GUI session to start the job
  in, ExitedAtOnce if the program exited before its pid was read."""
  import plistlib
  import subprocess
  import tempfile
  import time
  if env is None:
    env = dict(os.environ)
  program = _find_program(argv[0], env)
  name = name or os.path.basename(str(argv[0]))
  try:
    unload_stale_jobs(label_prefix)
  except Exception:                    # housekeeping; never blocks the launch
    pass
  label = job_label(label_prefix, os.getpid())
  job = {
    "Label": label,
    "Program": program,
    "ProgramArguments": [str(a) for a in argv],
    "EnvironmentVariables": utf8_environment(env, warn),
    "Umask": _current_umask(),
    "RunAtLoad": True,
    "KeepAlive": False,
    "ProcessType": "Interactive",
    "AbandonProcessGroup": True,
    "StandardInPath": os.devnull,
    "StandardOutPath": str(log_path or os.devnull),
    "StandardErrorPath": str(stderr_path or log_path or os.devnull),
  }
  try:
    job["WorkingDirectory"] = os.getcwd()
  except OSError:                      # deleted; launchd starts it in /
    pass
  attempted = refused = False
  try:
    fd, plist_path = tempfile.mkstemp(prefix="launchd_process_",
                                      suffix=".plist")
    try:
      # Binary, because an XML plist cannot hold the control characters some
      # environment variables carry (escape sequences in LESS_TERMCAP_*).
      with os.fdopen(fd, "wb") as f:
        plistlib.dump(job, f, fmt=plistlib.FMT_BINARY)
      attempted = True
      status, _, err = launchctl("bootstrap", "gui/%d" % os.getuid(),
                                 plist_path)
    finally:
      try:
        os.remove(plist_path)          # launchd has read it
      except OSError:
        pass
    if status != 0:
      refused = True
      # launchctl's first line; the next advises re-running as root, which
      # would leave no GUI domain to start the job in at all.
      lines = err.strip().splitlines()
      message = "launchctl bootstrap exited with status %d: %s" % (
        status, lines[0] if lines else "")
      if not _gui_domain_exists():
        raise LaunchdRefused(message)
      raise LaunchError(message)
    deadline = time.monotonic() + PID_TIMEOUT
    while True:
      pid, last_status = job_state(label)
      if pid is not None:
        _await_own_group(pid)
        return LaunchdProcess(label, pid, args=list(argv))
      # The job table's last exit is 0 both before the first run and after
      # an exit 0; the run count tells them apart.
      if last_status or job_runs(label):
        pid, last_status = job_state(label)   # it may have exited just now
        if pid is None:
          raise ExitedAtOnce(name, last_status or 0)
        continue
      if time.monotonic() >= deadline:
        raise LaunchError("launchd did not report a running %s" % name)
      time.sleep(0.05)
  except BaseException as exc:        # Ctrl-C included
    if attempted and not refused:
      _bootout_quietly(label)          # it may be loaded, and running
    if isinstance(exc, subprocess.TimeoutExpired):
      raise LaunchError(str(exc)) from exc
    raise


def spawn(argv, popen, env=None, log_path=None, stderr_path=None,
          label_prefix="org.cctbx.launchd", warn=None, name=None):
  """Start a GUI program: as a launchd job (:func:`launch`, with the other
  arguments) where :func:`use_launchd` says so, and otherwise with
  ``popen()``, the caller's own ``Popen`` call.

  Where there is no GUI login session to start a job in (an SSH login with
  nobody logged in at the screen, or sudo), it warns and falls back to
  ``popen()``: the program then shares the caller's coalition, as before
  launchd was used, which matters only in a GUI session. Any other launch
  failure is raised, a refusal inside a GUI session included, since falling
  back there would bring back the hazard."""
  if use_launchd():
    try:
      return launch(argv, env=env, log_path=log_path, stderr_path=stderr_path,
                    label_prefix=label_prefix, warn=warn, name=name)
    except LaunchdRefused as exc:
      if warn is not None:
        warn(refused_warning(name or os.path.basename(str(argv[0])), exc))
  return popen()


def refused_warning(name, exc):
  """The warning for starting program ``name`` as a child process after
  LaunchdRefused ``exc`` (no GUI login session)."""
  return ("WARNING: no GUI login session to start %s in as a launchd job "
          "(%s); starting it as a child process instead" % (name, exc))
