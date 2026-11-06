"""Run jobs as local subprocesses (spec 7)."""

import ctypes
import os
import signal
import subprocess
import sys
import threading

from libtbx.jobs.backends.base import Backend, PollResult, SubmitError


# The exit code ``TerminateJobObject`` gives every process of a cancelled
# job on Windows; 1, as ``Popen.terminate`` and ``taskkill /F``.
WINDOWS_TERMINATED_EXIT_CODE = 1


class _LocalHandle(object):
  """The subprocess behind a job.

  Attributes
  ----------
  proc : subprocess.Popen
  pid : int
  job_object : int or None
      Windows only: the job object holding the job's processes, until
      ``poll`` sees the job end and closes it.
  signalled : bool
      True once this backend has signalled the process group, or on
      Windows terminated the job object.
  """

  __slots__ = ("proc", "pid", "job_object", "signalled")

  def __init__(self, proc, job_object=None):
    self.proc = proc
    self.pid = proc.pid
    self.job_object = job_object
    self.signalled = False


class _BasicLimitInformation(ctypes.Structure):
  """``JOBOBJECT_BASIC_LIMIT_INFORMATION`` with fixed-width fields."""
  _fields_ = [
    ("PerProcessUserTimeLimit", ctypes.c_int64),
    ("PerJobUserTimeLimit", ctypes.c_int64),
    ("LimitFlags", ctypes.c_uint32),
    ("MinimumWorkingSetSize", ctypes.c_size_t),
    ("MaximumWorkingSetSize", ctypes.c_size_t),
    ("ActiveProcessLimit", ctypes.c_uint32),
    ("Affinity", ctypes.c_size_t),
    ("PriorityClass", ctypes.c_uint32),
    ("SchedulingClass", ctypes.c_uint32),
  ]


class _IoCounters(ctypes.Structure):
  """``IO_COUNTERS``."""
  _fields_ = [(name, ctypes.c_uint64) for name in (
    "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
    "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class _ExtendedLimitInformation(ctypes.Structure):
  """``JOBOBJECT_EXTENDED_LIMIT_INFORMATION``."""
  _fields_ = [
    ("BasicLimitInformation", _BasicLimitInformation),
    ("IoInfo", _IoCounters),
    ("ProcessMemoryLimit", ctypes.c_size_t),
    ("JobMemoryLimit", ctypes.c_size_t),
    ("PeakProcessMemoryUsed", ctypes.c_size_t),
    ("PeakJobMemoryUsed", ctypes.c_size_t),
  ]


class _JobObjects(object):
  """One Win32 job object per job.

  Every process a job's process starts joins its job object, so
  ``TerminateJobObject`` ends the whole tree, as ``os.killpg`` does on
  POSIX. A job object kills its processes when its last handle closes,
  so jobs also end when this process exits without ``on_shutdown``, as
  under ``TerminateProcess``, which is how Windows MCP clients stop a
  server.
  """

  # CREATE_SUSPENDED (absent from the subprocess module): the child joins
  # its job object before it can start a process outside it.
  # CREATE_NO_WINDOW: a hidden console of its own, so console programs the
  # job starts open no window and Ctrl+C in the host's console does not
  # reach the job (POSIX: start_new_session).
  creation_flags = 0x00000004 | 0x08000000

  _EXTENDED_LIMIT_INFORMATION_CLASS = 9
  # KILL_ON_JOB_CLOSE, and DIE_ON_UNHANDLED_EXCEPTION: a crashed process
  # exits instead of waiting on an error-reporting dialog nobody sees.
  _LIMIT_FLAGS = 0x00002000 | 0x00000400
  # PROCESS_TERMINATE | PROCESS_SET_QUOTA (AssignProcessToJobObject) |
  # PROCESS_SUSPEND_RESUME (NtResumeProcess)
  _PROCESS_ACCESS = 0x0001 | 0x0100 | 0x0800

  def __init__(self):
    from ctypes import wintypes
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    for name, argtypes, restype in (
        ("CreateJobObjectW", (wintypes.LPVOID, wintypes.LPCWSTR),
         wintypes.HANDLE),
        ("SetInformationJobObject",
         (wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID, wintypes.DWORD),
         wintypes.BOOL),
        ("OpenProcess", (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD),
         wintypes.HANDLE),
        ("AssignProcessToJobObject", (wintypes.HANDLE, wintypes.HANDLE),
         wintypes.BOOL),
        ("TerminateJobObject", (wintypes.HANDLE, wintypes.UINT),
         wintypes.BOOL),
        ("CloseHandle", (wintypes.HANDLE,), wintypes.BOOL)):
      function = getattr(k32, name)
      function.argtypes = argtypes
      function.restype = restype
    # Undocumented but stable since Windows XP; psutil resumes with it.
    # ResumeThread needs the thread handle, which Popen closes.
    ntdll = ctypes.WinDLL("ntdll")
    ntdll.NtResumeProcess.argtypes = (wintypes.HANDLE,)
    ntdll.NtResumeProcess.restype = ctypes.c_long         # NTSTATUS
    self._k32 = k32
    self._ntdll = ntdll

  @staticmethod
  def _error():
    return ctypes.WinError(ctypes.get_last_error())

  def adopt(self, pid):
    """Put the suspended process ``pid`` in a new job object and resume it.

    Returns the job object's handle. Raises ``OSError`` with the process
    still suspended, or already ended by the job object's close.
    """
    k32 = self._k32
    job_object = k32.CreateJobObjectW(None, None)
    if not job_object:
      raise self._error()
    try:
      limits = _ExtendedLimitInformation()
      limits.BasicLimitInformation.LimitFlags = self._LIMIT_FLAGS
      if not k32.SetInformationJobObject(
          job_object, self._EXTENDED_LIMIT_INFORMATION_CLASS,
          ctypes.byref(limits), ctypes.sizeof(limits)):
        raise self._error()
      process = k32.OpenProcess(self._PROCESS_ACCESS, False, pid)
      if not process:
        raise self._error()
      try:
        if not k32.AssignProcessToJobObject(job_object, process):
          raise self._error()
        status = self._ntdll.NtResumeProcess(process)
        if status < 0:
          raise OSError("NtResumeProcess failed (NTSTATUS 0x%08X)"
                        % (status & 0xFFFFFFFF))
      finally:
        k32.CloseHandle(process)
    except BaseException:
      k32.CloseHandle(job_object)
      raise
    return job_object

  def terminate(self, job_object):
    """End every process in the job object; raises ``OSError``."""
    if not self._k32.TerminateJobObject(job_object,
                                        WINDOWS_TERMINATED_EXIT_CODE):
      raise self._error()

  def close(self, job_object):
    """Close the handle, which ends any process still in the job object."""
    self._k32.CloseHandle(job_object)


class LocalBackend(Backend):
  """Spawn each job on this machine, detached from the host's terminal.

  On POSIX each job is a session of its own; a cancel signals its
  process group, SIGTERM and then SIGKILL. On Windows each job is a job
  object; a cancel ends all its processes at once with exit code
  ``WINDOWS_TERMINATED_EXIT_CODE`` (a console program cannot be asked to
  stop), and processes a job leaves behind end with it.

  Parameters
  ----------
  max_concurrent : int, optional
      Slots the manager may fill (default 8, minimum 1).
  """

  name = "local"
  handle_keys = ("pid",)

  def __init__(self, max_concurrent=8):
    self.max_concurrent = max(1, int(max_concurrent))
    self._lock = threading.Lock()
    self._job_objects = _JobObjects() if sys.platform == "win32" else None

  def capacity(self):
    return self.max_concurrent

  def submit(self, job):
    spec = job.spec
    env = dict(os.environ)
    # Do not leak a surrounding batch job's attempt count to children.
    env.pop("LIBTBX_JOBS_ATTEMPT", None)
    env.update(spec.env)
    if self._job_objects is None:
      detach = {"start_new_session": True}
    else:
      detach = {"creationflags": self._job_objects.creation_flags}
    try:
      out = open(spec.log_path, "a", encoding="utf-8")
    except OSError:
      out = None
    try:
      proc = subprocess.Popen(
        list(spec.argv),
        cwd=spec.cwd,
        env=env,
        # not the manager's stdin: under an MCP server it is the JSON-RPC
        # pipe, which a child reading stdin would drain
        stdin=subprocess.DEVNULL,
        stdout=out if out is not None else subprocess.DEVNULL,
        stderr=subprocess.STDOUT if out is not None else subprocess.DEVNULL,
        **detach)
    except OSError as exc:
      message = str(exc)
      # Windows errors name no file; a POSIX bad-cwd error names the cwd
      if spec.argv[0] not in message:
        message = "%s: %s" % (spec.argv[0], message)
      raise SubmitError(message)
    finally:
      if out is not None:
        out.close()
    job_object = None
    if self._job_objects is not None:
      try:
        job_object = self._job_objects.adopt(proc.pid)
      except BaseException as exc:
        # whatever went wrong, do not leave a suspended child behind
        proc.kill()
        proc.wait()
        if not isinstance(exc, Exception):
          raise
        raise SubmitError("cannot put the job in a Windows job object: %s"
                          % exc)
    return "running", _LocalHandle(proc, job_object)

  def poll(self, job):
    handle = job.handle
    with self._lock:
      rc = handle.proc.poll()
      if rc is not None and handle.job_object is not None:
        job_object, handle.job_object = handle.job_object, None
        self._job_objects.close(job_object)
    if rc is None:
      return PollResult("running")
    if self._job_objects is None:
      stopped = rc < 0
      reason = "killed by signal %d" % -rc if stopped else None
    else:
      # Windows cannot tell a terminated process from one that exited 1.
      stopped = handle.signalled and rc == WINDOWS_TERMINATED_EXIT_CODE
      reason = "terminated" if stopped else None
    if not handle.signalled:
      cancelled = False
    elif stopped:
      cancelled = True
    else:
      cancelled = None
    return PollResult("exited", exit_code=rc, reason=reason,
                      cancelled=cancelled)

  def _stop(self, job, signum):
    handle = job.handle
    with self._lock:
      if handle.proc.poll() is not None:
        return
      if self._job_objects is not None:
        self._job_objects.terminate(handle.job_object)
        handle.signalled = True
        return
      try:
        os.killpg(os.getpgid(handle.pid), signum)
        handle.signalled = True
      except ProcessLookupError:
        pass

  def terminate(self, job):
    self._stop(job, signal.SIGTERM)

  def kill(self, job):
    # Windows has no SIGKILL; there terminate has already ended the job.
    self._stop(job, getattr(signal, "SIGKILL", None))

  def on_shutdown(self, jobs):
    for job in jobs:
      handle = job.handle
      if handle is None:
        continue
      try:
        if self._job_objects is None:
          os.killpg(os.getpgid(handle.pid), signal.SIGTERM)
        else:
          job_object = handle.job_object      # poll may clear it meanwhile
          if job_object is not None:
            self._job_objects.terminate(job_object)
      except OSError:
        pass

  def handle_info(self, job):
    return {"pid": job.handle.pid if job.handle is not None else None}
