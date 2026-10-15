"""Controls for the claims that auto_queue.py makes (stage1/AUTO_QUEUE_SPEC.md,
including its Revision 2 section, which overrides the earlier text).

Claims controlled: ACTIVE resolution, init refusals, legal and illegal
transitions, atomic-write leftovers, run records and exit status, the run
refusal order (stop, lock, floor, criterion/candidate, allowed roots) and the
Revision 2 run preconditions (job state, this job's worker), stale runs after a
criterion or candidate change, the Ready-for-approval and Approved conditions
(worker or live run, whole-token quote, re-check), interrupted and
group-left-running runs, stop / check-stop / reconcile with real processes
(also from a different TZ, locale and PATH, and during the probe), the worker
registry, the queue lock (queue_dir ownership, extra files), and a truthful report.
Revision 9: the packet's SCREENING.txt (form, line numbers, included test, outside
access against the network grant, re-check at Approved and in the report) and the
manual-resume Waiting step.

The tool is only ever run through subprocess (-I -B); the test process never
imports it or reads its source. Every queue, control home, lock directory, space path,
allowed root and git repository lives under a temporary directory; --home and
--lock-dir are always passed, so the default home is never touched.
Each refusal is paired with a positive control (the same call without the
forbidden condition is accepted).
Revision 12: `run` needs macOS. Off macOS (real, or simulated inside the test
processes with GC_TEST_SIMULATE_NON_MACOS=1) the helper is run through the
test-only launcher tests/auto_queue_standin.py, whose stand-in for sandbox-exec
runs commands unrestricted; every command run through it is local. Only the
tests of the real macOS mechanisms (and their duplicates) are skipped there, each
with its reason (GUIDED_CODING_VERIFICATION.md lists them). The refusal tests
run the real helper (a subprocess that imports it and sets sys.platform), never
the launcher. `init` without --lock-dir is refused. The launcher cannot run a
real queue job (P8): it refuses (exit 16) any queue outside a temporary
directory of this test run that holds the test-only marker, which setUp writes.
P9: one test reads the helper's own process table (auto_queue.snapshot(), in a
subprocess that imports the helper) and checks that it lists a test process with
its parent and group, and that cleanup ends that process.
P10: the helper allows one second on the earlier bound of a run's window, because
Linux ps can report a start time up to 1 s early. On the Mac the tests imitate that
inside a test-side launcher (LINUX_PS_SCRIPT: the helper's ps start times are taken
as floor(exact fork time - 0.8 s) and its `ps -E` markers are dropped, as on Linux);
four tests use it always, and GC_TEST_IMITATE_LINUX_PS=1 (read only by the tests) runs
every helper call of a selected test through it. Off the Mac the launcher changes
nothing: real Linux ps already reports early.
"""

import calendar
import errno
import fcntl
import hashlib
import json
import math
import os
import shlex
import shutil
import socket
import struct
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path

TOOL = Path(__file__).resolve().parent.parent / "auto_queue.py"
# Revision 12: the test-only launcher, used only where the platform is not macOS
# (real, or simulated with GC_TEST_SIMULATE_NON_MACOS=1, read only by the tests).
STANDIN = Path(__file__).resolve().parent / "auto_queue_standin.py"
NON_MACOS = sys.platform != "darwin" or os.environ.get("GC_TEST_SIMULATE_NON_MACOS") == "1"
# Revision 12 (P8): the launcher runs only queues inside a directory directly under the
# temporary directory that holds this marker (written by setUp); otherwise exit 16.
STANDIN_MARKER = "GC_KIT_TEST_QUEUE_ROOT"
STANDIN_MARKER_TEXT = "guided_coding kit tests: test-only queue root for this test run; never a real queue\n"
STANDIN_REFUSED = 16
STANDIN_MESSAGE = ("REFUSED: tests/auto_queue_standin.py is a test-only launcher and runs only "
                   "queues the kit's own tests created; nothing was started or recorded.")
MACOS_ONLY = "macOS mechanism, skipped off macOS: "
PLATFORM_REFUSED = 15
PLATFORM_MESSAGE = "REFUSED: auto-mode runs need macOS (sandbox-exec and `ps -E`); this system is "
# The real helper with non-macOS forced inside its own process (no option or variable).
FORCED_NON_MACOS = ("import sys; sys.path.insert(0, sys.argv[1]); import auto_queue; "
                    "sys.platform = 'linux'; sys.exit(auto_queue.main(sys.argv[2:]))")
QID = "q1"
ALPHA, BETA = "job-alpha", "job-beta"
JOBS = [{"id": ALPHA, "title": "Alpha", "requires": ["unit", "build"], "depends_on": []},
        {"id": BETA, "title": "Beta", "requires": ["unit"], "depends_on": []}]
REPORT_HEADER = "| Job | Result | Decision or next step |"
WAITING_STEP = ("resume with /gc auto resume after the reset "
                "(automatic continuation is not supported)")
TIMEOUT = 15
# A CMD (or grandchild) that signals readiness and exits 0 once the release file
# appears (safety deadline 30 s; 120 s with the flag `long`, revision 14). mode is
# `+`-joined flags: `ignore`: SIGTERM ignored; `exit0`: SIGTERM exits 0 (a test
# script that traps TERM); `note`: SIGTERM writes <ready>.term and is otherwise
# ignored; `pgrp`: move to a new process group (same session); `setsid`: new
# session; `default`: SIGTERM kills it.
WAIT_SCRIPT = """import os, signal, sys, time
mode, ready, release = sys.argv[1:4]
flags = set(mode.split("+"))
if "ignore" in flags:
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
elif "exit0" in flags:
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
elif "note" in flags:
    signal.signal(signal.SIGTERM, lambda *_: open(ready + ".term", "w").close())
if "pgrp" in flags:
    os.setpgid(0, 0)
if "setsid" in flags:
    os.setsid()
open(ready, "w").close()
deadline = time.time() + (120 if "long" in flags else 30)
while not os.path.exists(release) and time.time() < deadline:
    time.sleep(0.05)
"""
# A CMD that starts children (inheriting its group), writes their pids, then
# exits 0 (`exit`) or stays until the deadline (`hold`).
LAUNCH_SCRIPT = """import os, subprocess, sys, time
end, wait, pids = sys.argv[1:4]
specs = [s.split(",") for s in sys.argv[4:]]
children = [subprocess.Popen([sys.executable, wait, *spec], stdin=subprocess.DEVNULL,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for spec in specs]
deadline = time.time() + 10
while not all(os.path.exists(s[1]) for s in specs) and time.time() < deadline:
    time.sleep(0.05)
with open(pids + ".tmp", "w") as handle:
    handle.write(" ".join(str(child.pid) for child in children))
os.replace(pids + ".tmp", pids)
if end == "hold":
    time.sleep(30)
"""
# Revision 12 (P9): print the helper's own process table, pid -> [ppid, pgid, sid, start],
# from auto_queue.snapshot() in this subprocess (the test process never imports the helper).
SNAPSHOT_SCRIPT = """import json, sys
sys.path.insert(0, sys.argv[1])
import auto_queue
print(json.dumps({str(pid): entry[:4] for pid, entry in auto_queue.snapshot().items()}))
"""
# Revision 13 (item 14): a test-side launcher that runs the helper (argv[2]) with its ps path
# replaced, inside this one process, by argv[3] (a missing path or a stand-in ps script); the
# installed helper is never changed. With argv[4] == "linux" the platform is also set to linux
# in that process (only for snapshot(), never for `run`).
FAKE_PS_SCRIPT = """import os, runpy, sys
kit, helper, ps = sys.argv[1:4]
sys.path.insert(0, kit)
import auto_queue
auto_queue.PS = ps
sys.argv = [helper] + sys.argv[4:]
if os.path.realpath(helper) == os.path.realpath(auto_queue.__file__):
    sys.exit(auto_queue.main(sys.argv[1:]))
runpy.run_path(helper, run_name="__main__")
"""
# Stand-ins for ps (Revision 13): refused as inside sandbox-exec, empty, malformed, and one that
# fails only for the `ps -E` marker query (as Linux ps has no -E).
# Revision 15: Python scripts run by this test's own interpreter (no /bin/sh). Each takes one
# argument word (-E or -ax) as its trigger, as the earlier `case " $* "` did, prints the same
# lines with the same exit status, and otherwise replaces itself with the real /bin/ps.
FAKE_PS_HEAD = "#!%s\nimport os, sys\nargs = sys.argv[1:]\n" % sys.executable
FAKE_PS_PASS = 'os.execv("/bin/ps", ["/bin/ps"] + args)\n'
FAKE_PS = {
    "refused": FAKE_PS_HEAD + 'print("ps: Operation not permitted", file=sys.stderr)\nsys.exit(71)\n',
    "empty": FAKE_PS_HEAD + "sys.exit(0)\n",
    "malformed": FAKE_PS_HEAD + 'print("PID garbage line")\nprint("12 x")\nsys.exit(0)\n',
    "no-marker": (FAKE_PS_HEAD + 'if "-E" in args:\n'
                  '    print("ps: unsupported option -E", file=sys.stderr)\n    sys.exit(1)\n'
                  + FAKE_PS_PASS),
    # Checkpoint 3 (REVIEW_CP2 findings 1 and 3): the whole-table query fails while
    # `ps -p` works; the `ps -E` marker query exits 0 with nothing, or with garbage.
    "ax-fails": (FAKE_PS_HEAD + 'if "-ax" in args:\n'
                 '    print("ps: Operation not permitted", file=sys.stderr)\n    sys.exit(71)\n'
                 + FAKE_PS_PASS),
    "empty-marker": FAKE_PS_HEAD + 'if "-E" in args:\n    sys.exit(0)\n' + FAKE_PS_PASS,
    "garbage-marker": (FAKE_PS_HEAD + 'if "-E" in args:\n    print("garbage")\n    sys.exit(0)\n'
                       + FAKE_PS_PASS),
}
# snapshot() in a subprocess that imports the helper with PS = argv[2] (and, with argv[3]
# == "linux", sys.platform = "linux"): prints `table N` (its size) or `refused CODE MESSAGE`.
SNAPSHOT_PROBE = """import sys
sys.path.insert(0, sys.argv[1])
import auto_queue
auto_queue.PS = sys.argv[2]
if sys.argv[3] == "linux":
    sys.platform = "linux"
try:
    procs = auto_queue.snapshot()
except auto_queue.Refused as refusal:
    print("refused", refusal.code, refusal)
else:
    print("table", len(procs), int(any(procs[p][4] for p in procs)))
"""
# Revision 14: the helper's own view of each pid in argv[2:] (imported in this subprocess; the
# test process never imports it): pid -> [process_exists, group_alive (pid as a group id),
# listed by snapshot()].
# The third view is "listed as a live process": in the table and (revision 14) not flagged as a
# zombie there.
ZOMBIE_PROBE = """import json, sys
sys.path.insert(0, sys.argv[1])
import auto_queue
procs = auto_queue.snapshot()
def listed(pid):
    entry = procs.get(pid)
    return entry is not None and not (len(entry) > 5 and entry[5])
print(json.dumps({p: [auto_queue.process_exists(int(p)), auto_queue.group_alive(int(p)),
                      listed(int(p))] for p in sys.argv[2:]}))
"""
# Revision 14 (review finding 2): survivors() and leaders_of() of the run record in the JSON file
# argv[2], in the helper's own process.
SURVIVORS_PROBE = """import json, sys
sys.path.insert(0, sys.argv[1])
import auto_queue
with open(sys.argv[2]) as handle:
    record = json.load(handle)
print(json.dumps({"survivors": auto_queue.survivors(record),
                  "leaders": sorted(auto_queue.leaders_of(record, auto_queue.snapshot()))}))
"""
# Revision 14 (review finding 1): a live process whose main thread has exited while a worker
# thread runs. On Linux its leader task then shows state Z (ps `Zl`) until the whole process
# exits. The main thread ends with the raw exit system call (not exit_group), made through
# ctypes, which releases the interpreter lock for the call; the worker writes a heartbeat
# file every 0.1 s until the release file appears (at most 60 s), then exits the process.
# SIGTERM is ignored, so a stop's signal leaves it running.
MT_SCRIPT = """import ctypes, os, signal, sys, threading, time
signal.signal(signal.SIGTERM, signal.SIG_IGN)
beat, release = sys.argv[1:3]
SYS_EXIT = {"x86_64": 60, "aarch64": 93}[os.uname().machine]
def worker():
    deadline, n = time.time() + 60, 0
    while not os.path.exists(release) and time.time() < deadline:
        n += 1
        with open(beat + ".tmp", "w") as handle:
            handle.write(str(n))
        os.replace(beat + ".tmp", beat)
        time.sleep(0.1)
    os._exit(0)
threading.Thread(target=worker).start()
time.sleep(0.3)
ctypes.CDLL(None, use_errno=True).syscall(SYS_EXIT, 0)
"""
# A fingerprint script for `install check`.
FINGERPRINT_SCRIPT = """#!{python}
import sys
print("FINGERPRINT {fingerprint}")
sys.exit({code})
"""
# The libtbx.easy_run pattern: CMD starts a shell in a new session, waits 2 s, exits.
# Revision 15: the "shell" is a Python process running SHELL_STANDIN (the file argv[3]), started
# with setsid as easy_run starts /bin/sh, with the JSON spec argv[2].
EASY_RUN_SCRIPT = """import os, subprocess, sys, time
shell = subprocess.Popen([sys.executable, sys.argv[3], sys.argv[2]], preexec_fn=os.setsid)
with open(sys.argv[1] + ".tmp", "w") as handle:
    handle.write(str(shell.pid))
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
time.sleep(2)
"""
# Revision 15 (review finding 1): a platform /bin/sh (and /bin/sleep) shows no environment to
# `ps -E` on macOS, so the helper cannot see a run-id marker on it and must find it through the
# wrapper's tracked descendants. The stand-in therefore first re-executes itself (same pid,
# start time and parent) with the marker moved to another variable, so that its own visible
# environment has no marker; each command runs as its child with the marker restored, except a
# command flagged "platform" (standing in for a platform binary such as /bin/sleep), which runs
# without it. Review finding 3: with "exec_last" the last command replaces the stand-in
# (os.execve), as bash 5.1 and later and dash do with the last simple command of `sh -c 'a; b'`;
# otherwise it runs as a child and the stand-in waits, as macOS /bin/sh (bash 3.2) and the bash
# 4.2 of cci-gpu-00 do. The exit status is the last command's.
SHELL_STANDIN = """import json, os, subprocess, sys
MARKER, KEPT = "GC_AUTO_RUN_ID", "GC_TEST_STANDIN_KEPT_RUN_ID"
if MARKER in os.environ:
    env = dict(os.environ)
    env[KEPT] = env.pop(MARKER)
    os.execve(sys.executable, [sys.executable] + sys.argv, env)
spec = json.loads(sys.argv[1])


def environment(platform):
    env = {key: value for key, value in os.environ.items() if key != KEPT}
    if KEPT in os.environ and not platform:
        env[MARKER] = os.environ[KEPT]
    return env


code, commands = 0, spec["commands"]
for index, command in enumerate(commands):
    env = environment(command.get("platform", False))
    if spec.get("exec_last") and index == len(commands) - 1:
        os.execve(command["argv"][0], command["argv"], env)
    code = subprocess.call(command["argv"], env=env)
sys.exit(code)
"""
# Revision 8 CMD: a background child (`plain`, or `ignore` SIGTERM) and a 60 s sleep.
TIMEOUT_SCRIPT = """import os, subprocess, sys, time
pidfile, mode = sys.argv[1:3]
body = "import signal, time\\n"
if mode == "ignore":
    body += "signal.signal(signal.SIGTERM, signal.SIG_IGN)\\n"
child = subprocess.Popen([sys.executable, "-c", body + "time.sleep(60)"], stdin=subprocess.DEVNULL,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
with open(pidfile + ".tmp", "w") as handle:
    handle.write(str(child.pid))
os.replace(pidfile + ".tmp", pidfile)
time.sleep(60)
"""
# An HTTPS request to argv[2]; prints argv[1] (a probe path) on success unless it is "-".
# With argv[3] the script also writes its own error text to that file (re-read 10 (a)).
NET_SCRIPT = """import sys, urllib.request
try:
    status = urllib.request.urlopen(sys.argv[2], timeout=15).status
except Exception as error:
    print("network error:", repr(error), file=sys.stderr)
    if len(sys.argv) > 3:
        with open(sys.argv[3], "w") as handle:
            handle.write("network error: " + repr(error) + "\\n")
    sys.exit(1)
print("status", status, file=sys.stderr)
if sys.argv[1] != "-":
    print(sys.argv[1])
"""
# A localhost TCP round trip to the test's listener.
LOCAL_SCRIPT = """import socket, sys
with socket.create_connection(("127.0.0.1", int(sys.argv[1])), timeout=10) as connection:
    connection.sendall(b"hi")
    sys.exit(0 if connection.recv(16) == b"ok\\n" else 1)
"""
# A TCP connection to a literal address argv[1], port 443 (no name lookup involved).
CONNECT_SCRIPT = """import socket, sys
try:
    socket.create_connection((sys.argv[1], 443), timeout=10).close()
except OSError as error:
    print("connect failed:", repr(error), file=sys.stderr)
    sys.exit(1)
print("connected")
"""
# A probe that detaches a 25 s python into its own session with output inherited,
# then runs a 30 s child (revision 15: a python child in its group instead of /bin/sleep 30).
DETACH_PROBE = """import os, subprocess, sys
child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(25)"],
                         start_new_session=True)
with open(sys.argv[1] + ".tmp", "w") as handle:
    handle.write(str(child.pid))
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
subprocess.run([sys.executable, "-c", "import time; time.sleep(30)"])
"""


def sleep_cmd(seconds):
    """Revision 15: a CMD that stays for `seconds` (instead of `sleep N`)."""
    return [sys.executable, "-c", "import time; time.sleep(%s)" % seconds]


# Revision 15: a CMD that exits 0 at once (instead of `true`).
QUICK_CMD = [sys.executable, "-c", "pass"]
# Revision 12 (P10): a test-side launcher that runs the helper (argv[2]: auto_queue.py, or the
# stand-in launcher, which then imports this already patched module) with Linux procps start
# times imitated inside this one process: every ps start time the helper reads (snapshot() and
# process_start()) becomes floor(exact fork time - LAG), the way procps adds clock ticks to a
# boot time truncated to whole seconds, and the `ps -E` markers are dropped (Linux ps has no
# -E). The exact fork time comes from the macOS kernel (sysctl KERN_PROC_PID, p_starttime).
# Off macOS it changes nothing. The installed helper is never changed.
IMITATE_LINUX_PS = os.environ.get("GC_TEST_IMITATE_LINUX_PS") == "1"
LINUX_PS_LAG = 0.8
LINUX_PS_SCRIPT = """import ctypes, math, os, runpy, sys, time
kit, helper, lag = sys.argv[1], sys.argv[2], float(sys.argv[3])
sys.path.insert(0, kit)
import auto_queue
if sys.platform == "darwin":
    libc = ctypes.CDLL(None)
    real_snapshot, real_process_start = auto_queue.snapshot, auto_queue.process_start

    def lagged(pid):
        mib = (ctypes.c_int * 4)(1, 14, 1, int(pid))  # CTL_KERN, KERN_PROC, KERN_PROC_PID
        buf, size = ctypes.create_string_buffer(1024), ctypes.c_size_t(1024)
        if libc.sysctl(mib, 4, buf, ctypes.byref(size), None, 0) != 0 or size.value < 16:
            return None
        sec = int.from_bytes(buf.raw[0:8], "little", signed=True)
        usec = int.from_bytes(buf.raw[8:12], "little", signed=True)
        return int(math.floor(sec + usec / 1e6 - lag)) if sec > 0 else None

    def snapshot():
        # Revision 14: any further fields (the zombie flag) are kept as they are.
        return {pid: [ppid, pgid, sid, None if start is None else lagged(pid), set(), *rest]
                for pid, (ppid, pgid, sid, start, marks, *rest) in real_snapshot().items()}

    def process_start(pid):
        text = real_process_start(pid)
        if text in (None, auto_queue.UNREADABLE):
            return text
        start = lagged(pid)
        return text if start is None else time.strftime("%a %b %d %H:%M:%S %Y", time.gmtime(start))

    auto_queue.snapshot, auto_queue.process_start = snapshot, process_start
sys.argv = [helper] + sys.argv[4:]
if os.path.realpath(helper) == os.path.realpath(auto_queue.__file__):
    sys.exit(auto_queue.main(sys.argv[1:]))
runpy.run_path(helper, run_name="__main__")
"""
# A name lookup; exit 0 when the host resolves.
DNS_SCRIPT = """import socket, sys
try:
    socket.getaddrinfo(sys.argv[1], 80)
except OSError as error:
    print("lookup failed:", repr(error), file=sys.stderr)
    sys.exit(1)
print("resolved", sys.argv[1])
"""
# Revision 9 (Q1): the screening record every ready packet carries by default.
SCREENING = "test: unit | outside: none | included\n"
# Revision 10: the default suite contacts no outside host. Blocked-path tests use
# 192.0.2.1 (RFC 5737 documentation range, no host); tests that would contact a
# real outside host or name run only with GC_TEST_NETWORK=1.
BLOCKED_URL, BLOCKED_ADDRESS = "https://192.0.2.1/", "192.0.2.1"
OUTSIDE_URL, OUTSIDE_ADDRESS, OUTSIDE_NAME = "https://www.rcsb.org", "1.1.1.1", "example.com"
OPT_IN_PREFIX = "outside network test skipped (set GC_TEST_NETWORK=1 to run): contacts "
# Re-read 10 (b): removed (either case) from the environment of blocked-path runs.
PROXY_VARIABLES = ("HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "NO_PROXY")
SCREEN_CHECK = Path(__file__).resolve().parent.parent / "payload" / "tools" / "screen_check.py"
# A CMD that reports what it sees: run-record statuses, cwd, one stderr line.
PEEK_SCRIPT = """import glob, json, os, sys
for path in glob.glob(os.path.join(sys.argv[1], "*.json")):
    with open(path) as handle:
        print("seen", json.load(handle).get("status"))
print("cwd", os.getcwd())
print("stderr-line", file=sys.stderr)
sys.exit(int(sys.argv[2]))
"""


# Revision 13 (item 12): the kit's own suite needs a working /bin/ps. Inside sandbox-exec (for
# example a queue run's outer sandbox) the setuid ps is refused; then this module runs ONE
# failing test naming the cause and the supported route instead of many unrelated failures.
PS_ROUTE = ("run the kit's own suite attended, outside the queue's outer sandbox (sandbox-exec), "
            "from an ordinary shell; see AUTO.md, 'Running the kit's own tests'. A queue job "
            "that must run it needs that exception stated in advance; never disable the "
            "network protection to make it run.")


def ps_problem(ps="/bin/ps"):
    """None when `ps -p <this pid>` works, else a message naming the cause and the route."""
    try:
        out = subprocess.run([ps, "-p", str(os.getpid()), "-o", "pid="], capture_output=True,
                             text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as error:
        cause = "%s could not be run: %s" % (ps, error)
    else:
        if out.returncode == 0 and out.stdout.split() == [str(os.getpid())]:
            return None
        cause = "%s exited %d%s" % (ps, out.returncode, (
            " (%s)" % out.stderr.strip().splitlines()[0]) if out.stderr.strip() else
            " and did not list this process")
    return ("PROCESS TABLE UNREADABLE: the auto-queue tests cannot run here: %s. This happens "
            "inside sandbox-exec, where the setuid ps is refused. Supported route: %s"
            % (cause, PS_ROUTE))


class KitSelfTestRoute(unittest.TestCase):
    """Revision 13 (item 12): the one test that stands for the whole module when ps cannot run."""

    def test_process_table_readable(self):
        problem = ps_problem()
        if problem:
            self.fail(problem)

    def test_preflight_fires_on_refused_ps_and_is_silent_normally(self):
        """A ps stand-in refused like ps inside sandbox-exec (exit 71), one that prints
        nothing, and a missing ps each give the message with the cause and the route;
        the real ps gives None."""
        self.assertIsNone(ps_problem())
        with tempfile.TemporaryDirectory() as temp:
            for mode in ("refused", "empty"):
                fake = Path(os.path.realpath(temp)) / mode
                fake.write_text(FAKE_PS[mode])
                fake.chmod(0o755)
                problem = ps_problem(str(fake))
                self.assertIsNotNone(problem, mode)
                self.assertTrue(problem.startswith("PROCESS TABLE UNREADABLE"), problem)
                self.assertIn("sandbox-exec", problem)
                self.assertIn("Running the kit's own tests", problem)
            self.assertIn("exited 71 (ps: Operation not permitted)",
                          ps_problem(str(Path(temp) / "refused")))
            self.assertIn("could not be run", ps_problem(str(Path(temp) / "missing")))

    @unittest.skipUnless(sys.platform == "darwin" and os.path.exists("/usr/bin/sandbox-exec"),
                         MACOS_ONLY + "the real sandbox-exec refusal of ps")
    def test_module_under_sandbox_exec_fails_once(self):
        """Under `sandbox-exec -p '(version 1)(allow default)'` this module runs exactly one
        test, which fails with the message; nothing else runs."""
        kit = Path(__file__).resolve().parent.parent
        result = subprocess.run(["/usr/bin/sandbox-exec", "-p", "(version 1)(allow default)",
                                 sys.executable, "-I", "-B", "-m", "unittest", "discover", "-s",
                                 "tests", "-p", Path(__file__).name], cwd=kit, text=True,
                                capture_output=True, timeout=300)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("Ran 1 test", result.stderr)
        self.assertIn("FAILED (failures=1)", result.stderr)
        self.assertEqual(result.stderr.count("PROCESS TABLE UNREADABLE"), 1, result.stderr)
        self.assertIn("Supported route", result.stderr)


def load_tests(loader, standard_tests, pattern):
    """Revision 13 (item 12): when ps cannot run, only the route test (failing once)."""
    if ps_problem() is None:
        return standard_tests
    return unittest.TestSuite([KitSelfTestRoute("test_process_table_readable")])


def ps_state(pid):
    """Revision 14: the ps stat of one process ("" when ps lists nothing for it)."""
    return subprocess.run(["/bin/ps", "-o", "stat=", "-p", str(pid)], text=True,
                          capture_output=True).stdout.strip()


def proc_stat(pid):
    """Revision 14, Linux: (state, num_threads) from /proc/<pid>/stat, or None."""
    try:
        data = Path("/proc/%d/stat" % pid).read_bytes()
        fields = data[data.rindex(b")") + 1:].split()
        return fields[0].decode(), int(fields[17])
    except (OSError, ValueError, IndexError):
        return None


def alive(pid):
    """Revision 14: as in the helper, a zombie (exited, not yet collected by its parent, as
    under a container's process 1 that never reaps) counts as gone: on Linux /proc state Z
    with one thread (a Z leader with running threads is alive); elsewhere ps stat Z without
    the `l` flag. Anything unread counts as alive."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    if sys.platform.startswith("linux"):
        reading = proc_stat(pid)
        return not (reading is not None and reading[0] == "Z" and reading[1] == 1)
    state = ps_state(pid)
    return not (state.startswith("Z") and "l" not in state[1:])


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def job_entry(data, job):
    """The job's entry in status output (a list of objects with `id`, or a map by id)."""
    jobs = data["jobs"]
    if isinstance(jobs, dict):
        return jobs[job]
    matches = [entry for entry in jobs if entry.get("id") == job]
    assert len(matches) == 1, (job, jobs)
    return matches[0]


# ---- Revision 15: archive and Mach-O fixtures written by the tests (no ar, ranlib or cc) ----
AR_DATE, AR_REDATE = 1000000000, 1100000000


def ar_header(field, when, size, uid=0):
    """One 60-byte `!<arch>` member header; `when` None leaves the date, owner, group and mode
    fields blank, as GNU ar writes its `//` name table."""
    if when is None:
        return b"%-16s%-12s%-6s%-6s%-8s%-10d`\n" % (field, b"", b"", b"", b"", size)
    return b"%-16s%-12d%-6d%-6d%-8s%-10d`\n" % (field, when, uid, uid, b"100644", size)


def ar_archive(members, when, style="gnu", index_when=None, uid=0):
    """Revision 15: an `!<arch>` archive of `members` (name -> bytes) written by the test,
    every member header dated `when` and the index member dated `index_when` (default
    `when`; macOS ranlib dates its index with the current time), in the two real layouts:
    - "gnu" (GNU and SysV ar): a `/` symbol table first (big-endian count, member offsets,
      NUL-terminated names), a `//` name table for names longer than 15 bytes (each name
      followed by `/` and a newline, the table padded to an even size with a newline inside
      its size, as GNU ar writes it; those members are then named `/N`), other names as
      `name/`; a member whose size is odd is
      followed by one newline outside its size;
    - "bsd" (BSD and macOS ar with ranlib): a `__.SYMDEF SORTED` index first (`#1/20`;
      ranlib entries of string index and member offset, then the strings, NUL-padded to a
      multiple of 8), every name in `#1/N` form (the name NUL-padded to N = 4 more than a
      multiple of 8, at least one NUL, at the start of the member data and counted in its
      size), member data padded with newlines to a multiple of 8 inside the size, so that every
      member's size is 4 more than a multiple of 8 and every header starts 8-aligned: as Apple
      ar and ranlib (Xcode 26) write it (review finding 5; compared in the record).
    The index names one symbol per member (`sym_` and the member name's stem)."""
    index_when = when if index_when is None else index_when
    parts, longnames = [], b""  # (owner member or None, name field, stored bytes, date)
    for name, data in members.items():
        raw = name.encode()
        if style == "bsd":
            raw = raw.ljust((len(raw) + 1 + 7) // 8 * 8 + 4, b"\0")
            body = raw + data + b"\n" * (-len(data) % 8)
            parts.append((name, b"#1/%d" % len(raw), body, when))
        elif len(raw) > 15:
            parts.append((name, b"/%d" % len(longnames), data, when))
            longnames += raw + b"/\n"
        else:
            parts.append((name, raw + b"/", data, when))
    if longnames:
        # GNU ar pads the name table to an even size with a newline inside its size (found by
        # test_real_tool_conformance_archive against GNU ar 2.27 on cci-gpu-00).
        parts.insert(0, (None, b"//", longnames + b"\n" * (len(longnames) % 2), None))
    symbols = [("sym_" + name.split(".")[0], name) for name in members]

    def index(offsets):
        if style == "bsd":
            ranlibs, strings = b"", b""
            for symbol, member in symbols:
                ranlibs += struct.pack("<II", len(strings), offsets[member])
                strings += symbol.encode() + b"\0"
            strings += b"\0" * (-len(strings) % 8)
            raw = b"__.SYMDEF SORTED\0\0\0\0"
            return b"#1/%d" % len(raw), (raw + struct.pack("<I", len(ranlibs)) + ranlibs
                                         + struct.pack("<I", len(strings)) + strings)
        return b"/", (struct.pack(">I", len(symbols))
                      + b"".join(struct.pack(">I", offsets[member]) for _, member in symbols)
                      + b"".join(symbol.encode() + b"\0" for symbol, _ in symbols))

    def stored(body):
        return 60 + len(body) + len(body) % 2

    pos, offsets = 8 + stored(index({member: 0 for _, member in symbols})[1]), {}
    for owner, _, body, _ in parts:
        if owner:
            offsets[owner] = pos
        pos += stored(body)
    field, body = index(offsets)
    out = [b"!<arch>\n", ar_header(field, index_when, len(body), uid), body, b"\n" * (len(body) % 2)]
    for _, field, body, date in parts:
        out += [ar_header(field, date, len(body), uid), body, b"\n" * (len(body) % 2)]
    return b"".join(out)


def ar_corruptions(data):
    """Revision 15: damaged forms of a parsable archive (name -> bytes): cut inside its last
    member, a member header shorter than 60 bytes, bad member header magic, a member size
    past the end, a size that is not a number, bad archive magic."""
    size = data[8 + 48:8 + 58]
    return {
        "truncated archive": data[:-5],
        "truncated member header": data[:8 + 30],
        "bad header magic": data[:8 + 58] + b"!!" + data[8 + 60:],
        "size past end": data[:8 + 48] + b"%-10d" % len(data) + data[8 + 58:],
        "size not a number": data[:8 + 48] + size.replace(b" ", b"x", 1) + data[8 + 58:],
        "bad archive magic": b"!<arcx>\n" + data[8:],
    }


MACHO_DATE = "Sun Oct 11 08:54:45 2026"
MACHO_SIGNATURE_SIZE = 256


def macho_bundle(date=MACHO_DATE, extra="phaser", value=1):
    """Revision 15: the bytes of a 64-bit little-endian Mach-O bundle written by the test (no
    cc). Review finding 2: it carries the load commands of `cc -bundle` output on arm64 macOS
    (Apple clang 17), in that order, with the same command sizes (848 bytes in all): __TEXT (the
    header, the load commands, then __text with the code `mov w0, #value; ret`, __cstring with
    the asctime build date and `extra`, and __unwind_info), __DATA, __LINKEDIT,
    LC_DYLD_CHAINED_FIXUPS, LC_DYLD_EXPORTS_TRIE, LC_SYMTAB, LC_DYSYMTAB, LC_UUID,
    LC_BUILD_VERSION, LC_SOURCE_VERSION, LC_LOAD_DYLIB (/usr/lib/libSystem.B.dylib),
    LC_FUNCTION_STARTS, LC_DATA_IN_CODE and LC_CODE_SIGNATURE; the data of the fixups, the
    exports trie, the function starts, the symbol and string tables lie in __LINKEDIT before
    the signature range, which ends at the end of __LINKEDIT and of the file. As the linker
    does, the UUID and the signature range are derived from the file's other bytes (here:
    SHA-256 of the file with both zeroed; the signature range is a block of derived bytes
    standing in for the linker's ad-hoc signature, not a real code signature), so both change
    whenever the date, `extra` or `value` does. The contents of the fixups, trie, function
    starts and tables are fixed bytes of the right sizes, not a working dynamic library."""
    page, text_at, cstring_at, linkedit_at = 0x4000, 0x3f80, 0x3f88, 0x8000
    code = struct.pack("<II", 0x52800000 | (value << 5), 0xd65f03c0)
    cstring = date.encode() + b"\0" + extra.encode() + b"\0"
    unwind_at = (cstring_at + len(cstring) + 3) // 4 * 4
    unwind = bytes(range(88))
    fixups, trie, starts = bytes(range(1, 81)), bytes(range(101, 157)), b"\x80\x7f" + b"\0" * 6
    strings = b" \0_build_time\0_version\0_what\0".ljust(32, b"\0")
    symbols = b"".join(struct.pack("<IBBHQ", strx, 0x0f, 1 + (n == 2), 0, 0x3f80 + 8 * n)
                       for n, strx in enumerate((2, 14, 23)))
    at = {"fixups": linkedit_at, "trie": linkedit_at + 80, "starts": linkedit_at + 136,
          "symtab": linkedit_at + 144, "strtab": linkedit_at + 192, "signature": linkedit_at + 224}
    end = at["signature"] + MACHO_SIGNATURE_SIZE

    def segment(name, fileoff, filesize, sections):
        head = struct.pack("<II16sQQQQiiII", 0x19, 72 + 80 * len(sections), name,
                           0x100000000 + fileoff, filesize, fileoff, filesize, 7, 5,
                           len(sections), 0)
        return head + b"".join(struct.pack("<16s16sQQIIIIIIII", sect, name,
                                           0x100000000 + offset, size, offset, 2, 0, 0, 0,
                                           0, 0, 0) for sect, offset, size in sections)

    def linkedit_data(cmd, offset, size):
        return struct.pack("<IIII", cmd, 16, offset, size)
    dylib = b"/usr/lib/libSystem.B.dylib\0".ljust(32, b"\0")
    commands = [segment(b"__TEXT", 0, page, [(b"__text", text_at, len(code)),
                                             (b"__cstring", cstring_at, len(cstring)),
                                             (b"__unwind_info", unwind_at, len(unwind))]),
                segment(b"__DATA", page, page, [(b"__data", page, 16)]),
                segment(b"__LINKEDIT", linkedit_at, end - linkedit_at, []),
                linkedit_data(0x80000034, at["fixups"], len(fixups)),  # LC_DYLD_CHAINED_FIXUPS
                linkedit_data(0x80000033, at["trie"], len(trie)),  # LC_DYLD_EXPORTS_TRIE
                struct.pack("<IIIIII", 0x2, 24, at["symtab"], 3, at["strtab"], len(strings)),
                struct.pack("<II18I", 0xb, 80, 0, 0, 0, 3, 3, *([0] * 13)),  # LC_DYSYMTAB
                struct.pack("<II16s", 0x1b, 24, b"\0" * 16),  # LC_UUID
                struct.pack("<IIIIIIII", 0x32, 32, 1, 26 << 16, 26 << 16 | 1 << 8, 1, 3,
                            1230 << 16 | 1 << 8),  # LC_BUILD_VERSION
                struct.pack("<IIQ", 0x2a, 16, 0),  # LC_SOURCE_VERSION
                struct.pack("<IIIIII", 0xc, 56, 24, 2, 1356 << 16, 1 << 16) + dylib,  # LC_LOAD_DYLIB
                linkedit_data(0x26, at["starts"], len(starts)),  # LC_FUNCTION_STARTS
                linkedit_data(0x29, at["symtab"], 0),  # LC_DATA_IN_CODE
                linkedit_data(0x1d, at["signature"], MACHO_SIGNATURE_SIZE)]  # LC_CODE_SIGNATURE
    loads = b"".join(commands)
    data = bytearray(end)
    data[0:32] = struct.pack("<IiiIIIII", 0xfeedfacf, 0x0100000c, 0, 8, len(commands),
                             len(loads), 0x85, 0)
    data[32:32 + len(loads)] = loads
    data[text_at:text_at + len(code)] = code
    data[cstring_at:cstring_at + len(cstring)] = cstring
    data[unwind_at:unwind_at + len(unwind)] = unwind
    data[page:page + 16] = struct.pack("<QQ", 0x100000000 + cstring_at,
                                       0x100000000 + cstring_at + len(date) + 1)
    for key, blob in (("fixups", fixups), ("trie", trie), ("starts", starts),
                      ("symtab", symbols), ("strtab", strings)):
        data[at[key]:at[key] + len(blob)] = blob
    uuid_at = 32 + sum(len(command) for command in commands[:7]) + 8
    derived = hashlib.sha256(bytes(data)).digest()
    data[uuid_at:uuid_at + 16] = derived[:16]
    block = b"GC-TEST-STANDIN-SIGNATURE\0"
    while len(block) < MACHO_SIGNATURE_SIZE:
        block += hashlib.sha256(derived + len(block).to_bytes(4, "little")).digest()
    data[at["signature"]:end] = block[:MACHO_SIGNATURE_SIZE]
    return bytes(data)


# Revision 15 (review finding 2): the load command types of `cc -bundle` output (Apple clang 17,
# arm64, macOS 26 SDK), in order, as macho_bundle() writes them.
MACHO_CC_COMMANDS = [0x19, 0x19, 0x19, 0x80000034, 0x80000033, 0x2, 0xb, 0x1b, 0x32, 0x2a, 0xc,
                     0x26, 0x29, 0x1d]
MACHO_LINKEDIT_DATA = {0x80000034: "fixups", 0x80000033: "exports trie", 0x26: "function starts"}
MACHO_OTHER_COMMANDS = {0x32: "build version command", 0xc: "load dylib command",
                        0xb: "dysymtab command"}


def macho_load_commands(data):
    """Revision 15: the load command types of a 64-bit Mach-O file, in order."""
    ncmds = struct.unpack_from("<I", data, 16)[0]
    types, offset = [], 32
    for _ in range(ncmds):
        cmd, size = struct.unpack_from("<II", data, offset)
        types.append(cmd)
        offset += size
    return types


def macho_ranges(data):
    """Revision 15: named byte ranges of a 64-bit Mach-O file, read from its load commands:
    uuid, signature, the sections (text, cstring, unwind_info, data), linkedit (before the
    signature), the __LINKEDIT data of the fixups, exports trie, function starts, symbol table
    and string table, some other load commands, and the whole load command area (commands)."""
    ncmds, sizeofcmds = struct.unpack_from("<II", data, 16)
    ranges, offset = {"commands": (32, 32 + sizeofcmds)}, 32
    for _ in range(ncmds):
        cmd, size = struct.unpack_from("<II", data, offset)
        if cmd == 0x19:
            name = data[offset + 8:offset + 24].rstrip(b"\0")
            fileoff, filesize = struct.unpack_from("<QQ", data, offset + 40)
            for index in range(struct.unpack_from("<I", data, offset + 64)[0]):
                at = offset + 72 + 80 * index
                sect = data[at:at + 16].rstrip(b"\0").decode()
                sect_size, sect_offset = struct.unpack_from("<QI", data, at + 40)
                ranges[sect.strip("_")] = (sect_offset, sect_offset + sect_size)
            if name == b"__LINKEDIT":
                ranges["linkedit"] = (fileoff, None)
        elif cmd == 0x1b:
            ranges["uuid"] = (offset + 8, offset + 24)
        elif cmd == 0x1d:
            dataoff, datasize = struct.unpack_from("<II", data, offset + 8)
            ranges["signature"] = (dataoff, dataoff + datasize)
        elif cmd in MACHO_LINKEDIT_DATA:
            dataoff, datasize = struct.unpack_from("<II", data, offset + 8)
            if datasize:
                ranges[MACHO_LINKEDIT_DATA[cmd]] = (dataoff, dataoff + datasize)
        elif cmd == 0x2:
            symoff, nsyms, stroff, strsize = struct.unpack_from("<IIII", data, offset + 8)
            ranges["symbol table"] = (symoff, symoff + 16 * nsyms)
            ranges["string table"] = (stroff, stroff + strsize)
        if cmd in MACHO_OTHER_COMMANDS:
            ranges[MACHO_OTHER_COMMANDS[cmd]] = (offset + 8, offset + size)
        offset += size
    ranges["linkedit"] = (ranges["linkedit"][0], ranges["signature"][0])
    return ranges


class AutoQueueChecks(unittest.TestCase):
    linux_ps = IMITATE_LINUX_PS  # Revision 12 (P10); the P10 tests set it for themselves
    fake_ps = None  # Revision 13: a ps stand-in path for the helper (FAKE_PS_SCRIPT)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.tmp = Path(os.path.realpath(self.temp.name))
        (self.tmp / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)  # Revision 12 (P8)
        self.procs, self.stray, self.workers, self.ready_packets = [], [], {}, {}
        self.addCleanup(self.kill_leftovers)
        user = self.tmp / "user"
        user.mkdir()
        (user / ".gitconfig").write_text("")
        self.home = self.tmp / "home"
        self.home.mkdir()
        self.env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "GC_AUTO"))}
        self.env.update(HOME=str(user), GC_AUTO_HOME=str(self.home), GIT_CONFIG_NOSYSTEM="1",
                        GIT_CONFIG_GLOBAL=str(user / ".gitconfig"))
        self.records = self.tmp / "records"
        self.records.mkdir()
        self.lock = self.tmp / "RUN.lock"
        self.src = self.tmp / "allowed" / "src"
        (self.src / "build").mkdir(parents=True)
        self.outside = self.tmp / "outside"
        self.outside.mkdir()
        self.work = self.tmp / "work"
        self.work.mkdir()
        self.wait_script = self.work / "wait.py"
        self.wait_script.write_text(WAIT_SCRIPT)
        self.linux_ps_script = self.work / "linux_ps.py"  # Revision 12 (P10)
        self.linux_ps_script.write_text(LINUX_PS_SCRIPT)
        self.counter = 0
        self.grant = self.tmp / "GRANT.md"
        self.grant.write_text("# Grant\nSynthetic grant for tests.\n")
        self.jobs = self.tmp / "JOBS.json"
        self.jobs.write_text(json.dumps(JOBS))
        self.criterion = self.tmp / "criterion.txt"
        self.criterion.write_text("Criterion one: the tool exits 0.\n")
        self.repo = self.tmp / "repo"
        self.repo.mkdir()
        self.git("init", "-q")
        self.git("config", "user.name", "Test Writer")
        self.git("config", "user.email", "test@example.invalid")
        self.git("config", "commit.gpgsign", "false")
        self.commit("first")
        self.queue = self.init_queue(QID, self.home, self.lock)

    # ---- helpers -------------------------------------------------------------

    def kill_leftovers(self):
        for proc in self.procs:
            if proc.poll() is None:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    proc.kill()
                proc.wait(5)
        for pid in self.stray:
            try:
                os.kill(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
        for path in self.tmp.glob("records*/*/runs/*.json"):
            try:
                record = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            for pgid in (record.get("pgid"), record.get("probe_pgid")):
                if isinstance(pgid, int) and pgid > 1 and pgid != os.getpgid(0):
                    try:
                        os.killpg(pgid, signal.SIGKILL)
                    except (ProcessLookupError, PermissionError):
                        pass

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.repo), *args], env=self.env, check=True,
                              text=True, capture_output=True).stdout.strip()

    def commit(self, text):
        (self.repo / "file.txt").write_text(text + "\n")
        self.git("add", "file.txt")
        self.git("commit", "-q", "-m", text)

    def tree(self):
        return self.git("rev-parse", "HEAD^{tree}")

    def command(self, args, home, queue, real=False):
        """Revision 12: off macOS through the test-only launcher, unless `real`."""
        helper = STANDIN if NON_MACOS and not real else TOOL
        command = [sys.executable, "-I", "-B", str(helper), "--home", str(home or self.home)]
        if self.fake_ps and not real:  # Revision 13: the helper's ps replaced
            script = self.work / "fake_ps_launcher.py"
            script.write_text(FAKE_PS_SCRIPT)
            command[3:4] = [str(script), str(TOOL.parent), str(helper), str(self.fake_ps)]
        elif self.linux_ps and not real:  # Revision 12 (P10): Linux ps start times imitated
            command[3:4] = [str(self.linux_ps_script), str(TOOL.parent), str(helper),
                            str(LINUX_PS_LAG)]
        if queue is True:
            queue = self.queue
        if queue:
            command += ["--queue", str(queue)]
        return command + [str(a) for a in args]

    def tool(self, *args, home=None, queue=True, env=None, real=False):
        return subprocess.run(self.command(args, home, queue, real), cwd=self.tmp,
                              env=dict(self.env, **(env or {})), text=True,
                              capture_output=True, timeout=60)

    def popen_tool(self, *args, env=None, new_session=True):
        out = open(self.work / f"popen-{len(self.procs)}.out", "w")
        self.addCleanup(out.close)
        proc = subprocess.Popen(self.command(args, None, True), cwd=self.tmp,
                                env=dict(self.env, **(env or {})), stdout=out,
                                stderr=subprocess.STDOUT, start_new_session=new_session)
        self.procs.append(proc)
        return proc

    def popen_run(self, job, cmd, probe=None, env=None, new_session=True):
        self.counter += 1
        probe = probe or "echo " + shlex.quote(str(self.src))
        return self.popen_tool("run", job, "unit", "--log", self.work / f"bg-{self.counter}.log",
                               "--probe", probe, "--", *cmd, env=env, new_session=new_session)

    def wait_cmd(self, mode):
        self.counter += 1
        ready, release = self.work / f"ready-{self.counter}", self.work / f"release-{self.counter}"
        return [sys.executable, str(self.wait_script), mode, str(ready), str(release)], ready, release

    def init_args(self, qid, lock, floor=0, allowed=True, records=None, space=None,
                  jobs=None, extra=()):
        args = ["init", "--records", records or self.records, "--queue-id", qid, "--grant",
                self.grant, "--jobs", jobs or self.jobs, "--lock-dir", lock, "--min-free-gib",
                floor, "--space-path", space or self.tmp, *extra]
        return args + (["--allowed-root", self.src.parent] if allowed else [])

    def use_new_queue(self, qid, **options):
        """Init another queue in its own home and lock and make it the default for helpers."""
        home = self.tmp / f"home-{qid}"
        home.mkdir()
        self.home, self.lock = home, self.tmp / f"RUN-{qid}.lock"
        self.queue = self.init_queue(qid, home, self.lock, **options)

    def init_queue(self, qid, home, lock, **options):
        result = self.tool(*self.init_args(qid, lock, **options), home=home, queue=None)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        expected = (options.get("records") or self.records) / qid
        self.assertEqual(result.stdout.strip().splitlines()[-1], str(expected))
        return expected

    def ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def refused(self, result, *codes):
        self.assertIn(result.returncode, codes, result.stdout + result.stderr)
        self.assertTrue(any(line.startswith("REFUSED:") for line in result.stderr.splitlines()),
                        result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def status(self, **where):
        return json.loads(self.ok(self.tool("status", **where)).stdout)

    def state(self, job, **where):
        return job_entry(self.status(**where), job)["state"]

    def set(self, job, state, *extra, **where):
        result = self.tool("set", job, state, *extra, **where)
        if state == "Ready for approval" and "--packet" in extra and not result.returncode:
            self.ready_packets[job] = Path(extra[list(extra).index("--packet") + 1])
        return result

    def approve(self, job, quote, manifest=True):
        """set Approved with the quote and (Revision 4) the stored packet's manifest identity."""
        extra = ["--quote", quote] if quote is not None else []
        if manifest is True:
            manifest = sha256_file(self.ready_packets[job] / "MANIFEST.sha256")
        if manifest:
            extra += ["--manifest", manifest]
        return self.set(job, "Approved", *extra)

    def start_worker(self, job, **where):
        self.counter += 1
        worker = f"w{self.counter}"
        self.ok(self.tool("worker", "start", job, "--id", worker, **where))
        self.workers[str(where.get("queue", self.queue))] = (job, worker)

    def end_worker(self, **where):
        job, worker = self.workers.pop(str(where.get("queue", self.queue)))
        self.ok(self.tool("worker", "end", job, "--id", worker, "--outcome", "finished", **where))

    def prepare(self, job, enter=True, **where):
        """Lock, criterion, candidate; with `enter`, also Preparing and a worker for the job."""
        self.ok(self.tool("lock", "take", **where))
        self.ok(self.tool("criterion", job, "--file", self.criterion, **where))
        self.ok(self.tool("candidate", job, "--repo", self.repo, "--commit", "HEAD", **where))
        if enter:
            self.ok(self.set(job, "Preparing", **where))
            self.start_worker(job, **where)

    def ready(self, job, *extra):
        """Ready for approval with no worker registered (B4); after a refusal a
        worker is registered again so the test can go on running."""
        if str(self.queue) in self.workers:
            self.end_worker()
        result = self.set(job, "Ready for approval", *extra)
        if result.returncode:
            self.start_worker(job)
        return result

    def run_cmd(self, job, label, code=0, probe=True, cmd=None, extra=(), **where):
        self.counter += 1
        args = ["run", job, label, "--log", self.work / f"run-{self.counter}.log", *extra]
        if probe is True:
            probe = "echo " + shlex.quote(str(self.src / "build"))
        if probe:
            args += ["--probe", probe]
        cmd = cmd or [sys.executable, "-c", f"import sys; sys.exit({code})"]
        return self.tool(*args, "--", *cmd, **where)

    def held_probe(self, name, ignore_term=False):
        """Checkpoint 3: a probe that marks its start, then waits (at most 30 s) for a
        release file the test creates, then names the source; replaces fixed `sleep N`
        probe windows, so a test no longer depends on how fast helper subprocesses start.
        With ignore_term (checkpoint 5) the probe and its children ignore SIGTERM, so a
        stop (which signals the run's marked processes) cannot end it early.
        Revision 15: the waiting is done by a python child of the probe's shell (the helper
        runs every probe with /bin/sh -c), which ignores SIGTERM itself with ignore_term,
        instead of touch, a sleep loop and echo.
        Returns (probe, started, release)."""
        started, release = self.work / f"{name}-started", self.work / f"{name}-release"
        code = ("import os, signal, sys, time\n"
                "if sys.argv[4] == '1':\n"
                "    signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
                "open(sys.argv[1], 'w').close()\n"
                "deadline = time.time() + 30\n"
                "while not os.path.exists(sys.argv[2]) and time.time() < deadline:\n"
                "    time.sleep(0.1)\n"
                "print(sys.argv[3])\n")
        probe = (("trap '' TERM; " if ignore_term else "")
                 + " ".join(shlex.quote(str(word)) for word in (
                     sys.executable, "-c", code, started, release, self.src,
                     "1" if ignore_term else "0")))
        return probe, started, release

    def marker_cmd(self, name):
        marker = self.work / name
        return marker, [sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"]

    def attempt(self, job, code, **where):
        """A run that must be refused with `code` (CMD never starts) or accepted (0)."""
        self.counter += 1
        marker, cmd = self.marker_cmd(f"attempt-{self.counter}")
        result = self.run_cmd(job, "unit", cmd=cmd, **where)
        if code:
            self.refused(result, code)
        else:
            self.ok(result)
        self.assertEqual(marker.exists(), not code)

    def runs(self, queue=None):
        found = {}
        for path in sorted(Path(queue or self.queue, "runs").glob("*.json")):
            found[path.stem] = json.loads(path.read_text())
        return found

    def recorded_candidate_id(self, tree):
        """The candidate id recorded for a job whose candidate includes `tree`, else
        the first recorded candidate id, else None (from `status`)."""
        jobs = self.status()["jobs"]
        entries = list(jobs.values()) if isinstance(jobs, dict) else jobs
        for entry in entries:
            if any(repo.get("tree") == tree for repo in entry.get("candidate") or []):
                return entry["candidate_id"]
        return next((entry["candidate_id"] for entry in entries
                     if entry.get("candidate_id")), None)

    def packet(self, tree, manifest=True, log="unit: exit 0\n", screening=SCREENING,
               candidate_id=True, trees=None):
        """A packet frozen with screen_check.py freeze (Revision 4); unfrozen without `manifest`.
        Revision 9 (Q1): it carries SCREENING.txt with `screening` (none when None).
        Revision 11: CODE_IDENTITY.txt has `candidate_id:` (True: the recorded id for
        `tree`; None: no line) and a `tested_tree:` line for each of `trees` (default
        [tree])."""
        if candidate_id is True:
            candidate_id = self.recorded_candidate_id(tree)
        self.counter += 1
        path = self.work / f"packet-{self.counter}"
        path.mkdir()
        identity = "candidate: synthetic\n"
        if candidate_id is not None:
            identity += f"candidate_id: {candidate_id}\n"
        identity += "".join(f"tested_tree: {t}\n" for t in ([tree] if trees is None else trees))
        (path / "CODE_IDENTITY.txt").write_text(identity)
        (path / "RUN_LOG.txt").write_text(log)
        if screening is not None:
            (path / "SCREENING.txt").write_text(screening)
        if manifest:
            frozen = subprocess.run([sys.executable, "-I", "-B", str(SCREEN_CHECK), "freeze",
                                     str(path)], text=True, capture_output=True)
            self.assertEqual(frozen.returncode, 0, frozen.stdout + frozen.stderr)
        return path

    def wait_for(self, condition, what, timeout=TIMEOUT):
        deadline = time.time() + timeout
        while time.time() < deadline:
            if condition():
                return
            time.sleep(0.05)
        self.fail("timed out waiting for " + what)

    def wait_record(self, before, status="running"):
        found = {}

        def check():
            for stem, record in self.runs().items():
                pid = record.get("pid") if status == "running" else record.get("wrapper_pid")
                if (stem not in before and record.get("status") == status
                        and isinstance(pid, int) and alive(pid)):
                    found[stem] = record
                    return True
            return False
        self.wait_for(check, f"a {status} run record")
        return next(iter(found.items()))

    def stop_until_stopped(self):
        results = []

        def check():
            results.append(self.tool("stop"))
            self.assertIn(results[-1].returncode, (0, 10), results[-1].stderr)
            return results[-1].returncode == 0
        self.wait_for(check, "stop to print Stopped")
        self.assertTrue(results[-1].stdout.strip().startswith("Stopped"), results[-1].stdout)
        return results[-1]

    def stopping(self, result, *names):
        self.assertEqual(result.returncode, 10, result.stdout + result.stderr)
        self.assertTrue(result.stdout.strip().startswith("Stopping"), result.stdout)
        for name in names:
            self.assertIn(name, result.stdout)

    def report(self):
        result = self.ok(self.tool("report"))
        self.assertIn(REPORT_HEADER, result.stdout.splitlines())
        return result.stdout

    def first_line(self, text):
        return [line for line in text.splitlines() if line.strip()][0].strip()

    def result_cell(self, report, job, column=1):
        rows = []
        for line in report.splitlines():
            if line.startswith("|") and line.strip() != REPORT_HEADER:
                cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
                if job in cells[0]:
                    rows.append(cells)
        self.assertEqual(len(rows), 1, report)
        return rows[0][column]

    # ---- queue, ACTIVE, init -------------------------------------------------

    def test_active_resolution_with_and_without_queue(self):
        """Global options / init / M4: ACTIVE names the queue; without --queue it is
        used, with --queue it is not needed; missing, empty or naming a directory
        without QUEUE.json: exit 2, `no active queue`."""
        self.assertEqual((self.home / "ACTIVE").read_text().splitlines()[0], str(self.queue))
        for name in ("QUEUE.json", "GRANT.md", "JOBS.json"):
            self.assertTrue((self.queue / name).is_file(), name)
        via_active = json.loads(self.ok(self.tool("status", queue=None)).stdout)
        self.assertEqual(via_active, self.status())
        self.assertEqual(job_entry(via_active, ALPHA)["state"], "Queued")
        self.assertEqual(job_entry(via_active, BETA)["state"], "Queued")
        self.ok(self.set(ALPHA, "Preparing", queue=None))
        self.assertEqual(self.state(ALPHA), "Preparing")
        empty_home = self.tmp / "home-empty"
        empty_home.mkdir()
        for content in (None, "", str(self.outside) + "\n"):
            if content is not None:
                (empty_home / "ACTIVE").write_text(content)
            result = self.tool("status", home=empty_home, queue=None)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertIn("no active queue", result.stdout + result.stderr)
        self.assertEqual(self.state(ALPHA, home=empty_home), "Preparing")  # explicit --queue

    def test_init_refused_while_active_queue_unfinished(self):
        """init: exit 8 while ACTIVE names an unfinished queue; exit 2 if ROOT/ID
        exists or (M3) the --space-path does not exist."""
        result = self.tool(*self.init_args("q2", self.tmp / "L2"), queue=None)
        self.refused(result, 8)
        self.assertFalse((self.records / "q2").exists())
        self.assertEqual((self.home / "ACTIVE").read_text().splitlines()[0], str(self.queue))
        other_home = self.tmp / "home-other"
        other_home.mkdir()
        self.refused(self.tool(*self.init_args(QID, self.tmp / "L3"), home=other_home,
                               queue=None), 2)
        self.refused(self.tool(*self.init_args("q3", self.tmp / "L3", space=self.tmp / "nope"),
                               home=other_home, queue=None), 2)
        self.assertFalse((self.records / "q3").exists())
        self.init_queue("q3", other_home, self.tmp / "L3")  # control: existing space path
        for job in (ALPHA, BETA):
            self.ok(self.set(job, "Discarded"))
        second = self.init_queue("q2", self.home, self.tmp / "L2")
        self.assertEqual((self.home / "ACTIVE").read_text().splitlines()[0], str(second))

    # ---- set -----------------------------------------------------------------

    def test_illegal_transitions_refused_and_state_kept(self):
        """set: only the listed transitions; anything else exit 8, no state change,
        no EVENTS.log line."""
        events = self.queue / "EVENTS.log"

        def lines():
            return len(events.read_text().splitlines()) if events.exists() else 0
        steps = [(ALPHA, "Testing", 8), (ALPHA, "Ready for approval", 8), (ALPHA, "Approved", 8),
                 (ALPHA, "Waiting", 8), (ALPHA, "Published", 8), (ALPHA, "Preparing", 0),
                 (ALPHA, "Approved", 8), (ALPHA, "Queued", 8), (ALPHA, "Testing", 0),
                 (ALPHA, "Discarded", 8), (ALPHA, "Blocked", 0), (ALPHA, "Testing", 8),
                 (ALPHA, "Discarded", 0), (ALPHA, "Preparing", 8), (ALPHA, "Queued", 8),
                 (BETA, "Preparing", 0), (BETA, "Waiting", 0), (BETA, "Discarded", 8),
                 (BETA, "Testing", 0)]
        for job, target, code in steps:
            with self.subTest(job=job, target=target):
                before_state, before_lines = self.state(job), lines()
                result = self.set(job, target)
                if code:
                    self.refused(result, code)
                    self.assertEqual(self.state(job), before_state)
                    self.assertEqual(lines(), before_lines)
                else:
                    self.ok(result)
                    self.assertEqual(self.state(job), target)
                    self.assertGreater(lines(), before_lines)

    def test_set_cannot_write_run_candidate_criterion_fields(self):
        """set: no option writes candidate, criterion, run or worker fields; unknown
        options are a usage error (exit 2) and change nothing."""
        before = self.status()
        for option in ("--candidate", "--candidate-id", "--candidate_id", "--criterion",
                       "--criterion-sha256", "--criterion_sha256", "--run", "--runs",
                       "--stale", "--status", "--exit-code", "--tested-paths", "--worker",
                       "--commit", "--tree"):
            with self.subTest(option=option):
                result = self.set(ALPHA, "Preparing", option, "x")
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertEqual(self.status(), before)
        result = self.set(ALPHA, "Preparing", "candidate_id=x")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(self.status(), before)
        self.ok(self.set(ALPHA, "Preparing", "--note", "ordinary note"))
        self.assertEqual(self.state(ALPHA), "Preparing")

    def test_atomic_write_leftovers_are_never_read(self):
        """Files: a leftover temp file never changes what is read; writes leave none."""
        before_status = self.status()
        decoy = json.loads((self.queue / "QUEUE.json").read_text())
        job_entry(decoy, ALPHA)["state"] = "Published"
        (self.queue / "runs").mkdir(exist_ok=True)
        planted = {"QUEUE.json.tmp": "{broken", ".QUEUE.json.1234.tmp": json.dumps(decoy),
                   "runs/.r1.json.tmp": "{broken", "runs/r2.json.tmp": "{\"status\": \"running\""}
        before_entries = {p.relative_to(self.queue).as_posix() for p in self.queue.rglob("*")}
        for name, text in planted.items():
            (self.queue / name).write_text(text)
        self.assertEqual(self.status(), before_status)
        self.ok(self.set(ALPHA, "Preparing"))
        self.assertEqual(self.state(ALPHA), "Preparing")
        for command in ("reconcile", "report"):
            self.assertNotIn("Traceback", self.ok(self.tool(command)).stderr)
        after = {p.relative_to(self.queue).as_posix() for p in self.queue.rglob("*")}
        self.assertEqual(after - before_entries - set(planted) - {".queue.lock"}, set())

    # ---- run -----------------------------------------------------------------

    def test_run_record_and_exit_status(self):
        """run: record written (starting, wrapper pid) before CMD, CMD output to the
        log, then `completed` with exit_code; tool exits with CMD's code or 128+signal."""
        self.prepare(BETA)
        peek = self.work / "peek.py"
        peek.write_text(PEEK_SCRIPT)
        result = self.run_cmd(BETA, "unit", probe="pwd", extra=["--cwd", self.src],
                              cmd=[sys.executable, str(peek), str(self.queue / "runs"), "13"])
        self.assertEqual(result.returncode, 13, result.stdout + result.stderr)
        (stem, record), = self.runs().items()
        self.assertEqual((record["job"], record["label"]), (BETA, "unit"))
        self.assertEqual(record["status"], "completed")
        self.assertEqual(record["exit_code"], 13)
        self.assertIs(record["stale"], False)
        self.assertEqual(record["tested_paths"], [str(self.src)])
        self.assertEqual(record["criterion_sha256"], sha256_file(self.criterion))
        self.assertTrue(record["candidate_id"])
        self.assertIsInstance(record["pid"], int)
        self.assertIsInstance(record["wrapper_pid"], int)
        log = Path(record["log"]).read_text()
        self.assertRegex(log, r"seen (starting|running)")
        self.assertIn(f"cwd {self.src}", log)
        self.assertIn("stderr-line", log)
        killed = self.run_cmd(BETA, "unit", cmd=[sys.executable, "-c",
                                                  "import os, signal; os.kill(os.getpid(), 15)"])
        self.assertEqual(killed.returncode, 128 + signal.SIGTERM, killed.stderr)
        passed = self.run_cmd(BETA, "unit")
        self.assertEqual(passed.returncode, 0, passed.stderr)
        statuses = sorted((r["status"], r["exit_code"]) for s, r in self.runs().items() if s != stem)
        self.assertEqual(statuses[0][0], "completed")
        self.assertIn(("completed", 0), statuses)

    def test_run_refusal_order_lock_criterion_candidate(self):
        """run: lock not held 4 before missing criterion 6; no candidate 6; CMD never starts."""
        self.ok(self.set(BETA, "Preparing"))
        self.start_worker(BETA)
        self.attempt(BETA, 4)
        self.ok(self.tool("lock", "take"))
        self.attempt(BETA, 6)
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        self.attempt(BETA, 6)
        self.ok(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD"))
        self.attempt(BETA, 0)

    def test_run_requires_preparing_or_testing_and_own_worker(self):
        """Revision 2 M4: run exit 8 unless the job is Preparing/Testing, exit 9 unless
        the registered, unended worker is for this job."""
        self.prepare(BETA)
        self.attempt(BETA, 0)
        for state, back in (("Blocked", "Preparing"), ("Waiting", "Testing")):
            self.ok(self.set(BETA, state))
            self.attempt(BETA, 8)
            self.ok(self.set(BETA, back))
            self.attempt(BETA, 0)
        self.end_worker()
        self.attempt(BETA, 9)
        self.start_worker(ALPHA)
        self.attempt(BETA, 9)
        self.end_worker()
        self.start_worker(BETA)
        self.attempt(BETA, 0)
        self.end_worker()
        self.prepare(ALPHA, enter=False)
        self.start_worker(ALPHA)
        self.attempt(ALPHA, 8)  # Queued
        self.ok(self.set(ALPHA, "Preparing"))
        self.attempt(ALPHA, 0)

    def test_probe_path_outside_allowed_roots_refused(self):
        """run / S2: a probe path outside the allowed roots (after realpath, relative
        lines against --cwd) or missing exit 7."""
        self.prepare(BETA)
        link = self.src / "escape"
        link.symlink_to(self.outside)
        inside = shlex.quote(str(self.src / "build"))
        probes = [("echo " + shlex.quote(str(self.outside)), ()),
                  ("echo " + shlex.quote(str(link)), ()),
                  (f"echo {inside}; echo {shlex.quote(str(self.outside))}", ()),
                  ("echo " + shlex.quote(str(self.src / "missing")), ()),
                  ("echo allowed/src/build", ("--cwd", self.outside))]  # exists only from tmp
        for probe, extra in probes:
            with self.subTest(probe=probe):
                marker, cmd = self.marker_cmd("ran-outside")
                self.refused(self.run_cmd(BETA, "unit", probe=probe, cmd=cmd, extra=extra), 7)
                self.assertFalse(marker.exists())
        marker, cmd = self.marker_cmd("ran-inside")
        self.ok(self.run_cmd(BETA, "unit", probe=f"echo {inside}; echo", cmd=cmd))
        self.assertTrue(marker.exists())
        self.ok(self.run_cmd(BETA, "unit", probe="echo build", extra=("--cwd", self.src)))
        # Revision 2: refused probes may leave records (written before the probe);
        # none of them is completed.
        records = [r for r in self.runs().values() if r.get("status") == "completed"]
        self.assertEqual(len(records), 2)
        for record in records:
            self.assertEqual(record["tested_paths"], [str(self.src / "build")])

    def test_run_refused_below_floor(self):
        """run: free space at the space path below the floor exit 5 (after the lock
        check, before the criterion check); floor 0 accepted."""
        free_gib = shutil.disk_usage(self.tmp).free / 2 ** 30
        high_home, low_home = self.tmp / "home-high", self.tmp / "home-low"
        high_home.mkdir()
        low_home.mkdir()
        high = self.init_queue("q-high", high_home, self.tmp / "L-high",
                               floor=math.ceil(free_gib) + 1000)
        where = dict(home=high_home, queue=high)
        self.ok(self.set(BETA, "Preparing", **where))
        self.start_worker(BETA, **where)
        self.attempt(BETA, 4, **where)
        self.ok(self.tool("lock", "take", **where))
        self.attempt(BETA, 5, **where)
        self.prepare(BETA, enter=False, **where)
        self.attempt(BETA, 5, **where)
        # control: floor 0, and no allowed roots (no probe path can then be accepted)
        low = self.init_queue("q-low", low_home, self.tmp / "L-low", floor=0, allowed=False)
        where = dict(home=low_home, queue=low)
        self.prepare(BETA, **where)
        marker, cmd = self.marker_cmd("ran-low")
        self.ok(self.run_cmd(BETA, "unit", probe=None, cmd=cmd, **where))
        self.assertTrue(marker.exists())
        self.refused(self.run_cmd(BETA, "unit", **where), 7)

    # ---- readiness, approval, stale ---------------------------------------------

    def test_ready_for_approval_conditions_and_approval_quote(self):
        """set Ready for approval: criterion+candidate, a passing current run per
        label with non-empty tested_paths, a packet; Approved needs a quote;
        M4: criterion/candidate refused (8) in Published and Discarded."""
        self.ok(self.set(BETA, "Preparing"))
        self.ok(self.set(BETA, "Testing"))
        self.refused(self.set(BETA, "Ready for approval", "--packet", self.packet(self.tree())), 8)
        self.prepare(ALPHA)
        self.ok(self.set(ALPHA, "Testing"))
        good = self.packet(self.tree())
        self.refused(self.ready(ALPHA, "--packet", good), 8)  # no runs
        self.ok(self.run_cmd(ALPHA, "unit"))
        self.refused(self.ready(ALPHA, "--packet", good), 8)  # no build
        self.ok(self.run_cmd(ALPHA, "build", probe=None))
        self.refused(self.ready(ALPHA, "--packet", good), 8)  # empty paths
        self.assertEqual(self.run_cmd(ALPHA, "build", code=1).returncode, 1)
        self.refused(self.ready(ALPHA, "--packet", good), 8)  # failed
        self.ok(self.run_cmd(ALPHA, "build"))
        self.refused(self.ready(ALPHA), 8)
        self.refused(self.ready(ALPHA, "--packet", self.packet(self.tree(), manifest=False)), 8)
        self.refused(self.ready(ALPHA, "--packet", self.packet("0" * 40)), 8)
        self.assertEqual(self.state(ALPHA), "Testing")
        self.ok(self.ready(ALPHA, "--packet", good))
        self.assertEqual(self.state(ALPHA), "Ready for approval")
        self.assertIn("local checks passed", self.result_cell(self.report(), ALPHA))
        self.refused(self.approve(ALPHA, None), 8)
        self.refused(self.approve(ALPHA, "approved, go ahead"), 8)
        quote = f"I approve {ALPHA} as tested"
        self.ok(self.approve(ALPHA, quote))
        self.assertIn(quote, json.dumps(self.status()))
        self.ok(self.tool("criterion", ALPHA, "--file", self.criterion))  # control: same text
        self.ok(self.set(ALPHA, "Published"))
        for target in ("Preparing", "Discarded", "Approved"):
            self.refused(self.set(ALPHA, target), 8)
        self.assertEqual(self.state(ALPHA), "Published")
        self.ok(self.set(BETA, "Blocked"))
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))  # control: Blocked
        self.ok(self.set(BETA, "Discarded"))
        for job in (ALPHA, BETA):
            with self.subTest(final=job):
                self.refused(self.tool("criterion", job, "--file", self.criterion), 8)
                self.refused(self.tool("candidate", job, "--repo", self.repo, "--commit",
                                       "HEAD"), 8)

    def test_approval_quote_whole_token_and_readiness_recheck(self):
        """Revision 2 B5: the quote must contain the job id as a whole token; Approved
        re-checks readiness and is refused (8) once a run became stale."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        for quote in (f"approve {BETA}2", f"approve x{BETA}", f"approve {BETA}-old",
                      f"approve {BETA}_x", f"approve sub-{BETA}", "approve job-bet"):
            with self.subTest(quote=quote):
                self.refused(self.approve(BETA, quote), 8)
                self.assertEqual(self.state(BETA), "Ready for approval")
        self.criterion.write_text("Criterion two: changed after readiness.\n")
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        self.refused(self.approve(BETA, f"Approve {BETA}."), 8)
        self.assertEqual(self.state(BETA), "Ready for approval")
        self.ok(self.set(BETA, "Preparing"))
        self.start_worker(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        self.ok(self.approve(BETA, f"Approve {BETA}."))

    def test_ready_refused_while_worker_or_run_live(self):
        """Revision 2 B4: Ready for approval refused (8) while a worker is registered
        and not ended, and while a run of the job is running."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        packet = self.packet(self.tree())
        self.refused(self.set(BETA, "Ready for approval", "--packet", packet), 8)
        cmd, ready, release = self.wait_cmd("default")
        proc = self.popen_run(BETA, cmd)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the held run")
        self.end_worker()
        self.refused(self.set(BETA, "Ready for approval", "--packet", packet), 8)
        release.write_text("")
        self.assertEqual(proc.wait(TIMEOUT), 0)
        self.assertEqual(self.runs()[stem]["status"], "completed")
        self.ok(self.set(BETA, "Ready for approval", "--packet", packet))

    def test_stale_after_criterion_change(self):
        """criterion: same text no change; different text marks runs stale, logs old
        and new hashes; report never says passed; readiness refused until rerun."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        old_hash = sha256_file(self.criterion)
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        self.assertFalse(any(r["stale"] for r in self.runs().values()))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.criterion.write_text("Criterion two: the tool exits 0 and logs.\n")
        new_hash = sha256_file(self.criterion)
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        self.assertTrue(all(r["stale"] is True for r in self.runs().values()))
        events = (self.queue / "EVENTS.log").read_text()
        self.assertIn(old_hash, events)
        self.assertIn(new_hash, events)
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("stale", cell)
        self.assertNotIn("passed", cell.lower())
        self.ok(self.set(BETA, "Preparing"))
        self.ok(self.set(BETA, "Testing"))
        self.refused(self.ready(BETA, "--packet", self.packet(self.tree())), 8)
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))

    def test_stale_after_candidate_change(self):
        """candidate: same candidate no change; unresolvable commit exit 2; a new
        candidate marks runs stale and readiness needs a new run and its tree."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        (first,) = self.runs().values()
        self.ok(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD"))
        bad = self.tool("candidate", BETA, "--repo", self.repo, "--commit", "no-such-rev")
        self.assertEqual(bad.returncode, 2, bad.stdout + bad.stderr)
        self.assertFalse(any(r["stale"] for r in self.runs().values()))
        old_tree = self.tree()
        self.commit("second")
        new_tree = self.tree()
        self.ok(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD"))
        self.assertTrue(all(r["stale"] is True for r in self.runs().values()))
        self.refused(self.ready(BETA, "--packet", self.packet(new_tree)), 8)
        self.ok(self.run_cmd(BETA, "unit"))
        fresh = [r for r in self.runs().values() if not r["stale"]]
        self.assertEqual(len(fresh), 1)
        self.assertNotEqual(fresh[0]["candidate_id"], first["candidate_id"])
        self.refused(self.ready(BETA, "--packet", self.packet(old_tree)), 8)
        self.ok(self.ready(BETA, "--packet", self.packet(new_tree)))

    # ---- interrupted and group-left runs ------------------------------------------

    def test_wrapper_sigterm_records_interrupted(self):
        """Revision 2 B3: a wrapper that receives SIGTERM ends the run `interrupted`
        with exit_code, even when CMD traps TERM and exits 0; never readiness."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        cmd, ready, _ = self.wait_cmd("exit0")
        proc = self.popen_run(BETA, cmd)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the trap to be installed")
        os.kill(proc.pid, signal.SIGTERM)
        proc.wait(TIMEOUT)
        record = self.runs()[stem]
        self.assertEqual(record["status"], "interrupted")
        self.assertIn("exit_code", record)
        packet = self.packet(self.tree())
        self.refused(self.ready(BETA, "--packet", packet), 8)
        self.assertNotIn("passed", self.result_cell(self.report(), BETA).lower())
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", packet))

    def launch(self, *modes):
        """A run whose CMD leaves children (one per mode) in its group and exits 0."""
        launcher, pids = self.work / "launch.py", self.work / f"pids-{self.counter}"
        launcher.write_text(LAUNCH_SCRIPT)
        specs, releases = [], []
        for mode in modes:
            cmd, ready, release = self.wait_cmd(mode)
            specs.append(f"{mode},{ready},{release}")
            releases.append(release)
        self.ok(self.run_cmd(BETA, "unit", probe="echo " + shlex.quote(str(self.src)),
                             cmd=[sys.executable, str(launcher), "exit", str(self.wait_script),
                                  str(pids), *specs]))
        children = [int(pid) for pid in pids.read_text().split()]
        self.stray.extend(children)
        self.assertTrue(all(alive(pid) for pid in children))
        (stem, record), = self.runs().items()
        self.assertIs(record.get("group_left_running"), True, record)
        return stem, children, releases

    def test_group_left_running_never_satisfies_readiness(self):
        """Revision 2 S1/B4: a leader exiting with live group members records
        group_left_running; readiness refused while the group lives and after."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        _, (child,), (release,) = self.launch("default")
        packet = self.packet(self.tree())
        self.refused(self.ready(BETA, "--packet", packet), 8)
        release.write_text("")
        self.wait_for(lambda: not alive(child), "the grandchild to exit")
        self.refused(self.ready(BETA, "--packet", packet), 8)
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", packet))

    def test_stop_signals_group_left_running_never_stopped_while_alive(self):
        """Revision 2 S1: stop sends SIGTERM to a group_left_running group and is
        never Stopped while a member lives."""
        self.prepare(BETA)
        stem, (deaf, plain), (release, _) = self.launch("ignore", "default")
        self.end_worker()
        self.stopping(self.tool("stop"))
        self.wait_for(lambda: not alive(plain), "the signalled grandchild to exit", 5)
        self.assertTrue(alive(deaf))
        report = self.report()
        self.assertEqual(self.first_line(report), "**Stopping**")
        self.assertNotIn("**Stopped**", report)
        release.write_text("")
        self.wait_for(lambda: not alive(deaf), "the TERM-ignoring grandchild to exit")
        self.stop_until_stopped()
        self.assertEqual(self.first_line(self.report()), "**Stopped**")

    # ---- stop, check-stop, reconcile, worker ----------------------------------

    def test_check_stop(self):
        """check-stop: exit 0 without a stop; exit 3 and STOP REQUESTED after one."""
        result = self.ok(self.tool("check-stop"))
        self.assertNotIn("STOP REQUESTED", result.stdout)
        self.assertFalse((self.queue / "STOP_REQUESTED").exists())
        self.assertTrue(self.ok(self.tool("stop")).stdout.strip().startswith("Stopped"))
        result = self.tool("check-stop")
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertIn("STOP REQUESTED", result.stdout)
        self.assertIn(str(self.tmp), (self.queue / "STOP_REQUESTED").read_text())  # cwd
        self.assertEqual(self.first_line(self.report()), "**Stopped**")

    def test_run_and_transitions_refused_under_stop(self):
        """stop: afterwards run (even before the lock check), worker start,
        criterion, candidate and transitions other than to Blocked/Waiting are
        refused with exit 3."""
        self.prepare(BETA)
        self.attempt(BETA, 0)
        self.ok(self.set(ALPHA, "Preparing"))
        self.ok(self.set(ALPHA, "Testing"))  # control: Preparing->Testing is legal
        self.stopping(self.tool("stop"), BETA)  # the worker is still registered
        self.attempt(BETA, 3)
        self.ok(self.tool("lock", "release"))
        self.attempt(BETA, 3)
        self.refused(self.set(BETA, "Testing"), 3)
        self.assertEqual(self.state(BETA), "Preparing")
        self.ok(self.set(BETA, "Blocked", "--note", "stopped by the Developer"))
        self.refused(self.set(BETA, "Preparing"), 3)
        self.ok(self.set(ALPHA, "Waiting"))
        self.refused(self.tool("criterion", BETA, "--file", self.criterion), 3)
        self.refused(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD"), 3)
        self.end_worker()
        self.refused(self.tool("worker", "start", BETA, "--id", "w-late"), 3)

    def test_stop_from_separate_process_stopping_then_stopped(self):
        """stop / Revision 2 B3: while a tracked run is alive the result is Stopping
        (exit 10) naming it; a run that outlives the stop ends `interrupted`, not
        completed; once it is gone, Stopped."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("ignore")
        proc = self.popen_run(BETA, cmd)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the held run to ignore SIGTERM")
        self.end_worker()
        self.stopping(self.tool("stop"), stem)
        self.assertEqual(self.tool("check-stop").returncode, 3)
        report = self.report()
        self.assertEqual(self.first_line(report), "**Stopping**")
        self.assertNotIn("**Stopped**", report)
        self.assertIn("running", self.result_cell(report, BETA))
        release.write_text("")
        proc.wait(TIMEOUT)
        record = self.runs()[stem]
        self.assertEqual(record["status"], "interrupted")
        self.assertIn("exit_code", record)
        self.stop_until_stopped()
        report = self.report()
        self.assertEqual(self.first_line(report), "**Stopped**")
        self.assertNotIn("passed", self.result_cell(report, BETA).lower())

    def test_stop_during_probe_never_starts_cmd(self):
        """Revision 2 B2: a `starting` record exists during the probe; a stop then
        means CMD never starts (`not started`, exit 3) and no Stopped while the
        wrapper lives. Checkpoint 3: the probe waits for a release file (bounded,
        30 s) that the test creates only after the stop has been recorded
        (STOP_REQUESTED exists), instead of a fixed `sleep 2`; so the stop always
        falls inside the probe, however slowly helper subprocesses start.
        Checkpoint 5 (review CP3 finding 2): the probe stays held until the stop
        command itself has finished, so its result is read while the starting
        wrapper is certainly alive: it must be Stopping (exit 10), never Stopped.
        The probe ignores SIGTERM (stop signals the run's marked processes, which
        would otherwise end the probe and the wrapper at once). Only then is the
        probe released; the wrapper exits 3 and a later stop prints Stopped."""
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("cmd-ran")
        probe, started, release = self.held_probe("stop-probe", ignore_term=True)
        proc = self.popen_run(BETA, cmd, probe=probe)
        self.wait_for(started.exists, "the probe to start")
        stem, record = self.wait_record(set(), status="starting")
        self.assertEqual(record["wrapper_pid"], proc.pid)
        self.end_worker()
        out = self.work / "stop-during-probe.out"
        with open(out, "w") as handle:
            stopper = subprocess.Popen(self.command(["stop"], None, True), cwd=self.tmp,
                                       env=self.env, stdout=handle, stderr=subprocess.STDOUT)
            self.procs.append(stopper)
            self.wait_for((self.queue / "STOP_REQUESTED").exists, "the stop to be recorded")
            self.assertIsNone(proc.poll(), "the wrapper ended before the probe was released")
            self.assertEqual(self.runs()[stem]["status"], "starting")
            code = stopper.wait(TIMEOUT)  # the probe is still held: the wrapper lives
        self.assertIsNone(proc.poll(), "the wrapper ended while the probe was held")
        self.assertFalse(release.exists())
        self.assertNotEqual(code, 0, "Stopped while the starting wrapper lived: "
                            + out.read_text())
        self.assertEqual(code, 10, out.read_text())
        self.assertTrue(out.read_text().startswith("Stopping"), out.read_text())
        release.write_text("")
        self.assertEqual(proc.wait(TIMEOUT), 3)
        self.assertFalse(marker.exists())
        self.assertEqual(self.runs()[stem]["status"], "not started")
        self.stop_until_stopped()

    def test_stop_from_different_tz_locale_and_path(self):
        """Revision 2 B1: stop and reconcile run with another TZ, LC_ALL and PATH than
        the run still recognise the live run; stop signals it."""
        self.prepare(BETA)
        proc = self.popen_run(BETA, sleep_cmd(5), env={"TZ": "UTC", "LC_ALL": "C"})
        stem, record = self.wait_record(set())
        foreign = {"TZ": "America/Denver", "LC_ALL": "fr_FR.ISO8859-1", "PATH": "/nonexistent"}
        self.ok(self.tool("reconcile", env=foreign))
        self.assertEqual(self.runs()[stem]["status"], "running")
        self.end_worker()
        first = self.tool("stop", env=foreign)
        self.assertIn(first.returncode, (0, 10), first.stdout + first.stderr)
        self.assertEqual(proc.wait(4.5), 128 + signal.SIGTERM)
        self.assertFalse(alive(record["pid"]))
        self.stop_until_stopped()
        self.assertEqual(self.runs()[stem]["status"], "interrupted")

    def test_stop_terminates_tracked_sleep_not_a_reused_pid(self):
        """stop: SIGTERM to the tracked sleep's group (run exits 128+15); a record
        whose pid has a different start time is not signalled and reconcile marks
        it unknown."""
        self.prepare(BETA)
        proc = self.popen_run(BETA, sleep_cmd(5))
        stem, record = self.wait_record(set())
        self.end_worker()
        # The decoy must start in a later second than the tracked run (checkpoint 3:
        # bounded polling on ps start times instead of a fixed 2.1 s sleep).
        self.wait_past_second(self.start_epoch(record["pid"]))
        decoy = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(5)"],
                                 start_new_session=True)
        self.procs.append(decoy)
        fake = {k: (stem + "-reused" if v == stem else v) for k, v in record.items()}
        fake.update(pid=decoy.pid, pgid=decoy.pid)
        staging = self.work / "fake.json"
        staging.write_text(json.dumps(fake))
        os.replace(staging, self.queue / "runs" / (stem + "-reused.json"))
        first = self.tool("stop")
        self.assertIn(first.returncode, (0, 10), first.stdout + first.stderr)
        self.assertEqual(proc.wait(4), 128 + signal.SIGTERM)
        self.stop_until_stopped()
        self.assertIsNone(decoy.poll(), "the reused pid was signalled")
        self.ok(self.tool("reconcile"))
        runs = self.runs()
        self.assertEqual(runs[stem]["status"], "interrupted")
        self.assertEqual(runs[stem + "-reused"]["status"], "unknown")
        self.assertIsNone(decoy.poll())

    def test_killed_run_becomes_unknown_after_reconcile(self):
        """reconcile: a `running` record whose process is gone becomes unknown and
        is never relaunched; report shows unknown; it cannot satisfy readiness."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        proc = self.popen_run(BETA, sleep_cmd(5))
        stem, record = self.wait_record(set())
        self.ok(self.tool("reconcile"))  # control: alive run stays running
        self.assertEqual(self.runs()[stem]["status"], "running")
        self.assertIn("running", self.result_cell(self.report(), BETA))
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait(5)
        try:
            os.killpg(record["pgid"], signal.SIGKILL)
        except ProcessLookupError:
            pass
        self.wait_for(lambda: not alive(record["pid"]), "the killed run to disappear")
        self.assertEqual(self.runs()[stem]["status"], "running")
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("unknown", cell)
        self.assertNotIn("passed", cell.lower())
        result = self.ok(self.tool("reconcile"))
        self.assertIn(stem, result.stdout)
        runs = self.runs()
        self.assertEqual(list(runs), [stem])
        self.assertEqual(runs[stem]["status"], "unknown")
        self.ok(self.tool("reconcile"))
        self.assertEqual(list(self.runs()), [stem])
        packet = self.packet(self.tree())
        self.refused(self.ready(BETA, "--packet", packet), 8)
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", packet))

    def test_other_jobs_worker_does_not_block_readiness(self):
        """Revision 3 R1: only a worker registered for this job blocks readiness; a
        parked ready job stays `local checks passed` and can be Approved while
        another job's worker is registered."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        packet = self.packet(self.tree())
        self.refused(self.set(BETA, "Ready for approval", "--packet", packet), 8)  # own worker
        self.end_worker()
        self.ok(self.set(ALPHA, "Preparing"))
        self.start_worker(ALPHA)
        self.ok(self.set(BETA, "Ready for approval", "--packet", packet))
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("local checks passed", cell)
        self.assertNotIn("stale", cell.lower())
        self.ok(self.approve(BETA, f"approve {BETA}"))
        self.assertEqual(self.state(BETA), "Approved")
        self.assertIn(str(self.queue), self.workers)  # ALPHA's worker was registered throughout

    def test_stop_never_signals_unconfirmed_identity(self):
        """Revision 3 N4: a live `running` record whose start time cannot be confirmed
        is reported as remaining (Stopping, EVENTS.log line) but never signalled."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        self.end_worker()
        (stem, record), = self.runs().items()
        cmd, ready, _ = self.wait_cmd("default")  # lives longer than stop's wait
        decoy = subprocess.Popen(cmd, start_new_session=True)
        self.procs.append(decoy)
        self.wait_for(ready.exists, "the decoy to start")
        fake = {k: (stem + "-unconfirmed" if v == stem else v) for k, v in record.items()
                if k not in ("exit_code", "ended", "group_left_running")}
        fake.update(status="running", pid=decoy.pid, pgid=decoy.pid, wrapper_pid=decoy.pid,
                    process_start="unreadable", wrapper_start="unreadable")
        staging = self.work / "fake.json"
        staging.write_text(json.dumps(fake))
        os.replace(staging, self.queue / "runs" / (stem + "-unconfirmed.json"))
        self.stopping(self.tool("stop"), stem + "-unconfirmed")
        self.assertIsNone(decoy.poll(), "a process of unconfirmed identity was signalled")
        self.assertIn(stem + "-unconfirmed", (self.queue / "EVENTS.log").read_text())
        self.assertEqual(self.first_line(self.report()), "**Stopping**")
        decoy.kill()
        decoy.wait(5)
        self.stop_until_stopped()  # control: once the process is gone, Stopped

    def test_stop_with_unended_worker(self):
        """worker/stop: one active worker (exit 9 for a second or a wrong id); stop
        with a registered unended worker is Stopping, never Stopped, and worker
        start is refused; after worker end, Stopped."""
        self.ok(self.tool("worker", "start", ALPHA, "--id", "w0"))
        self.ok(self.tool("worker", "end", ALPHA, "--id", "w0", "--outcome", "finished"))
        self.ok(self.tool("worker", "start", ALPHA, "--id", "w1"))  # ended worker is no bar
        self.refused(self.tool("worker", "start", BETA, "--id", "w2"), 9)
        self.refused(self.tool("worker", "end", ALPHA, "--id", "w2", "--outcome", "finished"), 9)
        self.stopping(self.tool("stop"), ALPHA)
        report = self.report()
        self.assertEqual(self.first_line(report), "**Stopping**")
        self.assertNotIn("**Stopped**", report)
        self.refused(self.tool("worker", "start", BETA, "--id", "w2"), 3, 9)
        self.assertEqual(self.tool("stop").returncode, 10)
        self.ok(self.tool("worker", "end", ALPHA, "--id", "w1", "--outcome", "cancelled"))
        self.stop_until_stopped()
        self.assertEqual(self.first_line(self.report()), "**Stopped**")
        self.refused(self.tool("worker", "start", BETA, "--id", "w2"), 3)

    # ---- lock, report ---------------------------------------------------------

    def test_foreign_lock_never_removed(self):
        """lock / S3: a foreign or unreadable lock is never taken or removed (exit 4)
        and run is refused; our own lock is taken, re-taken, shown and released;
        release refused (4, nothing removed) with extra files in the lock dir."""
        self.lock.mkdir()
        owner = "queue_id: other-queue\npid: 1\nhost: elsewhere\ntime: 2026-10-01T00:00:00Z\n"
        (self.lock / "owner.txt").write_text(owner)
        for action in ("take", "release"):
            self.refused(self.tool("lock", action), 4)
            self.assertEqual((self.lock / "owner.txt").read_text(), owner)
        self.assertIn("other-queue", self.ok(self.tool("lock", "show")).stdout)
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        self.ok(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD"))
        self.ok(self.set(BETA, "Preparing"))
        self.start_worker(BETA)
        self.attempt(BETA, 4)
        (self.lock / "owner.txt").chmod(0)  # unreadable owner file
        for action in ("take", "release"):
            self.refused(self.tool("lock", action), 4)
            self.assertTrue((self.lock / "owner.txt").exists())
        (self.lock / "owner.txt").unlink()  # unreadable: no owner.txt
        for action in ("take", "release"):
            self.refused(self.tool("lock", action), 4)
            self.assertTrue(self.lock.is_dir())
            self.assertFalse((self.lock / "owner.txt").exists())
        self.lock.rmdir()
        self.assertIn("free", self.ok(self.tool("lock", "show")).stdout)
        self.ok(self.tool("lock", "take"))
        lines = (self.lock / "owner.txt").read_text().splitlines()
        for key in ("queue_id:", "queue_dir:", "pid:", "host:", "time:"):
            self.assertTrue(any(line.startswith(key) for line in lines), (key, lines))
        self.assertIn(f"queue_id: {QID}", lines)
        self.assertIn(f"queue_dir: {self.queue}", lines)
        self.assertIn("already held", self.ok(self.tool("lock", "take")).stdout)
        self.assertIn(f"queue_id: {QID}", self.ok(self.tool("lock", "show")).stdout)
        self.attempt(BETA, 0)
        extra = self.lock / ".owner.tmp"
        extra.write_text("leftover\n")
        self.refused(self.tool("lock", "release"), 4)
        self.assertTrue(extra.exists())
        self.assertTrue((self.lock / "owner.txt").exists())
        extra.unlink()
        self.ok(self.tool("lock", "release"))
        self.assertFalse(self.lock.exists())
        self.assertIn("free", self.ok(self.tool("lock", "show")).stdout)

    def test_same_id_queue_elsewhere_is_not_owner(self):
        """Revision 2 S3: ownership needs queue_id and queue_dir; a queue with the same
        id in another directory can neither re-take nor release our lock."""
        other_records, other_home = self.tmp / "records-other", self.tmp / "home-other"
        other_records.mkdir()
        other_home.mkdir()
        other = self.init_queue(QID, other_home, self.lock, records=other_records)
        there = dict(home=other_home, queue=other)
        self.ok(self.tool("lock", "take"))
        owner = (self.lock / "owner.txt").read_text()
        taken = self.tool("lock", "take", **there)
        self.refused(taken, 4)
        self.assertNotIn("already held", taken.stdout)
        self.refused(self.tool("lock", "release", **there), 4)
        self.assertEqual((self.lock / "owner.txt").read_text(), owner)
        self.ok(self.tool("lock", "release"))
        self.ok(self.tool("lock", "take", **there))  # control: free lock is taken
        self.assertIn(f"queue_dir: {other}", (self.lock / "owner.txt").read_text())
        self.ok(self.tool("lock", "release", **there))

    def test_report_blocked_and_waiting_rows(self):
        """report/wait: Blocked shows `blocked` and its note; wait obeys the
        transition rules and shows `waiting`; no stop line without a stop.
        Revision 9 (Q2): a Waiting row's step is exactly the manual-resume text and
        never says the queue resumes or continues automatically."""
        self.refused(self.tool("wait", ALPHA, "--until", "2026-10-09T06:00:00Z"), 8)
        self.ok(self.set(ALPHA, "Preparing"))
        self.ok(self.tool("wait", ALPHA, "--until", "2026-10-09T06:00:00Z"))
        self.assertEqual(self.state(ALPHA), "Waiting")
        self.ok(self.set(BETA, "Preparing"))
        self.ok(self.set(BETA, "Blocked", "--note", "needs disk space"))
        report = self.report()
        self.assertNotIn("**Stopp", report)
        self.assertIn("waiting", self.result_cell(report, ALPHA).lower())
        row = [line for line in report.splitlines()
               if line.startswith("|") and ALPHA in line.split("|")[1]]
        self.assertEqual(len(row), 1, report)
        # Revision 9 (Q2): the step is exactly the manual-resume text.
        self.assertEqual(self.result_cell(report, ALPHA, column=2), WAITING_STEP)
        self.assertNotIn("resumes automatically", report)
        self.assertNotIn("continues if the app resumes", report)
        cell = self.result_cell(report, BETA)
        self.assertIn("blocked", cell.lower())
        self.assertIn("needs disk space", cell)
        self.refused(self.tool("wait", BETA, "--unknown", "--max-retries", "3"), 8)
        self.ok(self.set(BETA, "Preparing"))
        self.ok(self.tool("wait", BETA, "--unknown", "--max-retries", "3"))
        self.assertEqual(self.state(BETA), "Waiting")
        report = self.report()
        self.assertIn("waiting", self.result_cell(report, BETA).lower())
        self.assertEqual(self.result_cell(report, BETA, column=2), WAITING_STEP)

    # ---- Revision 4: surviving processes ---------------------------------------

    def survivor_after_wrapper_kill(self, flags):
        """CMD starts one child (flags) and holds; the wrapper and then the CMD leader
        are SIGKILLed, leaving the child. Returns (stem, child pid, ready, release)."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        launcher, pids = self.work / "launch.py", self.work / "pids-survivor"
        launcher.write_text(LAUNCH_SCRIPT)
        _, ready, release = self.wait_cmd(flags)
        proc = self.popen_run(BETA, [sys.executable, str(launcher), "hold", str(self.wait_script),
                                     str(pids), f"{flags},{ready},{release}"])
        stem, record = self.wait_record(set())
        self.wait_for(pids.exists, "the child pid")
        (child,) = [int(pid) for pid in pids.read_text().split()]
        self.stray.append(child)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(5)
        os.kill(record["pid"], signal.SIGKILL)
        self.wait_for(lambda: not alive(record["pid"]), "the CMD leader to exit")
        self.assertTrue(alive(child))
        self.end_worker()
        return stem, child, ready, release

    def check_survivor_counted(self, stem, child, ready, release):
        """A live survivor: never unknown, never passed, readiness refused, stop
        Stopping and signals it; once it is gone: unknown and Stopped (control)."""
        self.ok(self.tool("reconcile"))
        self.assertNotEqual(self.runs()[stem]["status"], "unknown")
        self.assertNotIn("passed", self.result_cell(self.report(), BETA).lower())
        self.refused(self.set(BETA, "Ready for approval", "--packet", self.packet(self.tree())), 8)
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(ready) + ".term").exists, "the survivor to get SIGTERM", 5)
        self.assertTrue(alive(child))
        report = self.report()
        self.assertEqual(self.first_line(report), "**Stopping**")
        self.assertNotIn("passed", self.result_cell(report, BETA).lower())
        release.write_text("")
        self.wait_for(lambda: not alive(child), "the survivor to exit")
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[stem]["status"], "unknown")
        self.stop_until_stopped()

    def test_survivor_group_member_after_wrapper_kill(self):
        """Revision 4 fix 1 (a): a CMD group member outliving its leader is found."""
        self.check_survivor_counted(*self.survivor_after_wrapper_kill("note"))

    def test_survivor_new_process_group_same_session(self):
        """Revision 4 fix 1 (b): a child in a new process group, same session, is found."""
        self.check_survivor_counted(*self.survivor_after_wrapper_kill("note+pgrp"))

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "finds the survivor by its GC_AUTO_RUN_ID "
                     "environment marker, read with `ps -E` (Linux ps has no -E)")
    def test_survivor_setsid_found_by_run_id_marker(self):
        """Revision 4 fix 1 (c): a non-system program calling os.setsid() is found by
        its GC_AUTO_RUN_ID environment marker."""
        self.check_survivor_counted(*self.survivor_after_wrapper_kill("note+setsid"))

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "finds the surviving probe by its GC_AUTO_RUN_ID "
                     "environment marker, read with `ps -E` (Linux ps has no -E)")
    def test_probe_surviving_killed_wrapper_is_found(self):
        """Revision 4 fix 1: a probe that outlives a SIGKILLed wrapper is counted;
        CMD never starts; once it is gone, unknown and Stopped."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("note")
        probe = " ".join(shlex.quote(part) for part in cmd) + "; echo " + shlex.quote(str(self.src))
        marker, run_cmd = self.marker_cmd("cmd-after-probe")
        proc = self.popen_run(BETA, run_cmd, probe=probe)
        self.wait_for(ready.exists, "the probe to start")
        stem, _ = self.wait_record(set(), status="starting")
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(5)
        self.end_worker()
        self.ok(self.tool("reconcile"))
        self.assertNotEqual(self.runs()[stem]["status"], "unknown")
        self.assertNotIn("passed", self.result_cell(self.report(), BETA).lower())
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(ready) + ".term").exists, "the probe to get SIGTERM", 5)
        release.write_text("")
        self.stop_until_stopped()
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[stem]["status"], "unknown")
        self.assertFalse(marker.exists())

    def test_unrelated_process_in_same_session_not_counted(self):
        """Revision 4 control: a process in the test's own session (the wrapper's
        too) that the run did not start is neither counted nor signalled."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("default")
        proc = self.popen_run(BETA, cmd, new_session=False)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the run")
        outsider_cmd, outsider_ready, _ = self.wait_cmd("note")
        outsider = subprocess.Popen(outsider_cmd)  # same session, started by the test
        self.procs.append(outsider)
        self.wait_for(outsider_ready.exists, "the outsider")
        release.write_text("")
        self.assertEqual(proc.wait(TIMEOUT), 0)
        self.end_worker()
        self.assertEqual(self.runs()[stem]["status"], "completed")
        result = self.ok(self.tool("stop"))
        self.assertTrue(result.stdout.strip().startswith("Stopped"), result.stdout)
        self.assertFalse(Path(str(outsider_ready) + ".term").exists())
        self.assertIsNone(outsider.poll())

    # ---- Revision 4: frozen packet, installation, probe failure, retries ---------

    def test_packet_must_verify_for_ready_and_approved(self):
        """Revision 4 fix 2: Ready needs a packet that verifies (tampered, extra or
        missing file: 8); a change after Ready drops `passed` and refuses Approved;
        Approved needs the right --manifest."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        tree = self.tree()
        for damage in ("tampered", "extra", "missing"):
            with self.subTest(damage=damage):
                packet = self.packet(tree)
                if damage == "tampered":
                    (packet / "RUN_LOG.txt").write_text("unit: exit 1\n")
                elif damage == "extra":
                    (packet / "EXTRA.txt").write_text("not listed\n")
                else:
                    (packet / "RUN_LOG.txt").unlink()
                self.refused(self.ready(BETA, "--packet", packet), 8)
        good = self.packet(tree)
        self.ok(self.ready(BETA, "--packet", good))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        log = good / "RUN_LOG.txt"
        original = log.read_bytes()
        log.write_text("changed after Ready\n")
        self.assertNotIn("passed", self.result_cell(self.report(), BETA).lower())
        quote = f"approve {BETA}"
        self.refused(self.approve(BETA, quote), 8)
        log.write_bytes(original)
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.refused(self.approve(BETA, quote, manifest="0" * 64), 8)
        self.refused(self.approve(BETA, quote, manifest=None), 8)
        self.assertEqual(self.state(BETA), "Ready for approval")
        self.ok(self.approve(BETA, quote))
        self.assertIn(sha256_file(good / "MANIFEST.sha256"), json.dumps(self.status()))

    def fingerprint_script(self, fingerprint, code=0):
        self.counter += 1
        path = self.work / f"fingerprint-{self.counter}.py"
        path.write_text(FINGERPRINT_SCRIPT.format(python=sys.executable,
                                                  fingerprint=fingerprint, code=code))
        path.chmod(0o755)
        return path

    def test_install_check_availability(self):
        """Revision 4: a mismatching or failing fingerprint makes the installation
        unavailable (12), run refused (12), report says so; a later match restores it."""
        self.prepare(BETA)
        good = self.fingerprint_script("abc123")

        def available():
            self.ok(self.tool("install", "check", "--script", good, "--expect", "abc123"))
            shown = self.ok(self.tool("install", "show")).stdout.lower()
            self.assertIn("available", shown)
            self.assertNotIn("unavailable", shown)
            self.attempt(BETA, 0)
            self.assertNotIn("unavailable", self.report().lower())
        available()
        for bad in (self.fingerprint_script("other999"), self.fingerprint_script("abc123", 1)):
            with self.subTest(script=bad.name):
                result = self.tool("install", "check", "--script", bad, "--expect", "abc123")
                self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
                self.assertIn("unavailable", self.ok(self.tool("install", "show")).stdout.lower())
                self.attempt(BETA, 12)
                report = self.report().lower()
                self.assertIn("installation", report)
                self.assertIn("unavailable", report)
                available()

    def test_failing_probe_refused_not_started(self):
        """A probe that fails: run exit 7, record `not started`, CMD never starts."""
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("ran-after-failed-probe")
        probe = "echo " + shlex.quote(str(self.src)) + "; exit 4"
        self.refused(self.run_cmd(BETA, "unit", probe=probe, cmd=cmd), 7)
        self.assertFalse(marker.exists())
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "not started")
        self.attempt(BETA, 0)

    def test_wait_unknown_default_retry_budget(self):
        """wait --unknown: the retry budget defaults to 8 (visible in status)."""
        def budgets(value):
            if isinstance(value, dict):
                return [v for k, v in value.items() if "retr" in k.lower()] + [
                    x for v in value.values() for x in budgets(v)]
            if isinstance(value, list):
                return [x for v in value for x in budgets(v)]
            return []
        for job, extra, expected in ((ALPHA, (), 8), (BETA, ("--max-retries", "3"), 3)):
            with self.subTest(job=job):
                self.ok(self.set(job, "Preparing"))
                self.ok(self.tool("wait", job, "--unknown", *extra))
                self.assertIn(expected, budgets(job_entry(self.status(), job)))

    # ---- Revision 5 ---------------------------------------------------------------

    def easy_run(self, commands, exec_last=False):
        """Run CMD in the libtbx.easy_run pattern; returns the shell's pid. Revision 15: the
        shell is SHELL_STANDIN (a python process started with setsid, as easy_run starts
        /bin/sh) running `commands`, (argv list, platform) pairs, one after another; see
        SHELL_STANDIN for `platform` and `exec_last`."""
        script, pidfile = self.work / "easy_run.py", self.work / f"shell-{self.counter}.pid"
        script.write_text(EASY_RUN_SCRIPT)
        standin = self.work / "shell_standin.py"
        standin.write_text(SHELL_STANDIN)
        spec = {"commands": [{"argv": argv, "platform": platform} for argv, platform in commands],
                "exec_last": exec_last}
        self.ok(self.run_cmd(BETA, "unit", cmd=[sys.executable, str(script), str(pidfile),
                                                 json.dumps(spec), str(standin)]))
        shell = int(pidfile.read_text())
        self.stray.append(shell)
        return shell

    def easy_run_shape(self, shell, child, exec_last=False):
        """Revision 15 (positive control of the stand-in): the shell leads its own session and
        group, which is not the run's; `child` is its child in that session, or (exec_last)
        the shell itself after it replaced itself with its last command. On macOS neither
        shows a run-id marker to `ps -E`, as a platform /bin/sh and /bin/sleep do not."""
        self.assertEqual((os.getsid(shell), os.getpgid(shell)), (shell, shell))
        self.assertNotIn(shell, {record.get("pgid") for record in self.runs().values()})
        if exec_last:
            self.assertEqual(child, shell, "the shell did not replace itself with its last command")
        else:
            self.assertNotEqual(child, shell)
            self.assertEqual(os.getsid(child), shell)
            self.assertEqual(subprocess.run(["/bin/ps", "-o", "ppid=", "-p", str(child)],
                                            text=True, capture_output=True).stdout.strip(),
                             str(shell))
        if sys.platform == "darwin":
            for pid in {shell, child}:
                shown = subprocess.run(["/bin/ps", "-E", "-ww", "-o", "command=", "-p", str(pid)],
                                       text=True, capture_output=True).stdout
                self.assertTrue(shown.strip(), "ps -E listed nothing for %d" % pid)
                self.assertNotIn("GC_AUTO_RUN_ID=", shown, "a run-id marker is visible on %d" % pid)

    def held_child(self, seconds, name):
        """Revision 15: an argv list for a python that writes its pid to a file (returned)
        and stays `seconds` seconds (instead of /bin/sleep)."""
        pidfile = self.work / f"{name}-{self.counter}.pid"
        code = ("import os, sys, time\n"
                "open(sys.argv[1] + '.tmp', 'w').write(str(os.getpid()))\n"
                "os.replace(sys.argv[1] + '.tmp', sys.argv[1])\n"
                "time.sleep(%s)\n" % seconds)
        return [sys.executable, "-c", code, str(pidfile)], pidfile

    def test_easy_run_pattern_shell_survivor(self):
        """Revision 5 S2: a shell started with setsid by CMD (python, then a 41 s child)
        keeps the run from passing; report says running; stop signals the shell and
        the child; neither remains after Stopped. Revision 15: the shell is a python stand-in
        and the 41 s child (formerly /bin/sleep 41) a python that names its own pid; neither
        shows a run-id marker (as the platform /bin/sh and /bin/sleep did not), so on macOS
        too they are found through the wrapper's tracked descendants; the shell's session and
        group and the child's parentage are checked (easy_run_shape). The shell waits for
        its last command, as macOS /bin/sh does; see the exec_last variant below."""
        self.check_easy_run_survivor(exec_last=False)

    def test_easy_run_pattern_exec_last_shell_survivor(self):
        """Revision 15 (review finding 3): the test above with a shell that replaces itself
        with its last command (the 41 s child), as bash 5.1 and later and dash do: that one
        process, in the shell's new session, must keep the run from passing and be ended by
        stop."""
        self.check_easy_run_survivor(exec_last=True)

    def check_easy_run_survivor(self, exec_last):
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        child_argv, child_pidfile = self.held_child(41, "child41")
        shell = self.easy_run([(sleep_cmd(4), False), (child_argv, True)], exec_last)
        self.assertTrue(alive(shell))
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("running", cell.lower())
        self.assertNotIn("passed", cell.lower())
        self.refused(self.ready(BETA, "--packet", self.packet(self.tree())), 8)
        self.wait_for(child_pidfile.exists, "the shell's 41 s child")
        child = int(child_pidfile.read_text())
        self.stray.append(child)
        self.easy_run_shape(shell, child, exec_last)
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("running", cell.lower())
        self.assertNotIn("passed", cell.lower())
        self.end_worker()
        first = self.tool("stop")
        self.assertIn(first.returncode, (0, 10), first.stdout + first.stderr)
        self.stop_until_stopped()
        self.assertFalse(alive(shell), "the shell survived Stopped")
        self.assertFalse(alive(child), "the 41 s child survived Stopped")

    def test_easy_run_pattern_quick_shell_passes(self):
        """Revision 5 S2 control: the same pattern with a shell that ends quickly
        completes and passes."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        shell = self.easy_run([(QUICK_CMD, False), (QUICK_CMD, True)])
        self.wait_for(lambda: not alive(shell), "the quick shell to exit")
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "completed")
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))

    def forge(self, record, name, **changes):
        fake = {k: (name if v == record.get("run_id") else v) for k, v in record.items()}
        fake.update(changes)
        for key in [k for k, v in fake.items() if v is None]:
            del fake[key]
        staging = self.work / "forged.json"
        staging.write_text(json.dumps(fake))
        os.replace(staging, self.queue / "runs" / f"{name}.json")
        return self.queue / "runs" / f"{name}.json"

    def test_ended_run_does_not_claim_later_process(self):
        """Revision 5 S1, isolated; the orphan starts 2.5 s after `ended`."""
        self.check_ended_bound(lambda ended: time.time() + 2.5)

    def test_ended_run_bound_exact_to_the_second(self):
        """Revision 5 S1 (window exact to whole UTC seconds); the orphan starts just
        after the next whole UTC second following `ended` (Revision 12: waits until
        ps reports that second for a new process)."""
        self.check_ended_bound(lambda ended: ended + 1.05)

    def utc(self, iso):
        return calendar.timegm(time.strptime(iso, "%Y-%m-%dT%H:%M:%SZ"))

    def start_epoch(self, pid):
        lstart = subprocess.run(["/bin/ps", "-o", "lstart=", "-p", str(pid)],
                                env={"LC_ALL": "C", "TZ": "UTC"}, text=True,
                                capture_output=True).stdout.strip()
        return calendar.timegm(time.strptime(lstart, "%a %b %d %H:%M:%S %Y"))

    def wait_past_second(self, bound):
        """Revision 12: wait until a newly started process's ps start time is a later
        whole UTC second than `bound`. Linux ps derives start times from the boot time
        and clock ticks, so the wall clock passing bound + 1 is not enough."""
        def later():
            probe = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(5)"])
            try:
                return self.start_epoch(probe.pid) > bound
            finally:
                probe.kill()
                probe.wait(5)
        self.wait_for(later, f"a process start time later than {bound}")

    def orphan_of_dead_leader(self):
        """An exited, unrelated session leader P and its orphan (sid = pgid = P,
        reparented away from P) that notes SIGTERM. Returns (P, orphan pid, ready, orphan
        start). Revision 14: the new parent may be any process but P (launchd or init is
        process 1; a container or a child subreaper can be another process).
        Revision 15: P is a python process (formerly `/bin/sh -c 'CMD & echo $! > FILE; exit
        0'`) that starts the child in P's group and session, with stdin from /dev/null as
        sh gives a background job, writes the child's pid and exits 0 without waiting."""
        cmd, ready, _ = self.wait_cmd("note")
        pidfile = self.work / f"orphan-{self.counter}.pid"
        code = ("import json, os, subprocess, sys\n"
                "child = subprocess.Popen(json.loads(sys.argv[1]), stdin=subprocess.DEVNULL)\n"
                "with open(sys.argv[2], 'w') as handle:\n"
                "    handle.write('%d\\n' % child.pid)\n"
                "os._exit(0)\n")
        shell = subprocess.Popen([sys.executable, "-c", code, json.dumps(cmd), str(pidfile)],
                                 start_new_session=True)
        self.assertEqual(shell.wait(TIMEOUT), 0)
        leader = shell.pid
        orphan = int(pidfile.read_text())
        self.stray.append(orphan)
        start = self.start_epoch(orphan)
        self.wait_for(ready.exists, "the orphaned child")
        self.assertFalse(alive(leader))
        self.assertEqual((os.getpgid(orphan), os.getsid(orphan)), (leader, leader))
        self.wait_for(lambda: subprocess.run(["/bin/ps", "-o", "ppid=", "-p", str(orphan)],
                                             text=True, capture_output=True).stdout.strip()
                      not in ("", str(leader)),
                      "the orphan to be reparented away from its dead leader")
        return leader, orphan, ready, start

    def check_counted_only_when(self, record, excluded, included, orphan, ready):
        """The forged record with `excluded` changes claims nothing (no 13, Stopped,
        orphan not signalled); with `included` it counts the orphan (13, Stopping,
        SIGTERM)."""
        name = record["run_id"] + "-forged"
        self.forge(record, name, **excluded)
        self.assertNotIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 0)  # nothing live: no 13
        self.forge(record, name, **included)  # control
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 13)
        self.end_worker()
        self.forge(record, name, **excluded)
        result = self.ok(self.tool("stop"))
        self.assertTrue(result.stdout.strip().startswith("Stopped"), result.stdout)
        self.assertFalse(Path(str(ready) + ".term").exists(), "an unclaimed process was signalled")
        self.assertTrue(alive(orphan))
        self.forge(record, name, **included)  # control
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(ready) + ".term").exists, "the counted orphan to get SIGTERM", 5)

    def check_ended_bound(self, start_at):
        """Revision 5 S1, isolated: the record's leader P is dead (so its group and
        session are trusted); P's orphaned child (sid = pgid = P, reparented) started
        after the record's `ended` is not counted or signalled, and a new run is
        accepted. Control: the same record without `ended` counts it (13, Stopping,
        SIGTERM)."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        ended = self.utc(record["ended"])
        time.sleep(max(0.0, start_at(ended) - time.time()))
        self.wait_past_second(ended)
        leader, orphan, ready, start = self.orphan_of_dead_leader()
        self.assertGreater(start, ended)
        excluded = dict(pid=leader, pgid=leader, process_start=None)
        self.check_counted_only_when(record, excluded, dict(excluded, ended=None,
                                                            status="running"), orphan, ready)

    def test_probe_window_bounds_probe_leader(self):
        """Revision 6 S1 variant: a running record's probe_pgid (an exited, unrelated
        session leader) claims its orphan only if the orphan started within
        [probe_started, probe_ended]. Revision 12: the orphan starts once ps reports
        a later second than probe_ended for a new process."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        probe_ended = self.utc(record["probe_ended"])
        self.assertLessEqual(self.utc(record["probe_started"]), probe_ended)
        time.sleep(max(0.0, probe_ended + 1.05 - time.time()))
        self.wait_past_second(probe_ended)
        leader, orphan, ready, start = self.orphan_of_dead_leader()
        self.assertGreater(start, probe_ended)
        excluded = dict(status="running", ended=None, probe_pgid=leader, probe_start=None)
        included = dict(excluded, probe_ended=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                            time.gmtime(start)))
        self.check_counted_only_when(record, excluded, included, orphan, ready)

    def test_shell_parent_of_marked_process_not_signalled(self):
        """Revision 6 M3: a bash started by the user during a run, whose child carries
        the marker text only in its arguments, is never signalled by stop. Revision 15: the
        user's "bash" is a python parent (the test's child) that notes SIGTERM in a file, as
        the bash trap did, runs the marked child and waits for it; its own command line holds
        the marker text only inside a JSON word (`...=<run id>"]`), as the bash command line
        held it only as `...=<run id>;`."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("default")
        proc = self.popen_run(BETA, cmd)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the run")
        marked, marked_ready, marked_release = self.wait_cmd("default")
        trapped, done = self.work / "bash-got-term", self.work / "bash-done"
        code = ("import json, signal, subprocess, sys\n"
                "signal.signal(signal.SIGTERM, lambda *_: open(sys.argv[1], 'w').close())\n"
                "subprocess.call(json.loads(sys.argv[3]))\n"
                "open(sys.argv[2], 'w').close()\n")
        bash = subprocess.Popen([sys.executable, "-c", code, str(trapped), str(done),
                                 json.dumps(marked + [f"GC_AUTO_RUN_ID={stem}"])])  # parent: the test
        self.procs.append(bash)
        self.wait_for(marked_ready.exists, "the marked python")
        self.end_worker()
        first = self.tool("stop")
        self.assertIn(first.returncode, (0, 10), first.stdout + first.stderr)
        release.write_text("")
        proc.wait(TIMEOUT)
        marked_release.write_text("")
        self.assertEqual(bash.wait(TIMEOUT), 0)
        self.assertTrue(done.exists())
        self.assertFalse(trapped.exists(), "stop signalled the user's bash")
        self.stop_until_stopped()

    def check_term_while_lock_held(self, send_term):
        """A held probe; the test takes .queue.lock, releases the probe and holds the
        lock until the probe process has ended (bounded polling; checkpoint 3, before:
        a 1 s probe and a fixed 1.5 s hold) and 0.5 s more; optionally SIGTERM the
        wrapper while it waits for the lock."""
        self.prepare(BETA)
        probe, started, release = self.held_probe(f"n1-probe-{self.counter}")
        marker, cmd = self.marker_cmd(f"n1-cmd-{self.counter}")
        proc = self.popen_run(BETA, cmd, probe=probe)
        output = self.work / f"popen-{len(self.procs) - 1}.out"
        self.wait_for(started.exists, "the probe")
        stem, record = self.wait_record(set(), status="starting")
        self.wait_for(lambda: self.runs()[stem].get("probe_pgid"), "the probe to be recorded")
        probe_pid = self.runs()[stem]["probe_pgid"]
        with open(self.queue / ".queue.lock", "a+") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                release.write_text("")
                # the probe ends; the wrapper (which reaps it) then waits for the lock
                self.wait_for(lambda: not alive(probe_pid), "the probe to end")
                if send_term:
                    os.kill(proc.pid, signal.SIGTERM)
                time.sleep(0.5)
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)
        code = proc.wait(TIMEOUT)
        self.assertNotIn("Traceback", output.read_text())
        (record,) = self.runs().values()
        return code, record, marker

    def test_sigterm_while_waiting_for_queue_lock(self):
        """Revision 6 N1: SIGTERM while the wrapper waits for the queue lock after the
        probe: no traceback, exit not 1, record not started or interrupted.
        Control: the same lock hold without SIGTERM completes."""
        code, record, marker = self.check_term_while_lock_held(True)
        self.assertNotEqual(code, 1)
        self.assertIn(record["status"], ("not started", "interrupted"))
        if record["status"] == "not started":
            self.assertFalse(marker.exists())
        self.ok(self.tool("reconcile"))
        self.assertNotIn(self.runs()[record["run_id"]]["status"], ("running", "unknown"))

    def test_queue_lock_hold_without_signal_completes(self):
        """Revision 6 N1 control: a held queue lock only delays the run."""
        code, record, marker = self.check_term_while_lock_held(False)
        self.assertEqual(code, 0)
        self.assertEqual(record["status"], "completed")
        self.assertTrue(marker.exists())

    def test_failed_new_baseline_keeps_expected_value(self):
        """Revision 6 N4: a failing --new-baseline check gives 12 and keeps the earlier
        expected value; a check with that value is then accepted without it."""
        good = self.fingerprint_script("abc123")
        self.ok(self.tool("install", "check", "--script", good, "--expect", "abc123"))
        wrong = self.fingerprint_script("xyz000")
        result = self.tool("install", "check", "--script", wrong, "--expect", "zzz999",
                           "--new-baseline")
        self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.ok(self.tool("install", "check", "--script", good, "--expect", "abc123"))
        self.assertNotIn("unavailable", self.tool("install", "show").stdout.lower())

    def test_child_of_tracked_pid_is_counted(self):
        """Revision 5: a running record's tracked, live pid makes its children
        survivors (report running, run refused 13, stop signals them).
        Revision 14: the tracked parent ignores SIGTERM, so it stays live through the stop;
        before, the parent killed by the stop stayed counted only as an uncollected zombie
        (a child of this test process), which now counts as gone."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        self.end_worker()
        (stem, record), = self.runs().items()
        launcher, pids = self.work / "launch.py", self.work / "pids-tracked"
        launcher.write_text(LAUNCH_SCRIPT)
        cmd, child_ready, child_release = self.wait_cmd("note")
        ignore_term = ("import runpy, signal, sys; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
                       "sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name='__main__')")
        parent = subprocess.Popen([sys.executable, "-c", ignore_term, str(launcher), "hold",
                                   str(self.wait_script), str(pids),
                                   f"note,{child_ready},{child_release}"],
                                  start_new_session=True)
        self.procs.append(parent)
        self.wait_for(pids.exists, "the tracked parent's child")
        (child,) = [int(pid) for pid in pids.read_text().split()]
        self.stray.append(child)
        lstart = subprocess.run(["/bin/ps", "-o", "lstart=", "-p", str(parent.pid)],
                                env={"LC_ALL": "C", "TZ": "UTC"}, text=True,
                                capture_output=True).stdout.strip()
        start = calendar.timegm(time.strptime(lstart, "%a %b %d %H:%M:%S %Y"))
        self.forge(record, stem + "-tracked", ended=None, status="running", pgid=child,
                   tracked=[{"pid": parent.pid, "start": start}])
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        self.start_worker(BETA)
        self.attempt(BETA, 13)
        self.end_worker()
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(child_ready) + ".term").exists, "the counted child to get SIGTERM", 5)
        self.assertIsNone(parent.poll(), "the tracked parent ignores SIGTERM")

    def test_no_relaunch_while_earlier_run_live(self):
        """Revision 5 S4: run of a job refused (13) while an earlier run of it has a
        live survivor; accepted once the survivor is gone."""
        self.prepare(BETA)
        _, (child,), (release,) = self.launch("default")
        self.attempt(BETA, 13)
        release.write_text("")
        self.wait_for(lambda: not alive(child), "the survivor to exit")
        self.attempt(BETA, 0)

    def test_install_check_bad_scripts_and_baseline(self):
        """Revision 5 S3: missing, non-executable or value-less fingerprint scripts give
        12 (unavailable, no traceback); a different --expect is refused (2) unless
        --new-baseline."""
        good = self.fingerprint_script("abc123")
        self.ok(self.tool("install", "check", "--script", good, "--expect", "abc123"))
        plain = self.fingerprint_script("abc123")
        plain.chmod(0o644)
        for bad in (self.work / "missing-script", plain, self.fingerprint_script("")):
            with self.subTest(script=bad.name):
                result = self.tool("install", "check", "--script", bad, "--expect", "abc123")
                self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertIn("unavailable", self.ok(self.tool("install", "show")).stdout.lower())
                self.ok(self.tool("install", "check", "--script", good, "--expect", "abc123"))
                self.assertNotIn("unavailable", self.tool("install", "show").stdout.lower())
        other = self.fingerprint_script("zzz999")
        result = self.tool("install", "check", "--script", other, "--expect", "zzz999")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.ok(self.tool("install", "check", "--script", other, "--expect", "zzz999",
                          "--new-baseline"))
        self.assertNotIn("unavailable", self.tool("install", "show").stdout.lower())

    def test_leaving_ready_clears_packet_and_approval(self):
        """Revision 5 M1 / Revision 6 N3: Approved -> Preparing clears the current
        packet, manifest and approval but keeps them in binding_history (quote also
        in EVENTS.log); a corrected packet in a new directory can be marked ready and
        approved; Approved -> Discarded keeps that approval in the history too."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        first = self.packet(self.tree())
        self.ok(self.ready(BETA, "--packet", first))
        old_manifest = sha256_file(first / "MANIFEST.sha256")
        old_quote = f"first approval of {BETA}"
        self.ok(self.approve(BETA, old_quote))
        self.ok(self.set(BETA, "Preparing"))
        self.check_history(BETA, old_quote, old_manifest)
        self.start_worker(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        second = self.packet(self.tree(), log="unit: exit 0 (corrected packet)\n")
        self.assertNotEqual(sha256_file(second / "MANIFEST.sha256"), old_manifest)
        self.ok(self.ready(BETA, "--packet", second))
        self.refused(self.approve(BETA, f"approve {BETA}", manifest=old_manifest), 8)
        self.ok(self.approve(BETA, f"approve {BETA}"))
        new_manifest = sha256_file(second / "MANIFEST.sha256")
        self.assertIn(new_manifest, json.dumps(self.status()))
        self.ok(self.set(BETA, "Discarded"))
        self.check_history(BETA, f"approve {BETA}", new_manifest)
        self.check_history(BETA, old_quote, old_manifest)

    def check_history(self, job, quote, manifest):
        """The approval is out of the current binding but in binding_history and the
        quote is in EVENTS.log."""
        entry = job_entry(self.status(), job)
        history = json.dumps(entry["binding_history"])
        self.assertIn(quote, history)
        self.assertIn(manifest, history)
        current = json.dumps({k: v for k, v in entry.items() if k != "binding_history"})
        self.assertNotIn(quote, current)
        self.assertNotIn(manifest, current)
        self.assertIn(quote, (self.queue / "EVENTS.log").read_text())

    def test_command_start_failure_not_started(self):
        """Revision 5: a --log that cannot be opened (a directory): exit 2 REFUSED,
        record `not started`, CMD never starts; a normal run is the control."""
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("ran-bad-log")
        result = self.tool("run", BETA, "unit", "--log", self.work, "--probe",
                           "echo " + shlex.quote(str(self.src)), "--", *cmd)
        self.refused(result, 2)
        self.assertFalse(marker.exists())
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "not started")
        self.attempt(BETA, 0)

    # ---- Revision 7: provisional choices, cleanup --------------------------------

    def test_provisional_choices_block_approval_until_decided(self):
        """Revision 7: provisional choices P1, P2; Approved refused while one is
        pending (step says `decide P<N> first`); decide needs the job id and P<N>
        as whole tokens; a rejected choice does not block approval."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        notes = ["used a smaller test map", "skipped the slow regression"]
        quote = f"approve {BETA}"
        for number, note in enumerate(notes, 1):
            result = self.ok(self.tool("provisional", BETA, "--note", note))
            self.assertIn(f"P{number}", result.stdout)
            self.refused(self.approve(BETA, quote), 8)
            step = self.result_cell(self.report(), BETA, column=2)
            if number == 1:
                self.assertIn("decide P1 first", step)
            else:  # two pending: both named before `first`
                self.assertRegex(step, r"decide P1\b.*\bP2 first")
        for bad in ("adopt P1", f"adopt {BETA}", f"{BETA} adopt P12 only",
                    f"x{BETA} adopt P1"):
            with self.subTest(quote=bad):
                self.refused(self.tool("decide", BETA, "1", "--adopt", "--quote", bad), 8)
        adopt = f"I adopt P1 for {BETA}"
        self.ok(self.tool("decide", BETA, "1", "--adopt", "--quote", adopt))
        self.refused(self.tool("decide", BETA, "1", "--reject", "--quote", adopt), 8)
        missing = self.tool("decide", BETA, "9", "--adopt", "--quote", f"P9 {BETA}")
        self.assertEqual(missing.returncode, 2, missing.stdout + missing.stderr)
        self.refused(self.approve(BETA, quote), 8)  # P2 still pending
        self.assertIn("decide P2 first", self.result_cell(self.report(), BETA, column=2))
        reject = f"reject P2 of {BETA}"
        self.ok(self.tool("decide", BETA, "2", "--reject", "--quote", reject))
        self.ok(self.approve(BETA, quote))
        self.assertEqual(self.state(BETA), "Approved")
        events = (self.queue / "EVENTS.log").read_text()
        for text in (*notes, adopt, reject):
            self.assertIn(text, events)

    def cleanup(self, *codes):
        result = self.tool("cleanup", BETA)
        self.assertIn(result.returncode, codes, result.stdout + result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        return result

    def test_cleanup_ends_group_survivor_not_outsider(self):
        """Revision 7: cleanup signals a run's leftover group member (Cleaned, exit 0,
        EVENTS line); report no longer running; an unrelated process is untouched."""
        self.prepare(BETA)
        _, (child,), _ = self.launch("default")
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        outsider_cmd, outsider_ready, _ = self.wait_cmd("note")
        outsider = subprocess.Popen(outsider_cmd)  # same session as the test, not the run's
        self.procs.append(outsider)
        self.wait_for(outsider_ready.exists, "the outsider")
        result = self.cleanup(0)
        self.assertTrue(result.stdout.strip().startswith("Cleaned"), result.stdout)
        self.assertFalse(alive(child))
        self.assertIn(str(child), (self.queue / "EVENTS.log").read_text())
        self.assertNotIn("running", self.result_cell(self.report(), BETA).lower())
        self.assertFalse(Path(str(outsider_ready) + ".term").exists(), "cleanup hit an outsider")
        self.assertIsNone(outsider.poll())
        self.attempt(BETA, 0)  # no live earlier run any more

    def test_cleanup_ends_easy_run_shell(self):
        """Revision 7: cleanup ends the setsid shell of the easy_run pattern. Revision 15:
        python stand-ins for the shell and its 41 s child, without visible run-id markers (see
        test_easy_run_pattern_shell_survivor); the child is ended too."""
        self.check_cleanup_easy_run(exec_last=False)

    def test_cleanup_ends_easy_run_exec_last_shell(self):
        """Revision 15 (review finding 3): the test above with a shell that replaces itself
        with its last command, as bash 5.1 and later and dash do."""
        self.check_cleanup_easy_run(exec_last=True)

    def check_cleanup_easy_run(self, exec_last):
        self.prepare(BETA)
        child_argv, child_pidfile = self.held_child(41, "cleanup41")
        shell = self.easy_run([(sleep_cmd(4), False), (child_argv, True)], exec_last)
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        self.wait_for(child_pidfile.exists, "the shell's 41 s child")  # revision 15
        child = int(child_pidfile.read_text())
        self.stray.append(child)
        self.easy_run_shape(shell, child, exec_last)
        result = self.cleanup(0)
        self.assertTrue(result.stdout.strip().startswith("Cleaned"), result.stdout)
        self.assertFalse(alive(shell))
        self.assertFalse(alive(child))  # revision 15
        self.assertNotIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 0)

    def helper_snapshot(self):
        """Revision 12 (P9): the helper's own process table, pid -> [ppid, pgid, sid, start]."""
        script = self.work / "snapshot.py"
        script.write_text(SNAPSHOT_SCRIPT)
        result = subprocess.run([sys.executable, "-I", "-B", str(script), str(TOOL.parent)],
                                cwd=self.tmp, env=self.env, text=True, capture_output=True,
                                timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return {int(pid): entry for pid, entry in json.loads(result.stdout).items()}

    def test_helper_process_table_lists_test_process_and_cleanup_ends_it(self):
        """Revision 12 (P9): a test process started by a run's CMD appears in the helper's
        own process table with its pid, its parent (the CMD leader) and its group (the
        run's); after the wrapper and the leader are killed, cleanup ends it. An empty
        or wrong table fails here."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        launcher, pids = self.work / "launch.py", self.work / "pids-table"
        launcher.write_text(LAUNCH_SCRIPT)
        _, ready, release = self.wait_cmd("default")
        proc = self.popen_run(BETA, [sys.executable, str(launcher), "hold", str(self.wait_script),
                                     str(pids), f"default,{ready},{release}"])
        stem, record = self.wait_record(set())
        self.wait_for(pids.exists, "the test process's pid")
        (child,) = [int(pid) for pid in pids.read_text().split()]
        self.stray.append(child)
        self.wait_for(ready.exists, "the test process to start")
        leader = record["pid"]
        self.assertEqual(os.getpgid(child), record["pgid"])
        procs = self.helper_snapshot()
        self.assertIn(child, procs, "the helper's process table lacks the test process %d "
                      "(%d processes listed)" % (child, len(procs)))
        self.assertIn(leader, procs, "the helper's process table lacks the CMD leader")
        ppid, pgid, sid, start = procs[child]
        self.assertEqual(ppid, leader, "parent of the test process in the helper's table")
        self.assertEqual(pgid, record["pgid"], "group of the test process in the helper's table")
        self.assertEqual(sid, os.getsid(child))
        self.assertIsInstance(start, int)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(5)
        os.kill(leader, signal.SIGKILL)
        self.wait_for(lambda: not alive(leader), "the CMD leader to exit")
        self.assertTrue(alive(child))
        self.end_worker()
        result = self.cleanup(0)
        self.assertTrue(result.stdout.strip().startswith("Cleaned"), result.stdout)
        self.assertFalse(alive(child), "cleanup left the test process running")
        self.assertIn("process %d of run %s" % (child, record["run_id"]),
                      (self.queue / "EVENTS.log").read_text())

    def test_cleanup_refused_while_run_leader_alive(self):
        """Revision 7: cleanup gives 13 and signals nothing while a run of the job is
        running with its leader alive; afterwards it is accepted."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("note")
        proc = self.popen_run(BETA, cmd)
        self.wait_record(set())
        self.wait_for(ready.exists, "the run")
        self.refused(self.tool("cleanup", BETA), 13)
        self.assertFalse(Path(str(ready) + ".term").exists(), "cleanup signalled a live run")
        release.write_text("")
        self.assertEqual(proc.wait(TIMEOUT), 0)
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))

    def test_cleanup_still_alive_and_allowed_under_stop(self):
        """Revision 7: a survivor ignoring SIGTERM gives `Still alive` (exit 10);
        cleanup runs under a stop; Cleaned once the survivor is gone."""
        self.prepare(BETA)
        _, (child,), (release,) = self.launch("ignore")
        self.end_worker()
        result = self.cleanup(10)
        self.assertTrue(result.stdout.strip().startswith("Still alive"), result.stdout)
        self.assertTrue(alive(child))
        self.stopping(self.tool("stop"))
        result = self.cleanup(10)  # allowed under a stop: not 3
        self.assertTrue(result.stdout.strip().startswith("Still alive"), result.stdout)
        release.write_text("")
        self.wait_for(lambda: not alive(child), "the survivor to exit")
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))
        self.stop_until_stopped()

    def test_lock_release_refused_while_survivor_lives(self):
        """Revision 7: lock release gives 13 naming the run while a run of the queue
        has live processes; the lock stays; after cleanup it is released."""
        self.prepare(BETA)
        stem, (child,), _ = self.launch("default")
        result = self.tool("lock", "release")
        self.refused(result, 13)
        self.assertIn(stem, result.stdout + result.stderr)
        self.assertTrue((self.lock / "owner.txt").exists())
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))
        self.assertFalse(alive(child))
        self.ok(self.tool("lock", "release"))
        self.assertFalse(self.lock.exists())

    def test_blocked_rows_say_whether_a_decision_is_needed(self):
        """Revision 7: Blocked --decision TEXT puts `decision needed: TEXT` in the step;
        without it the step says `nothing to decide`; --decision on another
        transition is a usage error (2)."""
        self.ok(self.set(ALPHA, "Preparing"))
        self.ok(self.set(BETA, "Preparing"))
        result = self.set(BETA, "Testing", "--decision", "choose option a or b")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(self.state(BETA), "Preparing")
        self.ok(self.set(ALPHA, "Blocked", "--note", "two fixes possible",
                         "--decision", "choose option a or b"))
        self.ok(self.set(BETA, "Blocked", "--note", "disk was full"))
        report = self.report()
        self.assertIn("decision needed: choose option a or b",
                      self.result_cell(report, ALPHA, column=2))
        step = self.result_cell(report, BETA, column=2)
        self.assertIn("nothing to decide", step)
        self.assertNotIn("decision needed", step)

    def test_cleanup_refused_during_probe(self):
        """Re-read 7: cleanup gives 13 while a run of the job is in its probe; the
        probe is not killed and the run completes normally."""
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("ran-after-probe-cleanup")
        # Checkpoint 3: the probe is held until cleanup has answered (before: a 3 s probe
        # and cleanup after a fixed 1 s, which a slow start could push past the probe).
        probe, started, release = self.held_probe("cleanup-probe")
        proc = self.popen_run(BETA, cmd, probe=probe)
        stem, _ = self.wait_record(set(), status="starting")
        self.wait_for(started.exists, "the probe to start")
        self.refused(self.tool("cleanup", BETA), 13)
        self.assertIsNone(proc.poll(), "the run ended before the probe was released")
        release.write_text("")
        self.assertEqual(proc.wait(TIMEOUT), 0)
        self.assertEqual(self.runs()[stem]["status"], "completed")
        self.assertTrue(marker.exists())

    def test_published_refused_while_provisional_pending(self):
        """Re-read 7: an Approved job with a pending provisional choice cannot be
        Published (8); after the choice is decided it can."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        self.ok(self.ready(BETA, "--packet", self.packet(self.tree())))
        self.ok(self.approve(BETA, f"approve {BETA}"))
        self.assertIn("P1", self.ok(self.tool("provisional", BETA, "--note",
                                                "renamed an output file")).stdout)
        self.refused(self.set(BETA, "Published"), 8)
        self.assertEqual(self.state(BETA), "Approved")
        self.ok(self.tool("decide", BETA, "1", "--adopt", "--quote", f"adopt P1 for {BETA}"))
        self.ok(self.set(BETA, "Published"))
        self.assertEqual(self.state(BETA), "Published")

    # ---- Revision 8: time limit, outside network, worker periods ----------------

    def timeout_cmd(self, child_mode):
        """CMD: start a background child (`plain` or `ignore` SIGTERM), sleep 60."""
        script, pidfile = self.work / "timeout_cmd.py", self.work / f"child-{self.counter}.pid"
        script.write_text(TIMEOUT_SCRIPT)
        return [sys.executable, str(script), str(pidfile), child_mode], pidfile

    def check_timed_out(self, child_mode, limit):
        self.use_new_queue("q-timeout", extra=("--default-timeout", "3"))
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        cmd, pidfile = self.timeout_cmd(child_mode)
        began = time.time()
        result = self.run_cmd(BETA, "unit", cmd=cmd)
        elapsed = time.time() - began
        self.assertEqual(result.returncode, 124, result.stdout + result.stderr)
        self.assertLess(elapsed, limit)
        child = int(pidfile.read_text())
        self.stray.append(child)
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "timed out")
        self.assertEqual(record["timeout_seconds"], 3)
        self.assertFalse(record.get("left_after_timeout"), record)
        self.assertFalse(alive(child), "the run's background child was left")
        self.assertFalse(alive(record["pid"]))
        self.refused(self.ready(BETA, "--packet", self.packet(self.tree())), 8)
        self.assertNotIn("passed", self.result_cell(self.report(), BETA).lower())
        self.ok(self.tool("lock", "release"))  # nothing is left
        self.assertFalse(self.lock.exists())

    def test_default_timeout_ends_run_and_its_child(self):
        """Revision 8: with --default-timeout 3 a 60 s CMD with a background child exits
        124 within about 15 s, `timed out`, nothing left, never ready; lock release ok."""
        self.check_timed_out("plain", 15)

    def test_timeout_kills_child_ignoring_sigterm(self):
        """Revision 8: a child that ignores SIGTERM is ended by SIGKILL (about 10 s later)."""
        self.check_timed_out("ignore", 25)

    def test_timeout_option_overrides_default(self):
        """Revision 8: --timeout 30 overrides the default 3 (a 2 s CMD completes);
        --timeout 0 is a usage error (2)."""
        self.use_new_queue("q-timeout", extra=("--default-timeout", "3"))
        self.prepare(BETA)
        sleeper = [sys.executable, "-c", "import time; time.sleep(4)"]
        marker, cmd = self.marker_cmd("ran-timeout-zero")
        result = self.run_cmd(BETA, "unit", cmd=cmd, extra=("--timeout", "0"))
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertFalse(marker.exists())
        self.ok(self.run_cmd(BETA, "unit", cmd=[sys.executable, "-c",
                                                 "import time; time.sleep(2)"],
                             extra=("--timeout", "30")))
        self.assertEqual(self.run_cmd(BETA, "unit", cmd=sleeper).returncode, 124)  # control
        statuses = sorted(r["status"] for r in self.runs().values() if r["status"] != "not started")
        self.assertEqual(statuses, ["completed", "timed out"])

    def net_script(self):
        script = self.work / "net.py"
        script.write_text(NET_SCRIPT)
        return script

    # Re-read 10 (c): None until measured, then the outcome of the outside control.
    outside_control = None

    def outside_tcp_control(self):
        """Outside the block, a 5 s TCP connection to 192.0.2.1:443 (documentation
        range, no host; only a SYN is sent). Returns (eperm, description)."""
        if AutoQueueChecks.outside_control is None:
            try:
                socket.create_connection((BLOCKED_ADDRESS, 443), timeout=5).close()
                outcome = (False, "connected")
            except OSError as error:
                eperm = (getattr(error, "errno", None) == errno.EPERM
                         or "not permitted" in str(error).lower())
                outcome = (eperm, repr(error))
            AutoQueueChecks.outside_control = outcome
        return AutoQueueChecks.outside_control

    def blocked_path_setup(self):
        """Re-read 10 (b, c): skip when the surrounding environment already refuses
        the outside address with EPERM (then the block is not what is shown); remove
        the proxy variables (either case) from the runs' environment."""
        eperm, description = self.outside_tcp_control()
        if eperm:
            self.skipTest("blocked-path test skipped: the surrounding environment already "
                          f"blocks outbound traffic (TCP {BLOCKED_ADDRESS}:443 outside the "
                          f"block failed with 'not permitted': {description})")
        for key in list(self.env):
            if key.upper() in PROXY_VARIABLES:
                del self.env[key]
        self.assertFalse([k for k in self.env if k.upper() in PROXY_VARIABLES])

    def test_outside_block_control_is_not_eperm(self):
        """Re-read 10 (c): outside the block, TCP to 192.0.2.1:443 for 5 s does not
        fail with "not permitted" (it times out or is unreachable); otherwise the
        environment already blocks outbound traffic and this (like the blocked-path
        tests) is skipped saying so."""
        eperm, description = self.outside_tcp_control()
        if eperm:
            self.skipTest("the surrounding environment already blocks outbound traffic: "
                          + description)
        self.assertNotIn("not permitted", description.lower())

    def outside_network_opt_in(self, destinations):
        """Revision 10: skip (before any pre-check or contact) unless GC_TEST_NETWORK
        is exactly `1`; the reason names the destinations."""
        if os.environ.get("GC_TEST_NETWORK") != "1":
            self.skipTest(OPT_IN_PREFIX + destinations)

    def internet_or_skip(self, script):
        check = subprocess.run([sys.executable, str(script), "-", OUTSIDE_URL], text=True,
                               capture_output=True, timeout=30)
        if check.returncode:
            self.skipTest("no outside network from this machine (outside the sandbox): "
                          + (check.stdout + check.stderr).strip()[-200:])

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "checks the real sandbox-exec refusal of outside "
                     "connections (the test-only stand-in does not block the network)")
    def test_outside_network_blocked_by_default(self):
        """Revision 8 / Revision 10: a job without "network" runs CMD with outside
        network blocked: an HTTPS request to https://192.0.2.1/ and a TCP connection
        to 192.0.2.1:443 (documentation range, no host) fail with "not permitted";
        record network `blocked`; a localhost connection in the same kind of run works.
        Re-read 10 (b, c): proxy variables removed; skipped if EPERM already happens
        outside the block."""
        self.blocked_path_setup()
        self.prepare(BETA)
        script = self.net_script()
        connect = self.work / "connect.py"
        connect.write_text(CONNECT_SCRIPT)
        for cmd, pattern in (([sys.executable, str(script), "-", BLOCKED_URL], r"not permitted"),
                             ([sys.executable, str(connect), BLOCKED_ADDRESS], r"not permitted")):
            with self.subTest(cmd=cmd[1]):
                before = set(self.runs())
                result = self.run_cmd(BETA, "unit", cmd=cmd)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                (record,) = [r for stem, r in self.runs().items() if stem not in before]
                self.assertEqual(record["network"], "blocked")
                self.assertEqual(record["status"], "completed")
                self.assertNotEqual(record["exit_code"], 0)
                self.assertRegex(Path(record["log"]).read_text().lower(), pattern)
        # control: localhost is not blocked
        server = socket.socket()
        self.addCleanup(server.close)
        server.bind(("127.0.0.1", 0))
        server.listen(1)

        def serve():
            connection, _ = server.accept()
            with connection:
                connection.sendall(b"ok\n" if connection.recv(16) == b"hi" else b"no\n")
        thread = threading.Thread(target=serve, daemon=True)
        thread.start()
        local = self.work / "local.py"
        local.write_text(LOCAL_SCRIPT)
        self.ok(self.run_cmd(BETA, "unit", cmd=[sys.executable, str(local),
                                                 str(server.getsockname()[1])]))
        thread.join(TIMEOUT)
        self.assertTrue(all(r["network"] == "blocked" for r in self.runs().values()))

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "checks the real sandbox-exec refusal of the "
                     "probe's outside connection (the test-only stand-in does not block "
                     "the network)")
    def test_probe_network_blocked_by_default(self):
        """Revision 8 / Revision 10: the probe runs under the same block (its HTTPS
        request to https://192.0.2.1/ fails: exit 7, not started). Re-read 10 (a):
        the probe's own error file says "not permitted" (a probe that merely failed,
        for example by timing out, would not show the block); (b, c) as above."""
        self.blocked_path_setup()
        self.prepare(BETA)
        script = self.net_script()
        marker, cmd = self.marker_cmd("ran-after-net-probe")
        errors = self.work / "probe-error.txt"
        probe = (f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} "
                 f"{shlex.quote(str(self.src))} {shlex.quote(BLOCKED_URL)} "
                 f"{shlex.quote(str(errors))}")
        self.refused(self.run_cmd(BETA, "unit", probe=probe, cmd=cmd), 7)
        self.assertFalse(marker.exists())
        self.assertTrue(errors.is_file(), "the probe wrote no error file")
        self.assertIn("not permitted", errors.read_text().lower())
        (record,) = self.runs().values()
        self.assertEqual((record["status"], record["network"]), ("not started", "blocked"))

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "the allowed-job control is the companion of the "
                     "sandbox-exec block (off macOS the stand-in never blocks, so success there "
                     "says nothing about the production control)")
    def test_network_allowed_job(self):
        """Revision 8 control: a job with "network": true records `allowed`; the same
        request succeeds in CMD and in the probe (skipped without internet).
        Revision 10: opt-in (GC_TEST_NETWORK=1)."""
        self.outside_network_opt_in(
            f"HTTPS GET {OUTSIDE_URL} (CMD and probe), TCP {OUTSIDE_ADDRESS}:443")
        script = self.net_script()
        self.internet_or_skip(script)
        jobs = self.tmp / "JOBS-net.json"
        jobs.write_text(json.dumps([dict(JOBS[1], network=True)]))
        self.use_new_queue("q-net", jobs=jobs)
        self.prepare(BETA)
        result = self.run_cmd(BETA, "unit", cmd=[sys.executable, str(script), "-", OUTSIDE_URL])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        probe = (f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} "
                 f"{shlex.quote(str(self.src))} {shlex.quote(OUTSIDE_URL)}")
        self.ok(self.run_cmd(BETA, "unit", probe=probe))
        connect = self.work / "connect.py"
        connect.write_text(CONNECT_SCRIPT)
        self.ok(self.run_cmd(BETA, "unit", cmd=[sys.executable, str(connect), OUTSIDE_ADDRESS]))
        records = list(self.runs().values())
        self.assertEqual(len(records), 3)
        for record in records:
            self.assertEqual((record["status"], record["network"]), ("completed", "allowed"))

    def test_workers_lists_periods_without_overlap(self):
        """Revision 8: `workers` lists each registration with start and end in order,
        `running` while one is registered, and a final `overlap: none`."""
        self.start_worker(BETA)
        first = self.workers[str(self.queue)][1]
        shown = self.ok(self.tool("workers")).stdout.strip().splitlines()
        self.assertIn("running", [line for line in shown if first in line][0])
        self.end_worker()
        self.start_worker(ALPHA)
        second = self.workers[str(self.queue)][1]
        self.end_worker()
        lines = self.ok(self.tool("workers")).stdout.strip().splitlines()
        self.assertEqual(lines[-1].strip(), "overlap: none")
        rows = lines[:-1]
        self.assertEqual(len(rows), 2, lines)
        for row, (worker, job) in zip(rows, ((first, BETA), (second, ALPHA))):
            self.assertIn(worker, row)
            self.assertIn(job, row)
            self.assertIn("finished", row)
            self.assertNotIn("running", row)
            self.assertGreaterEqual(len(re.findall(r"\d{4}-\d\d-\d\dT\d\d:\d\d", row)), 2, row)

    def dns_cmd(self, host):
        script = self.work / "dns.py"
        script.write_text(DNS_SCRIPT)
        return [sys.executable, str(script), host]

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "checks the real sandbox-exec refusal of an outside "
                     "name lookup (off macOS the stand-in does not block, so a failed lookup "
                     "would only show that this host is offline)")
    def test_outside_name_lookup_blocked(self):
        """Re-read 8: without "network", an outside name lookup fails inside CMD while
        localhost resolves. Revision 10: opt-in (it would query if the block failed)."""
        self.outside_network_opt_in(f"DNS {OUTSIDE_NAME} (under the block)")
        self.prepare(BETA)
        result = self.run_cmd(BETA, "unit", cmd=self.dns_cmd(OUTSIDE_NAME))
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.ok(self.run_cmd(BETA, "unit", cmd=self.dns_cmd("localhost")))
        records = sorted(self.runs().values(), key=lambda r: r["exit_code"])
        self.assertEqual([r["exit_code"] == 0 for r in records], [True, False])
        self.assertTrue(all(r["network"] == "blocked" for r in records))

    @unittest.skipIf(NON_MACOS, MACOS_ONLY + "the allowed-job control is the companion of the "
                     "sandbox-exec block (off macOS the stand-in never blocks, so success there "
                     "says nothing about the production control)")
    def test_outside_name_lookup_allowed_job(self):
        """Re-read 8 control: a job with "network": true resolves example.com (skipped
        without internet). Revision 10: opt-in (GC_TEST_NETWORK=1)."""
        self.outside_network_opt_in(f"DNS {OUTSIDE_NAME}")
        check = subprocess.run(self.dns_cmd(OUTSIDE_NAME), text=True, capture_output=True,
                               timeout=30)
        if check.returncode:
            self.skipTest("example.com does not resolve outside the sandbox: "
                          + check.stderr.strip()[-200:])
        jobs = self.tmp / "JOBS-net.json"
        jobs.write_text(json.dumps([dict(JOBS[1], network=True)]))
        self.use_new_queue("q-net", jobs=jobs)
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit", cmd=self.dns_cmd(OUTSIDE_NAME)))
        (record,) = self.runs().values()
        self.assertEqual(record["network"], "allowed")

    def test_probe_timeout(self):
        """Re-read 8: a probe longer than the limit (3 s against a 30 s probe): exit
        124, `not started` with reason `probe timed out`, CMD never runs, no probe
        process left. Control: a short probe runs CMD. Checkpoint 3: the time bound is
        measured from the run record's own timestamps (probe_started to ended, whole
        UTC seconds): at most 3 s limit + 2.5 s of the helper's TERM/KILL waits + 2 s
        for rounding and scheduling = 7 s, so helper start-up under load no longer
        counts; the whole call must still end well before the probe's 30 s (25 s).
        Checkpoint 5 (review CP3 finding 3): the lower bound is 5 s: the 3 s limit
        plus the helper's unconditional 2.5 s of TERM/KILL waits is 5.5 s, which in
        whole seconds is never below 5; a probe ended early (0.2 s) gives 2-3 s.
        The leftover check waits at most 5 s for the killed group to be reaped, and
        uses separate -o flags (procps reads `-o pgid=,command=` as one column)."""
        self.use_new_queue("q-probe-timeout", extra=("--default-timeout", "3"))
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("ran-after-probe-timeout")
        began = time.time()
        # Revision 15: the shell's 30 s child is a python process instead of `sleep 30`.
        result = self.run_cmd(BETA, "unit", cmd=cmd,
                              probe=" ".join(shlex.quote(word) for word in sleep_cmd(30))
                              + "; echo " + shlex.quote(str(self.src)))
        self.assertEqual(result.returncode, 124, result.stdout + result.stderr)
        self.assertLess(time.time() - began, 25)
        self.assertFalse(marker.exists())
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "not started")
        self.assertIn("probe timed out", record.get("reason", ""))
        took = self.utc(record["ended"]) - self.utc(record["probe_started"])
        self.assertGreaterEqual(took, 5, record)  # 3 s limit + 2.5 s waits, less rounding
        self.assertLessEqual(took, 7, record)

        def leftovers():
            # Revision 14: a zombie member (exited, not collected) is not a leftover.
            table = subprocess.run(["/bin/ps", "-A", "-o", "pgid=", "-o", "stat=", "-o",
                                    "command="], text=True, capture_output=True)
            self.assertEqual(table.returncode, 0, table.stderr)
            return [line for line in table.stdout.splitlines()
                    if line.split(None, 2)[:1] == [str(record["probe_pgid"])]
                    and not line.split(None, 2)[1].startswith("Z")]
        self.wait_for(lambda: not leftovers(), "the probe's group to be gone", 5)
        self.attempt(BETA, 0)

    @unittest.skipIf(NON_MACOS, "duplicate off macOS: there "
                     "test_linux_ps_lag_probe_timeout_finds_leader_and_detached_child runs the "
                     "same body with the real ps (the imitation changes nothing off the Mac)")
    def test_probe_timeout_with_detached_child_holding_output(self):
        """Probe time-out (confirmation of re-read 7): a probe whose detached python
        (own session, output inherited) outlives it does not hold the wrapper: exit
        124 within about 15 s, `not started` / `probe timed out`, the detached python
        is ended (by the run, or by cleanup before lock release is accepted).
        Checkpoint 3: on macOS this finds the detached python by its `ps -E` run-id
        marker, while the P10 variant drops the markers and finds it as the probe
        leader's child under imitated early Linux start times, so both run here. Off
        macOS the P10 variant is this same test with the real (early-reporting) ps,
        aligned to a whole second, so this one is skipped there."""
        self.check_probe_timeout_detached()

    def check_probe_timeout_detached(self, align=False):
        """The body of the probe time-out test; with `align` (Revision 12, P10) the run
        starts just after a whole UTC second and the record is returned."""
        self.prepare(BETA)
        script, pidfile = self.work / "detach_probe.py", self.work / "detached.pid"
        script.write_text(DETACH_PROBE)
        marker, cmd = self.marker_cmd("ran-after-detached-probe")
        probe = f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} {shlex.quote(str(pidfile))}"
        if align:
            time.sleep(1.02 - time.time() % 1)
        began = time.time()
        result = self.run_cmd(BETA, "unit", probe=probe, cmd=cmd, extra=("--timeout", "3"))
        elapsed = time.time() - began
        self.assertEqual(result.returncode, 124, result.stdout + result.stderr)
        self.assertLess(elapsed, 15)
        self.assertFalse(marker.exists())
        detached = int(pidfile.read_text())
        self.stray.append(detached)
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "not started")
        self.assertIn("probe timed out", record.get("reason", ""))
        release = self.tool("lock", "release")
        if release.returncode:
            self.refused(release, 13)
            self.assertTrue(self.lock.exists())
            self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))
            self.ok(self.tool("lock", "release"))
        self.assertFalse(alive(detached), "the detached python was left")
        self.assertFalse(self.lock.exists())
        return record

    # ---- Revision 12 (P10): one second allowed for early Linux ps start times ------

    def test_linux_ps_lag_probe_timeout_finds_leader_and_detached_child(self):
        """P10 (i): with ps start times reported up to 1 s early (Linux procps; imitated on
        the Mac by LINUX_PS_SCRIPT, markers dropped), a run started just after a whole
        second has its probe leader reported in the second before `started`. The helper
        still counts the leader (one second allowed), so it finds the leader's detached
        python (own session, reachable off macOS only as the leader's child) and ends it,
        at the time-out or by cleanup before lock release is accepted."""
        self.linux_ps = True
        record = self.check_probe_timeout_detached(align=True)
        if sys.platform == "darwin":  # the imitated case: the leader is reported early
            leader_start = calendar.timegm(time.strptime(record["probe_start"],
                                                         "%a %b %d %H:%M:%S %Y"))
            self.assertLess(leader_start, self.utc(record["started"]),
                            "the probe leader was not reported before `started`: %r" % record)

    def test_linux_ps_lag_earlier_process_not_claimed(self):
        """P10 (ii): with Linux ps start times imitated, an orphan of a dead group and
        session leader (sid = pgid = leader, reparented) whose record says the run started
        2 s after the orphan's start second (it started more than the allowance before
        the run) is not counted and not signalled by stop. Control: the same record with
        `started` at the orphan's start second counts it (13, Stopping, SIGTERM)."""
        self.linux_ps = True
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        leader, orphan, ready, start = self.orphan_of_dead_leader()
        iso = lambda second: time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(second))
        excluded = dict(pid=leader, pgid=leader, process_start=None, ended=None,
                        status="running", started=iso(start + 2))
        self.check_counted_only_when(record, excluded, dict(excluded, started=iso(start)),
                                     orphan, ready)

    def test_linux_ps_lag_cleanup_ends_survivor_not_outsider(self):
        """P10 (ii): with Linux ps start times imitated, cleanup still ends the run's
        group survivor and leaves an unrelated process (the test's own child, in the
        test's group and session) unsignalled."""
        self.linux_ps = True
        self.test_cleanup_ends_group_survivor_not_outsider()

    def test_linux_ps_lag_stop_spares_unrelated_process(self):
        """P10 (ii): with Linux ps start times imitated, a process in the test's own
        session that the run did not start, started while the run ran, is neither
        counted nor signalled by stop."""
        self.linux_ps = True
        self.test_unrelated_process_in_same_session_not_counted()

    # ---- Revision 9: screening record --------------------------------------------

    def screened_beta(self, network=False):
        """BETA Testing with a passing current unit run (a network job in its own
        queue when `network`); returns a function that freezes a packet with the
        given SCREENING.txt text (None: no file)."""
        if network:
            jobs = self.tmp / "JOBS-screen-net.json"
            jobs.write_text(json.dumps([dict(JOBS[1], network=True)]))
            self.use_new_queue("q-screen-net", jobs=jobs)
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        tree = self.tree()
        return lambda text: self.packet(tree, screening=text)

    def not_ready(self, packet, reason):
        """Ready for approval refused (8) with `not ready: ` and `reason`; state kept."""
        result = self.ready(BETA, "--packet", packet)
        self.refused(result, 8)
        self.assertIn("not ready: ", result.stderr)
        self.assertIn(reason, result.stderr)
        self.assertEqual(self.state(BETA), "Testing")

    def test_screening_record_missing_refused(self):
        """Revision 9 (Q1): a frozen packet without SCREENING.txt is not ready
        (`packet has no SCREENING.txt`); control: the same packet with a valid record."""
        frozen = self.screened_beta()
        self.not_ready(frozen(None), "packet has no SCREENING.txt")
        self.ok(self.ready(BETA, "--packet", frozen(SCREENING)))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))

    def test_screening_malformed_line_refused_with_line_number(self):
        """Revision 9 (Q1): a line that is not `test: NAME | outside: WHAT | DECISION`
        (NAME, WHAT non-empty; DECISION included or excluded) is refused as
        `SCREENING.txt line N malformed`, N counting every line from 1 (comments and
        blank lines included)."""
        frozen = self.screened_beta()
        prefix = "# screening of the tests\n\ntest: unit | outside: none | included\n"
        for bad in ("unit | outside: none | included",
                    "test: unit | none | included",
                    "test: unit | outside: none",
                    "test: unit | outside: none | included | extra",
                    "test:  | outside: none | included",
                    "test: unit | outside:  | included",
                    "test: unit | outside: none | ",
                    "test: unit | outside: none | maybe",
                    "test: unit, outside: none, included",
                    "just some words"):
            with self.subTest(line=bad):
                self.not_ready(frozen(prefix + bad + "\n"), "SCREENING.txt line 4 malformed")
        with self.subTest(position="first line"):
            self.not_ready(frozen("bogus\n" + SCREENING), "SCREENING.txt line 1 malformed")
        with self.subTest(position="after several comments and blank lines"):
            text = "# a\n\n# b\n   \n" + SCREENING + "\n# c\nbogus\n"
            self.not_ready(frozen(text), "SCREENING.txt line 8 malformed")
        self.ok(self.ready(BETA, "--packet", frozen(prefix)))  # control: same prefix only

    def test_screening_needs_an_included_test(self):
        """Revision 9 (Q1): a record with no `included` line is refused (`SCREENING.txt
        lists no included test`), also when it holds only comments and blank lines
        or is empty; control: one included line among excluded ones."""
        frozen = self.screened_beta()
        excluded = "test: slow | outside: none | excluded\n"
        for label, text in (("only excluded", excluded),
                            ("only comments and blanks", "# nothing screened\n\n   \n"),
                            ("empty", "")):
            with self.subTest(record=label):
                self.not_ready(frozen(text), "SCREENING.txt lists no included test")
        self.ok(self.ready(BETA, "--packet", frozen(excluded + SCREENING)))

    def test_screening_comments_and_blank_lines_ignored(self):
        """Revision 9 (Q1): blank lines and lines starting `#` are ignored (even when
        they would not have the record form)."""
        frozen = self.screened_beta()
        text = ("# test: broken | | line that would be malformed\n\n   \n"
                "#no space after the hash\n" + SCREENING + "\n# trailing comment\n")
        self.ok(self.ready(BETA, "--packet", frozen(text)))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))

    def test_screening_included_outside_access_needs_network_grant(self):
        """Revision 9 (Q1): an included test with outside access other than `none`
        is refused when the job has no "network": true (`SCREENING.txt includes NAME
        with outside access but the job has no network grant`); the same test
        excluded is accepted."""
        frozen = self.screened_beta()
        outside = "test: fetch_pdb | outside: HTTPS download from rcsb.org | {}\n"
        self.not_ready(frozen(SCREENING + outside.format("included")),
                       "SCREENING.txt includes fetch_pdb with outside access but the job "
                       "has no network grant")
        self.ok(self.ready(BETA, "--packet", frozen(SCREENING + outside.format("excluded"))))

    def test_screening_included_outside_access_accepted_with_network_grant(self):
        """Revision 9 (Q1) control: the same included outside-access line is accepted
        for a job with "network": true, even as the only included test."""
        frozen = self.screened_beta(network=True)
        outside = "test: fetch_pdb | outside: HTTPS download from rcsb.org | included\n"
        self.ok(self.ready(BETA, "--packet", frozen(SCREENING + outside)))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.ok(self.approve(BETA, f"approve {BETA}"))

    def test_screening_rechecked_at_approved_and_in_report(self):
        """Revision 4/9: readiness is re-run at Approved and at report time. Deleting
        SCREENING.txt after Ready makes the frozen packet fail verification first (a
        listed file is missing), so this shows the packet re-check, not the screening
        rule alone (see test_screening_rechecked_at_approved_with_verifying_packet):
        Approved is refused (8) and the row shows `results stale or incomplete
        (<reason>)`, never `passed`; restored, the packet verifies again and Approved
        is accepted."""
        frozen = self.screened_beta()
        packet = frozen(SCREENING)
        self.ok(self.ready(BETA, "--packet", packet))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        screening = packet / "SCREENING.txt"
        screening.unlink()
        cell = self.result_cell(self.report(), BETA)
        self.assertNotIn("passed", cell.lower())
        self.assertRegex(cell, r"results stale or incomplete \(\S.*\)")
        quote = f"approve {BETA}"
        result = self.approve(BETA, quote)
        self.refused(result, 8)
        self.assertEqual(self.state(BETA), "Ready for approval")
        screening.write_text(SCREENING)
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.ok(self.approve(BETA, quote))
        self.assertEqual(self.state(BETA), "Approved")

    def set_network_in_record(self, job, value):
        """Rewrite the job's `network` field in QUEUE.json by hand (atomically), as a
        record the tool must re-check; the checker trusts the record files."""
        path = self.queue / "QUEUE.json"
        data = json.loads(path.read_text())
        entry = job_entry(data, job)
        self.assertIs(entry.get("network"), not value, entry)
        entry["network"] = value
        temp = path.with_name("QUEUE.json.hand")
        temp.write_text(json.dumps(data))
        os.replace(temp, path)

    def test_screening_rechecked_at_approved_with_verifying_packet(self):
        """Revision 9 (Q1): the screening check is re-run at Approved and in the
        report while the packet still verifies. A network job is Ready with an
        included outside-access test; with its grant removed from the record the
        row shows `results stale or incomplete (SCREENING.txt includes fetch_pdb
        with outside access but the job has no network grant)` and Approved is
        refused (8) with that reason; with the grant restored Approved is accepted."""
        frozen = self.screened_beta(network=True)
        packet = frozen(SCREENING + "test: fetch_pdb | outside: HTTPS to rcsb.org | included\n")
        self.ok(self.ready(BETA, "--packet", packet))
        verify = subprocess.run([sys.executable, "-I", "-B", str(SCREEN_CHECK), "verify",
                                 str(packet)], text=True, capture_output=True)
        self.assertEqual(verify.returncode, 0, verify.stdout + verify.stderr)
        reason = ("SCREENING.txt includes fetch_pdb with outside access but the job has no "
                  "network grant")
        self.set_network_in_record(BETA, False)
        cell = self.result_cell(self.report(), BETA)
        self.assertIn(f"results stale or incomplete ({reason})", cell)
        self.assertNotIn("passed", cell.lower())
        quote = f"approve {BETA}"
        result = self.approve(BETA, quote)
        self.refused(result, 8)
        self.assertIn(reason, result.stderr)
        self.assertEqual(self.state(BETA), "Ready for approval")
        self.set_network_in_record(BETA, True)
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.ok(self.approve(BETA, quote))
        self.assertEqual(self.state(BETA), "Approved")

    def test_report_row_stale_or_incomplete_shows_reason(self):
        """Revision 9 (Q1): a Ready row whose readiness no longer holds shows
        `results stale or incomplete (<reason>)`, never `passed` (here a run became
        stale after a criterion change)."""
        frozen = self.screened_beta()
        self.ok(self.ready(BETA, "--packet", frozen(SCREENING)))
        self.criterion.write_text("Criterion two: changed after readiness.\n")
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        cell = self.result_cell(self.report(), BETA)
        self.assertRegex(cell, r"results stale or incomplete \(\S.*\)")
        self.assertNotIn("passed", cell.lower())

    # ---- Revision 9 clarification and re-read 8 ----------------------------------

    def test_screening_keywords_are_case_sensitive(self):
        """Clarification: `test:`, `outside:`, `included`, `excluded` are case-sensitive
        (another case is malformed); `None` is not `none`, so it is outside access."""
        frozen = self.screened_beta()
        for bad in ("Test: unit | outside: none | included",
                    "TEST: unit | outside: none | included",
                    "test: unit | Outside: none | included",
                    "test: unit | outside: none | Included",
                    "test: unit | outside: none | EXCLUDED"):
            with self.subTest(line=bad):
                self.not_ready(frozen(SCREENING + bad + "\n"), "SCREENING.txt line 2 malformed")
        for what in ("None", "NONE"):
            with self.subTest(what=what):
                self.not_ready(frozen(f"test: other | outside: {what} | included\n"),
                               "SCREENING.txt includes other with outside access but the job "
                               "has no network grant")
        self.ok(self.ready(BETA, "--packet", frozen(
            SCREENING + "test: other | outside: None | excluded\n")))

    def screening_refused_either(self, packet, *reasons):
        """Fail-safe: Ready refused (8, `not ready: `) with one of `reasons`."""
        result = self.ready(BETA, "--packet", packet)
        self.refused(result, 8)
        self.assertIn("not ready: ", result.stderr)
        self.assertTrue(any(reason in result.stderr for reason in reasons), result.stderr)
        self.assertEqual(self.state(BETA), "Testing")

    def test_screening_whitespace_and_encoding_edges(self):
        """Revision 9 (Q1) edges: a trailing space after the decision is malformed; a
        UTF-8 BOM before line 1 is malformed; an indented `#` line is not a comment
        (malformed); `none ` (extra space before ` |`) is not `none` (fail-safe:
        refused as outside access or malformed); a non-UTF-8 byte inside NAME is still
        a NAME (accepted)."""
        frozen = self.screened_beta()
        with self.subTest(edge="trailing space"):
            self.not_ready(frozen("test: unit | outside: none | included \n"),
                           "SCREENING.txt line 1 malformed")
        with self.subTest(edge="BOM on line 1"):
            self.not_ready(frozen("﻿" + SCREENING), "SCREENING.txt line 1 malformed")
        with self.subTest(edge="indented comment"):
            self.not_ready(frozen(SCREENING + "  # indented note\n"),
                           "SCREENING.txt line 2 malformed")
        with self.subTest(edge="none with an extra space"):
            self.screening_refused_either(
                frozen("test: unit | outside: none  | included\n"),
                "SCREENING.txt includes unit with outside access but the job has no network grant",
                "SCREENING.txt line 1 malformed")
        with self.subTest(edge="non-UTF-8 byte in NAME"):
            packet = frozen(None)
            (packet / "SCREENING.txt").write_bytes(b"test: unit\xff | outside: none | included\n")
            self.refreeze(packet)
            self.ok(self.ready(BETA, "--packet", packet))

    def refreeze(self, packet):
        manifest = packet / "MANIFEST.sha256"
        if manifest.exists():
            manifest.unlink()
        frozen = subprocess.run([sys.executable, "-I", "-B", str(SCREEN_CHECK), "freeze",
                                 str(packet)], text=True, capture_output=True)
        self.assertEqual(frozen.returncode, 0, frozen.stdout + frozen.stderr)

    def test_screening_crlf_line_endings_accepted(self):
        """Revision 9 (Q1): a record written with CRLF line endings is read line by
        line as usual (accepted, as the test writer was asked to check)."""
        frozen = self.screened_beta()
        packet = frozen(None)
        (packet / "SCREENING.txt").write_bytes(
            b"# screening\r\n\r\ntest: unit | outside: none | included\r\n")
        self.refreeze(packet)
        self.ok(self.ready(BETA, "--packet", packet))

    def test_screening_separator_inside_name_or_what(self):
        """Revision 9 (Q1): ` | ` inside NAME (`test: a | b | outside: none | included`)
        is part of NAME under the pattern reading (accepted); ` | ` inside WHAT
        (`outside: web | proxy | included`) is fail-safe: refused, as outside access
        without a grant or as malformed."""
        frozen = self.screened_beta()
        with self.subTest(where="WHAT"):
            self.screening_refused_either(
                frozen("test: unit | outside: web | proxy | included\n"),
                "with outside access but the job has no network grant",
                "SCREENING.txt line 1 malformed")
        with self.subTest(where="NAME"):
            self.ok(self.ready(BETA, "--packet", frozen("test: a | b | outside: none | included\n")))

    def test_screening_txt_as_empty_directory_refused(self):
        """Revision 9 (Q1): SCREENING.txt that is an empty directory is refused (8):
        the packet does not verify or has no SCREENING.txt file."""
        frozen = self.screened_beta()
        packet = frozen(None)
        (packet / "SCREENING.txt").mkdir()
        result = self.ready(BETA, "--packet", packet)
        self.refused(result, 8)
        self.assertIn("not ready: ", result.stderr)
        self.assertTrue(any(reason in result.stderr for reason in
                            ("packet does not verify", "packet has no SCREENING.txt")),
                        result.stderr)
        self.ok(self.ready(BETA, "--packet", frozen(SCREENING)))  # control

    def test_report_step_for_packet_and_screening_failures(self):
        """After re-read 8 (a): a Ready row failing with a reason that begins `packet`
        or `SCREENING.txt` has the step "freeze a corrected packet in a new directory
        and mark it ready again" (Revision 13: with the clerical-slip note, same evidence,
        ERRATA.md, no rerun); another failure (a stale run) keeps "rerun the affected checks"."""
        freeze_step = ("freeze a corrected packet in a new directory and mark it ready again "
                       "(clerical slip: same evidence, note in ERRATA.md, no rerun)")
        rerun_step = "rerun the affected checks"
        frozen = self.screened_beta(network=True)
        packet = frozen(SCREENING + "test: fetch_pdb | outside: HTTPS to rcsb.org | included\n")
        self.ok(self.ready(BETA, "--packet", packet))
        log = packet / "RUN_LOG.txt"
        original = log.read_bytes()
        log.write_text("changed after Ready\n")  # reason begins `packet does not verify`
        report = self.report()
        self.assertIn("results stale or incomplete (packet", self.result_cell(report, BETA))
        self.assertEqual(self.result_cell(report, BETA, column=2), freeze_step)
        log.write_bytes(original)
        self.set_network_in_record(BETA, False)  # reason begins `SCREENING.txt`
        report = self.report()
        self.assertIn("results stale or incomplete (SCREENING.txt", self.result_cell(report, BETA))
        self.assertEqual(self.result_cell(report, BETA, column=2), freeze_step)
        self.set_network_in_record(BETA, True)
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.criterion.write_text("Criterion two: changed after readiness.\n")
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        report = self.report()
        cell = self.result_cell(report, BETA)
        self.assertIn("results stale or incomplete (", cell)
        self.assertNotIn("(packet", cell)
        self.assertNotIn("(SCREENING.txt", cell)
        step = self.result_cell(report, BETA, column=2)
        self.assertIn(rerun_step, step)
        self.assertNotIn("freeze a corrected packet", step)

    def test_init_refuses_non_boolean_network(self):
        """After re-read 8 (b): init refuses (exit 2, `job ID: network must be true or
        false`) a `network` that is present and not a JSON boolean; nothing is
        created; true, false and an absent field are accepted."""
        for index, value in enumerate(("false", "true", 0, 1, None, [], {}, "yes")):
            with self.subTest(network=value):
                jobs = self.tmp / f"JOBS-bad-net-{index}.json"
                jobs.write_text(json.dumps([JOBS[0], dict(JOBS[1], network=value)]))
                home = self.tmp / f"home-bad-net-{index}"
                home.mkdir()
                qid = f"q-bad-net-{index}"
                result = self.tool(*self.init_args(qid, self.tmp / f"RUN-bad-{index}.lock",
                                                   jobs=jobs), home=home, queue=None)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertIn(f"job {BETA}: network must be true or false", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertFalse((self.records / qid).exists())
                self.assertFalse((home / "ACTIVE").exists())
        for index, value in enumerate((True, False)):
            with self.subTest(accepted=value):
                jobs = self.tmp / f"JOBS-good-net-{index}.json"
                jobs.write_text(json.dumps([JOBS[0], dict(JOBS[1], network=value)]))
                self.use_new_queue(f"q-good-net-{index}", jobs=jobs)
                self.assertIs(job_entry(self.status(), BETA).get("network"), value)


    def test_init_refuses_bad_job_ids(self):
        """After test round 2: init refuses (exit 2, "job id '<id>' must not start with
        '-' or contain spaces or '|'") a job id that is not a string, starts with `-`,
        or contains whitespace or `|`; nothing is created. Ids such as `a.b+`, `Job`
        and `a-b` are accepted."""
        tail = "must not start with '-' or contain spaces or '|'"
        for index, bad in enumerate(("---", "-a", "my job", "a|b", "a\tb", " a", "a\n",
                                     7, None, ["a"])):
            with self.subTest(job_id=bad):
                jobs = self.tmp / f"JOBS-bad-id-{index}.json"
                jobs.write_text(json.dumps([JOBS[1], dict(JOBS[0], id=bad)]))
                home = self.tmp / f"home-bad-id-{index}"
                home.mkdir()
                qid = f"q-bad-id-{index}"
                result = self.tool(*self.init_args(qid, self.tmp / f"RUN-bad-id-{index}.lock",
                                                   jobs=jobs), home=home, queue=None)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertIn(tail, result.stderr)
                self.assertIn("job id ", result.stderr)
                # How a non-string or unprintable id is shown is not stated; a printable
                # string id is shown in quotes as in the spec's message.
                if isinstance(bad, str) and bad.isprintable():
                    self.assertIn(f"job id '{bad}' {tail}", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertFalse((self.records / qid).exists())
                self.assertFalse((home / "ACTIVE").exists())
        for index, good in enumerate(("a.b+", "Job", "a-b", "x_1")):
            with self.subTest(accepted=good):
                jobs = self.tmp / f"JOBS-good-id-{index}.json"
                jobs.write_text(json.dumps([dict(JOBS[0], id=good)]))
                self.use_new_queue(f"q-good-id-{index}", jobs=jobs)
                self.assertEqual(job_entry(self.status(), good)["state"], "Queued")


    def test_outside_tests_skipped_without_opt_in(self):
        """Revision 10: with GC_TEST_NETWORK unset (or not exactly `1`) each opt-in
        test is skipped with the stated prefix and its destinations, before any
        pre-check or contact: setUp is replaced by a no-op and every way to start a
        process, connect, resolve or fetch raises, so anything run before the skip
        decision would show as an error instead of a skip. Re-read 10 (d) control:
        with the value `1` each opt-in test goes on and reaches an intercepted call
        (an error or failure naming the interception, not a skip). Checkpoint 5: off
        macOS (real or simulated) a test that is also macOS-only (checkpoint 3) is
        skipped for that reason first, for every value including `1`, and is checked
        for the macOS-mechanism reason instead; the opt-in rule is checked for every
        test that is not skipped as macOS-only (on macOS: all three)."""
        expected = {"test_network_allowed_job": (OUTSIDE_URL, f"{OUTSIDE_ADDRESS}:443"),
                    "test_outside_name_lookup_allowed_job": (OUTSIDE_NAME,),
                    "test_outside_name_lookup_blocked": (OUTSIDE_NAME,)}

        macos_only = {name for name in expected
                      if getattr(getattr(AutoQueueChecks, name), "__unittest_skip__", False)}
        self.assertEqual(macos_only, set(expected) if NON_MACOS else set())
        opt_in = {name: where for name, where in expected.items() if name not in macos_only}
        scratch = self.tmp / "gated"
        scratch.mkdir()

        class Gated(AutoQueueChecks):
            def setUp(self):
                # Only what the tests touch before their first process or contact.
                self.tmp = self.work = scratch
                self.home, self.queue, self.lock = scratch / "home", scratch / "queue", scratch / "L"
                self.env, self.counter, self.workers, self.procs = {}, 0, {}, []
                self.criterion = scratch / "criterion.txt"

        marker = "intercepted call"

        def forbidden(*args, **kwargs):
            raise AssertionError(f"{marker}: {args!r}")
        import unittest.mock
        import urllib.request
        for value in (None, "0", "true", "1 ", "", "1"):
            with self.subTest(GC_TEST_NETWORK=value):
                environment = {k: v for k, v in os.environ.items() if k != "GC_TEST_NETWORK"}
                if value is not None:
                    environment["GC_TEST_NETWORK"] = value
                suite = unittest.TestSuite(Gated(name) for name in expected)
                result = unittest.TestResult()
                with unittest.mock.patch.dict(os.environ, environment, clear=True), \
                        unittest.mock.patch.object(subprocess, "run", forbidden), \
                        unittest.mock.patch.object(subprocess, "Popen", forbidden), \
                        unittest.mock.patch.object(socket, "create_connection", forbidden), \
                        unittest.mock.patch.object(socket, "getaddrinfo", forbidden), \
                        unittest.mock.patch.object(urllib.request, "urlopen", forbidden):
                    suite.run(result)
                reasons = {test._testMethodName: reason for test, reason in result.skipped}
                self.assertEqual(result.testsRun, len(expected))
                for name in macos_only:
                    self.assertTrue(reasons.get(name, "").startswith(MACOS_ONLY), reasons.get(name))
                if value == "1":
                    self.assertEqual(sorted(reasons), sorted(macos_only))
                    reached = {test._testMethodName: trace
                               for test, trace in result.errors + result.failures}
                    self.assertEqual(sorted(reached), sorted(opt_in))
                    for name, trace in reached.items():
                        self.assertIn(marker, trace, name)
                    continue
                self.assertEqual((result.errors, result.failures), ([], []))
                self.assertEqual(sorted(reasons), sorted(expected))
                for name, destinations in opt_in.items():
                    self.assertTrue(reasons[name].startswith(OPT_IN_PREFIX), reasons[name])
                    for destination in destinations:
                        self.assertIn(destination, reasons[name])


    # ---- Revision 11: every repository in the packet, queue id ------------------

    def two_repository_beta(self):
        """BETA Testing with a passing current run on a candidate of two repositories;
        returns (candidate id as printed by `candidate`, tree 1, tree 2)."""
        second = self.tmp / "repo2"
        second.mkdir()
        for args in (("init", "-q"), ("config", "user.name", "Test Writer"),
                     ("config", "user.email", "test@example.invalid"),
                     ("config", "commit.gpgsign", "false")):
            subprocess.run(["git", "-C", str(second), *args], env=self.env, check=True,
                           capture_output=True)
        (second / "other.txt").write_text("other\n")
        for args in (("add", "other.txt"), ("commit", "-q", "-m", "other")):
            subprocess.run(["git", "-C", str(second), *args], env=self.env, check=True,
                           capture_output=True)
        tree2 = subprocess.run(["git", "-C", str(second), "rev-parse", "HEAD^{tree}"],
                               env=self.env, check=True, text=True,
                               capture_output=True).stdout.strip()
        self.ok(self.tool("lock", "take"))
        self.ok(self.tool("criterion", BETA, "--file", self.criterion))
        printed = self.ok(self.tool("candidate", BETA, "--repo", self.repo, "--commit", "HEAD",
                                    "--repo", second, "--commit", "HEAD")).stdout
        entry = job_entry(self.status(), BETA)
        self.assertEqual(sorted(r["tree"] for r in entry["candidate"]),
                         sorted([self.tree(), tree2]))
        self.assertIn(entry["candidate_id"], printed)
        self.ok(self.set(BETA, "Preparing"))
        self.start_worker(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        return entry["candidate_id"], self.tree(), tree2

    def identity_refused(self, packet, reason):
        result = self.ready(BETA, "--packet", packet)
        self.refused(result, 8)
        self.assertIn(reason, result.stderr)
        self.assertEqual(self.state(BETA), "Testing")

    EVERY_TREE = ("packet CODE_IDENTITY.txt does not name the tested tree of every candidate "
                  "repository")
    CANDIDATE_ID = "packet CODE_IDENTITY.txt does not name the candidate id"

    def test_packet_must_name_candidate_id(self):
        """Revision 11: CODE_IDENTITY.txt needs `candidate_id: ID` equal to the recorded
        candidate id; missing, different or a prefix of it is refused (`packet
        CODE_IDENTITY.txt does not name the candidate id`); control: the right id."""
        identity, tree1, tree2 = self.two_repository_beta()
        both = [tree1, tree2]
        for label, value in (("missing", None), ("different", "0" * 64),
                             ("prefix", identity[:32]), ("upper case", identity.upper())):
            with self.subTest(candidate_id=label):
                self.identity_refused(self.packet(tree1, candidate_id=value, trees=both),
                                      self.CANDIDATE_ID)
        self.ok(self.ready(BETA, "--packet", self.packet(tree1, candidate_id=identity,
                                                         trees=both)))

    def test_packet_must_name_every_candidate_tree(self):
        """Revision 11 (finding 1): with a two-repository candidate, a packet naming
        only one tree, or neither, is refused at Ready (`packet CODE_IDENTITY.txt does
        not name the tested tree of every candidate repository`); both trees are
        accepted, and extra tested_tree lines are allowed."""
        identity, tree1, tree2 = self.two_repository_beta()
        for label, trees in (("only the first", [tree1]), ("only the second", [tree2]),
                             ("neither", []), ("neither, another tree", ["f" * 40])):
            with self.subTest(trees=label):
                self.identity_refused(self.packet(tree1, candidate_id=identity, trees=trees),
                                      self.EVERY_TREE)
        self.ok(self.ready(BETA, "--packet", self.packet(
            tree1, candidate_id=identity, trees=[tree2, "e" * 40, tree1, "d" * 40])))
        self.assertIn("local checks passed", self.result_cell(self.report(), BETA))
        self.ok(self.approve(BETA, f"approve {BETA}"))

    def test_packet_both_trees_accepted_exactly(self):
        """Revision 11 control: the packet naming both trees (and the id) is accepted
        at Ready and Approved."""
        identity, tree1, tree2 = self.two_repository_beta()
        self.ok(self.ready(BETA, "--packet", self.packet(tree1, candidate_id=identity,
                                                         trees=[tree1, tree2])))
        self.ok(self.approve(BETA, f"approve {BETA}"))
        self.assertEqual(self.state(BETA), "Approved")

    def test_every_candidate_tree_rechecked_at_approved(self):
        """Revision 11 (finding 1): the rule is re-checked at Approved. A job Ready with
        both trees named gets a third repository added to its recorded candidate by
        hand (same candidate id, so runs stay current): the row is no longer passed
        and Approved is refused (8) with the every-repository reason; with the record
        restored, Approved is accepted."""
        identity, tree1, tree2 = self.two_repository_beta()
        self.ok(self.ready(BETA, "--packet", self.packet(tree1, candidate_id=identity,
                                                         trees=[tree1, tree2])))
        path = self.queue / "QUEUE.json"
        original = path.read_text()
        data = json.loads(original)
        entry = job_entry(data, BETA)
        entry["candidate"].append({"repo": str(self.tmp / "repo3"), "commit": "c" * 40,
                                   "tree": "a" * 40})
        temp = path.with_name("QUEUE.json.hand")
        temp.write_text(json.dumps(data))
        os.replace(temp, path)
        cell = self.result_cell(self.report(), BETA)
        self.assertNotIn("passed", cell.lower())
        self.assertIn(self.EVERY_TREE, cell)
        quote = f"approve {BETA}"
        result = self.approve(BETA, quote)
        self.refused(result, 8)
        self.assertIn(self.EVERY_TREE, result.stderr)
        self.assertEqual(self.state(BETA), "Ready for approval")
        temp.write_text(original)
        os.replace(temp, path)
        self.ok(self.approve(BETA, quote))

    def snapshot(self):
        return sorted(str(p) for p in self.tmp.rglob("*"))

    def test_init_queue_id_must_be_one_directory_name(self):
        """Revision 11 (finding 2): init refuses (exit 2, `queue id must be one directory
        name`) an id that is empty, `.` or `..`, contains `/`, or does not match
        `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; nothing is created anywhere (no directory
        inside or outside the records root, no ACTIVE file). Ordinary names, up to 64
        characters, are accepted."""
        for index, qid in enumerate(("..", ".", "../escaped-queue", str(self.tmp / "abs-queue"),
                                     "a/b", "", "-a", ".hidden", "a" * 65, "a b", "a|b",
                                     "../../escaped", "a/..", "qé")):
            with self.subTest(queue_id=qid):
                home = self.tmp / f"home-qid-{index}"
                home.mkdir()
                before = self.snapshot()
                args = self.init_args(qid, self.tmp / f"RUN-qid-{index}.lock")
                # `--queue-id=VALUE` so a value starting with `-` reaches the tool's check.
                at = args.index("--queue-id")
                args[at:at + 2] = [f"--queue-id={qid}"]
                result = self.tool(*args, home=home, queue=None)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertIn("queue id must be one directory name", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertEqual(self.snapshot(), before)
                self.assertFalse((home / "ACTIVE").exists())
                self.assertFalse((self.tmp / "escaped-queue").exists())
                self.assertFalse((self.tmp.parent / "escaped").exists())
        for qid in ("A", "night-1", "Night_2.b", "0q", "a" * 64, "a.b-c_d"):
            with self.subTest(accepted=qid):
                self.use_new_queue(qid)
                self.assertEqual(self.queue, self.records / qid)
                self.assertTrue((self.records / qid / "QUEUE.json").is_file())


    # ---- After re-read 11: containment check, candidate_id line forms -------------

    def identity_packet(self, identity_text):
        """A frozen packet whose CODE_IDENTITY.txt is exactly `identity_text`."""
        packet = self.packet(self.tree(), manifest=False, candidate_id=None)
        (packet / "CODE_IDENTITY.txt").write_text(identity_text)
        self.refreeze(packet)
        return packet

    def test_candidate_id_line_forms(self):
        """Revision 11 / after re-read 11: two candidate_id lines, one right and one
        wrong, are refused; a repeated identical line and spaces after the colon are
        accepted; a line with leading spaces is not a candidate_id line (refused when
        it is the only one)."""
        identity, tree1, tree2 = self.two_repository_beta()
        trees = f"tested_tree: {tree1}\ntested_tree: {tree2}\n"
        for label, text in (
                ("right then wrong", f"candidate_id: {identity}\ncandidate_id: {'0' * 64}\n"),
                ("wrong then right", f"candidate_id: {'0' * 64}\ncandidate_id: {identity}\n"),
                ("only line indented", f"  candidate_id: {identity}\n")):
            with self.subTest(refused=label):
                self.identity_refused(self.identity_packet(text + trees), self.CANDIDATE_ID)
        for label, text in (
                ("repeated identical", f"candidate_id: {identity}\ncandidate_id: {identity}\n"),
                ("spaces after the colon", f"candidate_id:    {identity}\n")):
            with self.subTest(accepted=label):
                self.ok(self.ready(BETA, "--packet", self.identity_packet(text + trees)))
                self.ok(self.set(BETA, "Preparing"))
                self.ok(self.set(BETA, "Testing"))

    def init_into(self, records, qid, index):
        """init with `records` (as given) in a fresh home (made by contain_home);
        returns (result, home)."""
        home = self.contain_home(index)
        result = self.tool(*self.init_args(qid, self.tmp / f"RUN-contain-{index}.lock",
                                           records=records), home=home, queue=None)
        self.assertNotIn("Traceback", result.stderr)
        return result, home

    def contain_home(self, index):
        home = self.tmp / f"home-contain-{index}"
        home.mkdir(exist_ok=True)
        return home

    def test_init_queue_path_link_outside_records_root_refused(self):
        """After re-read 11 (b): an existing link records/ID to a real directory outside
        the records root is refused (exit 2, `queue id must be one directory name
        inside the records root`); a dangling link pointing outside is refused (exit
        2); nothing is created anywhere, the target stays empty, no ACTIVE file."""
        target = self.tmp / "outside-target"
        target.mkdir()
        dangling = self.tmp / "outside-missing"
        for index, (qid, destination, message) in enumerate((
                ("linked-out", target,
                 "queue id must be one directory name inside the records root"),
                ("dangling-out", dangling, None))):
            with self.subTest(queue_id=qid):
                (self.records / qid).symlink_to(destination)
                self.contain_home(index)
                before = self.snapshot()
                result, home = self.init_into(self.records, qid, index)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                if message:
                    self.assertIn(message, result.stderr)
                self.assertEqual(self.snapshot(), before)
                self.assertFalse((home / "ACTIVE").exists())
                self.assertEqual(list(target.iterdir()), [])
                self.assertFalse(dangling.exists())
                self.assertTrue((self.records / qid).is_symlink())

    def test_init_queue_path_link_inside_records_root_refused(self):
        """After re-read 11 (a): a queue path that is a link inside the records root,
        dangling or to a real directory, is refused (exit 2, `queue directory …
        exists`), no traceback, nothing created."""
        real = self.records / "real-inside"
        real.mkdir()
        for index, (qid, destination) in enumerate((("dangling-in", self.records / "missing-in"),
                                                    ("linked-in", real))):
            with self.subTest(queue_id=qid):
                (self.records / qid).symlink_to(destination)
                self.contain_home(10 + index)
                before = self.snapshot()
                result, home = self.init_into(self.records, qid, 10 + index)
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertRegex(result.stderr, r"queue directory .* exists")
                self.assertEqual(self.snapshot(), before)
                self.assertFalse((home / "ACTIVE").exists())
                self.assertEqual(list(real.iterdir()), [])
                self.assertFalse((self.records / "missing-in").exists())

    def test_init_records_root_relative_or_through_link_is_resolved(self):
        """After re-read 11: a records root given relatively (to the tool's cwd) or
        through a symbolic link is resolved first; the queue is created in the real
        root, ACTIVE names it, and nothing appears under the link's own name."""
        relative = self.tmp / "records-rel"
        relative.mkdir()
        result, home = self.init_into("records-rel", "rel-1", 20)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.strip().splitlines()[-1], str(relative / "rel-1"))
        self.assertTrue((relative / "rel-1" / "QUEUE.json").is_file())
        self.assertEqual((home / "ACTIVE").read_text().splitlines()[0], str(relative / "rel-1"))
        real = self.tmp / "records-real"
        real.mkdir()
        link = self.tmp / "records-link"
        link.symlink_to(real)
        result, home = self.init_into(link, "linked-1", 21)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.strip().splitlines()[-1], str(real / "linked-1"))
        self.assertTrue((real / "linked-1" / "QUEUE.json").is_file())
        self.assertEqual((home / "ACTIVE").read_text().splitlines()[0], str(real / "linked-1"))
        self.assertTrue(link.is_symlink())
        self.assertEqual(sorted(p.name for p in self.tmp.glob("records-link*")), ["records-link"])


    # ---- Revision 12: macOS only, lock directory required -----------------------

    def forced_non_macos(self, *args, queue=True, env=None):
        """The real helper (auto_queue.py imported by a plain subprocess, never the
        test-only launcher) with sys.platform set to `linux` in that process only."""
        command = self.command(args, None, queue, real=True)
        command[1:4] = ["-I", "-B", "-c", FORCED_NON_MACOS, str(TOOL.parent)]
        return subprocess.run(command, cwd=self.tmp, env=dict(self.env, **(env or {})),
                              text=True, capture_output=True, timeout=60)

    def check_platform_refusal(self, result, cmd_marker, probe_marker):
        """Exit 15 with the one message; neither the probe nor CMD started; no run
        record or run event written."""
        self.refused(result, PLATFORM_REFUSED)
        self.assertEqual(result.stderr.strip().splitlines(), [PLATFORM_MESSAGE + "linux"])
        self.assertFalse(cmd_marker.exists(), "CMD started")
        self.assertFalse(probe_marker.exists(), "the probe started")
        self.assertEqual(self.runs(), {})
        events = (self.queue / "EVENTS.log").read_text()
        self.assertNotRegex(events, r"(?m)^\S+ run ")

    def refusal_attempt(self, job, name, env=None):
        cmd_marker, cmd = self.marker_cmd(f"{name}-cmd")
        probe_marker = self.work / f"{name}-probe"
        probe = f"touch {shlex.quote(str(probe_marker))}; echo {shlex.quote(str(self.src))}"
        result = self.forced_non_macos("run", job, "unit", "--log", self.work / f"{name}.log",
                                       "--probe", probe, "--", *cmd, env=env)
        self.check_platform_refusal(result, cmd_marker, probe_marker)
        return cmd, probe, cmd_marker, probe_marker

    def test_run_refused_off_macos_with_and_without_network_grant(self):
        """Revision 12: off macOS (forced inside the real helper's own process) `run`
        exits 15 with `auto-mode runs need macOS (sandbox-exec and `ps -E`)`, for a
        job without and a job with "network": true, before the probe or CMD starts and
        before any run record; the job stays Preparing and `status` still works.
        Control (real macOS only, where the real helper may run): the same call is
        accepted and runs the probe and CMD."""
        jobs = self.tmp / "JOBS-net.json"
        jobs.write_text(json.dumps([dict(JOBS[1], network=True)]))
        for network in (False, True):
            with self.subTest(network=network):
                if network:
                    home = self.tmp / "home-q-net"
                    home.mkdir()
                    self.home, self.lock = home, self.tmp / "RUN-q-net.lock"
                    self.ok(self.tool(*self.init_args("q-net", self.lock, jobs=jobs), home=home,
                                      queue=None, real=True))
                    self.queue = self.records / "q-net"
                self.assertEqual(job_entry(self.status(), BETA).get("network", False), network)
                self.prepare(BETA, real=True)
                cmd, probe, cmd_marker, probe_marker = self.refusal_attempt(
                    BETA, f"refused-{network}")
                self.assertEqual(self.state(BETA), "Preparing")
                shown = self.forced_non_macos("status")
                self.assertEqual(shown.returncode, 0, shown.stderr)
                self.assertEqual(job_entry(json.loads(shown.stdout), BETA)["state"], "Preparing")
                if sys.platform == "darwin" and not NON_MACOS:
                    self.ok(self.tool("run", BETA, "unit", "--log", self.work / "control.log",
                                      "--probe", probe, "--", *cmd, real=True))
                    self.assertTrue(cmd_marker.exists())
                    self.assertTrue(probe_marker.exists())
                    (record,) = self.runs().values()
                    self.assertEqual(record["network"], "allowed" if network else "blocked")

    def test_run_refusal_has_no_switch(self):
        """Revision 12: no environment variable or option lets `run` go ahead off
        macOS: with plausible variables set (including the tests' own simulation
        variable) and with --timeout and --cwd, the real helper still exits 15 with
        the same message and starts nothing."""
        self.prepare(BETA, real=True)
        plausible = {"GC_TEST_SIMULATE_NON_MACOS": "1", "GC_AUTO_ALLOW_NON_MACOS": "1",
                     "GC_AUTO_PLATFORM": "darwin", "GC_AUTO_NO_SANDBOX": "1",
                     "GC_AUTO_SANDBOX": "/usr/bin/true", "GC_ALLOW_UNSANDBOXED": "1",
                     "SANDBOX": "/bin/sh", "PLATFORM": "darwin", "OSTYPE": "darwin"}
        self.refusal_attempt(BETA, "switch-env", env=plausible)
        cmd_marker, cmd = self.marker_cmd("switch-options-cmd")
        result = self.forced_non_macos("run", BETA, "unit", "--log", self.work / "o.log",
                                       "--timeout", "30", "--cwd", self.work, "--", *cmd,
                                       env=plausible)
        self.check_platform_refusal(result, cmd_marker, self.work / "no-probe")

    def test_init_without_lock_dir_refused(self):
        """Revision 12: init without --lock-dir (or with an empty one) exits 2, saying
        the lock directory comes from the project method's lock directory; it creates
        nothing (no queue folder, no ACTIVE, ACTIVE unchanged where it exists). Control:
        the same call with --lock-dir creates the queue."""
        for index, home in enumerate((self.contain_home(40), self.home)):
            for lock in (None, ""):
                with self.subTest(home=home.name, lock=lock):
                    args = self.init_args("q-nolock", self.tmp / "unused.lock")
                    at = args.index("--lock-dir")
                    args[at:at + 2] = [] if lock is None else ["--lock-dir", lock]
                    active = home / "ACTIVE"
                    before_active = active.read_text() if active.exists() else None
                    before = self.snapshot()
                    result = self.tool(*args, home=home, queue=None, real=True)
                    self.refused(result, 2)
                    self.assertIn("--lock-dir", result.stderr)
                    self.assertIn("project method's lock directory", result.stderr)
                    self.assertEqual(self.snapshot(), before)
                    self.assertFalse((self.records / "q-nolock").exists())
                    self.assertEqual(active.read_text() if active.exists() else None,
                                     before_active)
        home = self.contain_home(41)
        self.ok(self.tool(*self.init_args("q-nolock", self.tmp / "RUN-nolock.lock"), home=home,
                          queue=None, real=True))
        self.assertTrue((self.records / "q-nolock" / "QUEUE.json").is_file())
        self.assertEqual(json.loads((self.records / "q-nolock" / "QUEUE.json").read_text())
                         ["lock_dir"], str(self.tmp / "RUN-nolock.lock"))

    # ---- Revision 12 (P8): the test-only launcher cannot run a real queue job ----

    def standin(self, *args, env=None, drop=(), cwd=None):
        """The test-only launcher itself, on any platform."""
        environment = {k: v for k, v in dict(self.env, **(env or {})).items() if k not in drop}
        return subprocess.run([sys.executable, "-I", "-B", str(STANDIN), *map(str, args)],
                              cwd=cwd or self.tmp, env=environment, text=True,
                              capture_output=True, timeout=60)

    def files_and_bytes(self, *directories):
        return {str(p): (p.read_bytes() if p.is_file() else None)
                for d in directories for p in sorted(Path(d).rglob("*"))}

    def other_temp_queue(self, marked):
        """A queue made with the real helper (init, lock, criterion, candidate, Preparing
        and a worker: no `run`, so on any platform) in a new temporary directory of its
        own, with or without the test-only marker."""
        other = tempfile.TemporaryDirectory()
        self.addCleanup(other.cleanup)
        root = Path(os.path.realpath(other.name))
        if marked:
            (root / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)
        home, records = root / "home", root / "records"
        home.mkdir()
        records.mkdir()
        self.ok(self.tool(*self.init_args("q-other", root / "RUN.lock", records=records),
                          home=home, queue=None, real=True))
        queue = records / "q-other"
        self.prepare(BETA, real=True, home=home, queue=queue)
        return root, home, queue

    def standin_refused(self, name, *where, reason, env=None, drop=(), watch=(), cwd=None):
        """`run` through the launcher is refused with exit 16 and one message naming
        `reason`; neither the probe nor CMD starts; nothing in `watch` changes."""
        cmd_marker, cmd = self.marker_cmd(f"{name}-cmd")
        probe_marker = self.work / f"{name}-probe"
        probe = f"touch {shlex.quote(str(probe_marker))}; echo {shlex.quote(str(self.src))}"
        before = self.files_and_bytes(*watch)
        result = self.standin(*where, "run", BETA, "unit", "--log", self.work / f"{name}.log",
                              "--probe", probe, "--", *cmd, env=env, drop=drop, cwd=cwd)
        self.refused(result, STANDIN_REFUSED)
        lines = result.stderr.strip().splitlines()
        self.assertEqual(len(lines), 1, result.stderr)
        self.assertTrue(lines[0].startswith(STANDIN_MESSAGE), lines[0])
        self.assertIn(reason, lines[0])
        self.assertFalse(cmd_marker.exists(), "CMD started")
        self.assertFalse(probe_marker.exists(), "the probe started")
        self.assertFalse((self.work / f"{name}.log").exists(), "run log written")
        self.assertEqual(self.files_and_bytes(*watch), before)
        return cmd, probe, cmd_marker

    def test_standin_refuses_queue_without_marker(self):
        """Revision 12 (P8): the launcher refuses (exit 16, before the probe, CMD or any
        record) a queue whose temporary directory has no test-only marker: by --queue,
        through the home's ACTIVE file, from the marked home of this test, with the
        default home, with the marker only in a subdirectory, as a link, or with other
        text; `status` too. Control: with the marker written, the same launcher runs the
        same queue's job (a local command)."""
        root, home, queue = self.other_temp_queue(marked=False)
        watch = (root,)
        no_marker = f"has no test-only marker {STANDIN_MARKER}"
        self.standin_refused("p8-queue", "--home", home, "--queue", queue,
                             reason=no_marker, watch=watch)
        self.standin_refused("p8-active", "--home", home, reason=no_marker, watch=watch)
        self.standin_refused("p8-env-home", "--queue", queue, env={"GC_AUTO_HOME": str(home)},
                             reason=no_marker, watch=watch)
        self.standin_refused("p8-marked-home", "--home", self.home, "--queue", queue,
                             reason=no_marker, watch=watch)
        self.standin_refused("p8-default-home", "--queue", queue, drop=("GC_AUTO_HOME",),
                             reason="no helper home given", watch=watch)
        result = self.standin("--home", home, "--queue", queue, "status")
        self.refused(result, STANDIN_REFUSED)
        self.assertEqual(result.stdout, "")
        (root / "records" / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)
        self.standin_refused("p8-deeper", "--home", home, "--queue", queue,
                             reason=no_marker, watch=watch)
        (root / STANDIN_MARKER).symlink_to(root / "records" / STANDIN_MARKER)
        self.standin_refused("p8-link", "--home", home, "--queue", queue,
                             reason=no_marker, watch=watch)
        (root / STANDIN_MARKER).unlink()
        (root / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT + "x")
        self.standin_refused("p8-text", "--home", home, "--queue", queue,
                             reason=no_marker, watch=watch)
        (root / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)
        shown = self.standin("--home", home, "--queue", queue, "status")
        self.assertEqual(shown.returncode, 0, shown.stderr)
        self.assertEqual(job_entry(json.loads(shown.stdout), BETA)["state"], "Preparing")
        cmd_marker, cmd = self.marker_cmd("p8-control-cmd")
        self.ok(self.standin("--home", home, "--queue", queue, "run", BETA, "unit", "--log",
                             self.work / "p8-control.log", "--probe",
                             "echo " + shlex.quote(str(self.src)), "--", *cmd))
        self.assertTrue(cmd_marker.exists())

    def test_standin_refuses_marked_queue_outside_temporary_directory(self):
        """Revision 12 (P8): a queue that has the marker but is not inside the launcher's
        temporary directory (TMPDIR moved to a folder inside this test's directory, so
        nothing is written elsewhere), or whose home is in another marked directory, is
        refused with exit 16 before the probe, CMD or any record; init through the
        launcher with its records root or lock directory outside a marked directory is
        refused and creates nothing. Control: the same queue with the real temporary
        directory is accepted (`status`)."""
        self.prepare(BETA, real=True)
        moved = self.tmp / "other-temp"
        moved.mkdir()
        watch = (self.queue, self.home)
        self.standin_refused("p8-outside", "--home", self.home, "--queue", self.queue,
                             env={"TMPDIR": str(moved)},
                             reason="is not inside the temporary directory", watch=watch)
        root, home, queue = self.other_temp_queue(marked=True)
        self.standin_refused("p8-split", "--home", home, "--queue", self.queue,
                             reason="different test directories", watch=watch + (root,))
        unmarked = tempfile.TemporaryDirectory()
        self.addCleanup(unmarked.cleanup)
        outside = Path(os.path.realpath(unmarked.name))
        for name, args in (("records", self.init_args("q-p8", self.tmp / "RUN-p8.lock",
                                                       records=outside)),
                           ("lock", self.init_args("q-p8", outside / "RUN.lock"))):
            with self.subTest(init=name):
                before = self.files_and_bytes(self.tmp, outside)
                result = self.standin("--home", self.home, *args)
                self.refused(result, STANDIN_REFUSED)
                self.assertIn(f"has no test-only marker {STANDIN_MARKER}", result.stderr)
                self.assertEqual(self.files_and_bytes(self.tmp, outside), before)
        shown = self.standin("--home", self.home, "--queue", self.queue, "status")
        self.assertEqual(shown.returncode, 0, shown.stderr)

    def test_standin_refuses_empty_home_that_helper_ignores(self):
        """Revision 12 (P8, REVIEW_4 F1): an empty or blank --home (`--home=`, `--home ""`,
        `--home " "`), or an empty GC_AUTO_HOME, is refused with exit 16 even when the
        current directory is a marked test directory. auto_queue.py treats an empty home
        as absent and falls back to GC_AUTO_HOME or the default home, whose ACTIVE names
        a "real" queue outside the launcher's temporary directory (TMPDIR is moved to a
        folder inside this test's directory, so nothing is written elsewhere). Includes
        the reviewer's reproduction (`--home "" run`, `--home= status`). Nothing in the
        real queue changes. Control: the real helper with `--home=` does open that queue."""
        real = self.tmp / "real"
        home, records, user = real / "home", real / "records", real / "user"
        for directory in (home, records, user):
            directory.mkdir(parents=True)
        self.ok(self.tool(*self.init_args("q-real", real / "RUN.lock", records=records),
                          home=home, queue=None, real=True))
        queue = records / "q-real"
        self.prepare(BETA, real=True, home=home, queue=queue)
        default_home = user / "GuidedCoding" / "auto"
        default_home.mkdir(parents=True)
        (default_home / "ACTIVE").write_text((home / "ACTIVE").read_text())
        temp = self.tmp / "launcher-temp"
        marked_cwd = temp / "test-run"
        marked_cwd.mkdir(parents=True)
        (marked_cwd / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)
        watch = (real,)
        via_env = {"TMPDIR": str(temp), "GC_AUTO_HOME": str(home)}
        via_default = {"TMPDIR": str(temp), "HOME": str(user)}
        empty = "empty or blank --home"
        for name, where, env, drop, reason in (
                ("p8-home-equals", ("--home=",), via_env, (), empty),
                ("p8-home-empty", ("--home", ""), via_env, (), empty),
                ("p8-home-blank", ("--home", " "), via_env, (), empty),
                ("p8-home-empty-default", ("--home", ""), via_default, ("GC_AUTO_HOME",), empty),
                ("p8-home-equals-queue", ("--home=", "--queue", queue), via_env, (), empty),
                ("p8-env-home-empty", (), dict(via_default, GC_AUTO_HOME=""), (),
                 "no helper home given"),
                ("p8-env-home-blank", (), dict(via_default, GC_AUTO_HOME=" "), (),
                 "no helper home given")):
            with self.subTest(case=name):
                self.standin_refused(name, *where, reason=reason, env=env, drop=drop,
                                     watch=watch, cwd=marked_cwd)
        before = self.files_and_bytes(real)
        result = self.standin("--home=", "status", env=via_default, drop=("GC_AUTO_HOME",),
                              cwd=marked_cwd)
        self.refused(result, STANDIN_REFUSED)
        self.assertEqual(result.stdout, "")
        self.assertEqual(self.files_and_bytes(real), before)
        shown = subprocess.run([sys.executable, "-I", "-B", str(TOOL), "--home=", "status"],
                               cwd=marked_cwd, env=dict(self.env, **via_default, GC_AUTO_HOME=""),
                               text=True, capture_output=True, timeout=60)
        self.assertEqual(shown.returncode, 0, shown.stderr)
        self.assertEqual(json.loads(shown.stdout)["queue_id"], "q-real")

    def test_standin_refuses_queue_named_by_marked_home_active(self):
        """Revision 12 (P8, REVIEW_4 F2): a MARKED home (this test's) whose ACTIVE file
        names a queue in an unmarked directory is refused with exit 16 when no --queue is
        given, by --home and by GC_AUTO_HOME; only the ACTIVE check can refuse here, since
        the home itself is marked. Control: with ACTIVE naming this test's own queue the
        same call runs the job (a local command)."""
        self.prepare(BETA, real=True)
        root, home, queue = self.other_temp_queue(marked=False)
        own_active = (self.home / "ACTIVE").read_text()
        (self.home / "ACTIVE").write_text(str(queue) + "\n")
        watch = (root, self.queue)
        reason = f"queue named by ACTIVE {queue}: "
        self.standin_refused("p8-active-arg", "--home", self.home, reason=reason, watch=watch)
        self.standin_refused("p8-active-env", reason=reason, watch=watch,
                             env={"GC_AUTO_HOME": str(self.home)})
        result = self.standin("--home", self.home, "status")
        self.refused(result, STANDIN_REFUSED)
        self.assertIn(reason, result.stderr)
        self.assertEqual(result.stdout, "")
        (self.home / "ACTIVE").write_text(own_active)
        cmd_marker, cmd = self.marker_cmd("p8-active-control-cmd")
        self.ok(self.standin("--home", self.home, "run", BETA, "unit", "--log",
                             self.work / "p8-active-control.log", "--probe",
                             "echo " + shlex.quote(str(self.src)), "--", *cmd))
        self.assertTrue(cmd_marker.exists())


    # ---- Revision 13 (item 14): an unreadable process table fails closed ---------

    def fake_ps_path(self, mode):
        if mode == "missing":
            return self.work / "no-such-ps"
        path = self.work / f"ps-{mode}"
        path.write_text(FAKE_PS[mode])
        path.chmod(0o755)
        return path

    def snapshot_probe(self, ps, platform="native"):
        script = self.work / "snapshot_probe.py"
        script.write_text(SNAPSHOT_PROBE)
        result = subprocess.run([sys.executable, "-I", "-B", str(script), str(TOOL.parent), str(ps),
                                 platform], cwd=self.tmp, env=self.env, text=True,
                                capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout.strip()

    def test_snapshot_fails_closed_when_ps_fails(self):
        """Revision 13 (item 14): snapshot() raises a refusal (17, "cannot read the process
        table") when ps cannot start, exits nonzero, prints nothing or prints a line that
        does not parse; it never returns an empty table. Control: the real ps gives a table.
        Linux: a failing `ps -E` marker query alone is tolerated off macOS (markers dropped,
        as before) and refused on macOS, where the marker is needed."""
        table = self.snapshot_probe("/bin/ps")
        self.assertTrue(table.startswith("table "), table)
        self.assertGreater(int(table.split()[1]), 1)
        for mode in ("missing", "refused", "empty", "malformed"):
            with self.subTest(ps=mode):
                out = self.snapshot_probe(self.fake_ps_path(mode))
                self.assertTrue(out.startswith("refused 17 cannot read the process table"), out)
                self.assertIn("whether queue processes are alive is unknown", out)
        no_marker = self.fake_ps_path("no-marker")
        out = self.snapshot_probe(no_marker, "linux")
        self.assertTrue(out.startswith("table "), out)
        self.assertEqual(out.split()[2], "0", "markers should be dropped, not invented")
        if sys.platform == "darwin":
            out = self.snapshot_probe(no_marker)
            self.assertTrue(out.startswith("refused 17 "), out)
            self.assertIn("exited 1", out)

    def test_snapshot_refuses_empty_or_unparsable_marker_output_on_macos(self):
        """Checkpoint 3 (REVIEW_CP2 finding 3): on macOS a `ps -E` marker query that exits
        0 with no output, or with a line that does not parse, refuses (17) like the main
        table, since a survivor found only by its marker would be invisible. Off macOS
        (simulated in the probe) the same output only drops the markers, as before."""
        for mode in ("empty-marker", "garbage-marker"):
            with self.subTest(ps=mode):
                out = self.snapshot_probe(self.fake_ps_path(mode), "linux")
                self.assertTrue(out.startswith("table "), out)
                if sys.platform != "darwin":
                    continue
                out = self.snapshot_probe(self.fake_ps_path(mode))
                self.assertTrue(out.startswith("refused 17 cannot read the process table"), out)
                self.assertIn("ps -E listed no environment line for this process", out)

    def test_stop_signals_confirmed_run_when_another_record_needs_the_table(self):
        """Checkpoint 3 (REVIEW_CP2 finding 1): a completed run's record sorts before a
        running one; with a ps that fails only for the whole table (`ps -p` works), stop
        still sends SIGTERM to the running run's group (confirmed by `ps -p`), then
        refuses 17 with Stopping, never Stopped. The run ends 128+15."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit", cmd=QUICK_CMD))
        (first,) = self.runs()
        time.sleep(1.02 - time.time() % 1)  # the next record's name sorts after the first
        proc = self.popen_run(BETA, sleep_cmd(30))
        stem, record = self.wait_record({first})
        self.assertLess(first, stem)
        self.end_worker()
        self.fake_ps = self.fake_ps_path("ax-fails")
        try:
            result = self.tool("stop")
        finally:
            self.fake_ps = None
        self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
        self.assertTrue(result.stdout.strip().startswith("Stopping:"), result.stdout)
        self.assertEqual(proc.wait(TIMEOUT), 128 + signal.SIGTERM)
        self.assertIn(f"stop sent SIGTERM to group {record['pgid']} of run {stem}",
                      (self.queue / "EVENTS.log").read_text())
        self.stop_until_stopped()

    def test_unreadable_table_refuses_stop_cleanup_release_transfer_without_runs(self):
        """Checkpoint 3 (REVIEW_CP2 finding 2): with no run record at all and the helper's
        ps refused, stop prints Stopping (17, never Stopped) and cleanup, lock release and
        transfer refuse 17 (lock kept, job not transferred), consistent with report.
        Control: with the real ps each succeeds."""
        self.prepare(BETA)
        self.end_worker()
        successor = self.tmp / "attended-none"
        successor.mkdir()
        self.fake_ps = self.fake_ps_path("refused")
        try:
            result = self.tool("cleanup", BETA)
            self.refused(result, 17)
            self.assertNotIn("Cleaned", result.stdout)
            self.refused(self.tool("lock", "release"), 17)
            self.assertTrue((self.lock / "owner.txt").exists())
            self.refused(self.tool("transfer", BETA, "--to", successor, "--note", "attended"), 17)
            self.assertNotEqual(self.state(BETA), "Transferred")
            result = self.tool("stop")
            self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
            self.assertTrue(result.stdout.strip().startswith("Stopping:"), result.stdout)
            self.assertNotIn("Stopped", result.stdout)
            report = self.tool("report")
            self.assertEqual(report.returncode, 17, report.stdout + report.stderr)
            self.assertEqual(self.first_line(report.stdout), "**Stopping**")
        finally:
            self.fake_ps = None
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))  # control
        self.stop_until_stopped()
        self.ok(self.tool("transfer", BETA, "--to", successor, "--note", "attended"))
        self.ok(self.tool("lock", "release"))

    def test_unreadable_process_table_never_reports_stopped_cleaned_or_released(self):
        """Revision 13 (item 14): with a run's leftover process alive and the helper's ps
        missing, refused (as inside sandbox-exec), empty or malformed: `run` refuses (17) and
        its CMD never starts; Ready, cleanup and lock release refuse (17) and change nothing
        (the leftover lives, the lock stays); reconcile marks no run unknown; report prints no
        readiness and shows the job unknown, exits 17; stop records the stop but prints
        Stopping, never Stopped (17), and report then says **Stopping**. Control: with the
        real ps cleanup prints Cleaned, the lock is released and stop prints Stopped."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        stem, (child,), _ = self.launch("default")
        packet = self.packet(self.tree())
        modes = ("missing", "refused", "empty", "malformed")
        for mode in modes:
            with self.subTest(ps=mode):
                self.fake_ps = self.fake_ps_path(mode)
                try:
                    marker, cmd = self.marker_cmd(f"fake-ps-{mode}")
                    result = self.run_cmd(BETA, "unit", cmd=cmd)
                    self.refused(result, 17)
                    self.assertIn("cannot read the process table", result.stderr)
                    self.assertFalse(marker.exists(), "CMD started without a process table")
                    self.end_worker()
                    self.refused(self.set(BETA, "Ready for approval", "--packet", packet), 17)
                    self.start_worker(BETA)
                    self.assertIn(self.tool("reconcile").returncode, (0, 17))
                    result = self.tool("cleanup", BETA)
                    self.refused(result, 17)
                    self.assertNotIn("Cleaned", result.stdout)
                    self.assertTrue(alive(child))
                    self.refused(self.tool("lock", "release"), 17)
                    self.assertTrue((self.lock / "owner.txt").exists())
                    result = self.tool("report")
                    self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
                    self.assertIn(REPORT_HEADER, result.stdout.splitlines())
                    self.assertIn("**Process table unreadable:**", result.stdout)
                    self.assertNotIn("**Stopped**", result.stdout)
                    self.assertIn("unknown", self.result_cell(result.stdout, BETA))
                    self.assertNotIn("passed", self.result_cell(result.stdout, BETA))
                finally:
                    self.fake_ps = None
        self.assertEqual(self.state(BETA), "Testing")
        self.assertNotEqual(self.runs()[stem]["status"], "unknown")
        for mode in modes:
            with self.subTest(stop=mode):
                self.fake_ps = self.fake_ps_path(mode)
                try:
                    result = self.tool("stop")
                    self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
                    self.assertTrue(result.stdout.strip().startswith("Stopping:"), result.stdout)
                    self.assertNotIn("Stopped", result.stdout)
                    self.assertTrue((self.queue / "STOP_REQUESTED").exists())
                    report = self.tool("report")
                    self.assertEqual(report.returncode, 17, report.stdout + report.stderr)
                    self.assertEqual(self.first_line(report.stdout), "**Stopping**")
                finally:
                    self.fake_ps = None
        self.assertTrue(alive(child))
        self.end_worker()
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))  # control
        self.assertFalse(alive(child))
        self.ok(self.tool("lock", "release"))
        self.stop_until_stopped()

    # ---- Revision 14: a zombie counts as gone ----------------------------------------

    def zombie(self, cmd=None):
        """A real zombie: a child of this test process, in its own session, that exits and
        is not waited for (kill_leftovers collects it at the end). With `cmd` it is that
        command (which may leave children in its group); otherwise it exits at once."""
        proc = subprocess.Popen(cmd or [sys.executable, "-c", "pass"], start_new_session=True,
                                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
        self.procs.append(proc)
        self.wait_for(lambda: ps_state(proc.pid).startswith("Z"), "a zombie (ps stat Z)")
        self.assertFalse(alive(proc.pid))
        return proc.pid

    def zombie_leader_with_member(self, mode="note"):
        """A zombie session and group leader Z whose child M (mode `note`: notes SIGTERM and
        stays until released) lives on in Z's group and session. Returns (Z, M, ready,
        release)."""
        launcher, pids = self.work / "launch.py", self.work / f"pids-zombie-{self.counter}"
        launcher.write_text(LAUNCH_SCRIPT)
        _, ready, release = self.wait_cmd(mode)
        leader = self.zombie([sys.executable, str(launcher), "exit", str(self.wait_script),
                              str(pids), f"{mode},{ready},{release}"])
        (member,) = [int(pid) for pid in pids.read_text().split()]
        self.stray.append(member)
        self.wait_for(ready.exists, "the zombie leader's live member")
        self.assertTrue(alive(member))
        self.assertEqual((os.getpgid(member), os.getsid(member)), (leader, leader))
        return leader, member, ready, release

    def zombie_probe(self, *pids):
        script = self.work / "zombie_probe.py"
        script.write_text(ZOMBIE_PROBE)
        result = subprocess.run([sys.executable, "-I", "-B", str(script), str(TOOL.parent),
                                 *[str(pid) for pid in pids]], cwd=self.tmp, env=self.env,
                                text=True, capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return {int(pid): views for pid, views in json.loads(result.stdout).items()}

    def lstart(self, pid):
        return subprocess.run(["/bin/ps", "-o", "lstart=", "-p", str(pid)],
                              env={"LC_ALL": "C", "TZ": "UTC", "PATH": "/bin:/usr/bin"},
                              text=True, capture_output=True).stdout.strip()

    def survivors_probe(self, record):
        path = self.work / f"probe-record-{self.counter}.json"
        self.counter += 1
        path.write_text(json.dumps({k: v for k, v in record.items() if v is not None}))
        script = self.work / "survivors_probe.py"
        script.write_text(SURVIVORS_PROBE)
        result = subprocess.run([sys.executable, "-I", "-B", str(script), str(TOOL.parent),
                                 str(path)], cwd=self.tmp, env=self.env, text=True,
                                capture_output=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads(result.stdout)

    def test_zombie_holding_recorded_leader_id_with_other_start_is_not_ours(self):
        """Revision 14 (review finding 2): a record names as its leader an id now held by a
        zombie Z with another start time (the id was reused), and Z's group has a live member:
        the leader is not the run's, the member is no survivor, and stop prints Stopped
        without signalling it. Control: with Z's real start time Z is the leader, the member
        a survivor, stop Stopping and the member gets SIGTERM. Revision 13 gives [] (the
        zombie's start time is compared); the first revision 14 candidate (ae6e63ff), which
        left zombies out of the table, gave [member]."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        leader, member, ready, release = self.zombie_leader_with_member()
        wrong = time.strftime("%a %b %d %H:%M:%S %Y",
                              time.gmtime(self.utc(record["started"]) - 3600))
        right = self.lstart(leader)
        foreign = dict(record, status="running", pid=leader, pgid=leader, process_start=wrong,
                       ended=None, exit_code=None)
        view = self.survivors_probe(foreign)  # (the run's own dead probe group stays a leader)
        self.assertNotIn(leader, view["leaders"], "reused id counted as ours")
        self.assertEqual(view["survivors"], [], "an outsider counted as a survivor")
        view = self.survivors_probe(dict(foreign, process_start=right))
        self.assertIn(leader, view["leaders"])
        self.assertEqual(view["survivors"], [member])
        name = record["run_id"] + "-reused-zombie"
        self.forge(record, name, status="running", pid=leader, pgid=leader,
                   process_start=wrong, ended=None, exit_code=None)
        self.end_worker()
        result = self.ok(self.tool("stop"))
        self.assertTrue(result.stdout.strip().startswith("Stopped"), result.stdout)
        self.assertFalse(Path(str(ready) + ".term").exists(), "stop signalled an outsider")
        self.assertTrue(alive(member))
        self.forge(record, name, status="running", pid=leader, pgid=leader,
                   process_start=right, ended=None, exit_code=None)  # control
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(ready) + ".term").exists, "the member to get SIGTERM", 5)
        release.write_text("")

    @unittest.skipUnless(sys.platform.startswith("linux")
                         and os.uname().machine in ("x86_64", "aarch64"),
                         "Linux only (x86_64 or aarch64): a process whose main thread has "
                         "exited shows state Z only on Linux (macOS shows it running, S)")
    def test_multithreaded_process_with_exited_main_thread_is_alive(self):
        """Revision 14 (review finding 1): a live process whose main thread has exited while
        another thread runs (its leader shows /proc state Z with 2 threads, ps `Zl`, and its
        heartbeat still advances) is alive in the helper's views (process_exists, group_alive,
        listed live in snapshot()) and in the test's own check; a `running` record whose
        leader it is stays running after reconcile, and stop prints Stopping. Fails on the
        first revision 14 candidate (ae6e63ff), which read any Z as gone."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        script, beat, release = self.work / "mt.py", self.work / "mt-beat", self.work / "mt-release"
        script.write_text(MT_SCRIPT)
        proc = subprocess.Popen([sys.executable, str(script), str(beat), str(release)],
                                start_new_session=True)
        self.procs.append(proc)
        self.wait_for(lambda: (proc_stat(proc.pid) or ("", 0))[0] == "Z",
                      "the main thread to exit (leader state Z)")
        state, threads = proc_stat(proc.pid)
        self.assertGreater(threads, 1, "the worker thread runs")
        self.assertIn("l", ps_state(proc.pid), "ps shows the multi-threaded flag")
        first = beat.read_text()
        self.wait_for(lambda: beat.read_text() != first, "a heartbeat after the main thread exited")
        self.assertTrue(alive(proc.pid))
        self.assertEqual(self.zombie_probe(proc.pid)[proc.pid], [True, True, True])
        name = record["run_id"] + "-threads"
        self.forge(record, name, status="running", pid=proc.pid, pgid=proc.pid,
                   process_start=self.lstart(proc.pid), ended=None, exit_code=None)
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[name]["status"], "running")
        self.end_worker()
        self.stopping(self.tool("stop"))
        release.write_text("")
        proc.wait(TIMEOUT)

    def test_zombie_is_gone_live_sibling_is_not_in_helper_views(self):
        """Revision 14: in the helper's own process (process_exists, group_alive with the
        pid as group id, and the snapshot() table) a zombie that is alone in its group is
        gone; a live process beside it, and a zombie leader's group while a live member
        remains in it, are alive and listed; once that member has exited, the group with
        only its zombie leader is gone. Fails on revision 13, which reads a zombie as
        alive."""
        lone = self.zombie()
        cmd, ready, release = self.wait_cmd("default")
        sibling = subprocess.Popen(cmd, start_new_session=True)
        self.procs.append(sibling)
        self.wait_for(ready.exists, "the live sibling")
        leader, member, _, member_release = self.zombie_leader_with_member()
        views = self.zombie_probe(lone, sibling.pid, leader, member)
        self.assertEqual(views[lone], [False, False, False], "zombie alone in its group")
        self.assertEqual(views[sibling.pid], [True, True, True], "live sibling")
        self.assertEqual(views[leader][0], False, "zombie leader as a process")
        self.assertEqual(views[leader][2], False, "zombie leader in the table")
        self.assertEqual(views[leader][1], True, "group with a live member")
        self.assertEqual(views[member][0], True, "live member as a process")
        self.assertEqual(views[member][2], True, "live member in the table")
        member_release.write_text("")
        self.wait_for(lambda: not alive(member), "the member to exit")
        self.assertTrue(ps_state(leader).startswith("Z"), "the leader is still a zombie")
        views = self.zombie_probe(leader, sibling.pid)
        self.assertEqual(views[leader], [False, False, False], "group of only a zombie")
        self.assertEqual(views[sibling.pid], [True, True, True], "live sibling")
        release.write_text("")
        self.assertEqual(sibling.wait(TIMEOUT), 0)

    def test_zombie_run_leader_is_gone_and_unreadable_table_still_refuses(self):
        """Revision 14: a `running` record whose CMD leader (pid and group, with its real start
        time) is a zombie: the report does not show it running, a new run of the job is not
        refused (13), reconcile marks it unknown; with the helper's ps refused, cleanup, lock
        release and stop still refuse (17: never Cleaned, released or Stopped); with the real
        ps stop prints Stopped at once (it does not wait for the zombie), cleanup prints
        Cleaned and the lock is released. Fails on revision 13 (running, 13, Stopping)."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        leader = self.zombie()
        start = self.lstart(leader)
        self.assertTrue(start, "ps gives the zombie's start time")
        name = record["run_id"] + "-zombie"
        self.forge(record, name, status="running", pid=leader, pgid=leader,
                   process_start=start, ended=None, exit_code=None)
        self.assertNotIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 0)
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[name]["status"], "unknown")
        self.end_worker()
        self.fake_ps = self.fake_ps_path("refused")
        try:
            result = self.tool("cleanup", BETA)
            self.refused(result, 17)
            self.assertNotIn("Cleaned", result.stdout)
            self.refused(self.tool("lock", "release"), 17)
            self.assertTrue((self.lock / "owner.txt").exists())
            result = self.tool("stop")
            self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
            self.assertNotIn("Stopped", result.stdout)
        finally:
            self.fake_ps = None
        result = self.ok(self.tool("stop"))
        self.assertTrue(result.stdout.strip().startswith("Stopped"), result.stdout)
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))
        self.ok(self.tool("lock", "release"))
        self.assertTrue(ps_state(leader).startswith("Z"), "the zombie was never collected")

    def test_zombie_leader_group_counted_while_member_lives(self):
        """Revision 14: a `running` record whose leader is a zombie but whose group still
        has a live member (the member's parent, the leader, has exited): the member is a
        survivor, so reconcile leaves the run running, readiness and a new run are refused,
        stop prints Stopping and sends it SIGTERM, cleanup gives Still alive (it ignores the
        signal). Once the member has exited, the zombie alone keeps nothing alive: reconcile
        marks the run unknown, stop prints Stopped and cleanup Cleaned. Fails on revision 13
        (cleanup refuses 13: the zombie leader reads as a running run).
        The member's safety deadline is 120 s (`long`): on a loaded Linux host the test took
        up to 80 s, and with 30 s the member could exit by itself before cleanup."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        leader, member, ready, release = self.zombie_leader_with_member("note+long")
        name = record["run_id"] + "-zombie"
        self.forge(record, name, status="running", pid=leader, pgid=leader,
                   process_start=self.lstart(leader), ended=None, exit_code=None)
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[name]["status"], "running")
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 13)
        self.refused(self.ready(BETA, "--packet", self.packet(self.tree())), 8)
        self.end_worker()
        self.stopping(self.tool("stop"))
        self.wait_for(Path(str(ready) + ".term").exists, "the member to get SIGTERM", 5)
        result = self.cleanup(10)
        self.assertTrue(result.stdout.strip().startswith("Still alive"), result.stdout)
        self.assertTrue(alive(member))
        release.write_text("")
        self.wait_for(lambda: not alive(member), "the member to exit")
        self.assertTrue(ps_state(leader).startswith("Z"), "the leader is still a zombie")
        self.ok(self.tool("reconcile"))
        self.assertEqual(self.runs()[name]["status"], "unknown")
        result = self.ok(self.tool("stop"))
        self.assertTrue(result.stdout.strip().startswith("Stopped"), result.stdout)
        self.assertTrue(self.cleanup(0).stdout.strip().startswith("Cleaned"))

    # ---- Revision 13 (item 10): transfer from auto to attended -------------------

    def successor(self, name="successor"):
        path = self.tmp / "attended" / name
        path.mkdir(parents=True)
        return path

    def test_transfer_from_blocked_is_final_and_offers_no_decision(self):
        """Revision 13 (item 10): `transfer JOB --to DIR --note TEXT` from Blocked makes the
        job Transferred: the report shows "transferred to DIR" and no decision; the decision
        and the note are kept in the record; the job cannot change afterwards (set, wait,
        criterion, candidate, worker start, transfer again); `set JOB Transferred` is not a
        way in; a missing successor directory is refused."""
        self.ok(self.set(ALPHA, "Preparing"))
        self.ok(self.set(ALPHA, "Blocked", "--note", "two fixes possible",
                         "--decision", "choose option a or b"))
        self.assertEqual(self.tool("set", BETA, "Transferred").returncode, 2)
        target = self.successor()
        result = self.tool("transfer", ALPHA, "--to", self.tmp / "no-such-dir", "--note", "x")
        self.refused(result, 2)
        self.assertEqual(self.state(ALPHA), "Blocked")
        result = self.ok(self.tool("transfer", ALPHA, "--to", target, "--note",
                                   "continued attended in the successor record"))
        self.assertIn(f"{ALPHA} Transferred to {target}", result.stdout)
        entry = job_entry(self.status(), ALPHA)
        self.assertEqual(entry["state"], "Transferred")
        self.assertEqual(entry["transfer"]["to"], str(target))
        self.assertEqual(entry["transfer"]["from_state"], "Blocked")
        self.assertEqual(entry["transfer"]["decision_needed_at_transfer"], "choose option a or b")
        self.assertNotIn("decision_needed", entry)
        report = self.report()
        self.assertEqual(self.result_cell(report, ALPHA), f"transferred to {target}")
        step = self.result_cell(report, ALPHA, column=2)
        self.assertNotIn("decision needed", step)
        self.assertNotIn("choose option", report)
        self.assertIn("Transferred to", (self.queue / "EVENTS.log").read_text())
        for args, code in ((("set", ALPHA, "Preparing"), 8), (("set", ALPHA, "Discarded"), 8),
                           (("wait", ALPHA, "--until", "2026-10-13T06:00:00Z"), 8),
                           (("criterion", ALPHA, "--file", self.criterion), 8),
                           (("candidate", ALPHA, "--repo", self.repo, "--commit", "HEAD"), 8),
                           (("worker", "start", ALPHA, "--id", "w-after"), 8),
                           (("transfer", ALPHA, "--to", target, "--note", "again"), 8)):
            with self.subTest(args=args[:2]):
                self.refused(self.tool(*args), code)
                self.assertEqual(self.state(ALPHA), "Transferred")

    def test_transfer_from_ready_keeps_the_packet_record(self):
        """Revision 13 (item 10): a Ready job can be transferred; its packet and manifest
        identity move into the transfer record (kept, no longer offered for approval) and
        the report offers no approval for it."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        self.ok(self.run_cmd(BETA, "unit"))
        packet = self.packet(self.tree())
        self.ok(self.ready(BETA, "--packet", packet))
        target = self.successor()
        self.ok(self.tool("transfer", BETA, "--to", target, "--note", "attended from here"))
        entry = job_entry(self.status(), BETA)
        self.assertEqual(entry["state"], "Transferred")
        self.assertEqual(entry["transfer"]["from_state"], "Ready for approval")
        self.assertEqual(entry["transfer"]["packet"], str(packet.resolve()))
        self.assertEqual(entry["transfer"]["packet_manifest_sha256"],
                         sha256_file(packet / "MANIFEST.sha256"))
        self.assertNotIn("packet", entry)
        report = self.report()
        self.assertNotIn("approve for integration", self.result_cell(report, BETA, column=2))
        self.refused(self.approve(BETA, f"approve {BETA}"), 8)
        self.assertTrue(packet.is_dir())

    def test_transfer_refused_while_worker_registered_and_all_finishes_queue(self):
        """Revision 13 (item 10): transfer from Preparing is refused (9) while a worker for
        the job is registered, and the job stays Preparing; after `worker end` it is accepted.
        `transfer --all` moves every job that is not final, and a queue whose jobs are all
        transferred is finished (another queue can be started in the same home)."""
        self.prepare(BETA)
        target = self.successor()
        result = self.tool("transfer", BETA, "--to", target, "--note", "handover")
        self.refused(result, 9)
        self.assertIn("still registered", result.stderr)
        self.assertEqual(self.state(BETA), "Preparing")
        self.refused(self.tool("init", *self.init_args("q-next", self.tmp / "RUN-next.lock")[1:],
                               queue=None), 8)
        self.end_worker()
        self.ok(self.tool("transfer", "--all", "--to", target, "--note", "handover"))
        self.assertEqual(self.state(BETA), "Transferred")
        self.assertEqual(self.state(ALPHA), "Transferred")
        self.refused(self.tool("transfer", "--all", "--to", target, "--note", "again"), 8)
        self.ok(self.tool("lock", "release"))
        self.ok(self.tool(*self.init_args("q-next", self.tmp / "RUN-next.lock"), queue=None))

    # ---- Revision 13 (item 13): timestamp-only installation differences ---------

    def write_fixture(self, path, data, writer=None):
        """Revision 15: write fixture bytes with `writer` (default: Path.write_bytes) and
        return what is then on disk. A step that reports success but leaves no file, or other
        bytes, is a test failure naming the file (revision 14 review finding 1: never
        FileNotFoundError)."""
        if path.exists():
            path.unlink()
        (writer or (lambda target, payload: target.write_bytes(payload)))(path, data)
        if not path.is_file():
            self.fail("the archive fixture step reported success but wrote no archive at %s"
                      % path)
        written = path.read_bytes()
        if written != data:
            self.fail("the archive fixture step wrote other bytes than intended at %s" % path)
        return written

    @staticmethod
    def archive_entries(data):
        """Revision 14: the members of an `!<arch>` archive as (header, name, content), the
        name read from BSD `#1/N` or GNU `name/` and `/N` forms; None if it does not parse."""
        if not data.startswith(b"!<arch>\n"):
            return None
        entries, pos, names = [], 8, b""
        while pos < len(data):
            header = data[pos:pos + 60]
            if len(header) < 60 or header[58:60] != b"`\n" or not header[48:58].strip().isdigit():
                return None
            size = int(header[48:58])
            body = data[pos + 60:pos + 60 + size]
            if len(body) < size:
                return None
            raw = header[:16].rstrip(b" ")
            if raw.startswith(b"#1/") and raw[3:].isdigit():
                length = int(raw[3:])
                name, content = body[:length].rstrip(b"\0"), body[length:]
            elif raw == b"//":
                names, name, content = body, raw, body
            elif raw.startswith(b"/") and raw[1:].isdigit() and names:
                name = names[int(raw[1:]):].split(b"/\n")[0]
                content = body
            else:
                name = raw[:-1] if raw.endswith(b"/") and raw != b"/" else raw
                content = body
            entries.append((header, name.decode("latin-1"), content))
            pos += 60 + size + (size % 2)
        return entries

    def archive_members_ok(self, data, members):
        """The archive's own members (index and name-table members left aside) are exactly
        `members`, in order, with these contents. macOS ar pads a member's content with up to
        seven newlines (inside its size) to keep 8-byte alignment; such padding is allowed."""
        entries = self.archive_entries(data)
        if entries is None:
            return False
        own = [(name, content) for _, name, content in entries
               if not (name.startswith("__.SYMDEF") or name in ("/", "//", "/SYM64"))]
        if [name for name, _ in own] != list(members):
            return False
        for (_, content), intended in zip(own, members.values()):
            rest = content[len(intended):]
            if not content.startswith(intended) or len(rest) > 7 or rest.strip(b"\n"):
                return False
        return True

    def dates_only_differ(self, first, second):
        """The two archives differ, and only in member header date fields (bytes 16-28)."""
        a, b = self.archive_entries(first), self.archive_entries(second)
        if first == second or a is None or b is None or len(a) != len(b):
            return False
        return all(x[0][:16] + x[0][28:] == y[0][:16] + y[0][28:] and x[2] == y[2]
                   for x, y in zip(a, b))

    def timestamp_archives(self, path, members, changed, same_size, style="gnu"):
        """Revision 15: (how, original, redated, changed, same_size) archive bytes for the
        timestamp test, written by the test (ar_archive, in the `style` layout; no ar or
        ranlib): the members dated 1000000000 (index 1000000007), the same members dated
        1100000000 (index 1100000007), and `changed` and `same_size` (a member changed
        without changing its size) dated 1100000000. Each is written to `path` and read back
        (write_fixture) and checked byte by byte: it holds exactly the intended members and
        contents, its index member comes first, and the re-dated one differs from the
        original in member header dates only. Any other result is a test failure that says
        which archive is wrong."""
        archives = []
        for contents, when in ((members, AR_DATE), (members, AR_REDATE), (changed, AR_REDATE),
                               (same_size, AR_REDATE)):
            archives.append(self.write_fixture(path, ar_archive(contents, when, style,
                                                                index_when=when + 7)))
        original, redated, after, same = archives
        for label, data, contents in (("original", original, members),
                                      ("re-dated", redated, members),
                                      ("content-changed", after, changed),
                                      ("same-size changed", same, same_size)):
            entries = self.archive_entries(data)
            if not self.archive_members_ok(data, contents):
                self.fail("the %s %s archive does not hold exactly the intended members"
                          % (label, style))
            index = "__.SYMDEF SORTED" if style == "bsd" else "/"
            if entries[0][1] != index:
                self.fail("the %s %s archive does not start with its index %r"
                          % (label, style, index))
        if not self.dates_only_differ(original, redated):
            self.fail("the re-dated %s archive differs from the original in more than the "
                      "member dates, or not at all" % style)
        return ("archives written by the test (%s layout)" % style, original, redated, after,
                same)

    def normalized_fingerprint_script(self, archive):
        """A fingerprint script printing the archive's exact digest and its normalized
        digest from payload/tools/fingerprint_normalize.py."""
        self.counter += 1
        path = self.work / f"fingerprint-ar-{self.counter}.py"
        normalizer = TOOL.parent / "payload" / "tools" / "fingerprint_normalize.py"
        path.write_text(f"""#!{sys.executable}
import hashlib, subprocess, sys
data = open({str(archive)!r}, "rb").read()
print("FINGERPRINT " + hashlib.sha256(data).hexdigest())
out = subprocess.run([sys.executable, "-I", "-B", {str(normalizer)!r}, "file", "--kind",
                      "archive", {str(archive)!r}], capture_output=True, text=True, check=True)
print("NORMALIZED " + out.stdout.split()[0])
""")
        path.chmod(0o755)
        return path

    def exact(self, script):
        out = subprocess.run([str(script)], capture_output=True, text=True, check=True).stdout
        return dict(line.split() for line in out.splitlines())

    def normalize(self, path, kind="archive"):
        """Revision 15: (digest, stderr) of payload/tools/fingerprint_normalize.py for one file."""
        normalizer = TOOL.parent / "payload" / "tools" / "fingerprint_normalize.py"
        out = subprocess.run([sys.executable, "-I", "-B", str(normalizer), "file", "--kind", kind,
                              str(path)], capture_output=True, text=True, check=True)
        return out.stdout.split()[0], out.stderr

    def test_install_check_classifies_timestamp_only_archive_change(self):
        """Revision 13 (item 13): for an archive, an exact match is available
        and records the normalized value; the same members re-archived with other dates
        change the exact fingerprint, and install check then says "timestamp-only
        differences: attended validation needed" but keeps the installation unavailable
        (12; runs refused; the recorded value is not replaced); a member whose content
        changed is not classified; restoring the exact bytes makes it available again.
        Revision 14 (review): also not a member whose content changed but not its size.
        Revision 15: the archives are written by the test (timestamp_archives) in both real
        layouts, GNU/SysV (with a `/` symbol table, a `//` long-name table and odd-size
        padding) and BSD/macOS (with a dated `__.SYMDEF SORTED` index, `#1/N` names and
        alignment padding), each on its own queue, on every platform; ar and ranlib are not
        used."""
        members = {"one.txt": b"first member\n", "two.txt": b"second member\n",
                   "a_member_with_a_long_name.o": b"long-named member\n"}
        changed_members = dict(members, **{"two.txt": b"second member, changed\n"})
        same_size_members = dict(members, **{"two.txt": b"SECOND member\n"})
        self.assertEqual(len(same_size_members["two.txt"]), len(members["two.txt"]))
        for style in ("gnu", "bsd"):
            with self.subTest(layout=style):
                if style != "gnu":
                    self.use_new_queue("q-archive-" + style)
                self.check_timestamp_only_archive(style, members, changed_members,
                                                  same_size_members)

    def check_timestamp_only_archive(self, style, members, changed_members, same_size_members):
        """The body of the test above for one archive layout, on the current queue."""
        self.prepare(BETA)
        archive = self.work / f"libsynthetic-{style}.a"
        how, original, redated, changed, same_size = self.timestamp_archives(
            archive, members, changed_members, same_size_members, style)
        archive.write_bytes(original)
        script = self.normalized_fingerprint_script(archive)
        values = self.exact(script)
        expect = values["FINGERPRINT"]
        self.ok(self.tool("install", "check", "--script", script, "--expect", expect))
        shown = json.loads(self.ok(self.tool("install", "show")).stdout)
        self.assertEqual(shown["state"], "available")
        self.assertEqual(shown["expected_normalized"], values["NORMALIZED"])
        archive.write_bytes(redated)
        self.assertNotEqual(self.exact(script)["FINGERPRINT"], expect, how)
        self.assertEqual(self.exact(script)["NORMALIZED"], values["NORMALIZED"], how)
        result = self.tool("install", "check", "--script", script, "--expect", expect)
        self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
        shown = json.loads(self.ok(self.tool("install", "show")).stdout)
        self.assertEqual(shown["state"], "unavailable")
        self.assertTrue(shown["reason"].startswith(
            "timestamp-only differences: attended validation needed"), shown["reason"])
        self.assertEqual(shown["expected"], expect)
        self.attempt(BETA, 12)
        self.assertIn("timestamp-only differences", self.report())
        archive.write_bytes(same_size)
        self.assertNotEqual(self.exact(script)["NORMALIZED"], values["NORMALIZED"], how)
        archive.write_bytes(changed)
        self.assertNotEqual(self.exact(script)["NORMALIZED"], values["NORMALIZED"], how)
        result = self.tool("install", "check", "--script", script, "--expect", expect)
        self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
        shown = json.loads(self.ok(self.tool("install", "show")).stdout)
        self.assertEqual(shown["reason"], "it no longer matches its recorded starting state")
        archive.write_bytes(original)
        self.ok(self.tool("install", "check", "--script", script, "--expect", expect))
        self.assertEqual(json.loads(self.ok(self.tool("install", "show")).stdout)["state"],
                         "available")
        self.attempt(BETA, 0)

    def test_normalizer_archive_fields_and_corruption(self):
        """Revision 15: the normalizer on archives written by the test, in both layouts. Only
        member header dates are ignored: the same members re-dated (the index too) give the
        same value; a changed member name, owner, index entry or member content (also of the
        same size), and a reordered member list, give another. Damaged archives (cut inside
        the last member, a short member header, bad header magic, a member size past the end,
        a size that is not a number, bad archive magic) are compared byte for byte, with
        "not a parsable archive" on stderr, so a damaged archive re-dated is never
        timestamp-only. Control: each undamaged pair has the same value."""
        members = {"one.txt": b"first member\n", "two.txt": b"second member\n",
                   "a_member_with_a_long_name.o": b"long-named member\n"}
        path = self.work / "fields.a"
        for style in ("gnu", "bsd"):
            with self.subTest(layout=style):
                base = ar_archive(members, AR_DATE, style, index_when=AR_DATE + 7)
                redated = ar_archive(members, AR_REDATE, style, index_when=AR_REDATE + 9)
                self.assertTrue(self.archive_members_ok(base, members))
                self.assertTrue(self.dates_only_differ(base, redated))
                self.write_fixture(path, base)
                reference, note = self.normalize(path)
                self.assertNotIn("not a parsable archive", note)
                self.assertNotEqual(reference, sha256_file(path))
                self.write_fixture(path, redated)
                self.assertEqual(self.normalize(path)[0], reference, "re-dated")
                reordered = dict(reversed(list(members.items())))
                variants = {
                    "member name": base.replace(b"two.txt", b"twX.txt", 1),
                    "owner": ar_archive(members, AR_DATE, style, index_when=AR_DATE + 7, uid=501),
                    "member order": ar_archive(reordered, AR_DATE, style, index_when=AR_DATE + 7),
                    "same-size content": ar_archive(dict(members, **{"two.txt": b"SECOND member\n"}),
                                                    AR_DATE, style, index_when=AR_DATE + 7),
                    "index entry": base.replace(b"sym_one", b"sym_ONE", 1)}
                for label in ("member name", "index entry"):  # the same bytes but one
                    self.assertEqual(len(variants[label]), len(base), label)
                for label, data in variants.items():
                    self.assertNotEqual(data, base, label)
                    self.write_fixture(path, data)
                    digest, note = self.normalize(path)
                    self.assertNotIn("not a parsable archive", note, label)
                    self.assertNotEqual(digest, reference, label)
                damaged_redated = ar_corruptions(redated)
                for label, data in ar_corruptions(base).items():
                    self.write_fixture(path, data)
                    digest, note = self.normalize(path)
                    self.assertIn("not a parsable archive; compared byte for byte", note, label)
                    self.assertEqual(digest, hashlib.sha256(data).hexdigest(), label)
                    self.write_fixture(path, damaged_redated[label])
                    self.assertNotEqual(self.normalize(path)[0], digest,
                                        "a damaged archive re-dated counted as timestamp-only: "
                                        + label)

    def test_archive_step_without_output_never_ends_in_file_not_found(self):
        """Revision 15 (revision 14 review, finding 1, the outside reader): a tool or step
        that reports success without writing the archive must lead to the stated fallback or
        a clear failure, never FileNotFoundError. (a) The reader's case: ar and ranlib
        stand-ins first on PATH that exit 0 and write nothing. With 1664dd1f's test code
        this ended in FileNotFoundError; now the archives are the intended ones (no outside
        tool is used). (b) The fixture's own write step replaced by one that returns without
        writing: a test failure that names the missing archive."""
        fakebin = self.work / "fake-ar-bin"
        fakebin.mkdir()
        for name in ("ar", "ranlib"):
            (fakebin / name).write_text("#!%s\nimport sys\nsys.exit(0)\n" % sys.executable)
            (fakebin / name).chmod(0o755)
        self.env["PATH"] = str(fakebin) + os.pathsep + self.env.get("PATH", os.defpath)
        members = {"one.txt": b"first member\n", "two.txt": b"second member\n"}
        changed = dict(members, **{"two.txt": b"second member, changed\n"})
        same_size = dict(members, **{"two.txt": b"SECOND member\n"})
        try:
            how, original, redated, after, same = self.timestamp_archives(
                self.work / "lib-no-output.a", members, changed, same_size)
        except FileNotFoundError as error:
            self.fail("an archive step that wrote nothing ended in FileNotFoundError: %s" % error)
        self.assertTrue(self.archive_members_ok(original, members), how)
        self.assertTrue(self.archive_members_ok(after, changed), how)
        self.assertTrue(self.archive_members_ok(same, same_size), how)
        self.assertTrue(self.dates_only_differ(original, redated), how)
        target = self.work / "never-written.a"
        with self.assertRaises(self.failureException) as caught:
            self.write_fixture(target, original, writer=lambda path, data: None)
        self.assertIn("wrote no archive at %s" % target, str(caught.exception))

    def test_install_check_expect_normalized_and_known_exit_codes(self):
        """Revision 13 (items 13, 14): on a queue with no stored normalized value,
        --expect-normalized gives the value to compare; the classification still leaves the
        installation unavailable. Exit 3 and exit 4 of a fingerprint script are named."""
        script = self.fingerprint_script("abc123")
        script.write_text(script.read_text().replace(
            'print("FINGERPRINT abc123")', 'print("FINGERPRINT abc123")\nprint("NORMALIZED n1")'))
        result = self.tool("install", "check", "--script", script, "--expect", "zzz999",
                           "--expect-normalized", "n1")
        self.assertEqual(result.returncode, 12, result.stdout + result.stderr)
        shown = json.loads(self.ok(self.tool("install", "show")).stdout)
        self.assertEqual(shown["state"], "unavailable")
        self.assertIn("timestamp-only differences", shown["reason"])
        result = self.tool("install", "check", "--script", script, "--expect", "zzz999",
                           "--expect-normalized", "n2")
        self.assertEqual(result.returncode, 12)
        self.assertNotIn("timestamp-only", self.ok(self.tool("install", "show")).stdout)
        for code, words in ((3, "a process runs from the installation"),
                            (4, "the process table could not be read")):
            with self.subTest(code=code):
                bad = self.fingerprint_script("zzz999", code)
                self.assertEqual(self.tool("install", "check", "--script", bad, "--expect",
                                           "zzz999").returncode, 12)
                reason = json.loads(self.ok(self.tool("install", "show")).stdout)["reason"]
                self.assertEqual(reason, f"its fingerprint script failed (exit {code}: {words})")

    def test_normalizer_build_time_masks_date_and_derived_fields(self):
        """Revision 13 (item 13): a bundle built twice from source differing only in an
        embedded asctime build date has a different exact digest and the same build-time
        digest (date, LC_UUID and code signature masked); a different other string changes
        it. A non-archive given as an archive is compared byte for byte.
        Revision 15: the bundles are written by the test (macho_bundle; formerly built with
        cc on macOS only), so this runs on every platform. Controls: the UUID and the
        signature range really differ between the two dates (so the masking is needed); a
        changed code byte, or one changed byte in the load commands (also LC_DYSYMTAB,
        LC_BUILD_VERSION and LC_LOAD_DYLIB, review finding 2), __cstring, __unwind_info, __DATA,
        or __LINKEDIT before the signature (also the fixups, exports trie, function starts,
        symbol and string tables), changes the digest; a changed byte inside the UUID or the
        signature range does not. The fixture's load commands are those of cc output, in
        order (MACHO_CC_COMMANDS)."""
        a = self.work / "a.so"
        a.write_bytes(macho_bundle())
        b = self.work / "b.so"
        b.write_bytes(macho_bundle(date="Tue Oct 13 01:02:03 2026"))
        c = self.work / "c.so"
        c.write_bytes(macho_bundle(extra="phasex"))
        seven = self.work / "seven.so"
        seven.write_bytes(macho_bundle(value=7))
        ranges = macho_ranges(a.read_bytes())
        for name in ("uuid", "signature"):
            begin, end = ranges[name]
            self.assertNotEqual(a.read_bytes()[begin:end], b.read_bytes()[begin:end], name)
        self.assertNotEqual(sha256_file(a), sha256_file(b))
        digest = self.normalize(a, "build-time")[0]
        self.assertEqual(self.normalize(b, "build-time")[0], digest)
        self.assertNotEqual(self.normalize(c, "build-time")[0], digest)
        self.assertNotEqual(self.normalize(seven, "build-time")[0], digest, "a code change")
        self.assertEqual(self.normalize(a, "archive")[0], sha256_file(a))
        original = a.read_bytes()
        changed = self.work / "changed.so"
        date_at = original.index(MACHO_DATE.encode())
        self.assertEqual(macho_load_commands(original), MACHO_CC_COMMANDS)
        for name, masked in (("text", False), ("cstring", False), ("unwind_info", False),
                             ("data", False), ("linkedit", False), ("fixups", False),
                             ("exports trie", False), ("function starts", False),
                             ("symbol table", False), ("string table", False),
                             ("dysymtab command", False), ("build version command", False),
                             ("load dylib command", False), ("commands", False),
                             ("uuid", True), ("signature", True)):
            begin, end = ranges[name]
            at = date_at + len(MACHO_DATE) + 2 if name == "cstring" else (begin + end) // 2
            with self.subTest(byte=name):
                data = bytearray(original)
                data[at] ^= 0x01
                changed.write_bytes(bytes(data))
                if masked:
                    self.assertEqual(self.normalize(changed, "build-time")[0], digest)
                else:
                    self.assertNotEqual(self.normalize(changed, "build-time")[0], digest)

    # ---- Revision 15 (review finding 4): real-tool conformance, where the tools exist ----------

    def real_tool(self, *names):
        """The paths of `names` on the test's PATH, or None when one is missing."""
        tools = [shutil.which(name, path=self.env.get("PATH")) for name in names]
        return tools if all(tools) else None

    @staticmethod
    def ar_family(ar):
        """Revision 15 (second reading, finding 1): "apple" for Apple (cctools) ar, whose
        `--version` prints its usage, "gnu" for GNU binutils ar, None for any other ar."""
        try:
            out = subprocess.run([ar, "--version"], capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.TimeoutExpired):
            return None
        text = out.stdout + out.stderr
        if "GNU ar" in text or "GNU Binutils" in text:
            return "gnu"
        if text.lstrip().startswith("usage:") and "ar -d" in text:
            return "apple"
        return None

    def test_real_tool_conformance_archive(self):
        """Revision 15 (review finding 4): where `ar` and `ranlib` are on PATH, an archive they
        write must have what the normalizer and the fixtures rely on: it parses as `!<arch>`
        (no "not a parsable archive"); its members hold exactly the intended contents; whether
        the fixture of the matching layout (BSD `#1/N` names: ar_archive "bsd"; otherwise
        "gnu") has the same member name fields, sizes and contents (the index member aside: its
        contents depend on the symbols) is reported on stderr (revision 15, second reading,
        finding 1: another valid `ar`, for example llvm-ar, pads names differently); for every
        `ar`, re-dating every member header
        of the real archive keeps the normalized value and changes the exact digest, and one
        changed content byte changes the normalized value. These are the only uses of real
        tools in the kit's tests; without the tools the test is skipped with that reason, and
        the fixture tests above still run. A tool that fails or writes no archive is reported
        as a skip with its reason, never as FileNotFoundError. Revision 15 follow-up (GitHub's
        macOS runner): tools whose archive does not hold exactly the intended members (an `ar`
        or `ranlib` that drops or reorders members that are not object files) are reported as
        a skip naming the members found, and a layout that differs from the fixture's is only
        noted on stderr, for every `ar`; the normalizer checks on the real archive stay strict."""
        tools = self.real_tool("ar", "ranlib")
        if not tools:
            self.skipTest("real-tool conformance check skipped: ar and ranlib are not both on "
                          "PATH (the archive fixture tests do not need them)")
        ar, ranlib = tools
        members = {"one.txt": b"first member\n", "two.txt": b"second member\n",
                   "a_member_with_a_long_name.o": b"long-named member\n"}
        build = self.work / "real-ar"
        build.mkdir()
        for name, data in members.items():
            (build / name).write_bytes(data)
            os.utime(build / name, (AR_DATE, AR_DATE))
        archive = build / "libreal.a"
        for argv in ([ar, "rc", str(archive), *members], [ranlib, str(archive)]):
            done = subprocess.run(argv, cwd=build, capture_output=True, env=self.env)
            if done.returncode != 0 or not archive.is_file():
                self.skipTest("real-tool conformance check skipped: %s exited %d%s" % (
                    argv[0], done.returncode, "" if archive.is_file() else " and wrote no archive"))
        real = archive.read_bytes()
        print("\n[conformance] real archive from %s and %s (%d bytes)" % (ar, ranlib, len(real)),
              file=sys.stderr)
        entries = self.archive_entries(real)
        self.assertIsNotNone(entries, "the real archive does not parse as !<arch>")
        if not self.archive_members_ok(real, members):
            self.skipTest("real-tool conformance check skipped: %s and %s wrote an archive "
                          "without exactly the intended members (found %s)" % (
                              ar, ranlib, [(e[1], len(e[2])) for e in entries]))
        digest, note = self.normalize(archive)
        self.assertNotIn("not a parsable archive", note)
        style = "bsd" if any(e[0][:3] == b"#1/" for e in entries) else "gnu"
        fixture = self.archive_entries(ar_archive(members, AR_DATE, style))

        def own(items):
            return [(e[0][:16], e[0][48:58], e[2]) for e in items
                    if not (e[1].startswith("__.SYMDEF") or e[1] in ("/", "/SYM64"))]
        print("[conformance] %s ar (%s): %s with the %s fixture" % (
            self.ar_family(ar) or "other", ar, "same member layout" if own(fixture) == own(entries)
            else "member layout DIFFERS", style), file=sys.stderr)
        redated, pos = bytearray(real), 8
        while pos < len(redated):
            size = int(redated[pos + 48:pos + 58])
            redated[pos + 16:pos + 28] = b"%-12d" % (AR_REDATE + pos)
            pos += 60 + size + size % 2
        self.write_fixture(archive, bytes(redated))
        self.assertNotEqual(sha256_file(archive), hashlib.sha256(real).hexdigest())
        self.assertEqual(self.normalize(archive)[0], digest, "re-dated real archive")
        changed = bytearray(real)
        changed[real.index(b"second member")] ^= 0x01
        self.write_fixture(archive, bytes(changed))
        self.assertNotEqual(self.normalize(archive)[0], digest, "changed content byte")

    def test_real_tool_conformance_macho_bundle(self):
        """Revision 15 (review finding 4): on macOS where `cc` is on PATH, bundles it builds
        (`cc -bundle`) must have what the normalizer and the fixture rely on: an LC_UUID and
        one LC_CODE_SIGNATURE whose range ends at the end of __LINKEDIT and of the file, after
        the load commands; two builds differing only in the asctime date differ exactly and
        have the same build-time digest, with nothing reported unmasked; a code change
        (`return 7`) changes it. A load command sequence other than the fixture's
        (MACHO_CC_COMMANDS) is reported on stderr, so the fixture can be brought up to date; it
        does not fail the test. Skipped, with the reason, off macOS, without cc, or when cc
        cannot build a bundle."""
        cc = self.real_tool("cc")
        if sys.platform != "darwin" or not cc:
            self.skipTest("real-tool conformance check skipped: needs macOS and cc on PATH "
                          "(the Mach-O fixture tests do not need them)")

        def bundle(name, date, value=1):
            source = self.work / f"{name}.c"
            source.write_text(f'const char *build_time = "{date}";\nconst char *what = "phaser";\n'
                              f"int version(void) {{ return {value}; }}\n")
            out = self.work / f"{name}.so"
            done = subprocess.run([cc[0], "-bundle", "-o", str(out), str(source)],
                                  capture_output=True, text=True, env=self.env)
            if done.returncode != 0 or not out.is_file():
                self.skipTest("real-tool conformance check skipped: cc -bundle exited %d: %s"
                              % (done.returncode, done.stderr.strip()[:200]))
            return out
        a = bundle("ra", MACHO_DATE)
        b = bundle("rb", "Tue Oct 13 01:02:03 2026")
        seven = bundle("r7", MACHO_DATE, 7)
        data = a.read_bytes()
        types = macho_load_commands(data)
        print("\n[conformance] cc bundle load commands: %s" % " ".join(hex(t) for t in types),
              file=sys.stderr)
        if types != MACHO_CC_COMMANDS:
            print("[conformance] NOTE: differs from the fixture's MACHO_CC_COMMANDS %s"
                  % " ".join(hex(t) for t in MACHO_CC_COMMANDS), file=sys.stderr)
        self.assertEqual(types.count(0x1b), 1, "LC_UUID")
        self.assertEqual(types.count(0x1d), 1, "LC_CODE_SIGNATURE")
        ranges = macho_ranges(data)
        sizeofcmds = struct.unpack_from("<I", data, 20)[0]
        begin, end = ranges["signature"]
        self.assertEqual(end, len(data), "the signature does not end at the end of the file")
        self.assertGreaterEqual(begin, 32 + sizeofcmds)
        self.assertGreaterEqual(begin, ranges["linkedit"][0])
        self.assertNotEqual(sha256_file(a), sha256_file(b))
        (first, note), (second, other) = self.normalize(a, "build-time"), self.normalize(b, "build-time")
        self.assertNotIn("not masked", note + other)
        self.assertNotIn("no build-time date found", note + other)
        self.assertEqual(first, second, "two cc builds differing only in the date")
        self.assertNotEqual(self.normalize(seven, "build-time")[0], first, "a code change")

    def test_normalizer_masks_code_signature_only_at_the_end(self):
        """Checkpoint 3 (REVIEW_CP2 finding 4): two bundles that differ in code, whose
        LC_CODE_SIGNATURE load command is rewritten to declare the whole file
        (dataoff 0, datasize = file length), keep different build-time digests (the
        declared range is not masked, and stderr says so). Control: unmodified, the same
        two bundles also differ, and a signature at the end is still masked (see the
        test above). Checkpoint 5 (review CP3 finding 6): the same holds for a range
        from the end of the load commands to the end of the file, and for that range
        with __LINKEDIT also rewritten to cover it (it would overlap __TEXT); such a
        file is compared byte for byte. Revision 15: the bundles are written by the test
        (macho_bundle, code `mov w0, #1` and `mov w0, #7`; formerly built with cc on macOS
        only), so this runs on every platform."""

        def bundle(name, value):
            out = self.work / f"{name}.so"
            out.write_bytes(macho_bundle(value=value))
            return out

        def widen(path, how):
            data = bytearray(path.read_bytes())
            ncmds, sizeofcmds = struct.unpack_from("<II", data, 16)
            begin = 0 if how == "whole" else 32 + sizeofcmds
            offset, found = 32, False
            for _ in range(ncmds):
                cmd, size = struct.unpack_from("<II", data, offset)
                if cmd == 0x1d:  # LC_CODE_SIGNATURE
                    struct.pack_into("<II", data, offset + 8, begin, len(data) - begin)
                    found = True
                elif (cmd == 0x19 and how == "linkedit"
                      and data[offset + 8:offset + 24].rstrip(b"\0") == b"__LINKEDIT"):
                    struct.pack_into("<QQ", data, offset + 40, begin, len(data) - begin)
                offset += size
            self.assertTrue(found, "the bundle has no LC_CODE_SIGNATURE")
            path.write_bytes(bytes(data))

        def digest(path):
            return self.normalize(path, "build-time")
        for how in ("whole", "after-load-commands", "linkedit"):
            with self.subTest(range=how):
                one, seven = bundle(f"b1-{how}", 1), bundle(f"b7-{how}", 7)
                self.assertNotEqual(digest(one)[0], digest(seven)[0])  # control
                widen(one, how)
                widen(seven, how)
                self.assertNotEqual(sha256_file(one), sha256_file(seven))
                (first, note), (second, _) = digest(one), digest(seven)
                self.assertNotEqual(first, second, "a widened signature range hid a code change")
                self.assertIn("LC_CODE_SIGNATURE not masked", note)
                self.assertEqual(first, sha256_file(one))  # compared byte for byte


if __name__ == "__main__":
    unittest.main()
