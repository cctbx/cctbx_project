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

The tool is only ever run through subprocess (-I -B); nothing here imports it
or reads its source. Every queue, control home, lock directory, space path,
allowed root and git repository lives under a temporary directory; --home and
--lock-dir are always passed, so the default home and lock are never touched.
Each refusal is paired with a positive control (the same call without the
forbidden condition is accepted).
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
QID = "q1"
ALPHA, BETA = "job-alpha", "job-beta"
JOBS = [{"id": ALPHA, "title": "Alpha", "requires": ["unit", "build"], "depends_on": []},
        {"id": BETA, "title": "Beta", "requires": ["unit"], "depends_on": []}]
REPORT_HEADER = "| Job | Result | Decision or next step |"
WAITING_STEP = ("resume with /gc auto resume after the reset "
                "(automatic continuation is not supported)")
TIMEOUT = 15
# A CMD (or grandchild) that signals readiness and exits 0 once the release file
# appears (safety deadline 30 s). mode is `+`-joined flags: `ignore`: SIGTERM
# ignored; `exit0`: SIGTERM exits 0 (a test script that traps TERM); `note`:
# SIGTERM writes <ready>.term and is otherwise ignored; `pgrp`: move to a new
# process group (same session); `setsid`: new session; `default`: SIGTERM kills it.
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
deadline = time.time() + 30
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
# A fingerprint script for `install check`.
FINGERPRINT_SCRIPT = """#!{python}
import sys
print("FINGERPRINT {fingerprint}")
sys.exit({code})
"""
# The libtbx.easy_run pattern: CMD starts a shell in a new session, waits 2 s, exits.
EASY_RUN_SCRIPT = """import os, subprocess, sys, time
shell = subprocess.Popen(["/bin/sh", "-c", sys.argv[2]], preexec_fn=os.setsid)
with open(sys.argv[1] + ".tmp", "w") as handle:
    handle.write(str(shell.pid))
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
time.sleep(2)
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
# then runs /bin/sleep 30.
DETACH_PROBE = """import os, subprocess, sys
child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(25)"],
                         start_new_session=True)
with open(sys.argv[1] + ".tmp", "w") as handle:
    handle.write(str(child.pid))
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
subprocess.run(["/bin/sleep", "30"])
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


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


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


class AutoQueueChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.tmp = Path(os.path.realpath(self.temp.name))
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

    def command(self, args, home, queue):
        command = [sys.executable, "-I", "-B", str(TOOL), "--home", str(home or self.home)]
        if queue is True:
            queue = self.queue
        if queue:
            command += ["--queue", str(queue)]
        return command + [str(a) for a in args]

    def tool(self, *args, home=None, queue=True, env=None):
        return subprocess.run(self.command(args, home, queue), cwd=self.tmp,
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
        wrapper lives."""
        self.prepare(BETA)
        started = self.work / "probe-started"
        marker, cmd = self.marker_cmd("cmd-ran")
        probe = f"touch {shlex.quote(str(started))}; sleep 2; echo {shlex.quote(str(self.src))}"
        proc = self.popen_run(BETA, cmd, probe=probe)
        self.wait_for(started.exists, "the probe to start")
        stem, record = self.wait_record(set(), status="starting")
        self.assertEqual(record["wrapper_pid"], proc.pid)
        self.end_worker()
        first = self.tool("stop")
        if first.returncode == 0:
            self.assertIsNotNone(proc.poll(), "Stopped while the starting wrapper lived")
        else:
            self.stopping(first)
        self.assertEqual(proc.wait(TIMEOUT), 3)
        self.assertFalse(marker.exists())
        self.assertEqual(self.runs()[stem]["status"], "not started")
        self.stop_until_stopped()

    def test_stop_from_different_tz_locale_and_path(self):
        """Revision 2 B1: stop and reconcile run with another TZ, LC_ALL and PATH than
        the run still recognise the live run; stop signals it."""
        self.prepare(BETA)
        proc = self.popen_run(BETA, ["sleep", "5"], env={"TZ": "UTC", "LC_ALL": "C"})
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
        proc = self.popen_run(BETA, ["sleep", "5"])
        stem, record = self.wait_record(set())
        self.end_worker()
        time.sleep(2.1)  # the decoy must start in a later second than the tracked run
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
        proc = self.popen_run(BETA, ["sleep", "5"])
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

    def test_survivor_setsid_found_by_run_id_marker(self):
        """Revision 4 fix 1 (c): a non-system program calling os.setsid() is found by
        its GC_AUTO_RUN_ID environment marker."""
        self.check_survivor_counted(*self.survivor_after_wrapper_kill("note+setsid"))

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

    def easy_run(self, shell_command):
        """Run CMD in the libtbx.easy_run pattern; returns the shell's pid."""
        script, pidfile = self.work / "easy_run.py", self.work / f"shell-{self.counter}.pid"
        script.write_text(EASY_RUN_SCRIPT)
        self.ok(self.run_cmd(BETA, "unit", cmd=[sys.executable, str(script), str(pidfile),
                                                 shell_command]))
        shell = int(pidfile.read_text())
        self.stray.append(shell)
        return shell

    def test_easy_run_pattern_shell_survivor(self):
        """Revision 5 S2: a shell started with setsid by CMD (python, then /bin/sleep 41)
        keeps the run from passing; report says running; stop signals the shell and
        the sleep; no sleep 41 remains after Stopped."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        python = shlex.quote(sys.executable)
        shell = self.easy_run(f"{python} -c 'import time; time.sleep(4)'; /bin/sleep 41")
        self.assertTrue(alive(shell))
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("running", cell.lower())
        self.assertNotIn("passed", cell.lower())
        self.refused(self.ready(BETA, "--packet", self.packet(self.tree())), 8)
        sleeper = []

        def sleep_started():
            table = subprocess.run(["/bin/ps", "-A", "-o", "pid=,ppid=,command="], text=True,
                                   capture_output=True).stdout
            for line in table.splitlines():
                pid, ppid, command = line.split(None, 2)
                if int(ppid) == shell and command.strip().endswith("sleep 41"):
                    sleeper.append(int(pid))
                    return True
            return False
        self.wait_for(sleep_started, "the shell's /bin/sleep 41")
        self.stray.extend(sleeper)
        cell = self.result_cell(self.report(), BETA)
        self.assertIn("running", cell.lower())
        self.assertNotIn("passed", cell.lower())
        self.end_worker()
        first = self.tool("stop")
        self.assertIn(first.returncode, (0, 10), first.stdout + first.stderr)
        self.stop_until_stopped()
        self.assertFalse(alive(shell), "the shell survived Stopped")
        self.assertFalse(alive(sleeper[0]), "sleep 41 survived Stopped")

    def test_easy_run_pattern_quick_shell_passes(self):
        """Revision 5 S2 control: the same pattern with a shell that ends quickly
        completes and passes."""
        self.prepare(BETA)
        self.ok(self.set(BETA, "Testing"))
        shell = self.easy_run(f"{shlex.quote(sys.executable)} -c 'pass'; /bin/sleep 0")
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
        after the next whole UTC second following `ended`."""
        self.check_ended_bound(lambda ended: ended + 1.05)

    def utc(self, iso):
        return calendar.timegm(time.strptime(iso, "%Y-%m-%dT%H:%M:%SZ"))

    def start_epoch(self, pid):
        lstart = subprocess.run(["/bin/ps", "-o", "lstart=", "-p", str(pid)],
                                env={"LC_ALL": "C", "TZ": "UTC"}, text=True,
                                capture_output=True).stdout.strip()
        return calendar.timegm(time.strptime(lstart, "%a %b %d %H:%M:%S %Y"))

    def orphan_of_dead_leader(self):
        """An exited, unrelated session leader P and its orphan (sid = pgid = P,
        ppid 1) that notes SIGTERM. Returns (P, orphan pid, ready, orphan start)."""
        cmd, ready, _ = self.wait_cmd("note")
        pidfile = self.work / f"orphan-{self.counter}.pid"
        shell = subprocess.Popen(["/bin/sh", "-c", " ".join(shlex.quote(c) for c in cmd)
                                  + f" & echo $! > {shlex.quote(str(pidfile))}; exit 0"],
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
                                             text=True, capture_output=True).stdout.strip() == "1",
                      "the orphan to be reparented to launchd")
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
        session are trusted); P's orphaned child (sid = pgid = P, ppid 1) started
        after the record's `ended` is not counted or signalled, and a new run is
        accepted. Control: the same record without `ended` counts it (13, Stopping,
        SIGTERM)."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        ended = self.utc(record["ended"])
        time.sleep(max(0.0, start_at(ended) - time.time()))
        leader, orphan, ready, start = self.orphan_of_dead_leader()
        self.assertGreater(start, ended)
        excluded = dict(pid=leader, pgid=leader, process_start=None)
        self.check_counted_only_when(record, excluded, dict(excluded, ended=None,
                                                            status="running"), orphan, ready)

    def test_probe_window_bounds_probe_leader(self):
        """Revision 6 S1 variant: a running record's probe_pgid (an exited, unrelated
        session leader) claims its orphan only if the orphan started within
        [probe_started, probe_ended]."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        (record,) = self.runs().values()
        probe_ended = self.utc(record["probe_ended"])
        self.assertLessEqual(self.utc(record["probe_started"]), probe_ended)
        time.sleep(max(0.0, probe_ended + 1.05 - time.time()))
        leader, orphan, ready, start = self.orphan_of_dead_leader()
        self.assertGreater(start, probe_ended)
        excluded = dict(status="running", ended=None, probe_pgid=leader, probe_start=None)
        included = dict(excluded, probe_ended=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                            time.gmtime(start)))
        self.check_counted_only_when(record, excluded, included, orphan, ready)

    def test_shell_parent_of_marked_process_not_signalled(self):
        """Revision 6 M3: a bash started by the user during a run, whose child carries
        the marker text only in its arguments, is never signalled by stop."""
        self.prepare(BETA)
        cmd, ready, release = self.wait_cmd("default")
        proc = self.popen_run(BETA, cmd)
        stem, _ = self.wait_record(set())
        self.wait_for(ready.exists, "the run")
        marked, marked_ready, marked_release = self.wait_cmd("default")
        trapped, done = self.work / "bash-got-term", self.work / "bash-done"
        script = (f"trap 'touch {shlex.quote(str(trapped))}' TERM; "
                  + " ".join(shlex.quote(c) for c in marked) + f" GC_AUTO_RUN_ID={stem}; "
                  + f"touch {shlex.quote(str(done))}")
        bash = subprocess.Popen(["/bin/bash", "-c", script])  # parent: the test
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
        """A 1 s probe; the test holds .queue.lock for about 2 s from the end of the
        probe; optionally SIGTERM the wrapper while it waits for the lock."""
        self.prepare(BETA)
        started = self.work / f"probe-started-{self.counter}"
        probe = f"touch {shlex.quote(str(started))}; sleep 1; echo {shlex.quote(str(self.src))}"
        marker, cmd = self.marker_cmd(f"n1-cmd-{self.counter}")
        proc = self.popen_run(BETA, cmd, probe=probe)
        output = self.work / f"popen-{len(self.procs) - 1}.out"
        self.wait_for(started.exists, "the probe")
        with open(self.queue / ".queue.lock", "a+") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                time.sleep(1.5)  # the probe has ended; the wrapper waits for the lock
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
        survivors (report running, run refused 13, stop signals them)."""
        self.prepare(BETA)
        self.ok(self.run_cmd(BETA, "unit"))
        self.end_worker()
        (stem, record), = self.runs().items()
        launcher, pids = self.work / "launch.py", self.work / "pids-tracked"
        launcher.write_text(LAUNCH_SCRIPT)
        cmd, child_ready, child_release = self.wait_cmd("note")
        parent = subprocess.Popen([sys.executable, str(launcher), "hold", str(self.wait_script),
                                   str(pids), f"note,{child_ready},{child_release}"],
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
        """Revision 7: cleanup ends the setsid shell of the easy_run pattern."""
        self.prepare(BETA)
        shell = self.easy_run(f"{shlex.quote(sys.executable)} -c 'import time; time.sleep(4)';"
                              " /bin/sleep 41")
        self.assertIn("running", self.result_cell(self.report(), BETA).lower())
        result = self.cleanup(0)
        self.assertTrue(result.stdout.strip().startswith("Cleaned"), result.stdout)
        self.assertFalse(alive(shell))
        self.assertNotIn("running", self.result_cell(self.report(), BETA).lower())
        self.attempt(BETA, 0)

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
        probe = "sleep 3; echo " + shlex.quote(str(self.src))
        proc = self.popen_run(BETA, cmd, probe=probe)
        stem, _ = self.wait_record(set(), status="starting")
        time.sleep(1.0)  # cleanup at about 1 s into the 3 s probe
        self.refused(self.tool("cleanup", BETA), 13)
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
        """Re-read 8: a probe longer than the limit: exit 124 within about 10 s,
        `not started` with reason `probe timed out`, CMD never runs, no probe
        process left. Control: a short probe runs CMD."""
        self.use_new_queue("q-probe-timeout", extra=("--default-timeout", "3"))
        self.prepare(BETA)
        marker, cmd = self.marker_cmd("ran-after-probe-timeout")
        began = time.time()
        result = self.run_cmd(BETA, "unit", cmd=cmd,
                              probe="sleep 30; echo " + shlex.quote(str(self.src)))
        self.assertEqual(result.returncode, 124, result.stdout + result.stderr)
        self.assertLess(time.time() - began, 10)
        self.assertFalse(marker.exists())
        (record,) = self.runs().values()
        self.assertEqual(record["status"], "not started")
        self.assertIn("probe timed out", record.get("reason", ""))
        table = subprocess.run(["/bin/ps", "-A", "-o", "pgid=,command="], text=True,
                               capture_output=True).stdout
        leftovers = [line for line in table.splitlines()
                     if line.split(None, 1)[0] == str(record["probe_pgid"])]
        self.assertEqual(leftovers, [], "a probe process was left")
        self.attempt(BETA, 0)

    def test_probe_timeout_with_detached_child_holding_output(self):
        """Probe time-out (confirmation of re-read 7): a probe whose detached python
        (own session, output inherited) outlives it does not hold the wrapper: exit
        124 within about 15 s, `not started` / `probe timed out`, the detached python
        is ended (by the run, or by cleanup before lock release is accepted)."""
        self.prepare(BETA)
        script, pidfile = self.work / "detach_probe.py", self.work / "detached.pid"
        script.write_text(DETACH_PROBE)
        marker, cmd = self.marker_cmd("ran-after-detached-probe")
        probe = f"{shlex.quote(sys.executable)} {shlex.quote(str(script))} {shlex.quote(str(pidfile))}"
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
        and mark it ready again"; another failure (a stale run) keeps "rerun the
        affected checks"."""
        freeze_step = "freeze a corrected packet in a new directory and mark it ready again"
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
        (an error or failure naming the interception, not a skip)."""
        expected = {"test_network_allowed_job": (OUTSIDE_URL, f"{OUTSIDE_ADDRESS}:443"),
                    "test_outside_name_lookup_allowed_job": (OUTSIDE_NAME,),
                    "test_outside_name_lookup_blocked": (OUTSIDE_NAME,)}

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
                if value == "1":
                    self.assertEqual(result.skipped, [])
                    reached = {test._testMethodName: trace
                               for test, trace in result.errors + result.failures}
                    self.assertEqual(sorted(reached), sorted(expected))
                    for name, trace in reached.items():
                        self.assertIn(marker, trace, name)
                    continue
                self.assertEqual((result.errors, result.failures), ([], []))
                self.assertEqual(result.testsRun, len(expected))
                reasons = {test._testMethodName: reason for test, reason in result.skipped}
                self.assertEqual(sorted(reasons), sorted(expected))
                for name, destinations in expected.items():
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


if __name__ == "__main__":
    unittest.main()
