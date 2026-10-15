"""GuidedCoding auto queue: one controller, one worker, one queue record.

State, stop flag, tracked command runs, the local run lock, reconciliation
and the morning report for `/gc auto`. Standard library only. No scheduler
and no worker pool. Interface: AUTO_QUEUE_SPEC, including Revisions 2-8
(record 2026-10-08-gc-auto-release-plan/stage1).
Revision 12: `run` refuses on any system other than macOS; `init` needs --lock-dir.
Revision 13: an unreadable process table refuses (exit 17) instead of reading as empty;
`transfer` hands jobs to an attended record (final state Transferred); `install check` can
classify timestamp-only differences from a NORMALIZED value (report-only, never available).
Revision 14: a zombie (an exited process its parent has not collected, as in a container
whose process 1 never reaps) counts as gone everywhere liveness is decided.

Stop from any Terminal, without the coding conversation:
    python3 <procedure root>/auto_queue.py stop
"""

import argparse
import calendar
import contextlib
import fcntl
import hashlib
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path

STATES = ("Queued", "Preparing", "Testing", "Ready for approval", "Approved",
          "Waiting", "Blocked", "Discarded", "Published", "Transferred")
TRANSITIONS = {
    "Queued": ("Preparing", "Blocked", "Discarded"),
    "Preparing": ("Testing", "Blocked", "Waiting", "Discarded"),
    "Testing": ("Preparing", "Ready for approval", "Blocked", "Waiting"),
    "Waiting": ("Preparing", "Testing", "Blocked"),
    "Blocked": ("Preparing", "Discarded"),
    "Ready for approval": ("Approved", "Preparing", "Discarded"),
    "Approved": ("Published", "Preparing", "Discarded"),
    "Published": (),
    "Discarded": (),
    # Revision 13: reached only through `transfer` (auto -> attended); final, never reopened:
    # the successor record owns the work from then on.
    "Transferred": (),
}
ACTIVE_STATES = ("Queued", "Preparing", "Testing", "Waiting")
FINAL_STATES = ("Published", "Discarded", "Transferred")
# Revision 12: `run` needs macOS (sandbox-exec for the network block, `ps -E` for the run-id
# marker); elsewhere it refuses with this exit code. Nothing but sys.platform decides it.
PLATFORM_REFUSED = 15
# Revision 13: the process table could not be read (ps failed, printed nothing usable, or
# left out this process): every command that needs it refuses with this code.
PROCESS_TABLE_UNREADABLE = 17
# Start times are read the same way whatever the caller's TZ, locale or PATH.
PS = "/bin/ps"
PS_ENV = {"LC_ALL": "C", "TZ": "UTC", "PATH": "/bin:/usr/bin"}
UNREADABLE = "unreadable"
SCREEN_CHECK = Path(__file__).resolve().parent / "payload" / "tools" / "screen_check.py"
MARKER = "GC_AUTO_RUN_ID"
SANDBOX = "/usr/bin/sandbox-exec"
# Outbound network denied except to this machine, and name lookups through the system resolver
# denied (localhost still resolves); everything else allowed.
NO_NETWORK = ('(version 1)(allow default)(deny network-outbound)'
              '(allow network-outbound (remote ip "localhost:*"))'
              '(allow network-outbound (remote unix-socket))'
              '(deny network-outbound (remote unix-socket (path-literal "/private/var/run/mDNSResponder")))')


class Refused(Exception):
    def __init__(self, code, message):
        Exception.__init__(self, message)
        self.code = code


def platform_supported():
    """Revision 12: auto-mode runs need macOS. Reads sys.platform only; there is no switch."""
    return sys.platform == "darwin"


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def tmp_name(path):
    return path.with_name(".%s.tmp.%d.%d" % (path.name, os.getpid(), time.monotonic_ns()))


def write_json(path, data):
    """Atomic write: unique temp file in the same directory, fsync, replace."""
    path = Path(path)
    tmp = tmp_name(path)
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, sort_keys=True)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def write_text(path, text):
    path = Path(path)
    tmp = tmp_name(path)
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------- processes

def proc_reading(pid):
    """Linux: (state, number of threads) of PID from /proc/<pid>/stat, or None if unreadable.
    The command name in parentheses may contain spaces or ')': the fields follow the last ')'
    (state is field 3, num_threads field 20)."""
    try:
        with open("/proc/%d/stat" % pid, "rb") as handle:
            data = handle.read()
        fields = data[data.rindex(b")") + 1:].split()
        return fields[0].decode("ascii", "replace"), int(fields[17])
    except (OSError, ValueError, IndexError):
        return None


def ps_reading(pid):
    """(ps stat, number of threads or None) of PID; ("", None) when ps cannot read it. The
    thread count (nlwp) is asked for on Linux only; macOS ps has no such column."""
    linux = sys.platform.startswith("linux")
    try:
        out = subprocess.run([PS, "-o", "stat="] + (["-o", "nlwp="] if linux else [])
                             + ["-p", str(pid)], env=PS_ENV, capture_output=True, text=True)
    except OSError:
        return "", None
    fields = out.stdout.split() if out.returncode == 0 else []
    if not fields:
        return "", None
    threads = int(fields[1]) if linux and len(fields) > 1 and fields[1].isdigit() else None
    return fields[0], threads


def zombie_gone(pid, state=None):
    """Revision 14: True only for a positive reading that PID is a zombie: a process that has
    exited and only waits for its parent to collect it (in a container whose process 1 never
    does, forever). Such a process runs nothing and counts as gone wherever liveness is
    decided. Everything else counts as alive, including any state that cannot be read.
    A state beginning with Z is not enough on Linux: there the leader of a multi-threaded
    process whose main thread has exited also shows Z (ps `Zl`) while its other threads run.
    So on Linux the process must also be single-threaded, read from /proc/<pid>/stat (state Z,
    num_threads 1) or else from ps (stat Z without the `l` flag, nlwp 1); an unknown thread
    count counts as alive. On macOS (no thread count in ps, and no `l` flag) a stat of Z
    without `l` is enough: there a process with a running thread never shows Z. STATE is the
    ps stat when the caller has already read it (macOS uses it as is)."""
    if sys.platform == "darwin":
        if state is None:
            state = ps_reading(pid)[0]
        return state.startswith("Z") and "l" not in state[1:]
    if state is not None and (not state.startswith("Z") or "l" in state[1:]):
        return False
    reading = proc_reading(pid) if sys.platform.startswith("linux") else None
    if reading is None:
        reading = ps_reading(pid)
    state, threads = reading
    return state.startswith("Z") and "l" not in state[1:] and threads == 1


def group_states(pgid):
    """(pid, ps stat) of the members of process group PGID from `ps -ax`, or None when the
    table cannot be read."""
    try:
        out = read_ps(["-ax", "-o", "pid=", "-o", "pgid=", "-o", "stat="])
    except ProcessTableUnreadable:
        return None
    states = []
    for line in out.splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[0].isdigit() and fields[1].isdigit() \
                and int(fields[1]) == pgid:
            states.append((int(fields[0]), fields[2]))
    return states


def process_exists(pid):
    """True while PID is a process that has not exited. Revision 14: a zombie is gone."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    return not zombie_gone(pid)


def process_start(pid):
    """Start time of a live process, None if it is gone, UNREADABLE if it exists unread."""
    if not pid or not process_exists(pid):
        return None
    try:
        out = subprocess.run([PS, "-o", "lstart=", "-p", str(pid)], env=PS_ENV,
                             capture_output=True, text=True)
        text = out.stdout.strip()
    except OSError:
        text = ""
    # Revision 13: only a start time that parses counts; anything else is unreadable.
    if text and epoch(text, "%a %b %d %H:%M:%S %Y") is not None:
        return text
    return UNREADABLE if process_exists(pid) else None


def same_process(pid, recorded_start):
    """True while the recorded process is alive; an unreadable start time counts as alive."""
    start = process_start(pid)
    if start is None:
        return False
    if start == UNREADABLE or recorded_start in (None, UNREADABLE):
        return True
    return start == recorded_start


def group_alive(pgid):
    """True while process group PGID has a member that has not exited. Revision 14: a group
    whose only members are zombies is gone (macOS then answers killpg with EPERM, Linux with
    success). An unreadable table keeps the group alive."""
    if not pgid:
        return False
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    states = group_states(pgid)
    if states is None:
        return True
    if not states:
        # Every member left between the two reads, or ps does not list them: ask again.
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True
    return not all(zombie_gone(pid, state) for pid, state in states)


_SNAPSHOT = {"time": 0.0, "procs": None}


def epoch(text, fmt):
    try:
        return calendar.timegm(time.strptime(text.strip(), fmt))
    except (ValueError, TypeError):
        return None


class ProcessTableUnreadable(Refused):
    """Revision 13: the process table could not be read, so whether a process is alive is
    unknown. Every caller refuses or reports unknown; none treats it as an empty table."""

    def __init__(self, why):
        Refused.__init__(self, PROCESS_TABLE_UNREADABLE,
                         "cannot read the process table (%s); whether queue processes are alive "
                         "is unknown" % why)


def read_ps(args):
    """stdout of `ps ARGS`; raises ProcessTableUnreadable if ps cannot start or exits nonzero."""
    try:
        out = subprocess.run([PS] + args, env=PS_ENV, capture_output=True, text=True)
    except OSError as error:
        raise ProcessTableUnreadable("%s could not be started: %s" % (PS, error))
    if out.returncode != 0:
        detail = out.stderr.strip().splitlines()[0] if out.stderr.strip() else ""
        raise ProcessTableUnreadable("%s exited %d%s" % (PS, out.returncode,
                                                         (": " + detail) if detail else ""))
    return out.stdout


def snapshot():
    """All processes: pid -> [ppid, pgid, sid, start epoch, run ids in the environment, zombie].
    Revision 14: a zombie (zombie_gone) stays in the table with its flag set, so that a
    recorded leader id is still compared with its start time (leaders_of); it counts as gone
    for liveness (live_table, survivors). A Z line whose start time does not parse is kept
    with start None.
    Revision 13: fails closed. A ps that cannot start, exits nonzero, prints a line that does
    not parse, or prints a table without this process raises ProcessTableUnreadable; nothing
    is cached then. The `ps -E` marker query is held to the same rule on macOS; elsewhere
    (Linux ps has no -E) its failure only drops the markers, as before. Checkpoint 3: on
    macOS the marker output must also contain a parsable line for this process (empty or
    unparsable output with exit 0 refuses, as for the main table); continuation lines of
    multi-line environment values are still skipped."""
    if _SNAPSHOT["procs"] is not None and time.time() - _SNAPSHOT["time"] < 0.3:
        return _SNAPSHOT["procs"]
    procs = {}
    out = read_ps(["-ax", "-o", "pid=", "-o", "ppid=", "-o", "pgid=", "-o", "stat=",
                   "-o", "lstart="])
    marker_args = ["-E", "-ww", "-ax", "-o", "pid=", "-o", "command="]
    if sys.platform == "darwin":
        marks = read_ps(marker_args)
    else:
        try:
            marks = read_ps(marker_args)
        except ProcessTableUnreadable:
            marks = ""
    for number, line in enumerate(out.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split(None, 4)
        ids = len(fields) >= 4 and all(f.isdigit() for f in fields[:3])
        zombie = ids and fields[3].startswith("Z") and zombie_gone(int(fields[0]), fields[3])
        start = epoch(fields[4], "%a %b %d %H:%M:%S %Y") if len(fields) == 5 else None
        if not ids or (start is None and not zombie):
            raise ProcessTableUnreadable("unparsable ps line %d: %r" % (number, line[:80]))
        pid, ppid, pgid = int(fields[0]), int(fields[1]), int(fields[2])
        try:
            sid = os.getsid(pid)
        except OSError:
            sid = None
        procs[pid] = [ppid, pgid, sid, start, set(), bool(zombie)]
    if os.getpid() not in procs:
        raise ProcessTableUnreadable("ps listed %d process(es) but not this one (%d)"
                                     % (len(procs), os.getpid()))
    own = False
    for line in marks.splitlines():
        # An environment value may contain a newline: such continuation lines are skipped.
        fields = line.split(None, 1)
        if len(fields) == 2 and fields[0].isdigit() and int(fields[0]) in procs:
            procs[int(fields[0])][4].update(re.findall(r"\b%s=(\S+)" % MARKER, fields[1]))
            own = own or int(fields[0]) == os.getpid()
    if sys.platform == "darwin" and not own:
        raise ProcessTableUnreadable("ps -E listed no environment line for this process (%d)"
                                     % os.getpid())
    _SNAPSHOT.update(time=time.time(), procs=procs)
    return procs


def live_table(procs):
    """Revision 14: the table without zombies (None stays None), for every liveness decision."""
    if procs is None:
        return None
    return {pid: entry for pid, entry in procs.items() if not (len(entry) > 5 and entry[5])}


def snapshot_or_none():
    """For a run's own bookkeeping while its command runs: the table, or None when it is
    unreadable (the caller then signals by group only and records the gap; it never
    concludes that nothing is left)."""
    try:
        return snapshot()
    except ProcessTableUnreadable:
        return None


def leaders_of(record, procs):
    """Group/session ids of the run's probe and CMD, unless the id may have been reused."""
    leaders = set()
    for key, start_key in (("pgid", "process_start"), ("probe_pgid", "probe_start")):
        leader = record.get(key)
        if not leader:
            continue
        recorded = epoch(record.get(start_key) or "", "%a %b %d %H:%M:%S %Y")
        live = procs.get(leader)
        # A leader with another start time means the id was reused: not ours. Revision 14:
        # this includes a zombie (still in the table) that holds the id.
        if live is not None and (recorded is None or live[3] is None or live[3] != recorded):
            continue
        leaders.add(leader)
    return leaders


def survivors(record):
    """Live processes of this run (spec revisions 4 and 5): descendants the wrapper recorded;
    processes in the run's group or session, or carrying its marker, started while it ran and
    descended from it, launchd or the wrapper; parents of marked processes started while it
    ran; and every later descendant of any of these."""
    procs = snapshot()
    me, wrapper = os.getpid(), record.get("wrapper_pid")
    started = epoch(record.get("started", ""), "%Y-%m-%dT%H:%M:%SZ")
    ended = epoch(record.get("ended", ""), "%Y-%m-%dT%H:%M:%SZ")

    def during(start):
        # Both times are whole UTC seconds: a process that started while the run ran cannot
        # have a later second than `ended`. Revision 12 (P10): the earlier bound allows one
        # second (`started` - 1), as the descendant rule below does, because Linux ps can
        # report a start time up to 1 s early; every other identity check still applies.
        return (start is not None and (started is None or start >= started - 1)
                and (ended is None or start <= ended))

    probe_start = epoch(record.get("probe_started", ""), "%Y-%m-%dT%H:%M:%SZ")
    probe_end = epoch(record.get("probe_ended", ""), "%Y-%m-%dT%H:%M:%SZ")

    def during_probe(start):
        # The same one-second allowance on the earlier bound (Revision 12, P10).
        return (start is not None and (probe_start is None or start >= probe_start - 1)
                and (probe_end is None or start <= probe_end))

    leaders = leaders_of(record, procs)  # Revision 14: with zombies, for the start-time check
    procs = live_table(procs)  # then only processes that have not exited
    found = set()
    for item in record.get("tracked", []):
        live = procs.get(item.get("pid"))
        if live is not None and live[3] is not None and live[3] == item.get("start"):
            found.add(item["pid"])
    probe_leader = record.get("probe_pgid") if record.get("probe_pgid") in leaders else None
    roots = {1, wrapper} | leaders
    marked = set()
    for pid, (ppid, pgid, sid, start, marks) in ((p, e[:5]) for p, e in procs.items()):
        if pid in (me, wrapper) or not during(start):
            continue
        if record["run_id"] in marks:
            marked.add(pid)
            found.add(pid)
            continue
        cmd_match = (pgid in leaders - {probe_leader} or sid in leaders - {probe_leader})
        probe_match = probe_leader is not None and (pgid == probe_leader or sid == probe_leader) \
            and during_probe(start)
        if (cmd_match or probe_match) and (pid in leaders or ppid in roots):
            found.add(pid)
    # A system shell running a marked program has no visible marker itself: count the chain of
    # parents of a marked process only when the chain connects to the run (its top is launchd,
    # the wrapper, a leader or an already counted process), never a user's own shell.
    for pid in marked:
        chain, parent = [], procs[pid][0]
        while parent in procs and parent not in (0, 1, me, wrapper) and parent not in found \
                and during(procs[parent][3]):
            chain.append(parent)
            parent = procs[parent][0]
        if chain and (parent in (1, wrapper) or parent in leaders or parent in found):
            found.update(chain)
    grew = True
    while grew:
        more = {pid for pid, info in procs.items()
                if pid not in found and pid not in (me, wrapper) and info[0] in found
                and (started is None or info[3] is None or info[3] >= started - 1)}
        found |= more
        grew = bool(more)
    return sorted(found)


def survivors_or_unknown(record):
    """survivors(record), or None when the process table is unreadable (Revision 13)."""
    try:
        return survivors(record)
    except ProcessTableUnreadable:
        return None


UNKNOWN_LEFT = "unknown (process table unreadable)"


def run_live(record):
    """Why a run record still has live processes of ours, or None."""
    status = record.get("status")
    if status == "starting" and same_process(record.get("wrapper_pid"), record.get("wrapper_start")):
        return "run %s starting" % record["run_id"]
    if status == "running" and same_process(record.get("pid"), record.get("process_start")):
        return "run %s still alive" % record["run_id"]
    left = survivors(record)
    if left:
        return "run %s has %d live process(es)" % (record["run_id"], len(left))
    return None


def verify_packet(packet):
    """(manifest identity, None) if the frozen packet verifies completely, else (None, reason)."""
    packet = Path(packet)
    env = dict(os.environ, GC_PAYLOAD_ROOT=str(SCREEN_CHECK.parent.parent))
    out = subprocess.run([sys.executable, "-I", "-B", str(SCREEN_CHECK), "verify", str(packet)],
                         env=env, capture_output=True, text=True)
    if out.returncode != 0:
        return None, "packet does not verify: %s" % (out.stderr.strip() or out.stdout.strip())
    return hashlib.sha256((packet / "MANIFEST.sha256").read_bytes()).hexdigest(), None


# ---------------------------------------------------------------- queue

class Queue:
    def __init__(self, directory):
        self.dir = Path(directory).resolve()
        self.file = self.dir / "QUEUE.json"
        if not self.file.is_file():
            raise Refused(2, "no queue record in %s" % self.dir)
        self.data = read_json(self.file)

    def reload(self):
        self.data = read_json(self.file)

    @contextlib.contextmanager
    def locked(self):
        """Exclusive lock for one read-modify-write of this queue's state."""
        with open(self.dir / ".queue.lock", "a") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                self.reload()
                yield self
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def save(self):
        write_json(self.file, self.data)

    def event(self, text):
        with open(self.dir / "EVENTS.log", "a", encoding="utf-8") as f:
            f.write("%s %s\n" % (now(), text))
            f.flush()
            os.fsync(f.fileno())

    def job(self, job_id):
        for job in self.data["jobs"]:
            if job["id"] == job_id:
                return job
        raise Refused(2, "no job %s" % job_id)

    def stop_requested(self):
        return (self.dir / "STOP_REQUESTED").exists()

    def runs(self, job_id=None):
        result = []
        run_dir = self.dir / "runs"
        if run_dir.is_dir():
            for path in sorted(run_dir.glob("*.json")):
                record = read_json(path)
                if job_id is None or record.get("job") == job_id:
                    result.append((path, record))
        return result

    def mark_stale(self, job_id, why):
        for path, record in self.runs(job_id):
            if not record.get("stale"):
                record["stale"] = True
                record["stale_reason"] = why
                write_json(path, record)

    def worker_active(self):
        worker = self.data.get("worker")
        return worker if worker and not worker.get("ended") else None

    def live_runs(self, job_id=None):
        return [(p, r, why) for p, r in self.runs(job_id) for why in [run_live(r)] if why]

    def remaining(self):
        """What keeps a stop from being complete."""
        left = []
        worker = self.worker_active()
        if worker:
            left.append("worker for job %s registered, not yet ended" % worker["job"])
        left.extend(why for _, _, why in self.live_runs())
        return left

    def finished(self):
        if any(j["state"] in ACTIVE_STATES for j in self.data["jobs"]):
            return False
        return not self.worker_active() and not (self.stop_requested() and self.remaining())

    def lock_dir(self):
        return Path(self.data["lock_dir"])

    def lock_owner(self):
        """Owner fields, None if there is no owner file; raises OSError if unreadable."""
        owner = self.lock_dir() / "owner.txt"
        if not owner.exists():
            return None
        fields = {}
        for line in owner.read_text(encoding="utf-8", errors="replace").splitlines():
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
        return fields

    def lock_held(self):
        try:
            owner = self.lock_owner()
        except OSError:
            return False
        return (owner is not None and owner.get("queue_id") == self.data["queue_id"]
                and owner.get("queue_dir") == str(self.dir))

    def readiness(self, job, packet=None):
        """None if the job may be Ready for approval (or Approved), else the failed condition."""
        if not job.get("criterion_sha256"):
            return "no criterion recorded"
        if not job.get("candidate_id"):
            return "no candidate recorded"
        worker = self.worker_active()
        if worker and worker["job"] == job["id"]:
            return "worker %s for this job is still registered" % worker["id"]
        live = self.live_runs(job["id"])
        if live:
            return live[0][2]
        roots = [os.path.realpath(r) for r in self.data.get("allowed_roots", [])]
        runs = [r for _, r in self.runs(job["id"])]
        for label in job.get("requires", []):
            good = [r for r in runs
                    if r.get("label") == label and r.get("status") == "completed"
                    and r.get("exit_code") == 0 and not r.get("stale")
                    and not r.get("group_left_running")
                    and not r.get("process_table_unreadable")
                    and r.get("candidate_id") == job["candidate_id"]
                    and r.get("criterion_sha256") == job["criterion_sha256"]
                    and r.get("tested_paths")
                    and all(under(p, roots) for p in r["tested_paths"])]
            if not good:
                return "no current passing run for required check %s" % label
        if packet is None:
            packet = job.get("packet")
        if not packet:
            return "no frozen packet given"
        packet = Path(packet)
        identity_sha, failed = verify_packet(packet)
        if failed:
            return failed
        stored = job.get("packet_manifest_sha256")
        if stored and job.get("packet") == str(packet.resolve()) and stored != identity_sha:
            return "packet changed since it was marked ready"
        trees = {r["tree"] for r in job.get("candidate", [])}
        identity = packet / "CODE_IDENTITY.txt"
        named, ids = set(), set()
        if identity.is_file():
            for line in identity.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("tested_tree:"):
                    named.add(line.split(":", 1)[1].strip())
                elif line.startswith("candidate_id:"):
                    ids.add(line.split(":", 1)[1].strip())
        if ids != {job["candidate_id"]}:
            return "packet CODE_IDENTITY.txt does not name the candidate id"
        if not trees or not trees <= named:
            return "packet CODE_IDENTITY.txt does not name the tested tree of every candidate repository"
        return screening_problem(packet / "SCREENING.txt", job)


SCREEN_LINE = re.compile(r"test: (\S.*?) \| outside: (\S.*?) \| (included|excluded)")


def screening_problem(path, job):
    """None if the packet's outside-access screening is recorded and fits the job's grant."""
    if not path.is_file():
        return "packet has no SCREENING.txt"
    included = False
    for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        match = SCREEN_LINE.fullmatch(line)
        if not match:
            return "SCREENING.txt line %d malformed" % number
        name, outside, decision = match.groups()
        if decision == "included":
            included = True
            if outside != "none" and not job.get("network"):
                return ("SCREENING.txt includes %s with outside access but the job has no "
                        "network grant" % name)
    if not included:
        return "SCREENING.txt lists no included test"
    return None


def under(path, roots):
    path = os.path.realpath(path)
    return any(path == r or path.startswith(r.rstrip(os.sep) + os.sep) for r in roots)


def names_job(quote, job_id):
    return re.search(r"(?<![A-Za-z0-9_-])%s(?![A-Za-z0-9_-])" % re.escape(job_id), quote) is not None


def home_dir(args):
    if args.home:
        return Path(args.home)
    if os.environ.get("GC_AUTO_HOME"):
        return Path(os.environ["GC_AUTO_HOME"])
    return Path.home() / "GuidedCoding" / "auto"


def active_dir(args):
    active = home_dir(args) / "ACTIVE"
    if active.is_file():
        lines = active.read_text(encoding="utf-8").splitlines()
        if lines and lines[0].strip() and Path(lines[0].strip(), "QUEUE.json").is_file():
            return lines[0].strip()
    return None


def open_queue(args):
    directory = args.queue or active_dir(args)
    if not directory:
        raise Refused(2, "no active queue (give --queue DIR)")
    return Queue(directory)


def free_gib(path):
    return shutil.disk_usage(path).free / float(1 << 30)


# ---------------------------------------------------------------- commands

def cmd_init(args):
    if not args.lock_dir:
        # Revision 12: no default lock; nothing is created or changed without one.
        raise Refused(2, "init needs --lock-dir: give the lock directory from the project "
                         "method's lock directory")
    previous = active_dir(args)
    if previous and not Queue(previous).finished():
        raise Refused(8, "active queue %s is not finished" % previous)
    if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", args.queue_id)
            or args.queue_id in (".", "..")):
        raise Refused(2, "queue id must be one directory name")
    records = Path(args.records).resolve()
    directory = records / args.queue_id
    if directory.resolve().parent != records:
        raise Refused(2, "queue id must be one directory name inside the records root")
    if directory.is_symlink() or directory.exists():
        raise Refused(2, "queue directory %s exists" % directory)
    if args.space_path and not Path(args.space_path).exists():
        raise Refused(2, "space path %s does not exist" % args.space_path)
    jobs = read_json(args.jobs)
    seen = set()
    for job in jobs:
        bad = "job id %r must not start with '-' or contain spaces or '|'"
        if isinstance(job, dict) and "id" in job and not isinstance(job["id"], str):
            raise Refused(2, bad % (job["id"],))
        if not isinstance(job, dict) or not job.get("id") or job["id"] in seen:
            raise Refused(2, "jobs file needs unique ids")
        if job["id"].startswith("-") or re.search(r"[\s|]", job["id"]):
            raise Refused(2, bad % (job["id"],))
        seen.add(job["id"])
        if not isinstance(job.get("network", False), bool):
            raise Refused(2, "job %s: network must be true or false" % job["id"])
    directory.mkdir(parents=True)
    (directory / "runs").mkdir()
    shutil.copyfile(args.grant, directory / "GRANT.md")
    shutil.copyfile(args.jobs, directory / "JOBS.json")
    data = {
        "queue_id": args.queue_id,
        "created": now(),
        "lock_dir": str(Path(args.lock_dir).resolve()),
        "min_free_gib": args.min_free_gib,
        "space_path": str(Path(args.space_path).resolve()) if args.space_path else str(directory),
        "allowed_roots": [str(Path(r).resolve()) for r in (args.allowed_root or [])],
        "default_timeout": args.default_timeout,
        "worker": None,
        "jobs": [{"id": j["id"], "title": j.get("title", ""),
                  "requires": list(j.get("requires", [])),
                  "depends_on": list(j.get("depends_on", [])),
                  "network": bool(j.get("network", False)),
                  "state": "Queued", "notes": []} for j in jobs],
    }
    write_json(directory / "QUEUE.json", data)
    queue = Queue(directory)
    queue.event("init queue %s jobs %s" % (args.queue_id, ",".join(sorted(seen))))
    home = home_dir(args)
    home.mkdir(parents=True, exist_ok=True)
    write_text(home / "ACTIVE", str(queue.dir) + "\n")
    print(queue.dir)


def cmd_status(args):
    queue = open_queue(args)
    data = dict(queue.data)
    data["stop_requested"] = queue.stop_requested()
    data["runs"] = [r for _, r in queue.runs()]
    print(json.dumps(data, indent=1, sort_keys=True))


def cmd_set(args):
    queue = open_queue(args)
    if args.state not in STATES:
        raise Refused(2, "unknown state %r" % args.state)
    if args.state == "Transferred":
        raise Refused(2, "a job is transferred only with `transfer JOB --to DIR --note TEXT`")
    with queue.locked():
        job = queue.job(args.job)
        if queue.stop_requested() and args.state not in ("Blocked", "Waiting"):
            raise Refused(3, "stop requested")
        if args.state not in TRANSITIONS[job["state"]]:
            raise Refused(8, "illegal transition %s -> %s" % (job["state"], args.state))
        if args.state == "Ready for approval":
            failed = queue.readiness(job, args.packet)
            if failed:
                raise Refused(8, "not ready: %s" % failed)
            job["packet"] = str(Path(args.packet).resolve())
            job["packet_manifest_sha256"] = verify_packet(job["packet"])[0]
        if args.state == "Published" and any(
                p["state"] == "pending" for p in job.get("provisional", [])):
            raise Refused(8, "job %s has an undecided provisional choice" % job["id"])
        if args.state == "Approved":
            if not args.quote or not names_job(args.quote, job["id"]):
                raise Refused(8, "approval needs the Developer's quote naming %s" % job["id"])
            if not args.manifest or args.manifest != job.get("packet_manifest_sha256"):
                raise Refused(8, "approval must name the manifest identity of the packet marked ready")
            pending = [p for p in job.get("provisional", []) if p["state"] == "pending"]
            if pending:
                raise Refused(8, "provisional choice P%d of job %s needs its own decision first"
                              % (pending[0]["number"], job["id"]))
            failed = queue.readiness(job)
            if failed:
                raise Refused(8, "no longer ready: %s" % failed)
            job["approval_quote"] = args.quote
            job["approval_manifest_sha256"] = args.manifest
            queue.event("job %s approval quote: %s" % (job["id"], args.quote.replace("\n", " ")))
        old = job["state"]
        if old in ("Ready for approval", "Approved") and args.state in ("Preparing", "Discarded"):
            job.setdefault("binding_history", []).append(
                {"time": now(), "left": old, "to": args.state,
                 **{key: job.get(key) for key in ("packet", "packet_manifest_sha256",
                                                  "approval_quote", "approval_manifest_sha256")
                    if job.get(key)}})
            for key in ("packet", "packet_manifest_sha256", "approval_quote", "approval_manifest_sha256"):
                job.pop(key, None)
        job["state"] = args.state
        if args.note:
            job["notes"].append({"time": now(), "state": args.state, "note": args.note})
        if args.state == "Blocked":
            job["decision_needed"] = args.decision or ""
        elif args.decision:
            raise Refused(2, "--decision applies only to Blocked")
        queue.save()
        queue.event("job %s %s -> %s%s" % (job["id"], old, args.state,
                                           (" note: " + args.note) if args.note else ""))
    print("%s %s" % (job["id"], args.state))


def changeable(queue, job):
    if job["state"] in FINAL_STATES:
        raise Refused(8, "job %s is %s" % (job["id"], job["state"]))
    if queue.stop_requested():
        raise Refused(3, "stop requested")


def cmd_criterion(args):
    queue = open_queue(args)
    text = Path(args.file).read_text(encoding="utf-8")
    digest = sha256_text(text)
    with queue.locked():
        job = queue.job(args.job)
        changeable(queue, job)
        old = job.get("criterion_sha256")
        if old == digest:
            print("criterion unchanged")
            return
        if old:
            queue.mark_stale(job["id"], "criterion changed")
            job.setdefault("criterion_history", []).append(
                {"time": now(), "old": old, "old_text": job.get("criterion", ""), "new": digest})
        job["criterion"] = text
        job["criterion_sha256"] = digest
        queue.save()
        queue.event("job %s criterion %s -> %s" % (job["id"], old or "-", digest))
    print(digest)


def cmd_candidate(args):
    queue = open_queue(args)
    if not args.repo or len(args.repo) != len(args.commit or []):
        raise Refused(2, "give --repo and --commit in pairs")
    parts = []
    for repo, rev in zip(args.repo, args.commit):
        commit = git_rev(repo, rev + "^{commit}")
        tree = git_rev(repo, rev + "^{tree}")
        parts.append({"repo": str(Path(repo).resolve()), "commit": commit, "tree": tree})
    parts.sort(key=lambda p: p["repo"])
    digest = sha256_text("\n".join("%s:%s" % (p["repo"], p["commit"]) for p in parts))
    with queue.locked():
        job = queue.job(args.job)
        changeable(queue, job)
        old = job.get("candidate_id")
        if old == digest:
            print("candidate unchanged")
            return
        if old:
            queue.mark_stale(job["id"], "candidate changed")
        job["candidate"] = parts
        job["candidate_id"] = digest
        queue.save()
        queue.event("job %s candidate %s -> %s" % (job["id"], old or "-", digest))
    print(digest)


def git_rev(repo, rev):
    out = subprocess.run(["git", "-C", repo, "rev-parse", "--verify", "--quiet", rev],
                         capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        raise Refused(2, "cannot resolve %s in %s" % (rev, repo))
    return out.stdout.strip()


def cmd_worker(args):
    queue = open_queue(args)
    with queue.locked():
        job = queue.job(args.job)
        active = queue.worker_active()
        if args.action == "start":
            if queue.stop_requested():
                raise Refused(3, "stop requested")
            if job["state"] in FINAL_STATES:
                raise Refused(8, "job %s is %s" % (job["id"], job["state"]))
            if active:
                raise Refused(9, "worker %s for job %s is still registered"
                              % (active["id"], active["job"]))
            queue.data["worker"] = {"job": args.job, "id": args.id, "started": now()}
            queue.save()
            queue.event("worker start %s job %s" % (args.id, args.job))
        else:
            if not active or active["id"] != args.id or active["job"] != args.job:
                raise Refused(9, "worker %s is not the registered worker" % args.id)
            active["ended"] = now()
            active["outcome"] = args.outcome
            queue.save()
            queue.event("worker end %s job %s outcome %s" % (args.id, args.job, args.outcome))
    print("worker %s %s" % (args.id, args.action))


def cmd_check_stop(args):
    queue = open_queue(args)
    if queue.stop_requested():
        print("STOP REQUESTED")
        raise SystemExit(3)
    print("no stop requested")


def cmd_lock(args):
    queue = open_queue(args)
    lock = queue.lock_dir()
    if args.action == "show":
        try:
            owner = lock / "owner.txt"
            print(owner.read_text(encoding="utf-8") if owner.is_file() else
                  ("held, owner unreadable" if lock.exists() else "free"))
        except OSError:
            print("held, owner unreadable")
        return
    if args.action == "take":
        lock.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.mkdir(lock)
        except FileExistsError:
            if queue.lock_held():
                print("already held")
                return
            raise Refused(4, "lock %s held by another owner" % lock)
        write_text(lock / "owner.txt",
                   "queue_id: %s\nqueue_dir: %s\npid: %d\nhost: %s\ntime: %s\n"
                   % (queue.data["queue_id"], queue.dir, os.getpid(), socket.gethostname(), now()))
        queue.event("lock taken %s" % lock)
        print("taken")
        return
    with queue.locked():
        if not queue.lock_held():
            raise Refused(4, "lock %s is not held by this queue" % lock)
        _SNAPSHOT["procs"] = None
        snapshot()  # Checkpoint 3: unreadable refuses (17) even with no run record
        live = queue.live_runs()
        if live:
            raise Refused(13, "queue-owned processes are still alive (%s); run cleanup first"
                          % "; ".join(why for _, _, why in live))
        extra = sorted(set(os.listdir(lock)) - {"owner.txt"})
        if extra:
            raise Refused(4, "lock %s contains other files (%s); not removed"
                          % (lock, ", ".join(extra)))
        os.remove(lock / "owner.txt")
        os.rmdir(lock)
        queue.event("lock released %s" % lock)
    print("released")


def cmd_run(args):
    if not platform_supported():
        # Revision 12: before the queue is read, any record is written or anything is started.
        raise Refused(PLATFORM_REFUSED,
                      "auto-mode runs need macOS (sandbox-exec and `ps -E`); this system is %s"
                      % sys.platform)
    queue = open_queue(args)
    command = list(args.command)
    if not command:
        raise Refused(2, "no command after --")
    cwd = os.path.abspath(args.cwd or os.getcwd())
    received = []

    def note_signal(signum, frame):
        received.append(signum)
    signal.signal(signal.SIGTERM, note_signal)
    signal.signal(signal.SIGINT, note_signal)

    run_id = "%s-%s-%d" % (time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()), args.job, os.getpid())
    record_path = queue.dir / "runs" / ("%s.json" % run_id)
    with queue.locked():
        job = queue.job(args.job)
        if queue.stop_requested():
            raise Refused(3, "stop requested")
        if job["state"] not in ("Preparing", "Testing"):
            raise Refused(8, "job %s is %s, not Preparing or Testing" % (job["id"], job["state"]))
        worker = queue.worker_active()
        if not worker or worker["job"] != job["id"]:
            raise Refused(9, "no registered worker for job %s" % job["id"])
        if not queue.lock_held():
            raise Refused(4, "lock not held by this queue")
        installation = queue.data.get("installation") or {}
        if installation.get("state") == "unavailable":
            raise Refused(12, "test installation unavailable: %s" % installation.get("reason", "-"))
        free = free_gib(queue.data["space_path"])
        if free < float(queue.data["min_free_gib"]):
            raise Refused(5, "free space %.2f GiB below floor %s GiB"
                          % (free, queue.data["min_free_gib"]))
        if not job.get("criterion_sha256") or not job.get("candidate_id"):
            raise Refused(6, "record the criterion and candidate before running checks")
        # Revision 13: a run starts only while the process table can be read (it is needed to
        # find earlier runs' processes and to track this one); unreadable refuses with 17.
        _SNAPSHOT["procs"] = None
        snapshot()
        live = queue.live_runs(job["id"])
        if live:
            raise Refused(13, "an earlier run of job %s still has live processes (%s)"
                          % (job["id"], live[0][2]))
        timeout = args.timeout if args.timeout is not None else int(queue.data.get("default_timeout", 3600))
        if timeout <= 0:
            raise Refused(2, "--timeout must be positive")
        network = "allowed" if job.get("network") else "blocked"
        if network == "blocked" and not os.path.exists(SANDBOX):
            raise Refused(14, "outside network must be blocked but %s is not available" % SANDBOX)
        record = {
            "run_id": run_id, "job": args.job, "label": args.label, "status": "starting",
            "candidate_id": job["candidate_id"], "criterion_sha256": job["criterion_sha256"],
            "command": command, "cwd": cwd, "log": str(Path(args.log).resolve()),
            "tested_paths": [], "wrapper_pid": os.getpid(),
            "wrapper_start": process_start(os.getpid()), "started": now(), "stale": False,
            "timeout_seconds": timeout, "network": network,
        }
        write_json(record_path, record)
        queue.event("run %s starting job %s label %s" % (run_id, args.job, args.label))

    def finish(status, **fields):
        with queue.locked():
            current = read_json(record_path)
            current.update(fields)
            current["status"] = status
            current["ended"] = now()
            write_json(record_path, current)
            queue.event("run %s %s%s" % (run_id, status,
                                         (" exit %s" % fields["exit_code"]) if "exit_code" in fields else ""))

    marked_env = dict(os.environ, **{MARKER: run_id, "GC_AUTO_QUEUE": str(queue.dir)})
    tested = []
    if args.probe:
        with queue.locked():
            probe_cmd = ["/bin/sh", "-c", args.probe]
            if network == "blocked":
                probe_cmd = [SANDBOX, "-p", NO_NETWORK] + probe_cmd
            probe = subprocess.Popen(probe_cmd, cwd=cwd, env=marked_env,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                     stdin=subprocess.DEVNULL, start_new_session=True)
            current = read_json(record_path)
            current["probe_pgid"] = probe.pid
            current["probe_start"] = process_start(probe.pid)
            current["probe_started"] = now()
            write_json(record_path, current)
        try:
            probe_out, _ = probe.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            # End the probe's group and every process identified as the run's, then stop
            # waiting on its output: a detached child may hold the pipe open.
            for sig in (signal.SIGTERM, signal.SIGKILL):
                _SNAPSHOT["procs"] = None
                targets = set(survivors_or_unknown(read_json(record_path)) or ())
                try:
                    os.killpg(probe.pid, sig)
                except OSError:
                    pass
                for pid in targets:
                    try:
                        os.kill(pid, sig)
                    except OSError:
                        pass
                time.sleep(2 if sig == signal.SIGTERM else 0.5)
            try:
                probe.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                for stream in (probe.stdout, probe.stderr):
                    try:
                        stream.close()
                    except OSError:
                        pass
                try:
                    probe.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
            _SNAPSHOT["procs"] = None
            left = survivors_or_unknown(read_json(record_path))
            finish("not started", reason="probe timed out after %d s" % timeout,
                   left_after_timeout=UNKNOWN_LEFT if left is None else left)
            raise SystemExit(124)
        with queue.locked():
            current = read_json(record_path)
            current["probe_ended"] = now()
            write_json(record_path, current)
        if received or queue.stop_requested():
            reason = "stop requested" if queue.stop_requested() else "signal during probe"
            finish("not started", reason=reason)
            if queue.stop_requested():
                raise Refused(3, "stop requested")
            raise SystemExit(128 + received[0])
        if probe.returncode != 0:
            finish("not started", reason="probe failed with exit %d" % probe.returncode)
            raise Refused(7, "probe failed with exit %d" % probe.returncode)
        tested = [os.path.realpath(os.path.join(cwd, line.strip()))
                  for line in probe_out.splitlines() if line.strip()]
        roots = [os.path.realpath(r) for r in queue.data.get("allowed_roots", [])]
        missing = [p for p in tested if not os.path.exists(p)]
        outside = [p for p in tested if not under(p, roots)]
        if missing or outside:
            why = ("tested path does not exist: %s" % missing[0] if missing
                   else "tested path outside allowed roots: %s" % outside[0])
            finish("not started", reason=why)
            raise Refused(7, why)
    if received:
        if queue.stop_requested():
            finish("not started", reason="stop requested")
            raise Refused(3, "stop requested")
        finish("not started", reason="signal before start")
        raise SystemExit(128 + received[0])

    log = Path(args.log).resolve()
    log.parent.mkdir(parents=True, exist_ok=True)
    with queue.locked():
        if queue.stop_requested():
            current = read_json(record_path)
            current.update({"status": "not started", "reason": "stop requested", "ended": now()})
            write_json(record_path, current)
            queue.event("run %s not started (stop requested)" % run_id)
            raise Refused(3, "stop requested")
        try:
            with open(log, "ab") as out:
                launch = ([SANDBOX, "-p", NO_NETWORK] + command) if network == "blocked" else command
                proc = subprocess.Popen(launch, cwd=cwd, stdout=out, stderr=subprocess.STDOUT,
                                        stdin=subprocess.DEVNULL, start_new_session=True,
                                        env=marked_env)
        except OSError as error:
            current = read_json(record_path)
            current.update({"status": "not started", "reason": "could not start: %s" % error,
                            "ended": now()})
            write_json(record_path, current)
            queue.event("run %s not started (%s)" % (run_id, error))
            raise Refused(2, "could not start the command: %s" % error)
        current = read_json(record_path)
        current.update({"status": "running", "tested_paths": tested, "pid": proc.pid,
                        "pgid": proc.pid, "process_start": process_start(proc.pid)})
        write_json(record_path, current)
        queue.event("run %s running pid %d" % (run_id, proc.pid))

    tracked = {}

    def forward(signum, frame):
        received.append(signum)
        try:
            os.killpg(proc.pid, signum)
        except OSError:
            pass
        procs = live_table(snapshot_or_none()) or {}  # unreadable: the group signal above only
        for pid, start in list(tracked.items()):
            live = procs.get(pid)
            if live is not None and live[3] == start:
                try:
                    os.kill(pid, signum)
                except OSError:
                    pass
    signal.signal(signal.SIGTERM, forward)
    signal.signal(signal.SIGINT, forward)
    for signum in list(received):
        forward(signum, None)
    def track():
        """Record every descendant of CMD (pid and start time), so a process that later leaves
        the group and the session is still known."""
        _SNAPSHOT["procs"] = None
        procs = live_table(snapshot_or_none())
        if procs is None:
            # Revision 13: nothing new can be learned now; the gap is recorded, and the end of
            # the run cannot then establish that nothing was left behind.
            gaps.append(now())
            return
        known = {proc.pid} | {pid for pid, start in tracked.items()
                              if pid in procs and procs[pid][3] == start}
        grew = True
        while grew:
            more = {pid for pid, info in procs.items() if pid not in known and info[0] in known}
            known |= more
            grew = bool(more)
        new = {pid: procs[pid][3] for pid in known
               if pid in procs and pid not in tracked and procs[pid][3] is not None}
        if new:
            tracked.update(new)
            with queue.locked():
                current = read_json(record_path)
                current["tracked"] = [{"pid": p, "start": t} for p, t in sorted(tracked.items())]
                write_json(record_path, current)

    gaps = []
    deadline = time.time() + timeout
    timed_out = False
    while True:
        try:
            code = proc.wait(timeout=1.0)
            break
        except subprocess.TimeoutExpired:
            track()
            if time.time() > deadline:
                timed_out = True
                break
    if timed_out:
        # End the run and everything identified as its own: TERM, up to 10 s, then KILL.
        record_now = read_json(record_path)
        procs = live_table(snapshot_or_none()) or {}
        targets = set(survivors_or_unknown(record_now) or ()) | {
            pid for pid, start in tracked.items() if procs.get(pid, [None] * 4)[3] == start}
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(proc.pid, sig)
            except OSError:
                pass
            for pid in targets:
                try:
                    os.kill(pid, sig)
                except OSError:
                    pass
            end = time.time() + (10 if sig == signal.SIGTERM else 2)
            while time.time() < end:
                _SNAPSHOT["procs"] = None
                procs = live_table(snapshot_or_none())
                if proc.poll() is not None and procs is not None and not [
                        p for p in targets if p in procs]:
                    break
                time.sleep(0.2)
            else:
                continue
            break
        try:
            code = proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            code = -signal.SIGKILL
        _SNAPSHOT["procs"] = None
        left = survivors_or_unknown(read_json(record_path))
        if left is None:
            left = UNKNOWN_LEFT
        finish("timed out", exit_code=124, left_after_timeout=left, signals=received,
               **({"process_table_unreadable": gaps} if gaps else {}))
        queue.event("run %s timed out after %d s%s" % (run_id, timeout,
                                                       (", still alive: %s" % left) if left else ""))
        raise SystemExit(124)
    track()
    exit_code = code if code >= 0 else 128 - code
    _SNAPSHOT["procs"] = None
    remaining = survivors_or_unknown(read_json(record_path))
    if remaining is None:
        gaps.append(now())
    # Revision 13: when the process table could not be read, leftovers cannot be excluded;
    # the gaps are recorded with the run, and readiness never accepts such a run.
    left = group_alive(proc.pid) or bool(remaining)
    interrupted = bool(received) or queue.stop_requested()
    finish("interrupted" if interrupted else "completed", exit_code=exit_code,
           group_left_running=left, signals=received,
           **({"process_table_unreadable": gaps} if gaps else {}))
    raise SystemExit(exit_code)


def identity_known(pid, recorded_start):
    """Signal only a process whose recorded start time is known and still matches."""
    if recorded_start in (None, UNREADABLE):
        return False
    return process_start(pid) == recorded_start


def signal_live(queue):
    """Send SIGTERM to every run's confirmed target and survivors. Checkpoint 3: records are
    taken one at a time, so a record whose liveness needs the process table, when the table
    is unreadable, no longer prevents signalling another run whose wrapper or group is
    confirmed by `ps -p` alone; the unreadable table is raised after the loop (stop then
    still ends with exit 17, never Stopped)."""
    unreadable = None
    for _, record in queue.runs():
        try:
            if not run_live(record):
                continue
        except ProcessTableUnreadable as error:
            unreadable = unreadable or error
            continue
        if record.get("status") == "starting":
            target = ("pid", record.get("wrapper_pid"))
            known = identity_known(record.get("wrapper_pid"), record.get("wrapper_start"))
        else:
            target = ("group", record.get("pgid"))
            known = identity_known(record.get("pid"), record.get("process_start")) or (
                record.get("group_left_running") and record.get("status") != "running")
        if known:
            try:
                if target[0] == "pid":
                    os.kill(target[1], signal.SIGTERM)
                else:
                    os.killpg(target[1], signal.SIGTERM)
                queue.event("stop sent SIGTERM to %s %s of run %s"
                            % (target[0], target[1], record["run_id"]))
            except OSError:
                pass
        else:
            queue.event("stop did not signal %s %s of run %s: identity not confirmed"
                        % (target[0], target[1], record["run_id"]))
        try:
            left = survivors(record)
        except ProcessTableUnreadable as error:
            unreadable = unreadable or error
            continue
        for pid in left:
            try:
                os.kill(pid, signal.SIGTERM)
                queue.event("stop sent SIGTERM to surviving process %d of run %s" % (pid, record["run_id"]))
            except OSError:
                pass
    if unreadable is not None:
        raise unreadable


def cmd_stop(args):
    queue = open_queue(args)
    with queue.locked():
        flag = queue.dir / "STOP_REQUESTED"
        if not flag.exists():
            write_text(flag, "time: %s\npid: %d\ncwd: %s\n" % (now(), os.getpid(), os.getcwd()))
            queue.event("stop requested")
        try:
            signal_live(queue)
            # Checkpoint 3: read the table even when no run record needs it, so that stop
            # never prints Stopped under an unreadable table (report would say Stopping).
            _SNAPSHOT["procs"] = None
            snapshot()
            left = queue.remaining()
        except ProcessTableUnreadable as error:
            stop_unknown(queue, error)
    deadline = time.time() + 5
    try:
        while left and time.time() < deadline:
            time.sleep(0.2)
            _SNAPSHOT["procs"] = None
            queue.reload()
            if queue.worker_active():
                left = queue.remaining()
                break
            left = queue.remaining()
    except ProcessTableUnreadable as error:
        stop_unknown(queue, error)
    if left:
        print("Stopping: " + "; ".join(left))
        raise SystemExit(10)
    print("Stopped")


def stop_unknown(queue, error):
    """Revision 13: the stop is recorded, but whether tracked processes remain is unknown:
    never Stopped. Prints Stopping with the reason and refuses (exit 17)."""
    queue.event("stop: %s" % error)
    worker = queue.worker_active()
    print("Stopping: %s%s" % (error, ("; worker for job %s registered, not yet ended" % worker["job"])
                              if worker else ""))
    raise error


def cmd_transfer(args):
    """Revision 13: hand jobs over from the auto queue to an attended successor record.
    The job (or every job that is not final, with --all) becomes Transferred, a final state:
    its decision_needed is no longer offered, the successor directory is recorded, and its
    history, runs and evidence are kept. Refused while a worker is registered for the job or
    a run of it may be alive (an unreadable process table counts as may be alive)."""
    queue = open_queue(args)
    if bool(args.job) == bool(args.all):
        raise Refused(2, "give one job or --all")
    if not args.note or not args.note.strip():
        raise Refused(2, "transfer needs --note saying why and what continues attended")
    successor = Path(args.to)
    if not successor.is_dir():
        raise Refused(2, "successor record directory %s does not exist" % args.to)
    successor = successor.resolve()
    if successor == queue.dir or queue.dir in successor.parents:
        raise Refused(2, "the successor record must be outside the queue directory")
    with queue.locked():
        if args.all:
            jobs = [j for j in queue.data["jobs"] if j["state"] not in FINAL_STATES]
            if not jobs:
                raise Refused(8, "no job left to transfer")
        else:
            jobs = [queue.job(args.job)]
        worker = queue.worker_active()
        _SNAPSHOT["procs"] = None
        snapshot()  # Checkpoint 3: unreadable refuses (17) even with no run record
        for job in jobs:
            if job["state"] in FINAL_STATES:
                raise Refused(8, "job %s is %s; a final job is not transferred" % (job["id"], job["state"]))
            if worker and worker["job"] == job["id"]:
                raise Refused(9, "worker %s for job %s is still registered; end it first"
                              % (worker["id"], job["id"]))
            live = queue.live_runs(job["id"])
            if live:
                raise Refused(13, "a run of job %s still has live processes (%s); run cleanup first"
                              % (job["id"], live[0][2]))
        for job in jobs:
            old = job["state"]
            entry = {"time": now(), "from_state": old, "to": str(successor), "note": args.note}
            if job.get("decision_needed"):
                entry["decision_needed_at_transfer"] = job.pop("decision_needed")
            for key in ("packet", "packet_manifest_sha256", "approval_quote",
                        "approval_manifest_sha256"):
                if job.get(key):
                    entry[key] = job.pop(key)
            job.pop("waiting", None)
            job["transfer"] = entry
            job["state"] = "Transferred"
            job["notes"].append({"time": entry["time"], "state": "Transferred", "note": args.note})
            queue.event("job %s %s -> Transferred to %s note: %s"
                        % (job["id"], old, successor, args.note.replace("\n", " ")))
        queue.save()
    for job in jobs:
        print("%s Transferred to %s" % (job["id"], successor))


def cmd_reconcile(args):
    queue = open_queue(args)
    found = 0
    with queue.locked():
        for path, record in queue.runs():
            if record.get("status") in ("starting", "running") and not run_live(record):
                record["status"] = "unknown"
                record["reconciled"] = now()
                write_json(path, record)
                queue.event("reconcile run %s unknown (process gone)" % record["run_id"])
                print("run %s: process gone, outcome unknown" % record["run_id"])
                found += 1
        worker = queue.worker_active()
    if worker:
        print("worker %s for job %s registered, not ended" % (worker["id"], worker["job"]))
    print("reconcile: %d run(s) marked unknown" % found)


def cmd_wait(args):
    queue = open_queue(args)
    with queue.locked():
        job = queue.job(args.job)
        if "Waiting" not in TRANSITIONS[job["state"]]:
            raise Refused(8, "illegal transition %s -> Waiting" % job["state"])
        if args.until:
            job["waiting"] = {"until": args.until, "since": now()}
            reason = "until %s" % args.until
        else:
            job["waiting"] = {"until": None, "retries_left": args.max_retries, "since": now()}
            reason = "reset time unknown, %d retries" % args.max_retries
        old = job["state"]
        job["state"] = "Waiting"
        queue.save()
        queue.event("job %s %s -> Waiting %s" % (job["id"], old, reason))
    print("%s Waiting %s" % (job["id"], reason))


def cmd_install(args):
    queue = open_queue(args)
    if args.action == "show":
        print(json.dumps(queue.data.get("installation") or {"state": "unchecked"}, sort_keys=True))
        return
    if not args.script or not args.expect:
        raise Refused(2, "install check needs --script and --expect")
    baseline = (queue.data.get("installation") or {}).get("expected")
    if baseline and args.expect != baseline and not args.new_baseline:
        raise Refused(2, "expected value differs from the recorded starting state; a new starting "
                         "state is an attended step (--new-baseline)")
    try:
        out = subprocess.run([args.script], capture_output=True, text=True)
        returncode, stdout = out.returncode, out.stdout
    except OSError as error:
        returncode, stdout = 126, ""
        start_error = str(error)
    else:
        start_error = None
    found = normalized = None
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == "FINGERPRINT":
            found = fields[1]
        elif len(fields) == 2 and fields[0] == "NORMALIZED":
            normalized = fields[1]
    # Availability needs the exact fingerprint, always (Revision 13 changes only the reason).
    ok = returncode == 0 and found is not None and found == args.expect
    known_codes = {3: "a process runs from the installation",
                   4: "the process table could not be read"}
    previous = queue.data.get("installation") or {}
    recorded_normalized = args.expect_normalized or previous.get("expected_normalized")
    if start_error:
        reason = "its fingerprint script could not be started"
    elif returncode:
        reason = "its fingerprint script failed (exit %d%s)" % (
            returncode, (": " + known_codes[returncode]) if returncode in known_codes else "")
    elif found is None:
        reason = "its fingerprint script printed no fingerprint"
    elif normalized and recorded_normalized and normalized == recorded_normalized:
        # Revision 13: report-only classification. The installation stays unavailable; only
        # an attended validation may record a new exact value.
        reason = ("timestamp-only differences: attended validation needed (the exact fingerprint "
                  "differs from the recorded starting state; the normalized value, which ignores "
                  "only known timestamp fields, matches)")
    else:
        reason = "it no longer matches its recorded starting state"
    with queue.locked():
        kept = baseline if (baseline and not ok) else args.expect
        if ok:
            kept_normalized = normalized
        else:
            kept_normalized = previous.get("expected_normalized") or args.expect_normalized
        queue.data["installation"] = {
            "state": "available" if ok else "unavailable", "checked": now(), "script": args.script,
            "expected": kept, "found": found, "reason": "" if ok else reason}
        if kept_normalized:
            queue.data["installation"]["expected_normalized"] = kept_normalized
        if normalized:
            queue.data["installation"]["found_normalized"] = normalized
        queue.save()
        queue.event("installation %s%s%s" % (queue.data["installation"]["state"],
                                             "" if ok else ": " + reason,
                                             " (new starting state)" if args.new_baseline else ""))
    print("installation %s" % queue.data["installation"]["state"])
    if not ok:
        raise SystemExit(12)


def cmd_provisional(args):
    queue = open_queue(args)
    with queue.locked():
        job = queue.job(args.job)
        choices = job.setdefault("provisional", [])
        number = len(choices) + 1
        choices.append({"number": number, "state": "pending", "time": now(), "note": args.note})
        queue.save()
        queue.event("job %s provisional P%d: %s" % (job["id"], number, args.note.replace("\n", " ")))
    print("P%d" % number)


def cmd_decide(args):
    queue = open_queue(args)
    with queue.locked():
        job = queue.job(args.job)
        choices = {p["number"]: p for p in job.get("provisional", [])}
        choice = choices.get(args.number)
        if choice is None:
            raise Refused(2, "job %s has no provisional choice P%d" % (job["id"], args.number))
        if choice["state"] != "pending":
            raise Refused(8, "P%d of job %s is already %s" % (args.number, job["id"], choice["state"]))
        if not names_job(args.quote, job["id"]) or not names_job(args.quote, "P%d" % args.number):
            raise Refused(8, "the decision must quote the Developer naming %s and P%d"
                          % (job["id"], args.number))
        choice.update({"state": "adopted" if args.adopt else "rejected", "decided": now(),
                       "quote": args.quote})
        queue.save()
        queue.event("job %s P%d %s; quote: %s" % (job["id"], args.number, choice["state"],
                                                 args.quote.replace("\n", " ")))
    print("P%d %s" % (args.number, choice["state"]))


def cmd_cleanup(args):
    queue = open_queue(args)
    with queue.locked():
        job = queue.job(args.job)
        # Checkpoint 3: the table is read even when the job has no run record; unreadable
        # refuses (17) rather than printing Cleaned.
        _SNAPSHOT["procs"] = None
        snapshot()
        for _, record in queue.runs(job["id"]):
            if (record.get("status") == "running" and same_process(record.get("pid"),
                                                                   record.get("process_start"))) or (
                    record.get("status") == "starting" and same_process(record.get("wrapper_pid"),
                                                                        record.get("wrapper_start"))):
                raise Refused(13, "run %s of job %s is still %s" % (record["run_id"], job["id"],
                                                                     record["status"]))
        for _, record in queue.runs(job["id"]):
            for pid in survivors(record):
                try:
                    os.kill(pid, signal.SIGTERM)
                    queue.event("cleanup sent SIGTERM to process %d of run %s" % (pid, record["run_id"]))
                except OSError:
                    pass
    deadline = time.time() + 5
    while True:
        _SNAPSHOT["procs"] = None
        left = [why for _, _, why in queue.live_runs(args.job)]
        if not left or time.time() > deadline:
            break
        time.sleep(0.2)
    if left:
        print("Still alive: " + "; ".join(left))
        raise SystemExit(10)
    print("Cleaned")


def cmd_workers(args):
    queue = open_queue(args)
    periods, open_ = [], {}
    for line in (queue.dir / "EVENTS.log").read_text(encoding="utf-8").splitlines():
        m = re.match(r"(\S+) worker (start|end) (\S+) job (\S+)(?: outcome (.*))?$", line)
        if not m:
            continue
        time_, kind, wid, job, outcome = m.groups()
        if kind == "start":
            open_[wid] = len(periods)
            periods.append([wid, job, time_, None, ""])
        elif wid in open_:
            periods[open_.pop(wid)][3:5] = [time_, outcome or ""]
    overlap = []
    for i, a in enumerate(periods):
        for b in periods[i + 1:]:
            if b[2] < (a[3] or "9999"):
                overlap.append("%s/%s" % (a[0], b[0]))
    for wid, job, start, end, outcome in periods:
        print("%s job %s %s -> %s %s" % (wid, job, start, end or "running", outcome))
    print("overlap: %s" % (", ".join(overlap) if overlap else "none"))


def cmd_report(args):
    queue = open_queue(args)
    lines = []
    # Revision 13: with an unreadable process table nothing about live processes is known:
    # never Stopped, no readiness, and every job with runs shows unknown; exit 17 after printing.
    table_error = None
    try:
        snapshot()
    except ProcessTableUnreadable as error:
        table_error = error
    if queue.stop_requested():
        lines.append("**Stopping**" if table_error or queue.remaining() else "**Stopped**")
        lines.append("")
    if table_error:
        lines.append("**Process table unreadable:** %s" % table_error)
        lines.append("")
    installation = queue.data.get("installation") or {}
    if installation.get("state") == "unavailable":
        lines.append("**Test installation unavailable:** %s" % installation.get("reason", "-"))
        lines.append("")
    lines.append("| Job | Result | Decision or next step |")
    lines.append("| --- | --- | --- |")
    for job in queue.data["jobs"]:
        runs = [r for _, r in queue.runs(job["id"])]
        state = job["state"]
        open_runs = [] if table_error else [
            r for r in runs if r.get("status") in ("starting", "running") or run_live(r)]
        if table_error and runs:
            result = "unknown (%s; process table unreadable)" % state.lower()
            step = "report again once the process table can be read"
        elif open_runs:
            alive = [r for r in open_runs if run_live(r)]
            result = ("running %s" % alive[0]["label"] if alive
                      else "unknown (a run ended without a completion record)")
            step = "wait for the run" if alive else "reconcile before deciding"
        elif state == "Ready for approval":
            failed = queue.readiness(job)
            pending = [p for p in job.get("provisional", []) if p["state"] == "pending"]
            if failed is None and pending:
                result = "local checks passed; candidate and ticket saved"
                step = "decide %s first, then approve / inspect / revise / discard" % ", ".join(
                    "P%d" % p["number"] for p in pending)
            elif failed is None:
                result = "local checks passed; candidate and ticket saved"
                step = "approve for integration / inspect / revise / discard"
            else:
                result = "results stale or incomplete (%s)" % failed
                # Revision 13: the helper refuses only what is in the frozen bytes, which an
                # ERRATA.md note cannot change; for a clerical slip the new freeze reuses the
                # same evidence with no rerun (WORKER.md).
                step = ("freeze a corrected packet in a new directory and mark it ready again "
                        "(clerical slip: same evidence, note in ERRATA.md, no rerun)"
                        if failed.startswith(("packet", "SCREENING.txt")) else "rerun the affected checks")
        elif state == "Transferred":
            result = "transferred to %s" % job["transfer"]["to"]
            step = "continued attended there; nothing to decide in this queue"
        elif state == "Blocked":
            result = "blocked: %s" % last_note(job)
            step = (("decision needed: " + job["decision_needed"]) if job.get("decision_needed")
                    else "nothing to decide; revise or discard if you wish")
        elif state == "Waiting":
            wait = job.get("waiting") or {}
            result = "waiting %s" % (("until " + wait["until"]) if wait.get("until")
                                     else "for an unknown reset time")
            step = "resume with /gc auto resume after the reset (automatic continuation is not supported)"
        else:
            result = state.lower()
            step = {"Queued": "not started", "Preparing": "in preparation", "Testing": "testing",
                    "Approved": "approved; integration on your instruction",
                    "Published": "published", "Discarded": "discarded"}.get(state, "-")
        title = job["id"] + ((" " + job["title"]) if job.get("title") else "")
        lines.append("| %s | %s | %s |" % (cell(title), cell(result), cell(step)))
    print("\n".join(lines))
    if table_error:
        raise table_error


def last_note(job):
    return job["notes"][-1]["note"] if job.get("notes") else "no reason recorded"


def cell(text):
    return str(text).replace("|", "/").replace("\n", " ")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="auto_queue.py", allow_abbrev=False)
    parser.add_argument("--queue")
    parser.add_argument("--home")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("init", allow_abbrev=False)
    p.add_argument("--records", required=True)
    p.add_argument("--queue-id", required=True)
    p.add_argument("--grant", required=True)
    p.add_argument("--jobs", required=True)
    p.add_argument("--lock-dir")
    p.add_argument("--min-free-gib", type=float, default=8.0)
    p.add_argument("--space-path")
    p.add_argument("--allowed-root", action="append")
    p.add_argument("--default-timeout", type=int, default=3600)
    p.set_defaults(func=cmd_init)

    sub.add_parser("status", allow_abbrev=False).set_defaults(func=cmd_status)

    p = sub.add_parser("set", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("state")
    p.add_argument("--note")
    p.add_argument("--packet")
    p.add_argument("--quote")
    p.add_argument("--manifest")
    p.add_argument("--decision")
    p.set_defaults(func=cmd_set)

    p = sub.add_parser("criterion", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("--file", required=True)
    p.set_defaults(func=cmd_criterion)

    p = sub.add_parser("candidate", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("--repo", action="append")
    p.add_argument("--commit", action="append")
    p.set_defaults(func=cmd_candidate)

    p = sub.add_parser("worker", allow_abbrev=False)
    p.add_argument("action", choices=("start", "end"))
    p.add_argument("job")
    p.add_argument("--id", required=True)
    p.add_argument("--outcome", default="finished")
    p.set_defaults(func=cmd_worker)

    sub.add_parser("check-stop", allow_abbrev=False).set_defaults(func=cmd_check_stop)

    p = sub.add_parser("lock", allow_abbrev=False)
    p.add_argument("action", choices=("take", "release", "show"))
    p.set_defaults(func=cmd_lock)

    p = sub.add_parser("run", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("label")
    p.add_argument("--log", required=True)
    p.add_argument("--probe")
    p.add_argument("--cwd")
    p.add_argument("--timeout", type=int)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("transfer", allow_abbrev=False)
    p.add_argument("job", nargs="?")
    p.add_argument("--all", action="store_true")
    p.add_argument("--to", required=True)
    p.add_argument("--note", required=True)
    p.set_defaults(func=cmd_transfer)

    sub.add_parser("stop", allow_abbrev=False).set_defaults(func=cmd_stop)
    sub.add_parser("reconcile", allow_abbrev=False).set_defaults(func=cmd_reconcile)

    p = sub.add_parser("wait", allow_abbrev=False)
    p.add_argument("job")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--until")
    group.add_argument("--unknown", action="store_true")
    p.add_argument("--max-retries", type=int, default=8)
    p.set_defaults(func=cmd_wait)

    sub.add_parser("report", allow_abbrev=False).set_defaults(func=cmd_report)

    sub.add_parser("workers", allow_abbrev=False).set_defaults(func=cmd_workers)

    p = sub.add_parser("provisional", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("--note", required=True)
    p.set_defaults(func=cmd_provisional)

    p = sub.add_parser("decide", allow_abbrev=False)
    p.add_argument("job")
    p.add_argument("number", type=int)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--adopt", action="store_true")
    group.add_argument("--reject", action="store_true")
    p.add_argument("--quote", required=True)
    p.set_defaults(func=cmd_decide)

    p = sub.add_parser("cleanup", allow_abbrev=False)
    p.add_argument("job")
    p.set_defaults(func=cmd_cleanup)

    p = sub.add_parser("install", allow_abbrev=False)
    p.add_argument("action", choices=("check", "show"))
    p.add_argument("--script")
    p.add_argument("--expect")
    p.add_argument("--expect-normalized")
    p.add_argument("--new-baseline", action="store_true")
    p.set_defaults(func=cmd_install)

    argv = list(sys.argv[1:] if argv is None else argv)
    command = []
    if "--" in argv:
        split = argv.index("--")
        argv, command = argv[:split], argv[split + 1:]
    args = parser.parse_args(argv)
    args.command = command
    try:
        args.func(args)
    except Refused as refusal:
        print("REFUSED: %s" % refusal, file=sys.stderr)
        return refusal.code
    return 0


if __name__ == "__main__":
    sys.exit(main())
