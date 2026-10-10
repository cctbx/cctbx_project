#!/usr/bin/env python3
"""Check screen form, frozen bytes, source identity, and local registration."""

import sys
# When run as a script, keep this directory out of module resolution before
# importing the standard library. The documented -I invocation also ignores
# PYTHONPATH and the current directory.
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
import hashlib
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

HEADINGS = {
    "plan": ("GOAL", "SCOPE", "CHECK", "LIMITS", "DECISION"),
    "result": ("CHANGED", "CHECKED", "LIMITS", "GATE", "DECISION"),
    "stop": ("WHY", "PRESERVED", "NEEDED"),
    "publication": ("BATCH", "SERVER CHECK", "LIMITS", "GATE", "DECISION"),
    "auto_start": ("JOBS", "CHECKS", "INSTALLATION", "RESULTS", "STOP", "NOT AUTHORIZED"),
    "auto_report": ("JOBS", "IN SIMPLE WORDS", "APPROVAL SUMMARIES", "NEXT"),
}
ACTIONS = {
    "plan": "ACTION: APPROVE / REVISE / STOP",
    "result": "ACTION: INTEGRATE / REVISE / DISCARD",
    "stop": "ACTION: DECIDE / CANCEL",
    "publication": "ACTION: PUBLISH / HOLD",
    "auto_start": "ACTION: NONE NEEDED (reply stop to stop the queue)",
    "auto_report": "ACTION: APPROVE / INSPECT / REVISE / DISCARD",
}
ATTENTION = {
    "auto_start": "**PLEASE READ — NO ACTION NEEDED**",
    "auto_report": "**PLEASE READ — ACTION NEEDED**",
}
BUDGET = {"auto_report": (90, 900)}
AUTO_QUEUE = Path(__file__).resolve().parent.parent.parent / "auto_queue.py"
HEX = r"[0-9a-f]{40}(?:[0-9a-f]{24})?"
SHA = r"[0-9a-f]{64}"
NAME = r"[a-z0-9][a-z0-9-]{0,39}"
DISPLAY_ID = re.compile(r"(?<![A-Za-z0-9])[0-9a-f]{8,64}(?:…[0-9a-f]{4,64})?(?![A-Za-z0-9])")
MIN_CLAUDE_VERSION = (2, 1, 281)
# The CLAUDE_CODE_ENTRYPOINT value observed in Claude app (Code tab) sessions on
# 2026-10-06. Any other or missing value takes the Terminal path below.
APP_ENTRYPOINT = "claude-desktop"


def fail(reason):
    raise ValueError(reason)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def regular_files(root):
    absolute = root.absolute()
    if any(component.is_symlink() for component in (absolute, *absolute.parents)):
        fail("symbolic link in evidence path")
    if root.is_symlink() or not root.is_dir():
        fail("evidence must be an ordinary directory")
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            fail(f"symbolic link in evidence: {path.relative_to(root)}")
        mode = path.stat().st_mode
        if stat.S_ISDIR(mode):
            continue
        if not stat.S_ISREG(mode) or path.stat().st_nlink != 1:
            fail(f"non-regular or linked file: {path.relative_to(root)}")
        name = path.relative_to(root).as_posix()
        if any(ord(char) < 32 or ord(char) == 127 for char in name):
            fail("control character in evidence path")
        files.append((name, path))
    return files


def freeze(root):
    manifest = root / "MANIFEST.sha256"
    if manifest.exists() or manifest.is_symlink():
        fail("manifest already exists; make a fresh evidence directory")
    files = regular_files(root)
    if not files:
        fail("empty evidence directory")
    lines = [f"{digest(path.read_bytes())}  {name}\n" for name, path in files]
    with manifest.open("x", encoding="utf-8", newline="\n") as stream:
        stream.writelines(lines)
    print(f"FROZEN ({len(files)} files)")


def verify(root):
    files = regular_files(root)
    manifest = root / "MANIFEST.sha256"
    if ("MANIFEST.sha256", manifest) not in files:
        fail("missing manifest")
    raw = manifest.read_bytes()
    try:
        contents = raw.decode("utf-8")
    except UnicodeDecodeError:
        fail("manifest is not UTF-8")
    expected = {}
    for line in contents.splitlines(keepends=True):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^\n\r]+)\n", line)
        if not match or match[2] in expected or match[2] == "MANIFEST.sha256":
            fail("malformed, duplicate, or self-listed manifest entry")
        expected[match[2]] = match[1]
    actual = {name: path for name, path in files if name != "MANIFEST.sha256"}
    if not expected or expected.keys() != actual.keys():
        fail("missing or extra evidence file")
    for name, path in actual.items():
        if digest(path.read_bytes()) != expected[name]:
            fail(f"changed evidence file: {name}")
    return digest(raw)


def verify_source(root):
    """Check the complete installed kit, including files absent from its manifest.

    Invoke this script with Python's -I -B flags so an unlisted module next
    to this script cannot run while the source is being checked.
    """
    files = regular_files(root)
    manifest = root / "SOURCE_MANIFEST.sha256"
    if ("SOURCE_MANIFEST.sha256", manifest) not in files:
        fail("missing source manifest")
    raw = manifest.read_bytes()
    try:
        contents = raw.decode("utf-8")
    except UnicodeDecodeError:
        fail("source manifest is not UTF-8")
    expected = {}
    for line in contents.splitlines(keepends=True):
        match = re.fullmatch(r"([0-9a-f]{64})  \./([^\n\r]+)\n", line)
        if (not match or match[2] in expected or
                match[2] == "SOURCE_MANIFEST.sha256" or
                any(part in ("", ".", "..") for part in match[2].split("/"))):
            fail("malformed, duplicate, or self-listed source entry")
        expected[match[2]] = match[1]
    actual = {name: path for name, path in files
              if name != "SOURCE_MANIFEST.sha256"}
    precompiled = [name for name in actual
                   if name not in expected
                   and precompiled_bytecode(name, expected, actual[name])]
    for name in precompiled:
        del actual[name]
    if not expected or expected.keys() != actual.keys():
        fail("missing or extra source file")
    for name, path in actual.items():
        if digest(path.read_bytes()) != expected[name]:
            fail(f"changed source file: {name}")
    if precompiled:
        print(f"NOTE ignored {len(precompiled)} precompiled bytecode file(s) beside "
              "listed modules; these tools never load them")
    return digest(raw)


def precompiled_bytecode(name, expected, path):
    """True only for the file compileall (libtbx.py_compile_all) writes for a
    listed module: <dir>/__pycache__/<stem>.<tag>[.opt-N].pyc beside a listed
    <dir>/<stem>.py, in a form Python checks against that source (timestamp
    or checked hash). Any other unlisted file, including bytecode for an
    unlisted module, a sourceless .pyc, or unchecked-hash bytecode (which
    Python would run without consulting the source), is still refused. The
    contents of accepted bytecode are not verified; these tools never load
    it."""
    parts = name.split("/")
    if len(parts) < 2 or parts[-2] != "__pycache__":
        return False
    match = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)\.[a-z]+-?[0-9]+(?:\.opt-[12])?\.pyc",
                         parts[-1])
    if not match or "/".join(parts[:-2] + [match[1] + ".py"]) not in expected:
        return False
    with path.open("rb") as stream:
        header = stream.read(8)
    # PEP 552 flags word: bit 0 = hash-based, bit 1 = check against source.
    flags = int.from_bytes(header[4:8], "little") if len(header) == 8 else -1
    return flags in (0, 3)


def claude_version(executable):
    """Return (version tuple, reported text) from `executable --version`, or
    (None, problem) when no version can be established: ("run", error text),
    ("timeout", seconds), ("exit", code) or ("format", output). A failed
    command yields no version even when its output contains a plausible
    number."""
    try:
        result = subprocess.run([executable, "--version"], capture_output=True,
                                text=True, timeout=10)
    except OSError as error:
        return None, ("run", error.strerror or str(error))
    except subprocess.TimeoutExpired:
        return None, ("timeout", 10)
    if result.returncode != 0:
        return None, ("exit", result.returncode)
    reported = result.stdout.strip()
    match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)(?: \(Claude Code\))?", reported)
    if match is None:
        return None, ("format", reported[:100])
    return tuple(int(part) for part in match.groups()), reported


def check_claude_version():
    """Gate the Claude Code that runs this session (minimum 2.1.281).

    Terminal session: the `claude` command on PATH, as before; any problem
    reading its version is a failure. Claude app session, recognized by
    CLAUDE_CODE_ENTRYPOINT=claude-desktop (the marker observed in app
    sessions on 2026-10-06; any other or missing value takes the Terminal
    path): the app's own Claude Code engine named by CLAUDE_CODE_EXECPATH.
    A separately installed `claude` never decides for the app. When the app
    engine's version cannot be established, print one NOT CHECKED line with
    the reason and continue; a version that is read and is below the
    minimum still fails.
    """
    if os.environ.get("CLAUDE_CODE_ENTRYPOINT") == APP_ENTRYPOINT:
        check_app_engine_version(os.environ.get("CLAUDE_CODE_EXECPATH"))
        return
    executable = shutil.which("claude")
    if executable is None:
        fail("Claude Code CLI not found on PATH; cannot verify minimum version 2.1.281")
    version, detail = claude_version(executable)
    if version is None:
        kind, value = detail
        fail({"run": f"could not run Claude Code CLI at {executable}: {value}",
              "timeout": "Claude Code CLI version check timed out",
              "exit": f"Claude Code CLI version check failed (exit {value})",
              "format": f"unrecognized Claude Code CLI version: {value!r}"}[kind])
    if version < MIN_CLAUDE_VERSION:
        fail(f"Claude Code CLI {detail} is below 2.1.281; run `claude update`, "
             "restart Claude Code, and check again")
    print(f"VERIFIED Claude Code CLI {detail} at {executable}, the claude command on "
          "this shell's PATH; 2.1.281 or newer is required")


def check_app_engine_version(executable):
    """The Claude app branch of check_claude_version (see its docstring)."""
    if not executable:
        reason = ("the app did not say where its Claude Code engine is "
                  "(CLAUDE_CODE_EXECPATH is not set)")
    else:
        version, detail = claude_version(executable)
        if version is not None:
            if version < MIN_CLAUDE_VERSION:
                fail(f"Claude app's Claude Code engine {detail} is below 2.1.281; "
                     "update the Claude app, restart it, and check again")
            print(f"VERIFIED Claude app's Claude Code engine {detail} at {executable}; "
                  "2.1.281 or newer is required")
            return
        kind, value = detail
        engine = f"the app's Claude Code engine at {executable}"
        reason = {"run": f"{engine} could not be run ({value})",
                  "timeout": f"{engine} did not answer a version request within {value} seconds",
                  "exit": f"the version command of {engine} failed (exit {value})",
                  "format": f"{engine} gave an unrecognized version answer: {value!r}"}[kind]
    print(f"NOT CHECKED Claude Code version: this is a Claude app session and {reason}. "
          "Continuing without confirming 2.1.281 or newer; a separately installed "
          "claude command does not decide for the app.")


def register_skill(source):
    """Register the verified source once, without following an occupied link.

    os.symlink creates the *exact* destination exclusively. Unlike `ln -s`,
    it cannot create a child inside an existing directory or directory link.
    """
    verify_source(source)
    check_claude_version()
    config = Path(os.environ.get("CLAUDE_CONFIG_DIR") or Path.home() / ".claude")
    if not config.is_absolute():
        fail("CLAUDE_CONFIG_DIR must be an absolute path")
    # Do not let mkdir create a missing prefix that changes how a later `..`
    # resolves. The preflight and the write must see the same parent path.
    if ".." in config.parts:
        fail("CLAUDE_CONFIG_DIR must not contain '..' components")
    skills = config / "skills"
    destination = skills / "guided_coding"
    if config.is_symlink() or skills.is_symlink():
        fail("personal configuration or skills directory is a symbolic link")
    if destination.exists() or destination.is_symlink():
        fail(f"skill destination already exists; inspect it without overwriting: {destination}")
    # Do not place the personal link inside the source, even through an
    # existing alias in a configuration path.
    if source.resolve() in destination.resolve(strict=False).parents:
        fail("personal skill destination is inside the verified source")
    skills.mkdir(parents=True, exist_ok=True)
    os.symlink(source.absolute(), destination, target_is_directory=True)
    print(f"REGISTERED {destination} -> {source.absolute()}")
    print("NEXT: In the project you choose, use /gc setup to recover saved settings "
          "and configure relevant work. Registration has not adopted a project.")


def blocks(kind, lines):
    headings = HEADINGS[kind]
    positions = []
    for heading in headings:
        occurrences = [i for i, value in enumerate(lines) if value == heading]
        if len(occurrences) != 1:
            fail(f"expected exactly one {heading} heading")
        positions.append(occurrences[0])
    if positions != sorted(positions) or positions[0] != 1:
        fail("headings out of order or content before first heading")
    if lines[-1] != ACTIONS[kind] or sum(x.startswith("ACTION:") for x in lines) != 1:
        fail("last line must be the one exact action line")
    result = {}
    for index, heading in enumerate(headings):
        end = positions[index + 1] if index + 1 < len(positions) else len(lines) - 1
        body = lines[positions[index] + 1 : end]
        if not body or any(not entry.strip() for entry in body):
            fail(f"empty or blank content in {heading}")
        result[heading] = body
    return result


def developer_view(kind, body):
    """Require the brief human decision fields; their truth remains a judgment."""
    if kind not in ("plan", "result"):
        return
    event = body["GOAL" if kind == "plan" else "CHANGED"]
    if not any(re.fullmatch(r"Why it matters: \S.*", line) for line in event):
        fail("missing Why it matters in Developer View")
    decision = body["DECISION"]
    choices = "APPROVE|REVISE|STOP" if kind == "plan" else "INTEGRATE|REVISE|DISCARD"
    if not any(re.fullmatch(rf"Recommendation: (?:{choices}) because \S.*", line)
               for line in decision):
        fail("missing recommendation and reason in Developer View")
    if not any(re.fullmatch(r"Next: \S.*", line) for line in decision):
        fail("missing next action in Developer View")


def reading_matches(kind, path, evidence_id, gate, disposition=None):
    """Bind the screen's GATE to one verbatim reading of this exact packet.

    The reading states its scope once: integration, publication, or both.
    A reading scoped to both may serve the RESULT and the PUBLICATION screen
    of the same frozen packet. A changed packet, commit, base or reading
    byte needs a new reading bound to the new identity (a correction is a
    new directory); the earlier verdict is never inherited.
    """
    if not path or not path.is_file() or path.is_symlink():
        fail("verbatim outside reading file is required")
    data = path.read_bytes()
    try:
        reading = data.decode("utf-8")
    except UnicodeDecodeError:
        fail("outside reading is not UTF-8")
    identities = re.findall(rf"^Packet identity:\s*({SHA})\s*$", reading, re.M)
    verdicts = re.findall(r"^Verdict:\s*(.+?)\s*$", reading, re.M)
    if identities != [evidence_id] or len(verdicts) != 1 or not re.fullmatch(r"PROCEED(?: IF .+)?", verdicts[0]):
        fail("outside reading does not name this packet with a proceeding verdict")
    declared = re.findall(r"^Scope:(.*)$", reading, re.M)
    if len(declared) != 1:
        fail("outside reading needs exactly one Scope line: integration, publication, "
             "or integration and publication")
    scopes = [declared[0].strip()]
    if scopes[0] not in ("integration", "publication", "integration and publication"):
        fail(f"reading Scope value not recognized: {scopes[0]!r}")
    if kind == "result" and scopes[0] == "publication":
        fail("reading scope does not cover integration")
    if kind == "publication" and scopes[0] == "integration":
        fail("reading scope does not cover publication")
    expected = [f"Reading SHA256: {digest(data)}", f"Verdict: {verdicts[0]}"]
    status = None
    if verdicts[0] != "PROCEED":
        if not disposition or not disposition.is_file() or disposition.is_symlink():
            fail("conditional verdict requires a separate disposition")
        note = disposition.read_text(encoding="utf-8")
        condition = verdicts[0].removeprefix("PROCEED IF ")
        statuses = re.findall(r"^Status: (.+)$", note, re.M)
        if (not note.startswith(f"Condition: {condition}\n") or
                statuses not in (["SATISFIED"], ["WAIVED"], ["PENDING"]) or
                not re.search(r"^Evidence: \S.*$", note, re.M)):
            fail("disposition does not name and address the exact condition")
        status = statuses[0]
        if status == "WAIVED" and not re.search(
                r'^Developer authorization: "\S.*"$', note, re.M):
            fail("waiver requires the Developer's quoted authorization")
        if status == "PENDING" and not re.search(r"^Pending: Developer \S.*$", note, re.M):
            fail("pending disposition requires the Developer decision")
        expected.append(f"Disposition SHA256: {digest(disposition.read_bytes())}")
    elif disposition:
        fail("unconditioned verdict has no disposition")
    if gate != expected:
        fail("screen gate does not identify the outside reading")
    return status


OUTGOING_KEYS = ("repository", "remote", "remote_url", "base", "commit", "tree", "refspec")
COMMIT = r"[0-9a-f]{40}"
BRANCH = r"[A-Za-z0-9._/-]+"


def parse_outgoing(path):
    """Read OUTGOING.txt: one block per repository, naming the exact proposed
    outgoing commit, its tree and base, the remote and the branch refspec.

    The checker binds these values to the PUBLICATION screen;
    publication_precheck.py compares them with the live repository before
    a push. Neither establishes that the push is authorized.
    """
    if not path.is_file() or path.is_symlink() or not path.stat().st_size:
        fail("missing or empty publication evidence: OUTGOING.txt")
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        fail("outgoing record is not UTF-8")
    needs = "outgoing record needs repository, remote, remote_url, base, commit, tree and refspec"
    blocks = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ": " not in line:
            fail("malformed outgoing record line")
        key, value = (part.strip() for part in line.split(": ", 1))
        if key == "repository":
            blocks.append({})
        if not blocks:
            fail(needs)
        if key in blocks[-1]:
            fail(f"duplicate key in outgoing record: {key}")
        blocks[-1][key] = value
    if not blocks:
        fail(needs)
    names = set()
    for block in blocks:
        if any(not block.get(key) for key in OUTGOING_KEYS):
            fail(needs)
        if block["repository"] in names:
            fail("duplicate repository in outgoing record")
        names.add(block["repository"])
        if any(not re.fullmatch(COMMIT, block[key]) for key in ("base", "commit", "tree")):
            fail("invalid Git identity in outgoing record")
        match = re.fullmatch(rf"({COMMIT}):refs/heads/({BRANCH})", block["refspec"])
        if not match or match[1] != block["commit"]:
            fail("outgoing refspec must be <commit>:refs/heads/<branch>")
        block["branch"] = match[2]
    return blocks


def suite_waiver(path):
    """Form check only. When the suite was not run, SERVER_SUITE.txt must
    carry the Developer's own words as a quoted block under a line
    'Waiver (Developer ...):', so that a publication-only permission cannot
    stand in for the suite waiver; a suffix after NOT RUN or a byte-order
    mark before it does not evade the check. The form establishes no authorization:
    whether those words waive the suite for this batch and scope is judged
    by the Developer and the Outside Reviewer, not by this check.
    """
    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeDecodeError:
        fail("publication evidence is not UTF-8: SERVER_SUITE.txt")
    if not lines or not lines[0].strip().startswith("SERVER_SUITE: NOT RUN"):
        return
    for index, line in enumerate(lines):
        if re.fullmatch(r"Waiver \(Developer[^)]*\):\s*", line):
            quoted = 0
            for following in lines[index + 1:]:
                if following.startswith("> ") and following[2:].strip():
                    quoted += 1
                else:
                    break
            if quoted:
                return
    fail("suite not run without the Developer's quoted waiver")


def code_identity(root, checked, changed):
    path = root / "CODE_IDENTITY.txt"
    for required in (path, root / "RUN_LOG.txt", root / "CHANGE.diff"):
        if not required.is_file() or not required.stat().st_size:
            fail(f"missing or empty change evidence: {required.name}")
    report = root / "APPROVAL_REPORT.md"
    if not report.is_file() or report.is_symlink():
        fail("missing approval report")
    try:
        body = report.read_text(encoding="utf-8")
        patch = (root / "CHANGE.diff").read_text(encoding="utf-8")
    except UnicodeDecodeError:
        fail("approval report and change diff must be UTF-8")
    if ("## Exact tested changes" not in body or
            "## New tests" not in body or patch not in body):
        fail("approval report must contain the complete exact CHANGE.diff and new-test section")
    if not any(x.startswith("Full diff and exact new tests: ") and
               "APPROVAL_REPORT.md" in x for x in changed):
        fail("result must point to the complete approval report")
    try:
        entries = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        fail("code identity is not UTF-8")
    values = {}
    for line in entries:
        if ": " not in line:
            fail("malformed code identity")
        key, value = line.split(": ", 1)
        if key in values or not value.strip():
            fail("duplicate or empty code identity field")
        values[key] = value
    for key in ("base_commit", "tested_tree", "check", "outcome", "criterion_by",
                "test_written_by", "run_by"):
        if key not in values:
            fail(f"missing code identity field: {key}")
    if not re.fullmatch(HEX, values["base_commit"]) or not re.fullmatch(HEX, values["tested_tree"]):
        fail("invalid Git identity")
    outcome = values["outcome"]
    # A complete run count is as explicit as PASS; do not require a new
    # frozen packet merely to prefix an otherwise truthful outcome line.
    counts = [re.fullmatch(r"(?:(?:before|after|baseline|candidate) )?([1-9]\d*)/([1-9]\d*) exit 0",
                           part.strip(), re.I) for part in outcome.split(";")]
    counted_pass = all(match and match[1] == match[2] for match in counts)
    if not (re.match(r"^PASS(?:$|[ :])", outcome) or counted_pass):
        fail("result requires a passing recorded check")
    tested = [x for x in checked if x.startswith("Tested tree: ")]
    if tested != [f'Tested tree: {values["tested_tree"]}']:
        fail("screen tested tree does not match the evidence")


def auto_screen(kind, body, queue):
    """Auto screens: the start notice states its limits; the report repeats the helper's rows."""
    if kind == "auto_start":
        if queue:
            fail("auto_start has no queue to compare")
        if f"python3 {AUTO_QUEUE} stop" not in body["STOP"] and f"python3 -I -B {AUTO_QUEUE} stop" not in body["STOP"]:
            fail("auto_start needs the stop command of this procedure's auto_queue.py")
        limits = " ".join(body["NOT AUTHORIZED"]).lower()
        if "master" not in limits or "publish" not in limits:
            fail("auto_start must say master is unchanged and nothing is published")
        text = " ".join(line for lines in body.values() for line in lines)
        if ("usage limit" not in text.lower() or not re.search(r"\bmanual(?:ly)?\b", text, re.IGNORECASE)
                or "/gc auto resume" not in text):
            fail("auto_start must say a usage limit stops the queue and resume is manual")
        return
    if queue is None:
        fail("auto_report requires --queue to compare with the queue record")
    helper = subprocess.run([sys.executable, "-I", "-B", str(AUTO_QUEUE), "--queue", str(queue),
                             "report"], capture_output=True, text=True)
    if helper.returncode != 0:
        fail("auto_queue.py report failed")
    expected = [line for line in helper.stdout.splitlines() if line.strip()]
    if body["JOBS"] != expected:
        fail("JOBS does not repeat the queue report exactly")
    words = " ".join(body["IN SIMPLE WORDS"] + body["NEXT"])
    if re.search(r"\bstopped\b", words, re.IGNORECASE) and expected[0] != "**Stopped**":
        fail("says stopped but the queue report does not")
    if re.search(r"\bpassed\b", words, re.IGNORECASE) and not any("passed" in row for row in expected):
        fail("says passed but no job row does")
    approval_summaries(body["APPROVAL SUMMARIES"], ready_jobs(expected))


SUMMARY_LABELS = ("Bug: ", "Fix: ", "Test: ", "Criterion: ", "Limits: ", "Approval: ")


def ready_jobs(rows):
    """Job ids of report rows whose step offers approval (or a provisional decision first)."""
    ready = []
    for row in rows:
        if row in ("| Job | Result | Decision or next step |", "| --- | --- | --- |"):
            continue
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        if len(cells) == 3 and (
                cells[2].startswith("approve for integration") or cells[2].startswith("decide P")):
            ready.append(cells[0].split()[0])
    return ready


def approval_summaries(section, ready):
    """One plain-words summary block per ready job, in report order; 'None ready.' otherwise."""
    if not ready:
        if section != ["None ready."]:
            fail("APPROVAL SUMMARIES must be 'None ready.' when no job is ready")
        return
    if len(section) != 7 * len(ready):
        fail("APPROVAL SUMMARIES needs exactly one seven-line block per ready job")
    for index, job in enumerate(ready):
        block = section[7 * index : 7 * index + 7]
        if not re.fullmatch(rf"Job {re.escape(job)}:.*", block[0]):
            fail(f"APPROVAL SUMMARIES block {index + 1} must start 'Job {job}:'")
        for label, line in zip(SUMMARY_LABELS, block[1:]):
            if not line.startswith(label) or not line[len(label):].strip():
                fail(f"summary for job {job} needs a non-empty '{label.strip()}' line in order")
        approval = block[6].lower()
        if "exact packet" not in approval or "merges nothing" not in approval:
            fail(f"summary for job {job}: Approval must say it accepts the exact packet and merges nothing")


def check(kind, screen, root=None, reading=None, disposition=None, report=True, queue=None):
    if not screen.is_file() or screen.is_symlink():
        fail("screen must be an ordinary file")
    try:
        # Bytes, so that a CR is seen rather than translated by universal newlines.
        source = screen.read_bytes().decode("utf-8")
    except UnicodeDecodeError:
        fail("screen is not UTF-8")
    lines = source.splitlines()
    max_lines, max_words = BUDGET.get(kind, (28, 220))
    if len(lines) > max_lines or len(source.split()) > max_words:
        fail("screen exceeds one-page budget")
    if not lines or not source.endswith("\n") or "\r" in source:
        fail("screen needs LF lines and a final newline")
    if re.search(r"<[^>]+>|\{\{[^}]+\}\}", source):
        fail("unfilled placeholder")
    parts = lines[0].split(" | ")
    if kind in ("plan", "result"):
        expected = 3 if kind == "plan" else 4
        if len(parts) != expected or parts[0] != kind.upper() or not re.fullmatch(NAME, parts[1]) or parts[2] not in ("light", "full"):
            fail("invalid plan/result header")
    elif kind in ATTENTION:
        if len(parts) != 2 or parts[0] != kind.upper().replace("_", " ") or not re.fullmatch(NAME, parts[1]):
            fail("invalid auto header")
    elif kind == "publication":
        if len(parts) != 3 or parts[0] != "PUBLICATION" or not re.fullmatch(NAME, parts[1]):
            fail("invalid publication header")
    elif len(parts) != 2 or parts[0] != "STOP" or not re.fullmatch(NAME, parts[1]):
        fail("invalid stop header")
    body = blocks(kind, lines)
    developer_view(kind, body)
    if kind in ATTENTION:
        if root or reading or disposition:
            fail("auto screens have no frozen evidence or outside reading")
        auto_screen(kind, body, queue)
    elif queue:
        fail("only auto_report compares with a queue")
    if kind in ("result", "publication"):
        if root is None or not re.fullmatch(SHA, parts[-1]):
            fail("evidence directory and manifest identity required")
        identity = verify(root)
        if parts[-1] != identity:
            fail("screen does not identify the current evidence")
        if kind == "result":
            code_identity(root, body["CHECKED"], body["CHANGED"])
            if parts[2] == "light":
                if reading or disposition or body["GATE"] != ["Reading: deferred to publication", "Verdict: deferred to publication"]:
                    fail("light-path outside reading belongs to publication batch")
            else:
                status = reading_matches("result", reading, identity, body["GATE"], disposition)
                if status in ("WAIVED", "PENDING") and not any(
                        status in line for heading in ("LIMITS", "DECISION")
                        for line in body[heading]):
                    fail("waived or pending condition must be visible in Developer View")
        else:
            for filename in ("SERVER_SUITE.txt", "ROSTER_COMPARISON.txt"):
                item = root / filename
                if not item.is_file() or not item.stat().st_size:
                    fail(f"missing or empty publication evidence: {filename}")
            batch = "\n".join(body["BATCH"])
            for block in parse_outgoing(root / "OUTGOING.txt"):
                if not re.search(rf"(?<![0-9A-Za-z]){block['commit']}(?![0-9A-Za-z])", batch):
                    fail("outgoing commit is not named in BATCH")
            suite_waiver(root / "SERVER_SUITE.txt")
            status = reading_matches("publication", reading, identity, body["GATE"], disposition)
            if status in ("WAIVED", "PENDING") and not any(
                    status in line for heading in ("LIMITS", "DECISION")
                    for line in body[heading]):
                fail("waived or pending condition must be visible in Developer View")
    elif kind not in ATTENTION and (root or reading or disposition):
        fail("plan and stop have no frozen evidence or outside reading")
    if report:
        print(f"CHECK OK {kind}")


def present(kind, screen, root=None, reading=None, disposition=None, queue=None):
    """Check the saved decision, then render a conspicuous Developer view."""
    check(kind, screen, root, reading, disposition, report=False, queue=queue)
    lines = screen.read_text(encoding="utf-8").splitlines()
    if kind == "plan" or kind in ATTENTION:
        spaced = []
        for line in lines:
            if line in HEADINGS[kind]:
                spaced.append("")
            spaced.append(line)
        lines = spaced
    if kind in ("result", "publication"):
        lines[0] = " | ".join(lines[0].split(" | ")[:-1])
    lines = [line for line in lines if not line.startswith(
        ("Tested tree: ", "Reading SHA256: ", "Disposition SHA256: "))]
    lines = [DISPLAY_ID.sub(
        lambda match: "[recorded]" if any(c in "abcdef" for c in match.group())
        else match.group(), line) for line in lines]
    body = "\n".join(lines)
    # A longer outer fence keeps any backticks in a screen inside its body.
    longest = max((len(match.group()) for match in re.finditer(r"`+", body)), default=0)
    fence = "`" * max(3, longest + 1)
    print(ATTENTION.get(kind, "**PLEASE READ — YOUR DECISION IS NEEDED**") + "\n\n" +
          fence + "\n" + body + "\n" + fence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("freeze", "verify"):
        commands.add_parser(name).add_argument("evidence", type=Path)
    commands.add_parser("verify-source").add_argument("source", type=Path)
    commands.add_parser("check-claude-version")
    commands.add_parser("register-skill").add_argument("source", type=Path)
    command = commands.add_parser("check")
    presenter = commands.add_parser("present")
    for command in (command, presenter):
        command.add_argument("kind", choices=tuple(HEADINGS))
        command.add_argument("screen", type=Path)
        command.add_argument("--evidence", type=Path)
        command.add_argument("--reading", type=Path)
        command.add_argument("--disposition", type=Path)
        command.add_argument("--queue", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "freeze":
            freeze(args.evidence)
        elif args.command == "verify":
            verify(args.evidence)
            print("VERIFIED evidence")
        elif args.command == "verify-source":
            verify_source(args.source)
            print("VERIFIED complete source")
        elif args.command == "check-claude-version":
            check_claude_version()
        elif args.command == "register-skill":
            register_skill(args.source)
        elif args.command == "check":
            check(args.kind, args.screen, args.evidence, args.reading, args.disposition,
                  queue=args.queue)
        else:
            present(args.kind, args.screen, args.evidence, args.reading, args.disposition,
                    queue=args.queue)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
