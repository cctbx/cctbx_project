#!/usr/bin/env python3
"""Check screen form, frozen bytes, source identity, and local registration."""

import sys
# Revision 18: the kit needs Python 3.10 or newer; an older Python refuses here, before
# anything else runs, with exit 19 (the same message and code in every tool).
if sys.version_info < (3, 10):
    sys.stderr.write("REFUSED: requires Python 3.10 or newer; this is Python %d.%d.%d (%s)\n"
                     % (tuple(sys.version_info[:3]) + (sys.executable,)))
    sys.exit(19)
# When run as a script, keep this directory out of module resolution before
# importing the standard library. The documented -I invocation also ignores
# PYTHONPATH and the current directory.
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
import hashlib
import importlib.util
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

HEADINGS = {
    "plan": ("GOAL", "SCOPE", "CHECK", "LIMITS", "DECISION"),
    "result": ("CHANGED", "CHECKED", "LIMITS", "GATE", "DECISION"),
    "accept": ("CANDIDATE", "CHECKED", "PENDING", "LIMITS", "DECISION"),
    "stop": ("WHY", "PRESERVED", "NEEDED"),
    "publication": ("BATCH", "SERVER CHECK", "LIMITS", "GATE", "DECISION"),
    "auto_start": ("JOBS", "CHECKS", "INSTALLATION", "RESULTS", "STOP", "NOT AUTHORIZED"),
    "auto_report": ("JOBS", "IN SIMPLE WORDS", "APPROVAL SUMMARIES", "NEXT"),
}
ACTIONS = {
    "plan": "ACTION: APPROVE / REVISE / STOP",
    "result": "ACTION: INTEGRATE / REVISE / DISCARD",
    "accept": "ACTION: ACCEPT / REVISE / DISCARD",
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
    """Require the brief human decision fields; their truth remains a judgment.

    Revision 17 (the reporting rule): the STOP and PUBLICATION screens also carry
    a `Recommendation: ... because ...` line and a `Next: ...` line, in NEEDED and
    in DECISION. A STOP recommends one of the options it names in NEEDED (free
    text, which may be the smallest useful check); a PUBLICATION recommends
    PUBLISH or HOLD."""
    if kind in ("stop", "publication"):
        section = body["NEEDED" if kind == "stop" else "DECISION"]
        choice = r"\S.*?" if kind == "stop" else r"(?:PUBLISH|HOLD)"
        if not any(re.fullmatch(rf"Recommendation: {choice} because \S.*", line)
                   for line in section):
            fail("missing recommendation and reason in Developer View")
        if not any(re.fullmatch(r"Next: \S.*", line) for line in section):
            fail("missing next action in Developer View")
        return
    if kind not in ("plan", "result", "accept"):
        return
    event = body[{"plan": "GOAL", "result": "CHANGED", "accept": "CANDIDATE"}[kind]]
    if not any(re.fullmatch(r"Why it matters: \S.*", line) for line in event):
        fail("missing Why it matters in Developer View")
    decision = body["DECISION"]
    choices = {"plan": "APPROVE|REVISE|STOP", "result": "INTEGRATE|REVISE|DISCARD",
               "accept": "ACCEPT|REVISE|DISCARD"}[kind]
    if not any(re.fullmatch(rf"Recommendation: (?:{choices}) because \S.*", line)
               for line in decision):
        fail("missing recommendation and reason in Developer View")
    if not any(re.fullmatch(r"Next: \S.*", line) for line in decision):
        fail("missing next action in Developer View")


def reading_matches(kind, path, evidence_id, gate, disposition=None, carried=None):
    """Bind the screen's GATE to one verbatim reading of this exact packet.

    The reading states its scope: integration, publication, or both.
    A reading scoped to both may serve the RESULT and the PUBLICATION screen
    of the same frozen packet. A changed packet, commit, base or reading
    byte needs a new reading, or a short addendum, bound to the new
    identity; the earlier verdict is never inherited.

    `Packet identity:`, `Verdict:` and `Scope:` lines that repeat the same
    value (for example a closing block repeated in a summary) count as one;
    repeats with different values are refused. Checkpoint 3: so are
    such lines written indented, bulleted, quoted or in bold (for example
    `- **Verdict:** DO NOT PROCEED`) whose value differs from the column-0
    value; before, only column-0 lines were read and these were ignored.
    The reading file is only read, and the GATE names the SHA-256 of its
    exact bytes.

    Revision 19 (item 2.2): an indented, bulleted, quoted, bold or table copy
    whose value is a pass word (`| Packet identity | Matches |`, `Verified`,
    `Pass`, `Passed`, ...) or the plain value followed by one (`- Packet
    identity: `<id>` — matches.`) is a reviewer's integrity note, not a
    conflicting value (passing_note). Conflicting values still block, and so
    does a second column-0 `Scope:` line with another value, even a prose
    one such as `Scope: Restage of ...`: the checker cannot tell a
    description from a narrower or wider scope.

    Revision 19 (carry-over): `carried` is carryover_evidence() of the packet:
    the identity of the earlier packet that the outside reader read, from
    CARRYOVER.txt and the packet's evidence/reviewed/ copies, and the mapping of
    the reviewed identities to the outgoing ones. A reading that names it
    instead of this packet's identity is accepted for publication only; check()
    then also requires the LIMITS line (CARRYOVER_LIMITS). A carried `PROCEED
    IF` keeps its condition: the disposition is still required, addresses the
    condition as worded and states the mapping of every reviewed identity the
    condition names (carried_condition).
    """
    if not path or not path.is_file() or path.is_symlink():
        fail("verbatim outside reading file is required")
    data = path.read_bytes()
    try:
        reading = data.decode("utf-8")
    except UnicodeDecodeError:
        fail("outside reading is not UTF-8")
    identities = set(re.findall(rf"^Packet identity:\s*(\S.*?)\s*$", reading, re.M))
    verdicts = sorted(set(re.findall(r"^Verdict:\s*(.+?)\s*$", reading, re.M)))
    named = evidence_id
    if carried is not None and kind == "publication" and identities == {carried["reviewed_id"]}:
        named = carried["reviewed_id"]
    if identities != {named} or len(verdicts) != 1 or not re.fullmatch(r"PROCEED(?: IF .+)?", verdicts[0]):
        if carried is None and len(identities) == 1 and kind == "publication":
            fail("outside reading does not name this packet with a proceeding verdict (a "
                 "reading of an earlier packet carries over only with CARRYOVER.txt and "
                 "evidence/reviewed/ in this packet)")
        fail("outside reading does not name this packet with a proceeding verdict")
    declared = set(value.strip() for value in re.findall(r"^Scope:(.*)$", reading, re.M))
    if len(declared) != 1:
        fail("outside reading needs exactly one Scope line (identical repeats count as "
             "one): integration, publication, or integration and publication")
    scopes = [declared.pop()]
    plain = {"Packet identity": {named}, "Verdict": set(verdicts), "Scope": set(scopes)}
    for label, value in decorated_fields(reading):
        if value not in plain[label] and not passing_note(value, plain[label]):
            fail(f"outside reading has an indented, bulleted, quoted or bold {label} line "
                 f"that differs from its {label} line: {value[:60]!r}")
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
        if named != evidence_id:
            carried_condition(condition, note, carried, evidence_id)
        expected.append(f"Disposition SHA256: {digest(disposition.read_bytes())}")
    elif disposition:
        fail("unconditioned verdict has no disposition")
    if gate != expected:
        fail("screen gate does not identify the outside reading")
    return status


def reading_identities(path):
    """The column-0 Packet identity values of a reading (as reading_matches reads them)."""
    return set(re.findall(r"^Packet identity:\s*(\S.*?)\s*$",
                          path.read_bytes().decode("utf-8"), re.M))


WAIVER_GATE = "Reading: WAIVED by the Developer"


def waiver_matches(path, evidence_id, gate, root):
    """Revision 17: the Developer's waiver of the outside reading for one full-path
    integration, in place of a reading (`--waiver`, never with `--reading`).

    OUTSIDE_WAIVER.md (WORKER.md) lies outside the frozen packet and has the
    lines `Outside reading: WAIVED`, `Packet identity: <sha256 of this packet's
    MANIFEST.sha256>`, `Scope: integration` and `Developer authorization:`,
    followed at once by one or more `> ` lines quoting the Developer's own words.
    Each of the first three appears once (identical repeats count as one; an
    indented, bulleted, quoted or bold copy with another value is refused, as for
    a reading). The screen's GATE is then exactly `Reading: WAIVED by the
    Developer` and `Waiver SHA256: <sha256 of the waiver file>`; check() also
    requires WAIVED in LIMITS. A waiver covers only the packet identity it names
    and never publication. The form establishes no authorization: whether the
    quoted words waive the reading for this change is judged by the Developer
    and the Outside Reviewer, not by this check. The waiver file is only read."""
    if not path or not path.is_file() or path.is_symlink():
        fail("waiver file is required and must be an ordinary file")
    inside, packet = path.resolve(), root.resolve()
    if packet in inside.parents:
        fail("waiver must be outside the frozen packet (in the companions directory)")
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        fail("waiver is not UTF-8")
    lines = text.splitlines()
    # Repair 1: fields are read within their own line ([ \t]*, not \s*, which crossed lines).
    if set(re.findall(r"^Outside reading:[ \t]*(\S.*?)[ \t]*$", text, re.M)) != {"WAIVED"}:
        fail("waiver needs the line 'Outside reading: WAIVED'")
    for line in lines:
        if line.startswith("Outside reading:"):
            continue
        bare = re.sub(rf"^(?:[\s>+#|{DECORATION}-]|\d+[.)])*", "", line)
        match = re.match(rf"Outside reading[{DECORATION}]*\s*[:|][{DECORATION}]*(.*)$", bare)
        if match and match[1].strip().rstrip("|").strip().strip(DECORATION).strip() != "WAIVED":
            fail(f"waiver has an indented, bulleted, quoted or bold Outside reading line that "
                 f"differs from 'WAIVED': {line[:60]!r}")
    if any(label == "Verdict" for label, _ in decorated_fields(text)) or re.search(
            r"^Verdict:", text, re.M):
        fail("a waiver has no Verdict line (it replaces the reading; it is not one)")
    identities = set(re.findall(r"^Packet identity:[ \t]*(\S.*?)[ \t]*$", text, re.M))
    if identities != {evidence_id}:
        fail("waiver does not name this packet identity")
    scopes = set(value.strip() for value in re.findall(r"^Scope:(.*)$", text, re.M))
    if scopes != {"integration"}:
        fail("waiver needs exactly one Scope line, 'Scope: integration' (a waiver never "
             "covers publication)")
    plain = {"Packet identity": identities, "Scope": scopes}
    for label, value in decorated_fields(text):
        if label in plain and value not in plain[label]:
            fail(f"waiver has an indented, bulleted, quoted or bold {label} line that "
                 f"differs from its {label} line: {value[:60]!r}")
    quoted = 0
    for index, line in enumerate(lines):
        if re.fullmatch(r"Developer authorization:\s*", line):
            for following in lines[index + 1:]:
                if following.startswith("> ") and following[2:].strip():
                    quoted += 1
                else:
                    break
            break
    if not quoted:
        fail("waiver needs 'Developer authorization:' followed by the Developer's own words "
             "as '> ' lines")
    if gate != [WAIVER_GATE, f"Waiver SHA256: {digest(data)}"]:
        fail("screen gate does not identify the waiver: it must be "
             f"'{WAIVER_GATE}' and 'Waiver SHA256: <sha256 of the waiver file>'")


DECORATION = "*_`\"'"
# Revision 19 (item 2.2): words with which a reviewer reports that a check passed.
PASS_WORDS = r"(?:verified|matches|match|matched|pass|passed|ok|confirmed|identical|correct)"


def passing_note(value, plain):
    """Revision 19 (item 2.2): True when a decorated field's value is a pass word alone
    (`Matches`, `Verified`, `Pass`, `Passed`, ...) or one of the plain values followed by a
    pass word (`<id> — matches.`), with decoration and closing punctuation ignored. Such a
    row reports a check, it does not state another value."""
    bare = re.sub(rf"[{re.escape(DECORATION)}]", "", value).strip()
    if re.fullmatch(rf"{PASS_WORDS}[.!]?", bare, re.I):
        return True
    return any(re.fullmatch(rf"{re.escape(item)}\s*(?:[-–—:,;(]\s*)?{PASS_WORDS}[.!)]*", bare, re.I)
               for item in plain)


def decorated_fields(reading):
    """(label, value) of every Packet identity, Verdict or Scope line that is
    not in the plain column-0 form: indented, bulleted (-, *, +, a number),
    quoted (>, quote marks) or bold/italic/code, with that decoration removed.
    Checkpoint 5 (review CP3 note 8): also a Markdown heading (`## Verdict:
    ...`) and a table row (`| Verdict | ... |`, the cell separator standing
    for the colon); a table header row naming these labels fails closed."""
    found = []
    for line in reading.splitlines():
        if re.match(r"(?:Packet identity|Verdict|Scope):", line):
            continue
        bare = re.sub(rf"^(?:[\s>+#|{DECORATION}-]|\d+[.)])*", "", line)
        match = re.match(rf"(Packet identity|Verdict|Scope)[{DECORATION}]*\s*[:|][{DECORATION}]*(.*)$",
                         bare)
        if match:
            value = match[2].strip()
            if line.lstrip().startswith("|"):
                value = value.rstrip("|").strip()
            found.append((match[1], value.strip(DECORATION).strip()))
    return found


OUTGOING_KEYS = ("repository", "remote", "remote_url", "base", "commit", "tree", "refspec")
COMMIT = r"[0-9a-f]{40}"
BRANCH = r"[A-Za-z0-9._/-]+"


def parse_outgoing(path):
    """Read OUTGOING.txt: one block per repository, naming the exact proposed
    outgoing commit, its tree and base, the remote and the branch refspec.

    The checker binds these values to the PUBLICATION screen;
    publication_precheck.py compares them with the live repository before
    a push. Neither establishes that the push is authorized.

    Checkpoint 3: the same strict line rule as check-draft (column-0
    lowercase `key: value` lines of the seven keys only; blank lines and
    `#` comments ignored), so the draft and the final check agree. A packet
    frozen earlier with prose or another key in OUTGOING.txt is refused and
    needs a corrected packet. `publication_precheck.py outgoing` writes a
    block in this form.
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
        # Checkpoint 3: one rule for the draft and the final check (replaces the
        # final check's looser `key: value` split): column-0 lowercase key lines
        # only, and only the seven keys, so prose is refused even with ': '.
        if not re.fullmatch(r"[a-z][a-z0-9_]*: \S.*", line):
            fail(f"OUTGOING.txt has a prose or non-key line (put prose in "
                 f"PUBLICATION_NOTES.md): {line[:60]!r}")
        key, value = (part.strip() for part in line.split(": ", 1))
        if key not in OUTGOING_KEYS:
            fail(f"OUTGOING.txt key is not one of the seven outgoing keys (put descriptions "
                 f"in PUBLICATION_NOTES.md): {key}")
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

    Checkpoint 3: the first line (after a byte-order mark and surrounding
    spaces) must be `SERVER_SUITE: <outcome>`, as PUBLICATION_PACKET.md
    says, in the draft and the final check; before, any other first line
    passed unchecked. What the outcome says is not judged here.

    Revision 19 (plan revision 5): the server suite is never carried over to a
    restaged commit. A first line `SERVER_SUITE: CARRIED OVER ...` is refused: run
    the suite, or quote the Developer's waiver; CARRYOVER.txt only informs that
    decision.
    """
    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeDecodeError:
        fail("publication evidence is not UTF-8: SERVER_SUITE.txt")
    if not lines or not re.match(r"SERVER_SUITE: \S", lines[0].strip()):
        fail("SERVER_SUITE.txt first line must be 'SERVER_SUITE: <outcome>'")
    if re.match(r"SERVER_SUITE: CARRIED\s+OVER\b", lines[0].strip(), re.I):
        fail("the server suite is never carried over to a restaged commit: run it, or quote "
             "the Developer's waiver (SERVER_SUITE: NOT RUN)")
    if not lines[0].strip().startswith("SERVER_SUITE: NOT RUN"):
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


# Revision 19: the reading carry-over after a restage (publication_precheck.py carry-over).
CARRYOVER_HEADER = "CARRYOVER (publication_precheck.py carry-over, revision 19)"
CARRYOVER_OK = "CARRYOVER OK reading"
REVIEWED_MANIFEST = "evidence/reviewed/MANIFEST.sha256"
REVIEWED_OUTGOING = "evidence/reviewed/OUTGOING.txt"
# The LIMITS line a carried-over reading needs (the Developer's condition, 2026-10-18), and
# the key words the checker enforces on the one line that contains "reading carried over".
CARRYOVER_LIMITS = ("Reading carried over: the check shows the changed files are identical and "
                    "untouched upstream; it does not show that upstream changes elsewhere "
                    "cannot interact with them.")
CARRYOVER_WORDS = ("reading carried over", "identical", "untouched upstream", "does not show",
                   "upstream changes elsewhere", "cannot interact")
# Upstream files CARRYOVER.txt flags for the Developer's server-suite decision (information
# only; the suite is never carried over).
FLAG_KINDS = ("compiled", "test", "configuration", "data", "checkout")
FLAGGED_INTRO = ("flagged upstream files (information for the server-suite decision; the suite "
                 "is not carried over): ")
# Plan revision 6 (a): upstream files a changed Python file imports or names by a literal path.
IMPORTED_LABEL = "upstream files imported by a changed file"
NAMED_LABEL = "upstream files named by a changed file"
COMPILED_SUFFIXES = (".c", ".cc", ".cpp", ".cxx", ".h", ".hpp", ".ipp", ".pyx", ".cu", ".so")
# Repair 1 (review finding 7): files that change how the identical files are checked out or
# archived (attributes, submodules).
CHECKOUT_NAMES = (".gitattributes", ".gitmodules")
COMPILED_NAMES = ("SConscript", "SConstruct", "libtbx_config")
CONFIG_SUFFIXES = (".cfg", ".ini", ".toml", ".yaml", ".yml", ".conf", ".phil")
CONFIG_NAMES = ("SConscript", "SConstruct", "libtbx_config", "libtbx_refresh.py", "setup.py",
                "setup.cfg", "pyproject.toml", "CMakeLists.txt", "Makefile", "requirements.txt",
                "environment.yml", "conftest.py", "run_tests.py")
DATA_SUFFIXES = (".pdb", ".cif", ".mmcif", ".ent", ".mtz", ".hkl", ".sca", ".ccp4", ".mrc",
                 ".map", ".dat", ".csv", ".tsv", ".npy", ".npz", ".pkl", ".pickle", ".gz",
                 ".bz2", ".zip", ".tar", ".json", ".xml", ".seq", ".fa", ".fasta", ".pir")


def upstream_flags(path):
    """The kinds (FLAG_KINDS order) of an upstream file: compiled (C/C++ sources and headers,
    SCons files, libtbx_config), test (tst_*, test_*, *_test.py, or under tests/, test/ or
    regression/), configuration (build, setup and parameter files) and data (structure,
    reflection, map and other data files, or under data/)."""
    parts = path.split("/")
    name = parts[-1]
    lowered = name.lower()
    kinds = []
    if lowered.endswith(COMPILED_SUFFIXES) or name in COMPILED_NAMES:
        kinds.append("compiled")
    if (re.match(r"(?:tst|test)_", name) or lowered.endswith("_test.py")
            or any(part in ("tests", "test", "regression") for part in parts[:-1])):
        kinds.append("test")
    if lowered.endswith(CONFIG_SUFFIXES) or name in CONFIG_NAMES:
        kinds.append("configuration")
    if lowered.endswith(DATA_SUFFIXES) or "data" in parts[:-1]:
        kinds.append("data")
    if name in CHECKOUT_NAMES:
        kinds.append("checkout")
    return tuple(kinds)


def carryover_section(lines, name):
    """The lines of one repository's section of CARRYOVER.txt, from `== repository: NAME`
    through its `repository NAME: ...` result line (empty when absent)."""
    try:
        first = lines.index(f"== repository: {name}")
    except ValueError:
        return []
    for index in range(first + 1, len(lines)):
        if lines[index].startswith(f"repository {name}: "):
            return lines[first:index + 1]
        if lines[index].startswith("== repository: ") or lines[index].startswith("CARRYOVER "):
            break
    return []


# Repair 2 (review 2, findings 1 and 5): the content differences, one shared list used by
# publication_precheck.py to write them and by this checker to decide what may be waived.
DIFF_BLOBS = "changed file content differs: "
DIFF_DIFF = "git diff of the outgoing commit is not byte-identical to the reviewed one"
DIFF_TOUCH = "upstream commits touch changed paths: "
DIFF_MESSAGE = "the commit message differs from the reviewed one"
WAIVABLE_REFUSALS = (DIFF_BLOBS, DIFF_DIFF, DIFF_TOUCH, DIFF_MESSAGE)


def content_result(count):
    """The result of a repository section whose only refusals are content differences."""
    return f"{count} content difference(s)"


def consistent_section(section, name):
    """Revision 19: rebuild one repository's section from the values it states and refuse any
    difference, so a hand edit that keeps the identities and the final line (a blob marked
    DIFFERS or same against its values, a count that does not match its rows, an upstream row
    on a changed path not reported, wrong flags, differences that the rows do not show) is
    refused. Repair 2: the non-waivable facts (ancestry, the same changed paths) come first;
    then every content comparison is recorded, and a section with content differences ends
    with one `difference (NAME): ...` line per difference and `repository NAME: REFUSED: N
    content difference(s)`. Facts that only Git can establish are rechecked by
    `publication_precheck.py carry-over --expect` and `check --packet` at the PUBLISH precheck.
    Returns ({imported, named}, differences)."""
    rows = iter(section[3:])
    rebuilt = list(section[:3])

    def take(pattern, what):
        line = next(rows, None)
        match = re.fullmatch(pattern, line or "")
        if not match:
            fail(f"CARRYOVER.txt is not internally consistent ({name}: {what}): {(line or '')[:80]!r}")
        return line, match
    rebuilt.append(take(r"reviewed base is an ancestor of the outgoing base: yes", "ancestor")[0])
    line, match = take(r"changed paths: ([1-9]\d*) \(the same set in both commits\)", "changed paths")
    rebuilt.append(line)
    changed, python, upstream, differing, order = set(), set(), set(), [], []
    for _ in range(int(match[1])):
        line, blob = take(r"blob (same|DIFFERS) reviewed (deleted|\d{6} [0-9a-f]{40}) outgoing "
                          r"(deleted|\d{6} [0-9a-f]{40}) (.+)", "blob rows")
        if (blob[1] == "same") != (blob[2] == blob[3]) or blob[4] in changed:
            fail(f"CARRYOVER.txt is not internally consistent ({name}: blob row {blob[4]!r})")
        changed.add(blob[4])
        order.append(blob[4])
        if blob[1] == "DIFFERS":
            differing.append(blob[4])
        if blob[4].endswith(".py") and blob[3] != "deleted":
            python.add(blob[4])
        rebuilt.append(line)
    if order != sorted(order):
        fail(f"CARRYOVER.txt is not internally consistent ({name}: blob rows out of order)")
    differences = [DIFF_BLOBS + " ".join(differing)] if differing else []
    line, match = take(r"diff sha256 reviewed ([0-9a-f]{64}) outgoing ([0-9a-f]{64}) \(\d+ bytes\)",
                       "diff hashes")
    rebuilt.append(line)
    if match[1] == match[2]:
        rebuilt.append(take(r"diff: byte-identical", "diff")[0])
    else:
        rebuilt.append(take(r"diff: DIFFERS", "diff")[0])
        differences.append(DIFF_DIFF)
    line, match = take(r"upstream commits: (\d+)", "upstream count")
    rebuilt.append(line)
    flagged = {kind: set() for kind in FLAG_KINDS}
    touching = []
    for _ in range(int(match[1])):
        line, commit = take(r"upstream commit ([0-9a-f]{40}): (\d+) file\(s\)", "upstream commits")
        rebuilt.append(line)
        for _ in range(int(commit[2])):
            line, row = take(r"  ([A-Z][0-9]*) (.+?)(?: \[([a-z,]+)\])?", "upstream files")
            path = row[2]
            flags = upstream_flags(path)
            if tuple((row[3] or "").split(",") if row[3] else ()) != flags:
                fail(f"CARRYOVER.txt is not internally consistent ({name}: upstream file {path!r})")
            if path in changed and f"{commit[1][:12]} {path}" not in touching:
                touching.append(f"{commit[1][:12]} {path}")
            for kind in flags:
                flagged[kind].add(path)
            upstream.add(path)
            rebuilt.append(line)
    if touching:
        rebuilt.append(take(re.escape("upstream commits touching a changed path: " + "; ".join(touching)),
                            "touching")[0])
        differences.append(DIFF_TOUCH + "; ".join(touching))
    else:
        rebuilt.append(take(r"upstream commits touching a changed path: none", "touching")[0])
    rebuilt.append(FLAGGED_INTRO + ", ".join(f"{kind} {len(flagged[kind])}" for kind in FLAG_KINDS))
    rebuilt += [f"flagged {kind}: {' '.join(sorted(flagged[kind])) or 'none'}" for kind in FLAG_KINDS]
    for _ in FLAG_KINDS:
        next(rows, None)
    next(rows, None)
    line, read = take(r"changed Python files read for imports: (\d+) \(not parsed: (.+)\)", "imports read")
    unparsed = [] if read[2] == "none" else read[2].split(" ")
    if int(read[1]) + len(unparsed) != len(python) or not set(unparsed) <= python:
        fail(f"CARRYOVER.txt is not internally consistent ({name}: changed Python files read)")
    rebuilt.append(line)
    flagged_imports = {}
    for label, verb in ((IMPORTED_LABEL, "imported"), (NAMED_LABEL, "named")):
        line, found = take(rf"{re.escape(label)}: (.+)", label)
        entries = [] if found[1] == "none" else found[1].split("; ")
        for entry in entries:
            item = re.fullmatch(rf"(\S+) \({verb} by (\S+)\)", entry)
            if not item or item[1] not in upstream or item[2] not in python or item[2] in unparsed:
                fail(f"CARRYOVER.txt is not internally consistent ({name}: {label}: {entry[:60]!r})")
        flagged_imports[verb] = [entry.split(" ")[0] for entry in entries]
        rebuilt.append(line)
    line, message = take(r"commit message sha256 reviewed ([0-9a-f]{64}) outgoing ([0-9a-f]{64})",
                         "commit message hashes")
    rebuilt.append(line)
    if message[1] == message[2]:
        rebuilt.append(take(r"commit message: identical", "commit message")[0])
    else:
        rebuilt.append(take(r"commit message: DIFFERS", "commit message")[0])
        differences.append(DIFF_MESSAGE)
    rebuilt += [f"difference ({name}): {text}" for text in differences]
    rebuilt.append(f"repository {name}: " + (f"REFUSED: {content_result(len(differences))}"
                                              if differences else "OK"))
    if rebuilt != section:
        fail(f"CARRYOVER.txt is not internally consistent ({name}: its flagged lists, its "
             "differences, its result line or extra lines)")
    return flagged_imports, differences


def waivable(lines):
    """True when CARRYOVER.txt refused only in repository sections whose result is `N
    content difference(s)` and every `difference (...)` line is one of WAIVABLE_REFUSALS: a
    missing or invalid commit, an ancestry failure, other changed paths or a repository
    mismatch compared no content and cannot be waived (repairs 1 and 2)."""
    refused = [line.split(": REFUSED: ", 1)[1] for line in lines
               if re.match(r"repository \S+: REFUSED: ", line)]
    found = [line.split("): ", 1)[1] for line in lines if re.match(r"difference \(\S+\): ", line)]
    return (bool(refused) and all(re.fullmatch(r"[1-9]\d* content difference\(s\)", reason)
                                  for reason in refused)
            and bool(found) and all(text.startswith(WAIVABLE_REFUSALS) for text in found))


def normal(text):
    return " ".join(text.split()).lower()


# Repairs 1-3 (reviews 1-3, finding 2; the Developer's choice A on STOP_1): one LIMITS line
# must contain the Developer's sentence itself, which is the real requirement. As a backstop,
# no LIMITS line of a carried or waived reading may make one of these explicit claims (any
# case): "cannot interact", "can not interact", "can't interact", "no interaction", "does not
# affect", "do not affect", "cannot affect", "will not affect". Honest limit lines such as
# "possible interaction ... is not excluded" or "may interact" pass; a paraphrase is the
# reviewer's matter.
INTERACTION_CLAIM = re.compile(
    r"\b(?:cannot|can not|can['’]t)\s+interact\w*|\bno\s+interaction\b"
    r"|\b(?:does|do|will)\s+not\s+affect\b|\bcannot\s+affect\b", re.I)


def interaction_claim(line):
    """True when a LIMITS line, once the Developer's sentence is removed, claims that nothing
    upstream can interact with or affect the change."""
    text = normal(line).replace(normal(CARRYOVER_LIMITS), " ")
    return bool(INTERACTION_CLAIM.search(text))


def carryover_limits_line(line):
    text = normal(line)
    sentence = normal(CARRYOVER_LIMITS)
    if sentence not in text:
        return False
    return not interaction_claim(line)


def ordinary(path, what):
    if path.is_symlink() or not path.is_file() or not path.stat().st_size:
        fail(f"carry-over evidence missing, empty or not an ordinary file: {what}")
    return path.read_bytes()


def carryover_evidence(root, top=None):
    """Revision 19: the carry-over evidence of a publication packet, or None when the packet
    has no CARRYOVER.txt. When present it must hold, or the check refuses:
    - CARRYOVER.txt in the fixed form `publication_precheck.py carry-over` writes, with no
      refusal, ending with the one line `CARRYOVER OK reading`, and every repository section
      internally consistent (consistent_section);
    - evidence/reviewed/MANIFEST.sha256, a copy of the MANIFEST.sha256 of the packet the
      outside reader read (its SHA-256 is that packet's identity), with an `OUTGOING.txt`
      entry, and evidence/reviewed/OUTGOING.txt whose SHA-256 is that entry;
    - CARRYOVER.txt's reviewed_outgoing_sha256 and outgoing_sha256 equal to the SHA-256 of
      evidence/reviewed/OUTGOING.txt and of this packet's OUTGOING.txt, and its reviewed: and
      outgoing: lines equal to each repository block of the two records, which name the same
      repositories, remotes, destinations and branches.
    The Git facts are established by the tool that wrote CARRYOVER.txt; this checker has no
    repository, so the PUBLISH precheck recomputes them (`carry-over --expect`, and `check`).
    Returns {reviewed_id, reviewed_commits, mapping}; mapping pairs each reviewed base, commit
    and tree with the outgoing one."""
    top = top or {}
    path = root / "CARRYOVER.txt"
    if not path.exists() and not path.is_symlink():
        # Repair 1 (review finding 1): a CARRYOVER.txt elsewhere in the packet, or the reviewed
        # copies without one, is refused rather than ignored.
        if any(found.name == "CARRYOVER.txt" for found in root.rglob("CARRYOVER.txt")) or (
                root / "evidence" / "reviewed").exists():
            fail("the packet holds a CARRYOVER.txt or evidence/reviewed/ but no CARRYOVER.txt at "
                 "its top level")
        return None
    try:
        lines = ordinary(path, "CARRYOVER.txt").decode("utf-8").splitlines()
    except UnicodeDecodeError:
        fail("CARRYOVER.txt is not UTF-8")
    if not lines or lines[0] != CARRYOVER_HEADER:
        fail("CARRYOVER.txt is not in the fixed form written by publication_precheck.py "
             "carry-over")
    refused = lines[-1].startswith("CARRYOVER REFUSED: ")
    finals = [line for line in lines[1:] if line.startswith("CARRYOVER ")]
    if refused:
        if finals != [lines[-1]] or CARRYOVER_OK in lines:
            fail("CARRYOVER.txt must end with exactly one final CARRYOVER line")
    else:
        if any(line.startswith("CARRYOVER REFUSED") or ": REFUSED: " in line for line in lines):
            fail("CARRYOVER.txt records a refusal: the reading does not carry over (a real "
                 "addendum is needed)")
        if lines[-1] != CARRYOVER_OK or lines.count(CARRYOVER_OK) != 1 or finals != [CARRYOVER_OK]:
            fail(f"CARRYOVER.txt must end with the one line '{CARRYOVER_OK}'")

    def field(name, index):
        if not re.fullmatch(rf"{name}: {SHA}", lines[index] if len(lines) > index else ""):
            fail(f"CARRYOVER.txt needs its {name} line in place")
        return lines[index][len(name) + 2:]
    reviewed_sha, outgoing_sha = field("reviewed_outgoing_sha256", 1), field("outgoing_sha256", 2)
    manifest = ordinary(root / REVIEWED_MANIFEST, REVIEWED_MANIFEST)
    reviewed_outgoing = ordinary(root / REVIEWED_OUTGOING, REVIEWED_OUTGOING)
    try:
        entries = [re.fullmatch(r"([0-9a-f]{64})  ([^\n\r]+)", line)
                   for line in manifest.decode("utf-8").splitlines()]
    except UnicodeDecodeError:
        fail(f"{REVIEWED_MANIFEST} is not UTF-8")
    if not entries or any(entry is None for entry in entries):
        fail(f"{REVIEWED_MANIFEST} is not a packet manifest")
    listed = [entry[1] for entry in entries if entry[2] == "OUTGOING.txt"]
    if len(listed) != 1:
        fail(f"{REVIEWED_MANIFEST} lists no single OUTGOING.txt")
    if digest(reviewed_outgoing) != listed[0]:
        fail(f"{REVIEWED_OUTGOING} is not the OUTGOING.txt that {REVIEWED_MANIFEST} lists")
    if reviewed_sha != digest(reviewed_outgoing):
        fail(f"CARRYOVER.txt names another reviewed OUTGOING.txt than {REVIEWED_OUTGOING}")
    outgoing_path = top.get("OUTGOING.txt", root / "OUTGOING.txt")
    if outgoing_sha != digest(ordinary(outgoing_path, "OUTGOING.txt")):
        fail("CARRYOVER.txt names another OUTGOING.txt than this packet's (a stale carry-over "
             "from an earlier restage?)")
    old = {block["repository"]: block for block in parse_outgoing(root / REVIEWED_OUTGOING)}
    new = {block["repository"]: block for block in parse_outgoing(outgoing_path)}
    if set(old) != set(new):
        fail("the reviewed and the outgoing OUTGOING.txt name different repositories")
    names = [line[len("== repository: "):] for line in lines if line.startswith("== repository: ")]
    if sorted(names) != sorted(new) and not (refused and not names):
        fail("CARRYOVER.txt does not cover exactly the repositories of both OUTGOING.txt files")
    covered = 3
    mapping, flagged, refusals = [], {"imported": [], "named": []}, []
    for name in names:
        before, after = old[name], new[name]
        for key in ("remote", "remote_url", "branch"):
            if before[key] != after[key]:
                fail(f"the reviewed and the outgoing {name} blocks name another {key}")
        section = carryover_section(lines, name)
        expected = [f"== repository: {name}",
                    f"reviewed: base {before['base']} commit {before['commit']} tree {before['tree']}",
                    f"outgoing: base {after['base']} commit {after['commit']} tree {after['tree']}"]
        if section[:3] != expected:
            fail(f"CARRYOVER.txt does not record the reviewed and outgoing {name} blocks of the "
                 "OUTGOING.txt files")
        result = section[-1][len(f"repository {name}: "):] if section else ""
        if result.startswith("REFUSED: ") and not refused:
            fail("CARRYOVER.txt records a refusal: the reading does not carry over")
        if result.startswith("REFUSED: ") and not re.fullmatch(
                r"REFUSED: [1-9]\d* content difference\(s\)", result):
            # A refusal that compared no content: the section stops where the tool stopped.
            refusals.append(f"{name}: " + result[len("REFUSED: "):])
        else:
            found, differences = consistent_section(section, name)
            for verb, paths in found.items():
                flagged[verb] += paths
            if differences:
                refusals.append(f"{name}: {content_result(len(differences))}")
        covered += len(section)
        mapping += [(before[key], after[key]) for key in ("base", "commit", "tree")]
    if covered + 1 != len(lines):
        fail("CARRYOVER.txt has lines outside its repository sections")
    if refused and names and lines[-1] != "CARRYOVER REFUSED: " + "; ".join(refusals):
        fail("CARRYOVER.txt's final line does not repeat the refusals of its sections")
    if refused and not names and len(lines) != 4:
        fail("CARRYOVER.txt has lines outside its repository sections")
    return {"reviewed_id": digest(manifest), "mapping": mapping, "refused": refused,
            "lines": lines, "flagged": flagged, "outgoing": new,
            "reviewed_commits": sorted(block["commit"] for block in old.values())}


ADDENDUM_WAIVER = "ADDENDUM_WAIVER.md"


def addendum_waiver(root, carry):
    """Plan revision 6 (b): when carry-over refused, the earlier reading is used only with the
    Developer's waiver of the addendum, in the packet as ADDENDUM_WAIVER.md, fixed form:

        Addendum: WAIVED
        Reviewed packet identity: <sha256 of evidence/reviewed/MANIFEST.sha256>
        Outgoing: <repository> base <40 hex> commit <40 hex> tree <40 hex>   (one per block)
        Differences (from CARRYOVER.txt):
        > <each 'difference (NAME): ...' line of CARRYOVER.txt, copied, and nothing else>
        Developer decision: "<his own words>"

    The identities must be those of this packet and of the read packet; the quoted lines must
    be exactly the content differences CARRYOVER.txt records (repair 2). Form only: whether the words waive the addendum for these
    differences is the Developer's decision, judged with the reviewer."""
    if not waivable(carry["lines"]):
        fail("this carry-over refusal cannot be waived: only a content difference (changed file "
             "content, the diff, an upstream commit touching a changed file, the commit message) "
             "can; a missing or invalid commit, an ancestry failure or a repository mismatch "
             "needs a new reading")
    path = root / ADDENDUM_WAIVER
    try:
        text = ordinary(path, ADDENDUM_WAIVER).decode("utf-8")
    except UnicodeDecodeError:
        fail(f"{ADDENDUM_WAIVER} is not UTF-8")
    lines = text.splitlines()

    def once(pattern, what):
        found = [line for line in lines if re.fullmatch(pattern, line)]
        if len(found) != 1:
            fail(f"{ADDENDUM_WAIVER} needs exactly one line {what}")
        return found[0]
    once(r"Addendum: WAIVED", "'Addendum: WAIVED'")
    if once(rf"Reviewed packet identity: {SHA}", "'Reviewed packet identity: <sha256>'")[26:] != carry["reviewed_id"]:
        fail(f"{ADDENDUM_WAIVER} names another reviewed packet identity than evidence/reviewed/")
    outgoing = sorted(line for line in lines if line.startswith("Outgoing: "))
    expected = sorted(f"Outgoing: {b['repository']} base {b['base']} commit {b['commit']} tree {b['tree']}"
                      for b in carry["outgoing"].values())
    if outgoing != expected:
        fail(f"{ADDENDUM_WAIVER} does not name this packet's outgoing identities (one 'Outgoing:' "
             "line per repository block)")
    once(r'Developer decision: "\S.*"', 'Developer decision: "<his own words>"')
    start = [index for index, line in enumerate(lines) if line == "Differences (from CARRYOVER.txt):"]
    if len(start) != 1:
        fail(f"{ADDENDUM_WAIVER} needs one line 'Differences (from CARRYOVER.txt):'")
    quoted = []
    for line in lines[start[0] + 1:]:
        if not line.startswith("> "):
            break
        quoted.append(line[2:])
    # Repair 2 (review 2, finding 1): exactly the content differences CARRYOVER.txt records,
    # each once, and nothing else.
    required = [line for line in carry["lines"] if re.match(r"difference \(\S+\): ", line)]
    if sorted(quoted) != sorted(required):
        fail(f"{ADDENDUM_WAIVER} must quote, as '> ' lines, exactly the 'difference (...)' lines "
             "of CARRYOVER.txt, each once, and nothing else")


def carried_condition(condition, note, carry, evidence_id=None):
    """Revision 19 (plan revision 5): a carried `PROCEED IF` keeps its condition. The
    disposition addresses it as worded; for every reviewed base, commit or tree that the
    condition names (in full or by a prefix of 7 or more hex digits) it states the mapping on a
    line `Mapping: <reviewed id> -> <outgoing id>`. Nothing marks it satisfied automatically."""
    # Repair 1 (review finding 6): any case, and 64-hex packet identities as well as 40-hex
    # Git identities; the reviewed packet identity maps to this packet's identity.
    tokens = [token.lower() for token in
              re.findall(r"(?<![0-9a-fA-F])[0-9a-fA-F]{7,64}(?![0-9a-fA-F])", condition)]
    pairs = list(carry["mapping"])
    if evidence_id:
        pairs.append((carry["reviewed_id"], evidence_id))
    for before, after in pairs:
        if any(before.startswith(token) for token in tokens):
            if not re.search(rf"^Mapping: {before} -> {after}\s*$", note, re.M | re.I):
                fail(f"the carried condition names the reviewed identity {before[:12]}: the "
                     f"disposition needs the line 'Mapping: {before} -> {after}'")


def carried_reading(body, carry, top, root):
    """Revision 19: a reading of the earlier packet is used. LIMITS needs one line containing
    "reading carried over" with the key words of CARRYOVER_LIMITS (case-insensitive), or,
    when carry-over refused and ADDENDUM_WAIVER.md is used, one line with "addendum waived"
    and "because" (the reason); every upstream file a changed file imports or names
    (plan revision 6) is named in LIMITS with "imported by a changed file" or "named by a
    changed file"; and CODE_IDENTITY.txt has a line `reading_carried_over_from: <reviewed
    commit> [...]` naming the reviewed commits of every repository block."""
    limits = [line.lower() for line in body["LIMITS"]]
    # Repair 2 (review 2, finding 2): no LIMITS line may claim that nothing upstream can
    # interact with or affect the change.
    for line in body["LIMITS"]:
        if interaction_claim(line):
            fail("a LIMITS line claims that upstream changes cannot interact with or affect the "
                 f"change, which the carry-over check does not show: {line[:80]!r}")
    if carry["refused"]:
        if not any("addendum waived" in line and "because" in line for line in limits):
            fail("a reading used with ADDENDUM_WAIVER.md needs a LIMITS line stating that the "
                 "addendum was waived and why ('addendum waived because ...')")
    else:
        if not any(carryover_limits_line(line) for line in body["LIMITS"]):
            fail("a carried-over reading must be stated in LIMITS on one line with the exact "
                 f"sentence, and without claiming that nothing can interact: '{CARRYOVER_LIMITS}'")
    for verb, label in (("imported", "imported by a changed file"), ("named", "named by a changed file")):
        for path in carry["flagged"][verb]:
            if not any(label in line and path.lower() in line for line in limits):
                fail(f"CARRYOVER.txt flags {path} as {label}: LIMITS must state it on a line "
                     f"with '{label}' and the path")
    identity = top.get("CODE_IDENTITY.txt", root / "CODE_IDENTITY.txt")
    try:
        text = identity.read_text(encoding="utf-8") if identity.is_file() else ""
    except UnicodeDecodeError:
        fail("code identity is not UTF-8")
    named = re.findall(r"^reading_carried_over_from: (.+)$", text, re.M)
    if len(named) != 1 or sorted(named[0].split()) != carry["reviewed_commits"]:
        fail("a carried-over reading needs CODE_IDENTITY.txt 'reading_carried_over_from: "
             "<reviewed commit>' naming the reviewed commit of each repository")


# Checkpoint 3 (1b): the top-level files the checks read, which a correction
# companion may supply when a verified copy was frozen in the wrong place.
CORRECTABLE = ("CODE_IDENTITY.txt", "RUN_LOG.txt", "CHANGE.diff", "APPROVAL_REPORT.md",
               "OUTGOING.txt", "SERVER_SUITE.txt", "ROSTER_COMPARISON.txt")
# Repair 1 (review finding 1): never supplied by a correction companion, because the PUBLISH
# precheck recomputes CARRYOVER.txt beside OUTGOING.txt and must not be steered around.
NOT_CORRECTABLE = ("CARRYOVER.txt", "ADDENDUM_WAIVER.md")


def corrections_map(root, corrections):
    """Checkpoint 3 (1b): misplaced verified evidence named by a correction
    companion. `corrections` is a directory outside the frozen packet holding
    ERRATA.md, whose lines `Relocated: NAME from PATH` say that the top-level
    file NAME was frozen at the packet-relative PATH instead, and a copy of
    it as corrections/NAME. The copy is used for NAME only when NAME is
    absent from the packet's top level, PATH is inside the packet and listed
    in its MANIFEST.sha256, PATH's own file name is NAME (checkpoint 5: a
    file frozen under another name, such as notes/draft.txt, cannot stand in
    for SERVER_SUITE.txt), no other line names the same PATH, and the copy's
    SHA-256 equals that listed hash (verify() has already matched the packet
    file to it). Anything else is refused, so no content enters that was not
    frozen and verified under that name. Other ERRATA.md lines are not read.
    Returns {NAME: copy path}."""
    if corrections is None:
        return {}
    if corrections.is_symlink() or not corrections.is_dir():
        fail("corrections must be an ordinary directory outside the packet")
    inside, packet = corrections.resolve(), root.resolve()
    if inside == packet or packet in inside.parents:
        fail("corrections directory must be outside the frozen packet")
    errata = corrections / "ERRATA.md"
    if errata.is_symlink() or not errata.is_file():
        fail("corrections directory needs an ordinary ERRATA.md")
    try:
        text = errata.read_bytes().decode("utf-8")
    except UnicodeDecodeError:
        fail("ERRATA.md is not UTF-8")
    listed = dict(reversed(line.split("  ", 1)) for line in
                  (root / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines())
    found, sources = {}, set()
    for line in text.splitlines():
        match = re.fullmatch(r"Relocated: (\S+) from (\S+)\s*", line)
        if not match:
            continue
        name, source = match[1], match[2]
        if name in NOT_CORRECTABLE:
            fail(f"ERRATA.md relocates {name}, which a correction companion may never supply "
                 "(it must be frozen at the packet top level)")
        if name not in CORRECTABLE:
            fail(f"ERRATA.md relocates a file the checker does not read: {name}")
        if name in found:
            fail(f"ERRATA.md relocates {name} twice")
        if (source.startswith("/") or "\\" in source
                or any(part in ("", ".", "..") for part in source.split("/"))):
            fail(f"ERRATA.md names a path outside the packet: {source}")
        if source not in listed:
            fail(f"ERRATA.md names a file not listed in the frozen manifest: {source}")
        if (root / name).exists() or (root / name).is_symlink():
            fail(f"ERRATA.md relocates {name}, but the packet has it at the top level")
        if source.split("/")[-1] != name:
            fail(f"ERRATA.md relocates {name} from a file of another name: {source}")
        if source in sources:
            fail(f"ERRATA.md names {source} for more than one file")
        sources.add(source)
        copy = corrections / name
        if copy.is_symlink() or not copy.is_file():
            fail(f"corrections directory has no ordinary copy of {name}")
        if digest(copy.read_bytes()) != listed[source]:
            fail(f"corrected copy of {name} differs from the frozen {source}")
        found[name] = copy
    return found


REVIEW_HEADING = "## Plain-English review"
REVIEW_LABELS = ("Purpose", "Code", "Tests", "Assessment")


def plain_english_review(body):
    """Revision 19 (item 0.2): the approval report opens with `## Plain-English review` (it
    replaces `## Summary`), one block per proposed commit under the labels Purpose, Code
    (short before and after snippets), Tests (key assertions, whether they cover the fix,
    the gaps, or why none were added) and Assessment. Form only: the heading once and each
    label at the start of a line in that section (optionally bulleted or bold); what they say
    is the Developer's and the reviewer's to judge. A report frozen before revision 19 has
    `## Summary` instead and is refused here; frozen records are not rewritten."""
    lines = body.splitlines()
    starts = [index for index, line in enumerate(lines) if line.strip() == REVIEW_HEADING]
    if len(starts) != 1:
        fail("approval report needs one '## Plain-English review' section (it replaces "
             "'## Summary'; templates/APPROVAL_REPORT.md)")
    section = []
    for line in lines[starts[0] + 1:]:
        if line.startswith("## "):
            break
        section.append(line)
    for label in REVIEW_LABELS:
        if not any(re.match(rf"\s*(?:[-*+]\s+)?(?:\*\*)?{label}:(?:\*\*)?(?:\s|$)", line)
                   for line in section):
            fail(f"the Plain-English review needs a '{label}:' part (Purpose, Code, Tests, "
                 "Assessment)")


def code_identity(root, checked, changed, top=None):
    top = top or {}
    path = top.get("CODE_IDENTITY.txt", root / "CODE_IDENTITY.txt")
    log = top.get("RUN_LOG.txt", root / "RUN_LOG.txt")
    diff = top.get("CHANGE.diff", root / "CHANGE.diff")
    for required in (path, log, diff):
        if not required.is_file() or not required.stat().st_size:
            fail(f"missing or empty change evidence: {required.name}")
    report = top.get("APPROVAL_REPORT.md", root / "APPROVAL_REPORT.md")
    if not report.is_file() or report.is_symlink():
        fail("missing approval report")
    try:
        body = report.read_text(encoding="utf-8")
        patch = diff.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        fail("approval report and change diff must be UTF-8")
    if ("## Exact tested changes" not in body or
            "## New tests" not in body or patch not in body):
        fail("approval report must contain the complete exact CHANGE.diff and new-test section")
    plain_english_review(body)
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


# Revision 17: a shell word for an absolute path, bare or quoted. Unquoted, only characters
# no supported shell treats specially (letters, digits and _ . / + @ % : , = -), which is
# tighter than revision 16's any non-space; anything else must be quoted. Inside single
# quotes nothing is special in sh, bash, zsh, csh or tcsh except the closing quote and csh's
# history character `!` (refused in every form); inside double quotes `$`, a backquote and a
# backslash are also refused. Revision 18 (repair 1, reading 1, finding 2; repair 1b): inside
# either kind of quotes every character must also be printable (str.isprintable: letters with
# accents and other scripts pass; tab, DEL and other control characters such as Ctrl-U, CR and ESC,
# format characters such as the bidi overrides U+202E and U+2066 and U+200B, line and paragraph
# separators, other non-space separators and surrogates do not), see printable_path(), so that a
# terminal or line editor has nothing to act on when the line is pasted; before, anything but the
# excluded characters and a newline.
SAFE_BARE = r"[A-Za-z0-9_./+@%:,=-]+"
QUOTED_SINGLE = r"[^'!]"  # with printable_path(): printable, except ' and !
QUOTED_DOUBLE = r"[^\"$`\\!]"  # also except " $ \ and a backquote
PATH_WORD = (r"(?:(?P<{0}_bare>/" + SAFE_BARE + r")|'(?P<{0}_single>/" + QUOTED_SINGLE + r"+)'"
             r"|\"(?P<{0}_double>/" + QUOTED_DOUBLE + r"+)\")")
# Revision 18: any absolute path to a Python 3 interpreter, by its file name: python, python3,
# python3.N or libtbx.python (a venv's bin/python, a system or Homebrew python3.N, the cctbx
# conda Python and the python that libtbx.python starts on macOS). Before: only
# /.../libtbx.python or /.../conda_base/bin/python3[.N].
INTERPRETER = re.compile(r"/.+/(?:python(?:3(?:\.\d+)?)?|libtbx\.python)")


def path_word(match, name):
    """The path a PATH_WORD group matched, without its quotes."""
    return next(match.group(f"{name}_{form}") for form in ("bare", "single", "double")
                if match.group(f"{name}_{form}") is not None)


def printable_path(path):
    """Revision 18 (repair 1b): every character printable (str.isprintable), as a quoted stop path
    needs; refuses control, format and separator characters other than the space, and surrogates."""
    return path.isprintable()


def stop_line(line):
    """The auto_start stop command: `INTERPRETER -I -B AUTO_QUEUE stop`, both absolute,
    each bare or quoted (a path with a space must be quoted); the interpreter is an
    existing executable regular file named python, python3, python3.N or libtbx.python
    (revision 18; before, only libtbx.python or conda_base/bin/python3[.N])."""
    match = re.fullmatch(PATH_WORD.format("python") + " -I -B " + PATH_WORD.format("queue")
                         + " stop", line)
    if not match or path_word(match, "queue") != str(AUTO_QUEUE):
        return False
    interpreter = path_word(match, "python")
    # Repair 1 (reading 1, finding 5): the checker accepts only what stop_command can write,
    # so a single quote (which stop_command cannot put inside single quotes) is refused in
    # the double-quoted form too.
    if "'" in interpreter or "'" in path_word(match, "queue"):
        return False
    if not (printable_path(interpreter) and printable_path(path_word(match, "queue"))):
        return False
    if not INTERPRETER.fullmatch(interpreter):
        return False
    resolved = os.path.realpath(interpreter)
    return os.path.isfile(resolved) and os.access(resolved, os.X_OK)


def stop_command(interpreter=None, queue=AUTO_QUEUE):
    """Revision 17: the stop line for the auto_start screen, each path bare when it has no
    space or quote, otherwise in single quotes; refused unless stop_line() accepts it (an
    absolute, existing, executable python, python3, python3.N or libtbx.python, revision 18).
    Revision 18: without an interpreter, the absolute sys.executable of the Python running
    this checker, which must then have psutil (the helper the stop line runs needs it)."""
    def word(path):
        path = str(path)
        return path if re.fullmatch(SAFE_BARE, path) else "'%s'" % path
    if interpreter is None:
        interpreter = sys.executable
        if not interpreter or not os.path.isabs(interpreter):
            fail("this Python reports no absolute interpreter path (sys.executable); "
                 "name the interpreter: stop-command INTERPRETER")
        if importlib.util.find_spec("psutil") is None:
            fail(f"this Python ({word(interpreter)}) has no psutil, which auto_queue.py needs: "
                 "make a venv and install it there (python3 -m venv DIR, then DIR/bin/python -m "
                 "pip install psutil; not --user, which the kit's -I ignores), then run this "
                 "command with DIR/bin/python")
    line = f"{word(interpreter)} -I -B {word(queue)} stop"
    if not stop_line(line):
        fail("the interpreter must be an absolute path to an existing executable Python 3 "
             "interpreter named python, python3, python3.N or libtbx.python, without a ' or ! "
             f"character or a character that is not printable: {interpreter!r}")
    return line


def auto_screen(kind, body, queue):
    """Auto screens: the start notice states its limits; the report repeats the helper's rows."""
    if kind == "auto_start":
        if queue:
            fail("auto_start has no queue to compare")
        # Revision 16: the stop command runs the helper with a Python that has psutil, by its
        # absolute path so that it works from any Terminal whatever PATH is (repair 1, finding
        # 3). Revision 18: any absolute python, python3, python3.N or libtbx.python (before,
        # only libtbx.python or conda_base/bin/python3[.N]); stop-command without an argument
        # prints the absolute sys.executable of the Python that runs it. A bare or relative
        # name is refused. Repair 2 (reading 2, finding 2): the interpreter must be an
        # existing, executable regular file (symbolic links resolved). Revision 17 (outside
        # reading of revision 16): either path may be written in single or double quotes,
        # which a path containing a space needs ('/Users/a b/.../libtbx.python' -I -B ...);
        # an unquoted path with a space, and a quoted relative or bare name, are refused.
        if not any(stop_line(line) for line in body["STOP"]):
            fail("auto_start needs the stop command of this procedure's auto_queue.py, run by an "
                 "existing executable interpreter")
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
    if helper.returncode == 18:
        # Revision 16: the helper needs psutil.
        fail("auto_queue.py report failed: " + (helper.stderr.strip().splitlines() or [""])[-1])
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


# Checkpoint 3: every verb form of integrate, merge, push and publish (the nouns
# "integration" and "publication" stay allowed, for PENDING and the follow-up note).
ACCEPT_CLAIMS = re.compile(r"\b(?:integrat(?:e|ed|es|ing)|merg(?:e|ed|es|ing)|"
                           r"push(?:ed|es|ing)?|publish(?:ed|es|ing)?)\b", re.IGNORECASE)


def accept_screen(body):
    """ACCEPT accepts one exact frozen candidate while the outside reading is
    pending. It names what is still pending, says that master and origin
    stay unchanged, and neither offers nor claims integration or
    publication (those are the RESULT and PUBLICATION decisions): any verb
    form of integrate, merge, push or publish is refused. Checkpoint 3: the
    DECISION's required words are now "master and origin stay unchanged"
    (replacing "merges nothing and publishes nothing", which used the
    refused verbs)."""
    text = "\n".join(line for lines in body.values() for line in lines)
    claim = ACCEPT_CLAIMS.search(text)
    if claim:
        fail(f"accept screen must not offer or claim integration or publication: {claim.group()!r}")
    if not any(re.fullmatch(r"Pending: \S.*", line) for line in body["PENDING"]):
        fail("accept screen needs a 'Pending:' line naming what is still pending")
    decision = " ".join(body["DECISION"]).lower()
    if "master and origin stay unchanged" not in decision:
        fail("accept DECISION must say master and origin stay unchanged")


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


def check(kind, screen, root=None, reading=None, disposition=None, report=True, queue=None,
          corrections=None, waiver=None):
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
    elif kind in ("publication", "accept"):
        if len(parts) != 3 or parts[0] != kind.upper() or not re.fullmatch(NAME, parts[1]):
            fail(f"invalid {kind} header")
    elif len(parts) != 2 or parts[0] != "STOP" or not re.fullmatch(NAME, parts[1]):
        fail("invalid stop header")
    body = blocks(kind, lines)
    developer_view(kind, body)
    # Revision 17: the outside-reading waiver serves only a full-path RESULT, in place of a
    # reading; never with one, never with a disposition and never for publication.
    if waiver is not None:
        if reading is not None:
            fail("give either --reading or --waiver, not both")
        if disposition is not None:
            fail("a waiver replaces the reading; it takes no disposition")
        if kind == "result" and parts[2:3] == ["light"]:
            fail("--waiver serves only a full-path result; a light result defers its reading "
                 "to publication")
        if kind != "result" or parts[2:3] != ["full"]:
            fail("--waiver serves only a full-path result (integration); it never covers "
                 "publication")
    if kind in ATTENTION:
        if root or reading or disposition:
            fail("auto screens have no frozen evidence or outside reading")
        auto_screen(kind, body, queue)
    elif queue:
        fail("only auto_report compares with a queue")
    if kind in ("result", "publication", "accept"):
        if root is None or not re.fullmatch(SHA, parts[-1]):
            fail("evidence directory and manifest identity required")
        identity = verify(root)
        if parts[-1] != identity:
            fail("screen does not identify the current evidence")
        top = corrections_map(root, corrections)
        # Checkpoint 5 (review CP3 finding 5): post-reading errata are named in LIMITS.
        if top and not any("ERRATA.md" in line for line in body["LIMITS"]):
            fail("a correction companion was used: LIMITS must name ERRATA.md")
        if report:
            for name, copy in sorted(top.items()):
                print(f"NOTE {name} read from the correction companion {copy} "
                      "(byte-identical to the frozen file ERRATA.md names)")
        if kind == "accept":
            if reading or disposition:
                fail("accept is used while the outside reading is pending; it takes no reading")
            code_identity(root, body["CHECKED"], body["CANDIDATE"], top)
            accept_screen(body)
        elif kind == "result":
            code_identity(root, body["CHECKED"], body["CHANGED"], top)
            if parts[2] == "light":
                if reading or disposition or body["GATE"] != ["Reading: deferred to publication", "Verdict: deferred to publication"]:
                    fail("light-path outside reading belongs to publication batch")
            elif waiver is not None:
                waiver_matches(waiver, identity, body["GATE"], root)
                # The word WAIVED on a line naming the outside reading or review; "not WAIVED",
                # in any case, is refused (repair 1, reading 1 finding 10).
                if not any(re.search(r"\bWAIVED\b", line)
                           and not re.search(r"\bnot\s+WAIVED\b", line, re.I)
                           and re.search(r"\boutside (?:reading|review)\b", line, re.I)
                           for line in body["LIMITS"]):
                    fail("a waived outside reading must be named in LIMITS (WAIVED, on a line "
                         "naming the outside reading)")
            else:
                status = reading_matches("result", reading, identity, body["GATE"], disposition)
                if status in ("WAIVED", "PENDING") and not any(
                        status in line for heading in ("LIMITS", "DECISION")
                        for line in body[heading]):
                    fail("waived or pending condition must be visible in Developer View")
        else:
            # Checkpoint 3: companions inside the packet are refused here as in
            # check-draft (before, only the draft refused them).
            for name, _ in regular_files(root):
                if Path(name).name in COMPANION_NAMES:
                    fail(f"companion file inside the packet (keep it outside): {name}")
            for filename in ("SERVER_SUITE.txt", "ROSTER_COMPARISON.txt"):
                item = top.get(filename, root / filename)
                if not item.is_file() or not item.stat().st_size:
                    fail(f"missing or empty publication evidence: {filename}")
            batch = "\n".join(body["BATCH"])
            for block in parse_outgoing(top.get("OUTGOING.txt", root / "OUTGOING.txt")):
                if not re.search(rf"(?<![0-9A-Za-z]){block['commit']}(?![0-9A-Za-z])", batch):
                    fail("outgoing commit is not named in BATCH")
            carry = carryover_evidence(root, top)
            suite_waiver(top.get("SERVER_SUITE.txt", root / "SERVER_SUITE.txt"))
            waived = (root / ADDENDUM_WAIVER).exists() or (root / ADDENDUM_WAIVER).is_symlink()
            names_reviewed = bool(carry and reading and reading.is_file()
                                  and reading_identities(reading) == {carry["reviewed_id"]} != {identity})
            if waived and not (carry and carry["refused"]):
                fail(f"{ADDENDUM_WAIVER} is used only when CARRYOVER.txt records a refusal")
            if carry and carry["refused"] and names_reviewed:
                if not waived:
                    fail("CARRYOVER.txt records a refusal: the reading does not carry over; ask "
                         f"for an addendum, or record the Developer's waiver in {ADDENDUM_WAIVER}")
                addendum_waiver(root, carry)
            elif carry and carry["refused"]:
                carry = None  # an addendum reading of this packet; the refusal is only evidence
            status = reading_matches("publication", reading, identity, body["GATE"], disposition,
                                     carried=carry)
            if carry and reading_identities(reading) == {carry["reviewed_id"]} != {identity}:
                carried_reading(body, carry, top, root)
            if status in ("WAIVED", "PENDING") and not any(
                    status in line for heading in ("LIMITS", "DECISION")
                    for line in body[heading]):
                fail("waived or pending condition must be visible in Developer View")
    elif kind not in ATTENTION and (root or reading or disposition):
        fail("plan and stop have no frozen evidence or outside reading")
    if corrections is not None and kind not in ("result", "publication", "accept"):
        fail("only result, publication and accept read a correction companion")
    if report:
        print(f"CHECK OK {kind}")
    return top if kind in ("result", "publication", "accept") else {}


PUBLICATION_FILES = ("OUTGOING.txt", "SERVER_SUITE.txt", "ROSTER_COMPARISON.txt")
PUBLICATION_LISTED = ("README_FIRST.md", "PROOF_SUMMARY.md", "APPROVAL_REPORT.md",
                      "PUBLICATION_NOTES.md")
COMPANION_NAMES = ("APPROVALS.md", "CHECKER.md", "ERRATA.md")


def check_draft(kind, root):
    """Run, on an unfrozen publication packet, every publication check that
    needs neither the freeze, the screen nor the outside reading (see
    payload/PUBLICATION_PACKET.md), and say which checks were skipped.

    Checkpoint 3: the OUTGOING.txt line rule (column-0 key lines of the
    seven keys, prose refused), the SERVER_SUITE.txt first-line form and
    the companion refusal are the final check's own rules (before, the draft
    alone applied the first and the last). A passing draft is not a frozen
    packet, a reading or a PUBLISH decision.
    """
    if kind != "publication":
        fail("check-draft supports only publication")
    files = dict(regular_files(root))
    for name in PUBLICATION_FILES:
        item = root / name
        if name not in files or not item.stat().st_size:
            fail(f"missing or empty publication evidence at the packet top level: {name}")
    for name in files:
        if Path(name).name in COMPANION_NAMES:
            fail(f"companion file inside the packet (keep it outside): {name}")
    blocks = parse_outgoing(root / "OUTGOING.txt")
    carry = carryover_evidence(root)
    suite_waiver(root / "SERVER_SUITE.txt")
    if (root / ADDENDUM_WAIVER).exists():
        if not (carry and carry["refused"]):
            fail(f"{ADDENDUM_WAIVER} is used only when CARRYOVER.txt records a refusal")
        addendum_waiver(root, carry)
        print(f"checked: {ADDENDUM_WAIVER} form against the refused CARRYOVER.txt")
    if carry:
        print("checked: CARRYOVER.txt is consistent and agrees with evidence/reviewed/ and "
              "OUTGOING.txt (the PUBLISH precheck recomputes it)")
    print(f"DRAFT OK publication ({len(blocks)} outgoing repository block(s))")
    print("checked: OUTGOING.txt key lines and blocks, SERVER_SUITE.txt and "
          "ROSTER_COMPARISON.txt at top level, SERVER_SUITE first line and NOT RUN waiver "
          "form, no companions inside")
    if "MANIFEST.sha256" in files:
        verify(root)
        print("checked: existing MANIFEST.sha256 matches the files")
    else:
        print("skipped: freeze and packet identity (not frozen)")
    print("skipped: outgoing commits named in BATCH (needs the PUBLICATION screen)")
    print("skipped: outside reading, its scope and any disposition (needs the reading)")
    for name in PUBLICATION_LISTED:
        if name not in files:
            print(f"NOTE not present (listed in PUBLICATION_PACKET.md, not checked): {name}")


HIDDEN_PREFIXES = ("Tested tree: ", "Reading SHA256: ", "Disposition SHA256: ", "Waiver SHA256: ")
# Revision 19 (item 0.3): characters Markdown would interpret, escaped with a backslash.
MARKDOWN_SPECIAL = re.compile(r"([\\`*_\[\]<>#|~&!{}$])")
COMMAND_LINE = re.compile(r"""['"]?/\S""")


def escape(text, table=False):
    """Text as Markdown shows it literally: backslash before each special character (not
    `|` inside a table row), before each character of a leading run of -, + or =, and after a
    leading number followed by . or )."""
    escaped = MARKDOWN_SPECIAL.sub(lambda m: m.group() if table and m.group() == "|" else "\\" + m.group(),
                                   text)
    return re.sub(r"^(\s*)([-+=]+|\d+(?=[.)]))", lambda m: m[1] + (
        "".join("\\" + char for char in m[2]) if not m[2][0].isdigit() else m[2] + "\\"), escaped)


def code_span(text):
    """An inline code span that shows `text` exactly (a command to copy)."""
    longest = max((len(m.group()) for m in re.finditer(r"`+", text)), default=0)
    fence = "`" * (longest + 1)
    pad = " " if text.startswith("`") or text.endswith("`") else ""
    return f"{fence}{pad}{text}{pad}{fence}"


def title(words):
    return words.capitalize()


def render(kind, lines):
    """Revision 19 (item 0.3): the checked screen as a Markdown block quote, in the normal
    typeface. The header becomes **Kind** · id · level; each section name is bold title case
    (the saved file keeps CAPITALS, and every check reads the file); every source line is
    its own paragraph, which the window wraps (repair 2; before, lines were joined), table rows
    stay rows, and a command line (an absolute path) is shown as code;
    a blank quoted line separates paragraphs and sections; the action line becomes **Your
    choice:** with the options in bold. Everything else is escaped, so it shows literally."""
    headings = HEADINGS[kind]
    parts = lines[0].split(" | ")
    out = [" · ".join([f"**{title(parts[0].replace('_', ' '))}**"] + [escape(p) for p in parts[1:]])]
    paragraphs = []  # (kind, text) with kind "text", "table", "code", "heading"
    for line in lines[1:-1]:
        if line in headings:
            paragraphs.append(("heading", f"**{title(line)}**"))
        elif line.lstrip().startswith("|"):
            if paragraphs and paragraphs[-1][0] == "table":
                paragraphs[-1] = ("table", paragraphs[-1][1] + "\n" + escape(line, table=True))
            else:
                paragraphs.append(("table", escape(line, table=True)))
        elif COMMAND_LINE.match(line):
            paragraphs.append(("code", code_span(line)))
        else:
            # Repair 2 (review 2, note 7): every source line is its own paragraph, so nothing
            # merges into the line before it; the reader's window still wraps each one.
            paragraphs.append(("text", escape(line.strip())))
    for _, text in paragraphs:
        out.append("")
        out += text.split("\n")
    action = lines[-1].removeprefix("ACTION: ")
    options = []
    for option in action.split(" / "):
        match = re.fullmatch(r"([A-Z][A-Z ]*[A-Z]|[A-Z])( .*)?", option)
        options.append(f"**{title(match[1])}**" + escape(match[2] or "") if match else escape(option))
    out += ["", "**Your choice:** " + " / ".join(options)]
    return [("> " + line) if line else ">" for line in out]


def present(kind, screen, root=None, reading=None, disposition=None, queue=None,
            corrections=None, waiver=None):
    """Check the saved decision, then render a conspicuous Developer view."""
    top = check(kind, screen, root, reading, disposition, report=False, queue=queue,
                corrections=corrections, waiver=waiver)
    lines = screen.read_text(encoding="utf-8").splitlines()
    if kind in ("result", "publication", "accept"):
        lines[0] = " | ".join(lines[0].split(" | ")[:-1])
    lines = [line for line in lines if not line.startswith(HIDDEN_PREFIXES)]
    lines = [DISPLAY_ID.sub(
        lambda match: "[recorded]" if any(c in "abcdef" for c in match.group())
        else match.group(), line) for line in lines]
    print(ATTENTION.get(kind, "**PLEASE READ — YOUR DECISION IS NEEDED**") + "\n\n" +
          "\n".join(render(kind, lines)))
    # Revision 19: the rendered screen says that the reading was carried over.
    if kind == "publication" and root is not None:
        carry = carryover_evidence(root, top)
        if carry and reading_identities(reading) == {carry["reviewed_id"]}:
            print("NOTE: the outside reading names the earlier packet it read; " +
                  ("the carry-over check refused and the Developer waived the addendum "
                   f"({ADDENDUM_WAIVER})" if carry["refused"] else
                   "it is carried over by CARRYOVER.txt, which the PUBLISH precheck recomputes"))
    # Checkpoint 5 (review CP3 finding 5): the rendered screen says that a
    # correction companion supplied evidence, and names it.
    for name, copy in sorted(top.items()):
        print(f"NOTE: correction companion used: {name} was read from {copy} "
              f"(ERRATA.md in {copy.parent}; byte-identical to the frozen file it names)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("freeze", "verify"):
        commands.add_parser(name).add_argument("evidence", type=Path)
    commands.add_parser("verify-source").add_argument("source", type=Path)
    commands.add_parser("check-claude-version")
    commands.add_parser("register-skill").add_argument("source", type=Path)
    commands.add_parser("stop-command", help="revision 17: print the auto_start stop line, "
                        "quoting a path that needs it; revision 18: without INTERPRETER, the "
                        "absolute path of the Python running this command").add_argument(
                            "interpreter", nargs="?")
    draft = commands.add_parser("check-draft")
    draft.add_argument("kind", choices=("publication",))
    draft.add_argument("evidence", type=Path)
    command = commands.add_parser("check")
    presenter = commands.add_parser("present")
    for command in (command, presenter):
        command.add_argument("kind", choices=tuple(HEADINGS))
        command.add_argument("screen", type=Path)
        command.add_argument("--evidence", type=Path)
        command.add_argument("--reading", type=Path)
        command.add_argument("--disposition", type=Path)
        command.add_argument("--queue", type=Path)
        command.add_argument("--corrections", type=Path,
                             help="correction companion directory (ERRATA.md and byte-identical "
                                  "copies of misplaced frozen files), outside the packet")
        command.add_argument("--waiver", type=Path,
                             help="revision 17: OUTSIDE_WAIVER.md, the Developer's waiver of the "
                                  "outside reading for a full-path result; not with --reading")
    args = parser.parse_args()
    try:
        if args.command == "freeze":
            freeze(args.evidence)
        elif args.command == "verify":
            identity = verify(args.evidence)
            print("VERIFIED evidence")
            # Revision 20: the packet identity, the SHA-256 of MANIFEST.sha256.
            print(f"packet identity: {identity}")
        elif args.command == "verify-source":
            verify_source(args.source)
            print("VERIFIED complete source")
        elif args.command == "check-claude-version":
            check_claude_version()
        elif args.command == "register-skill":
            register_skill(args.source)
        elif args.command == "stop-command":
            print(stop_command(args.interpreter))
        elif args.command == "check-draft":
            check_draft(args.kind, args.evidence)
        elif args.command == "check":
            check(args.kind, args.screen, args.evidence, args.reading, args.disposition,
                  queue=args.queue, corrections=args.corrections, waiver=args.waiver)
        else:
            present(args.kind, args.screen, args.evidence, args.reading, args.disposition,
                    queue=args.queue, corrections=args.corrections, waiver=args.waiver)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
