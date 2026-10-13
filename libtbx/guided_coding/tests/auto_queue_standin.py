"""Test-only launcher for auto_queue.py (Revision 12). Not part of the helper.

The kit tests use it, instead of running auto_queue.py itself, only where the
platform is not macOS (real, or simulated inside the test processes with
GC_TEST_SIMULATE_NON_MACOS=1, a variable read only by the tests). It imports
auto_queue and, inside this one process only, (a) accepts the platform and
(b) replaces /usr/bin/sandbox-exec with a stand-in that checks the
`-p PROFILE CMD...` form, ignores the profile and executes CMD unrestricted
(in the same process, as sandbox-exec does). The installed auto_queue.py has
no hook, variable or option for this; run directly, auto_queue.py refuses
`run` off macOS.

This launcher refuses real queues (Revision 12, P8). Before it imports the
helper, starts anything or writes anything, it exits 16 unless every queue
path it checks (the helper home, every --queue, the queue named by the home's
ACTIVE file, and init's records root and lock directory) lies inside one
directory directly under the temporary
directory (tempfile.gettempdir()) that holds the test-only marker file
GC_KIT_TEST_QUEUE_ROOT with exactly MARKER_TEXT. It resolves those paths as
auto_queue.py does: the helper home is --home, else GC_AUTO_HOME, else the
default home, which is refused. An empty or blank --home, --queue, --records
or --lock-dir, or an empty or blank GC_AUTO_HOME used as the home, is refused
too, so that the paths checked are the ones the helper opens: the helper
treats a blank value (for example --home " ") as a relative path, falls back
to another home for an empty --home or GC_AUTO_HOME (and to the ACTIVE queue
for an empty --queue), resolves an empty --records to the current directory,
and refuses an empty --lock-dir itself. Checked: every --home or GC_AUTO_HOME, every
--queue, the queue named by each home's ACTIVE file, and init's --records and
--lock-dir. The kit tests write the marker into the temporary directory they
create for each test; real queues (under .claude/records, ~/GuidedCoding/auto
or anywhere else without it) are refused. Limits: a marker written by hand
inside the temporary directory is accepted, and so is a symlink planted by
hand inside a marked test queue (for example its runs/ pointing to a real
queue's runs/, which then receives the run record), because the check
resolves the queue directory, not what is inside it. This guards against
mistaken use, not deliberate misuse. Only those queue paths are checked:
other paths the helper opens or writes (--log, --cwd, --file, --packet,
--script, --grant, --jobs) are not, so a run started through the launcher
can write its log anywhere.

Because the stand-in does not block the network, the tests run through it
only controlled local commands; the tests of the real network block and of
the `ps -E` run-id marker are skipped off macOS.

    python3 -I -B auto_queue_standin.py [auto_queue.py arguments]
"""

import os
import subprocess
import sys
import tempfile
import types

STANDIN = os.path.realpath(__file__)
FLAG = "--sandbox-standin"
MARKER = "GC_KIT_TEST_QUEUE_ROOT"
MARKER_TEXT = "guided_coding kit tests: test-only queue root for this test run; never a real queue\n"
REFUSED = 16


def sandbox_standin(argv):
    """Stand-in for `sandbox-exec -p PROFILE CMD [ARG...]`."""
    if len(argv) < 3 or argv[0] != "-p":
        print("sandbox stand-in: expected -p PROFILE CMD [ARG...]", file=sys.stderr)
        return 2
    try:
        os.execvp(argv[2], argv[2:])
    except OSError as error:
        print("sandbox stand-in: cannot run %s: %s" % (argv[2], error), file=sys.stderr)
        return 127


def option_values(argv, name):
    """Every value given as `NAME VALUE` or `NAME=VALUE` before `--` (the helper's
    parsers use allow_abbrev=False, so no other spelling reaches it)."""
    tokens = argv[:argv.index("--")] if "--" in argv else argv
    values = []
    for index, token in enumerate(tokens):
        if token == name and index + 1 < len(tokens):
            values.append(tokens[index + 1])
        elif token.startswith(name + "="):
            values.append(token[len(name) + 1:])
    return values


def test_root(path, temp):
    """The directory directly under `temp` that contains `path` (resolved), or None."""
    real = os.path.realpath(path)
    if os.path.commonpath([real, temp]) != temp or real == temp:
        return None
    return os.path.join(temp, os.path.relpath(real, temp).split(os.sep)[0])


def marked(root):
    marker = os.path.join(root, MARKER)
    if os.path.islink(marker) or not os.path.isfile(marker):
        return False
    try:
        with open(marker, encoding="utf-8") as handle:
            return handle.read() == MARKER_TEXT
    except (OSError, UnicodeDecodeError):
        return False


def not_a_test_queue(argv, environ):
    """None when every checked queue path (helper home, --queue, the ACTIVE queue,
    init's --records and --lock-dir) is inside one marked test root; else the
    reason. Other paths (--log, --cwd, --file, ...) are not checked. Mirrors auto_queue.home_dir/open_queue. An empty or blank
    value is refused here: the helper would read it as absent, as the current
    directory or as a relative path, not as the value checked here."""
    temp = os.path.realpath(tempfile.gettempdir())
    for name in ("--home", "--queue", "--records", "--lock-dir"):
        if any(not value.strip() for value in option_values(argv, name)):
            return ("empty or blank %s (refused so that the paths checked are the ones "
                    "the helper opens)" % name)
    homes = option_values(argv, "--home")
    if not homes:
        if not environ.get("GC_AUTO_HOME", "").strip():
            return ("no helper home given (a non-empty --home or GC_AUTO_HOME); "
                    "the default home is refused")
        homes = [environ["GC_AUTO_HOME"]]
    paths = [("helper home", home) for home in homes]
    paths += [("queue", queue) for queue in option_values(argv, "--queue")]
    for home in homes:
        active = os.path.join(home, "ACTIVE")
        if os.path.isfile(active):
            with open(active, encoding="utf-8", errors="replace") as handle:
                lines = handle.read().splitlines()
            if lines and lines[0].strip():
                paths.append(("queue named by ACTIVE", lines[0].strip()))
    paths += [("records root", value) for value in option_values(argv, "--records")]
    paths += [("lock directory", value) for value in option_values(argv, "--lock-dir")]
    roots = set()
    for what, path in paths:
        root = test_root(path, temp)
        if root is None:
            return "%s %s is not inside the temporary directory %s" % (what, path, temp)
        if not marked(root):
            return "%s %s: %s has no test-only marker %s" % (what, path, root, MARKER)
        roots.add(root)
    if len(roots) != 1:
        return "the queue paths are in different test directories: %s" % ", ".join(sorted(roots))
    return None


def launch(argv):
    reason = not_a_test_queue(argv, os.environ)
    if reason:
        print("REFUSED: tests/auto_queue_standin.py is a test-only launcher and runs only "
              "queues the kit's own tests created; nothing was started or recorded. Use "
              "auto_queue.py for real queues. Not a test queue: %s" % reason, file=sys.stderr)
        return REFUSED
    sys.path.insert(0, os.path.dirname(os.path.dirname(STANDIN)))
    import auto_queue
    auto_queue.platform_supported = lambda: True
    auto_queue.SANDBOX = STANDIN

    def popen(args, *more, **options):
        if isinstance(args, (list, tuple)) and args and args[0] == STANDIN:
            args = [sys.executable, "-I", "-B", STANDIN, FLAG] + list(args[1:])
        return subprocess.Popen(args, *more, **options)

    names = {name: getattr(subprocess, name) for name in dir(subprocess)
             if not name.startswith("__")}
    names["Popen"] = popen
    auto_queue.subprocess = types.SimpleNamespace(**names)
    return auto_queue.main(argv)


if __name__ == "__main__":
    if sys.argv[1:2] == [FLAG]:
        sys.exit(sandbox_standin(sys.argv[2:]))
    sys.exit(launch(sys.argv[1:]))
