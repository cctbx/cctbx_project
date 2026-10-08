#!/usr/bin/env python3
"""Read-only Git checks before an authorized push: effective destinations,
forbidden push settings, outgoing bindings, an optional dry run, and the
form of a proposed push command. It never pushes; it fetches only with
--fetch. Passing establishes that the recorded proposal matches the live
repository, not that the Developer has decided to publish."""

import sys
# When run as a script, keep this directory out of module resolution before
# importing the standard library (the documented -I invocation also does).
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
from pathlib import Path
import re
import subprocess
import types

# Load the reviewed sibling explicitly from its source text, as
# review_bundle.py does, so cached bytecode is never executed.
_checker_path = Path(__file__).with_name("screen_check.py")
_checker = types.ModuleType("screen_check")
_checker.__file__ = str(_checker_path)
exec(compile(_checker_path.read_bytes(), str(_checker_path), "exec"), _checker.__dict__)
parse_outgoing, fail = _checker.parse_outgoing, _checker.fail
COMMIT, BRANCH = _checker.COMMIT, _checker.BRANCH

TRUE_VALUES = ("true", "yes", "on", "1")


def git(repo, *args):
    """Run a read-only Git query and return its stdout only when Git exited 0.
    A nonzero exit refuses, whatever was printed: output that looks right
    after a failed command is not evidence."""
    result = subprocess.run(["git", "-C", str(repo), *args], text=True, capture_output=True)
    if result.returncode != 0:
        fail(f"git {' '.join(args)} failed (exit {result.returncode}): "
             f"{result.stderr.strip() or result.stdout.strip() or 'no output'}")
    return result.stdout


def git_query(repo, *args):
    """Run a Git query whose nonzero exit may be a legitimate answer; the
    caller decides from (exit, stdout, stderr)."""
    result = subprocess.run(["git", "-C", str(repo), *args], text=True, capture_output=True)
    return result.returncode, result.stdout, result.stderr


def select_block(blocks, name):
    if name is None:
        if len(blocks) != 1:
            fail("several repositories; name one with --repository")
        return blocks[0]
    for block in blocks:
        if block["repository"] == name:
            return block
    fail(f"repository not in outgoing record: {name}")


def check_destinations(repo, block):
    """Exactly one effective push URL, equal to the record, and no setting
    that rewrites, mirrors, adds refspecs or follows tags."""
    remote = block["remote"]
    urls = git(repo, "remote", "get-url", "--all", "--push", remote).splitlines()
    if urls != [block["remote_url"]]:
        fail(f"push destination differs from outgoing record: {' '.join(urls) or 'none'}")
    pattern = (r"^(url\..*\.(insteadof|pushinsteadof)|remote\." + re.escape(remote)
               + r"\.(pushurl|mirror|push|tagopt)|push\..*)$")
    code, listing, errors = git_query(repo, "config", "--show-origin", "--get-regexp", pattern)
    if code == 1 and not listing.strip() and not errors.strip():
        listing = ""  # Git's documented "no match" result
    elif code != 0:
        fail(f"git configuration query failed (exit {code}): {errors.strip() or listing.strip() or 'no output'}")
    for line in listing.splitlines():
        origin, _, setting = line.partition("\t")
        key, _, value = setting.partition(" ")
        lowered = key.lower()
        if lowered.startswith("push."):
            if lowered == "push.followtags" and value.strip().lower() in TRUE_VALUES:
                fail(f"forbidden push setting: {key}={value} ({origin})")
            print(f"NOTE push setting present: {key}={value} ({origin})")
            continue
        fail(f"forbidden push setting: {key}={value} ({origin})")


def check_bindings(repo, block):
    """The recorded commit exists, has the recorded tree, its first parent is
    the recorded base, and the remote-tracking branch is at that base."""
    commit = block["commit"]

    def resolve(what, ref):
        code, out, err = git_query(repo, "rev-parse", "--verify", "-q", ref)
        if code != 0:
            fail(f"outgoing binding mismatch: {what} {ref} not found (exit {code}): "
                 f"{err.strip() or out.strip() or 'no output'}")
        return out.strip()

    if resolve("commit", f"{commit}^{{commit}}") != commit:
        fail("outgoing binding mismatch: commit not in repository")
    tree = git(repo, "rev-parse", f"{commit}^{{tree}}").strip()
    if tree != block["tree"]:
        fail(f"outgoing binding mismatch: tree is {tree}")
    parent = resolve("parent of", f"{commit}^")
    if parent != block["base"]:
        fail(f"outgoing binding mismatch: parent {parent} is not the base")
    tracking = f"refs/remotes/{block['remote']}/{block['branch']}"
    remote_head = resolve("remote-tracking branch", tracking)
    if remote_head != block["base"]:
        fail(f"outgoing binding mismatch: {tracking} is {remote_head}, not the base")


def dry_run(repo, block):
    """One fast-forward update of the approved branch and nothing else."""
    result = subprocess.run(
        ["git", "-C", str(repo), "-c", "push.followTags=false", "push", "--no-follow-tags",
         "--dry-run", "--porcelain", block["remote"], block["refspec"]],
        text=True, capture_output=True)
    # Git's own result first: output that looks right means nothing after a failed command.
    if result.returncode != 0:
        fail(f"dry run failed (exit {result.returncode}): "
             f"{result.stderr.strip() or result.stdout.strip() or 'no output'}")
    updates = [line for line in result.stdout.splitlines() if "\t" in line]
    if not updates:
        fail(f"dry run refused: {result.stderr.strip() or result.stdout.strip() or 'no output'}")
    if len(updates) != 1:
        fail(f"dry run would update more than the approved branch: {updates}")
    flag, _, rest = updates[0].partition("\t")
    refs, _, summary = rest.partition("\t")
    if flag != " " or not refs.endswith(f":refs/heads/{block['branch']}"):
        fail(f"dry run refused: {flag.strip() or 'flag'} {refs} {summary}".rstrip())
    return summary or "fast-forward"


def check_command(args):
    blocks = parse_outgoing(args.outgoing)
    block = select_block(blocks, args.repository)
    if not args.repo.is_dir():
        fail("repository directory not found")
    if args.fetch:
        git(args.repo, "fetch", block["remote"])
    check_destinations(args.repo, block)
    check_bindings(args.repo, block)
    if args.dry_run:
        print(f"DRY RUN one ref update: {dry_run(args.repo, block)}")
    print(f"PRECHECK OK {block['repository']} {block['commit']} -> {block['remote']} {block['refspec']}")
    print(f"git -C {args.repo} -c push.followTags=false push --no-follow-tags "
          f"{block['remote']} {block['refspec']}")


def vet_command(words):
    """Accept only the documented push command shape and refuse everything
    else. Git accepts abbreviated options (--ta for --tags), combined short
    options (-vf) and environment-driven settings (--config-env=...), so a
    list of forbidden words cannot be complete; the shape is a whitelist,
    in this exact order:
    [git] [-C <path>] [-c push.followTags=false] push <options> <remote> <refspec>
    where options are drawn only from --no-follow-tags (required), --dry-run
    and --porcelain, each at most once and all before the operands; the
    remote is a configured name matching [A-Za-z0-9_][A-Za-z0-9._-]* and
    not . or ..; the refspec is <40hex>:refs/heads/<branch>. The validator
    checks the command text submitted to it; passing authorizes nothing."""
    words = words[1:] if words[:1] == ["--"] else list(words)
    if not words:
        fail("no push command given")
    index = 0
    if words[index] == "git":
        index += 1
    if index + 1 < len(words) and words[index] == "-C":
        if words[index + 1].startswith("-"):
            fail(f"forbidden push form: -C {words[index + 1]}")
        index += 2
    if index + 1 < len(words) and words[index] == "-c":
        key, separator, value = words[index + 1].partition("=")
        if not (key.lower() == "push.followtags" and separator and value == "false"):
            fail(f"forbidden push form: -c {words[index + 1]}")
        index += 2
    if index >= len(words) or words[index] != "push":
        fail("push command must have the documented shape: expected push, found "
             f"{words[index] if index < len(words) else 'nothing'}")
    index += 1
    allowed = ("--no-follow-tags", "--dry-run", "--porcelain")
    options, operands = [], []
    for word in words[index:]:
        if word in allowed:
            if operands:
                fail("push command must have the documented shape: options come before the operands")
            if word in options:
                fail(f"push command must have the documented shape: repeated option {word}")
            options.append(word)
        elif word.startswith("-"):
            fail(f"forbidden push form: {word}")
        else:
            operands.append(word)
    if "--no-follow-tags" not in options:
        fail("push command must disable tag following")
    if len(operands) != 2:
        fail("push command must have the documented shape: one remote name and exactly one "
             f"refspec, found {len(operands)} operands")
    remote, spec = operands
    if remote in (".", "..") or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9._-]*", remote):
        fail(f"forbidden push form: the remote must be a configured remote name: {remote}")
    if not re.fullmatch(rf"{COMMIT}:refs/heads/{BRANCH}", spec):
        fail(f"forbidden push form: {spec}")
    print("VET OK")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="compare OUTGOING.txt with the live repository")
    check.add_argument("repo", type=Path)
    check.add_argument("outgoing", type=Path)
    check.add_argument("--repository", help="block to check when the record names several")
    check.add_argument("--fetch", action="store_true", help="git fetch the remote first (network)")
    check.add_argument("--dry-run", action="store_true",
                       help="git push --dry-run and require one fast-forward of the branch (network)")
    vet = commands.add_parser("vet", help="refuse forbidden forms in a proposed push command")
    vet.add_argument("words", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    try:
        if args.command == "check":
            check_command(args)
        else:
            vet_command(args.words)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
