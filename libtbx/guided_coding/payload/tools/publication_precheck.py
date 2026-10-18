#!/usr/bin/env python3
"""Read-only Git checks before an authorized push: effective destinations,
forbidden push settings, outgoing bindings, an optional dry run, and the
form of a proposed push command. It never pushes; it fetches only with
--fetch. Passing establishes that the recorded proposal matches the live
repository, not that the Developer has decided to publish.

`outgoing` (checkpoint 3) writes the strict one-commit OUTGOING.txt block
for one repository from the live repository, so it is generated, not
typed; it applies the same checks to what it wrote.

`confirm` (revision 17) reads the remote branch and its tag state with
`git ls-remote`, before and after a push the Developer runs himself, for
projects whose rules forbid the agent any `git push`."""

import sys
# Revision 18: the kit needs Python 3.10 or newer; an older Python refuses here, before
# anything else runs, with exit 19 (the same message and code in every tool).
if sys.version_info < (3, 10):
    sys.stderr.write("REFUSED: requires Python 3.10 or newer; this is Python %d.%d.%d (%s)\n"
                     % (tuple(sys.version_info[:3]) + (sys.executable,)))
    sys.exit(19)
# When run as a script, keep this directory out of module resolution before
# importing the standard library (the documented -I invocation also does).
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
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


def check_destinations(repo, block, notes=None):
    """Exactly one effective push URL, equal to the record, and no setting
    that rewrites, mirrors, adds refspecs or follows tags. Notes about
    harmless push settings go to `notes` (standard output by default;
    `outgoing` sends them to standard error so its output is only the block)."""
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
            print(f"NOTE push setting present: {key}={value} ({origin})", file=notes or sys.stdout)
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


def generate_outgoing(args):
    """Write the OUTGOING.txt block for one commit whose parent is the remote
    tip: repository, remote, remote_url, base, commit, tree, refspec, in the
    strict form parse_outgoing accepts. Refused when the working tree is
    dirty, the remote is missing or has no single push URL, the commit is a
    merge, its parent is not the (fetched, with --fetch) remote-tracking tip,
    or it is a chain of several commits on that tip. Multi-commit ranges are
    not supported (deferred). The block is printed, or written to a new file
    with --output (an existing file is refused); before that it is parsed
    back with parse_outgoing and checked against the repository exactly as
    `check` does (destinations, push settings, bindings)."""
    repo = args.repo
    if not repo.is_dir():
        fail("repository directory not found")
    top = Path(git(repo, "rev-parse", "--show-toplevel").strip())
    name = args.repository or top.name
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9._-]*", name):
        fail(f"repository name must be a plain name: {name}")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9._-]*", args.remote) or args.remote in (".", ".."):
        fail(f"the remote must be a configured remote name: {args.remote}")
    if not re.fullmatch(BRANCH, args.branch) or args.branch.startswith("-"):
        fail(f"invalid branch name: {args.branch}")
    dirty = git(repo, "status", "--porcelain", "--untracked-files=normal")
    if dirty.strip():
        fail("working tree is not clean; commit or remove the changes first:\n" + dirty.rstrip())
    code, out, err = git_query(repo, "remote", "get-url", "--all", "--push", args.remote)
    if code != 0:
        fail(f"remote not found: {args.remote} ({err.strip() or out.strip() or 'no output'})")
    urls = out.splitlines()
    if len(urls) != 1 or not urls[0].strip():
        fail(f"remote {args.remote} has {len(urls)} push URLs; exactly one is required: "
             f"{' '.join(urls) or 'none'}")
    if args.fetch:
        git(repo, "fetch", args.remote)
    tracking = f"refs/remotes/{args.remote}/{args.branch}"
    code, out, err = git_query(repo, "rev-parse", "--verify", "-q", tracking + "^{commit}")
    if code != 0:
        fail(f"remote-tracking branch {tracking} not found; fetch it first (--fetch)")
    tip = out.strip()
    code, out, err = git_query(repo, "rev-parse", "--verify", "-q", args.commit + "^{commit}")
    if code != 0:
        fail(f"commit not found: {args.commit}")
    commit = out.strip()
    parents = git(repo, "rev-list", "--parents", "-n", "1", commit).split()[1:]
    if len(parents) != 1:
        fail(f"commit {commit} has {len(parents)} parents; one commit on the remote tip is required")
    if parents[0] != tip:
        code, _, _ = git_query(repo, "merge-base", "--is-ancestor", tip, parents[0])
        if code == 0:
            count = git(repo, "rev-list", "--count", f"{tip}..{commit}").strip()
            fail(f"commit {commit} is a chain of {count} commits on {tracking}; one commit is "
                 "required (multi-commit ranges are not supported)")
        fail(f"parent {parents[0]} of commit {commit} is not the remote tip {tracking} {tip}")
    tree = git(repo, "rev-parse", f"{commit}^{{tree}}").strip()
    text = (f"repository: {name}\nremote: {args.remote}\nremote_url: {urls[0]}\n"
            f"base: {tip}\ncommit: {commit}\ntree: {tree}\n"
            f"refspec: {commit}:refs/heads/{args.branch}\n")
    with tempfile.TemporaryDirectory() as scratch:
        draft = Path(scratch) / "OUTGOING.txt"
        draft.write_text(text, encoding="utf-8")
        (block,) = parse_outgoing(draft)
    check_destinations(repo, block, notes=sys.stderr)
    check_bindings(repo, block)
    if args.output:
        with open(args.output, "x", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
        print(f"OUTGOING written: {args.output}")
    else:
        sys.stdout.write(text)
    print(f"OUTGOING OK {name} {commit} on {tracking} {tip}", file=sys.stderr)


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


def no_prompt_environment(repo):
    """Revision 17 (repair 1, reading 1 finding 3): the environment for `confirm`'s git calls.
    GIT_TERMINAL_PROMPT=0 stops git asking for an HTTPS user name or password; ssh gets
    `-o BatchMode=yes`, so it never asks for a key passphrase, a password or a host-key
    confirmation and fails instead (an agent's keys still work). It is appended to the ssh
    command git would otherwise use, keeping that command's options: GIT_SSH_COMMAND if set,
    else the repository's core.sshCommand, else GIT_SSH (a program name; one that is not
    OpenSSH may refuse -o, which fails closed), else `ssh`."""
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")
    command = (env.get("GIT_SSH_COMMAND") or "").strip()
    if not command:
        code, out, _ = git_query(repo, "config", "--get", "core.sshCommand")
        command = out.strip() if code == 0 else ""
    if not command and env.get("GIT_SSH"):
        command = shlex.quote(env["GIT_SSH"])
    env["GIT_SSH_COMMAND"] = (command or "ssh") + " -o BatchMode=yes"
    return env


def remote_heads(repo, block):
    """Revision 17: (branch head on the remote or None, tag state), read with `git ls-remote`
    from the recorded push URL (read-only; network; with no_prompt_environment, so neither
    git nor ssh waits for typed input). The tag state is `<count>/<sha256 of the sorted tag
    ref lines>`, so a tag added and another removed is seen as well as a changed count. No
    other ref is read: a branch pushed besides the outgoing one is not seen."""
    url, branch = block["remote_url"], block["branch"]
    env = no_prompt_environment(repo)

    def ls_remote(*words):
        result = subprocess.run(["git", "-C", str(repo), "ls-remote", *words, url],
                                text=True, capture_output=True, env=env)
        if result.returncode != 0:
            fail(f"git ls-remote {' '.join(words)} failed (exit {result.returncode}): "
                 f"{result.stderr.strip() or result.stdout.strip() or 'no output'}")
        return result.stdout.splitlines()
    lines = ls_remote("--heads")
    heads = [line.split("\t")[0] for line in lines
             if line.split("\t")[1:] == [f"refs/heads/{branch}"]]
    if len(heads) > 1:
        fail(f"remote lists refs/heads/{branch} more than once")
    tags = sorted(line for line in ls_remote("--tags") if line.strip())
    digest = hashlib.sha256("".join(line + "\n" for line in tags).encode()).hexdigest()
    count = sum(1 for line in tags if not line.endswith("^{}"))
    return (heads[0] if heads else None), f"{count}/{digest}"


def confirm_command(args):
    """Revision 17, the no-push publication variant (a project whose rules forbid the agent
    any `git push`, read as including dry runs): the agent runs `check` without --dry-run and
    `vet`s the exact commands, the Developer runs them himself, and this command then confirms
    the remote. Read-only: it reads the remote with `git ls-remote` on the recorded push URL
    (no HTTPS prompt, ssh in BatchMode: no_prompt_environment) and never pushes or fetches. With --before it confirms that
    the remote branch is still at the recorded base (nothing pushed yet) and prints the
    remote's tag state; without it, --tags STATE (the state printed before the push) is
    required, and it confirms that the branch is now at the outgoing commit and the tag
    state is unchanged (refused when the branch is still at the base, or anywhere else, or
    any tag was added, removed or moved)."""
    blocks = parse_outgoing(args.outgoing)
    block = select_block(blocks, args.repository)
    if not args.repo.is_dir():
        fail("repository directory not found")
    if not args.before and args.tags is None:
        fail("after the push, give --tags with the tag state that confirm --before printed")
    check_destinations(args.repo, block)
    head, tags = remote_heads(args.repo, block)
    branch = f"refs/heads/{block['branch']}"
    if args.before:
        if head != block["base"]:
            fail(f"remote {branch} is {head or 'absent'}, not the recorded base {block['base']}")
    elif head != block["commit"]:
        if head == block["base"]:
            fail(f"remote {branch} is still at the base {block['base']}: the push has not "
                 "happened")
        fail(f"remote {branch} is {head or 'absent'}, not the outgoing commit {block['commit']}")
    if args.tags is not None and tags != args.tags:
        fail(f"remote tag state is {tags}, not {args.tags}: a tag was added, removed or moved")
    state = "still at the base" if args.before else "at the outgoing commit"
    print(f"REMOTE {block['repository']} {block['remote']} {branch} {head} ({state})")
    print(f"REMOTE TAGS {tags}")
    print("CONFIRMED before the push" if args.before else
          f"CONFIRMED {block['repository']} {block['commit']} is {block['remote']} {branch}")


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
    out = commands.add_parser("outgoing", help="write the OUTGOING.txt block for one commit")
    out.add_argument("repo", type=Path)
    out.add_argument("--commit", required=True, help="the one outgoing commit")
    out.add_argument("--remote", default="origin")
    out.add_argument("--branch", default="master")
    out.add_argument("--repository", help="block name (default: the repository directory name)")
    out.add_argument("--fetch", action="store_true", help="git fetch the remote first (network)")
    out.add_argument("--output", type=Path, help="new file to write (default: standard output)")
    confirm = commands.add_parser("confirm", help="revision 17: read the remote branch and tag "
                                  "count after (or, with --before, before) the Developer's push")
    confirm.add_argument("repo", type=Path)
    confirm.add_argument("outgoing", type=Path)
    confirm.add_argument("--repository", help="block to confirm when the record names several")
    confirm.add_argument("--before", action="store_true",
                         help="confirm the branch is still at the base (before the push)")
    confirm.add_argument("--tags", help="after the push (required): the tag state COUNT/SHA256 "
                         "that confirm --before printed")
    vet = commands.add_parser("vet", help="refuse forbidden forms in a proposed push command")
    vet.add_argument("words", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command in ("check", "outgoing", "confirm") and shutil.which("git") is None:
        # Revision 16: these read the live repository with git; vet needs no git.
        print(f"ERROR: requires git on PATH: `{args.command}` reads the repository with git "
              "(`vet` needs no git)", file=sys.stderr)
        return 2
    try:
        if args.command == "check":
            check_command(args)
        elif args.command == "outgoing":
            generate_outgoing(args)
        elif args.command == "confirm":
            confirm_command(args)
        else:
            vet_command(args.words)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
