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
projects whose rules forbid the agent any `git push`.

`carry-over` (revision 19) compares the OUTGOING.txt the outside reader
read with the OUTGOING.txt of a restaged packet and writes CARRYOVER.txt:
the same changed files with the same blobs, a byte-identical diff, and an
upstream range that touches none of them."""

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
import ast
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
        return check_record_remote(blocks[0])
    for block in blocks:
        if block["repository"] == name:
            return check_record_remote(block)
    fail(f"repository not in outgoing record: {name}")


def record_remote_problem(block):
    """The record's remote name and URL are passed to git (`fetch`, `ls-remote`), so a value
    starting with "-" would be read as a git option (for example --upload-pack, which runs a
    command). They are held to the rule `outgoing` applies when it writes the record; the
    reason a block fails it, or None."""
    remote, url = block["remote"], block["remote_url"]
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9._-]*", remote) or remote in (".", ".."):
        return f"the outgoing record's remote must be a configured remote name: {remote}"
    if url.startswith("-") or any(c.isspace() for c in url):
        return f"the outgoing record's remote_url is not a plain URL: {url}"
    # A remote-helper URL (<transport>::<address>) runs a helper program; ext:: runs any
    # command when git allows it. A publication remote never needs one.
    if re.match(r"[A-Za-z][A-Za-z0-9+.-]*::", url):
        return f"the outgoing record's remote_url names a remote helper (<transport>::), not a plain URL: {url}"
    return None


def check_record_remote(block):
    problem = record_remote_problem(block)
    if problem:
        fail(problem)
    return block


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


def check_packet_identity(args):
    """Repair 3 (review 3, finding 1): `check` is tied to the packet the Developer approved.
    --packet DIR and --identity SHA256 (the packet identity on the PUBLICATION screen) are
    required; DIR is the directory of the OUTGOING.txt given, its MANIFEST.sha256 must verify,
    and its SHA-256 must be that identity. A re-frozen packet, or a copy of OUTGOING.txt
    elsewhere, is refused."""
    if args.packet is None or args.identity is None:
        fail("check needs --packet DIR and --identity SHA256: the frozen publication packet that "
             "holds this OUTGOING.txt and its identity as shown on the PUBLICATION screen")
    if not re.fullmatch(_checker.SHA, args.identity):
        fail("--identity must be the 64-hex packet identity (SHA-256 of the packet's MANIFEST.sha256)")
    packet = args.packet
    if packet.is_symlink() or not packet.is_dir() or args.outgoing.resolve() != (packet / "OUTGOING.txt").resolve():
        fail("--packet must be the packet directory that holds this OUTGOING.txt")
    identity = _checker.verify(packet)
    if identity != args.identity:
        fail(f"the packet's identity is {identity}, not the approved {args.identity}: it is not the "
             "packet on the PUBLICATION screen")
    print(f"PACKET IDENTITY {identity} (as approved)")


def check_command(args):
    blocks = parse_outgoing(args.outgoing)
    block = select_block(blocks, args.repository)
    check_packet_identity(args)
    if not args.repo.is_dir():
        fail("repository directory not found")
    if args.fetch:
        git(args.repo, "fetch", block["remote"])
    check_destinations(args.repo, block)
    check_bindings(args.repo, block)
    check_carryover(args, block)
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


# Revision 19: the carry-over check. The header, the OK line, the flags and the consistency
# rules are screen_check.py's, which reads CARRYOVER.txt.
CARRYOVER_OK, CARRYOVER_HEADER = _checker.CARRYOVER_OK, _checker.CARRYOVER_HEADER
upstream_flags, FLAG_KINDS = _checker.upstream_flags, _checker.FLAG_KINDS


class CarryOverRefused(Exception):
    """One carry-over condition failed; the reason is recorded in CARRYOVER.txt."""


def raw_git(repo, *args):
    """A read-only git query returning bytes (diffs are compared byte for byte)."""
    result = subprocess.run(["git", "-C", str(repo), *args], capture_output=True)
    if result.returncode != 0:
        fail(f"git {' '.join(args)} failed (exit {result.returncode}): "
             f"{result.stderr.decode('utf-8', 'replace').strip() or 'no output'}")
    return result.stdout


# Fixed diff options, so that no user or repository setting (renames, colour, external diff
# drivers, text conversion, quoting) changes the bytes that are compared.
DIFF = ("-c", "core.quotePath=true", "diff", "--no-color", "--no-ext-diff", "--no-textconv",
        "--binary", "--full-index", "--no-renames")


def changed_entries(repo, base, commit):
    """{path: (mode, blob) or None for a deletion} of the paths `commit` changes against `base`."""
    out = raw_git(repo, "-c", "core.quotePath=true", "diff-tree", "-r", "-z", "--no-renames",
                  base, commit)
    fields = out.split(b"\0")
    entries = {}
    index = 0
    while index + 1 < len(fields) and fields[index]:
        meta, path = fields[index].decode(), fields[index + 1].decode("utf-8", "surrogateescape")
        index += 2
        _, new_mode, _, new_blob, status = meta.lstrip(":").split(" ")
        entries[path] = None if status == "D" else (new_mode, new_blob)
    return entries


def module_paths(module, files):
    """In-repository files a module name can be: a/b/c.py or a/b/c/__init__.py, at the
    repository root (cctbx_project and phenix keep their packages there)."""
    if not module:
        return set()
    base = module.replace(".", "/")
    return {path for path in (base + ".py", base + "/__init__.py") if path in files}


def imports_and_names(source, path, files):
    """Revision 19 (plan revision 6): read one changed Python file with ast, never running
    it, and return (in-repository files it imports, files it names by a literal path).
    `import a.b` gives a/b(.py|/__init__.py) and the packages a and a.b; `from a import b`
    gives a and a/b; a relative import is resolved from the file's own package. A literal
    string naming an existing file of the repository, relative to the file's directory or to
    the repository root, counts as named. Dynamic imports (importlib, __import__, exec) and
    computed paths are not seen. Returns None when the file cannot be parsed."""
    try:
        tree = ast.parse(source, filename=path)
    except (SyntaxError, ValueError):
        return None
    package = path.split("/")[:-1]
    imported, named = set(), set()

    def add(module):
        parts = module.split(".")
        for index in range(1, len(parts) + 1):
            imported.update(module_paths(".".join(parts[:index]), files))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if node.level - 1 > len(package):
                    continue
                prefix = package[:len(package) - (node.level - 1)]
                base = ".".join(prefix + ([node.module] if node.module else []))
            else:
                base = node.module or ""
            if base:
                add(base)
            for alias in node.names:
                if alias.name != "*":
                    imported.update(module_paths(f"{base}.{alias.name}" if base else alias.name, files))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            text = node.value.strip()
            if not text or "\n" in text or len(text) > 300 or text.startswith("/"):
                continue
            for candidate in ("/".join(package + [text]), text):
                normal = [part for part in candidate.split("/") if part not in ("", ".")]
                while ".." in normal:
                    index = normal.index("..")
                    if index == 0:
                        normal = None
                        break
                    del normal[index - 1:index + 1]
                if normal and "/".join(normal) in files:
                    named.add("/".join(normal))
    imported.discard(path)
    named.discard(path)
    return imported, named


def resolve_commit(repo, ref, what):
    """Repair 1 (review finding 3): the refusal names no path, so a recompute from another
    spelling of the repository path gives the same text."""
    code, out, err = git_query(repo, "rev-parse", "--verify", "-q", ref + "^{commit}")
    if code != 0 or out.strip() != ref:
        raise CarryOverRefused(f"{what} {ref} is not a commit in the repository")
    return out.strip()


def commit_message(repo, commit):
    """The commit's message bytes exactly (what is published), without author, committer
    and dates."""
    raw = raw_git(repo, "cat-file", "commit", commit)
    header, separator, message = raw.partition(b"\n\n")
    return message if separator else b""


def carry_block(repo, old, new, lines):
    """The carry-over checks for one repository block, appending the report to `lines`.

    Repair 2 (review 2, finding 1): first every check that cannot be waived (the same remote,
    destination and branch; every commit resolves, with its recorded tree, as one commit on its
    recorded base; the old base is an ancestor of the new base; the same set of changed paths),
    any of which raises CarryOverRefused at once. Then every content comparison is made and
    recorded (blobs, the diff, upstream commits touching changed paths, the commit message), and
    all content differences together end the section as `difference (NAME): ...` lines and a
    `repository NAME: REFUSED: N content difference(s)` result, which a Developer may waive.
    Returns the list of content differences (empty when the block passes)."""
    name = new["repository"]
    lines += [f"== repository: {name}",
              f"reviewed: base {old['base']} commit {old['commit']} tree {old['tree']}",
              f"outgoing: base {new['base']} commit {new['commit']} tree {new['tree']}"]
    for key in ("remote", "remote_url", "branch"):
        if old[key] != new[key]:
            raise CarryOverRefused(f"the {key} differs between the reviewed and the outgoing record")
    for block, what in ((old, "reviewed"), (new, "outgoing")):
        resolve_commit(repo, block["base"], f"{what} base")
        resolve_commit(repo, block["commit"], f"{what} commit")
        if git(repo, "rev-parse", block["commit"] + "^{tree}").strip() != block["tree"]:
            raise CarryOverRefused(f"{what} commit {block['commit']} does not have the recorded "
                                   f"tree {block['tree']}")
        parents = git(repo, "rev-list", "--parents", "-n", "1", block["commit"]).split()[1:]
        if parents != [block["base"]]:
            raise CarryOverRefused(f"{what} commit {block['commit']} is not one commit on its "
                                   f"recorded base {block['base']}")
    code, _, err = git_query(repo, "merge-base", "--is-ancestor", old["base"], new["base"])
    if code != 0:
        raise CarryOverRefused(f"the reviewed base {old['base']} is not an ancestor of the "
                               f"outgoing base {new['base']}"
                               + (f" ({err.strip()})" if code != 1 and err.strip() else ""))
    lines.append("reviewed base is an ancestor of the outgoing base: yes")
    reviewed = changed_entries(repo, old["base"], old["commit"])
    outgoing = changed_entries(repo, new["base"], new["commit"])
    if not outgoing:
        raise CarryOverRefused("the outgoing commit changes no file")
    if set(reviewed) != set(outgoing):
        only_r = sorted(set(reviewed) - set(outgoing))
        only_o = sorted(set(outgoing) - set(reviewed))
        raise CarryOverRefused(f"the commits change different paths "
                               f"(reviewed only: {' '.join(only_r) or 'none'}; "
                               f"outgoing only: {' '.join(only_o) or 'none'})")
    lines.append(f"changed paths: {len(outgoing)} (the same set in both commits)")
    # Content comparisons: each is made and recorded, whatever the others find.
    differences = []
    differing = []
    for path in sorted(outgoing):
        before, after = reviewed[path], outgoing[path]
        same = before == after
        describe = (lambda entry: "deleted" if entry is None else f"{entry[0]} {entry[1]}")
        lines.append(f"blob {'same' if same else 'DIFFERS'} reviewed {describe(before)} "
                     f"outgoing {describe(after)} {path}")
        if not same:
            differing.append(path)
    if differing:
        differences.append(_checker.DIFF_BLOBS + " ".join(differing))
    reviewed_diff = raw_git(repo, *DIFF, old["base"], old["commit"])
    outgoing_diff = raw_git(repo, *DIFF, new["base"], new["commit"])
    lines.append(f"diff sha256 reviewed {hashlib.sha256(reviewed_diff).hexdigest()} "
                 f"outgoing {hashlib.sha256(outgoing_diff).hexdigest()} "
                 f"({len(outgoing_diff)} bytes)")
    if reviewed_diff != outgoing_diff:
        lines.append("diff: DIFFERS")
        differences.append(_checker.DIFF_DIFF)
    else:
        lines.append("diff: byte-identical")
    upstream = git(repo, "rev-list", "--reverse", f"{old['base']}..{new['base']}").split()
    lines.append(f"upstream commits: {len(upstream)}")
    touched_changed, flagged = [], {kind: set() for kind in FLAG_KINDS}
    upstream_paths = {}
    for commit in upstream:
        # -m: a merge is compared with each parent, so a path brought in by either side counts.
        out = raw_git(repo, "-c", "core.quotePath=true", "diff-tree", "-r", "-m", "-z",
                      "--no-renames", "--no-commit-id", "--name-status", commit)
        fields = [field.decode("utf-8", "surrogateescape") for field in out.split(b"\0") if field]
        files = sorted(set(zip(fields[0::2], fields[1::2])), key=lambda item: (item[1], item[0]))
        lines.append(f"upstream commit {commit}: {len(files)} file(s)")
        upstream_paths[commit] = [path for _, path in files]
        for status, path in files:
            flags = upstream_flags(path)
            lines.append(f"  {status} {path}" + (f" [{','.join(flags)}]" if flags else ""))
            for kind in flags:
                flagged[kind].add(path)
            if path in outgoing and f"{commit[:12]} {path}" not in touched_changed:
                touched_changed.append(f"{commit[:12]} {path}")
    # The net base-to-base change is covered by the per-commit rows; it is kept as a guard.
    net = set(git(repo, "-c", "core.quotePath=true", "diff", "--name-only", "--no-renames",
                  old["base"], new["base"]).splitlines())
    reported = {entry.split(" ", 1)[1] for entry in touched_changed}
    if (net & set(outgoing)) - reported:
        raise CarryOverRefused("upstream changed paths that no upstream commit lists: "
                               + " ".join(sorted((net & set(outgoing)) - reported)))
    if touched_changed:
        lines.append("upstream commits touching a changed path: " + "; ".join(touched_changed))
        differences.append(_checker.DIFF_TOUCH + "; ".join(touched_changed))
    else:
        lines.append("upstream commits touching a changed path: none")
    lines.append(_checker.FLAGGED_INTRO + ", ".join(f"{kind} {len(flagged[kind])}"
                                                     for kind in FLAG_KINDS))
    for kind in FLAG_KINDS:
        lines.append(f"flagged {kind}: {' '.join(sorted(flagged[kind])) or 'none'}")
    # Plan revision 6 (a): upstream files that a changed Python file imports or names
    # (information only; dynamic imports are not seen).
    touched = set().union(*(set(paths) for paths in upstream_paths.values())) if upstream_paths else set()
    files = set(path.decode("utf-8", "surrogateescape") for path in raw_git(
        repo, "-c", "core.quotePath=true", "ls-tree", "-r", "-z", "--name-only", new["commit"]).split(b"\0") if path)
    readable, unreadable, imported, named = 0, [], [], []
    for path in sorted(outgoing):
        if outgoing[path] is None or not path.endswith(".py"):
            continue
        source = raw_git(repo, "cat-file", "blob", outgoing[path][1])
        found = imports_and_names(source, path, files)
        if found is None:
            unreadable.append(path)
            continue
        readable += 1
        imported += [f"{target} (imported by {path})" for target in sorted(found[0] & touched)]
        named += [f"{target} (named by {path})" for target in sorted(found[1] & touched)]
    lines.append(f"changed Python files read for imports: {readable} (not parsed: "
                 f"{' '.join(unreadable) or 'none'})")
    lines.append(f"{_checker.IMPORTED_LABEL}: {'; '.join(imported) or 'none'}")
    lines.append(f"{_checker.NAMED_LABEL}: {'; '.join(named) or 'none'}")
    # Repair 1 (review finding 4): the published message is part of what was read.
    reviewed_message = commit_message(repo, old["commit"])
    outgoing_message = commit_message(repo, new["commit"])
    lines.append(f"commit message sha256 reviewed {hashlib.sha256(reviewed_message).hexdigest()} "
                 f"outgoing {hashlib.sha256(outgoing_message).hexdigest()}")
    if reviewed_message != outgoing_message:
        lines.append("commit message: DIFFERS")
        differences.append(_checker.DIFF_MESSAGE)
    else:
        lines.append("commit message: identical")
    lines += [f"difference ({name}): {text}" for text in differences]
    return differences


def repository_paths(words, names):
    """REPO arguments, in the pattern of `check --repository NAME`: NAME=PATH for each
    repository block (NAME is the block's `repository:` value). With one block, one bare PATH
    is that block's repository."""
    if len(names) == 1 and len(words) == 1 and "=" not in words[0]:
        return {next(iter(names)): Path(words[0])}
    paths = {}
    for word in words:
        name, separator, path = word.partition("=")
        if not separator or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9._-]*", name) or not path:
            raise CarryOverRefused(f"with several repository blocks give each repository as "
                                   f"NAME=PATH (NAME as in the block's repository: line): {word}")
        if name in paths:
            raise CarryOverRefused(f"repository named twice: {name}")
        paths[name] = Path(path)
    return paths


def path_free(refusal, repo):
    """A refusal's text with the repository path replaced, so that it depends only on the
    repository's content and its NAME (repair 1, review finding 3)."""
    text = " ".join(str(refusal).split())
    for spelling in {str(repo), str(Path(repo).absolute()), str(Path(repo).resolve())}:
        if spelling and spelling != ".":
            text = text.replace(spelling, "<repository>")
    return text


def carry_over_lines(repos_words, reviewed, outgoing, fetch=False):
    """Revision 19: CARRYOVER.txt in its fixed form, and whether the check passed. Every
    repository block is checked on its own, and the result is OK only when all pass."""
    lines = [CARRYOVER_HEADER]
    try:
        old_bytes, new_bytes = reviewed.read_bytes(), outgoing.read_bytes()
        lines.append(f"reviewed_outgoing_sha256: {hashlib.sha256(old_bytes).hexdigest()}")
        lines.append(f"outgoing_sha256: {hashlib.sha256(new_bytes).hexdigest()}")
        old_blocks = {block["repository"]: block for block in parse_outgoing(reviewed)}
        new_blocks = {block["repository"]: block for block in parse_outgoing(outgoing)}
        if set(old_blocks) != set(new_blocks):
            raise CarryOverRefused("the records name different repositories (reviewed: "
                                   f"{' '.join(sorted(old_blocks))}; outgoing: "
                                   f"{' '.join(sorted(new_blocks))})")
        repos = repository_paths(repos_words, new_blocks)
        missing = sorted(set(new_blocks) - set(repos))
        extra = sorted(set(repos) - set(new_blocks))
        if missing or extra:
            raise CarryOverRefused("each repository block needs exactly one REPO (missing: "
                                   f"{' '.join(missing) or 'none'}; not in the records: "
                                   f"{' '.join(extra) or 'none'})")
    except CarryOverRefused as refusal:
        return lines + [f"CARRYOVER REFUSED: {refusal}"], False
    except (ValueError, OSError) as error:
        return lines + [f"CARRYOVER REFUSED: {error}"], False
    problems = [problem for block in list(old_blocks.values()) + list(new_blocks.values())
                for problem in [record_remote_problem(block)] if problem]
    if problems:   # before any git call, as in check and confirm
        return lines + [f"CARRYOVER REFUSED: {problems[0]}"], False
    refusals = []
    for name in sorted(new_blocks):
        try:
            if not repos[name].is_dir():
                raise CarryOverRefused("repository directory not found")
            if fetch:
                git(repos[name], "fetch", new_blocks[name]["remote"])
            differences = carry_block(repos[name], old_blocks[name], new_blocks[name], lines)
            if differences:
                result = _checker.content_result(len(differences))
                lines.append(f"repository {name}: REFUSED: {result}")
                refusals.append(f"{name}: {result}")
            else:
                lines.append(f"repository {name}: OK")
        except (CarryOverRefused, ValueError, OSError) as refusal:
            if not any(line == f"== repository: {name}" for line in lines):
                lines.append(f"== repository: {name}")
            reason = path_free(refusal, repos[name])
            lines.append(f"repository {name}: REFUSED: {reason}")
            refusals.append(f"{name}: {reason}")
    if refusals:
        return lines + ["CARRYOVER REFUSED: " + "; ".join(refusals)], False
    return lines + [CARRYOVER_OK], True


def carry_over_command(args):
    lines, passed = carry_over_lines(args.repo, args.reviewed, args.outgoing, args.fetch)
    text = "".join(line + "\n" for line in lines)
    if args.expect:
        # The PUBLISH precheck: the live repositories must still give the frozen file, and it
        # must be an OK result. A result that differs or cannot be computed is refused.
        # Plan revision 6 (b): a refused result is accepted only when it is the frozen file and
        # the packet beside it holds a valid ADDENDUM_WAIVER.md (repair 2: its form and
        # identities are checked, and it must list exactly the recorded differences).
        waived = False
        if (not passed and args.expect.is_file()
                and args.expect.read_bytes() == text.encode("utf-8")
                and (args.expect.parent / _checker.ADDENDUM_WAIVER).exists()):
            try:
                validate_waiver(args.expect.parent)
                waived = True
            except ValueError as error:
                print(f"ADDENDUM_WAIVER.md refused: {error}", file=sys.stderr)
        if (not args.expect.is_file() or args.expect.read_bytes() != text.encode("utf-8")
                or (not passed and not waived)):
            sys.stdout.write(text)
            print(f"CARRYOVER DIFFERS from {args.expect}: the recomputed carry-over is not the "
                  "frozen CARRYOVER.txt (or it refused without ADDENDUM_WAIVER.md)", file=sys.stderr)
            return 1
        print(f"CARRYOVER SAME as {args.expect}" + ("" if passed else
              f" (refused; the addendum is waived in {_checker.ADDENDUM_WAIVER})"), file=sys.stderr)
        if not passed:
            sys.stdout.write(text)
            return 0
    if args.output:
        with open(args.output, "x", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
        print(f"CARRYOVER written: {args.output}", file=sys.stderr)
    sys.stdout.write(text)
    return 0 if passed else 1


def packet_has_carryover(packet):
    """Whether a directory holds, or its MANIFEST.sha256 lists, a CARRYOVER.txt, or holds
    evidence/reviewed/."""
    manifest = packet / "MANIFEST.sha256"
    listed = []
    if manifest.is_file() and not manifest.is_symlink():
        listed = [line.split("  ", 1)[-1] for line in
                  manifest.read_text(encoding="utf-8", errors="replace").splitlines()]
    return (any(name.split("/")[-1] == "CARRYOVER.txt" for name in listed)
            or (packet / "CARRYOVER.txt").exists() or (packet / "CARRYOVER.txt").is_symlink()
            or any(True for _ in packet.rglob("CARRYOVER.txt"))
            or (packet / "evidence" / "reviewed").exists())


def validate_waiver(packet):
    """Repair 2 (review 2, finding 3): the tools check ADDENDUM_WAIVER.md's fixed form and
    identities with the checker's own rules (an ordinary file, the reviewed and outgoing
    identities, exactly the recorded differences, waivable kinds only)."""
    carry = _checker.carryover_evidence(packet)
    if carry is None or not carry["refused"]:
        fail("ADDENDUM_WAIVER.md without a refused CARRYOVER.txt")
    _checker.addendum_waiver(packet, carry)


def check_carryover(args, block):
    """`check` and the carry-over (revisions 19, repairs 1 and 2).

    Without --packet: refused when the OUTGOING.txt sits in a directory that holds a
    CARRYOVER.txt, a MANIFEST.sha256 listing one, or evidence/reviewed/ (use --packet).
    With --packet DIR: DIR is the directory of OUTGOING.txt; its MANIFEST.sha256 must verify;
    when it holds CARRYOVER.txt, every repository block needs a repository (NAME=PATH with
    --repos, or the one REPO for a one-block record), the whole carry-over is recomputed and
    must equal the frozen file byte for byte, and a refused result is accepted only with a
    valid ADDENDUM_WAIVER.md. A packet that lists or holds CARRYOVER.txt anywhere else, or holds
    evidence/reviewed/, without one at its top level is refused."""
    outgoing = args.outgoing
    packet = args.packet  # bound to the approved identity by check_packet_identity()
    frozen = packet / "CARRYOVER.txt"
    if not frozen.exists() and not frozen.is_symlink():
        if packet_has_carryover(packet):
            fail("the packet lists or holds a CARRYOVER.txt (or evidence/reviewed/) but has no "
                 "CARRYOVER.txt beside OUTGOING.txt: the carry-over cannot be recomputed")
        print("PACKET VERIFIED (no carry-over)")
        return
    reviewed = packet / _checker.REVIEWED_OUTGOING
    if frozen.is_symlink() or not frozen.is_file() or reviewed.is_symlink() or not reviewed.is_file():
        fail("CARRYOVER.txt is present, but it or evidence/reviewed/OUTGOING.txt is not an "
             "ordinary file: the carry-over cannot be recomputed")
    blocks = parse_outgoing(outgoing)
    if args.repos:
        words = list(args.repos)
    elif len(blocks) == 1:
        words = [str(args.repo)]
    else:
        fail("the packet holds CARRYOVER.txt and several repository blocks: give every "
             "repository as NAME=PATH with --repos, so that all of them are recomputed")
    lines, passed = carry_over_lines(words, reviewed, outgoing)
    text = "".join(line + "\n" for line in lines)
    if frozen.read_bytes() != text.encode("utf-8"):
        fail("the carry-over recomputed for every repository differs from the frozen CARRYOVER.txt")
    if not passed:
        if not (packet / _checker.ADDENDUM_WAIVER).exists():
            fail("CARRYOVER.txt records a refusal and the packet has no ADDENDUM_WAIVER.md")
        validate_waiver(packet)
    print(f"CARRYOVER RECOMPUTED {len(blocks)} repository block(s): same as the frozen "
          "CARRYOVER.txt" + ("" if passed else " (refused; the addendum is waived)"))


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
    check.add_argument("--packet", type=Path, help="revision 19 (required): the frozen packet "
                       "directory holding OUTGOING.txt; its manifest is verified and any "
                       "carry-over is recomputed for every repository")
    check.add_argument("--identity", help="revision 19 (required): the packet identity shown on "
                       "the PUBLICATION screen; refused unless it is the SHA-256 of DIR/MANIFEST.sha256")
    check.add_argument("--repos", nargs="+", metavar="NAME=PATH", help="with --packet: every "
                       "repository block's repository, for the carry-over recompute")
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
    carry = commands.add_parser("carry-over", help="revision 19: check that a completed outside "
                                "reading carries over to a restaged commit; writes CARRYOVER.txt")
    carry.add_argument("repo", nargs="+", help="the repository: PATH when the records have one "
                       "block, otherwise NAME=PATH for every block (NAME as in its repository: line)")
    carry.add_argument("--reviewed", type=Path, required=True,
                       help="the OUTGOING.txt of the packet the outside reader read")
    carry.add_argument("--outgoing", type=Path, required=True,
                       help="the OUTGOING.txt of the restaged packet")
    carry.add_argument("--fetch", action="store_true", help="git fetch each remote first (network)")
    carry.add_argument("--output", type=Path, help="new file to write CARRYOVER.txt to (an "
                       "existing file is refused); the text is also printed")
    carry.add_argument("--expect", type=Path, help="a frozen CARRYOVER.txt that the output must "
                       "equal byte for byte (at the PUBLISH precheck)")
    vet = commands.add_parser("vet", help="refuse forbidden forms in a proposed push command")
    vet.add_argument("words", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command in ("check", "outgoing", "confirm", "carry-over") and shutil.which("git") is None:
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
        elif args.command == "carry-over":
            return carry_over_command(args)
        else:
            vet_command(args.words)
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
