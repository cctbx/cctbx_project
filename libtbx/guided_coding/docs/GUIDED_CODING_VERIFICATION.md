# GuidedCoding r10 rev11 candidate verification

This describes a proposed kit. It has not been installed in the Developer's
Claude Code Desktop app, integrated into `cctbx_project`, or approved for
Worker use.

## Release and source checks

From the canonical `libtbx/guided_coding` root, first compare the downloaded
archive with the independently supplied SHA-256 in the review cover. Then run:

```bash
shasum -a 256 -c SOURCE_MANIFEST.sha256
python3 -I -B payload/tools/screen_check.py verify-source .
cat payload/RELEASE
```

Require every listed hash to pass, `VERIFIED complete source`, and the
`r10 rev11` release label. The first command establishes the listed bytes of
the checker before running it. The second compares the entire directory
with the manifest, including unlisted files, and refuses links, hardlinks,
missing files and changed bytes. `-I` prevents a module in the tools
folder or working directory from shadowing Python's standard library;
`-B` avoids writing bytecode into the source tree. Repeat both checks before
running central tools, integration, or publication. Run the tools with
`python3 -I -B`. These checks establish internal consistency with the
supplied manifest; they do not authenticate the supplier. A separately
stated archive hash and an independent review are still needed.

## CLI version gate and rev11 registration

After verifying the source, and before creating the personal skill link or
using `/gc setup` or `/gc <task>`, run:

```bash
python3 -I -B payload/tools/screen_check.py check-claude-version
```

The command looks up `claude` on this shell's PATH, calls `--version`, and
requires a readable version at least 2.1.281. It refuses missing,
unrecognized, failing, and older executables without installing updates.
The cutoff accounts for the rev8 Mac trial's 2.1.268 `AGENTS.md` loading
failure and [Claude Code's documented limitations before 2.1.281](https://code.claude.com/docs/en/memory#when-agentsmd-support-is-unavailable).
Tests use disposable
fake executables: 2.1.268 and 2.1.280 fail; 2.1.281, the Developer's
reported 2.1.284, and 2.2.0 pass. They also cover missing, malformed and
nonzero results. These are author-side tests on Linux, not a live Mac
installation or proof that Claude Code Desktop loads `AGENTS.md`.
The gate trusts the stdout of whichever `claude` executable is found on
PATH. A misleading executable can report a passing version. Inspect the
reported executable path in the live setup; this gate does not authenticate
it or the Desktop app.
The rev9 user-guide command had a defect: when `skills/guided_coding` already
existed as a directory or directory symlink, `ln -s` created a nested link
inside it and returned success. With an existing link to the source, that
nested link invalidated the complete-source inventory. This was reproduced
in a disposable setup. Rev10 replaces the shell `ln` with
`screen_check.py register-skill .`, which rechecks source and CLI, creates
the exact link exclusively, and refuses an occupied path without writing a
nested link. Its unit cases cover fresh registration, a second run, a link
to another package, a real directory, a dangling link, a symlinked skills
parent, and an older CLI. The rev10 author-side Linux suite ran **41 tests:
40 passed, one Mac case-insensitive test skipped**. Its reviewer found a
further case: an absolute configuration path with a missing prefix followed
by `..` let recursive directory creation change the effective path after
preflight. It could create an empty directory in the source before refusing
an occupied destination, or write through a hidden symlinked parent. Rev11
rejects any `..` component before directory creation. Its new unit case
snapshots source files and directories and the other directory across fresh,
occupied, configuration-link and skills-link endpoints. The rev11 author-side
Linux suite ran **42 tests: 41 passed, one Mac case-insensitive test skipped**.
No real Claude Code CLI or Desktop run was made for this rev11 candidate.
Read-only `/gc help` and `/gc status`, and the personal-link removal
`/gc uninstall`, remain usable if the version gate fails. The one-time
conversational setup is instructed to run the
same command before it creates a link. The registration helper enforces its
own checks even when called directly; there is still no per-project installer.

## Rev7 review findings addressed in rev8

| Reviewer probe | Rev7 result | Rev8 check and result |
| --- | --- | --- |
| Add unlisted `payload/tools/hashlib.py` | Listed-hash command passed and the module executed on `screen_check.py --help` | Complete inventory refuses the extra file before guided work. A disposable CLI test confirms neither tool's direct `--help` loads it. |
| Put companions in `packet/companions/` | Bundle helper accepted them in both domains | Bundle helper refuses before writing; packet verification still passes. It also refuses reserved companion filenames already included in a packet, even if another companion directory is outside. |
| Send source archive without release proof | Reviewer could not establish gate provenance or packet identity | Rev8 review handoff includes source archive, stated archive SHA-256, complete diff from rev7, frozen packet, separate companions, and author-side criterion/runner record. This is evidence to review, not an approval. |

The rev8 suite reported 36 tests, OK, one skipped on Linux with Python 3.12.
The skipped case needs a case-insensitive Mac volume. The new tests cover
clean and changed source, extra code, extra ordinary file and symlink, direct
CLI loading, and nested companions. An extraction check compares the source
archive's files to the candidate and reruns the source verifier. The review
packet records the exact commands, outputs and exit statuses. Author-side
checks are not an independent audit.

## Existing behavior and limits

The personal `~/.claude/skills/guided_coding` symlink is intended to expose
one central, manual `/gc` skill to local projects. A target's current
instructions must adopt the exact Developer–Guide Contract identity and
identify its own project method before GuidedCoding governs that target.
`/gc help`, `status`, `setup` and `uninstall` are control requests; setup and
adoption are conversational instructions, not an executable installer.
An ordinary task in a fresh conversation remains ordinary after older
always-on instructions are retired. These are design claims checked against
the files; they are not live Desktop observations.

The screen checker verifies a frozen evidence packet's listed bytes and
exact inventory. For that packet's identity, calculate SHA-256 of its
`MANIFEST.sha256`; `screen_check.py verify` does not print that hash. The
bundle helper requires its output in the packet's own parent directory and
companions outside the packet. It does not authenticate authors or test
truth. Source verification and evidence verification are distinct.

Before Worker use on the intended Mac, independently test `/skills`, `/gc`,
`/guided_coding`, an ordinary fresh conversation, target selection, adopted
and unadopted projects, a case-insensitive path probe, and the relevant
project instruction loading (especially an `AGENTS.md` project). The
central source and its complete release packet need independent review and
explicit adoption. No live Mac, server, or Claude Code run is claimed here.

The earlier r09 documentation review and r10 rev3–rev10 review findings
informed this candidate. They do not approve rev11. The kit still has 24
regular source files including its manifest; it does not add a per-project
installer or another permanent layer of prompts.
