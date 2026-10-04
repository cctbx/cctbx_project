# GuidedCoding — an opt-in coding procedure

GuidedCoding helps a developer direct an LLM through a bounded change:
agree on the problem and a check, build and test, inspect the evidence,
then decide whether to integrate and publish. It combines instructions
for the people and agents involved with small Python tools that check
source inventories, evidence packets and decision-screen formats.

Start a guided conversation with **`/guided_coding`**. The shorter `/gc`
remains available; the longer name makes the purpose clearer. Neither
command is a separate Claude Code permission mode.

## Status

GuidedCoding 2.0 r10 rev17 (`followups-20261004`) is a **candidate
revision of the limited opt-in pilot**. The pilot implementation
`enumcheck-20261002T194333Z` was published to `cctbx_project` on
2026-10-02 (commit `c36887c7f489018af4f91773246ecde32d1b4e24`), and its
documentation revision `docs-20261003` on 2026-10-03 (commit
`b0747a4a55f29db3abe04358480d5867e94cb792`). This revision changes the
checker, adds two read-only tools and revises the procedure texts; the
Developer–Guide Contract, the adoption rule and the setup recipe are
unchanged. Its own commit, source manifest and review identities belong
in its task and publication records, not in this file.

The pilot includes:

- One centrally stored, manually invoked skill, with help, status, setup
  and uninstall instructions.
- Project-specific setup that recovers saved settings, asks for relevant
  missing information and prepares complete changes for authorization.
- A checked setup-save sequence: preserve approved text and prior files,
  verify the recovery records, then apply and read back authorized edits.
- Source and evidence inventory checks, concise decision screens, and a
  bundle format for an Outside Reviewer.

This revision adds:

- Publication conventions without tags: the destination's convention is
  recorded with its source, a push uses one explicit `<commit>:refs/heads/<branch>`
  refspec with tag following disabled, and `publication_precheck.py`
  compares the effective push destinations, settings, outgoing bindings
  and a dry run with the frozen record before the authorized push, and
  accepts only the documented command shape when a command text is
  submitted to it (it does not intercept Git run any other way).
- Explicit waivers: a suite marked `NOT RUN` must carry the Developer's
  quoted words for that batch and scope. The checker checks the form;
  people judge the authorization. "Commit and publish" waives nothing.
- One outside reading for integration and publication of the same frozen
  packet, bound by `OUTGOING.txt` and a `Scope:` line. A changed proposal
  gets a new identity and a reassessment; the earlier verdict is not inherited.
- Publication plans that name each installation to update, with permitted
  operations, recovery boundary and stop conditions. Publication authorizes
  no installation update by itself.
- A per-job `JOB_SUMMARY.txt` and the read-only `/guided_coding history`
  listing of a project's records.
- Session naming through the client's rename tool with a manual fallback,
  recorded search roots for discovery, and a rule for denied privacy access.

The pilot has material limits:

- The save sequence depends on the Guide following it. Failed verification
  must hold the protected edits. A failure during apply or readback can
  leave **PARTIAL** state; the Guide must stop and preserve recovery records.
- Command controls demonstrate specified failures stopping in those
  commands. They do not establish reliable compliance by every future Guide.
  Native handling of genuinely damaged or missing recovery records before
  verification remains untested.
- Recorded final-candidate client observations are on a Mac: CLI
  help/status, live CLI help, and Desktop Code-tab **help only**. They do
  not establish all setup/recovery behavior or client support on other
  platforms. See the [verification record](GUIDED_CODING_VERIFICATION.md).
- The full PHENIX server suite was **not run for this GC-only publication**;
  the Developer explicitly waived it for that pilot. Local package tests
  are not an equivalent, and the waiver is not a standing exemption.
- There is no parallel-job coordinator, automatic shared-installation
  reservation or new installation-restoration engine. Client transcripts,
  settings bookkeeping and auto-memory can exist outside project records.
- The new checker rules and tools check form and Git state. They do not
  authorize a push, update an installation, reserve a resource or
  establish a reviewer's independence. The native observations of this
  revision (a history listing, the CLI title instruction, damaged and
  missing recovery records, fresh-session recovery) were made with one Mac
  CLI client in non-interactive auto mode; Desktop title behaviour was
  observed only as far as the verification record states.

The source carries the literal `r10 rev17 candidate` label used by its
startup check. That label is not evidence of publication or approval; use
the recorded commit, source manifest and decision instead.
The development tag `guided_coding-r10-rev16-pilot` was removed under the
repository's release-only tag convention. It is not an installation target
and should not be recreated. Pilot status belongs in these docs and the
publication record, not in a development tag.

## Start using it

Keep one reviewed source copy at
`cctbx_project/libtbx/guided_coding/`. Register it once per local Claude
Code configuration at `~/.claude/skills/guided_coding`, then invoke it in
the project you want to work on:

```text
/guided_coding setup
/guided_coding Fix the regression described in this issue.
```

Registration and project adoption are separate. Before guided work, the
project's current authority must adopt the exact contract identity and
name its own project method. Setup helps prepare that decision; it does
not assume it. Existing settings are proposed defaults, not new permission.

The package supplies general setup defaults and a method template. It
contains no personal PHENIX profile, account, server requirement or `t96`
definition. An optional project defaults card can be supplied separately;
another developer adapts its paths and permissions to their own environment.

Ordinary work remains ordinary in a fresh conversation when no always-on
project instruction activates GuidedCoding. Starting a new conversation
does not erase persistent instructions or auto-memory.

## Documentation map

| Read this | For |
| --- | --- |
| [User guide](GUIDED_CODING_USER_GUIDE.md) | Register, set up a project, choose commands, update and uninstall |
| [Architecture](GUIDED_CODING_ARCHITECTURE.md) | What lives where, who decides, what the tools actually enforce |
| [Verification](GUIDED_CODING_VERIFICATION.md) | Source-to-claim map, recorded observations, limits and check commands |
| [Skill entry](../SKILL.md) | The instructions loaded by an explicit invocation |
| [Contract](../payload/DEVELOPER_GUIDE_CONTRACT.md) and [roles](../payload/ROLES.md) | Authority, responsibilities and review boundaries |
| [Setup workflow](../payload/SETUP.md) and [general defaults](../payload/SETUP_DEFAULTS.md) | Recover settings and prepare a checked save |
| [Guide](../payload/GUIDE.md) and [Worker](../payload/WORKER.md) | The governing work procedure |
| [Reviewer brief](../payload/OUTSIDE_REVIEWER_BRIEF.md) and [transport](../payload/REVIEW_TRANSPORT.md) | Prepare and read a review bundle |

This pilot supplies a structured process and inspectable evidence. It
makes no general claim that an LLM following it cannot make mistakes.
