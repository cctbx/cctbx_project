# GuidedCoding verification and limits

This file separates three kinds of evidence: what the supplied source
implements, what earlier execution records report, and what remains
unestablished. A successful checksum, test suite or client-loading check
is not approval to integrate or publish.

## Status and identities

The published implementation baseline is r10 rev16
`enumcheck-20261002T194333Z`, published on 2026-10-02 as the limited opt-in
pilot. Its recorded identities are:

| Item | Identity |
| --- | --- |
| Published commit | `c36887c7f489018af4f91773246ecde32d1b4e24` |
| Source archive SHA-256 | `fce3627e46982f3b5fad1baaad7ada538984e9087a84b2147548b4564e844968` |
| Source manifest SHA-256 | `b86d0490df59f854c214e6771ad8bb44f03069a38612be198edff4ebb11238e0` |
| Contract version / SHA-256 | `2026-09-17` / `ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b` |
| Publication review bundle SHA-256 | `0e4c112d3726026bd397a4e50ac2a3cff79b146883f8ceab269e3a069bc2c6ca` |
| Publication packet identity | `9834e7a5bfdf6b8dd8401972265715f527432ba277b63c4ee5c50f130508f44b` |

The `docs-20261003` revision changes documentation and its source identity.
It retains the contract, skill entry, Python tools, tests and SETUP §5
shell recipe. Earlier reviews and test runs apply to their recorded source;
they do not approve this new documentation diff. Its new manifest must be
bound to its own check and publication records. The literal release-family
check still expects `r10 rev16 candidate`; it does not establish an exact
revision or approval.

The October 2–3 publication, client and installation facts below come from
the supplied controller and outside-review records. The documentation pass
checked the uploaded baseline against its complete manifest and read the
implementation and tests; it did not rerun those sessions, inspect live
remote refs or act as an Outside Reviewer. Raw historical evidence is kept
in the associated task/review records, not duplicated in the generic package.

## Release and source checks

First compare an archive with the independently supplied SHA-256 in its
review/publication handoff. From a trusted extraction, replace the path
below with the canonical source root and run this complete Bash block:

```bash
cd /absolute/path/to/libtbx/guided_coding &&
shasum -a 256 -c SOURCE_MANIFEST.sha256 &&
GC_PAYLOAD_ROOT="$(pwd -P)/payload" python3 -I -B payload/tools/screen_check.py verify-source . &&
cat payload/RELEASE &&
grep -Fq 'r10 rev16 candidate' payload/RELEASE
```

Every dependency uses `&&`. A failed directory change or listed checksum
must stop before executing the checker; a failed inventory check must stop
the later commands. Do not append unguarded preparation or registration.
After a successful block, separately compare the manifest's SHA-256 and
full release label with the expected handoff identity. Repeat the guarded
source check before central tools and consequential work as the skill
requires.

`verify-source` rejects missing/changed listed bytes, extra files, symlinks
and multiply linked files. The accepted installer-bytecode pattern is an
exception, not authentication of bytecode contents. `-I` isolates imports
from the current directory and `PYTHONPATH`; `-B` suppresses writing caches.
The tool invocation uses source rather than package caches. These controls
establish consistency, not supplier identity or model compliance.

## Trace a claim to the source

For example: **registration refuses an existing destination instead of
creating a nested link.**

1. Open [`register_skill` in `screen_check.py`](../payload/tools/screen_check.py).
   It computes the single `skills/guided_coding` destination, checks
   `exists()` or `is_symlink()`, and calls `fail` if occupied. The creation
   is a direct `os.symlink`, not a recursive `ln -s` command.
2. Read `SkillRegistrationChecks` in
   [`tst_screen_check_v20.py`](../tests/tst_screen_check_v20.py). The fresh/
   second-run case and the other-link/directory/dangling-link case assert
   refusal and absence of a nested link. Other cases cover version and
   parent-path refusals.
3. To establish that those tests ran for a particular revision, inspect
   its captured command, source identity, stdout/stderr and exit status.
   Reading the test source alone establishes its intended assertion,
   not an observed pass.

Other important mappings:

| Claim | Implementation / test source | What it does not establish |
| --- | --- | --- |
| Complete-source consistency | `verify_source`, `regular_files`; `SourceInventoryChecks` | Origin, semantic correctness or authenticated cache contents |
| CLI minimum 2.1.281 | `MIN_CLAUDE_VERSION`, `check_claude_version`; `ClaudeVersionChecks` | Desktop runtime, actual instruction loading or trustworthy PATH executable |
| Screen form and result identity | `blocks`, `developer_view`, `code_identity`, `check`; `ScreenChecks` | Correct criterion, true logs, complete tests or Developer authorization |
| Matching outside reading | `reading_matches`; conditional/waiver/pending controls | Independence or authenticity of the reviewer, satisfaction of a semantic condition |
| Two-domain review transport | `make_bundle` in `review_bundle.py`; `ReviewBundleChecks` | Privacy, evidence completeness or permission to publish |
| Checked setup ordering | `SKILL.md`, `SETUP.md` §5 and its embedded recipe | Mechanical enforcement by Claude Code or reliable obedience by every Guide |
| Test isolation in the repository | [`libtbx/tst_guided_coding.py`](../../tst_guided_coding.py) | Complete source verification of arbitrary unlisted files in the original installation |

There are **45 test methods** in the two shipped package test files.
Subtests exercise additional cases; a count of 45 is not 45 independent
behavioral claims. The setup shell recipe and native LLM scenarios are
not covered by those 45 unit methods. Their controls were supplied in
separate review evidence.

## Recorded checks of the published implementation

| Evidence | Recorded result and reach |
| --- | --- |
| Enumeration correction and follow-up reading | The two `find` producers' failures, including partial listings, stop verification. A/B/B2 did not print VERIFY_OK; the tested driver held apply; protected files were unchanged. Seven earlier recipe controls retained their results. R-F1 was closed. |
| Package and repository checks | Controller reported 45 package tests passing on the Mac. The publication reviewer ran the repository harness on Linux: unit success with one case-insensitive-filesystem skip and a libtbx precompile skip. Skips are not passes of those branches. |
| C1, disposable CLI configuration | Exact enumcheck candidate, CLI 2.1.284: help and status, source guard and command-to-expansion comparison. Reach is CLI loading/help/status, not setup/recovery. |
| 4b, live CLI | Help loaded the expected source and ran the guard. The successful attempt followed an expired-auth failure and user sign-in. |
| D1, Desktop Code tab | Help loaded the personal central link and ran the source guard; transcript reported runtime 2.1.286. Reach is **loading/help only**, not Desktop status or full setup/recovery. |
| Publication r2 | Outside reading: PROCEED for the limited pilot; P-F1/P-F2/P-F3/P-N1 closed, R-F1 retained closed. A clerical description of three test assertions was corrected externally. Genuine publication screen passed; Developer PUBLISH and verified remote pushes were recorded. |
| Mac and shared Linux installation synchronization | Controller reported clean checkouts at the publication commit on the Mac and on one shared Linux installation, source equality between them and repository harness success on both. Synchronized GC files on that installation are not native Claude Code client validation: no Claude client was installed or run there. This was not a PHENIX full-suite run. |
| Tag removal, October 3 | Controller reported the pilot tag absent remotely, its local copy removed after preservation, and the published commit still in upstream history. Repository tags are reserved for releases; use the commit/manifest, not the removed tag. |

The full PHENIX server suite and a baseline/candidate roster were **NOT RUN /
NOT PRODUCED for the GC-only pilot publication**, under the Developer's
specific waiver. The publication record states this publicly. Local tests
do not replace that suite, and the waiver does not apply automatically to
a later change.

## Remaining qualifications

- **Native damaged-record branch:** no recorded native session of the final
  enumcheck candidate demonstrates the Guide responding to genuinely
  damaged or missing recovery records before verification. Command controls
  cover stopping in the commands; the driver/Guide still determines what
  gets invoked next.
- **Exact proposal and recovery history:** earlier cases omitted complete
  proposals or approved recovery text. Assisted and later native cases
  showed improvement; they do not erase those failures. A relocated setup
  record was disclosed but not wholly literal to its displayed proposal.
- **Fresh-session context:** disposable auto-memory was written outside the
  target during earlier trials. Later sessions using that configuration
  cannot be called memory-free or dependent only on target records without
  checking those inputs. No broad fresh-recovery guarantee follows.
- **A1 provenance:** an untouched pre-redaction debug log was not located
  by the reported bounded search. Other transcripts and partial excerpts
  remain; the reviewer found no material behavioral claim unsupported solely
  by that loss. The complete original and “only alteration” claim remain
  unverified. This is a provenance limit, not proof of another behavioral
  failure.
- **Client/platform scope:** a passing CLI version check is not a Desktop
  version check. Earlier TUI observations had continuation/permission-mode
  qualifications. Native Windows, real csh and all-client setup/recovery
  behavior are not established by the final-candidate observations above.
  The repository wrapper explicitly limits the tools to macOS and Linux;
  its Windows skip is not evidence of Windows support.
- **Isolation:** observed file equality and absence of recorded tool calls
  have bounded reach. They do not establish whole-machine immutability,
  all syscalls or unchanged live authentication state.
- **General effectiveness:** this package does not establish a failure
  rate, superiority to unstructured LLM use, or reliability for every future
  Guide. Separate test authorship is recorded, not proof of independence.

## How to check a later documentation or source revision

Keep the previous evidence unchanged and bind the new source separately.
For a docs-only revision, check the diff against the implementation,
command syntax and preserved authority boundaries. Record unchanged
executable/contract bytes explicitly. If instructions or claims change,
assess those changes; a unit suite cannot establish their behavioral effect.

After the source guard succeeds, package tests can be run from a trusted
clean, manifest-listed copy with a resolved, symlink-free temporary path:

```bash
python3 -I -B -m unittest discover -s tests -p 'tst_*.py' -v
```

From an appropriate cctbx environment, the repository wrapper is:

```bash
libtbx.python /absolute/path/to/cctbx_project/libtbx/tst_guided_coding.py
```

These are separate authorized check steps, not commands to continue after
a failed source guard. The wrapper copies listed files into scratch, uses
a resolved temporary directory, runs the package tests, and checks that
listed installed bytes did not change. When libtbx precompilation is
available it exercises that branch on another copy. It skips Windows and
Python 2; it can also skip the real precompile branch. Report each skip.
The wrapper reports existing unlisted files rather than using them for the
test copy, so retain the complete-source check as a separate requirement.

For native validation, identify the exact source, client/runtime, selected
skill path, loaded project instructions and persistent context. Inspect
recorded command-to-expansion linkage, actual operations and resulting
files. A version string, registration link, model-reported path or procedural
simulation alone is insufficient. Matching decoded transcript text is not
proof of file-byte loading. Select observations for the changed behavior
and intended support scope; do not label an earlier candidate's session as
a fresh run of the new one.

The Outside Reviewer evaluates the new source and bound evidence; the
Developer's acceptance, integration, activation and publication decisions
remain distinct. A genuine final screen must cite a qualifying reading of
its own frozen packet. Use the existing contract/project rules for any
waiver; do not infer one from this document.

## Relevant history

Rev10–11 corrected exclusive registration and configuration-path handling.
Rev12 corrected user guidance from the Mac trial. Rev13 added the narrow
installer-bytecode exception and source-only tool loading. Rev14 introduced
recoverable, relevant setup. Failed rev15 parallel-job work was excluded
from the rev16 lineage.

Rev16 stabilized rev14, then corrected fail-fast checksum instructions,
approved-text preservation, explicit comparison failure propagation and
checked file enumeration. A failed `check && echo OK` inside a script was
not a sufficient stopping rule; unchecked producer exits were also not
sufficient. The final published implementation is enumcheck, not an earlier
F1/checked-save/failprop candidate. Earlier reviews belong to those earlier
identities; frozen packets and errata retain that history.

Current client reference material: [skills](https://code.claude.com/docs/en/skills)
and [instruction loading](https://code.claude.com/docs/en/memory).
Vendor documentation describes the client; it is not a GC native test run.
