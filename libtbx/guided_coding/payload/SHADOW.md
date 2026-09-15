# The shadow installation (remote verification without touching the host installation)

How a per-change workspace runs a project suite on a server whose runtime binds its own import paths at dispatch time. Established by the 2026-09-11 shakedown on this profile's server; generated per site from this template and PROVEN per site before use (see "Commissioning" below).

## Why a plain checkout is not enough

On PHENIX-family installations the dispatchers (`build/bin/*`) rebuild the import environment themselves: an inherited PYTHONPATH is erased unless PHENIX_TRUST_OTHER_ENV is set, and when kept it is appended LAST, after the installation's own module paths. So no environment variable can redirect `import <module>` to a workspace checkout. A workable workspace must supply ITS OWN dispatchers - its own LIBTBX_BUILD.

## The shadow root

Under the per-change workspace `<WORKSPACE>/`, build a root `<WORKSPACE>/root/` that looks like an installation but swaps the modules under test:

- `root/modules/<name>` -> symlink to the host installation's module, for EVERY module EXCEPT the swapped set.
- `root/modules/<swapped>` -> symlink to that repository's workspace checkout, one per repository the change touches.
- `root/conda_base` -> symlink to the host's conda_base (dispatchers reach it as `$LIBTBX_BUILD/../conda_base`).
- `root/build/<entry>` -> symlink to the host's build entries, EXCEPT `bin`.
- `root/build/bin` -> a REAL COPY of the host's dispatchers. This is the crux: a dispatcher resolves its own physical location to compute LIBTBX_BUILD, so copied dispatchers make the shadow root the build root, while a symlinked one would resolve back to the host.

**What this arrangement is and is not.** It is an EXECUTION arrangement: it controls which code the runtime imports and which tests it harvests. It is NOT write protection. Symbolic links are transparent - a test that writes through `root/modules/<unswapped>/...` writes into the host installation. Host safety is provided separately (next section), never by the shadow itself.

## Protecting the developer's working installation (a separate mechanism, chosen in the profile)

The shadow links into SOME installation; a test writing through those links writes into it. So the installation the shadow links into must be one of, recorded in the profile:
- **An administrator-supported read-only execution environment**: the test process runs under a restricted account with read-only access to the shared installation and write access only to designated output areas. The only option that also resists a misbehaving test.
- **A separately commissioned DISPOSABLE test installation** (the practical choice for a developer with a build routine): a second configured installation on the server, built once by the normal update/build procedure, used only by this procedure, and rebuilt when stale. The shadow links into IT - never into the working installation - so a test's stray writes land in something rebuildable. Its currency is checked against the intended version before each comparison, like any environment.
- **Normal-account execution against the working installation, with the exposure stated accurately**: any write through a link changes the developer's working installation; a post-run check (`git status --porcelain` in every host module repository plus hashes of tracked files) DETECTS changes to tracked files only - it cannot see changes inside untracked or ignored files, nor in `conda_base` or the build - and nothing is prevented. Recorded as the weakest posture; permitted only by explicit election.

No routine recursive permission changes are made to any installation: they cannot be reversed faithfully without preserving every object's original mode, they interfere with other work, and an interruption leaves the installation in an altered state.

Tests are launched from a fresh scratch working directory inside `<WORKSPACE>` for EACH execution, never from a host directory, so ordinary cwd-relative output lands in the workspace and cannot carry from one run to the next.

## The swapped set is the plan's repository scope

Whatever repositories the approved plan touches are the swapped ones - a phenix-module change swaps phenix; a cctbx_project change swaps cctbx_project; a test-only change swaps the test repository; a cross-repository change swaps both. Each swapped repository is staged by its own git bundle at its own named commit. Everything else symlinks to the frozen host installation.

## Proof before use, every time

Two properties must be shown by execution, per swapped module, not assumed:

1. **Import proof** - code executed through the shadow dispatcher reports the swapped module's `__file__` inside the shadow root's workspace checkout.
2. **Harvest proof** - the suite's module directory (for PHENIX: `PHENIX_DIST`) resolves to the shadow's swapped module, AND a uniquely named dummy `tst_` file placed only in the workspace checkout appears in the harvested roster (then is removed). Without this, a mixed mode - harvest from the host, import from the workspace - can pass a fail/pass check while running the host's tests.

A decisive fail/pass pair whose test exists in BOTH the host and the workspace cannot by itself distinguish these modes. When an already-recorded pair is reused as the instrument, run the harvest proof alongside it.

If no arrangement passes both proofs, remote verification is reported UNAVAILABLE and the decision returns to the developer. The developer's standing installation is never a silent fallback; using it as the runtime requires an explicit recorded election with the risks stated.

## Limits, stated

- **Compiled code is out of scope.** Built libraries come from the frozen host build via symlink, so a change to compiled sources cannot be shadow-tested without a rebuild - the same boundary as this procedure's pure-source scope.
- **One canonical path.** The workspace, the shadow root, and the server lock live under one canonical physical path, recorded in the profile; a lock guarding a different path than the work is not exclusivity.
- **Host activity during a run.** While a shadow run executes, the host installation must not be edited or rebuilt EXCEPT in the swapped repositories - everything else is read live through symlinks, and a rebuild between a comparison's two runs voids the comparison. The run's announcement names which repositories are swapped.
- **Site-specific by construction.** This template exploits how libtbx dispatchers compute their build root. A non-PHENIX profile needs its own arrangement established at its own commissioning.

## Commissioning (per site, once)

Generated from the profile's two answers - the host installation path and the workspace location - then proven: build the shadow, run both proofs against a known fail/pass pair, run the interruption drill, and record the result. Only the executed proof clears the arrangement for use; generating it from this template is setup, not evidence. The commissioning record is kept even if the developer later turns the remote leg off, so turning it back on is re-verification against the current build rather than starting over.
