# The incident procedure

Triggered by ANY of: an interrupted verification, a subagent writing outside its copy, a conflict mid-rebase, a manifest mismatch, a lost or ambiguous lock (local or server), evidence that two things used a tree at once, or anything the Worker cannot classify. When this document applies, it overrides everything except the developer's word.

Three distinct permissions, so the procedure never contradicts itself:
- (i) STOP: halt new source mutations and safely stop or suspend identified, Worker-owned work.
- (ii) PRESERVE: bounded writes are always permitted to the incident record, to tags, and to a designated backup location - never to source.
- (iii) RECOVER: only the specific actions the Developer subsequently authorizes.

The five steps, always in order:

1. FREEZE (permission i). No cleanup, no retry, no "quick fix".
2. INSPECT, read-only. Establish actual state: git status, fingerprints of plan-scope files AND of any path the incident may have touched outside scope, lock states, worktree list. Never reason from what should be true.
3. PRESERVE (permission ii). Copy volatile scratch to the backup location, pin unreachable commits with tags, record fingerprints.
4. REPORT. Four labeled parts to the Developer; a self-contained packet to the Guide (the incident analyst) for anything above "nothing lost". Severity vocabulary: "nothing lost" (only with evidence - never the default while inspection is incomplete), "IMPACT UNKNOWN", "evidence invalidated", "developer files at risk".
5. RECOVER (permission iii). Only the steps the Developer authorizes, smallest reversible first, verified by content afterward. The incident and its recovery go in the change record and the packet's incident log - never silently absorbed.
