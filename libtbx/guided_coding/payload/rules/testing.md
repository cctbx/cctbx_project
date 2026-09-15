# Testing rules (core, for this module)

Who answers what: `.claude/DECISIONS.md`. Why the rules are what they are: `.claude/RATIONALE.md`. Full verification procedure: `.claude/skills/workflow/VERIFY.md` - the Worker loads it at the start of every verification and follows it exactly; it is not summarized from memory.

Core rules that bind everyone, subagents included:
- Subagents construct on the NAMED starting commit, verify that themselves as their first action, and use ABSOLUTE paths derived from their own copy root in every command.
- The executed decisive test is always the approved bytes, verified by hash, on both code versions.
- No check is recorded as a bare "passed": preserve command, identities, context, exit status, and raw output (empty output recorded as empty).
- The tree lock (VERIFY.md) is taken atomically before any stage that executes code from the live tree.
- On anything unexpected or unclassifiable, .claude/INCIDENT.md governs.
