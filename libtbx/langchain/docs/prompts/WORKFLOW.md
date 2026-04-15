# LLM WORKFLOW RULES (CRITICAL)

This codebase is worked on incrementally with possible interruption at any time.

## PRIMARY PRINCIPLE

The system must ALWAYS be in a recoverable state.

At any moment, a new session must be able to resume work using ONLY:

* HANDOFF.json
* CODE_INDEX.md
* This file

---

## HANDOFF.json REQUIREMENTS

HANDOFF.json is the single source of truth.

It MUST always include:

* current_task
* status
* relevant_files
* completed_steps
* remaining_steps
* decisions
* assumptions
* open_questions
* last_stable_state
* next_step

It must be:

* accurate
* up to date
* sufficient for full restart

---

## MANDATORY CHECKPOINT RULES

### Micro-Step Checkpointing (NEW)

After **any reasoning, analysis, investigation, or code change**, immediately update HANDOFF.json with:

* current progress
* partial results or findings
* next planned step

Do not wait until the end of a step or task.
Even intermediate thoughts must be recorded.

Before proceeding to the next micro-step, confirm that HANDOFF.json accurately reflects the current truth.

### Standard Checkpoints

You MUST update HANDOFF.json:

1. After understanding the task
2. After identifying relevant files
3. After forming a plan
4. After EACH code change
5. BEFORE running tests
6. AFTER test results
7. BEFORE ending response

Never defer updating HANDOFF.json.

---

## STEP EXECUTION MODEL

Work must proceed in micro-steps:

1. Understand problem
2. Identify files
3. Create PLAN
4. Refine PLAN (optional but recommended)
5. Implement change
6. Validate (tests)
7. Finalize

After EACH micro-step:

* **Checkpoint immediately** in HANDOFF.json

---

## PLAN FORMAT (MANDATORY)

All plans MUST follow this exact format:

PLAN:

1. <step description>
   - goal: <what this step achieves>
   - output: <what is produced>

2. <step description>
   - goal: ...
   - output: ...

Rules:

* Maximum 6–8 steps
* Each step must be independently executable
* Each step should take no more than one short iteration to complete
* Each step must produce a clear, concrete output
* Do NOT use JSON
* Do NOT write long prose

The plan should be concise and easy to follow.

---

## ITERATIVE PLAN REFINEMENT (RECOMMENDED)

Plans should be improved through brief critique cycles when the task is non-trivial.

Process:

1. Generate initial PLAN
2. Pause for external review (e.g., Gemini or human)
3. Incorporate feedback
4. Repeat once more if needed (maximum 2 refinement cycles)

Rules:

* Limit to 1–2 refinement iterations
* Do NOT proceed to implementation until refinement is complete
* Summarize key changes after each refinement
* **Checkpoint HANDOFF.json after each revision**

If no external feedback is provided:

* Perform a self-critique:

  * identify weaknesses
  * improve the plan once

---

## PLAN BEFORE EXECUTION

You MUST create a PLAN before making any code changes.

Do NOT begin implementation until the plan is written, refined, and critically reviewed.
**Checkpoint HANDOFF.json immediately after PLAN creation and after any refinement.**

---

## TEST EXECUTION RULES

* Do NOT run full test suite without checkpointing
* Prefer:

  * targeted tests first
  * then full suite
* ALWAYS checkpoint HANDOFF.json BEFORE running tests

---

## INTERRUPTION PROTOCOL

If you suspect you may run out of time:

1. STOP immediately
2. Update HANDOFF.json with:

   * current progress
   * partial work
   * exact next step
3. Do NOT start new work

After any interruption, **restart from HANDOFF.json**, not prior memory.

---

## FILE CHANGE RULES

* Modify no more than 3 files per step unless necessary
* Prefer minimal, surgical diffs
* Do not refactor unrelated code

---

## OUTPUT REQUIREMENTS (UPDATED FOR FULL FILES)

When making changes, always provide:

1. **Full copies of all changed files**, preserving the same directory structure as the original code archive provided to you.
2. A compressed archive (tar.gz) of these files matching the original archive structure.
3. Updated HANDOFF.json (mandatory after every micro-step).
4. Short explanation of changes and reasoning for each modified file.

Notes:

* Do not send only diffs or partial files.
* File names, paths, and directory hierarchy must match the original codebase layout.
* The archive should be self-contained so it can be extracted and used directly without manual reassembly.

After each code-change micro-step:

* Stage the changed file(s) in the archive structure
* Checkpoint HANDOFF.json
* Only then continue to the next micro-step

Prompt snippet for Claude at output:

```text id="output_prompt"
Provide full copies of all changed files in a tar.gz archive that preserves the original directory structure.
Include updated HANDOFF.json and a short explanation for each file.
Do not provide patches or partial files alone.
```

---

## FAILURE HANDLING

If tests fail:

* Record failure in HANDOFF.json
* Include error output summary
* Propose next fix
* Do NOT proceed without updating state

---

## SUCCESS CRITERIA

A task is ONLY complete when:

* Code changes are implemented
* Tests pass
* HANDOFF.json reflects final state
* Full file copies have been provided in proper archive structure

---

## REMEMBER

You are not just solving problems.

You are maintaining a persistent, recoverable work state across sessions.

