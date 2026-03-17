# How to Program with an LLM

A practical guide for programmers who are new to using
large language models as coding partners. The first
half covers general principles that apply to any LLM
and any codebase. The second half describes the
specific prompt files and workflow we use on this
project.

---

## Part 1: General Principles

### 1.1 Mental Model

An LLM is not a compiler, not a search engine, and not
a junior developer. The closest analogy is a well-read
colleague who has seen a vast amount of code but has
never run any of it. This has practical consequences:

- It can produce fluent, plausible code very quickly.
- It cannot execute code in its head. It does not
  truly know whether its output works.
- It will sometimes write code that looks correct,
  passes a casual reading, and fails at runtime on an
  edge case it never considered.
- It has no persistent memory between sessions unless
  you explicitly provide context.

Your job is to be the **verifier and navigator**. The
LLM generates; you validate. This division of labor
is the foundation of everything else in this document.

### 1.2 Context is Everything

LLMs have a fixed context window — the total amount of
text they can see at once (the conversation so far,
any attached files, and their own response). Everything
the LLM needs to know must be inside that window or it
effectively does not exist.

**What to provide:**

- The specific files the LLM will read or modify.
  Don't give it the entire codebase if you only need
  changes to two files.
- Coding standards and style guides. The LLM will
  match whatever conventions you show it.
- Error messages, test output, and log snippets when
  debugging. Copy the actual text — don't paraphrase.
- A description of what you want, including constraints
  (performance requirements, compatibility, style).

**What to leave out:**

- Files that are not relevant to the current task.
  Extra context dilutes the LLM's attention.
- Redundant repetitions of the same information. State
  it once, clearly.

**Practical tips:**

- Attach files rather than pasting them into the chat
  when possible. Pasted text makes the conversation
  harder to follow.
- If the conversation is getting long and the LLM
  seems confused, start a fresh session with a clean
  summary rather than continuing.
- When providing code, include enough surrounding
  context that the LLM can see the function signatures,
  imports, and calling patterns — not just the single
  line you want changed.

### 1.3 How to Ask for Code

The quality of what you get is directly proportional to
the clarity of what you ask for. Vague requests produce
vague code.

**Good request:**

> Add a `_safe_float()` helper to `metric_evaluator.py`
> that converts string or None values to float, returning
> None on failure. Then use it at lines 352, 392, and in
> `calculate_improvement_rate()` to coerce `r_free`
> values before arithmetic. Follow the same pattern as
> `_safe_float()` in `metrics_analyzer.py`.

**Weak request:**

> Fix the type error in metric_evaluator.py.

The good request tells the LLM what to build, where to
put it, what pattern to follow, and which call sites to
update. The weak request forces the LLM to guess all of
that.

**Structure for effective requests:**

1. **What** — the change you want
2. **Where** — which files and locations
3. **Why** — the problem it solves (helps the LLM make
   reasonable choices in ambiguous situations)
4. **Constraints** — style, compatibility, patterns to
   follow or avoid
5. **Examples** — show existing code that does
   something similar, if applicable

### 1.4 The Plan–Review–Implement Cycle

Never let an LLM jump straight to writing code for a
non-trivial task. The cost of a bad plan executed
thoroughly is much higher than the cost of planning
time.

The recommended cycle:

1. **Plan.** Ask the LLM to produce a written plan —
   problem description, approach, implementation steps,
   risks. Review it yourself.
2. **Cross-review.** For important changes, give the
   plan to a second LLM (or a colleague) for critique.
   A different model catches different blind spots.
3. **Revise.** Feed the critique back. Limit to 1–2
   revision rounds — diminishing returns set in fast.
4. **Implement.** Only now does the LLM write code,
   following the agreed plan step by step.
5. **Verify.** You run the code, run the tests, and
   confirm correctness. The LLM cannot do this for you.

This cycle prevents the most expensive failure mode:
the LLM confidently building out an approach that was
wrong from the start.

### 1.5 Verification Discipline

This is the single most important habit to develop.

**Never trust LLM output without verification.**

This applies even when:

- The code looks obviously correct
- The LLM says "this is tested and working"
- The change is small and apparently trivial
- You asked for a trivial refactor

Concrete verification steps:

- **Parse-check** every changed file. At minimum:
  `python -c "import ast; ast.parse(open('file.py').read())"`
- **Run the tests** that exercise the changed code.
  Not "the LLM says the tests pass" — you run them.
- **Read the diff** carefully. LLMs sometimes make
  unrelated changes, remove code they thought was
  unnecessary, or subtly change logic while claiming
  to only reformat.
- **Check edge cases** the LLM may not have considered:
  None values, empty collections, missing keys,
  unexpected types, concurrent access.

### 1.6 Common LLM Failure Modes

Understanding how LLMs fail helps you catch problems
before they reach your codebase.

**Hallucinated APIs.** The LLM invents function names,
keyword arguments, or library features that don't
exist. It does this with complete confidence. Always
verify that the API calls it writes actually match the
real signatures.

**Confident wrongness.** When an LLM makes an error
and you point it out, it will sometimes apologize and
produce a new answer that is equally wrong but
differently worded. If a correction round doesn't
converge, re-state the problem from scratch rather
than asking the LLM to "try again."

**Sycophantic agreement.** If you suggest a wrong
approach, many LLMs will agree with you rather than
push back. Don't use leading questions like "shouldn't
we use X here?" if you're genuinely unsure. Instead
ask "what are the options for this?" and evaluate the
answer.

**Drift over long conversations.** As the conversation
gets longer, the LLM gradually loses track of earlier
instructions and constraints. If you notice the LLM
forgetting your coding standards or making mistakes it
avoided earlier, it's time to start a fresh session
with a clean summary.

**Unasked-for changes.** LLMs sometimes "improve"
code you didn't ask them to touch — renaming
variables, refactoring adjacent functions, removing
comments they consider redundant. Always diff the
output against the original and reject unrelated
changes.

**Incomplete changes.** The opposite problem: the LLM
changes the function but forgets to update the callers,
the tests, or the documentation. Ask explicitly:
"what else needs to change to keep everything
consistent?"

### 1.7 Session Management

LLM sessions are inherently fragile. They can be
interrupted by context length limits, timeouts,
network issues, or simply because you closed the
browser tab. Plan for this.

**Checkpointing.** After each meaningful step, have
the LLM write a checkpoint file (we use
`HANDOFF.json`) that captures: what was done, what
remains, what decisions were made, and what the next
step is. If the session dies, a new session can pick
up from the checkpoint with no loss.

**Session length.** Shorter sessions with clear
handoffs are more reliable than marathon sessions.
After about 15–20 back-and-forth exchanges, consider
wrapping up and starting fresh. Each new session gets
a clean context window and a fresh start on the
instructions you provide.

**Starting a new session.** Always provide:

1. The reference documents (coding guidelines, workflow
   rules)
2. The current checkpoint file
3. The relevant source files
4. A clear statement of what to do next

Do not rely on the LLM "remembering" anything from a
prior session. Even if the LLM platform has a memory
feature, treat each coding session as standalone.

### 1.8 Using Multiple LLMs

Different LLMs have different strengths and weaknesses.
Using more than one is genuinely useful in two
situations:

**Plan review.** Have one LLM write a plan and a
different LLM critique it. They make different kinds
of mistakes, so the reviewer will often catch things
the author missed. This is the most effective multi-LLM
pattern.

**Stuck debugging.** If one LLM has gone in circles on
a bug, describe the problem fresh to a different model.
A different set of biases often finds a different (and
sometimes correct) path.

There is no need to use multiple LLMs for routine
coding tasks. The overhead of managing multiple
sessions isn't worth it for straightforward changes.

### 1.9 What LLMs Are Good and Bad At

**Good at:**

- Boilerplate and repetitive code (serialization,
  CRUD, test scaffolding)
- Translating a clear specification into code
- Explaining unfamiliar code or error messages
- Suggesting approaches you hadn't considered
- Drafting documentation, docstrings, and comments
- Refactoring with a well-defined pattern (e.g.,
  "add null-safety to all `.get()` calls in this file")

**Bad at:**

- Reasoning about complex state across many files
- Performance optimization (it optimizes for
  readability, not speed, unless told otherwise)
- Knowing the actual runtime behavior of code it
  writes
- Maintaining consistency over long sessions
- Understanding your specific codebase's unwritten
  conventions (unless you tell it)
- Anything requiring information that isn't in its
  context window

### 1.10 Code Review for LLM-Generated Code

Review LLM-generated code the same way you would
review a pull request from a new team member who is
talented but unfamiliar with your project:

- Does it follow the project's style and conventions?
- Does it handle the edge cases that matter for this
  codebase?
- Did it introduce any unnecessary dependencies?
- Did it change anything beyond what was requested?
- Are the names clear and consistent with existing
  code?
- Would you be comfortable maintaining this code six
  months from now?

If the answer to any of these is "no," send it back
for revision — just as you would with a human author.

---

## Part 2: Project-Specific Prompts and Files

The `docs/prompts/` directory contains a set of
ready-made prompts and reference documents for driving
LLM coding sessions on this project. Each file has a
specific role in the workflow.

### 2.1 Reference Documents (attach at session start)

These files give the LLM the context it needs to
produce correct code. Attach them at the beginning of
any coding session.

| File | Purpose | When to provide |
|------|---------|-----------------|
| `CCTBX_LLM_PROGRAMMING_GUIDELINES.md` | Coding standards, style rules, cctbx/PHENIX-specific patterns, and common pitfalls. This is the base guideline that all LLM sessions need. | **Every session.** Attach before any coding request. |
| `WORKFLOW.md` | Defines the incremental work model: HANDOFF.json requirements, checkpoint rules, plan format, interruption protocol, and output requirements. | **Every session.** Attach alongside the coding guidelines so the LLM knows how to structure its work and maintain recoverable state. |

For agent-specific work, also attach
`AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md`, which covers
patterns and pitfalls specific to the `agent/`,
`knowledge/`, and `programs/ai_agent.py` code.

### 2.2 Session Prompts (use at specific moments)

These are short prompts you paste or send at specific
points during a session. They are not reference
documents — they are instructions that tell the LLM
what to do next.

**`WORKFLOW_PROMPT.txt` — Start a new session**

Use this as the opening message when starting a coding
session. It tells the LLM to read HANDOFF.json (if
one exists), summarize the current state, propose a
plan, and follow the checkpoint discipline defined in
WORKFLOW.md. Pair it with the two reference documents
above.

Typical start of session:

1. Attach: `CCTBX_LLM_PROGRAMMING_GUIDELINES.md`,
   `WORKFLOW.md`, code archive, and `HANDOFF.json`
   (if resuming)
2. Paste the prompt below as the first message

```text
You are working on a large codebase with strict workflow rules.

Your primary goal is NOT just to complete tasks, but to maintain a continuously
valid HANDOFF.json so work can resume at any moment if interrupted.

## Core Rules

1. HANDOFF.json must ALWAYS be complete and restartable
2. After EVERY meaningful step, update HANDOFF.json
3. NEVER batch large work without checkpointing
4. Before any long operation (like running tests), update HANDOFF.json
5. If you suspect time/resource limits, STOP and write HANDOFF.json immediately

## Required Workflow

At the start:

* Read HANDOFF.json (if present)
* Summarize current task and state
* Propose a step-by-step plan
* Ask for any missing files BEFORE coding

During work:

* Break work into small steps
* After each step:

  * Update HANDOFF.json
  * Ensure it reflects current truth

Before running tests:

* Fully update HANDOFF.json with:

  * changes made
  * expected outcomes
  * what success/failure means

If interrupted, HANDOFF.json must allow a new session to continue with NO
additional context.

## Constraints

* Modify at most 3 files unless justified
* Prefer minimal diffs over rewrites
* Preserve public APIs unless explicitly required

## Output Format

When you make changes:

1. Show patch/diff
2. Update HANDOFF.json
3. Brief explanation

Acknowledge these rules and begin by reading HANDOFF.json or requesting it.
```

**`PLAN_PROMPT.txt` — Request a plan**

Use this when you want the LLM to produce a detailed
implementation plan before writing code. Provide it
along with a description of the problem or the bug
report. The output should be a markdown plan document.

```text
Please make a plan for fixing these problems. Include a full discussion of the
problem, the overall approach, the details of the approach, implementation plan,
and risks involved and their mitigation. write as md
```

**`REVIEW.txt` — Cross-review a plan with a second LLM**

After the first LLM produces a plan, give the plan
and this prompt to a **different** LLM (the prompt
text names Gemini, but any second LLM works). The
reviewer is instructed to critique each step: find
flaws, unclear assumptions, and workflow violations.
It will not give generic praise — it focuses on
concrete issues and suggested improvements.

```text
Gemini, you are a senior engineer tasked with **critically reviewing this plan**.
Your goal is to find every potential flaw, unclear assumption, or workflow
violation. Be specific and constructive; do not give generic praise.

For each step:
1. Identify issues or risks
2. Explain why it is a problem
3. Suggest improvements or alternatives

Limit to 1–2 paragraphs per step. Do not provide general compliments.
```

The review cycle:

1. LLM A produces a plan (via `PLAN_PROMPT.txt`)
2. Give the plan + `REVIEW.txt` to LLM B
3. Give LLM B's critique back to LLM A to revise
4. Optionally repeat once (max 2 refinement rounds,
   per WORKFLOW.md)

**`CONTINUE_PROMPT.txt` — Resume after interruption**

Use this when an LLM session was interrupted (timeout,
context limit, crash) and you are starting a new
session to continue the same work. Attach the current
HANDOFF.json and the code archive alongside this
prompt.

```text
You were interrupted.

Do NOT continue from memory.

1. Read HANDOFF.json
2. Summarize current state
3. Confirm the next step
4. Continue execution from there

If anything is unclear or missing, ask before proceeding.
```

### 2.3 Typical Workflow Sequence

A full task from start to finish typically uses the
files in this order:

```
1. START SESSION
   Attach: CCTBX_LLM_PROGRAMMING_GUIDELINES.md
           WORKFLOW.md
           code archive
           HANDOFF.json (if resuming)
   Send:   WORKFLOW_PROMPT.txt

2. PLAN
   Send:   PLAN_PROMPT.txt + problem description
   LLM produces plan as markdown

3. REVIEW (recommended for non-trivial changes)
   Give plan + REVIEW.txt to a second LLM
   Feed critique back to first LLM
   Revise plan (1–2 rounds max)

4. IMPLEMENT
   LLM executes the plan step by step,
   checkpointing HANDOFF.json after each step

5. IF INTERRUPTED
   Start new session
   Attach: same reference docs + HANDOFF.json
   Send:   CONTINUE_PROMPT.txt
```

### 2.4 Notes

- The reference documents
  (`CCTBX_LLM_PROGRAMMING_GUIDELINES.md`,
  `WORKFLOW.md`) should be attached as files, not
  pasted inline, to avoid consuming the LLM's
  context window with boilerplate.
- The session prompts (`WORKFLOW_PROMPT.txt`,
  `PLAN_PROMPT.txt`, `REVIEW.txt`,
  `CONTINUE_PROMPT.txt`) are short enough to paste
  directly into the chat.
- `REVIEW.txt` references Gemini by name but works
  with any LLM. The important thing is to use a
  different model (or at minimum a fresh session)
  so the reviewer is not biased by having produced
  the plan.
- `HANDOFF.json` is the critical link between
  sessions. The workflow rules in `WORKFLOW.md`
  exist primarily to ensure this file stays current.
  If the LLM is not updating it after each step,
  remind it.
