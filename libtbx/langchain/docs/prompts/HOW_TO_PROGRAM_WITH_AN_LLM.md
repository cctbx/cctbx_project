# How to Program with an LLM

A practical guide for programmers who are new to using
large language models as coding partners. The quick-start
below gets you going immediately. Part 1 covers general
principles that apply to any LLM and any codebase, with
the project-specific prompts shown in context where each
principle is introduced. Part 2 covers the project
workflow: rules of engagement, which guideline files to
attach and when, and checklists for the human side of
the process.

All prompt files live in `docs/prompts/`. The small ones
are reproduced in full throughout this document so you
can read the guide without opening anything else.

---

## Quick-Start Checklist

If you want to start using these files right now, here
is the recipe. Everything referenced here is explained
in detail in Parts 1 and 2.

1. **Gather your materials.** You need the code files
   to be changed, a description of the task, and the
   prompt files from `docs/prompts/`.

2. **Open a new LLM session.** Attach as files:
   `CCTBX_LLM_PROGRAMMING_GUIDELINES.md`, `WORKFLOW.md`,
   `ARCHITECTURE.md` (if you have one — see §1.11),
   your code archive, and `HANDOFF.json` if resuming.
   (Also attach
   `AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md` if working
   on agent code. See §2.2 for when to use which file.)

3. **Paste `WORKFLOW_PROMPT.txt`** as the first message.
   The LLM will acknowledge the rules and ask for
   `HANDOFF.json` or begin reading it.

4. **Describe the task and paste `PLAN_PROMPT.txt`.**
   Wait for the plan. Read it carefully.

5. **Optionally review the plan.** Copy the plan to a
   second LLM with `REVIEW.txt`. Feed the critique
   back. Revise once or twice.

6. **Let the LLM implement.** Verify each step (see
   §1.5). Make sure it updates `HANDOFF.json` after
   each step.

7. **If the session dies,** start a new session with the
   same attachments and paste `CONTINUE_PROMPT.txt`.

8. **Ask "what did we miss?"** Before closing out, ask
   the LLM (or a second LLM) whether there are side
   effects, callers, edge cases, or downstream
   consequences that weren't addressed. See §1.12.

For what each file contains, see the **File Reference**
table at the end of this document.

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
text they can see at once (the conversation so far, any
attached files, and their own response). Everything the
LLM needs to know must be inside that window or it
effectively does not exist.

**What to provide:**

- The specific files the LLM will read or modify. Don't
  give it the entire codebase if you only need changes
  to two files.
- An architecture document that describes the overall
  system structure, module responsibilities, and data
  flow (see §1.11). This is much more effective than
  hoping the LLM will infer the architecture from
  scattered source files.
- Coding standards and style guides. The LLM will match
  whatever conventions you show it. (On this project,
  that means attaching
  `CCTBX_LLM_PROGRAMMING_GUIDELINES.md` — see §2.2.)
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
- If the conversation is getting long and the LLM seems
  confused, start a fresh session with a clean summary
  rather than continuing. (See §1.7 on session
  management and §1.6 on conversation drift.)
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

The good request tells the LLM **what** to build,
**where** to put it, **what pattern** to follow, and
**which call sites** to update. The weak request forces
the LLM to guess all of that.

**Structure for effective requests:**

1. **What** — the change you want
2. **Where** — which files and locations
3. **Why** — the problem it solves (helps the LLM make
   reasonable choices in ambiguous situations)
4. **Constraints** — style, compatibility, patterns to
   follow or avoid
5. **Examples** — show existing code that does something
   similar, if applicable

### 1.4 The Plan–Review–Implement Cycle

Never let an LLM jump straight to writing code for a
non-trivial task. The cost of a bad plan executed
thoroughly is much higher than the cost of planning
time.

The recommended cycle:

```
START ──► PLAN ──► REVIEW ──► IMPLEMENT ──► VERIFY ──► WHAT DID WE MISS?
  │          │        │           │            │              │
Attach:    Paste:   Give plan   LLM writes  You run        Ask LLM:
 CCTBX     PLAN_    + REVIEW    code,       tests,         any side
 Guide     PROMPT   .txt to a   updates     read diffs,    effects,
 WORKFLOW  + prob.  second LLM  HANDOFF     check edges    callers,
 ARCH.md   desc.                each step                  or gaps?
 code
Paste:
 WORKFLOW_PROMPT
```

**Step 1 — Plan.** Ask the LLM to produce a written
plan: problem description, approach, implementation
steps, risks. If you have an `ARCHITECTURE.md` (see
§1.11), provide it here so the plan accounts for the
system's structure. If you don't have one and the task
is substantial, ask the LLM to draft one first — review
and correct it, then use it as input to the plan.
When you are ready to request the plan, use this prompt:

```text
Please make a plan for fixing these problems. Include
a full discussion of the problem, the overall approach,
the details of the approach, implementation plan, and
risks involved and their mitigation. write as md
```

(This is `PLAN_PROMPT.txt`.)

**Step 2 — Cross-review.** For important changes, give
the plan to a second LLM for critique. A different model
catches different blind spots. Use this prompt with the
second LLM:

```text
Gemini, you are a senior engineer tasked with
**critically reviewing this plan**. Your goal is to find
every potential flaw, unclear assumption, or workflow
violation. Be specific and constructive; do not give
generic praise.

For each step:
1. Identify issues or risks
2. Explain why it is a problem
3. Suggest improvements or alternatives

Limit to 1–2 paragraphs per step. Do not provide
general compliments.
```

(This is `REVIEW.txt`. It names Gemini but works with
any LLM — the important thing is a different model from
the one that wrote the plan.)

**Step 3 — Revise.** Feed the critique back to the first
LLM. Limit to 1–2 revision rounds — diminishing returns
set in fast.

**Step 4 — Implement.** Only now does the LLM write
code, following the agreed plan step by step. The LLM
should checkpoint after each step so that work is
recoverable if interrupted (see §1.7).

**Step 5 — Verify.** You run the code, run the tests,
and confirm correctness. The LLM cannot do this for you
(see §1.5).

**Step 6 — Ask "what did we miss?"** Before closing out,
ask the LLM to consider whether the changes have
unaddressed side effects, untouched callers, missing
test coverage, or downstream consequences. This is
especially valuable with a second LLM that wasn't
involved in the implementation. See §1.12.

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
  unnecessary, or subtly change logic while claiming to
  only reformat.
- **Check edge cases** the LLM may not have considered:
  None values, empty collections, missing keys,
  unexpected types, concurrent access.

**Project-specific checks** (see §2.4 for the full
checklist):

- **Command check:** If the code builds PHENIX commands,
  verify every parameter against the program's PHIL
  definition or `programs.yaml` strategy flags. The LLM
  will hallucinate parameters that don't exist.
- **Import check:** If the code is in `agent/`, verify
  every `libtbx.langchain` import has a
  `try/except ImportError` fallback.
- **State check:** Did the LLM update `HANDOFF.json`
  before starting this step? If not, remind it now —
  before you lose the session.

### 1.6 Common LLM Failure Modes

Understanding how LLMs fail helps you catch problems
before they reach your codebase.

**Hallucinated APIs.** The LLM invents function names,
keyword arguments, or library features that don't exist.
It does this with complete confidence. Always verify
that the API calls it writes actually match the real
signatures. (This is especially common with PHENIX —
see "Don't invent parameters" in
`CCTBX_LLM_PROGRAMMING_GUIDELINES.md` §3.)

**Confident wrongness.** When an LLM makes an error and
you point it out, it will sometimes apologize and
produce a new answer that is equally wrong but
differently worded. If a correction round doesn't
converge, re-state the problem from scratch rather than
asking the LLM to "try again."

**Sycophantic agreement.** If you suggest a wrong
approach, many LLMs will agree with you rather than push
back. Don't use leading questions like "shouldn't we use
X here?" if you're genuinely unsure. Instead ask "what
are the options for this?" and evaluate the answer.

**Drift over long conversations.** As the conversation
gets longer, the LLM gradually loses track of earlier
instructions and constraints. If you notice the LLM
forgetting your coding standards or making mistakes it
avoided earlier, it's time to start a fresh session with
a clean summary. (See §1.7.)

**Unasked-for changes.** LLMs sometimes "improve" code
you didn't ask them to touch — renaming variables,
refactoring adjacent functions, removing comments they
consider redundant. Always diff the output against the
original and reject unrelated changes.

**Incomplete changes.** The opposite problem: the LLM
changes the function but forgets to update the callers,
the tests, or the documentation. Ask explicitly: "what
else needs to change to keep everything consistent?"

### 1.7 Session Management

LLM sessions are inherently fragile. They can be
interrupted by context length limits, timeouts, network
issues, or simply because you closed the browser tab.
Plan for this.

#### Checkpointing with HANDOFF.json

After each meaningful step, the LLM writes a checkpoint
file that captures: what was done, what remains, what
decisions were made, and what the next step is. On this
project we use `HANDOFF.json`. A blank template:

```json
{
  "current_task": "UNKNOWN",
  "status": "not_started",
  "relevant_files": [],
  "completed_steps": [],
  "remaining_steps": [],
  "decisions": [],
  "assumptions": [],
  "open_questions": [],
  "last_stable_state": null,
  "next_step": "analyze task"
}
```

The LLM fills in and updates every field as work
progresses. The critical property: **a new LLM session
must be able to resume using only `HANDOFF.json` and
the code, with no other context about what happened
before.** The full rules for when and how to update it
are in `WORKFLOW.md` (see §2.1).

If the LLM is not updating `HANDOFF.json` after each
step, remind it. This is the most common workflow
violation.

#### Session length

Shorter sessions with clear handoffs are more reliable
than marathon sessions. After about 15–20 back-and-forth
exchanges, consider wrapping up and starting fresh. Each
new session gets a clean context window and a fresh start
on the instructions you provide.

#### Resuming after interruption

When a session dies (timeout, context limit, crash),
here is the exact recovery procedure:

1. Open a new LLM session.
2. Attach the same reference documents as before
   (`CCTBX_LLM_PROGRAMMING_GUIDELINES.md`,
   `WORKFLOW.md`, and
   `AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md` if
   applicable).
3. Attach the current `HANDOFF.json` and the code
   archive.
4. Paste this prompt as the opening message:

```text
You were interrupted.

Do NOT continue from memory.

1. Read the attached HANDOFF.json to understand the
   state.
2. Read the code archive to see what was actually
   implemented.
3. Summarize where we are and what the very next
   micro-step is.

If anything is unclear or missing, ask before
proceeding.
```

(This is `CONTINUE_PROMPT.txt`.)

The LLM will read `HANDOFF.json`, tell you where things
stand, and propose the next step. Confirm before it
proceeds.

#### Starting a new session (not a resume)

Always provide:

1. The reference documents (coding guidelines, workflow
   rules — see §2.2 for which to attach)
2. The current `HANDOFF.json` (if any prior work exists)
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
different LLM critique it. They make different kinds of
mistakes, so the reviewer will often catch things the
author missed. This is the most effective multi-LLM
pattern. The `REVIEW.txt` prompt (shown in §1.4) is
designed for exactly this — give it plus the plan to
the second LLM.

**Stuck debugging.** If one LLM has gone in circles on
a bug, describe the problem fresh to a different model.
A different set of biases often finds a different (and
sometimes correct) path.

There is no need to use multiple LLMs for routine
coding tasks. The overhead of managing multiple sessions
isn't worth it for straightforward changes.

### 1.9 What LLMs Are Good and Bad At

**Good at:**

- Boilerplate and repetitive code (serialization, CRUD,
  test scaffolding)
- Translating a clear specification into code
- Explaining unfamiliar code or error messages
- Suggesting approaches you hadn't considered
- Drafting documentation, docstrings, and comments
- Refactoring with a well-defined pattern (e.g., "add
  null-safety to all `.get()` calls in this file")

**Bad at:**

- Reasoning about complex state across many files
- Performance optimization (it optimizes for
  readability, not speed, unless told otherwise)
- Knowing the actual runtime behavior of code it writes
- Maintaining consistency over long sessions
- Understanding your specific codebase's unwritten
  conventions (unless you tell it)
- Anything requiring information that isn't in its
  context window

### 1.10 Code Review for LLM-Generated Code

Review LLM-generated code the same way you would review
a pull request from a new team member who is talented
but unfamiliar with your project:

- Does it follow the project's style and conventions?
- Does it handle the edge cases that matter for this
  codebase?
- Did it introduce any unnecessary dependencies?
- Did it change anything beyond what was requested?
- Are the names clear and consistent with existing code?
- Would you be comfortable maintaining this code six
  months from now?

If the answer to any of these is "no," send it back for
revision — just as you would with a human author.

### 1.11 Supplying an Architecture Document

LLMs work from the files you give them. If you hand them
three source files with no explanation of how those files
fit into the larger system, they will guess — and they
will guess wrong. An architecture document eliminates the
most damaging category of guessing: guessing about
structure.

**What to put in `ARCHITECTURE.md`:**

- The major components or modules and what each one is
  responsible for.
- How data flows between them (e.g., "the PLAN node
  produces a command list that the BUILD node executes").
- Which parts are client-side vs server-side, if
  applicable.
- Key design constraints (e.g., "agent/ must not import
  server-only dependencies").
- Any non-obvious coupling (e.g., "changes to
  `session.data` fields must also appear in
  `create_initial_state()` and `contract.py`").

This does not need to be long. A one-page document with
a component list and a data flow description is far more
useful than nothing. The goal is to prevent the LLM from
making structurally wrong changes — like putting
server-only logic in client code, or modifying a file
without updating the three other files that must stay in
sync with it.

**When to provide it:**

- Attach `ARCHITECTURE.md` alongside the coding
  guidelines at session start, especially when the task
  involves multiple files or crosses module boundaries.
- If you don't have one yet and the project is
  substantial, **create it as part of the planning
  process** (§1.4 Step 1). Before asking for an
  implementation plan, ask the LLM to draft an
  architecture document from the source files you've
  provided. Review and correct it, then use it as input
  to the plan. This pays for itself immediately — the
  plan will be structurally sound instead of built on
  guesses about how the system fits together.
- For smaller or one-off tasks, you can skip this. But
  if you find yourself explaining the same system
  structure to the LLM across multiple sessions, that's
  a sign you need an `ARCHITECTURE.md`.

**When to update it:**

- After any task that changes the system's structure
  (new modules, new data flow paths, changed
  responsibilities). If the architecture document is
  stale, it will actively mislead future sessions.

### 1.12 The "What Did We Miss?" Check

After the code is written and the tests pass, there is
one more step before you close out: ask the LLM whether
anything was missed.

This catches a class of problems that verification
(§1.5) does not: side effects, downstream consequences,
and unstated assumptions that neither you nor the LLM
thought to check. It is cheap to do and occasionally
catches expensive mistakes.

**How to do it.** After the implementation is complete
and tests pass, ask:

> We just made these changes: [brief summary or diff].
> Are there any side effects, callers, downstream
> consumers, edge cases, documentation, or tests that
> we haven't addressed? Think carefully about what
> could break that we haven't considered.

**Why a second LLM is especially useful here.** The LLM
that wrote the code has a blind spot: it has already
convinced itself the approach is correct. A fresh LLM
(or a fresh session) has no such commitment. Give it
the changed files, the original files, and the question
above. It will often spot things the implementing LLM
overlooked — a caller that passes None, a serialization
path that expects the old field name, a test file that
hard-codes an assumption the change just invalidated.

**What to look for in the answer:**

- Callers or importers of changed functions that weren't
  updated
- Serialization round-trips (`to_dict` / `from_dict`)
  that need new fields
- Configuration or YAML files that reference changed
  names
- Documentation that describes the old behavior
- Test files that test the old behavior or hard-code
  values that just changed

You don't need to act on every suggestion — the LLM may
flag things that aren't actually problems. But reviewing
the list takes a minute and occasionally saves hours.

---

## Part 2: Project Workflow

### 2.1 Rules of Engagement

When you start a session using `WORKFLOW_PROMPT.txt`
(the full text is in the Appendix), the LLM is
instructed to follow these rules. Knowing them helps you
understand why the LLM behaves the way it does — and
when to correct it.

**Small steps.** The LLM must modify at most 3 files
per step unless it justifies more. Work is broken into
micro-steps: understand the problem, identify files,
create a plan, implement one change, validate, repeat.

**No placeholders.** The LLM must provide full copies of
all changed files, never diffs-only or
`# ... rest of code stays the same ...` comments. You
should be able to drop the output directly into the
codebase.

**Mandatory checkpoints.** The LLM must update
`HANDOFF.json` after every micro-step — after reasoning,
after each code change, before running tests, after test
results, and before ending its response. This is the
rule that gets violated most often. If you notice the
LLM batching multiple changes without checkpointing,
remind it.

**Plan before code.** The LLM must write a plan (max
6–8 steps, each independently executable) and have it
reviewed before writing any code. See §1.4.

**Archive output.** When delivering changes, the LLM
produces a tar.gz archive preserving the original
directory structure, plus the updated `HANDOFF.json` and
a short explanation for each modified file.

The full workflow rules are in `WORKFLOW.md`, which the
LLM reads at session start.

### 2.2 Which Guideline Files to Attach

Attaching the right guideline files gives the LLM the
conventions and pitfalls specific to the code it will
touch. Attaching the wrong ones wastes context window
space. Here is when to use each:

**Always attach (every session):**

- **`CCTBX_LLM_PROGRAMMING_GUIDELINES.md`** — The base
  layer. Covers coding standards (2-space indent,
  80-char lines, naming conventions), cctbx/PHENIX
  patterns (PHIL parameters, flex arrays, space group
  conventions, import patterns), common pitfalls
  (None-safety on `dict.get()`, R-free as fraction not
  percentage, hardcoded path separators), error handling
  rules ("never raises" pattern, no bare excepts), and
  a pre-submission checklist. Use for any code involving
  crystallographic logic, PHIL parameters, or standard
  PHENIX utilities.

- **`WORKFLOW.md`** — Defines the checkpoint discipline,
  plan format, interruption protocol, and output
  requirements. This is what makes sessions recoverable.

- **`ARCHITECTURE.md`** (if you have one) — Describes
  the system's component structure, data flow, and
  design constraints. Especially important when the task
  crosses module boundaries or involves multiple files.
  See §1.11.

**Attach when working on agent code:**

- **`AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md`** — The
  agent layer. Use only when working within `agent/`,
  `knowledge/`, or `programs/ai_agent.py`. Contains
  strict rules about parameter verification against
  `programs.yaml` (the command sanitizer strips unknown
  flags silently), `libtbx.langchain` import fallbacks
  (both paths must resolve — this is different from the
  general cctbx pattern), the three error classification
  systems (which must agree when you add a new error
  pattern), client vs server code path awareness,
  `session.data` persistence (every new state field must
  appear in three places), and agent-specific test
  patterns.

**Don't over-attach.** If you are fixing a bug in
`iotbx/pdb/`, you don't need the agent guidelines. If
you are modifying `agent/error_classifier.py`, you do.
Extra guideline files cost context window space that
could hold more of the actual code the LLM needs to see.

### 2.3 The Workflow Visualized

Here is the full cycle showing which prompt to use at
each stage and what the human does at each step:

```
  ┌─────────────────────────────────────────────────┐
  │               START SESSION                     │
  │                                                 │
  │  Attach as files (do NOT paste):                │
  │    CCTBX_LLM_PROGRAMMING_GUIDELINES.md          │
  │    WORKFLOW.md                                  │
  │    ARCHITECTURE.md (if you have one)            │
  │    AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md       │
  │      (only if working on agent code)            │
  │    code archive                                 │
  │    HANDOFF.json (if resuming)                   │
  │                                                 │
  │  Paste into chat:                               │
  │    WORKFLOW_PROMPT.txt                           │
  └────────────────────┬────────────────────────────┘
                       │
                       ▼
  ┌─────────────────────────────────────────────────┐
  │                  PLAN                           │
  │  Paste into chat:                               │
  │    PLAN_PROMPT.txt + problem description        │
  │  Human: Read the plan. Does it make sense?      │
  └────────────────────┬────────────────────────────┘
                       │
                       ▼
  ┌─────────────────────────────────────────────────┐
  │            REVIEW (recommended)                 │
  │  In a DIFFERENT LLM session:                    │
  │    Give the plan + REVIEW.txt                   │
  │  Feed critique back to first LLM.               │
  │  Revise 1–2 rounds.                             │
  └────────────────────┬────────────────────────────┘
                       │
                       ▼
  ┌─────────────────────────────────────────────────┐
  │               IMPLEMENT                         │
  │  LLM executes plan step by step.                │
  │  LLM updates HANDOFF.json after each step.      │
  │  Human: verify each step (§1.5, §2.4)           │
  └────────────────────┬────────────────────────────┘
                       │
                       ▼
  ┌─────────────────────────────────────────────────┐
  │          WHAT DID WE MISS? (§1.12)              │
  │  Ask the LLM (or a second LLM):                 │
  │    any side effects, callers, edge cases,        │
  │    or downstream consequences we missed?         │
  │  Review the answer. Act on real issues.          │
  └────────────────────┬────────────────────────────┘
                       │
              ┌────────┴────────┐
              ▼                 ▼
  ┌──────────────────┐  ┌──────────────────────────┐
  │      DONE        │  │     INTERRUPTED           │
  │  HANDOFF.json    │  │  New session              │
  │  reflects final  │  │  Same file attachments    │
  │  state.          │  │  Paste CONTINUE_PROMPT    │
  │  Update          │  │  (see §1.7)               │
  │  ARCHITECTURE.md │  │                           │
  │  if structure    │  │                           │
  │  changed.        │  │                           │
  └──────────────────┘  └──────────────────────────┘
```

### 2.4 The Human Navigator's Checklist

Your job is verification and navigation. The LLM cannot
do these things for itself. Use this checklist to make
sure you are holding up your end.

**At session start:**

- [ ] Did I attach `CCTBX_LLM_PROGRAMMING_GUIDELINES.md`
      and `WORKFLOW.md`?
- [ ] Did I attach `ARCHITECTURE.md` (if I have one)?
- [ ] Did I also attach
      `AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md` (if
      working on agent code)?
- [ ] Did I attach `HANDOFF.json` (if resuming)?
- [ ] Did I provide only the relevant code files, not
      the entire codebase?

**After the LLM produces a plan:**

- [ ] Did I read the plan and check it makes sense?
- [ ] For non-trivial changes, did I send the plan to
      a second LLM for review?

**After each implementation step:**

- [ ] Did I run `ast.parse()` (or equivalent) on every
      changed file?
- [ ] Did I run the tests that exercise the changed
      code?
- [ ] Did I read the diff and reject unrelated changes?
- [ ] Did the LLM hallucinate any PHENIX parameters?
      (Check against the actual PHIL definition or
      `programs.yaml` strategy flags.)
- [ ] Is the LLM keeping `HANDOFF.json` current? If
      not, remind it.

**Agent-specific checks (when working on agent code):**

- [ ] If a new state field was added, did the LLM
      update `create_initial_state()`, `session.data`,
      and `contract.py`?
- [ ] If a new error pattern was added, does it appear
      in all three classification systems?
- [ ] Do all `libtbx.langchain` imports have fallback
      `except ImportError` blocks?

**Before closing out:**

- [ ] Did I ask the LLM (or a second LLM) "what did we
      miss?" — any side effects, callers, edge cases,
      or downstream consequences? (§1.12)

**At session end:**

- [ ] Does `HANDOFF.json` reflect the final state?
- [ ] Could a new session resume from `HANDOFF.json`
      alone?
- [ ] Did the LLM provide full file copies (not just
      diffs)?
- [ ] If the system's structure changed, did I update
      `ARCHITECTURE.md`?

---

## Appendix A: File Reference

All files live in `docs/prompts/`.

**Attach as files** — these are large reference
documents. Do NOT paste their contents into the chat;
that crowds out the LLM's reasoning space and leaves
less room for your actual code.

| File | Size | Purpose | When to attach |
|------|------|---------|----------------|
| `CCTBX_LLM_PROGRAMMING_GUIDELINES.md` | ~650 lines | Coding standards, cctbx patterns, pitfalls, checklist | **Every** session |
| `WORKFLOW.md` | ~250 lines | Checkpoint rules, plan format, interruption protocol | **Every** session |
| `ARCHITECTURE.md` | varies | System structure, module responsibilities, data flow, design constraints | **Every** session (if you have one). See §1.11. |
| `AI_AGENT_LLM_PROGRAMMING_GUIDELINES.md` | ~310 lines | Agent-specific patterns: imports, state persistence, error systems | Only when working on **agent code** |

**Paste into chat** — these are short prompts you type
or paste at specific moments during a session. They are
small enough that pasting is fine.

| File | Size | Purpose | When to paste |
|------|------|---------|---------------|
| `WORKFLOW_PROMPT.txt` | ~55 lines | Start a session; sets checkpoint rules | As the **first message** (full text in Appendix B) |
| `PLAN_PROMPT.txt` | 2 lines | Ask the LLM for an implementation plan | When **requesting a plan** (full text in §1.4) |
| `REVIEW.txt` | ~8 lines | Ask a second LLM to critique a plan | Into a **different LLM** with the plan (full text in §1.4) |
| `CONTINUE_PROMPT.txt` | ~11 lines | Resume after interruption | When **resuming** a dead session (full text in §1.7) |

**Carried between sessions** — the LLM creates and
maintains this file; you save it and provide it to the
next session.

| File | Size | Purpose |
|------|------|---------|
| `HANDOFF.json` | ~12 lines (blank template) | Session state: what was done, what remains, what to do next. Template in §1.7. |

---

## Appendix B: WORKFLOW_PROMPT.txt (full text)

This is the opening message you paste to start a
session. It is the longest of the session prompts
(~55 lines), so it is kept here rather than inline
in Part 1.

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
