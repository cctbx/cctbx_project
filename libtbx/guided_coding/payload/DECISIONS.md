# The decision table

Who answers what, and what happens if nobody says anything. The Worker consults this before putting ANY question to the developer. A question the developer cannot answer is a defect in the question.

**Precedence.** This table ROUTES questions; it does not override a rule stated elsewhere. When it and another document differ, the specific rule wins and you say which you followed. In particular: `/wrap`'s removal rules govern anything of the DEVELOPER's (subagent copies from other changes, their untracked files); the workflow and VERIFY.md name decisions that are the developer's but are not listed in column B below - setting aside a developer file, accepting an environment update, accepting a weaker fallback - and those are asked exactly as column B items are. Column B is a floor, not a ceiling.

**The test, applied to every question before you ask it:** can this be answered from the traceback, the source, the profile in CLAUDE.local.md, the plan, or this table? If yes, answer it yourself and say what you concluded. Ask only what genuinely needs the developer's authority or knowledge.

## A. The Worker answers these. Never ask.

| Question | How you answer it | Say this |
|---|---|---|
| Which repository does this change belong in? | Start from the traceback's file path, then check the callers and callees: a function can fail because another repository handed it a bad value, so where it SURFACED is not always where it BELONGS. Read enough to say which, and preserve the uncertainty if the diagnosis is not yet settled | "The crash surfaces in cctbx_project; I have not yet established whether the defect belongs there or in the caller" - and if the session is in the wrong repository for the change as diagnosed, say so |
| Which tests exercise the changed code? | `phenix.find_program` per changed file, plus callers | Name them and their number |
| Where does a new test belong? | The three tiers in the workflow, using the project's own tools | State the tier and why the better tier did not fit |
| Is the tree clean / where is master / what moved upstream? | Read it | Report as fact, not as a question |
| Does a colleague's upstream commit affect this change? | Check overlap with files touched, called or tested | Speak only on overlap; otherwise log one line |
| Should I clean up my own scratch copies, worktrees, temporary files? | Yours from THIS change only: it is your housekeeping. Anything belonging to the developer or to another change follows /wrap's removal rules - ask | Do your own, report it; name the others and ask |
| Which of my own commands need approval? | The settings; construct commands so routine work does not prompt | A surprise prompt is a defect of yours |
| What does a term mean? | Use plain words; the term once in parentheses | Never assume shorthand is shared |

## B. The developer answers these. Always ask, in the four-line shape.

| Question | Why it is theirs | Default if they say nothing |
|---|---|---|
| Is this diagnosis right, and is this the change to make? | The plan gate | No default - wait |
| What SHOULD the code do here? (meaning, not mechanism) | Their code, their science | No default - wait |
| Integrate or discard? | Their tree | No default - wait |
| Publish, and with what conditions overridden? | Their colleagues, their users | No default - wait |
| Is anything else using this tree right now? | Only they know | No default - wait |
| May I expand scope beyond the approved plan? | Scope is theirs | No - proceed as planned |
| Should a failed condition be overridden? | Theirs to accept | No - hold, and say what is missing |

## C. Rules, not questions. Apply them; do not ask.

| Situation | The rule |
|---|---|
| Do the case-specific registered tests always run? | ALWAYS, both sides, before the integrate decision. WHERE they run is the next row's judgment: locally by default, on the server when the local cost is large |
| What if they will take more than ~10 minutes locally? | Say so, give the estimate, and propose the server instead. Do not silently start a long local run |
| When does the full server suite (t96) run? | For PUBLICATION only. Never for a local integration |
| When are the Guide and the Helper involved? | Read the arrangement from the profile; never ask "is a Helper active?" - a first-time developer cannot answer that. The Guide reads every report; a Helper reading of the proof packet is required before publication unless the developer records a waiver |
| May I tell the user their bug is fixed? | Not until it is PUBLISHED. On a saved-locally outcome, say the fix is in testing |
| Is a proof packet needed? | Assemble and have it checker-signed whenever publication is pending, even if publication is deferred |
| Something went wrong or cannot be classified | `.claude/INCIDENT.md` governs. Freeze, inspect, preserve, report, recover only as authorized |
| A procedure file needs changing | Never edit it. Record the finding; it ships in the next package revision |
| The developer asks a question mid-change | Answer it, then re-issue the complete current set of commands so the latest message alone says what to do |

## D. How anything in column B is asked

Four lines, always, in this order:

1. **What I want to do** - in plain words, the action and its purpose
2. **What you will notice** - the visible effect, on which machine
3. **What cannot go wrong, and what could** - the real risk, checked not guessed; uncertainty preserved
4. **What I need from you** - the decision, with concrete options, **the default if they say nothing**, and your recommendation with a one-line reason

And: you do the work and bring the decision. Never hand the developer a task list of things you could have executed.
