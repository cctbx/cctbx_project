# The roles

**2026-09-26, adopted.** Five roles. The first question is not what they are,
it is when you need one.

**Adoption history.** The 2026-09-26 version superseded 2026-09-17,
SHA-256 `a25195ded09cb6f23c29fb29e195978e81f432bc67968234ad437a419cb17e9c`,
by adding two sentences after the Helper paragraph in "When do you need
each one?" The 2026-10-03 documentation revision corrects the obsolete
installed-name note below and removes a historical review-session name
from the general role description. Role authority is unchanged; the
earlier adoption date is not an approval of this documentation revision.

**Parent.** This document and its three briefs specialize
`DEVELOPER_GUIDE_CONTRACT.md` version 2026-09-17, SHA-256
`ab4586810ef683702c8648267d5475fb50659b9fada4b2409293e93a87e6876b`.
They add constraints and mechanics. Where they touch the contract, the
contract governs, and none of them may move a decision the contract
reserves.

## When do you need each one?

**A Helper, when you do not understand.** It is *your* Helper. It reports
to no one else, it is briefed by no one else, and nobody else sees what
you ask it. Use its advice however you like — that is the point of having
it. What it is not: a gate result, and never an input that reaches the
work except through a decision you make.

Before approving a working window's result, try to restate what it did
and trace one central claim to its source. If you cannot, or you need the
options and their costs explained for a decision you own, a Helper can
work through the evidence with you.

You may open a separate Helper conversation at any time, including while
a Guide or Worker has work under way. Opening it does not transfer their
work to the Helper; commitments that depend on your decision still wait.

**A Guide, when you do not know how to do it, or you want to delegate
the work or the authority.** A Guide may decide **anything you have not
reserved** — with *your* authority. That is why it is a Guide and not a
Leader: what it decides, it decides as you, it records as yours, and you
can reverse.

Six decisions are reserved to you by the contract and are not delegated:
**meaning, value, risk, waiver, acceptance, and consequential action.**
Meaning is the conclusion you adopt — anyone may read the evidence and
recommend one. Consequential action covers integration, publication,
anything claimed in public, and anything else whose effects reach beyond
the working copy.

Risk and waiver are on that list. Deciding to accept a risk, or to waive
a requirement, is not delegable however convenient it would be. Carrying out a publication you have
specifically authorized is not the same as deciding to publish; the first
is delegable and the second is not. Everything else may be delegated by
naming it.

**A Worker, when the job is bounded.** It needs an outcome, a scope,
acceptance criteria, and a named range of implementation choices it may
make on its own. It does not need to be fully specified, because no job
is: new facts appear while the work is done.

**What a Worker does when it meets a choice nobody covered.** It may
proceed only inside the discretion it was given, and inside the agreed
limits on effects and recovery. A choice that changes the meaning, the
scope, an acceptance criterion, a public commitment, or those limits goes
back to you, unless an authorization already in force settles it. It
stops the affected action, says what the missing decision is and what
each option would cost, and asks the smallest question that lets it
continue. Unrelated safe work carries on. An isolated experiment may be
recorded as provisional, with its assumption, its owner and how to undo
it — and recording it does not make it adopted. Every provisional choice
is adopted or rejected by whoever has the authority, before the work that
rests on it is accepted.

**An Outside Reviewer, when you want something evaluated
independently.** Not explained, not improved — evaluated, by a reader who
did not make it and has no stake in it passing.

---

## The five roles

Five roles. Their authority and their deliverables
are distinct, even where some activities overlap: the Guide and the
Helper both explain, three different readers all review, and the
Developer and the Guide both reason about choices. What must never
overlap is who may decide and who owns which output.

**Names.** Records written before 2026-09-17 use "Helper" for what is
here called the Outside Reviewer. Records before 2026-09-14 use "Guide"
for the Worker and "Helper" for the Guide. This package uses **Outside
Reviewer** at the gate and **Helper** for the Developer's own explainer.
The old parenthetical "called Helper in the installed procedure" belongs
only to those historical records.

| Role | What it is | What it does | What it never does |
|---|---|---|---|
| **Developer** | The person. | Chooses scope and priority. Keeps six decisions: meaning, value, risk, waiver, acceptance, and consequential action. Others may interpret and recommend. Delegates anything else by naming it, and can take it back. Runs the pastes and the blocks. Owns his personal profile. Carries material between the other four, verbatim. | Edits procedure files by hand. Approves what he cannot restate. Summarizes what he carries. Delegates one of the six reserved decisions. |
| **Guide** | The session that plans and directs. It may be a Claude Code session with hands in the repositories and on the servers, or a chat with hands only in its own environment. Which it is depends on where the project runs, not on what the role is. | Keeps the ordered queue and recommends the next change. Establishes facts it can reach. Directs the work: writes the front door for a separate Worker session, or briefs a bounded subagent. Reads what comes back. Prepares what goes to the Outside Reviewer. Builds and tests releases. Keeps the handoff and the release notes. | Assesses at a gate, or stands in for a reviewer of anything it produced or planned. It checks delivery, which is not a gate. Repairs a defect it finds in someone else's change instead of sending it back. Decides any of the six reserved decisions, decides anything else the Developer has not delegated by name, or keeps a delegated decision out of the record. Integrates or publishes without the Developer's quoted word. Produces the evidence for its own plan without saying so. Briefs, messages or sees the Helper. |
| **Worker** | One bounded job with hands. Either a separate Claude Code session in one repository, or a subagent the Guide briefs. Its reading subagents belong to it. | Investigates, plans, builds in isolated copies, verifies in the live tree with save-and-restore, and stops. What it must produce depends on the path the change is assigned: the full path takes one review from its reviewer subagent, a frozen proof packet, and a checker subagent's reading; the light path takes an evidence archive instead. Integrates or pushes only on the Developer's quoted authorization, re-verifying that the authorization still applies. Works inside its named discretion; stops the affected action when a choice would change meaning, scope, an acceptance criterion, a public commitment or the agreed limits. A separate session's stops go to the Developer directly; a subagent's reach him through the Guide. | Integrates or publishes without the quoted word. Edits procedure files or the profile. Decides meaning, scope or acceptance. Upgrades the procedure. |
| **Outside Reviewer** | A reader outside the repositories and outside the Guide's conversation, with no hands on the Developer's machines; a different model where possible. | Reads at three gates: a full-path packet before local integration, a publication-batch packet before the push, and a release kit before it is installed for Workers to use. Verifies the bundle and the packet, forms its own view before reading the verdict, runs the shipped test, probes what the code refuses, says what it did not reach, and says whether the thing may proceed and on what condition. | Repairs anything. Changes the Developer's systems or the artifact under review — its own throwaway scratch is its to use freely. Accepts a description of evidence in place of the evidence. Takes the Developer's authority as its own. Reopens the architecture because it would have designed it differently. |
| **Helper** | The Developer's own window: a conversation he keeps for his own understanding, of whatever model he chooses. Started whenever he wants one, often mid-stream from the record alone. | Explains what is being done, in his words. Says what a decision actually is and what each option costs. Asks for his reading first when he has one. Tells him when his questions are one-sided. Available at any time about anything. | Project execution of any kind: code, plans, packets, releases. Any decision. Any gate. Reporting to, being briefed by, or speaking for anyone but him. Filling in what it was not sent. Holding the only copy of anything durable. |

## Authority, and taking it back

A delegation is a **grant**, defined in the contract at §9: what it
covers, where it applies, when it ends, how its effects would be
recovered, and whether it is still in force. It is **filed where the
project keeps its approvals** — not held in a session's memory — and a
session **re-reads it** before any consequential action. What this
document adds is where revocation has to reach.

Taking a grant back does not reach a session by itself. **A revocation
is an artifact**, and travels the way artifacts do: you choose that it
passes, it passes exactly, and where an exact transfer can be automated
it may be. It has to reach every session working under the grant,
including a subagent, which hears it through its Guide. Before integrating, publishing, or any other
consequential action outside its own workspace, a session checks that its
authorization still applies. If it cannot establish that, it stops that
action and keeps its work.

## Can it be undone?

The contract's rule at §9: recoverable is not "there is an undo command",
and a grant states how its effects would be recovered and what bound on
consequence and recovery cost is acceptable. The question when choosing
is whether the effects can be recovered inside that bound.

## To start each session

**Helper:** paste `HELPER.md`, then the record and the message in
front of you.
**Guide:** paste `GUIDE.md` and the current handoff.
**Worker:** the Guide's front door, which names the repository, the
bounded outcome, the discretion, and the endpoint.
**Outside Reviewer:** `OUTSIDE_REVIEWER_BRIEF.md` with the bundle and
the Guide's message.

## Who reviews what

Nobody reviews their own work. The Guide **checks delivery** — did the
decisive check fail and then pass, were the runs in order, was every
finding dispositioned, what was settled that nobody assigned — and
coordinates the corrections. That is not a gate and is never recorded as
one.

**Independent assessment at a gate belongs to the Outside Reviewer, and
to nothing the Guide produced or planned.** A Worker's change is reviewed
by its own reviewer subagent before the packet is frozen, by its checker
subagent after, and by the Outside Reviewer at the gate.

These are three chances to catch an error at three different moments.
They are not three independent votes: the first two are the same model
family as the Worker, and the third may be. Do not read agreement among
them as confirmation. A release the Guide built is reviewed by the Outside Reviewer. The
Guide reads what comes back and tells the Developer what it sees; that
reading is not a gate and is never recorded as one.

## Three readers, kept apart

The Worker's **reviewer** subagent reads the candidate before the packet
is frozen. The Worker's **checker** subagent reads the frozen packet
against a checklist. Both are inside the change and the same model as
the Worker. The **Outside Reviewer** reads at the three gates above, from
outside, and is a different model where possible.

The Helper is not a reader at any gate. It may read the same packet with
the Developer, and its opinion is never the gate result.

## Which changes reach the Outside Reviewer

A full-path change's packet, at its gate. A publication batch's packet,
before the push. Every release kit. A light-path change — one the
workflow sorts as low-consequence by what it touches, never by size —
carries an evidence archive instead of a packet, and reaches the Outside
Reviewer only inside the batch that publishes it.

## Who talks to whom

The Developer decides what passes between windows, and it passes
**verbatim**. What keeps the Guide, the Outside Reviewer and the Helper
independent is that they share no context: only the named artifact
crosses, and a summary on the way defeats that. His hands are not what
creates the independence, so where an exact transfer can be automated it
may be — the requirement is that nothing but the artifact passes, and
that he chooses what passes.

The Guide talks to the Developer; through him it starts Workers and
briefs the Outside Reviewer. A separate Worker session talks to the
Developer directly; a subagent cannot, so its questions reach him through
the Guide, which says which kind it was. The Outside Reviewer's reply
reaches both the Worker and the Guide through him, verbatim.
Nobody talks to the Helper but the Developer, and nobody recruits it but
him. He may recruit any role at any time.

## What follows from the Guide having hands

Claude Code is a Guide, and it has hands. The earlier definition — a
Guide as a chat with no hands — was taken from one project where the
Guide happens to be a chat window. That was a fact about that project,
not about the role.

Two things follow, and the second is a real cost.

**Facts stop travelling through the Developer.** A Guide that can look
does not ask him to paste. It states where and when it measured
anything, in the shell that will use it.

**The session that plans a change may also produce its evidence.** That
is the arrangement with a stake in its own approval. It is not forbidden,
because forbidding it would mean giving up hands, but it has to be
visible and it has to be bounded:

- The **decisive check and the evidence** are produced by something that
  did not plan the change: a subagent with a bounded brief, or a separate
  session. A script that reruns the same way for anyone gives
  reproducibility, not independence.
- The packet records **who chose the criterion, who wrote the test, and
  who ran it.** Where none of them was independent of the planning, the
  packet says so and the Outside Reviewer challenges the criterion rather
  than only rerunning the test.
- Whenever the Guide produced the evidence for its own plan, the packet
  and the message **say so**.
- When that is true, the **Outside Reviewer matters more, not less**, and
  the Guide never argues it out of a finding.
