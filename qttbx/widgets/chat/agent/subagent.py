"""Run a child conversation and return its output to the parent.

The child is a full :class:`AgentSession` at ``parent.depth + 1`` with its own
agent, tool registry, policy, conversation and event sink. Four properties are
load-bearing and each has a reason that is not obvious from the code:

- **Its own agent.** ``TokenUsage.context_tokens`` is a monotonic peak kept per
  agent, so reusing the parent's would permanently inflate the parent's context
  reading and fire a false handoff nudge.
- **Its own event sink.** A child ``TurnDone`` delivered to the parent's handler
  unlocks the parent's UI mid-turn, after which the user's next message is
  silently dropped by the busy guard. The sink is not a no-op, though: it is
  the ONLY place the child's terminal disposition is observable (see
  :class:`_TerminalOutcome`).
- **Its own persistence call.** ``AgentSession`` disables autosave and
  checkpointing at ``depth > 0``, so a child transcript never reaches disk
  unless this module stores it.
- **Its own self-answering approvals.** Nothing downstream of a private sink
  can ever answer a ``ToolApprovalRequest``, so a child tool resolving to
  ``ask`` would park the parent's worker forever. This module owns the sink, so
  it owns the answer -- see :class:`_HeadlessApprovals`.

And one property the whole module exists for: **a child that was cut off must
never be mistakable for a child that ran fine and found nothing.** Cut off
covers the turn cap, cancellation, a crash, a truncated answer, a provider
block, a turn that ended with NO terminal event at all, and being denied every
MEASUREMENT tool. Each lands on ``incomplete_reason`` (what code reads) AND as
prose in ``final_text`` (all the parent MODEL ever sees).

Whose cap fired matters, because the two backend families are capped in
different places. ``AgentSession.run_turn(max_turns=...)`` caps the API
backends and reports it through a transcript marker; ``claude_code`` runs its
agentic loop inside the SDK subprocess, which reports its own cap as
``TurnDone(finish=CAP)`` and never writes that marker. Reading only the marker
therefore made a capped child on the DEFAULT backend byte-identical to a clean
one. Both signals are read here, and the ``finish`` half is read off the
canonical, backend-agnostic enum rather than off ``stop_reason`` strings.

The same split runs through the other two zero-measurement signals, and each
had a half nothing was reading:

- **The terminal event can be ABSENT.** Every backend surfaces a provider
  failure as ``yield AgentError(...)`` followed by a return -- no ``TurnDone``
  -- so ``AgentSession`` leaves ``stop_reason`` empty, the sink latches nothing
  and a table of known-bad finishes has nothing to look up. Measured: a review
  killed by a rate limit came back byte-identical to one that reviewed
  everything and found nothing. ``_TerminalOutcome`` therefore records whether
  a terminal event arrived AT ALL, not merely what it said.
- **Whose permission layer denied.** The API backends route every approval
  through ``_HeadlessApprovals`` here; ``claude_code`` resolves permissions
  inside its own SDK callback, where nothing this module owns is consulted. So
  the denied-every-tool backstop was dead on the DEFAULT backend. ``run_child``
  hands such an agent a denial sink (``set_denial_sink``, duck-typed so no
  backend is imported) and reads the union.

What counts as a MEASUREMENT is the third half, and it is asked as ONE question
-- *did any tool that can actually measure return a real result?* -- rather
than per failure shape. Four shapes reach the same place and only one of them
was ever caught: denied, errored (an unknown tool, a server that never
started), never called at all, and answered successfully by a tool that
measures nothing. The old gate gated on a DENIAL record, so the other three
walked past it; the gate is now the outcome itself, with the denial kept only
to say WHICH tools were refused. Two structural exclusions decide what counts,
both read off the wiring the child was actually built from so neither can drift
(:func:`_is_measurement_tool`):

- the skill-introspection tools, pre-allowed on every child so an on_demand
  skill can deliver its operating instructions -- they read the child's own
  prompt material. Read off the registry (``ToolRegistry.source_of``), the same
  signal the phenix-side factory uses to choose what to pre-allow, never a
  hardcoded name list that a renamed or added skill tool walks past;
- any tool from an MCP server the child's PROFILE did not ask for. A backend
  adds servers of its own: ``claude_code`` builds an in-process ``phenix_chat``
  server on every agent, and its ``phenix_get_job_history`` -- auto-approved
  ahead of every policy check, and returning ``"(no jobs recorded)"`` as a
  SUCCESSFUL result when no provider is wired, which is every child -- made
  "did any tool succeed" true for a child that had measured nothing. The
  bundle's ``measurement_servers`` is what the profile really wired up, so the
  whole class is excluded at once without naming a single tool.

The flag itself is gated on the child having been EQUIPPED to measure
(:func:`_child_was_equipped`): a child that held measurement tools and got
nothing back from any of them measured nothing, however that came about, while
a child that was never equipped is a build-time fault the phenix-side factory
refuses to spawn at all.

One profile key, one meaning. ``subagents.default_max_turns`` sits in the
``subagents`` block beside ``enabled``, ``max_depth`` and ``default_profile``,
and like all three it means *the default I hand the children I spawn* --
nothing else, wherever it is read from. A profile's OWN cap is the separate
top-level ``Profile.max_turns``, because a profile in the middle of a spawn
chain is a child and a parent at once and would otherwise have to carry both
meanings in one number. Precedence is stated once, in
:func:`_effective_max_turns`: the child's declaration beats the caller's
default beats :data:`_DEFAULT_CHILD_MAX_TURNS`. Whichever wins is VALIDATED --
a cap below 1 stops the child before it dispatches a single tool, which is a
zero-measurement child produced by a configuration nobody would think to
suspect, so it is refused rather than obeyed or quietly clamped.

This module must not import anything phenix-side: the backend-specific work of
turning a profile name into a built agent arrives as the injected
``build_child`` callable, and the finish vocabulary is defined in ``events``
(which phenix's ``finish`` module re-exports) rather than imported from phenix.
"""

import concurrent.futures
import os
import re
import sys

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.approval import ApprovalCoordinator
from qttbx.widgets.chat.agent.base import ToolSpec, close_client
from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, SubagentRecord, TokenUsage, now)
from qttbx.widgets.chat.agent.events import (
  CANCELLED, CAP, CLEAN, ERROR, TOOL_USE, TRUNCATED, AgentError, TurnDone)
from qttbx.widgets.chat.agent.session import AgentSession
from qttbx.widgets.chat.agent.tools import ToolApprovalResponse


RUN_SUBAGENT_SCHEMA = {
  "type": "object",
  "properties": {
    "task": {
      "type": "string",
      "description": "The complete instructions for the child. It sees only "
                     "this -- not your conversation.",
    },
    "profile": {
      "type": "string",
      "description": "Profile to run the child under. Omit for the default.",
    },
    "model": {
      "type": "string",
      "description": "Model for the child. Omit to inherit yours.",
    },
    "backend": {
      "type": "string",
      "description": "Backend for the child. Omit to inherit yours.",
    },
    "subject_paths": {
      "type": "array",
      "items": {"type": "string"},
      "description": "Absolute paths to the files the child is being sent to "
                     "examine (e.g. the model under review). Hashed and "
                     "recorded so a later reader can tell whether the child "
                     "examined them as they now stand.",
    },
  },
  "required": ["task"],
}


# The wording the MODEL reads before deciding to delegate, so it must say that
# the child starts blind and cannot ask. ONE tool on TWO backend surfaces --
# the registry builtin below and claude_code's SDK tool -- and a description
# that drifted between them would describe a tool that does not exist on one of
# them. Shared as a constant so it cannot; a parity test pins the pair.
RUN_SUBAGENT_DESCRIPTION = (
  "Run a child conversation to completion and return its final text. The "
  "child sees ONLY the task you write -- not your conversation -- so state "
  "everything it needs. It cannot spawn further children, and it has no "
  "user: any tool call of its own that would need approval is denied.")


class _HeadlessApprovals(ApprovalCoordinator):
  """Approval coordinator for a session nobody can answer: denies immediately.

  A child has no user and no UI -- its ``on_event`` goes to a private sink that
  only records -- so a ``ToolApprovalRequest`` surfaced from it reaches nothing
  that could ever call :meth:`ApprovalCoordinator.submit`. The base coordinator
  would hand the worker a future resolvable only by that impossible ``submit``,
  and ``AgentSession._await_approval`` parks on it with NO timeout: the shared
  ``CancelToken`` is not polled while parked, so not even the parent's Stop
  breaks the wedge. Every ``ask`` in a child would therefore be a deadlock.

  This coordinator answers each request itself, with ``deny``:

  - An ``ask`` is a request for the user's consent. With the consenting party
    absent, "no" is the only truthful answer. Auto-approving instead would turn
    a subagent into a laundering channel for precisely the calls the parent's
    policy said the user must see first.
  - A denial is an ordinary tool outcome: ``_dispatch_and_build_results`` turns
    it into a ``Denied by user policy`` error result the child reads and can
    report on, so the child's turn finishes and the parent still gets a record.

  A profile that wants its children to use tools says so explicitly, via
  ``tool_policy_default`` or per-tool ``allow`` entries -- a decision the user
  owns rather than one this module makes on their behalf.

  It sits on the coordinator rather than on the policy because
  ``AgentSession._resolve_and_approve`` rewrites ``allow`` to ``ask`` for a tool
  registered ``allow_remember=False`` AFTER consulting the policy; the
  coordinator is the single point every path into ``_await_approval`` funnels
  through, so coercing the policy alone would leave that floor still parking.

  Every denial is also RECORDED, not merely logged: a child refused every tool
  runs, reports fluently and has measured nothing, and the session log is not
  where the parent looks. ``run_child`` turns the record into an
  ``incomplete_reason`` when nothing else succeeded.

  The record is deliberately WIDER than this coordinator's own decisions.
  ``claude_code`` resolves every permission inside its own SDK callback, so a
  child on the DEFAULT backend never reaches ``open`` at all and the backstop
  fired for nobody. :meth:`record_denial` is therefore public: ``run_child``
  hands it to any agent that denies on its own behalf.

  Parameters
  ----------
  log : file-like, optional
      Stream the denials are reported on. Defaults to ``sys.stdout``.

  Attributes
  ----------
  denied_tools : list of str
      Names of every tool denied, in order, with repeats -- whether this
      coordinator decided it or a backend reported it.
  """

  def __init__(self, log=None):
    ApprovalCoordinator.__init__(self)
    self.log = log if log is not None else sys.stdout
    self.denied_tools = []

  def record_denial(self, tool_name,
                    reason="approval required, no user available"):
    """Record (and report) one denied tool, whoever decided it.

    The single place a denial becomes part of the child's RECORD rather than
    just a log line. Best-effort on the log so an agent calling this from
    inside a provider callback -- where any raise is attributed to the tool
    call rather than to the logging -- cannot be broken by a closed stream.

    Parameters
    ----------
    tool_name : str
        The tool that was refused, as the deciding layer named it.
    reason : str, optional
        Short why, for the log line.
    """
    self.denied_tools.append(tool_name)
    try:
      print("chat: subagent denied '%s' (%s)" % (tool_name, reason),
            file=self.log)
    except Exception:
      pass

  def open(self, req, emit):
    """Deny ``req`` at once, returning an already-resolved future.

    ``emit`` is still called so the base contract (register-then-surface) is
    preserved for any sink a caller does supply, and the denial is logged AND
    recorded so a child that came back empty-handed is diagnosable from the
    record, not just from the log. Nothing is registered: there is no pending
    state, hence nothing for ``cancel_turn`` to abandon.

    Parameters
    ----------
    req : ToolApprovalRequest
        The request the child's session wants answered.
    emit : callable
        The child session's ``on_event`` (a recording sink by construction).

    Returns
    -------
    concurrent.futures.Future
        Already resolved to a ``deny`` ``ToolApprovalResponse``.
    """
    emit(req)
    self.record_denial(req.tool_name)
    fut = concurrent.futures.Future()
    fut.set_result(ToolApprovalResponse(request_id=req.request_id,
                                        decision="deny"))
    return fut


class _TerminalOutcome:
  """The child's private event sink: forwards nothing, remembers the end.

  Two jobs in one object, because the child needs exactly one sink and both
  jobs need it to be this module's:

  - Nothing may reach the PARENT's handler. A child ``TurnDone`` delivered
    there unlocks the parent's UI mid-turn, after which the busy guard silently
    drops the user's next message.
  - Something must still SEE the child's terminal disposition. ``AgentSession``
    records only ``stop_reason`` onto the message (``_accumulate``), and the
    canonical ``finish`` -- the one field that distinguishes a cap, a backend
    error, a truncation and a provider block from a clean end -- exists ONLY on
    the event. A literal no-op sink discarded it, which is what made a
    ``claude_code`` child capped inside its SDK subprocess byte-identical to
    one that finished and found nothing.

  ``TurnDone`` is per-round and ``AgentSession`` loops while ``finish`` stays
  ``TOOL_USE``, so the LAST one seen is the turn's real disposition -- the same
  latch-and-overwrite rule the headless runner uses.

  A terminal event can also be ABSENT, which is a third state and not a
  degenerate ``finish``. Every backend reports a provider failure by yielding
  an ``AgentError`` and returning -- anthropic's ``APIError`` /
  ``TransportError`` arms, openai's, gemini's ``_map_error``, claude_code's
  input-extraction guard -- so no ``TurnDone`` is ever emitted, ``_accumulate``
  leaves ``stop_reason`` empty and ``AgentSession`` exits its loop. ``finish``
  alone cannot tell that from a backend that reported a ``TurnDone`` carrying
  no canonical finish, and the two mean opposite things, so both are recorded.

  Attributes
  ----------
  finish : str
      Canonical finish of the last round, ``""`` when none was reported.
  saw_turn_done : bool
      Whether the turn ENDED on a terminal event -- not whether one ever
      arrived. A latch that only ever set was wrong for the multi-round API
      backends, whose every completed round reports a
      ``TurnDone(stop_reason='tool_use')`` before the next one begins: a
      provider failure on round two therefore found the latch already set, took
      the ``_outcome_for_finish`` path instead of the no-terminal one, and read
      round one's ``TOOL_USE`` as a TURN CAP. The parent was told to raise a cap
      the child never reached, for a rate limit or an outage the message never
      mentioned. So an ``AgentError`` -- which every backend yields immediately
      before returning without a ``TurnDone`` -- clears it again.
  error : str
      Message of the last ``AgentError``, ``""`` when none was reported. The
      only account of WHY a turn with no terminal event ended.
  """

  def __init__(self):
    self.finish = ""
    self.saw_turn_done = False
    self.error = ""

  def __call__(self, event):
    if isinstance(event, TurnDone):
      self.finish = event.finish or ""
      # stop_reason is deliberately NOT latched here. AgentSession already
      # records it onto the message (_accumulate), which is the argument this
      # class's docstring uses to justify capturing `finish` -- the field that
      # exists ONLY on the event. A second copy nothing read was the leftover.
      self.saw_turn_done = True
    elif isinstance(event, AgentError):
      self.error = getattr(event, "message", "") or ""
      # The turn did not end on a terminal event, whatever an earlier ROUND
      # reported. Every backend yields this and returns, so nothing follows it
      # in the same turn -- and if a later round does report a TurnDone, the
      # branch above sets the latch again and this error is simply the last one
      # the turn survived.
      self.saw_turn_done = False


# The cap the child runs under when nobody named one: neither the child's own
# profile nor the spawning caller. Was read off the CHILD profile's
# ``subagents_default_max_turns`` (with 25 as the getattr fallback), which is
# that profile's default for ITS children -- the second meaning of a key that
# may only have one. A module constant says the same 25 without borrowing
# another key's number.
_DEFAULT_CHILD_MAX_TURNS = 25

# The session writes its cap marker as a message of exactly this shape, and as
# nothing else in the block. Matched whole rather than by prefix so a child that
# QUOTES the marker -- a review agent asked to review this very module writes it
# verbatim -- cannot hand the parent a false "it was cut off". (The prefix-only
# constant this replaced is gone: keeping a looser matcher beside the real one
# is an invitation to reach for it.)
_CAP_MARKER_RE = re.compile(r"^\[Subagent stopped at turn cap \(-?\d+\)\]$")
# The role AgentSession appends that marker under. Checked, because the marker
# is only evidence of a cap when the SESSION wrote it; the same text in an
# assistant turn is the child talking about caps, not being stopped by one.
_CAP_MARKER_ROLE = "user"
_CANCEL_MARKER = "[Subagent cancelled before it finished]"
_SDK_CAP_MARKER = "[Subagent stopped at its backend's turn cap]"
_TRUNCATED_MARKER = (
  "[Subagent's answer was cut off at the model's output-token limit]")
_ERROR_MARKER = "[Subagent ended on a backend error]"
_NO_TEXT_MARKER = "[Subagent produced no text]"
_PROVIDER_ERROR_MARKER = "[Subagent ended on a backend error: %s]"
_NO_TERMINAL_MARKER = (
  "[Subagent's turn ended with no terminal event; it did not finish]")
_NO_MEASUREMENT_MARKER = (
  "[Subagent measured nothing: it held measurement tools and not one of them "
  "returned a result, so its report is unverified]")
_CRASH_MARKER = (
  "[Subagent crashed before it finished: %s. Anything above is partial.]")

# How many TEXT-FREE assistant messages :func:`_final_text` may walk back over
# before giving up. ONE, because that is exactly how many the session can
# append after the child's last real turn: ``_collect_one_response`` commits
# each completed tool iteration and opens ONE fresh assistant message, which is
# the residual a cap or a provider error leaves behind. Walking further reaches
# rounds the child itself went on past -- and a model that calls tools without
# narrating them leaves many -- so an unbounded walk hands the parent the
# child's OPENING line as its findings.
_FINAL_TEXT_WALK_BACK = 1

# Block types that carry a TOOL'S ANSWER. Two, because a tool can run in two
# places: ``tool_result`` for one this session (or a backend's own loop)
# dispatched, ``server_tool_result`` for one the PROVIDER executed -- web
# search and its siblings, which a profile opts into through ``server_tools``
# and which ``_accumulate`` records under its own type. Reading only the first
# meant a child that measured entirely through provider-executed tools was
# reported as having measured nothing: a false accusation, and the mirror of
# the failure this backstop exists to catch. It stayed hidden while the flag
# was gated on a denial record, since such a child is denied nothing.
_RESULT_BLOCK_TYPES = ("tool_result", "server_tool_result")

# Canonical finish -> (incomplete_reason, prose marker). CLEAN and "" (unset)
# are absent on purpose: they are the only two that mean "nothing to report",
# and everything NOT in this table still fails closed below rather than
# defaulting to clean.
_FINISH_OUTCOMES = {
  CANCELLED: ("cancelled", _CANCEL_MARKER),
  # Both cap shapes land on the same reason. CAP is the backend's own cap
  # (claude_code's SDK subprocess); a FINAL TOOL_USE is an API backend's loop
  # stopped by AgentSession's cap while the model still wanted tools -- the
  # same reading the headless runner's disposition makes.
  CAP: ("turn_cap", _SDK_CAP_MARKER),
  TOOL_USE: ("turn_cap", _SDK_CAP_MARKER),
  TRUNCATED: ("truncated", _TRUNCATED_MARKER),
  ERROR: ("error", _ERROR_MARKER),
}


def _outcome_for_finish(finish):
  """Map a canonical ``TurnDone.finish`` onto ``(incomplete_reason, marker)``.

  FAILS CLOSED, exactly as the per-backend mappers that produce ``finish`` do:
  a value that is neither clean nor recognized is a provider state this code
  has never seen (``content_filter``, ``pause_turn``, ``SAFETY``, a future SDK
  subtype) and is reported under its own name rather than passed off as a
  finished run.

  Parameters
  ----------
  finish : str
      Canonical finish from the child's terminal ``TurnDone``. ``""`` means the
      backend reported none, which carries no information either way -- the
      caller's other signals decide.

  Returns
  -------
  tuple of (str, str)
      ``incomplete_reason`` and the prose marker for ``final_text``; both
      ``""`` when the child finished cleanly.
  """
  if not finish or finish == CLEAN:
    return "", ""
  known = _FINISH_OUTCOMES.get(finish)
  if known is not None:
    return known
  return finish, "[Subagent ended abnormally (%s)]" % finish


def _cap_marker(conversation):
  """Return the session's turn-cap marker, or ``""`` if the cap never fired.

  Three narrowings, each closing a different way a child's own words could be
  mistaken for the session's signal. ``AgentSession`` appends
  ``[Subagent stopped at turn cap (N)]`` when the cap fires, unconditionally,
  as the final message, in a text block containing that and nothing else -- so
  a marker is only genuine when all three hold:

  - it is on the LAST message, since that is where the session puts it;
  - that message's role is the one the session writes it under
    (:data:`_CAP_MARKER_ROLE`, ``"user"``). Without this, an ASSISTANT turn is
    read as the session's own cap signal: a child whose closing answer opens by
    quoting the marker -- a review agent asked to review this very module
    writes it verbatim -- is reported to the parent as cut off when it finished
    fine, and its whole answer is then relabelled a truncated one;
  - the text is the marker WHOLE (:data:`_CAP_MARKER_RE`), not merely prefixed
    by it, so a user-role message that opens with the phrase and continues into
    other prose -- a task quoting it, a tool result echoing it -- is not read as
    a cap either.

  Should the session ever change the marker's text, role or position, the cap
  exercises fail loudly; a test that breaks beats a live false positive nobody
  can see.

  This covers the API backends only. ``claude_code``'s loop is capped inside
  its SDK subprocess, which reports through ``TurnDone(finish=CAP)`` and writes
  no marker at all; that half is read from the event sink.

  Parameters
  ----------
  conversation : Conversation
      The child's transcript.

  Returns
  -------
  str
      The marker text, or ``""`` when the child was not capped.
  """
  for message in conversation.messages[-1:]:
    if message.role != _CAP_MARKER_ROLE:
      continue
    for block in message.content:
      if block.type != "text":
        continue
      candidate = (block.data.get("text", "") or "").strip()
      if _CAP_MARKER_RE.match(candidate):
        return candidate
  return ""


def _final_text(conversation):
  """The child's last assistant turn that SAID something, plus any cap marker.

  "Last assistant turn" is not "last assistant message", and the difference is
  the whole point of this function. A child's transcript routinely ends in an
  assistant message carrying no text at all, on the two paths that matter most:

  - the turn cap. ``AgentSession`` commits each completed tool iteration and
    starts a fresh assistant message; when the cap fires on the next round that
    residual is appended with only its ``tool_use`` blocks, or with nothing at
    all, and THEN the cap marker goes on after it.
  - a backend error. The same residual is appended, empty, when the provider
    fails after an iteration was already committed.

  Stopping at the last assistant message therefore read the residual, found no
  text, and returned "" -- so a child that measured a structure, found real
  problems and then ran out of turns reported *nothing but the cap marker* to
  the parent, which is worse than the truncation it was reporting. Scanning
  back to the last assistant message that actually carries text returns the
  findings the child did produce; they are still the child's own words, and the
  marker appended below is what says they are partial.

  Only assistant messages are considered, in reverse, and the FIRST one with
  text wins -- never a concatenation of several, which would graft an earlier
  round's superseded reasoning onto the final answer.

  And the walk is BOUNDED at :data:`_FINAL_TEXT_WALK_BACK`. Unbounded it
  reaches back through however many silent tool rounds the child ran -- a model
  that calls tools without narrating them writes no text for round after round
  -- and returns its OPENING line as its report: "Let me start by looking at
  the model", handed to the parent as the findings. That is worse than the
  empty string this walk was added to avoid, because it reads like an answer.
  One step is the whole of what the residual costs and the whole of what is
  safe: ``_collect_one_response`` commits each completed iteration and starts
  ONE fresh assistant message, so at most one text-free residual can be
  appended after the last real turn. Anything further back is a round the child
  itself went on past.

  The cap marker is the difference between "the child finished" and "the child
  was cut off", so it must reach the parent rather than being trimmed away with
  the rest of the transcript. It is appended as a USER message, so the
  last-assistant scan alone never sees it.

  Parameters
  ----------
  conversation : Conversation
      The child's transcript.

  Returns
  -------
  str
      The child's closing text, with the cap marker appended when one fired.
  """
  text = ""
  examined = 0
  for message in reversed(conversation.messages):
    if message.role != "assistant":
      continue
    parts = [b.data.get("text", "") for b in message.content
             if b.type == "text"]
    candidate = "\n".join(p for p in parts if p)
    if candidate.strip():
      text = candidate
      break
    examined += 1
    if examined > _FINAL_TEXT_WALK_BACK:
      # Past the residual, so anything still earlier is a round this child ran
      # on past. It said nothing after that; "" is the honest answer, and the
      # caller's _NO_TEXT_MARKER is what says so.
      break
  marker = _cap_marker(conversation)
  if marker and marker not in text:
    text = (text + "\n" + marker).strip()
  return text


def _tool_names_by_use_id(conversation):
  """``{tool_use_id: tool_name}`` for every tool call in the transcript.

  A ``tool_result`` block carries only the id it answers, so the name has to
  come from the ``tool_use`` block that raised it. Both backend families write
  that block: the API loop from its own ``ToolUseRequested``, ``claude_code``
  from the tool_use it observed in its SDK subprocess.

  ``server_tool_use`` counts too. A provider-executed tool (web search, code
  execution -- whatever a profile's ``server_tools`` opts into) runs at the
  provider and is recorded by ``_accumulate`` under its own block type, so a
  scan for ``tool_use`` alone knew none of their names and
  :func:`_any_tool_succeeded` could not credit their results.
  """
  names = {}
  for message in conversation.messages:
    for block in message.content:
      if block.type in ("tool_use", "server_tool_use"):
        use_id = block.data.get("id")
        if use_id:
          names[use_id] = block.data.get("name", "")
  return names


def _is_skill_tool(name, tools):
  """Whether *name* is one of the child's skill-introspection tools.

  The structural signal, not a name list: ``ToolRegistry.source_of`` reports
  ``'skill'`` for exactly the pairs ``SkillLoader.tools`` produced, which is
  the same set the phenix-side factory pre-allows on a child's policy. A name
  list would have to be kept in step with that factory across two repos, and a
  renamed or newly added skill tool would silently start counting as a
  measurement again.

  Registered under the BARE name on both surfaces, so an SDK-prefixed
  ``mcp__<server>__<bare>`` is retried bare -- ``claude_code`` wraps the skill
  tools into an in-process MCP server and reports them that way. A name the
  registry does not know is NOT a skill tool: an MCP tool a claude_code child
  ran inside its own subprocess never enters the registry at all, and calling
  that "no measurement" would fire the flag on a child that measured plenty.

  Parameters
  ----------
  name : str
      Tool name as the transcript recorded it.
  tools : ToolRegistry or None
      The child's registry (``bundle.tools``). ``None`` -- or anything
      without ``source_of`` -- classifies nothing as a skill tool.

  Returns
  -------
  bool
  """
  source_of = getattr(tools, "source_of", None)
  if not callable(source_of) or not name:
    return False
  if source_of(name) == "skill":
    return True
  if name.startswith("mcp__"):
    parts = name.split("__", 2)                # ['mcp', server, bare]
    if len(parts) > 2 and source_of(parts[2]) == "skill":
      return True
  return False


def _tool_server(name, tools):
  """The MCP server a tool belongs to, or ``""`` when it belongs to none.

  Two surfaces, one question. ``claude_code`` names every MCP tool
  ``mcp__<server>__<bare>``, so the server is IN the name -- which is the only
  route available for it, because a tool its SDK subprocess ran never enters
  the registry. The API backends register theirs through
  ``ToolRegistry.register_mcp_tool``, which stamps the source ``mcp:<server>``
  (and renames a colliding tool ``<server>:<bare>``, leaving the source
  untouched), so the registry answers there.

  A built-in, a skill tool and a name nothing knows all return ``""``: they
  belong to no server, which is a different answer from "a server nobody
  configured".

  Parameters
  ----------
  name : str
      Tool name as the transcript recorded it.
  tools : ToolRegistry or None
      The child's registry.

  Returns
  -------
  str
  """
  if not name:
    return ""
  if name.startswith("mcp__"):
    parts = name.split("__", 2)                  # ['mcp', server, bare]
    if len(parts) > 2 and parts[1]:
      return parts[1]
  source_of = getattr(tools, "source_of", None)
  if callable(source_of):
    source = source_of(name) or ""
    if source.startswith("mcp:"):
      return source[len("mcp:"):]
  return ""


def _is_measurement_tool(name, tools, servers):
  """Whether a result from *name* says anything about the thing under review.

  ONE predicate, because the four ways a child ends up measuring nothing --
  denied, errored, never given the tool, or answered only by bookkeeping --
  differ in how they fail and not in what they leave behind. Both exclusions
  are STRUCTURAL, read off the same wiring the child was built from, so
  neither can drift out of step with a renamed or newly added tool:

  - **Skill-introspection tools** (:func:`_is_skill_tool`): registry source
    ``'skill'``, the same signal the phenix-side factory uses to choose what to
    pre-allow. They read the child's own prompt material.
  - **Tools from a server the child's PROFILE did not ask for.** *servers*
    names the MCP servers this build really wired up. A backend adds servers of
    its own on top of those -- ``claude_code`` builds an in-process
    ``phenix_chat`` server on EVERY agent, carrying the job-history and
    ask-user tools, and auto-approves the job-history one ahead of any policy
    check. With no provider wired (nothing wires one for a child)
    ``phenix_get_job_history`` returns ``"(no jobs recorded)"`` as a
    SUCCESSFUL result, which satisfied "did any tool succeed" for a child that
    had measured nothing at all. Asking whose server it came from excludes
    every such housekeeping tool at once, present and future, without naming
    any of them.

  Everything else counts, and that direction matters as much: a provider's own
  file/search built-ins belong to no MCP server and are not in the registry, so
  they stay measurements -- a child that read the files it was asked about
  measured plenty.

  Parameters
  ----------
  name : str
      Tool name as the transcript recorded it.
  tools : ToolRegistry
      REQUIRED. The child's registry (``bundle.tools``).
  servers : set of str
      Names of the MCP servers the child's profile really gave it
      (``bundle.measurement_servers``). REQUIRED. There is deliberately no
      "unknown" value: the two things this could do without an inventory are
      opposite -- count every tool as a measurement and never flag, or count
      none and flag every child -- so guessing is worse than refusing.
      ``_measurement_servers`` raises rather than inventing one.

  Returns
  -------
  bool
  """
  if not name:
    return False
  if _is_skill_tool(name, tools):
    return False
  server = _tool_server(name, tools)
  return not server or server in servers


def _result_is_error(block):
  """Whether a result block records a FAILURE rather than an answer.

  The two result types say so differently. A ``tool_result`` carries the
  boolean ``is_error`` this session (or a backend's own loop) set. A
  ``server_tool_result`` carries the provider's raw payload, opaque at this
  layer -- so the two shapes a provider actually uses to signal failure are
  recognized: a ``type`` ending ``_error`` (Anthropic's
  ``web_search_tool_result_error``) and an ``error_code`` key.

  Anything else in that opaque dict reads as an ANSWER. That direction is
  deliberate: an unrecognized payload treated as a failure would accuse a child
  that measured, which is the fault this predicate was widened to fix, and the
  accusation is the one a reviewer cannot argue with.

  Parameters
  ----------
  block : ContentBlock
      A block whose ``type`` is in :data:`_RESULT_BLOCK_TYPES`.

  Returns
  -------
  bool
  """
  data = block.data or {}
  if block.type == "tool_result":
    return bool(data.get("is_error"))
  content = data.get("content")
  if not isinstance(content, dict):
    return False
  if content.get("error_code"):
    return True
  return str(content.get("type", "") or "").endswith("_error")


def _any_tool_succeeded(conversation, tools, servers):
  """Whether the child got at least one non-error result from a MEASUREMENT
  tool.

  The question a denial record cannot answer on its own: a child denied one
  destructive tool while reading five files measured plenty, and flagging it
  would make the flag meaningless. Only a child that got nothing back from
  anything it could measure with is the "reports fluently, measured nothing"
  case.

  What counts is :func:`_is_measurement_tool`'s question, so a bookkeeping
  tool answering successfully is not a measurement however cleanly it
  answered.

  Parameters
  ----------
  conversation : Conversation
      The child's transcript.
  tools : ToolRegistry
      The child's registry, for classifying skill tools. REQUIRED, for the
      same reason as *servers* below.
  servers : set of str
      REQUIRED. The child's measurement servers; see
      :func:`_is_measurement_tool`.

  Returns
  -------
  bool
      ``True`` when any measurement ``tool_result`` block is not an error.
  """
  names = _tool_names_by_use_id(conversation)
  for message in conversation.messages:
    for block in message.content:
      if block.type not in _RESULT_BLOCK_TYPES or _result_is_error(block):
        continue
      if _is_measurement_tool(names.get(block.data.get("tool_use_id"), ""),
                              tools, servers):
        return True
  return False


def _measurement_servers(bundle):
  """The MCP servers *bundle* says this child can really measure with.

  Parameters
  ----------
  bundle : libtbx.group_args
      The child bundle.

  Returns
  -------
  set of str
      The child's measurement servers. Raises rather than returning a
      third "unknown" value; see :func:`_is_measurement_tool`.
  """
  names = getattr(bundle, "measurement_servers", None)
  if names is None:
    from libtbx.utils import Sorry
    raise Sorry(
      "child bundle declares no measurement_servers. The inventory is "
      "REQUIRED: without it the zero-measurement gate has to guess, and its "
      "two possible guesses are opposite -- treat every tool as a "
      "measurement and never flag, or treat none as one and flag every "
      "child. Build the bundle through the factory, which always sets it.")
  return set(names)


def _child_was_equipped(bundle):
  """Whether this child was built holding anything it could measure with.

  The gate on the zero-measurement flag, and the reason the flag can be
  unconditional rather than waiting for a denial: a child that HAD measurement
  tools and came back with no result from any of them measured nothing,
  whether it was denied, errored, or simply never called one. A child that was
  never equipped is a different fault with a different owner -- the phenix-side
  factory refuses to spawn one at all (``_assert_measurement_tools``) -- and a
  caller that hands this module no inventory has told it nothing, so it says
  nothing.

  Both routes are read, because the two backend families deliver tools
  differently and each is invisible on the other's route: ``claude_code``'s MCP
  tools live in its SDK subprocess and never enter the registry (only the
  server list on the bundle records them), while an API child's arrive as
  registry entries sourced ``mcp:<server>``.

  Parameters
  ----------
  bundle : libtbx.group_args
      The child bundle.

  Returns
  -------
  bool
  """
  if getattr(bundle, "measurement_servers", None):
    return True
  source_map = getattr(getattr(bundle, "tools", None), "tool_to_source_map",
                       None)
  if not callable(source_map):
    return False
  return any(str(source or "").startswith("mcp:")
             for source in source_map().values())


def _stamp_provenance(conversation, bundle):
  """Stamp the child's assistant messages with WHO actually produced them.

  ``AgentSession`` stamps each assistant message from what it can see --
  ``agent.model`` and ``profile.backend`` -- and for a child neither is
  authoritative. The factory resolves the pair before building the agent, and
  substitutes a different backend (with the parent's model) when the requested
  one has no usable credentials; the child's ``Profile`` still names the
  backend that was ASKED for. So a persisted child transcript said it came from
  the backend somebody wanted rather than the one that answered -- and a
  cross-model review whose transcript cannot say who reviewed is not a
  cross-model review.

  ``bundle.model`` / ``bundle.backend`` are the resolved pair, the same values
  ``SubagentRecord.model`` and the child ``ConversationMeta`` already carry, so
  this makes the per-message stamp agree with the record around it rather than
  introducing a third answer. A bundle that leaves either blank stamps nothing
  for it and the session's own guess stands, which is the only honest fallback.

  Parameters
  ----------
  conversation : Conversation
      The child's transcript, mutated in place.
  bundle : libtbx.group_args
      The child bundle, supplying the resolved ``model`` / ``backend``.
  """
  model = getattr(bundle, "model", None)
  backend = getattr(bundle, "backend", None)
  for message in conversation.messages:
    if message.role != "assistant":
      continue
    if model:
      message.model = model
    if backend:
      message.backend = backend


def _peak_usage(conversation):
  """Roll the child's usage up into one :class:`TokenUsage`.

  The four cost fields sum; ``context_tokens`` is a peak and takes the maximum.
  Mixing those two rules is the bug this function exists to prevent.
  """
  total = TokenUsage()
  for message in conversation.messages:
    usage = getattr(message, "usage", None)
    if usage is None:
      continue
    total.input += usage.input
    total.output += usage.output
    total.cache_read += usage.cache_read
    total.cache_creation += usage.cache_creation
    total.context_tokens = max(total.context_tokens, usage.context_tokens)
  return total


def _release_child(bundle, log):
  """Release everything ``build_child`` opened for one child.

  A child that FAILS releases through ``build_child``'s own guard, which covers
  only a failed build; a child that SUCCEEDS held its resources forever. Per
  spawn that is one provider client or ``claude`` CLI subprocess (plus the SDK
  subprocess's own MCP servers and a daemon asyncio loop thread) and one
  ``phenix.mcp_server`` subprocess per session-side server. The declared
  consumer is a review agent invoked repeatedly in a long-lived session, so
  "the process will exit eventually" is not a bound -- and
  ``McpServerConnection`` has no finalizer, so gc does not do it either.

  Deterministic, not best-effort-on-teardown: this is the same pair the
  parent's own ``closeEvent`` performs, done as soon as the child can no longer
  use them. Both halves swallow their errors -- a release that raised would
  convert a completed child into a failed tool call, losing its text.

  Parameters
  ----------
  bundle : libtbx.group_args
      The child bundle. ``agent`` is closed via
      :func:`qttbx.widgets.chat.agent.base.close_client`; ``connections`` (when
      the factory started any) are stopped. A bundle carrying neither is fine.
  log : file-like
      Session log for a failed stop.
  """
  close_client(getattr(bundle, "agent", None))
  for conn in (getattr(bundle, "connections", None) or []):
    try:
      conn.stop()
    except Exception as exc:
      print("chat: subagent MCP server stop failed: %s" % exc, file=log)


def _declared_child_max_turns(bundle):
  """The cap the CHILD's own profile declared for itself, or ``None``.

  ONE meaning, read from keys that carry only that meaning. The child's own cap
  is ``Profile.max_turns`` -- a top-level key whose whole purpose is "the cap
  that governs when this profile is run" -- delivered either straight off the
  child's ``Profile`` or, equivalently, pre-extracted by the factory onto
  ``bundle.max_turns``. The bundle field wins when both are set, because a
  factory that resolved the number saw the child's own JSON and this module did
  not.

  ``subagents.default_max_turns`` is deliberately NOT consulted here. That key
  lives in the ``subagents`` block alongside ``enabled``, ``max_depth`` and
  ``default_profile``, every one of which is about the children a profile
  SPAWNS; reading it off the child as "my own cap" gave one key two opposite
  meanings depending on who held it, and a profile in the middle of a spawn
  chain would have had to carry both at once.

  Parameters
  ----------
  bundle : libtbx.group_args
      The child bundle. ``max_turns`` (optional) is the declared cap the
      factory read from the child's own JSON; ``profile.max_turns`` is the same
      declaration read off the loaded ``Profile``.

  Returns
  -------
  int or None
      The declared cap, or ``None`` when the child's author wrote none -- which
      is what leaves the spawning caller's default in force.
  """
  declared = getattr(bundle, "max_turns", None)
  if declared is not None:
    return declared
  return getattr(getattr(bundle, "profile", None), "max_turns", None)


def _validate_max_turns(value, source):
  """Return *value* as a usable turn cap, or raise ``Sorry``.

  REFUSES rather than clamps, and refuses before the child runs. A cap of ``0``
  or a negative reaches ``AgentSession.run_turn``'s ``iterations >= max_turns``
  test on the FIRST pass, so the child is cut off before a single tool is
  dispatched: it measures nothing, and the only trace is a ``turn_cap`` record
  nobody would think to blame on a configuration file. Clamping to 1 instead
  would run the child under a number its author never wrote -- and one turn is
  itself a zero-measurement child on any real task -- so the misconfiguration
  would stay hidden while still costing a provider call. Refusal puts the fault
  in front of the only person who can fix it, the profile's author, and costs
  nothing scientifically: a child capped at zero turns was never going to
  produce a usable answer.

  A refusal here is also an outcome this system already has a vocabulary for --
  the same "the spawn was refused before anything ran" shape the measurement
  gate uses -- rather than a new kind of half-run child.

  Parameters
  ----------
  value : int
      The cap selected by :func:`_effective_max_turns`.
  source : str
      Where it came from, named in the message so the author knows which key to
      edit.

  Returns
  -------
  int
      *value*, when it is a usable cap.

  Raises
  ------
  libtbx.utils.Sorry
      When *value* is not an integer, or is below 1.
  """
  # bool is an int subclass; True would silently mean a 1-turn child.
  if isinstance(value, bool) or not isinstance(value, int):
    raise Sorry("subagent turn cap (%s) must be an integer, not %r"
                % (source, value))
  if value < 1:
    raise Sorry("subagent turn cap (%s) is %d; it must be at least 1. A cap "
                "below 1 stops the child before it can dispatch a single "
                "tool, so it would report on nothing it measured."
                % (source, value))
  return value


def _effective_max_turns(bundle, max_turns):
  """The cap the child really runs under, validated.

  The one place the precedence is stated, so no call site can disagree with
  another:

  1. the cap the CHILD's own profile declared (``Profile.max_turns``, via
     :func:`_declared_child_max_turns`) -- a specific decision about this
     child, made by the author who knows what the job needs;
  2. else the caller's number, which for both tool paths is the PARENT
     profile's ``subagents.default_max_turns`` -- a *default* for children, as
     its name says, so anything the child declared outranks it;
  3. else :data:`_DEFAULT_CHILD_MAX_TURNS`, for a direct caller that passed
     nothing and a bundle that declares nothing.

  Applied here rather than at the call sites because both of them pass the
  parent's number and neither holds the child's declaration. A reviewer profile
  asking for 40 turns otherwise ran at the parent's 25 on every backend but
  ``claude_code`` (whose separate SDK cap the factory writes directly), which
  truncates the review and blocks the readiness claim it was spawned to settle.

  Whatever wins is then validated (:func:`_validate_max_turns`): an
  out-of-range cap is refused, not silently obeyed.

  Parameters
  ----------
  bundle : libtbx.group_args
      The child bundle.
  max_turns : int or None
      The caller's default cap for the children it spawns.

  Returns
  -------
  int
      The cap to hand ``AgentSession.run_turn``.

  Raises
  ------
  libtbx.utils.Sorry
      When the winning cap is not a usable one.
  """
  declared = _declared_child_max_turns(bundle)
  if declared is not None:
    return _validate_max_turns(declared, "the child profile's max_turns")
  if max_turns is not None:
    return _validate_max_turns(
      max_turns, "the caller's subagents.default_max_turns")
  return _DEFAULT_CHILD_MAX_TURNS


def run_child(parent_session, task, bundle, cancel, tool_use_id,
              max_turns=None, subject_paths=None):
  """Run one child conversation to completion and persist its record.

  Parameters
  ----------
  parent_session : AgentSession
      The spawning session; supplies depth, storage, the parent id and the
      usage rollup.
  task : str
      The child's entire instructions.
  bundle : libtbx.group_args
      ``agent`` / ``tools`` / ``policy`` / ``profile`` / ``profile_name`` /
      ``model`` / ``backend``, as returned by the injected ``build_child``,
      plus four optional fields it may add: ``connections`` (session-side MCP
      servers to release afterwards), ``notice`` (prose the parent must see
      about how the child was actually built -- a substituted backend/model,
      say, which changes WHO wrote the answer), ``max_turns`` (the cap the
      child's own profile declared for ITSELF -- ``Profile.max_turns``, NOT
      its ``subagents.default_max_turns``, which is that profile's default for
      the children IT spawns -- outranking *max_turns* below; see
      :func:`_effective_max_turns`) and ``measurement_servers`` (the names
      of the MCP servers this build really wired up). The last two of those are
      what tells a MEASUREMENT tool from the pre-allowed skill tools and from
      the backend's own housekeeping servers (:func:`_is_measurement_tool`),
      and whether the child was equipped to measure at all
      (:func:`_child_was_equipped`); a bundle carrying neither is refused by
      :func:`_measurement_servers`, because the gate cannot guess an inventory
      without choosing between never flagging and always flagging. ``agent`` is offered the
      denial sink when it has a ``set_denial_sink`` hook, which is how a
      backend that decides permissions itself feeds the same record.
  cancel : CancelToken
      Shared with the parent's Stop; the child polls it between events.
  tool_use_id : str
      The parent tool call this child answers; recorded for the GUI.
  max_turns : int, optional
      The caller's DEFAULT child turn cap for ``AgentSession.run_turn``, which
      governs the API backends; a cap the child's own profile declared outranks
      it (:func:`_effective_max_turns`). Both tool paths pass the PARENT
      profile's ``subagents.default_max_turns``, so ``None`` arrives only when
      the parent profile omits the field -- a plain ``group_args`` parent, or a
      ``Profile`` predating it -- and then :data:`_DEFAULT_CHILD_MAX_TURNS`
      applies. A ``claude_code`` child is capped elsewhere entirely: its agentic
      loop runs in the SDK subprocess, reached only through
      ``extra_sdk_options['max_turns']``.

  Returns
  -------
  SubagentRecord
      Persisted under the parent conversation. ``incomplete_reason`` says
      whether the child actually finished, and how it did not; a child that was
      cut off ANY way also carries a marker in ``final_text``, since that
      string is all the parent model sees.

  Raises
  ------
  libtbx.utils.Sorry
      When the effective turn cap is out of range (:func:`_validate_max_turns`).
  Exception
      Whatever a crashed child raised, re-raised AFTER its record has been
      built, rolled up and persisted -- a crash is not a finish, so it must not
      be reported as one, but the transcript of the work it did before crashing
      is exactly what the parent needs and nothing else saves it.
  """
  started_at = now()
  # Hashed BEFORE the child runs, so the digest describes the
  # file the child was pointed at rather than whatever it looks
  # like once the child (or anything else) has finished.
  subject_digests = digest_subject_paths(subject_paths)
  child_conv = Conversation.new(profile_name=bundle.profile_name,
                                model=bundle.model,
                                backend=bundle.backend)
  outcome = _TerminalOutcome()
  approvals = _HeadlessApprovals(log=parent_session.log)
  # A backend that resolves permissions in its OWN layer (claude_code, inside
  # the SDK callback) never reaches the coordinator above, so its denials were
  # invisible to the record and the denied-every-tool backstop fired for nobody
  # on the DEFAULT backend. Duck-typed: this module must not import, or know
  # about, any particular backend, and an agent without the hook is unaffected.
  sink_setter = getattr(bundle.agent, "set_denial_sink", None)
  if callable(sink_setter):
    sink_setter(approvals.record_denial)
  child = AgentSession(
    agent=bundle.agent, conversation=child_conv,
    storage=parent_session.storage, tools=bundle.tools,
    policy=bundle.policy, profile=bundle.profile,
    depth=parent_session.depth + 1,
    on_event=outcome,                       # never the parent's sink
    approvals=approvals,
    log=parent_session.log,
    # The child conversation is NEVER saved -- its transcript lives inside the
    # SubagentRecord, under the parent. An image filed under child_conv's id
    # therefore lands in a conversations/<child_id>/attachments/ directory for
    # a conversation that does not exist, and the record's image blocks, which
    # the GUI resolves against the conversation the record is stored under,
    # point at nothing. Filing under the parent's id puts the bytes where the
    # record is.
    attachment_conv_id=parent_session.conv.meta.id)
  # The session mints its own sub_id at depth>0, so the record carries THAT id
  # rather than a second one: two ids for one child would leave the stored
  # record and the session that produced it with no name in common.
  sub_id = child.sub_id
  user_message = Message(
    role="user", content=[ContentBlock(type="text", data={"text": task})],
    timestamp=started_at)
  from qttbx.widgets.chat.agent.errors import TurnCancelled
  last_message = None
  cancelled = False
  crash = None
  try:
    # Inside the try, so a refused cap still releases what build_child opened.
    max_turns = _effective_max_turns(bundle, max_turns)
    try:
      last_message = child.run_turn(user_message, cancel, max_turns=max_turns)
    except TurnCancelled:
      # A cancelled child is a partial result, not an error: the parent asked
      # for it and needs whatever was produced, plus the record for the GUI.
      # An escape here would cost both -- the record is built BELOW, so the
      # transcript would never be persisted and the GUI would have no handle on
      # a child that did real work. The parent's own dispatch would turn the
      # raise into a bare "Cancelled during execution" tool_result.
      cancelled = True
    except Exception as exc:
      # A CRASH, which is neither a cancel nor a finish. Held, not swallowed:
      # it is re-raised at the very end, after the record is built, rolled up
      # and stored. Letting it propagate from here instead lost the entire
      # transcript of a child that had run twenty turns and measured plenty --
      # depth>0 disables autosave and checkpointing, so this module's own store
      # call is the ONLY thing that ever writes a child's messages to disk, and
      # it is BELOW. The parent still gets its error tool_result (the re-raise),
      # the GUI still gets the transcript, and `crashed` on the record keeps the
      # two from being confused: recording a crash as a clean, empty finish is
      # the exact confusion this module exists to prevent.
      # `Exception`, never `BaseException`: a KeyboardInterrupt or a
      # SystemExit must unwind immediately, not stop to write a record.
      crash = exc
    # The other cancel shape: the shared token trips mid-run and the session
    # short-circuits, stamping the returned message 'cancelled' instead of
    # raising.
    if getattr(last_message, "stop_reason", None) == "cancelled":
      cancelled = True
    # The resolved pair the FACTORY chose, stamped on every assistant message
    # the child produced. The session stamps its own guess -- agent.model and
    # profile.backend -- which is not the same thing: the factory substitutes a
    # backend (and its model) when the requested one has no usable credentials,
    # and the profile still names the one that was asked for. A persisted child
    # transcript that cannot say which backend really wrote it defeats the
    # point of commissioning a cross-model review, so the bundle's answer -- the
    # same one `record.model` already carries -- wins wherever it has one.
    _stamp_provenance(child_conv, bundle)
    final_text = _final_text(child_conv)
    # Every non-clean end lands on BOTH fields. incomplete_reason is what code
    # reads; the marker is what the parent MODEL reads, since the tool result
    # string is the only thing about the child it ever sees. Without the prose
    # half a child cut off mid-review is byte-identical to one that reviewed
    # everything and found nothing.
    incomplete_reason, marker = _outcome_for_finish(outcome.finish)
    if cancelled:
      # Wins over the sink: a session-synthesized cancel is reported through
      # run_turn's return, and a TurnCancelled unwind may leave the sink
      # holding an earlier round's tool_use finish.
      incomplete_reason, marker = "cancelled", _CANCEL_MARKER
    elif crash is not None:
      # Ahead of every other non-clean reading below, for the same reason
      # cancelled is: the unwind left the sink holding whatever the LAST
      # completed round said, which for a child that crashed after real work is
      # a perfectly clean-looking finish. What actually ended the child is the
      # exception, and that is what the record must say.
      incomplete_reason = "crashed"
      marker = _CRASH_MARKER % ("%s: %s" % (type(crash).__name__, crash))
    elif _cap_marker(child_conv):
      # AgentSession's own cap. The transcript marker is already in final_text
      # (see _final_text), so no second one is added.
      incomplete_reason, marker = "turn_cap", ""
    elif not outcome.saw_turn_done:
      # No terminal event AT ALL. Every backend reports a provider failure by
      # yielding an AgentError and returning, so there is no finish to look up
      # and _outcome_for_finish's fail-closed table never sees it -- the turn
      # simply stopped, and an untouched incomplete_reason reads as "it ran to
      # its own end". The two shapes are told apart because they diagnose
      # differently: an error has a provider message the parent can act on, a
      # silent stop has nothing to say beyond its own name.
      if outcome.error:
        incomplete_reason = "error"
        marker = _PROVIDER_ERROR_MARKER % outcome.error
      else:
        incomplete_reason = "no_terminal_event"
        marker = _NO_TERMINAL_MARKER
    # The zero-measurement reading is computed for EVERY non-cancelled,
    # non-crashed child, whatever ended it -- it is a fact about what the child
    # got back, not an alternative to the terminal cause. Ordered as one more
    # `elif` arm it was simply lost whenever another arm matched first: a child
    # refused every tool AND cut off by a provider failure reported only the
    # failure, so the parent was told to retry a spawn whose real fault was a
    # policy that granted it nothing -- and retrying reproduces it exactly.
    zero_reason, zero_marker = "", ""
    if not cancelled and crash is None and not _any_tool_succeeded(
        child_conv, bundle.tools, _measurement_servers(bundle)):
      # NOTHING it could measure with came back: it wrote a plausible report
      # about nothing. Gated on the OUTCOME, not on how the outcome came about
      # -- the previous gate was ``approvals.denied_tools``, which made the
      # whole backstop unreachable for a child whose calls came back as ERRORS
      # (an unknown tool, a server that never started, a backing program that
      # failed) and for one that never called a tool at all. Both report as
      # fluently as a denied one, and neither populates a denial record.
      #
      # Two triggers, one predicate, because the two say different things:
      if approvals.denied_tools:
        # Refused every tool it tried. Reported by name: which tools a child
        # was refused is what the user acts on.
        zero_reason = "tools_denied"
        zero_marker = (
          "[Subagent was denied every tool it tried (%s); it measured "
          "nothing]" % ", ".join(sorted(set(approvals.denied_tools))))
      elif _child_was_equipped(bundle):
        # No denial, so nothing was refused -- this child HELD measurement
        # tools and still came back with no result from any of them.
        zero_reason = "no_measurement"
        zero_marker = _NO_MEASUREMENT_MARKER
    if zero_reason and not incomplete_reason:
      # Nothing else ended it, so this IS what ended it, and it becomes the
      # reason code as well as the prose.
      incomplete_reason = zero_reason
    if marker and marker not in final_text:
      final_text = (final_text + "\n" + marker).strip()
    if zero_marker and zero_marker not in final_text:
      # Appended even when another cause owns `incomplete_reason`. The two
      # markers answer different questions -- why the child stopped, and
      # whether anything it says was measured -- and the parent MODEL, which
      # sees only this string, needs both to decide whether to retry or to fix
      # the configuration.
      final_text = (final_text + "\n" + zero_marker).strip()
    if not final_text.strip():
      # An empty tool result is a provider 400 on the parent's next request,
      # and reads as "ran fine, nothing to report" if it survives -- which is
      # what a claude_code finish=ERROR produced, its assistant turn carrying
      # no text at all. Checked BEFORE the notice below, or a notice would
      # stand in for an answer that was never written.
      final_text = _NO_TEXT_MARKER
    notice = getattr(bundle, "notice", "") or ""
    if notice:
      # FIRST, and always -- even on a clean run. The parent asked for a second
      # opinion from a named model; if it got one from somewhere else, that is
      # the first thing it needs to know about the answer below.
      final_text = (notice + "\n" + final_text).strip()
    usage = _peak_usage(child_conv)
    record = SubagentRecord(
      sub_id=sub_id,
      parent_conversation_id=parent_session.conv.meta.id,
      parent_tool_use_id=tool_use_id,
      task=task,
      profile_name=bundle.profile_name,
      model=bundle.model,
      started_at=started_at,
      finished_at=now(),
      final_text=final_text,
      token_usage=usage,
      messages=list(child_conv.messages),
      incomplete_reason=incomplete_reason,
      subject_digests=subject_digests)
    # The child's spend rides on the record (token_usage above) and nowhere
    # else. It used to be filed a second time on the parent session, which
    # nothing read: two sinks for one number, the redundant one kept alive
    # only by the tests that pinned it.
    if parent_session.storage is not None:
      try:
        parent_session.storage.store_subagent(
          parent_session.conv.meta.id, record)
      except Exception as exc:
        # Every sibling save in this subsystem is best-effort for the same
        # reason: a full disk or a read-only project must not convert a
        # completed child into a bare tool error that also throws away its
        # text. The record still returns; only the GUI's later handle on it is
        # lost, and that is the smaller loss.
        print("chat: subagent record not saved: %s" % exc,
              file=parent_session.log)
    if crash is not None:
      # Everything above has run: the record is built, the usage is rolled up
      # and the transcript is on disk. NOW the crash resumes its unwind, so the
      # parent's dispatch turns it into the error tool_result it always was. A
      # crashed child must never come back as a returned record -- that would
      # make it read as a child that finished with nothing to say.
      raise crash
    # Offer the bundle a last word on what the PARENT reads. Duck-typed for the
    # same reason as the denial sink above: this module must not know what any
    # particular caller wants appended. Phenix uses it to hand a review round
    # the loop state derived from the recorded digests, so the exit condition
    # arrives next to the decision instead of thousands of messages behind it.
    #
    # Applied AFTER store_subagent and only to the returned copy: the persisted
    # report stays the reviewer's own words, which is what later readers cite.
    # Best-effort -- an annotator that raises must not turn a completed child
    # into a tool error and throw away its report.
    annotate = getattr(bundle, "annotate_result", None)
    if callable(annotate):
      extra = None
      try:
        extra = annotate(record)
      except Exception as exc:
        print("chat: subagent result annotator failed: %s" % exc,
              file=parent_session.log)
      if extra:
        record.final_text = ("%s\n\n%s" % (record.final_text, extra)).strip()
    return record
  finally:
    _release_child(bundle, parent_session.log)


#: Longest derived title kept intact; longer ones are elided. The sidebar
#: elides for display anyway, but the stored title is also what search matches
#: and what an export header carries, so it is bounded here too.
_TITLE_MAX = 80


def subagent_conversation_title(record):
  """Derive a human title for a child's transcript from the child's own report.

  Deterministic, and free: no model call, so it works headless, offline, and
  after the fact on records written before this function existed. The source
  is the report's first markdown heading, which the review harness already
  requires reports to lead with.

  Preamble before the heading is skipped -- a child commonly narrates a line or
  two ("I now have coverage; compiling the report") before the document starts
  -- and fenced blocks are skipped so a ``#`` comment inside a shell sample
  cannot be mistaken for the title.

  Parameters
  ----------
  record : SubagentRecord
      The finished child. Only ``final_text``, ``profile_name`` and the
      timestamps are read.

  Returns
  -------
  str
      The heading text, elided to a bounded length; or, when the report has no
      heading (it failed, was capped, or simply answered in prose), a stable
      fallback naming the profile and when the child ran.
  """
  fenced = False
  for line in (getattr(record, "final_text", "") or "").splitlines():
    stripped = line.strip()
    if stripped.startswith("```"):
      fenced = not fenced
      continue
    if fenced or not stripped.startswith("#"):
      continue
    heading = stripped.lstrip("#").strip()
    if heading:
      if len(heading) > _TITLE_MAX:
        heading = heading[:_TITLE_MAX - 1].rstrip() + "…"
      return heading
  # No heading: name the run rather than inventing a description of it.
  #
  # "Subagent", NOT "Review". This layer knows nothing about what a child was
  # spawned to do, and calling every one a review made a generic child --
  # anything spawned through the default subagent profile -- read as an
  # adversarial review in the lists that consume these records. A caller that
  # counts reviews then counts children that never reviewed anything. The
  # profile name is the only thing here that says what the child actually was,
  # so it carries the distinction.
  when = getattr(record, "finished_at", None) or getattr(
    record, "started_at", None)
  stamp = when.strftime("%Y-%m-%d %H:%M") if when is not None else "unfinished"
  return "Subagent (%s) — %s" % (
    getattr(record, "profile_name", "") or "subagent", stamp)


def digest_subject_paths(paths):
  """Hash each path now, returning ``{absolute path -> sha256}``.

  Computed HERE, by code, rather than accepted as a number in the task text.
  That is the whole point: the agent supplying the paths is the one whose
  later claim ("nothing changed since that review") the digest exists to
  check, so a digest it wrote itself would verify nothing.

  A path that cannot be read maps to ``""`` -- recorded as *unverifiable*,
  which a reader must not confuse with *unchanged*. Raising instead would let
  one unreadable file abort a review that was otherwise fine.
  """
  from libtbx.utils import sha256_hexdigest
  out = {}
  for raw in (paths or []):
    try:
      path = os.path.abspath(os.path.expanduser(str(raw)))
    except (TypeError, ValueError):
      continue
    try:
      # libtbx's chunked reader, not a hand-rolled loop: it is the same
      # algorithm three copies of this had drifted into writing separately.
      out[path] = sha256_hexdigest(path)
    except OSError:
      out[path] = ""
  return out


def normalize_subject_paths(value):
  """Validate the ``subject_paths`` argument; return a list of strings.

  Same reasoning as ``normalize_subagent_task``: the schema is a request, not
  a guarantee, and a model emitting a bare string here is an ordinary
  mis-generation. A single string is accepted as a one-element list rather
  than refused -- it is unambiguous -- while anything else is dropped, since
  a malformed identity must not block the review.
  """
  if value is None:
    return []
  if isinstance(value, str):
    return [value] if value.strip() else []
  if isinstance(value, (list, tuple)):
    return [v for v in value if isinstance(v, str) and v.strip()]
  return []


def normalize_subagent_task(task):
  """Validate a ``run_subagent`` ``task`` argument and return it stripped.

  PUBLIC because ``task`` arrives on TWO tool surfaces -- this module's registry
  builtin and ``claude_code``'s SDK tool -- which must refuse exactly the same
  inputs, or an argument that spawns a child on one backend family is refused
  on the other.

  The schema says ``task`` is a string, and a schema is a request, not a
  guarantee: a model emitting an object or an array here is an ordinary
  mis-generation, and the old ``(args.get("task") or "").strip()`` met it with
  a bare ``AttributeError`` -- an unhandled crash inside a tool dispatch rather
  than the ``Sorry`` this layer contracts to raise for bad input, and one whose
  message ("'dict' object has no attribute 'strip'") tells the model nothing it
  can correct. The type is named in the refusal so the retry is informed.

  Whitespace-only is refused for the pre-existing reason: it makes an empty user
  message (a provider 400), and a child with no instructions reports on nothing.

  Parameters
  ----------
  task : object
      The raw ``task`` argument as the model emitted it.

  Returns
  -------
  str
      The task, stripped.

  Raises
  ------
  libtbx.utils.Sorry
      When *task* is missing, not a string, or blank.
  """
  if task is None:
    raise Sorry("run_subagent requires a non-empty 'task'")
  if not isinstance(task, str):
    raise Sorry("run_subagent 'task' must be a string holding the child's "
                "complete instructions, not %s. Write the instructions as one "
                "string." % type(task).__name__)
  task = task.strip()
  if not task:
    raise Sorry("run_subagent requires a non-empty 'task'")
  return task


def subagents_permitted(profile, depth=0):
  """Whether *profile* permits spawning a child at *depth*.

  The same two fields ``register_run_subagent``'s handler checks when it is
  CALLED. Exposed so an assembler can apply the policy BEFORE building a tool,
  which the handler's own check cannot cover: a backend that runs its own
  agentic loop (claude_code) pins its tool set at client construction and never
  dispatches through ``AgentSession``, so its copy of ``run_subagent`` has no
  backstop and a disabled profile has to be honored earlier.

  It lives here, beside the handler, because ``Profile`` owns these fields and
  a second reading of them in an embedder is a second place for the policy to
  drift. Gating before the tool is built also keeps a tool that could only ever
  refuse off the model's list, instead of costing a turn to discover.

  Parameters
  ----------
  profile : Profile or None
      Read for ``subagents_enabled`` and ``subagents_max_depth``.
  depth : int, optional
      The prospective PARENT's depth; 0 for a top-level session. A child is
      permitted when the parent's depth is below the profile's limit.

  Returns
  -------
  bool
  """
  if not getattr(profile, "subagents_enabled", True):
    return False
  return depth < getattr(profile, "subagents_max_depth", 1)


def register_run_subagent(registry, build_child):
  """Register the ``run_subagent`` builtin on *registry*.

  Parameters
  ----------
  registry : ToolRegistry
      Registry the builtin is added to.
  build_child : callable
      ``build_child(profile_name, model, backend) -> ChildBundle``. Injected so
      this module never imports the backend-specific factory.
  """
  def _handler(name, input, cancel, session, tool_use_id):
    args = input or {}
    profile = getattr(session, "profile", None)
    if not getattr(profile, "subagents_enabled", True):
      raise Sorry("subagents are disabled for this profile")
    max_depth = getattr(profile, "subagents_max_depth", 1)
    if session.depth >= max_depth:
      raise Sorry("subagent depth %d reached the limit of %d; a child may "
                  "not spawn its own child" % (session.depth, max_depth))
    task = normalize_subagent_task(args.get("task"))
    bundle = build_child(args.get("profile"), args.get("model"),
                         args.get("backend"))
    # The PARENT profile's subagents.default_max_turns is exactly what its name
    # says -- a DEFAULT for the children this profile spawns. A child profile
    # that declares its OWN max_turns overrides it inside run_child, the one
    # layer that sees both numbers.
    record = run_child(session, task, bundle, cancel, tool_use_id,
                       max_turns=getattr(
                         profile, "subagents_default_max_turns", None),
                       subject_paths=normalize_subject_paths(
                         args.get("subject_paths")))
    return record.final_text

  spec = ToolSpec(name="run_subagent", description=RUN_SUBAGENT_DESCRIPTION,
                  input_schema=RUN_SUBAGENT_SCHEMA)
  registry.register_builtin(spec=spec, handler=_handler, risk="write")
