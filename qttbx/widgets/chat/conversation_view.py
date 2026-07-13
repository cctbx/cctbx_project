"""Scrollable list of message bubbles, approval cards, and question cards.

The view is the streaming target for the runner:
  - add_message(m)              - append a finalized Message bubble
  - start_assistant_bubble()    - create + return an in-progress bubble
  - finalize_assistant_bubble() - close the in-progress bubble
  - add_approval_request(req)   - coalesce concurrent same-batch
                                  requests; fresh card once decided
  - add_question_request(req)   - render a multi-choice question card

It does NOT own the AgentSession or QtAgentRunner; the chat window wires
runner signals to view methods, and the view emits approval_decided /
question_answered so the window can push responses through the runner.
"""

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.agent.conversation import Message, now
from qttbx.widgets.chat.message_bubble import MessageBubble
from qttbx.widgets.chat.question_card import QuestionCard
from qttbx.widgets.chat.tool_approval import ToolApprovalCard


def _is_tool_result_answer(message):
  """True for a user message whose blocks are ALL tool_result -- the answering
  message the session appends after an assistant tool_use turn (a dispatched
  batch, or claude_code's observed results). On reload these fold into the
  bubble holding the matching tool_use cells; a text or mixed user message (the
  user's own input, the turn-cap marker) is left as its own bubble.
  """
  if getattr(message, "role", None) != "user":
    return False
  content = getattr(message, "content", None) or []
  return bool(content) and all(
    getattr(b, "type", None) == "tool_result" for b in content)


class ConversationView(QtWidgets.QScrollArea):
  """Scrollable streaming view of bubbles, approval cards, and questions.

  Parameters
  ----------
  parent : QtWidgets.QWidget, optional
      Parent widget.
  storage : object, optional
      Attachment store passed to each ``MessageBubble`` for image lookup.
  conv_id : str, optional
      Identifier of the conversation being displayed.
  """

  approval_decided = QtCore.Signal(list)         # list[ToolApprovalResponse]
  question_answered = QtCore.Signal(str, dict)   # request_id, answers
  image_clicked = QtCore.Signal(str, str)        # conv_id, sha256

  def __init__(self, parent=None, storage=None, conv_id=None,
               assistant_label="Assistant"):
    super().__init__(parent)
    self._storage = storage
    self._conv_id = conv_id
    # Current backend's assistant display name; the fallback for bubbles /
    # question cards with no per-message backend stamp. Fixed for the
    # ChatWindow session (backend/profile don't change) and set here at
    # construction; set_assistant_label() can update it later (used by tests).
    self._assistant_label = assistant_label or "Assistant"
    self.setWidgetResizable(True)
    self._container = QtWidgets.QWidget(self)
    self._layout = QtWidgets.QVBoxLayout(self._container)
    self._layout.setContentsMargins(12, 8, 12, 8)
    # Spacing 0: MessageBubble carries its own 12 px top margin for
    # turn separation, so adding layout spacing here would double-space
    # the bubbles.
    self._layout.setSpacing(0)
    self._layout.addStretch(1)                     # bubbles pushed to top
    self.setWidget(self._container)
    self._bubbles = []                             # type: list[MessageBubble]
    self._approval_cards = []                      # type: list[ToolApprovalCard]
    self._approval_by_batch = {}                   # batch_id -> ToolApprovalCard
    self._question_cards = []                      # type: list[QuestionCard]
    self._in_progress = None
    # Sticky "follow the bottom" flag. True by default and re-asserted
    # whenever the user takes an action that should bring the latest
    # content into view (sends a message, the assistant starts replying).
    # Only flipped off when the user MANUALLY scrolls away from the
    # bottom (we listen for the scrollbar's actionTriggered signal,
    # which fires only on real user interaction, not programmatic
    # setValue calls). With auto_height growing bubbles in place, the
    # naive 'bar.value() >= bar.maximum() - 4' check broke after the
    # first growth: maximum jumps up while value stays put, so we'd
    # decide we were 'no longer at bottom' and stop following.
    self._follow_bottom = True
    bar = self.verticalScrollBar()
    bar.actionTriggered.connect(self._on_user_scroll_action)
    # When a bubble grows (auto_height resize, new cell appended) the
    # scrollbar's range increases. _maybe_scroll_to_bottom's
    # singleShot(0) fires before Qt processes that layout pass, so
    # setValue(maximum) lands on a stale max. Hook rangeChanged as
    # the authoritative "content grew" signal -- if we're following,
    # snap to the new max once it lands.
    bar.rangeChanged.connect(self._on_range_changed)

  def set_conversation(self, storage, conv_id):
    self._storage = storage
    self._conv_id = conv_id

  # ---- bubble API ----------------------------------------------------------

  def set_assistant_label(self, name):
    """Set the fallback assistant display name for new bubbles / cards.

    Stamped messages still show their own backend; this name is used only
    when a message carries no backend stamp (live streaming, legacy data).
    """
    self._assistant_label = name or "Assistant"

  def add_message(self, message, fold_tool_results=False):
    # On reload (fold_tool_results=True), an answering tool_result message folds
    # into the bubble that holds its tool_use cells (the preceding assistant
    # message) instead of becoming its own bubble: each result transitions its
    # tool_use cell out of 'running' and renders inline, so a reloaded turn
    # shows one cell per tool -- the call and its result together -- rather than
    # a bank of stuck-'running' call cells followed by a detached bank of
    # 'result' cells. A result with no matching cell falls back to its own cell
    # (MessageBubble._add_block's orphan path), so nothing is dropped. Live
    # callers keep the default: a streamed tool_result batch is its own bubble.
    if fold_tool_results and self._bubbles and _is_tool_result_answer(message):
      target = self._bubbles[-1]
      target.fold_tool_results(message.content)
      self._maybe_scroll_to_bottom()
      return target
    bubble = MessageBubble(message, parent=self._container,
                           storage=self._storage, conv_id=self._conv_id,
                           assistant_label=self._assistant_label)
    bubble.image_clicked.connect(self.image_clicked)
    self._insert_widget(bubble)
    self._bubbles.append(bubble)
    # A new message is something the user wants to see -- either they
    # just sent it themselves, or the runner is appending a response
    # for them. Re-assert follow regardless of where the scroll was.
    self._follow_bottom = True
    self._maybe_scroll_to_bottom()
    return bubble

  def start_assistant_bubble(self):
    if self._in_progress is not None:
      return self._in_progress
    msg = Message(role="assistant", content=[], timestamp=now())
    bubble = MessageBubble(msg, parent=self._container,
                           storage=self._storage, conv_id=self._conv_id,
                           assistant_label=self._assistant_label)
    bubble.image_clicked.connect(self.image_clicked)
    self._insert_widget(bubble)
    self._bubbles.append(bubble)
    self._in_progress = bubble
    self._follow_bottom = True
    self._maybe_scroll_to_bottom()
    return bubble

  def finalize_assistant_bubble(self, stop_reason):
    if self._in_progress is not None:
      # A cancelled turn (the user hit Stop) ends without the runner
      # delivering tool_results for any in-flight tool calls, so transition
      # the bubble's still-running tool cells to the cancelled terminal state
      # rather than leaving them stuck spinning.
      if stop_reason == "cancelled":
        self._in_progress.cancel_running_tools()
      self._in_progress.message.stop_reason = stop_reason
      self._in_progress = None
    # An approval card the turn left undecided must not stay clickable into a
    # later turn (approval-misroute guard); the same holds for a question card
    # the turn left unanswered (a late Submit would misroute).
    self.finalize_pending_approvals()
    self.finalize_pending_questions()
    self._approval_by_batch.clear()                # batches don't carry over

  def append_text_delta_to_current(self, text):
    if self._in_progress is None:
      self.start_assistant_bubble()
    self._in_progress.append_text_delta(text)
    self._maybe_scroll_to_bottom()

  def append_thinking_delta_to_current(self, text):
    if self._in_progress is None:
      self.start_assistant_bubble()
    self._in_progress.append_thinking_delta(text)
    self._maybe_scroll_to_bottom()

  def append_block_to_current(self, block):
    """Append a content block to the in-progress assistant bubble.

    For ToolUseRequested or post-stream tool_use cells the agent
    materializes onto the current bubble.
    """
    if self._in_progress is None:
      self.start_assistant_bubble()
    self._in_progress.append_block(block)
    self._maybe_scroll_to_bottom()

  def finish_tool_cell(self, tool_use_id, is_error=False, result=None):
    """Mark a tool cell on the in-progress bubble as finished or failed.

    The claude_code backend runs its tools inside the SDK subprocess, bypassing
    the session's dispatch loop -- the path that finishes an API backend's cell
    by delivering its tool_result. So nothing else finishes the live tool_use
    cell; observing the result is the signal that the tool completed, and this
    call transitions the matching in-progress cell to a terminal state. Without
    it the cell stays ``running`` until the turn ends, and
    ``finalize_assistant_bubble``'s Stop-time sweep (``cancel_running_tools``)
    would mislabel a COMPLETED tool ``cancelled``. No-op when there is no
    in-progress bubble or no cell matches ``tool_use_id`` (e.g. an API backend
    that already finished the cell itself).

    ``is_error`` transitions the cell to the failed state (carrying ``result``
    as the error text) rather than finished, so a failed in-subprocess tool is
    reported live in the same red state the reloaded view shows -- otherwise the
    failure reads as a neutral success until the next reload / conversation
    switch rebuilds the view.

    Parameters
    ----------
    tool_use_id : str
        Identifier of the tool_use cell to finish (matches the originating
        ``ToolUseRequested.id``).
    is_error : bool, optional
        Whether the observed result was an error.
    result : str, optional
        Flattened result text; shown as the error detail when ``is_error``.
    """
    if self._in_progress is None:
      return
    if is_error:
      self._in_progress.set_tool_use_finished(
        tool_use_id, error=result or "error")
    else:
      self._in_progress.set_tool_use_finished(tool_use_id)

  # ---- approval API --------------------------------------------------------

  def add_approval_request(self, req):
    """Render an approval card for ``req``.

    Same-batch_id requests coalesce into one card only while that card is
    still awaiting a decision -- a backend that emits a whole batch at
    once gets a single Approve-all card. The Anthropic backend instead
    dispatches serially, blocking on the approval queue between tools, so
    request 2 arrives only after request 1's card has been decided and
    hidden. That card is skipped here and a fresh, visible card is
    created; appending to the hidden card would strand the request and
    hang the turn.

    Parameters
    ----------
    req : ToolApprovalRequest
        The approval request to render. Its ``batch_id`` drives the
        coalescing behaviour described above.
    """
    existing = (self._approval_by_batch.get(req.batch_id)
                if req.batch_id else None)
    if existing is not None and not existing.is_decided():
      existing.set_requests(list(existing._requests) + [req])
      return
    card = ToolApprovalCard(self._container)
    card.set_requests([req])
    card.decided.connect(self._on_card_decided)
    self._insert_widget(card)
    self._approval_cards.append(card)
    if req.batch_id:
      self._approval_by_batch[req.batch_id] = card
    self._maybe_scroll_to_bottom()

  def _on_card_decided(self, responses):
    self.approval_decided.emit(responses)

  def finalize_pending_approvals(self):
    """Finalize (disable) every still-undecided approval card.

    A card left undecided when its turn ends must not stay clickable: a later
    click would route a stale response into the NEXT turn (the approval-misroute
    bug), since the approval path matches on ``is_busy()`` alone, not the
    request id. Called from ``ChatWindow._on_stop`` (the parked-cancel path,
    where no terminal ``TurnDone`` reaches the view) and from
    ``finalize_assistant_bubble`` (the normal / error turn-end paths).
    """
    for card in self._approval_cards:
      if not card.is_decided():
        card.finalize()

  # ---- question API --------------------------------------------------------

  def add_question_request(self, req):
    """Render a QuestionCard for an ``AskUserQuestionRequested``.

    One card per request -- there's no batching since the model groups
    its questions into a single tool call already.

    Parameters
    ----------
    req : AskUserQuestionRequested
        The question request, carrying ``request_id`` and ``questions``.
    """
    card = QuestionCard(self._container,
                        assistant_label=self._assistant_label)
    card.set_request(req.request_id, req.questions)
    card.answered.connect(self._on_question_answered)
    self._insert_widget(card)
    self._question_cards.append(card)
    self._maybe_scroll_to_bottom()

  def _on_question_answered(self, request_id, answers):
    self.question_answered.emit(request_id, answers)

  def finalize_pending_questions(self):
    """Finalize (disable) every still-unanswered question card.

    The QuestionCard analogue of ``finalize_pending_approvals``: a card left
    unanswered when its turn ends must not stay clickable, because a late
    Submit carries a ``request_id`` the parked worker no longer waits on -- the
    answer is silently dropped while the stale card lingers with Submit
    enabled. Called from the same turn-end paths: ``ChatWindow._on_stop`` (the
    parked-cancel path) and ``finalize_assistant_bubble`` (normal / error
    turn-end).
    """
    for card in self._question_cards:
      if not card.is_resolved():
        card.finalize()

  # ---- clear ---------------------------------------------------------------

  def clear(self):                                       # noqa: A003
    for w in (list(self._bubbles) + list(self._approval_cards)
              + list(self._question_cards)):
      w.setParent(None)
      w.deleteLater()
    self._bubbles.clear()
    self._approval_cards.clear()
    self._approval_by_batch.clear()
    self._question_cards.clear()
    self._in_progress = None

  # ---- introspection (for tests) ------------------------------------------

  def bubbles(self):
    return list(self._bubbles)

  def bubble_count(self):
    return len(self._bubbles)

  def approval_cards(self):
    return list(self._approval_cards)

  def approval_card_count(self):
    return len(self._approval_cards)

  # ---- internals -----------------------------------------------------------

  def _insert_widget(self, w):
    # Insert before the trailing stretch so widgets stack from the top.
    count = self._layout.count()
    self._layout.insertWidget(count - 1, w)

  def _on_user_scroll_action(self, action):
    """Re-evaluate follow-mode after a real user scroll action.

    Qt fires actionTriggered for drag / wheel / pageup / etc., but NOT
    for programmatic setValue calls -- so layout-driven scrollbar
    updates (auto_height growing a bubble) never flip the flag.
    """
    if action == QtWidgets.QAbstractSlider.SliderNoAction:
      return
    bar = self.verticalScrollBar()
    # sliderPosition already reflects this action's target -- Qt updates it
    # before emitting actionTriggered and commits value() only after -- so
    # disengage follow SYNCHRONOUSLY once the user has moved away from the
    # bottom. Deferring this (all we used to do) loses the race with a
    # streaming delta: rangeChanged fires _on_range_changed while _follow_bottom
    # is still True, snapping the viewport back to the bottom and stranding the
    # user, and the deferred refresh then reads that snapped-back position and
    # re-asserts follow. Re-engaging is still left to the deferred
    # _refresh_follow_state (which reads the committed value), so a scroll that
    # lands back at the bottom resumes following. Same 24 px band as the refresh.
    if bar.sliderPosition() < bar.maximum() - 24:
      self._follow_bottom = False
    QtCore.QTimer.singleShot(0, self._refresh_follow_state)

  def _refresh_follow_state(self):
    bar = self.verticalScrollBar()
    # Bigger threshold than the literal-bottom check: auto_height
    # growth happens between deltas, so a 4-px window would frequently
    # flip the flag off mid-stream. 24 px (~ one line of body text)
    # is a comfortable 'within follow distance' band.
    self._follow_bottom = bar.value() >= bar.maximum() - 24

  def _maybe_scroll_to_bottom(self):
    if not self._follow_bottom:
      return
    # Two-arg singleShot(msec, callable) is the only form supported on
    # both PySide2 and PySide6 -- PySide6 added a 3-arg
    # (msec, QObject, bound_method) overload that PySide2 lacks. The
    # bound method already carries self, so the receiver argument is
    # redundant; the trade-off is no auto-cancel-on-destroy guard, but
    # with msec=0 the receiver is virtually certain to still exist.
    QtCore.QTimer.singleShot(0, self._do_scroll_to_bottom)

  def _do_scroll_to_bottom(self):
    bar = self.verticalScrollBar()
    bar.setValue(bar.maximum())

  def _on_range_changed(self, _min, _max):
    """Snap to the bottom when content grows, if following.

    Fires whenever the scrollbar range changes -- i.e. whenever the
    contained content grows or shrinks. If we're following the bottom,
    snap to the new maximum now that the layout has settled.
    """
    if self._follow_bottom:
      bar = self.verticalScrollBar()
      bar.setValue(bar.maximum())
