"""Scrollable list of MessageBubbles + inline ToolApprovalCards.

The view is the streaming target for the runner:
  - add_message(m)              - append a finalized Message bubble
  - start_assistant_bubble()    - create + return an in-progress bubble
  - finalize_assistant_bubble() - close the in-progress bubble
  - add_approval_request(req)   - coalesce by batch_id into one card

It does NOT own the AgentSession or QtAgentRunner; the chat window wires
runner signals to view methods, and the view emits approval_decided so the
window can push responses through the runner."""

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.agent.conversation import Message, now
from qttbx.widgets.chat.message_bubble import MessageBubble
from qttbx.widgets.chat.tool_approval import ToolApprovalCard


class ConversationView(QtWidgets.QScrollArea):

  approval_decided = QtCore.Signal(list)         # list[ToolApprovalResponse]
  image_clicked = QtCore.Signal(str, str)        # conv_id, sha256

  def __init__(self, parent=None, storage=None, conv_id=None):
    super().__init__(parent)
    self._storage = storage
    self._conv_id = conv_id
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

  def set_conversation(self, storage, conv_id):
    self._storage = storage
    self._conv_id = conv_id

  # ---- bubble API ----------------------------------------------------------

  def add_message(self, message):
    bubble = MessageBubble(message, parent=self._container,
                           storage=self._storage, conv_id=self._conv_id)
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
                           storage=self._storage, conv_id=self._conv_id)
    bubble.image_clicked.connect(self.image_clicked)
    self._insert_widget(bubble)
    self._bubbles.append(bubble)
    self._in_progress = bubble
    self._follow_bottom = True
    self._maybe_scroll_to_bottom()
    return bubble

  def finalize_assistant_bubble(self, stop_reason):
    if self._in_progress is not None:
      self._in_progress.message.stop_reason = stop_reason
      self._in_progress = None
    self._approval_by_batch.clear()                # batches don't carry over

  def in_progress_bubble(self):
    return self._in_progress

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
    """For ToolUseRequested or post-stream tool_use cells the agent
    materializes onto the current bubble."""
    if self._in_progress is None:
      self.start_assistant_bubble()
    self._in_progress.append_block(block)
    self._maybe_scroll_to_bottom()

  # ---- approval API --------------------------------------------------------

  def add_approval_request(self, req):
    """Append a request to an existing batched card if batch_id matches,
    otherwise create a new card."""
    if req.batch_id and req.batch_id in self._approval_by_batch:
      card = self._approval_by_batch[req.batch_id]
      card.set_requests(list(card._requests) + [req])
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

  # ---- clear ---------------------------------------------------------------

  def clear(self):                                       # noqa: A003
    for w in list(self._bubbles) + list(self._approval_cards):
      w.setParent(None)
      w.deleteLater()
    self._bubbles.clear()
    self._approval_cards.clear()
    self._approval_by_batch.clear()
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
    """Re-evaluate follow-mode after a real user scroll action. Qt
    fires actionTriggered for drag / wheel / pageup / etc., but NOT
    for programmatic setValue calls -- so layout-driven scrollbar
    updates (auto_height growing a bubble) never flip the flag."""
    if action == QtWidgets.QAbstractSlider.SliderNoAction:
      return
    # Defer until after Qt updates value() for the action.
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
