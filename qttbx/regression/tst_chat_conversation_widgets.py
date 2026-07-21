"""Conversation-area widgets: ConversationView (the scrolling stack of
MessageBubbles in the centre column) and ConversationList (the sidebar
of past conversations). Both compose MessageBubble / ToolApprovalCard
into the user-visible chat surface, just at different sides of the
splitter."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times, null_out

try:
  from qttbx.qt import QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, ConversationMeta, Message, now)
from qttbx.widgets.chat.agent.tools import ToolApprovalRequest


def _qapp():
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  return app


def _user(text="hi"):
  return Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": text})])


# The shape of error a force-killed Coot bridge produces: long, single-line,
# and routed to the failing tool's disclosure cell.
LONG_TOOL_ERROR = (
  "MCP error -32000: Connection closed. Failed to connect to the Coot RPC "
  "bridge at 127.0.0.1:44100: [Errno 61] Connection refused. The Coot "
  "process may have exited or been killed; restart Coot and retry.")


def _png_bytes(width=10, height=10):
  from qttbx.qt import QtCore, QtGui
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(255, 0, 0))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  return bytes(buf.data())


def _meta(i, title=None):
  ts = now()
  return ConversationMeta(
    id="id%d" % i, title="Conv %d" % i if title is None else title,
    profile_name="t", model="m", created_at=ts, updated_at=ts)


# ---- ConversationView ----------------------------------------------------


def exercise_add_bubble_then_streaming_then_finalize():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.add_message(_user("hello"))
  assert v.bubble_count() == 1
  bub = v.start_assistant_bubble()
  bub.append_text_delta("hi ")
  bub.append_text_delta("there")
  v.finalize_assistant_bubble("end_turn")
  assert v.bubble_count() == 2
  assert "there" in v.bubbles()[-1].combined_text()


def exercise_finalize_cancelled_marks_running_tool_cells_cancelled():
  """When a turn is cancelled (the user hit Stop), the in-progress bubble is
  finalized with stop_reason='cancelled'. Any tool cell still 'running' -- its
  result will never arrive -- must be transitioned to the cancelled terminal
  state rather than left stuck spinning."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  bub = v.start_assistant_bubble()
  v.append_block_to_current(ContentBlock(type="tool_use", data={
    "id": "t1", "name": "phenix_start_job", "input": {}}))
  cell = bub._tool_cells_by_id["t1"]
  assert cell.is_running()
  v.finalize_assistant_bubble("cancelled")
  assert not cell.is_running(), \
    "cancelled turn must terminate the still-running tool cell"
  assert "cancelled" in cell.header_button.text()


def exercise_finalize_non_cancel_leaves_running_cells_untouched():
  """A non-cancel finalize (end_turn) must NOT relabel a still-running tool
  cell -- the sweep is cancel-specific (only stop_reason='cancelled' calls
  cancel_running_tools). This pins that gate directly; it does NOT assume the
  cell was already finished -- on claude_code a live tool_use cell is only
  finished when its result is observed (see
  exercise_observed_result_finishes_cell_so_cancel_sweep_skips_it)."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  bub = v.start_assistant_bubble()
  v.append_block_to_current(ContentBlock(type="tool_use", data={
    "id": "t1", "name": "phenix_start_job", "input": {}}))
  cell = bub._tool_cells_by_id["t1"]
  v.finalize_assistant_bubble("end_turn")
  assert cell.is_running(), \
    "non-cancel finalize must leave running cells alone"


def exercise_observed_result_finishes_cell_so_cancel_sweep_skips_it():
  """[F#2] claude_code runs its tools in the SDK subprocess, bypassing the
  session dispatch that finishes an API backend's cell -- so nothing finishes
  the live tool_use cell and it stays 'running'. ConversationView's
  finish_tool_cell (driven by ChatWindow._on_tool_result_observed when a result
  is observed) must transition the matching in-progress cell to the finished
  terminal state. Otherwise the Stop-time sweep -- finalize_assistant_bubble(
  'cancelled') -> cancel_running_tools -- mislabels a tool that actually
  COMPLETED as 'cancelled'.

  Pins three things: (a) an observed cell becomes finished, (b) a later
  turn-cancel does NOT relabel that finished cell, and (c) a still-running
  sibling cell IS marked cancelled.

  Revert-proof: make finish_tool_cell a no-op and the cell stays running, so
  assertion (a) fails immediately -- and were it reached, the cancel sweep in
  (b) would relabel it 'cancelled'."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  bub = v.start_assistant_bubble()
  for tid in ("done", "running"):
    v.append_block_to_current(ContentBlock(type="tool_use", data={
      "id": tid, "name": "phenix_start_job", "input": {}}))
  done_cell = bub._tool_cells_by_id["done"]
  running_cell = bub._tool_cells_by_id["running"]
  assert done_cell.is_running() and running_cell.is_running()
  # (a) 'done's result is observed -> its cell transitions to finished.
  v.finish_tool_cell("done")
  assert not done_cell.is_running()
  assert "finished" in done_cell.header_button.text(), \
    done_cell.header_button.text()
  # The user hits Stop -> the cancel sweep runs.
  v.finalize_assistant_bubble("cancelled")
  # (b) the finished cell is left exactly as it was -- never 'cancelled'.
  assert "cancelled" not in done_cell.header_button.text(), \
    done_cell.header_button.text()
  assert "finished" in done_cell.header_button.text()
  # (c) the still-running sibling IS swept to cancelled.
  assert not running_cell.is_running()
  assert "cancelled" in running_cell.header_button.text()


def exercise_finish_tool_cell_is_a_no_op_without_in_progress_bubble():
  """finish_tool_cell must tolerate being called with no in-progress bubble
  (e.g. a stray late result after the turn already finalized) -- a guard so a
  mis-timed ToolResultObserved can't raise into the GUI thread."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.finish_tool_cell("nope")               # no in-progress bubble -> no raise
  v.start_assistant_bubble()
  v.finish_tool_cell("unknown_id")          # unknown id on the bubble -> no raise


def exercise_reload_folds_tool_results_into_matching_tool_cells():
  """[F#2] On reload, an answering tool_result message folds into the bubble
  holding its tool_use cells instead of becoming its own bubble: the reloaded
  turn shows one finished cell per tool (call + result), not a stuck-'running'
  call cell plus a detached 'result' cell.

  add_message(..., fold_tool_results=True) is the reload path (ChatWindow.
  _rebuild_view). Pins: (a) no extra bubble for the answer, (b) the matching
  cell leaves 'running' and shows 'finished', (c) an is_error result shows
  'failed', and (d) the default (live) path still gives the answer its own
  bubble and leaves the call cell running (the legacy shape this retires)."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()

  def _assistant_with_tool(tid):
    return Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="text", data={"text": "working"}),
      ContentBlock(type="tool_use", data={
        "id": tid, "name": "phenix_get_status", "input": {}})])

  def _answer(tid, text, is_error=False):
    return Message(role="user", timestamp=now(), content=[
      ContentBlock(type="tool_result", data={
        "tool_use_id": tid,
        "content": [ContentBlock(type="text", data={"text": text})],
        "is_error": is_error})])

  # Reload: an assistant tool_use turn, then its answering tool_result message.
  v = ConversationView()
  v.add_message(_assistant_with_tool("t1"), fold_tool_results=True)
  assert v.bubble_count() == 1
  cell = v.bubbles()[-1]._tool_cells_by_id["t1"]
  assert cell.is_running()                        # nothing has answered it yet
  v.add_message(_answer("t1", "Job 1: done"), fold_tool_results=True)
  # (a) folded -- no second bubble for the answer.
  assert v.bubble_count() == 1, \
    "answering tool_result must fold into the assistant bubble, not add one"
  # (b) the matching cell left 'running' and shows the finished result.
  assert not cell.is_running(), "folded result must finish the tool cell"
  assert "finished" in cell.header_button.text(), cell.header_button.text()

  # (c) a failed observed result folds as an error cell.
  v.add_message(_assistant_with_tool("t2"), fold_tool_results=True)
  err_cell = v.bubbles()[-1]._tool_cells_by_id["t2"]
  v.add_message(_answer("t2", "boom", is_error=True), fold_tool_results=True)
  assert not err_cell.is_running()
  assert "failed" in err_cell.header_button.text(), err_cell.header_button.text()

  # (d) default (live) path is unchanged: the answer is its own bubble and the
  # call cell stays 'running' -- the orphan/legacy reload shape this fix retires.
  v2 = ConversationView()
  v2.add_message(_assistant_with_tool("t3"))
  v2.add_message(_answer("t3", "done"))
  assert v2.bubble_count() == 2
  assert v2.bubbles()[0]._tool_cells_by_id["t3"].is_running()


def exercise_finish_tool_cell_error_marks_cell_failed():
  """A failed observed tool (is_error=True) finishes its LIVE cell as 'failed',
  not neutral 'finished' -- so a failure the user is watching is reported at
  once, matching the reloaded view (which renders the persisted is_error
  tool_result as a red cell). A successful finish stays 'finished'."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  bub = v.start_assistant_bubble()
  for tid in ("bad", "good"):
    v.append_block_to_current(ContentBlock(type="tool_use", data={
      "id": tid, "name": "Bash", "input": {}}))
  v.finish_tool_cell("bad", is_error=True, result="exit status 1")
  v.finish_tool_cell("good")
  bad = bub._tool_cells_by_id["bad"]
  good = bub._tool_cells_by_id["good"]
  assert not bad.is_running() and "failed" in bad.header_button.text(), \
    bad.header_button.text()
  assert not good.is_running() and "finished" in good.header_button.text(), \
    good.header_button.text()


def exercise_failed_tool_cell_does_not_floor_the_view_width():
  """A failed tool's error text must not pin the view's width.

  The bubbles track the window because ConversationView is a
  setWidgetResizable QScrollArea -- but only while nothing inside floors the
  container's minimum width. The disclosure header is a QToolButton, whose
  minimumSizeHint is its FULL text width and which never elides, so routing a
  long error into it floored the container: the bubbles stopped tracking the
  window, a horizontal scrollbar appeared, and every message was clipped.
  Seen for real when a force-killed Coot made mcp__coot__* calls fail with a
  long MCP connection error. One such cell anywhere in the conversation pins
  the whole view, even scrolled out of sight.
  """
  from qttbx.widgets.chat.conversation_view import ConversationView
  app = _qapp()
  v = ConversationView()
  bub = v.add_message(_user("hello"))
  bub.add_tool_use_cell(tool_id="t1", name="mcp__coot__run_python", args={})
  bub.set_tool_use_finished(tool_id="t1", error=LONG_TOOL_ERROR)
  v.resize(1100, 700)
  v.show()
  app.processEvents()
  v.resize(500, 700)
  app.processEvents()
  assert v.widget().width() <= v.viewport().width(), (
    "container %d exceeds viewport %d -- a failed tool cell floored the view"
    % (v.widget().width(), v.viewport().width()))
  assert v.horizontalScrollBar().maximum() == 0, (
    "horizontal scrollbar appeared (max=%d); bubbles no longer track the "
    "window width" % v.horizontalScrollBar().maximum())


def exercise_set_assistant_label_flows_to_new_bubbles():
  """set_assistant_label provides the fallback name for new assistant bubbles
  that carry no per-message backend stamp (the live streaming case)."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.set_assistant_label("Gemini")
  bubble = v.start_assistant_bubble()
  bubble.append_text_delta("hello")
  assert "Gemini" in bubble.first_text_cell_html()


def exercise_question_card_uses_assistant_label():
  """A question card rendered by the view names the active assistant."""
  from qttbx.widgets.chat.agent.events import AskUserQuestionRequested
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.set_assistant_label("Gemini")
  v.add_question_request(AskUserQuestionRequested(
    request_id="r1",
    questions=[{"question": "Which?",
                "options": [{"label": "a"}, {"label": "b"}]}]))
  card = v._question_cards[-1]
  labels = [w.text() for w in card.findChildren(QtWidgets.QLabel)]
  assert any("Gemini needs an answer" in t for t in labels), labels


def exercise_stop_finalizes_pending_question_cards():
  """A turn stopped while a QuestionCard is still awaiting an answer must
  finalize that card (disable Submit + hide it) -- the same sweep the pending
  approval cards get. Otherwise the stale card stays clickable and a late
  Submit carries a request_id the parked worker no longer waits on, so the
  answer is silently dropped. finalize_assistant_bubble is the normal/error
  turn-end path; on a parked cancel ChatWindow._on_stop sweeps immediately,
  then the session's synthesized terminal cancel TurnDone finalizes via this
  same path afterwards (queued)."""
  from qttbx.widgets.chat.agent.events import AskUserQuestionRequested
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.start_assistant_bubble()
  v.add_question_request(AskUserQuestionRequested(
    request_id="r1",
    questions=[{"question": "Which?",
                "options": [{"label": "a"}, {"label": "b"}]}]))
  card = v._question_cards[-1]
  assert card._buttons_widget.isEnabled()            # live: Submit clickable
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  # The user hits Stop -> the turn finalizes.
  v.finalize_assistant_bubble("cancelled")
  # The undecided question card must be swept: Submit disabled + card hidden.
  assert not card._buttons_widget.isEnabled(), \
    "stopped question card must have Submit disabled"
  assert card.isHidden(), "stopped question card must be hidden"
  # A late click on the finalized card must not emit a stale answer.
  card.click_submit()
  assert fired == [], "a finalized question card must not emit a late answer"


def exercise_batched_approval_coalesces_by_batch_id():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.add_message(_user("x"))
  v.start_assistant_bubble()
  reqs = [ToolApprovalRequest(
    request_id="r%d" % i, tool_name="t", tool_source="builtin",
    input={"i": i}, risk="write", batch_id="B")
    for i in range(3)]
  for r in reqs:
    v.add_approval_request(r)
  assert v.approval_card_count() == 1
  card = v.approval_cards()[0]
  assert len(card._requests) == 3


def exercise_two_batches_two_cards():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.start_assistant_bubble()
  for bid in ("B1", "B1", "B2"):
    v.add_approval_request(ToolApprovalRequest(
      request_id="r-" + bid + "-x", tool_name="t",
      tool_source="builtin", input={}, risk="write",
      batch_id=bid))
  assert v.approval_card_count() == 2


def exercise_solo_request_without_batch_id_gets_own_card():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.start_assistant_bubble()
  v.add_approval_request(ToolApprovalRequest(
    request_id="r1", tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id=None))
  v.add_approval_request(ToolApprovalRequest(
    request_id="r2", tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id=None))
  assert v.approval_card_count() == 2


def exercise_same_batch_after_decision_starts_a_fresh_card():
  """Multi-tool turns dispatch serially: the session emits approval
  request 2 only AFTER request 1 has been decided and its card hidden
  (it blocks on the approval queue between tools). A second same-batch
  request must therefore get its own fresh card -- appending it to the
  already-decided, hidden card strands the request and hangs the turn
  (the multi-tool approval deadlock)."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.start_assistant_bubble()
  v.add_approval_request(ToolApprovalRequest(
    request_id="t1", tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id="B"))
  assert v.approval_card_count() == 1
  card1 = v.approval_cards()[0]
  # The user approves the first tool; the card emits its decision + hides.
  card1.click_approve_all()
  # The session dispatches tool 1, then emits the second same-batch request.
  v.add_approval_request(ToolApprovalRequest(
    request_id="t2", tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id="B"))
  # A second, distinct card must carry t2 -- not an append to hidden card1.
  assert v.approval_card_count() == 2, v.approval_card_count()
  card2 = v.approval_cards()[1]
  assert card2 is not card1
  assert [r.request_id for r in card2._requests] == ["t2"]
  # The decided card keeps only its own request.
  assert [r.request_id for r in card1._requests] == ["t1"]


def exercise_clear():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.add_message(_user("hi"))
  v.add_message(_user("there"))
  v.clear()
  assert v.bubble_count() == 0


def exercise_image_click_propagates_to_view():
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    att = storage.store_attachment("c1", _png_bytes(), "image/png")
    v = ConversationView(storage=storage, conv_id="c1")
    m = Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="image", data={
        "attachment_sha256": att.sha256, "mime": "image/png"})])
    v.add_message(m)
    fired = []
    v.image_clicked.connect(lambda c, s: fired.append((c, s)))
    # Walk to the _ImageCell inside the bubble.
    bub = v.bubbles()[0]
    cells = [w for w in bub.findChildren(QtWidgets.QFrame)
             if hasattr(w, "click") and hasattr(w, "sha256")]
    assert cells
    cells[0].click()
    assert fired == [("c1", att.sha256)], fired
  finally:
    shutil.rmtree(tmp)


def exercise_add_message_always_scrolls_to_bottom():
  """A new message represents a new logical turn -- either the user
  just sent something (and wants to see it) or a complete assistant
  message just landed (and the user wants to see it). add_message
  re-asserts follow regardless of where the user had scrolled."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  app = _qapp()
  v = ConversationView()
  v.resize(300, 200)
  v.show()
  for i in range(20):
    v.add_message(_user("line %d" % i))
  app.processEvents()
  bar = v.verticalScrollBar()
  # Simulate the user scrolling up to re-read older content.
  v._follow_bottom = False
  bar.setValue(max(1, bar.maximum() // 3))
  app.processEvents()
  # User sends a new message. The view should scroll to show it.
  v.add_message(_user("new message after scrolling up"))
  app.processEvents()
  assert v._follow_bottom is True
  assert bar.value() >= bar.maximum() - 24, (bar.value(), bar.maximum())


def exercise_streaming_deltas_do_not_yank_user_who_scrolled_up():
  """During streaming, if the user has manually scrolled away from
  the bottom (follow flipped off), incoming text-deltas must NOT
  pull them back. The flag is the source of truth -- only a real
  user scroll action (or add_message / start_assistant_bubble)
  flips it."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  app = _qapp()
  v = ConversationView()
  v.resize(300, 200)
  v.show()
  for i in range(20):
    v.add_message(_user("line %d" % i))
  app.processEvents()
  # Start an assistant bubble (re-asserts follow); then 'user scrolls up'.
  v.start_assistant_bubble()
  app.processEvents()
  bar = v.verticalScrollBar()
  v._follow_bottom = False                     # explicit user-scrolled-away
  bar.setValue(max(1, bar.maximum() // 3))
  mid = bar.value()
  app.processEvents()
  # Stream deltas. The bubble grows, but the scrollbar must NOT jump.
  for chunk in ("Once upon a time ", "there was a ", "very long answer."):
    v.append_text_delta_to_current(chunk)
    app.processEvents()
  assert v._follow_bottom is False
  # Position should still be roughly where the user left it -- not at
  # the bottom. Allow some slack since the bubble's growth changes
  # bar.maximum() (value stays put but the gap from max widens).
  assert bar.value() < bar.maximum() - 24, (bar.value(), bar.maximum(), mid)


def exercise_scroll_up_disengages_follow_synchronously_during_growth():
  """During streaming, a manual scroll-up must disengage follow-mode
  SYNCHRONOUSLY. If _on_user_scroll_action only DEFERS the follow re-eval
  (singleShot), _follow_bottom is still True when the next delta's
  rangeChanged fires -> _on_range_changed snaps the viewport back to the
  bottom, and the deferred refresh (reading the snapped-back position) then
  re-asserts follow. The user could never scroll away to re-read mid-stream."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  bar = v.verticalScrollBar()
  bar.setRange(0, 500)
  v._follow_bottom = True
  bar.setValue(120)                    # user parked well above the bottom
  # A real user scroll-up action (wheel / arrow / drag). It must flip follow
  # off NOW, not on a later event-loop turn.
  v._on_user_scroll_action(QtWidgets.QAbstractSlider.SliderSingleStepSub)
  assert not v._follow_bottom, \
    "scroll-up did not disengage follow-mode synchronously"
  # A streaming delta grows the content before the deferred refresh runs;
  # rangeChanged must NOT yank the viewport back to the new bottom.
  bar.setRange(0, 600)
  v._on_range_changed(0, 600)
  assert bar.value() != bar.maximum(), (bar.value(), bar.maximum())


def exercise_range_change_while_following_snaps_to_new_max():
  """When the assistant bubble grows (auto_height resize, new cell
  appended), the scrollbar's range increases AFTER the synchronous
  delta handler returns. The follow-mode invariant is: at the end of
  the layout pass, the bar must be at the new maximum -- otherwise
  long streaming output disappears off the bottom of the viewport.

  Regression: _maybe_scroll_to_bottom's singleShot(0) fired before
  layout, so setValue(max) landed on the stale max. rangeChanged
  hook catches the post-layout grow."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  app = _qapp()
  v = ConversationView()
  v.resize(300, 200)
  v.show()
  v.start_assistant_bubble()
  app.processEvents()
  # Stream enough text that the bubble has to grow well past the
  # initial viewport height.
  for i in range(40):
    v.append_text_delta_to_current("line %02d of streamed output\n" % i)
    app.processEvents()
  bar = v.verticalScrollBar()
  assert v._follow_bottom is True
  assert bar.value() >= bar.maximum() - 24, (bar.value(), bar.maximum())


# ---- ConversationList ----------------------------------------------------


def exercise_list_populate_and_select():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1), _meta(2)])
  selected = []
  w.selected.connect(lambda cid: selected.append(cid))
  w.select_index(1)
  assert selected == ["id1"], selected


def exercise_list_set_active_marks_one_row_with_checkmark():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1), _meta(2)])
  w.set_active("id1")
  assert w._active_id == "id1", w._active_id
  # Exactly the active row carries a non-null (checkmark) icon; the others
  # carry a same-size transparent placeholder (null cacheKey differs).
  active_icon = w._list.item(1).icon()
  assert not active_icon.isNull(), "active row must show a checkmark icon"
  for i in (0, 2):
    assert w._list.item(i).icon().cacheKey() != active_icon.cacheKey(), i


def exercise_list_rebuild_preserves_active_checkmark():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id1")
  w.set_conversations([_meta(0), _meta(1)])          # rebuild
  assert not w._list.item(1).icon().isNull(), "rebuild must keep the checkmark"
  assert w._list.item(0).icon().cacheKey() != w._list.item(1).icon().cacheKey()


def exercise_list_checkmark_is_small_and_gutter_fixed():
  """The active-row checkmark is small (a pinned iconSize) and the inactive-row
  blank gutter is the SAME size, so the checkmark can't shift the conversation
  name relative to an inactive row."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  sz = w._list.iconSize()
  assert sz.width() == sz.height() and 0 < sz.width() <= 14, \
    (sz.width(), sz.height())                        # small + square
  # checkmark and blank render at the same footprint -> constant decoration
  # width -> the name never indents as active moves between rows
  assert w._check_icon.pixmap(sz).size() == w._blank_icon.pixmap(sz).size() == sz


def exercise_list_set_active_on_empty_title_emits_no_rename():
  """setIcon is setData(DecorationRole), which fires itemChanged. Without a
  blockSignals guard, an empty-title row (text 'Untitled', cached title '')
  would emit a spurious rename_requested(cid, 'Untitled')."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0, title=""), _meta(1, title="")])
  renames = []
  w.rename_requested.connect(lambda cid, t: renames.append((cid, t)))
  w.set_active("id0")
  assert renames == [], renames


def exercise_list_new_button_emits_signal():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  news = []
  w.new_requested.connect(lambda: news.append(True))
  w.click_new()
  assert news == [True]


def exercise_list_delete_button_emits_signal_for_selected():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.select_index(0)
  dels = []
  w.delete_requested.connect(lambda cid: dels.append(cid))
  w.click_delete()
  assert dels == ["id0"], dels


def exercise_list_items_are_editable_for_rename():
  """Every row must carry the Qt.ItemIsEditable flag so the in-place
  editor opens on F2 / the Rename button / double-click on the active
  row (a double-click on a non-active row activates it instead)."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  for i in range(w._list.count()):
    flags = w._list.item(i).flags()
    assert bool(flags & QtCore.Qt.ItemIsEditable), \
      "row %d must be editable for rename" % i


def exercise_list_rename_button_emits_signal_when_text_changes():
  """Clicking 'Rename' opens the editor; we simulate the user typing a
  new title by mutating the item text directly (Qt fires itemChanged
  whether the change came from the editor or a setText call). The
  signal must carry (conv_id, new_title)."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.select_index(0)
  renames = []
  w.rename_requested.connect(
    lambda cid, title: renames.append((cid, title)))
  # click_rename opens the editor; in headless tests the editor isn't
  # interactive, so we drive itemChanged the same way Qt would when
  # the user commits.
  w.click_rename()
  w._list.item(0).setText("Renamed conversation")
  assert renames == [("id0", "Renamed conversation")], renames
  # Cached meta is also updated so a subsequent identical rename is
  # treated as a no-op.
  w._list.item(0).setText("Renamed conversation")
  assert renames == [("id0", "Renamed conversation")], renames


def exercise_list_double_click_path_emits_rename_via_item_changed():
  """Double-clicking the ACTIVE row opens the editor; committing (simulated by
  a setText, as Qt fires itemChanged either way) emits rename_requested."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  renames = []
  w.rename_requested.connect(lambda cid, t: renames.append((cid, t)))
  w._on_item_double_clicked(w._list.item(0))         # active -> editor opens
  w._list.item(0).setText("Renamed conversation")
  assert renames == [("id0", "Renamed conversation")], renames


def exercise_list_double_click_non_active_emits_activated():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w._on_item_double_clicked(w._list.item(1))         # non-active
  assert got == ["id1"], got


def exercise_list_double_click_active_opens_editor_no_activated():
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w._on_item_double_clicked(w._list.item(0))         # active -> rename
  assert got == [], got
  assert w._list.state() == QtWidgets.QAbstractItemView.EditingState


def exercise_list_single_click_does_not_activate():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w.select_index(1)                                  # selection change only
  assert got == [], got


def exercise_list_empty_rename_is_rejected_and_reverted():
  """An empty new title must not emit and must restore the previous
  text in the list so the visible label doesn't end up blank."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0)])
  renames = []
  w.rename_requested.connect(
    lambda cid, title: renames.append((cid, title)))
  w._list.item(0).setText("")
  assert renames == [], renames
  assert w._list.item(0).text() == "Conv 0", w._list.item(0).text()


def exercise_list_rename_with_same_title_is_no_op():
  """Setting the same title (e.g., editor opened then committed
  unchanged) must NOT emit -- otherwise every editor-open cycle
  would churn storage."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0)])
  renames = []
  w.rename_requested.connect(
    lambda cid, title: renames.append((cid, title)))
  w._list.item(0).setText("Conv 0")
  assert renames == [], renames


def exercise_list_rename_normalises_padded_title_in_display():
  """A committed rename to a padded title emits and stores the whitespace-
  STRIPPED title, so the row's DISPLAYED text must be stripped to match --
  otherwise the visible row shows '  Notes  ' while the stored/emitted title is
  'Notes', a mismatch that lingers for a whole in-flight turn (the busy path
  skips the sidebar repopulate that would otherwise redraw it)."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0)])
  w.select_index(0)
  renames = []
  w.rename_requested.connect(
    lambda cid, title: renames.append((cid, title)))
  w.click_rename()                                   # opens the editor
  w._list.item(0).setText("  Notes  ")               # commit a padded title
  assert renames == [("id0", "Notes")], renames      # emitted stripped
  assert w._cached_title("id0") == "Notes", w._cached_title("id0")   # stored stripped
  assert w._list.item(0).text() == "Notes", w._list.item(0).text()  # displayed too


def exercise_list_rename_button_with_no_selection_is_no_op():
  """Pressing Rename with nothing selected must not crash; the
  editor open is just suppressed."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.click_rename()


def exercise_list_select_id_selects_without_activating():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1), _meta(2)])
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w.select_id("id2")
  assert w.selected_id() == "id2", w.selected_id()
  assert got == [], got                              # programmatic: no activate


def exercise_list_select_id_unknown_is_noop():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.select_index(0)
  w.select_id("does-not-exist")
  assert w.selected_id() == "id0", w.selected_id()   # selection unchanged


def _return_event():
  from qttbx.qt import QtCore, QtGui
  return QtGui.QKeyEvent(QtCore.QEvent.KeyPress, QtCore.Qt.Key_Return,
                         QtCore.Qt.NoModifier)


def exercise_list_return_on_non_active_activates():
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  w.select_index(1)                                  # non-active selected
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  handled = w.eventFilter(w._list, _return_event())
  assert handled is True, "Return must be consumed"
  assert got == ["id1"], got


def exercise_list_return_on_active_opens_no_editor():
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  w.select_index(0)                                  # active selected
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  handled = w.eventFilter(w._list, _return_event())
  assert handled is True, "Return must be consumed (no fall-through to editor)"
  assert got == [], got
  assert w._list.state() != QtWidgets.QAbstractItemView.EditingState


def exercise_list_return_after_rename_commit_does_not_activate():
  """Renaming a NON-active row (F2/Rename button) and pressing Return to commit
  must not then activate that row. The EditingState check is the primary guard;
  the closeEditor flag is backup. If this passes with the flag removed on both
  PySide2 and PySide6, the flag may be dropped."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w.set_active("id0")
  w.select_index(1)                                  # non-active row to rename
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w.click_rename()                                   # opens editor on row 1
  w._list.item(1).setText("Renamed")                 # commit (fires itemChanged)
  w._list.closePersistentEditor(w._list.item(1))     # editor closes on commit
  w._arm_return_suppression()                        # simulate closeEditor firing
  handled = w.eventFilter(w._list, _return_event())  # the propagated Return
  assert handled is True and got == [], (handled, got)


def exercise_list_locked_conversation_is_dimmed_but_selectable():
  """A conversation open in another process (locked_ids) renders DIMMED but stays
  selectable + enabled (so a click can re-check / offer Unlock); it is NOT
  editable, and double-clicking it emits no activated (can't switch to a locked
  one until it's unlocked)."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  it0, it1 = w._list.item(0), w._list.item(1)
  assert it0.flags() & QtCore.Qt.ItemIsEditable                # unlocked: editable
  assert it1.flags() & QtCore.Qt.ItemIsSelectable              # locked: selectable
  assert it1.flags() & QtCore.Qt.ItemIsEnabled                 #   + enabled
  assert not (it1.flags() & QtCore.Qt.ItemIsEditable)          #   but not editable
  assert it1.foreground() == w._dim_brush                      # dimmed
  assert it0.foreground() != w._dim_brush                      # unlocked: normal
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  w._on_item_double_clicked(it1)                               # locked -> inert
  assert got == [], got


def exercise_list_unlock_button_and_mark_unlocked():
  """set_unlock_button flips Delete -> Unlock (disabling Rename) and routes its
  click to unlock_requested; mark_unlocked re-styles a row normal (editable,
  undimmed) and drops it from the locked set; set_delete_button restores Delete +
  Rename."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  assert w._del_btn.text() == "Delete"                         # default
  got = []
  w.unlock_requested.connect(lambda cid: got.append(cid))
  w.set_unlock_button("id1")
  assert w._del_btn.text() == "Unlock" and w._unlock_target == "id1"
  assert not w._rename_btn.isEnabled()                         # locked: no rename
  w._on_del_button()                                           # routes to unlock
  assert got == ["id1"]
  w.mark_unlocked("id1")                                       # freed -> normal
  assert "id1" not in w._locked_ids
  it1 = w._list.item(1)
  assert (it1.flags() & QtCore.Qt.ItemIsEditable) \
    and it1.foreground() != w._dim_brush
  w.set_delete_button()
  assert w._del_btn.text() == "Delete" and w._unlock_target is None
  assert w._rename_btn.isEnabled()                             # restored
  # the button reserves the wider label's width, so the Delete<->Unlock toggle
  # can't resize it (which would nudge the sidebar)
  for _label in ("Delete", "Unlock"):
    w._del_btn.setText(_label)
    assert w._del_btn.minimumWidth() >= w._del_btn.sizeHint().width()


def exercise_list_mark_locked_dims_a_row():
  """mark_locked dims a row and adds it to the locked set (the inverse of
  mark_unlocked) -- used when an on-select re-check finds a conversation now open
  in another process. The _LOCKED_ROLE flag is set so it stays dimmed even when
  selected."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.conversation_list import ConversationList, _LOCKED_ROLE
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])          # none locked initially
  assert w.listed_ids() == ["id0", "id1"]            # rows in display order
  it1 = w._list.item(1)
  assert it1.foreground() != w._dim_brush
  assert it1.flags() & QtCore.Qt.ItemIsEditable
  w.mark_locked("id1")                               # a 2nd window opened
  assert "id1" in w._locked_ids
  assert it1.foreground() == w._dim_brush
  assert not (it1.flags() & QtCore.Qt.ItemIsEditable)
  assert it1.data(_LOCKED_ROLE) is True              # dimmed even when selected
  w.mark_locked("id1")                               # idempotent
  assert "id1" in w._locked_ids


def exercise_list_locked_row_stays_dimmed_when_selected():
  """A locked row stays dimmed even when it is the selected row: the item
  delegate forces HighlightedText (the selected row's text colour) to the dim
  colour, so the 'locked' cue survives selection. A normal row is unaffected."""
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.conversation_list import (
    ConversationList, _LockAwareDelegate)
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  d = w._list.itemDelegate()
  assert isinstance(d, _LockAwareDelegate)
  dim = w.palette().color(QtGui.QPalette.Disabled, QtGui.QPalette.Text)
  opt1 = QtWidgets.QStyleOptionViewItem()
  d.initStyleOption(opt1, w._list.model().index(1, 0))          # locked row
  assert opt1.palette.color(QtGui.QPalette.HighlightedText) == dim
  opt0 = QtWidgets.QStyleOptionViewItem()
  d.initStyleOption(opt0, w._list.model().index(0, 0))          # unlocked row
  assert opt0.palette.color(QtGui.QPalette.HighlightedText) != dim


def exercise_list_enter_on_locked_row_is_inert():
  """Enter/Return on a locked selected row does NOT activate it (parity with the
  double-click guard) -- it must be unlocked first. Deleting the eventFilter's
  locked-row clause would fail this."""
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  w.select_index(1)                                    # select the locked row
  got = []
  w.activated.connect(lambda cid: got.append(cid))
  ev = QtGui.QKeyEvent(QtCore.QEvent.KeyPress, QtCore.Qt.Key_Return,
                       QtCore.Qt.NoModifier)
  w.eventFilter(w._list, ev)
  assert got == [], "Enter on a locked row must not activate"


def exercise_list_non_left_double_click_is_ignored():
  """itemDoubleClicked is button-agnostic; the viewport filter swallows a
  non-left double-click so only a LEFT one activates/renames."""
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  vp = w._list.viewport()

  def _dbl(button):
    return QtGui.QMouseEvent(
      QtCore.QEvent.MouseButtonDblClick, QtCore.QPointF(5, 5),
      button, button, QtCore.Qt.NoModifier)

  assert w.eventFilter(vp, _dbl(QtCore.Qt.RightButton)) is True    # consumed
  assert w.eventFilter(vp, _dbl(QtCore.Qt.MiddleButton)) is True   # consumed
  # a LEFT double-click is NOT consumed here (falls through to itemDoubleClicked)
  assert w.eventFilter(vp, _dbl(QtCore.Qt.LeftButton)) is not True


def exercise_list_rebuild_resets_unlock_button_to_delete():
  """A rebuild (set_conversations) restores Delete + Rename even if a locked row
  had flipped the button to Unlock, so no stale Unlock target lingers."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  w.set_unlock_button("id1")
  assert w._del_btn.text() == "Unlock" and not w._rename_btn.isEnabled()
  w.set_conversations([_meta(0), _meta(1)])            # rebuild
  assert w._del_btn.text() == "Delete" and w._unlock_target is None
  assert w._rename_btn.isEnabled()


def exercise_list_lock_mid_rename_closes_editor_and_blocks_rename():
  """mark_locked on a row being renamed closes its editor, and even a forced
  commit on the now-locked row emits no rename_requested (the flag change alone
  would not retro-close the editor)."""
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  got = []
  w.rename_requested.connect(lambda cid, t: got.append((cid, t)))
  it1 = w._list.item(1)
  w._list.setCurrentItem(it1)
  w._list.editItem(it1)                                # open the editor
  assert w._list.state() == QtWidgets.QAbstractItemView.EditingState
  w.mark_locked("id1")                                 # a poll dims it mid-edit
  assert w._list.state() != QtWidgets.QAbstractItemView.EditingState, \
    "mark_locked must close the open editor"
  w._list.blockSignals(True)                           # force a stray commit
  it1.setText("hacked")
  w._list.blockSignals(False)
  w._on_item_changed(it1)
  assert got == [], "a locked row must not emit rename_requested"


def exercise_list_revert_last_rename_restores_pre_edit_title():
  """revert_last_rename restores a row's pre-edit title + cache after an
  optimistic in-place rename, without a rebuild -- backs the chat window's
  rename-refusal path (which can't rely on a rebuild a list failure would skip)."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  m1 = _meta(1)
  m1.title = "Orig"
  w.set_conversations([_meta(0), m1])
  got = []
  w.rename_requested.connect(lambda cid, t: got.append((cid, t)))
  it1 = w._list.item(1)
  w._list.blockSignals(True)
  it1.setText("Edited")
  w._list.blockSignals(False)
  w._on_item_changed(it1)                              # optimistic commit
  assert got == [("id1", "Edited")] and it1.text() == "Edited"
  w.revert_last_rename("id1")                          # refusal -> revert
  assert it1.text() == "Orig" and w._cached_title("id1") == "Orig"
  w.revert_last_rename("id1")                          # idempotent (already reverted)
  assert it1.text() == "Orig"


def exercise_list_action_button_debounces_on_label_flip():
  """The Delete/Unlock button disables itself briefly on a REAL label flip (so a
  click aimed at the old label drops), and does NOT re-arm when set to the same
  state twice (a steady poll must not keep it disabled)."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)], locked_ids={"id1"})
  assert w._del_btn.isEnabled()
  w.set_unlock_button("id1")                            # flip Delete -> Unlock
  assert not w._del_btn.isEnabled(), "a label flip debounces the button"
  w._enable_action_button()                             # (the timer would)
  w.set_unlock_button("id1")                            # same state -> no re-arm
  assert w._del_btn.isEnabled(), "a repeat in the same state must not debounce"


def exercise_list_selection_change_clears_stale_unlock_target():
  """Changing the selection must reset the shared Delete/Unlock button
  SYNCHRONOUSLY. The async lock re-check that re-arms Unlock for a locked row
  runs later (off-thread scan -> _sync_action_button_to_selection), so a stale
  'Unlock B' target left over from a previously-selected locked row must NOT
  survive into a newly-selected UNLOCKED row -- otherwise a click during the scan
  fires unlock_requested for B, the wrong, no-longer-selected conversation. The
  debounce only guards LABEL FLIPS, not selection changes, so the reset has to
  happen in _on_row_changed."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  # id1 (row 1) locked, id2 (row 2) not.
  w.set_conversations([_meta(0), _meta(1), _meta(2)], locked_ids={"id1"})
  w.select_index(1)                                     # select locked row B
  w.set_unlock_button("id1")                            # a scan armed Unlock B
  assert w._del_btn.text() == "Unlock" and w._unlock_target == "id1"
  w.select_index(2)                                     # select UNLOCKED row C
  # Before any async scan reconciles, the button must no longer be armed to
  # unlock B (else a click here unlocks the wrong conversation).
  assert w._unlock_target != "id1", w._unlock_target
  assert w._del_btn.text() == "Delete", w._del_btn.text()


def exercise_list_return_suppression_flag_auto_clears():
  """The post-rename Return guard must self-clear on the NEXT event-loop turn
  (QTimer.singleShot(0)), not linger. closeEditor also fires on a focus-out
  commit or an Escape revert -- cases where no Return follows -- so a
  never-cleared flag would silently swallow the user's next genuine Return.
  Arm it, pump the event loop, and assert it cleared: deleting the singleShot
  line leaves the flag armed and fails this."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.conversation_list import ConversationList
  app = _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0), _meta(1)])
  w._arm_return_suppression()
  assert w._suppress_next_return is True, "flag should arm"
  deadline = QtCore.QElapsedTimer()
  deadline.start()
  while w._suppress_next_return and deadline.elapsed() < 500:
    app.processEvents(QtCore.QEventLoop.AllEvents, 20)
  assert w._suppress_next_return is False, "flag must auto-clear next loop turn"


def exercise():
  exercise_add_bubble_then_streaming_then_finalize()
  exercise_finalize_cancelled_marks_running_tool_cells_cancelled()
  exercise_finalize_non_cancel_leaves_running_cells_untouched()
  exercise_observed_result_finishes_cell_so_cancel_sweep_skips_it()
  exercise_finish_tool_cell_is_a_no_op_without_in_progress_bubble()
  exercise_reload_folds_tool_results_into_matching_tool_cells()
  exercise_finish_tool_cell_error_marks_cell_failed()
  exercise_failed_tool_cell_does_not_floor_the_view_width()
  exercise_set_assistant_label_flows_to_new_bubbles()
  exercise_question_card_uses_assistant_label()
  exercise_stop_finalizes_pending_question_cards()
  exercise_batched_approval_coalesces_by_batch_id()
  exercise_two_batches_two_cards()
  exercise_solo_request_without_batch_id_gets_own_card()
  exercise_same_batch_after_decision_starts_a_fresh_card()
  exercise_clear()
  exercise_image_click_propagates_to_view()
  exercise_add_message_always_scrolls_to_bottom()
  exercise_streaming_deltas_do_not_yank_user_who_scrolled_up()
  exercise_scroll_up_disengages_follow_synchronously_during_growth()
  exercise_range_change_while_following_snaps_to_new_max()
  exercise_list_populate_and_select()
  exercise_list_set_active_marks_one_row_with_checkmark()
  exercise_list_rebuild_preserves_active_checkmark()
  exercise_list_checkmark_is_small_and_gutter_fixed()
  exercise_list_set_active_on_empty_title_emits_no_rename()
  exercise_list_new_button_emits_signal()
  exercise_list_delete_button_emits_signal_for_selected()
  exercise_list_items_are_editable_for_rename()
  exercise_list_rename_button_emits_signal_when_text_changes()
  exercise_list_double_click_path_emits_rename_via_item_changed()
  exercise_list_double_click_non_active_emits_activated()
  exercise_list_double_click_active_opens_editor_no_activated()
  exercise_list_single_click_does_not_activate()
  exercise_list_empty_rename_is_rejected_and_reverted()
  exercise_list_rename_with_same_title_is_no_op()
  exercise_list_rename_normalises_padded_title_in_display()
  exercise_list_rename_button_with_no_selection_is_no_op()
  exercise_list_select_id_selects_without_activating()
  exercise_list_select_id_unknown_is_noop()
  exercise_list_return_on_non_active_activates()
  exercise_list_return_on_active_opens_no_editor()
  exercise_list_return_after_rename_commit_does_not_activate()
  exercise_list_return_suppression_flag_auto_clears()
  exercise_list_locked_conversation_is_dimmed_but_selectable()
  exercise_list_unlock_button_and_mark_unlocked()
  exercise_list_mark_locked_dims_a_row()
  exercise_list_locked_row_stays_dimmed_when_selected()
  exercise_list_enter_on_locked_row_is_inert()
  exercise_list_non_left_double_click_is_ignored()
  exercise_list_rebuild_resets_unlock_button_to_delete()
  exercise_list_lock_mid_rename_closes_editor_and_blocks_rename()
  exercise_list_revert_last_rename_restores_pre_edit_title()
  exercise_list_action_button_debounces_on_label_flip()
  exercise_list_selection_change_clears_stale_unlock_target()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
