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


def _png_bytes(width=10, height=10):
  from qttbx.qt import QtCore, QtGui
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(255, 0, 0))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  return bytes(buf.data())


def _meta(i):
  ts = now()
  return ConversationMeta(
    id="id%d" % i, title="Conv %d" % i, profile_name="t", model="m",
    created_at=ts, updated_at=ts)


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


def exercise_batched_approval_coalesces_by_batch_id():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.add_message(_user("x"))
  v.start_assistant_bubble()
  reqs = [ToolApprovalRequest(
    request_id="r%d" % i, tool_name="t", tool_source="builtin",
    input={"i": i}, risk="write", summary=None, batch_id="B")
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
      summary=None, batch_id=bid))
  assert v.approval_card_count() == 2


def exercise_solo_request_without_batch_id_gets_own_card():
  from qttbx.widgets.chat.conversation_view import ConversationView
  _qapp()
  v = ConversationView()
  v.start_assistant_bubble()
  v.add_approval_request(ToolApprovalRequest(
    request_id="r1", tool_name="t", tool_source="builtin",
    input={}, risk="write", summary=None, batch_id=None))
  v.add_approval_request(ToolApprovalRequest(
    request_id="r2", tool_name="t", tool_source="builtin",
    input={}, risk="write", summary=None, batch_id=None))
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
    input={}, risk="write", summary=None, batch_id="B"))
  assert v.approval_card_count() == 1
  card1 = v.approval_cards()[0]
  # The user approves the first tool; the card emits its decision + hides.
  card1.click_approve_all()
  # The session dispatches tool 1, then emits the second same-batch request.
  v.add_approval_request(ToolApprovalRequest(
    request_id="t2", tool_name="t", tool_source="builtin",
    input={}, risk="write", summary=None, batch_id="B"))
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
  editor opens on double-click / F2 / the Rename button."""
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
  """The double-click rename path goes through the same itemChanged
  slot as the button-driven editor. Pin that the slot emits the
  signal whether the user got there via double-click or the button."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.set_conversations([_meta(0)])
  renames = []
  w.rename_requested.connect(
    lambda cid, title: renames.append((cid, title)))
  w._list.item(0).setText("Edited via double-click")
  assert renames == [("id0", "Edited via double-click")], renames


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


def exercise_list_rename_button_with_no_selection_is_no_op():
  """Pressing Rename with nothing selected must not crash; the
  editor open is just suppressed."""
  from qttbx.widgets.chat.conversation_list import ConversationList
  _qapp()
  w = ConversationList()
  w.click_rename()


def exercise():
  exercise_add_bubble_then_streaming_then_finalize()
  exercise_batched_approval_coalesces_by_batch_id()
  exercise_two_batches_two_cards()
  exercise_solo_request_without_batch_id_gets_own_card()
  exercise_same_batch_after_decision_starts_a_fresh_card()
  exercise_clear()
  exercise_image_click_propagates_to_view()
  exercise_add_message_always_scrolls_to_bottom()
  exercise_streaming_deltas_do_not_yank_user_who_scrolled_up()
  exercise_range_change_while_following_snaps_to_new_max()
  exercise_list_populate_and_select()
  exercise_list_new_button_emits_signal()
  exercise_list_delete_button_emits_signal_for_selected()
  exercise_list_items_are_editable_for_rename()
  exercise_list_rename_button_emits_signal_when_text_changes()
  exercise_list_double_click_path_emits_rename_via_item_changed()
  exercise_list_empty_rename_is_rejected_and_reverted()
  exercise_list_rename_with_same_title_is_no_op()
  exercise_list_rename_button_with_no_selection_is_no_op()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
