"""Find-in-conversation tests: cell enumeration, scroll/pin plumbing,
ToolCallDisclosure.expand, the SearchBar widget, and the
ConversationSearch controller."""

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtCore, QtGui, QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.chat.agent.conversation import ContentBlock, Message, now


def _app():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  return app


def _pump(app, ms=100):
  timer = QtCore.QElapsedTimer()
  timer.start()
  while timer.elapsed() < ms:
    app.processEvents(QtCore.QEventLoop.AllEvents, 20)


def _text_message(text, role="user"):
  return Message(role=role, timestamp=now(), content=[
    ContentBlock(type="text", data={"text": text})])


def _mixed_message():
  """Assistant message with one cell of each searchable kind: a text
  cell, a thinking cell, and a tool disclosure (args + folded result)."""
  return Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "alpha needle one"}),
    ContentBlock(type="thinking", data={"text": "needle in thought"}),
    ContentBlock(type="tool_use", data={
      "id": "t1", "name": "grep", "input": {"pattern": "needle"}}),
    ContentBlock(type="tool_result", data={
      "tool_use_id": "t1", "is_error": False,
      "content": [ContentBlock(
        type="text", data={"text": "needle result line"})]})])


def exercise_bubble_searchable_cells_order_and_kinds():
  app = _app()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  b = MessageBubble(_mixed_message())
  cells = b.searchable_cells()
  kinds = [k for k, _w in cells]
  assert kinds == ["text", "thinking", "tool", "tool"], kinds
  # The tool pair is (args_view, result_view) of the same disclosure.
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  disc = b.findChildren(ToolCallDisclosure)[0]
  assert cells[2][1] is disc.args_view
  assert cells[3][1] is disc.result_view
  # Widgets carry the expected text.
  assert "alpha needle one" in cells[0][1].toPlainText()
  assert "needle in thought" in cells[1][1].toPlainText()


def exercise_bubble_enumerates_every_thinking_cell():
  """A bubble can hold several thinking cells; the _thinking_cell
  attribute only points at the LAST one (each thinking block replaces
  it), so enumeration must walk the layout."""
  app = _app()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="thinking", data={"text": "first thought"}),
    ContentBlock(type="text", data={"text": "interlude"}),
    ContentBlock(type="thinking", data={"text": "second thought"})])
  b = MessageBubble(m)
  kinds = [k for k, _w in b.searchable_cells()]
  assert kinds == ["thinking", "text", "thinking"], kinds
  texts = [w.toPlainText() for k, w in b.searchable_cells()
           if k == "thinking"]
  assert "first thought" in texts[0]
  assert "second thought" in texts[1]


def exercise_view_searchable_cells_flattens_bubbles_in_order():
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  v = ConversationView()
  v.add_message(_text_message("one"))
  v.add_message(_mixed_message())
  cells = v.searchable_cells()
  kinds = [k for k, _w in cells]
  assert kinds == ["text", "text", "thinking", "tool", "tool"], kinds
  assert "one" in cells[0][1].toPlainText()


def _filled_view(app, n=30):
  from qttbx.widgets.chat.conversation_view import ConversationView
  v = ConversationView()
  v.resize(400, 300)
  v.show()
  for i in range(n):
    v.add_message(_text_message("filler line %d" % i))
  _pump(app)
  return v


def exercise_ensure_visible_disengages_follow_and_holds():
  app = _app()
  v = _filled_view(app)
  assert v._follow_bottom is True
  _kind, w = v.searchable_cells()[0]
  cursor = QtGui.QTextCursor(w.document())
  cursor.setPosition(0)
  rect = w.cursorRect(cursor)
  v.ensure_child_rect_visible(w, rect)
  _pump(app)
  assert v._follow_bottom is False
  assert v._autofollow_held is True
  # The viewport actually moved up toward the first message.
  bar = v.verticalScrollBar()
  assert bar.value() < bar.maximum(), (bar.value(), bar.maximum())


def exercise_add_message_respects_autofollow_hold():
  app = _app()
  v = _filled_view(app, n=3)
  v.set_autofollow_held(True)
  v._follow_bottom = False
  v.add_message(_text_message("mid-turn tool results"))
  assert v._follow_bottom is False
  v.set_autofollow_held(False)
  v.add_message(_text_message("after release"))
  assert v._follow_bottom is True


def exercise_expand_recalcs_collapsed_result_height():
  """expand() must leave the inner result view at its real multi-line
  height SYNCHRONOUSLY-after-layout, not the ~1-line height it was
  given while hidden (sized at viewport width 0). The manual-click path
  defers this by an event-loop tick; the search reveal path cannot
  scroll against stale geometry."""
  app = _app()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  d = ToolCallDisclosure(name="grep", status="finished", parent=host)
  layout.addWidget(d)
  host.resize(400, 600)
  host.show()
  d.set_result("\n".join("line %d" % i for i in range(40)))
  _pump(app)
  assert not d.body.isVisible()
  one_line = d.result_view.fontMetrics().lineSpacing()
  d.expand()
  _pump(app)
  assert d.body.isVisible()
  assert d.header_button.isChecked()
  assert d.result_view.height() > 5 * one_line, d.result_view.height()
  # Idempotent: expanding an expanded row is a no-op, not a toggle.
  d.expand()
  assert d.body.isVisible()


def exercise_expand_calls_refresh_inner_heights_synchronously():
  """expand() must run the inner-height refresh SYNCHRONOUSLY (before it
  returns): the scroll-to-match caller measures geometry in the same
  tick. Offscreen Qt corrects the height via setVisible's synchronous
  resize regardless -- a height assertion passes with or without the
  explicit refresh -- so the only way to pin the production contract
  (where that resize is async) is to observe the call itself. Shadow
  the refresh with a counter, expand once and assert it fired exactly
  once with NO pump, then expand again (already open) and assert the
  early-return no-op path does not fire it a second time."""
  app = _app()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  d = ToolCallDisclosure(name="grep", status="finished", parent=host)
  layout.addWidget(d)
  host.resize(400, 600)
  host.show()
  d.set_result("\n".join("line %d" % i for i in range(40)))
  _pump(app)
  counter = {"n": 0}
  original = d._refresh_inner_heights
  def counting_refresh():
    counter["n"] += 1
    original()
  d._refresh_inner_heights = counting_refresh
  d.expand()
  assert counter["n"] == 1, counter["n"]   # synchronous, before any pump
  d.expand()                               # already expanded -> no-op
  assert counter["n"] == 1, counter["n"]


def exercise_search_bar_signals_and_defaults():
  app = _app()
  from qttbx.widgets.chat.search_bar import SearchBar
  bar = SearchBar()
  assert bar.included_kinds() == {"text", "thinking"}
  got = {"q": [], "next": 0, "prev": 0, "scope": 0, "close": 0}
  bar.query_changed.connect(lambda t: got["q"].append(t))
  bar.next_requested.connect(
    lambda: got.__setitem__("next", got["next"] + 1))
  bar.prev_requested.connect(
    lambda: got.__setitem__("prev", got["prev"] + 1))
  bar.scope_changed.connect(
    lambda: got.__setitem__("scope", got["scope"] + 1))
  bar.close_requested.connect(
    lambda: got.__setitem__("close", got["close"] + 1))
  bar._edit.setText("needle")
  assert got["q"] == ["needle"], got["q"]
  bar._next_btn.click()
  bar._prev_btn.click()
  bar._scope_boxes["thinking"].setChecked(False)
  bar._scope_boxes["tool"].setChecked(True)
  bar._close_btn.click()
  assert (got["next"], got["prev"], got["scope"], got["close"]) \
      == (1, 1, 2, 1), got
  assert bar.included_kinds() == {"text", "tool"}
  assert bar.query() == "needle"


def exercise_search_bar_count_states():
  app = _app()
  from qttbx.widgets.chat.search_bar import SearchBar
  bar = SearchBar()
  bar.set_count(3, 17)
  assert bar._count.text() == "3/17", bar._count.text()
  bar.set_count(0, 0)
  assert bar._count.text() == "0/0", bar._count.text()
  bar.set_count(17, 500, capped=True)
  assert bar._count.text() == "17/500+", bar._count.text()
  bar.clear_count()
  assert bar._count.text() == ""


def exercise_search_bar_keys_from_all_children():
  """Enter -> next, Shift+Enter -> previous, Escape -> close, with focus
  on ANY interactive child of the bar -- the line edit, a checkbox that
  took click focus, or a button. Escape's ShortcutOverride is claimed
  everywhere too, so the window-level Escape=stop shortcut cannot fire
  from inside the bar and kill a running turn."""
  app = _app()
  from qttbx.widgets.chat.search_bar import SearchBar
  bar = SearchBar()
  got = {"next": 0, "prev": 0, "close": 0}
  bar.next_requested.connect(
    lambda: got.__setitem__("next", got["next"] + 1))
  bar.prev_requested.connect(
    lambda: got.__setitem__("prev", got["prev"] + 1))
  bar.close_requested.connect(
    lambda: got.__setitem__("close", got["close"] + 1))
  children = (bar._edit, bar._prev_btn, bar._next_btn,
              bar._scope_boxes["thinking"], bar._scope_boxes["tool"],
              bar._close_btn)
  for child in children:
    ov = QtGui.QKeyEvent(QtCore.QEvent.ShortcutOverride,
                         QtCore.Qt.Key_Escape, QtCore.Qt.NoModifier)
    app.sendEvent(child, ov)
    assert ov.isAccepted(), child
  def key(target, k, mods=QtCore.Qt.NoModifier):
    app.sendEvent(target, QtGui.QKeyEvent(QtCore.QEvent.KeyPress, k, mods))
  key(bar._edit, QtCore.Qt.Key_Return)
  key(bar._edit, QtCore.Qt.Key_Enter)                     # numpad Enter
  key(bar._edit, QtCore.Qt.Key_Return, QtCore.Qt.ShiftModifier)
  key(bar._edit, QtCore.Qt.Key_Escape)
  key(bar._scope_boxes["thinking"], QtCore.Qt.Key_Return)
  key(bar._scope_boxes["tool"], QtCore.Qt.Key_Escape)
  key(bar._close_btn, QtCore.Qt.Key_Return, QtCore.Qt.ShiftModifier)
  assert (got["next"], got["prev"], got["close"]) == (3, 2, 2), got


def _controller_fixture(app):
  """View with two mixed messages + a ConversationSearch bound to it."""
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.add_message(_text_message("first needle here"))
  v.add_message(_mixed_message())
  _pump(app)
  return v, ConversationSearch(v)


def exercise_controller_open_close_toggle():
  app = _app()
  v, cs = _controller_fixture(app)
  assert not cs.is_open()
  closed = []
  cs.closed.connect(lambda: closed.append(True))
  cs.toggle()
  assert cs.is_open() and cs.bar.isVisible()
  # Overlay sits inside the view, near the top-right.
  assert cs.bar.parent() is v
  assert cs.bar.y() == 8
  assert cs.bar.x() + cs.bar.width() <= v.width()
  cs.toggle()
  assert not cs.is_open()
  assert closed == [True]
  # Reopening retains the previous query, pre-selected for overtyping.
  cs.open()
  cs.bar._edit.setText("kept")
  cs.close()
  cs.open()
  assert cs.bar.query() == "kept"
  assert cs.bar._edit.selectedText() == "kept"


def exercise_scan_default_scope_and_counts():
  app = _app()
  v, cs = _controller_fixture(app)
  cs.open()
  cs.bar._edit.setText("needle")
  # Default scope: text + thinking (Thinking starts checked) --
  # "first needle here", "alpha needle one", "needle in thought".
  # Tool matches are opt-in.
  assert len(cs._matches) == 3, cs._matches
  assert cs._current == 0
  assert cs.bar._count.text() == "1/3", cs.bar._count.text()
  # Scope checkboxes narrow/widen the scan: -1 thinking, +2 tool (the
  # args JSON and the result line each contain 'needle').
  cs.bar._scope_boxes["thinking"].setChecked(False)
  assert len(cs._matches) == 2, len(cs._matches)
  cs.bar._scope_boxes["thinking"].setChecked(True)
  assert len(cs._matches) == 3, len(cs._matches)
  cs.bar._scope_boxes["tool"].setChecked(True)
  assert len(cs._matches) == 5, len(cs._matches)
  assert cs.bar._count.text() == "1/5", cs.bar._count.text()
  # Case-insensitive.
  cs.bar._edit.setText("NEEDLE")
  assert len(cs._matches) == 5
  # No matches -> 0/0; empty query -> blank label, no matches.
  cs.bar._edit.setText("zzz-not-there")
  assert cs.bar._count.text() == "0/0"
  cs.bar._edit.setText("")
  assert cs.bar._count.text() == ""
  assert cs._matches == []
  # Emptying the query releases the autofollow hold (spec: pin releases on
  # close, send, and query-goes-empty). Re-populate first so clearing it
  # transitions non-empty -> empty and actually fires the edit signal.
  cs.bar._edit.setText("needle")
  v.set_autofollow_held(True)
  cs.bar._edit.setText("")
  assert v._autofollow_held is False


def exercise_highlights_applied_and_cleared():
  app = _app()
  v, cs = _controller_fixture(app)
  cs.open()
  cs.bar._edit.setText("needle")
  first_widget = cs._matches[0][0]
  sels = first_widget.extraSelections()
  assert len(sels) == 1, len(sels)
  # Current match is orange; a non-current match in another cell is
  # yellow. Both force black text.
  assert sels[0].format.background().color().name().upper() == "#FF9800"
  assert sels[0].format.foreground().color().name().upper() == "#000000"
  second_widget = cs._matches[1][0]
  sels2 = second_widget.extraSelections()
  assert sels2[0].format.background().color().name().upper() == "#FFF176"
  # Close clears every highlighted cell and releases the pin.
  v.set_autofollow_held(True)
  cs.close()
  assert first_widget.extraSelections() == []
  assert second_widget.extraSelections() == []
  assert v._autofollow_held is False


def exercise_match_cap():
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import (
    ConversationSearch, MAX_MATCHES)
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.add_message(_text_message("x " * (MAX_MATCHES + 50)))
  _pump(app)
  cs = ConversationSearch(v)
  cs.open()
  cs.bar._edit.setText("x")
  assert len(cs._matches) == MAX_MATCHES
  assert cs._capped is True
  assert cs.bar._count.text() == "1/500+", cs.bar._count.text()


def exercise_navigation_steps_and_wraps():
  app = _app()
  v, cs = _controller_fixture(app)
  cs.open()
  cs.bar._edit.setText("needle")     # 2 text + 1 thinking match
  assert (cs._current, len(cs._matches)) == (0, 3)
  cs.bar.next_requested.emit()
  assert cs._current == 1
  assert cs.bar._count.text() == "2/3"
  cs.bar.next_requested.emit()
  assert cs._current == 2
  cs.bar.next_requested.emit()       # wraps forward
  assert cs._current == 0
  cs.bar.prev_requested.emit()       # wraps backward
  assert cs._current == 2
  _pump(app)                         # let the deferred scroll run
  assert v._autofollow_held is True
  assert v._follow_bottom is False


def exercise_navigation_relocates_after_offset_shift():
  """Streaming markdown re-sets its document, re-flowing offsets. The
  nav re-scan relocates by exact (widget, start), else nearest offset
  in the same widget -- never jumping back to the first match."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.add_message(_text_message("lead bo tail bo end"))
  _pump(app)
  cs = ConversationSearch(v)
  cs.open()
  cs.bar._edit.setText("bo")
  assert len(cs._matches) == 2
  cs.bar.next_requested.emit()       # current -> second match
  assert cs._current == 1
  # Simulate a markdown re-flow that shifts earlier text: bold markers
  # around 'lead' vanish from the plain-text projection on re-render.
  view_widget = cs._matches[1][0]
  view_widget.clear()
  view_widget.append_markdown("**lead** bo tail bo end")
  # Plain text is now "lead bo tail bo end" rendered from new raw; the
  # second 'bo' offset differs from the recorded one only if markers
  # shifted -- either way the relocation must stay on the SAME
  # occurrence (nearest offset), not reset to match 0.
  cs.bar.next_requested.emit()
  # From the relocated second match, +1 wraps to the first.
  assert cs._current == 0, cs._current
  assert len(cs._matches) == 2


def exercise_stale_deferred_scroll_noop_after_close():
  """A next/prev step schedules the scroll-to-match with singleShot(0);
  closing the bar before that deferred call runs must leave it a no-op --
  it must not re-engage the hold (via ensure_child_rect_visible) a search the user already
  closed. Without the is_open() guard in _scroll_to the deferred scroll
  fires against the still-alive controller and re-pins the view."""
  app = _app()
  v, cs = _controller_fixture(app)
  cs.open()
  cs.bar._edit.setText("needle")     # >= 1 match
  assert len(cs._matches) >= 1
  cs.bar.next_requested.emit()       # schedules the deferred scroll
  cs.close()                         # close BEFORE the deferred scroll runs
  _pump(app)                         # let the singleShot(0) fire
  assert v._autofollow_held is False


def exercise_same_cell_current_vs_other_highlight():
  """Two matches in ONE cell: the current match renders orange while its
  same-cell neighbour stays yellow (the cross-cell test can't catch a
  per-widget-uniform-color bug)."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.add_message(_text_message("needle then needle"))
  cs = ConversationSearch(v)
  cs.open()
  cs.bar._edit.setText("needle")
  assert len(cs._matches) == 2
  assert cs._matches[0][0] is cs._matches[1][0]      # same widget
  sels = cs._matches[0][0].extraSelections()
  assert len(sels) == 2, len(sels)
  assert sels[0].format.background().color().name().upper() == "#FF9800"
  assert sels[1].format.background().color().name().upper() == "#FFF176"


def exercise_overlay_repositions_when_scrollbar_appears():
  """Content growth can summon the vertical scrollbar without any view
  Resize; the overlay must shift left so the scrollbar never overlaps
  it."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.add_message(_text_message("just one line"))
  _pump(app)
  cs = ConversationSearch(v)
  cs.open()
  _pump(app)
  sb = v.verticalScrollBar()
  assert not sb.isVisible()
  x_without_sb = cs.bar.x()
  for i in range(30):
    v.add_message(_text_message("filler line %d" % i))
  _pump(app)
  assert sb.isVisible()
  assert cs.bar.x() < x_without_sb, (cs.bar.x(), x_without_sb)
  assert cs.bar.x() + cs.bar.width() <= v.width() - sb.width()


def exercise_scroll_to_top_and_bottom():
  """scroll_to_top parks the viewport at the first message and stops
  following; scroll_to_bottom returns to the latest message, resumes
  following, and releases the autofollow hold (an explicit jump to the end
  abandons the match the pin was protecting)."""
  app = _app()
  v = _filled_view(app)
  bar = v.verticalScrollBar()
  assert bar.maximum() > 0 and bar.value() == bar.maximum()
  v.scroll_to_top()
  assert bar.value() == 0, bar.value()
  assert v._follow_bottom is False
  v.set_autofollow_held(True)
  v.scroll_to_bottom()
  assert bar.value() == bar.maximum()
  assert v._follow_bottom is True
  assert v._autofollow_held is False
  # scroll_to_top leaves an engaged pin alone: follow is off either way,
  # and the pin keeps suppressing add_message's follow re-assert while
  # the search stays open on a match.
  v.set_autofollow_held(True)
  v.scroll_to_top()
  assert v._autofollow_held is True
  assert v._follow_bottom is False


def exercise_find_next_and_previous_controller_api():
  """find_next()/find_previous() drive navigation without the bar's own
  keys: with the bar closed they open it and reveal the retained
  query's current match WITHOUT stepping (reopening must not skip
  match 1); with it open they step exactly like Enter / Shift+Enter."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.activateWindow()
  cs = ConversationSearch(v)
  v.add_message(_text_message("needle one sits at the top"))
  v.add_message(_text_message("needle two just below"))
  for i in range(30):
    v.add_message(_text_message("filler line %d" % i))
  _pump(app)
  # Retain a query, then close: the shortcut path starts from a closed
  # bar with a remembered query.
  cs.open()
  cs.bar._edit.setText("needle")
  cs.close()
  bar = v.verticalScrollBar()
  assert bar.value() == bar.maximum()
  cs.find_next()                       # closed -> open + reveal, no step
  _pump(app)                           # deferred scroll-to-match
  assert cs.is_open()
  assert cs._current == 0, cs._current
  assert cs.bar._count.text() == "1/2", cs.bar._count.text()
  assert v._autofollow_held is True
  assert bar.value() < bar.maximum()   # left the bottom for match 1
  cs.find_next()                       # open -> step, like Enter
  assert cs._current == 1, cs._current
  cs.find_previous()                   # open -> step back
  assert cs._current == 0, cs._current
  # find_previous from closed also reveals the current match, no step.
  cs.close()
  cs.find_previous()
  _pump(app)
  assert cs.is_open()
  assert cs._current == 0, cs._current
  # With no retained query the shortcut just opens the (blank) bar.
  cs.bar._edit.setText("")
  cs.close()
  cs.find_next()
  assert cs.is_open()
  assert cs._matches == []
  assert cs.bar._count.text() == ""


def exercise_close_preserves_scroll_position():
  """Closing the bar must not move the viewport. Hiding the focused bar
  makes Qt hand focus down the tab chain; if that lands on a message
  cell, QScrollArea's focusNextPrevChild scrolls that cell into view --
  yanking a bottom-parked conversation back to its first bubble. The
  controller must park focus on the view BEFORE hiding the bar.

  The controller is created BEFORE the cells (the ChatWindow
  construction order) -- load-bearing: the cells then sit after the bar
  in the tab chain, which is the arrangement that reproduces the jump.
  """
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  v.activateWindow()
  cs = ConversationSearch(v)
  for i in range(40):
    v.add_message(_text_message("filler line %d" % i))
  _pump(app)
  bar = v.verticalScrollBar()
  assert bar.maximum() > 0 and bar.value() == bar.maximum()
  bottom = bar.value()
  cs.open()                            # focuses the query edit
  _pump(app, 50)
  assert bar.value() == bottom, (bar.value(), bottom)
  cs.close()
  _pump(app)
  assert bar.value() == bottom, (bar.value(), bottom)
  # Parked mid-conversation: open/close must leave the position alone
  # too (close must not re-engage follow either).
  mid = bar.maximum() // 2
  bar.setValue(mid)
  v._follow_bottom = False
  cs.open()
  _pump(app, 50)
  assert bar.value() == mid, (bar.value(), mid)
  cs.close()
  _pump(app)
  assert bar.value() == mid, (bar.value(), mid)
  assert v._follow_bottom is False


def exercise_bar_width_tracks_view():
  """The overlay scales with the conversation area: about half the
  view's width, never narrower than its content hint, and re-sized live
  when the view resizes."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(1200, 300)
  v.show()
  v.add_message(_text_message("one line"))
  _pump(app)
  cs = ConversationSearch(v)
  cs.open()
  _pump(app, 50)
  hint_w = cs.bar.sizeHint().width()
  assert hint_w < 600, hint_w          # precondition for the half checks
  assert cs.bar.width() == 600, cs.bar.width()
  # Growing the view grows the bar (the view Resize repositions it).
  v.resize(1600, 300)
  _pump(app, 50)
  assert cs.bar.width() == 800, cs.bar.width()
  assert cs.bar.x() + cs.bar.width() <= v.width()
  # A view too narrow for half-width falls back to the content hint
  # (capped to the view, as before).
  v.resize(500, 300)
  _pump(app, 50)
  assert cs.bar.width() == min(hint_w, 500 - 16), cs.bar.width()
  assert cs.bar.x() + cs.bar.width() <= v.width()


def exercise_navigation_pin_engages_before_deferred_scroll():
  """The pin must engage SYNCHRONOUSLY on navigation, not in the deferred
  scroll. Otherwise a message arriving between Enter and the deferred
  tick sees _follow_bottom still True, schedules its own deferred
  snap-to-bottom, and yanks the viewport off the match the user just
  navigated to (the snap fires after the match scroll)."""
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  v = ConversationView()
  v.resize(400, 300)
  v.show()
  v.add_message(_text_message("the needle sits at the top"))
  for i in range(30):
    v.add_message(_text_message("filler line %d" % i))
  _pump(app)
  cs = ConversationSearch(v)
  cs.open()
  cs.bar._edit.setText("needle")
  assert len(cs._matches) >= 1
  assert v._follow_bottom is True          # a live stream would follow
  cs.bar.next_requested.emit()             # schedules the deferred scroll
  # Pin engaged in the SAME tick -- before any interleaved message.
  assert v._autofollow_held is True
  assert v._follow_bottom is False
  # A mid-turn message lands before the deferred scroll fires; it must
  # neither re-assert follow nor schedule a snap that outruns the match
  # scroll.
  v.add_message(_text_message("streamed mid-turn tool results"))
  _pump(app)
  bar = v.verticalScrollBar()
  assert bar.value() < bar.maximum(), (bar.value(), bar.maximum())


def exercise_deferred_bottom_snap_rechecks_follow():
  """A deferred snap-to-bottom scheduled while following must no-op if
  follow was disengaged (search navigation) before it fires."""
  app = _app()
  v = _filled_view(app)
  v.verticalScrollBar().setValue(0)
  v._follow_bottom = True
  v._maybe_scroll_to_bottom()              # schedules the deferred snap
  v._follow_bottom = False                 # disengaged before it fires
  _pump(app)
  assert v.verticalScrollBar().value() == 0, v.verticalScrollBar().value()


def exercise_navigate_into_collapsed_disclosure_expands_and_scrolls():
  app = _app()
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.conversation_search import ConversationSearch
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  v = ConversationView()
  v.resize(500, 300)
  v.show()
  for i in range(10):
    v.add_message(_text_message("filler %d" % i))
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="tool_use", data={
      "id": "t9", "name": "phil", "input": {"q": "params"}}),
    ContentBlock(type="tool_result", data={
      "tool_use_id": "t9", "is_error": False,
      "content": [ContentBlock(type="text", data={
        "text": "\n".join("filler row %d" % i for i in range(60))
                + "\nthe hidden needle\n"
                + "\n".join("tail row %d" % i for i in range(60))})]})])
  v.add_message(m)
  _pump(app)
  disc = v.widget().findChildren(ToolCallDisclosure)[0]
  assert not disc.body.isVisible()
  cs = ConversationSearch(v)
  cs.open()
  cs.bar._scope_boxes["tool"].setChecked(True)
  cs.bar._edit.setText("hidden needle")
  assert len(cs._matches) == 1
  bar0 = v.verticalScrollBar().value()
  cs.bar.next_requested.emit()
  _pump(app, ms=250)                 # expansion + deferred scroll + layout
  assert disc.body.isVisible()
  # The final scroll position points into the EXPANDED result view: the
  # match row is inside the viewport band, which is only possible after
  # the auto-height refresh grew the cell (stale ~1-line geometry would
  # leave the scrollbar near bar0).
  w, start, _n = cs._matches[0]
  cursor = QtGui.QTextCursor(w.document())
  cursor.setPosition(start)
  pt = w.viewport().mapTo(v.widget(), w.cursorRect(cursor).topLeft())
  bar = v.verticalScrollBar()
  assert bar.value() != bar0, (bar.value(), bar0)
  assert bar.value() <= pt.y() <= bar.value() + v.viewport().height(), \
      (bar.value(), pt.y(), v.viewport().height())


def exercise():
  exercise_bubble_searchable_cells_order_and_kinds()
  exercise_bubble_enumerates_every_thinking_cell()
  exercise_view_searchable_cells_flattens_bubbles_in_order()
  exercise_ensure_visible_disengages_follow_and_holds()
  exercise_add_message_respects_autofollow_hold()
  exercise_expand_recalcs_collapsed_result_height()
  exercise_expand_calls_refresh_inner_heights_synchronously()
  exercise_search_bar_signals_and_defaults()
  exercise_search_bar_count_states()
  exercise_search_bar_keys_from_all_children()
  exercise_controller_open_close_toggle()
  exercise_scan_default_scope_and_counts()
  exercise_highlights_applied_and_cleared()
  exercise_match_cap()
  exercise_navigation_steps_and_wraps()
  exercise_navigation_relocates_after_offset_shift()
  exercise_stale_deferred_scroll_noop_after_close()
  exercise_same_cell_current_vs_other_highlight()
  exercise_overlay_repositions_when_scrollbar_appears()
  exercise_scroll_to_top_and_bottom()
  exercise_find_next_and_previous_controller_api()
  exercise_close_preserves_scroll_position()
  exercise_bar_width_tracks_view()
  exercise_navigation_pin_engages_before_deferred_scroll()
  exercise_deferred_bottom_snap_rechecks_follow()
  exercise_navigate_into_collapsed_disclosure_expands_and_scrolls()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
