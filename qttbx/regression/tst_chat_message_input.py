"""MessageInput widget test. Covers: send signal payload, button toggle,
stop signal, and that pressing send with empty text is a no-op."""

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


def exercise_send_signal_carries_text_and_empty_attachments():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  w.set_text("hello")
  received = []
  w.send.connect(lambda msg, atts: received.append((msg, atts)))
  w.click_send()
  assert received == [("hello", [])], received
  # After send, the input is cleared.
  assert w.text() == ""


def exercise_empty_send_is_no_op():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  received = []
  w.send.connect(lambda msg, atts: received.append(msg))
  w.click_send()
  assert received == [], received


def exercise_set_busy_toggles_to_stop():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  stops = []
  w.stop.connect(lambda: stops.append(True))
  w.set_busy(True)
  # Now the action button emits stop, not send.
  w.click_send()
  assert stops == [True]
  w.set_busy(False)


def exercise_attach_bytes_appears_in_send_payload():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  w.attach_bytes(b"\x89PNG\r\n\x1a\nfake", "image/png",
                 filename="img.png")
  received = []
  w.send.connect(lambda msg, atts: received.append((msg, atts)))
  w.set_text("look at this")
  w.click_send()
  assert len(received) == 1
  text, atts = received[0]
  assert text == "look at this"
  assert len(atts) == 1
  assert atts[0]["mime"] == "image/png"
  assert atts[0]["filename"] == "img.png"
  assert atts[0]["bytes"].startswith(b"\x89PNG")


def exercise_unsupported_mime_is_dropped_with_warning():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  warned = []
  w.attachment_rejected.connect(lambda msg: warned.append(msg))
  w.attach_bytes(b"hello", "text/plain", filename="x.txt")
  assert warned, "expected rejection signal for text/plain"
  assert w.attachment_count() == 0


def exercise_oversized_image_is_resampled():
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  img = QtGui.QImage(2000, 2000, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(50, 100, 150))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  big = bytes(buf.data())
  w = MessageInput()
  # Force a tiny cap so resampling kicks in.
  w._max_image_bytes = max(1024, len(big) // 4)
  w.attach_bytes(big, "image/png", filename="big.png")
  assert w.attachment_count() == 1
  atts = w._attachments
  assert len(atts[0]["bytes"]) <= w._max_image_bytes


def exercise_save_chat_button_emits_signal():
  """The 'Save chat' button in the lower-left of the button row emits
  the parameterless save_chat signal so the chat window can prompt for
  a path and write the markdown export. No coupling to a conversation
  here -- MessageInput stays UI-only."""
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  fired = []
  w.save_chat.connect(lambda: fired.append(True))
  w._save_chat_btn.click()
  assert fired == [True], fired


def exercise_placeholder_set_and_reset():
  """set_placeholder swaps the edit's placeholderText; reset returns
  to MessageInput.DEFAULT_PLACEHOLDER. Used by ChatWindow to cycle
  'Thinking...' verbs while a turn is in flight."""
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  assert w._edit.placeholderText() == MessageInput.DEFAULT_PLACEHOLDER
  w.set_placeholder("Refining...")
  assert w._edit.placeholderText() == "Refining..."
  w.reset_placeholder()
  assert w._edit.placeholderText() == MessageInput.DEFAULT_PLACEHOLDER


def exercise_placeholder_dim_flag_controls_palette_role():
  """set_placeholder(..., dim=False) renders the placeholder at the
  regular Text colour (full contrast) rather than the theme's dim
  PlaceholderText colour. reset_placeholder restores the dim. Used
  by ChatWindow so 'Thinking...' verbs aren't visually muted."""
  from qttbx.qt import QtGui
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  text_color = w._edit.palette().color(QtGui.QPalette.Text)
  dim_color = w._dim_placeholder_color
  # Idle: dim placeholder.
  assert w._edit.palette().color(QtGui.QPalette.PlaceholderText) == \
    dim_color
  # Thinking: full contrast.
  w.set_placeholder("Refining...", dim=False)
  assert w._edit.palette().color(QtGui.QPalette.PlaceholderText) == \
    text_color
  # Reset returns to dim.
  w.reset_placeholder()
  assert w._edit.palette().color(QtGui.QPalette.PlaceholderText) == \
    dim_color


def exercise_auto_approve_button_is_checkable_and_emits_signal():
  """The centred 'Auto-approve' button is a checkable QPushButton.
  Clicking it flips the checked state and emits
  auto_approve_changed(checked); the label flips so the user can tell
  the toggle is active. No colour override on either state so the
  text stays readable in both light and dark themes."""
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  w = MessageInput()
  assert w._auto_approve_btn.isCheckable()
  assert not w._auto_approve_btn.isChecked()
  assert w._auto_approve_btn.text() == "Auto-approve"
  changes = []
  w.auto_approve_changed.connect(lambda v: changes.append(v))
  w._auto_approve_btn.click()
  assert changes == [True], changes
  assert w._auto_approve_btn.isChecked()
  assert "ON" in w._auto_approve_btn.text()
  # No colour rule on the button: a hardcoded colour reads poorly on
  # one theme or the other, and Qt's native checked-state rendering
  # is enough of a visual cue alongside the label flip.
  assert "color:" not in (
    w._auto_approve_btn.styleSheet() or "").replace(" ", "").lower()
  w._auto_approve_btn.click()
  assert changes == [True, False], changes
  assert not w._auto_approve_btn.isChecked()
  assert w._auto_approve_btn.text() == "Auto-approve"


def _press(w, key):
  """Send a plain key press through MessageInput's event filter, the way
  Qt delivers it. If the filter does not consume it, let the edit handle
  it (so a non-history Up/Down still moves the cursor)."""
  from qttbx.qt import QtCore, QtGui
  ev = QtGui.QKeyEvent(QtCore.QEvent.KeyPress, key, QtCore.Qt.NoModifier)
  if not w.eventFilter(w._edit, ev):
    w._edit.keyPressEvent(ev)


def _new_input():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  return MessageInput()


def exercise_up_arrow_recalls_previous_inputs():
  """Up at the top line walks back through the supplied history, newest
  first, and stops at the oldest entry."""
  from qttbx.qt import QtCore
  w = _new_input()
  w.set_history(["first", "second", "third"])
  _press(w, QtCore.Qt.Key_Up)
  assert w.text() == "third", w.text()
  _press(w, QtCore.Qt.Key_Up)
  assert w.text() == "second", w.text()
  _press(w, QtCore.Qt.Key_Up)
  assert w.text() == "first", w.text()
  _press(w, QtCore.Qt.Key_Up)              # already oldest -> stays
  assert w.text() == "first", w.text()


def exercise_down_arrow_walks_forward_and_restores_draft():
  """Down walks toward newer entries and, past the newest, restores the
  in-progress draft that was being typed when navigation began."""
  from qttbx.qt import QtCore
  w = _new_input()
  w.set_history(["a", "b", "c"])
  w.set_text("draft")
  _press(w, QtCore.Qt.Key_Up)              # save draft, -> newest "c"
  assert w.text() == "c", w.text()
  _press(w, QtCore.Qt.Key_Up)              # -> "b"
  assert w.text() == "b", w.text()
  _press(w, QtCore.Qt.Key_Down)            # -> "c"
  assert w.text() == "c", w.text()
  _press(w, QtCore.Qt.Key_Down)            # past newest -> restore draft
  assert w.text() == "draft", w.text()
  _press(w, QtCore.Qt.Key_Down)            # not navigating -> no-op
  assert w.text() == "draft", w.text()


def exercise_sent_message_is_appended_to_history():
  """A sent message becomes the newest recall entry immediately."""
  from qttbx.qt import QtCore
  w = _new_input()
  w.set_history(["earlier"])
  w.set_text("hello world")
  w.click_send()
  assert w.text() == ""                    # cleared on send
  _press(w, QtCore.Qt.Key_Up)
  assert w.text() == "hello world", w.text()
  _press(w, QtCore.Qt.Key_Up)
  assert w.text() == "earlier", w.text()


def exercise_up_arrow_below_first_line_moves_cursor_not_history():
  """In a multi-line draft, Up only recalls history when the cursor is
  on the first line; elsewhere it is a normal cursor move."""
  from qttbx.qt import QtCore, QtGui
  w = _new_input()
  w.set_history(["old"])
  w.set_text("line1\nline2")
  cur = w._edit.textCursor()
  cur.movePosition(QtGui.QTextCursor.End)   # cursor on last line (line2)
  w._edit.setTextCursor(cur)
  _press(w, QtCore.Qt.Key_Up)               # not first line -> no recall
  assert w.text() == "line1\nline2", w.text()
  cur = w._edit.textCursor()
  cur.movePosition(QtGui.QTextCursor.Start)  # cursor on first line
  w._edit.setTextCursor(cur)
  _press(w, QtCore.Qt.Key_Up)               # first line -> recall
  assert w.text() == "old", w.text()


def exercise_set_history_swaps_recall_list_and_resets_navigation():
  """Switching conversations replaces the recall list AND resets the
  navigation pointer and draft. Uses a longer new list and a fresh draft
  so a stale _history_index would land on the wrong entry, and a stale
  _history_draft would resurface foreign text -- a non-reset
  implementation fails here rather than passing vacuously."""
  from qttbx.qt import QtCore
  w = _new_input()
  w.set_history(["a1", "a2"])
  w.set_text("draftA")
  _press(w, QtCore.Qt.Key_Up)               # save draftA, index 2->1 -> "a2"
  assert w.text() == "a2", w.text()
  # Switch to a longer conversation's history with a new draft.
  w.set_history(["n1", "n2", "n3"])
  w.set_text("draftB")
  _press(w, QtCore.Qt.Key_Up)               # must start at the newest, "n3"
  assert w.text() == "n3", w.text()
  _press(w, QtCore.Qt.Key_Down)             # past newest -> restore draftB
  assert w.text() == "draftB", w.text()


def exercise_set_history_drops_blank_entries():
  """set_history honors its contract: empty and whitespace-only entries
  are dropped so a blank, invisible line never surfaces during recall."""
  w = _new_input()
  w.set_history(["", "  ", "\n", "real"])
  assert w._history == ["real"], w._history


def exercise():
  exercise_send_signal_carries_text_and_empty_attachments()
  exercise_empty_send_is_no_op()
  exercise_set_busy_toggles_to_stop()
  exercise_attach_bytes_appears_in_send_payload()
  exercise_unsupported_mime_is_dropped_with_warning()
  exercise_oversized_image_is_resampled()
  exercise_save_chat_button_emits_signal()
  exercise_auto_approve_button_is_checkable_and_emits_signal()
  exercise_placeholder_set_and_reset()
  exercise_placeholder_dim_flag_controls_palette_role()
  exercise_up_arrow_recalls_previous_inputs()
  exercise_down_arrow_walks_forward_and_restores_draft()
  exercise_sent_message_is_appended_to_history()
  exercise_up_arrow_below_first_line_moves_cursor_not_history()
  exercise_set_history_swaps_recall_list_and_resets_navigation()
  exercise_set_history_drops_blank_entries()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
