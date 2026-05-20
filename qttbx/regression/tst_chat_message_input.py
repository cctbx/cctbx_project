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
  w = MessageInput()
  received = []
  w.send.connect(lambda msg, atts: received.append(msg))
  w.click_send()
  assert received == [], received


def exercise_set_busy_toggles_to_stop():
  from qttbx.widgets.chat.message_input import MessageInput
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
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


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
