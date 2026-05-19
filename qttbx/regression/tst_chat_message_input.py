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


def exercise():
  exercise_send_signal_carries_text_and_empty_attachments()
  exercise_empty_send_is_no_op()
  exercise_set_busy_toggles_to_stop()
  exercise_attach_bytes_appears_in_send_payload()
  exercise_unsupported_mime_is_dropped_with_warning()
  exercise_oversized_image_is_resampled()
  exercise_save_chat_button_emits_signal()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
