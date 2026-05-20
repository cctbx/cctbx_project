"""MarkdownView widget test. Smoke-only: construct, set markdown, assert the
underlying QTextDocument is non-empty."""

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


def exercise_set_and_append_markdown():
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  v = MarkdownView()
  v.set_markdown("**hello**")
  assert "hello" in v.toPlainText()
  v.append_markdown("\n\nworld")
  text = v.toPlainText()
  assert "hello" in text and "world" in text


def exercise_clear():
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  v = MarkdownView()
  v.set_markdown("x")
  v.clear()
  assert v.toPlainText().strip() == ""


def exercise_auto_height_no_scrollbar():
  """After the chat UI redesign, MarkdownView has no scrollbar of its
  own (auto_height). The outer ConversationView is the sole scroller.
  Pins the contract so a regression that re-enables the widget's own
  scrollbar surfaces here instead of as a stacked-scrolls UX bug."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  v = MarkdownView()
  assert v.verticalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff
  assert v.horizontalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff


def exercise():
  exercise_set_and_append_markdown()
  exercise_clear()
  exercise_auto_height_no_scrollbar()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
