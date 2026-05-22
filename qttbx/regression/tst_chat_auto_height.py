"""auto_height.set_auto_height(): a QTextEdit-family widget grows with
its content instead of scrolling. Both scrollbars are disabled and
the widget's height tracks document().size().height() + frame chrome.
Spec: ConversationView is the sole scroller; per-bubble scrollbars
stack visually badly."""

import os

from libtbx.utils import format_cpu_times


def _qapp():
  from qttbx.qt import QtWidgets
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
  init_default_app_font(app)
  return app


def exercise_helper_disables_both_scrollbars():
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.qt import QtCore, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  w = QtWidgets.QPlainTextEdit()
  set_auto_height(w)
  assert w.verticalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff
  assert w.horizontalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff


def exercise_height_grows_with_content():
  """A widget configured by set_auto_height(): setPlainText with more
  content must produce a larger fixed height than less content AND
  the height must actually fit the rendered text (rough check: at
  least ~N - 2 lines for an N-line plain-text document).

  Regression: QPlainTextDocumentLayout reports documentSize().height()
  as the BLOCK count, not pixels -- the long_h > short_h check passed
  even when long_h was off by ~16x (84 px for 80 lines that need
  ~1280 px). The 'fits the content' assertion below catches that.

  The widget must be shown so its viewport has a real width --
  without it block layout treats the document as zero-width."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  w = QtWidgets.QPlainTextEdit()
  w.resize(400, 800)
  w.show()
  set_auto_height(w)
  w.setPlainText("one short line")
  app.processEvents()
  short_h = w.height()
  n_lines = 80
  many = "\n".join("line %d" % i for i in range(n_lines))
  w.setPlainText(many)
  app.processEvents()
  long_h = w.height()
  assert long_h > short_h, (short_h, long_h)
  line_h = w.fontMetrics().lineSpacing()
  needed = (n_lines - 2) * line_h
  assert long_h >= needed, (long_h, needed, line_h, n_lines)


def exercise():
  exercise_helper_disables_both_scrollbars()
  exercise_height_grows_with_content()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
