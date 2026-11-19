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
  # And fits it whole: a view a few pixels short of its text -- the
  # top document margin was left out -- kept a one-line scroll range,
  # and a wheel notch, PageDown or Space over it shifted the text up a
  # line instead of scrolling the conversation.
  _assert_nothing_to_scroll(w)


def _assert_nothing_to_scroll(w):
  ranges = (w.verticalScrollBar().maximum(), w.horizontalScrollBar().maximum())
  assert ranges == (0, 0), (
    "%d lines / %d px left for the widget's own scrollbars" % ranges)


def exercise_viewport_and_scroll_range_follow_a_rewrap():
  """A height set from inside the viewport's resize reaches the viewport.

  QAbstractScrollArea delivers resizeEvent for its viewport's resize,
  from inside its own child layout (guarded against re-entry) or, on
  first show, before the viewport is visible. Either way the widget
  grew and the viewport did not follow: an argument preview that
  rewrapped from three lines to five at a narrower window showed the
  first three, with the rest behind a scroll range of two for a wheel
  notch to reveal; a card shown at its minimum width kept the whole
  preview but one line. Widened back, the widget shrinks and the
  viewport has to follow it down as well.
  """
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height

  def host_with(view):
    host = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(host)
    layout.addWidget(view)
    layout.addStretch(1)
    return host

  def settled(host, view):
    last = None
    for _ in range(50):
      app.processEvents()
      now = (view.height(), view.viewport().height())
      if now == last:
        break
      last = now
    assert view.viewport().height() == view.height() - 2 * view.frameWidth(), (
      "the viewport is %d px tall in a %d px widget"
      % (view.viewport().height(), view.height()))
    _assert_nothing_to_scroll(view)
    cursor = view.textCursor()
    cursor.movePosition(QtGui.QTextCursor.End)
    # bottom() is the last pixel row; row height() is the first outside.
    assert view.cursorRect(cursor).bottom() < view.viewport().height(), (
      "the last line is clipped", view.cursorRect(cursor).bottom(),
      view.viewport().height())

  # First show at a width the text needs several lines at.
  view = QtWidgets.QPlainTextEdit()
  view.setPlainText("ACGT" * 60)
  set_auto_height(view)
  host = host_with(view)
  host.resize(250, 600)
  host.show()
  settled(host, view)
  three_lines = view.height()
  # Then the squeeze: the text rewraps and the widget grows in the hook.
  host.resize(140, 600)
  settled(host, view)
  assert view.height() > three_lines, (view.height(), three_lines)
  # And back: fewer lines, a shorter widget, the viewport with it.
  host.resize(250, 600)
  settled(host, view)
  assert view.height() == three_lines, (view.height(), three_lines)
  host.close()


def exercise_height_follows_a_font_change():
  """A larger font relays out the document without a contentsChanged,
  and the view kept its height: a Ctrl+wheel zoom of a read-only
  QPlainTextEdit left three lines behind an inner scroll range that
  then swallowed plain wheel notches. The application's font reaches a
  view without one of its own the same way."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height

  def grown(change, cls):
    view = cls()
    view.setPlainText("ACGT" * 60)
    set_auto_height(view)
    host = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(host)
    layout.addWidget(view)
    layout.addStretch(1)
    host.resize(250, 600)
    host.show()
    app.processEvents()
    before = view.height()
    bigger = QtGui.QFont(view.font())
    bigger.setPointSize(max(24, bigger.pointSize() * 2))
    change(view, bigger)
    for _ in range(10):
      app.processEvents()
    assert view.height() > before, (view.height(), before)
    _assert_nothing_to_scroll(view)
    cursor = view.textCursor()
    cursor.movePosition(QtGui.QTextCursor.End)
    assert view.cursorRect(cursor).bottom() < view.viewport().height(), (
      "the last line is clipped", view.cursorRect(cursor).bottom(),
      view.viewport().height())
    host.close()

  # A tool result's QPlainTextEdit and a reply's QTextEdit alike.
  for cls in (QtWidgets.QPlainTextEdit, QtWidgets.QTextEdit):
    grown(lambda view, font: view.setFont(font), cls)
    saved = QtGui.QFont(app.font())
    try:
      grown(lambda view, font: app.setFont(font), cls)
    finally:
      app.setFont(saved)


def exercise_a_view_measured_hidden_reopens_with_nothing_to_scroll():
  """A folded view changed at an unchanged height shows with no inner
  scroll range.

  Qt's scrollbar pass for a hidden QPlainTextEdit leaves the whole
  document minus a line in range, and only a resize re-runs it: a
  letter-spacing change, or a same-height result, made while a tool
  row was folded reopened to a view that a wheel notch scrolled by
  three lines instead of the conversation."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height

  def spaced(view):
    font = QtGui.QFont(view.font())
    font.setLetterSpacing(QtGui.QFont.AbsoluteSpacing, 2)
    view.setFont(font)

  def rewritten(view):
    view.setPlainText("\n".join("LINE %d" % i for i in range(40)))

  for change in (spaced, rewritten):
    host = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(host)
    body = QtWidgets.QWidget(host)
    layout.addWidget(body)
    layout.addStretch(1)
    view = QtWidgets.QPlainTextEdit(body)
    QtWidgets.QVBoxLayout(body).addWidget(view)
    view.setPlainText("\n".join("line %d" % i for i in range(40)))
    set_auto_height(view)
    host.resize(400, 1200)
    host.show()
    app.processEvents()
    _assert_nothing_to_scroll(view)
    body.hide()
    app.processEvents()
    folded = view.height()
    change(view)
    app.processEvents()
    assert view.height() == folded, (change.__name__, view.height(), folded)
    body.show()
    app.processEvents()
    _assert_nothing_to_scroll(view)
    host.close()


def exercise_a_relaid_document_is_measured_shown_or_folded():
  """A wrap mode or a tab stop relays out the document as a font does,
  and the height follows: at once while shown, at the show while folded.

  Only the view's resize was hooked: a wrap mode set while shown left
  the height of the wrapped lines, and one set while folded opened at
  it."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  body = QtWidgets.QWidget(host)
  layout.addWidget(body)
  layout.addStretch(1)
  view = QtWidgets.QPlainTextEdit(body)
  QtWidgets.QVBoxLayout(body).addWidget(view)
  view.setPlainText("\n".join(
    "line %d\t" % i + "\tword " * 30 for i in range(40)))
  set_auto_height(view)
  host.resize(400, 1200)
  host.show()
  app.processEvents()
  wrapped = view.height()

  def settled():
    # A resize measures again; a height that moves was stale. Up and
    # back: one pixel more crosses a wrap step of the tabbed text on
    # Qt 6.8, which measures a different height for a different width.
    width = host.width()
    host.resize(width + 1, 1200)
    app.processEvents()
    host.resize(width, 1200)
    app.processEvents()
    return view.height()

  view.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
  app.processEvents()
  unwrapped = view.height()
  assert unwrapped < wrapped, (unwrapped, wrapped)
  assert settled() == unwrapped, (view.height(), unwrapped)
  body.hide()
  app.processEvents()
  view.setLineWrapMode(QtWidgets.QPlainTextEdit.WidgetWidth)
  assert view.height() == unwrapped, (view.height(), unwrapped)
  body.show()
  app.processEvents()
  opened = view.height()
  assert opened > unwrapped, (opened, unwrapped)
  assert settled() == opened, (view.height(), opened)
  view.setTabStopDistance(300.0)
  app.processEvents()
  tabbed = view.height()
  assert settled() == tabbed, (view.height(), tabbed)
  _assert_nothing_to_scroll(view)


def exercise_lines_appended_are_measured_shown_or_at_the_show():
  """Lines a result appends reach the height at once while shown and
  at the show while folded, and a folded view's text set whole opens
  with nothing to scroll.

  The layout says nothing for an append -- it speaks for a text set
  whole, a resize, a font -- so the document's own word is what
  measures one. A height set while folded leaves the viewport to the
  show: brought up folded, the viewport had the resize the show would
  have delivered taken from it, and opened with the whole document but
  a line in its scroll range."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  body = QtWidgets.QWidget(host)
  layout.addWidget(body)
  layout.addStretch(1)
  view = QtWidgets.QPlainTextEdit(body)
  QtWidgets.QVBoxLayout(body).addWidget(view)
  view.setPlainText("line 0")
  set_auto_height(view)
  host.resize(400, 1200)
  host.show()
  app.processEvents()
  one = view.height()
  for i in range(1, 8):
    view.appendPlainText("line %d" % i)
  app.processEvents()
  shown = view.height()
  assert shown > one, (shown, one)
  _assert_nothing_to_scroll(view)
  body.hide()
  app.processEvents()
  for i in range(8, 15):
    view.appendPlainText("line %d" % i)
  body.show()
  app.processEvents()
  assert view.height() > shown, (view.height(), shown)
  _assert_nothing_to_scroll(view)
  body.hide()
  app.processEvents()
  view.setPlainText("\n".join("other %d" % i for i in range(25)))
  body.show()
  app.processEvents()
  _assert_nothing_to_scroll(view)


def exercise_a_folded_view_is_measured_at_its_show_after_a_font_change():
  """A font change reaches a folded view, whose height nobody reads
  until it is shown: measuring it there laid its whole document out,
  and the application's font over fifty folded 1000-line results took
  0.27 s against 0.003 s. The show measures it."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  body = QtWidgets.QWidget(host)
  layout.addWidget(body)
  layout.addStretch(1)
  view = QtWidgets.QPlainTextEdit(body)
  QtWidgets.QVBoxLayout(body).addWidget(view)
  view.setPlainText("\n".join("line %d" % i for i in range(40)))
  set_auto_height(view)
  host.resize(400, 1200)
  host.show()
  app.processEvents()
  body.hide()
  app.processEvents()
  folded = view.height()
  saved = QtGui.QFont(app.font())
  bigger = QtGui.QFont(saved)
  bigger.setPointSize(max(24, saved.pointSize() * 2))
  try:
    app.setFont(bigger)
    app.processEvents()
    assert view.height() == folded, (view.height(), folded)
    body.show()
    app.processEvents()
    assert view.height() > folded, (view.height(), folded)
    _assert_nothing_to_scroll(view)
    cursor = view.textCursor()
    cursor.movePosition(QtGui.QTextCursor.End)
    assert view.cursorRect(cursor).bottom() < view.viewport().height(), (
      "the last line is clipped", view.cursorRect(cursor).bottom(),
      view.viewport().height())
  finally:
    app.setFont(saved)
  host.close()


def exercise_a_resize_measures_once():
  """A resize that leaves the height alone measures the document once.

  The layout's own signal, sent from the view's resize, is the measure
  for a new width, and the show of a folded row delivers the resize
  queued while it was folded: a second measure on either, from a hook
  on the resize beside the signal, was a walk of every block of a tool
  result for nothing."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  view = QtWidgets.QPlainTextEdit()
  view.setPlainText("\n".join("line %d" % i for i in range(40)))
  set_auto_height(view)
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  layout.addWidget(view)
  layout.addStretch(1)
  host.resize(400, 1200)
  host.show()
  app.processEvents()
  calls = []
  fixed = view.setFixedHeight
  view.setFixedHeight = lambda h: (calls.append(h), fixed(h))
  host.resize(350, 1200)
  app.processEvents()
  assert calls == [view.height()], (calls, view.height())
  host.close()


def exercise_the_first_show_lays_the_document_out_once():
  """A first show lays the document out at the viewport's width once.

  QPlainTextEdit lays its document out again on any resize whose old
  width differs from the new one. The viewport's catch-up resize, sent
  for the one its show left pending, names the width it has as the old
  one, so the document the measure just walked stands."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  view = QtWidgets.QPlainTextEdit()
  view.setPlainText("\n".join("line %d " % i + "word " * 30
                              for i in range(40)))
  set_auto_height(view)
  host = QtWidgets.QWidget()
  layout = QtWidgets.QVBoxLayout(host)
  layout.addWidget(view)
  layout.addStretch(1)
  host.resize(400, 3000)
  before = view.height()
  relayouts = []
  view.document().documentLayout().update.connect(
    lambda *args: relayouts.append(view.width()))
  host.show()
  assert view.height() != before, (view.height(), before)
  assert len(relayouts) == 1, relayouts
  host.close()


def exercise_a_content_change_does_not_lay_the_document_out_again():
  """Measuring after a content change leaves the document's layout to
  the view.

  The view wraps its document at the viewport's width from its own
  resize event, ahead of the hook, and QPlainTextEdit's layout wraps
  at a width of its own that the document's does not reach. Setting
  the document's width from the measure laid the whole document out
  again on every change: an append to a 200-paragraph markdown view
  took a second."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  app = _qapp()
  from qttbx.qt import QtGui, QtWidgets
  from qttbx.widgets.chat.auto_height import set_auto_height
  calls = []

  def counting(name, original):
    def method(self, *args):
      calls.append(name)
      return original(self, *args)
    return method

  originals = {}
  for name in ("setTextWidth", "setPageSize"):
    originals[name] = getattr(QtGui.QTextDocument, name)
    setattr(QtGui.QTextDocument, name, counting(name, originals[name]))
  try:
    for cls in (QtWidgets.QPlainTextEdit, QtWidgets.QTextEdit):
      view = cls()
      set_auto_height(view)
      host = QtWidgets.QWidget()
      layout = QtWidgets.QVBoxLayout(host)
      layout.addWidget(view)
      layout.addStretch(1)
      host.resize(250, 600)
      host.show()
      app.processEvents()
      del calls[:]
      view.setPlainText("ACGT" * 60)
      app.processEvents()
      assert calls == [], (cls.__name__, calls)
      host.close()
  finally:
    for name, original in originals.items():
      setattr(QtGui.QTextDocument, name, original)


def exercise():
  exercise_helper_disables_both_scrollbars()
  exercise_height_grows_with_content()
  exercise_viewport_and_scroll_range_follow_a_rewrap()
  exercise_height_follows_a_font_change()
  exercise_a_view_measured_hidden_reopens_with_nothing_to_scroll()
  exercise_a_folded_view_is_measured_at_its_show_after_a_font_change()
  exercise_a_relaid_document_is_measured_shown_or_folded()
  exercise_lines_appended_are_measured_shown_or_at_the_show()
  exercise_a_resize_measures_once()
  exercise_the_first_show_lays_the_document_out_once()
  exercise_a_content_change_does_not_lay_the_document_out_again()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
