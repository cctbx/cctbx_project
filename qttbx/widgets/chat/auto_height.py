"""Make a QTextEdit-family widget grow to fit its content.

The chat view's outer QScrollArea (ConversationView) is the sole
scrollable surface — per-bubble text widgets that grow their own
scrollbars stack visually badly and force the user to scroll twice
to read a long message. set_auto_height(widget) reconfigures the
widget so its height tracks the document height as content changes
and as the widget's width changes (the available wrap width).
"""

from qttbx.qt import QtCore, QtWidgets


def set_auto_height(widget, min_lines=1):
  """Configure a QTextEdit / QPlainTextEdit / QTextBrowser to grow
  with its content instead of scrolling.

  Effects:
    - Both scrollbars disabled (the outer scroll area handles it).
    - Vertical size policy = Fixed, horizontal = Expanding.
    - Any existing maximumHeight constraint is lifted.
    - On document changes, the widget's height is set to the
      document's current rendered height (plus the frame chrome).
    - On viewport resize (wrap width changes), the height is
      recomputed because long lines may re-wrap.

  min_lines clamps the minimum rendered height to N text lines so an
  empty widget doesn't collapse to zero pixels mid-stream.
  """
  widget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
  widget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
  widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                       QtWidgets.QSizePolicy.Fixed)
  widget.setMaximumHeight(16777215)   # Qt's QWIDGETSIZE_MAX -- "no cap"

  def _refresh():
    doc = widget.document()
    # Wrap at the current viewport width so height tracks reality.
    viewport_w = widget.viewport().width()
    if viewport_w > 0:
      doc.setTextWidth(viewport_w)
    frame = 2 * widget.frameWidth()
    line_h = widget.fontMetrics().lineSpacing()
    # QPlainTextEdit uses QPlainTextDocumentLayout, whose
    # documentSize().height() returns the BLOCK count, not pixels --
    # so a 30-line plain-text document reports height 30.0 instead of
    # ~30 * lineSpacing. Sum block bounding rects for an accurate,
    # wrap-aware pixel height. QTextEdit and QTextBrowser use the
    # default QTextDocumentLayout which reports pixel heights directly.
    if isinstance(widget, QtWidgets.QPlainTextEdit):
      layout = doc.documentLayout()
      doc_h = 0
      block = doc.firstBlock()
      while block.isValid():
        doc_h += layout.blockBoundingRect(block).height()
        block = block.next()
      doc_h = int(doc_h)
    else:
      doc_h = int(doc.size().height())
    min_h = max(min_lines, 1) * line_h + frame
    widget.setFixedHeight(max(doc_h + frame + 2, min_h))

  widget.document().contentsChanged.connect(_refresh)

  # Hook the resize event so width-driven rewrap recomputes the height.
  # Wrapping the bound method (rather than subclassing) keeps the
  # helper usable on existing widget subclasses without forcing a
  # new class per text-widget kind.
  _orig_resize = widget.resizeEvent

  def _resize_event(event):
    _orig_resize(event)
    _refresh()

  widget.resizeEvent = _resize_event

  # Initial sizing so the widget isn't a 0-px sliver before the first
  # contentsChanged fires.
  _refresh()

  # Expose the refresh callable so callers that toggle the widget's
  # visibility (e.g. ToolCallDisclosure expanding a collapsed body)
  # can force a recalc once Qt has assigned real geometry. Content set
  # while the widget was hidden gets sized for viewport().width() == 0
  # otherwise, and the inner view stays collapsed even after show.
  widget._auto_height_refresh = _refresh
