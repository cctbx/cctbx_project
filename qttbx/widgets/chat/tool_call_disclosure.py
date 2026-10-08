"""Collapsible rows used inside MessageBubble: tool calls and thinking.

``DisclosureRow`` is the fold: a single-line ``▸ <label>`` header that
click-expands to ``▾`` over a body hidden by default (or shown, for a
thinking cell). ``ToolCallDisclosure`` fills the body with pretty-printed
JSON args on top and a plain-text result below; status text and colour are
settable so the runner can transition running -> finished/failed/cancelled
without restructuring the widget. The thinking cell in message_bubble is
the other subclass.
"""

import json

from qttbx.qt import QtCore, QtGui, QtWidgets

from qttbx.widgets.chat.eliding import ElidingToolButton


# Header state styling. The (status) parenthetical in the header text
# already communicates state literally; the style here is just a hint.
# 'running' is italic (no colour override) so legibility doesn't depend
# on hue; 'error' uses a red that reads on both light and dark themes;
# 'muted' uses the theme's mid colour so it stays in palette. 'cancelled'
# is the terminal state for a tool aborted by Stop: a muted+italic look
# that is distinct from both 'default'/finished (no style) and 'error'
# (red), so a cancelled call doesn't read as a successful or failed one.
_COLORS = {
  None:        "",                       # default text color
  "default":   "",
  "running":   "font-style: italic;",
  "error":     "color: #c0392b;",
  "muted":     "color: palette(mid);",
  "cancelled": "color: palette(mid); font-style: italic;",
}

# Longest status rendered in the header, which carries a one-line
# '<name> (<status>)' summary. A status with newlines in it, or a whole tool
# result pasted into it, is a caller error; this keeps the rendered string
# short and single-line and hangs the rest on the tooltip.
#
# This is a READABILITY bound, not the layout's protection: a character budget
# cannot bound a pixel width -- 60 chars of a proportional font is still
# ~880 px -- and the tool name it sits next to is not bounded at all. The
# header not flooring the view is ElidingToolButton's job.
_MAX_HEADER_STATUS = 60


def set_transparent_background(view):
  """Make a text view paint no background of its own, theme-switch safe.

  The obvious recipe -- ``setStyleSheet("background: transparent;")`` --
  has a Qt 5.15 trap: once a stylesheet-styled view has been polished,
  it can keep the palette the stylesheet engine resolved for it, so an
  OS light/dark switch afterwards leaves its text at the OLD theme's
  colour. Setting the ``Base`` role transparent on the widget's own
  palette pins only that role; ``Text`` keeps inheriting, so the view
  follows the app palette like an unstyled widget.

  Parameters
  ----------
  view : QtWidgets.QAbstractScrollArea
      A QPlainTextEdit / QTextEdit style view.
  """
  palette = view.palette()
  palette.setColor(QtGui.QPalette.Base, QtCore.Qt.transparent)
  view.setPalette(palette)
  view.viewport().setAutoFillBackground(False)


class DisclosureRow(QtWidgets.QFrame):
  """A ``▸ label`` header over a foldable body.

  Subclasses populate ``self.body_layout`` and implement
  ``_inner_views()`` (the auto-height views to refresh on open) and
  ``_header_label()`` (the text after the arrow, plus a tooltip or
  ``None``). The base constructor already renders the header
  (``_refresh_header()`` -> ``_header_label()``), so a subclass must set
  whatever state those read -- including any ``_refresh_header``
  override's, e.g. the tool row's ``_color`` -- BEFORE it calls
  ``super().__init__()``.

  Parameters
  ----------
  parent : QtWidgets.QWidget, optional
      Parent widget.
  expanded : bool, optional
      Start with the body shown. Tool rows start folded.
  """

  def __init__(self, parent=None, expanded=False):
    super().__init__(parent)
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(2)

    # Eliding: the header renders text nobody bounds (an MCP server names
    # its own tools, and a caller can put anything in a status). A plain
    # QToolButton reports that whole string as its minimumSizeHint and
    # never elides, so it floors the bubble's minimum width and with it the
    # whole ConversationView's -- the bubbles stop tracking the window and
    # a horizontal scrollbar appears. The base sets no colour, so the
    # header takes the theme's text colour; a subclass may style it (the
    # tool row colours its header per status).
    self.header_button = ElidingToolButton(self)
    self.header_button.setAutoRaise(True)
    self.header_button.setCheckable(True)
    self.header_button.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
    self.header_button.setCursor(QtCore.Qt.PointingHandCursor)
    self.header_button.clicked.connect(self._on_toggled)
    layout.addWidget(self.header_button)

    # The body indents its content under the header so the fold reads as
    # one unit.
    self.body = QtWidgets.QWidget(self)
    self.body_layout = QtWidgets.QVBoxLayout(self.body)
    self.body_layout.setContentsMargins(16, 2, 0, 4)
    self.body_layout.setSpacing(4)
    layout.addWidget(self.body)
    self.header_button.setChecked(bool(expanded))
    self._sync_body_to_header()

  # ---- hooks --------------------------------------------------------------

  def _inner_views(self):
    """The auto-height views in the body, refreshed when the row opens."""
    raise NotImplementedError

  def _header_label(self):
    """``(text, tooltip)`` shown after the arrow; tooltip ``None`` leaves
    it to the eliding button (full text only when elided)."""
    raise NotImplementedError

  # ---- public API ---------------------------------------------------------

  def is_expanded(self):
    """True when the body is shown."""
    return self.header_button.isChecked()

  def set_expanded(self, on):
    """Show (``True``) or fold (``False``) the body programmatically.

    On open, the body's views are re-measured before this returns:
    ``ConversationView.set_thinking_expanded`` re-anchors the viewport in
    the same call and must read real geometry (the click path defers the
    refresh a tick; search's reveal defers its own scroll). Showing the
    body resizes a view that opens at a new width, and auto_height's
    resize hook re-measures it; the explicit refresh covers a view that
    opens at the width it already had (no resize event), e.g. after a
    font change made while it was folded. The refresh runs BEFORE the
    body's and the row's layouts are activated, so the height it sets is
    what the row's size hint carries out of this call. No-op when already
    in the requested state.
    """
    on = bool(on)
    if self.is_expanded() == on:
      return
    self.header_button.setChecked(on)
    self._sync_body_to_header()
    if on:
      self._refresh_inner_heights()
      self.body_layout.activate()
      self.layout().activate()

  def ensure_revealed(self):
    """Reveal hidden searchable content (duck-typed protocol).

    ConversationSearch walks a match's ancestors and calls this on any
    that expose it, so the controller needs no knowledge of this class
    or of whether the body is currently collapsed.
    """
    self.set_expanded(True)

  # ---- internals ----------------------------------------------------------

  def _sync_body_to_header(self):
    """Show/hide the body to match the header's checked state."""
    self.body.setVisible(self.header_button.isChecked())
    self._refresh_header()

  def _on_toggled(self):
    self._sync_body_to_header()
    if self.header_button.isChecked():
      # A view that opens at a new width is re-measured by auto_height's
      # resize hook when the layout pass assigns that width; one that
      # opens at the width it already had gets no resize event, so
      # re-measure it after that pass -- deferred a tick, since nothing
      # reads geometry in the click's own call.
      QtCore.QTimer.singleShot(0, self._refresh_inner_heights)

  def _refresh_inner_heights(self):
    for view in self._inner_views():
      # The deferred singleShot may fire after the widget's C++ object
      # has been destroyed (parent garbage-collected while the callback
      # was still in Qt's event queue). Guard so a stale callback no-ops.
      try:
        if not view.isHidden():
          view._auto_height_refresh()
      except RuntimeError:
        # 'Internal C++ object already deleted' -- widget is gone.
        return

  def _refresh_header(self):
    arrow = "▾" if self.header_button.isChecked() else "▸"
    text, tooltip = self._header_label()
    self.header_button.set_full_text("%s %s" % (arrow, text), tooltip=tooltip)


class ToolCallDisclosure(DisclosureRow):
  """Collapsible disclosure row for a single tool call.

  Parameters
  ----------
  name : str
      Tool name shown in the header.
  status : str
      Initial status text shown in the header parenthetical. A status
      of ``'running'`` starts the header in the italic running style.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  def __init__(self, name, status, parent=None):
    self._name = name
    self._status = status
    self._color = "running" if status == "running" else None
    super().__init__(parent)

    # Args + result views grow with their content; the outer
    # ConversationView is the sole scroller per the chat UI redesign.
    # No setMaximumHeight cap -- a 10K-line tool result produces a
    # 10K-line bubble that the outer view scrolls through.
    from qttbx.widgets.chat.auto_height import set_auto_height

    self.args_view = QtWidgets.QPlainTextEdit(self.body)
    self.args_view.setReadOnly(True)
    self.args_view.setFrameStyle(QtWidgets.QFrame.NoFrame)
    # Monospace for JSON args; no explicit color so the text follows
    # the active theme's palette (a hardcoded grey like '#555' rendered
    # as near-invisible grey-on-dark under dark themes).
    # Resolve the platform's actual fixed-pitch family via QFontDatabase
    # rather than asking for 'monospace' in the stylesheet -- Qt has no
    # font literally named "Monospace" on macOS / Windows, so the CSS
    # generic forces a one-time ~50 ms alias scan and prints a
    # qt.qpa.fonts warning. systemFont(FixedFont) returns Menlo on
    # macOS, Consolas on Windows, etc.
    self.args_view.setFont(QtGui.QFontDatabase.systemFont(
      QtGui.QFontDatabase.FixedFont))
    set_transparent_background(self.args_view)
    set_auto_height(self.args_view)
    self.args_view.hide()

    self.result_view = QtWidgets.QPlainTextEdit(self.body)
    self.result_view.setReadOnly(True)
    self.result_view.setFrameStyle(QtWidgets.QFrame.NoFrame)
    set_transparent_background(self.result_view)
    set_auto_height(self.result_view)
    self.result_view.hide()

    self.body_layout.addWidget(self.args_view)
    self.body_layout.addWidget(self.result_view)

  # ---- public API ---------------------------------------------------------

  def set_status(self, status, color=None):
    """Update the header status text.

    Parameters
    ----------
    status : str
        New status text shown in the header parenthetical. Coerced with
        ``str``, so an exception object renders rather than raising. Keep it
        short: the header is a one-line summary, so only the first line
        survives and anything past ``_MAX_HEADER_STATUS`` characters moves to
        the tooltip. Bulk text belongs in ``set_result``.
    color : str or None, optional
        One of ``None``, ``'default'``, ``'running'``, ``'error'``,
        ``'muted'``, or ``'cancelled'``. ``None`` preserves the current color.
    """
    self._status = status
    if color is not None:
      self._color = color
    self._refresh_header()

  def set_args(self, args):
    """Set the tool arguments, rendered as pretty-printed JSON.

    Parameters
    ----------
    args : dict or None
        Tool arguments (typically a dict for MCP tools). ``None`` clears
        and hides the args view.
    """
    if args is None:
      self.args_view.setPlainText("")
      self.args_view.hide()
      return
    try:
      text = json.dumps(args, indent=2, sort_keys=True)
    except (TypeError, ValueError):
      text = repr(args)
    self.args_view.setPlainText(text)
    self.args_view.show()

  def set_result(self, text):
    """Set the tool result text shown below the args.

    Parameters
    ----------
    text : str
        Plain-text result. An empty value clears and hides the result
        view. Coerced with ``str``, so an exception object or any other
        payload renders instead of raising out of the handler that is
        delivering the tool result.
    """
    if not text:
      self.result_view.setPlainText("")
      self.result_view.hide()
      return
    self.result_view.setPlainText(str(text))
    self.result_view.show()

  def searchable_cells(self):
    """This row's searchable text: args and result views, kind ``"tool"``.

    Both views are reported even while the body is collapsed -- hidden
    tool text is searchable, and navigation reveals it via
    ``ensure_revealed``.
    """
    return [("tool", self.args_view), ("tool", self.result_view)]

  def is_running(self):
    """Return True while the call is still in its initial ``running`` state.

    The predicate the turn-cancel sweep uses to find tool cells that never
    reached a terminal state (finished / failed / cancelled) -- their result
    will never arrive, so they would otherwise stay stuck spinning.
    """
    return self._status == "running"

  # ---- hooks / internals --------------------------------------------------

  def _inner_views(self):
    return (self.args_view, self.result_view)

  def _status_text(self):
    """Return the status as a string.

    Coerced rather than assumed: the header used to build its text with
    ``'%s'``, so a caller could pass an exception object -- or anything else
    -- and see it rendered. Reading the raw value to clamp it would raise
    ``AttributeError`` on everything that is not a str.
    """
    return "" if self._status is None else str(self._status)

  def _header_status(self):
    """Return the status reduced to one short line for the header.

    Returns ``(text, clamped)``; ``clamped`` is True when anything was
    dropped, so the caller can hang the full status on the tooltip.
    """
    full = self._status_text()
    first = full.split("\n", 1)[0].strip()
    if len(first) > _MAX_HEADER_STATUS:
      return first[:_MAX_HEADER_STATUS - 1].rstrip() + "…", True
    return first, first != full

  def _header_label(self):
    status, clamped = self._header_status()
    # A clamped status keeps the whole of itself on the tooltip. Otherwise
    # leave the tooltip to the button, which shows the full header text only
    # when it had to elide it -- so an ordinary 'finished' row at a normal
    # width sprouts no redundant hover label.
    tooltip = None
    if clamped:
      tooltip = "%s (%s)" % (self._name, self._status_text())
    return "%s (%s)" % (self._name, status), tooltip

  def _refresh_header(self):
    super()._refresh_header()
    self.header_button.setStyleSheet(_COLORS.get(self._color, ""))
