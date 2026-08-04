"""Single-line collapsible row used inside MessageBubble for tool calls.

Used for tool calls and tool-approval prompts. Header:
``▸ tool_name (status)`` -- click to expand to ``▾``. Body (hidden by
default): pretty-printed JSON args on top, plain-text result below.
Status text and color are settable so the runner can transition
running → finished/failed/cancelled without restructuring the widget.
"""

import json

from qttbx.qt import QtCore, QtWidgets

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


class ToolCallDisclosure(QtWidgets.QFrame):
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
    super().__init__(parent)
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self._name = name
    self._status = status
    self._color = "running" if status == "running" else None

    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(2)

    # Eliding: the header renders a tool name and a status, neither of them
    # bounded (an MCP server names its own tools, and a caller can put
    # anything in a status). A plain QToolButton reports that whole string as
    # its minimumSizeHint and never elides, so it floors the bubble's minimum
    # width and with it the whole ConversationView's -- the bubbles stop
    # tracking the window and a horizontal scrollbar appears.
    self.header_button = ElidingToolButton(self)
    self.header_button.setAutoRaise(True)
    self.header_button.setCheckable(True)
    self.header_button.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
    self.header_button.setCursor(QtCore.Qt.PointingHandCursor)
    self.header_button.clicked.connect(self._on_toggled)
    layout.addWidget(self.header_button)

    self.body = QtWidgets.QWidget(self)
    body_layout = QtWidgets.QVBoxLayout(self.body)
    body_layout.setContentsMargins(16, 2, 0, 4)
    body_layout.setSpacing(4)

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
    from qttbx.qt import QtGui
    self.args_view.setFont(QtGui.QFontDatabase.systemFont(
      QtGui.QFontDatabase.FixedFont))
    self.args_view.setStyleSheet("background: transparent;")
    set_auto_height(self.args_view)
    self.args_view.hide()

    self.result_view = QtWidgets.QPlainTextEdit(self.body)
    self.result_view.setReadOnly(True)
    self.result_view.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self.result_view.setStyleSheet("background: transparent;")
    set_auto_height(self.result_view)
    self.result_view.hide()

    body_layout.addWidget(self.args_view)
    body_layout.addWidget(self.result_view)
    layout.addWidget(self.body)
    self.body.hide()

    self._refresh_header()

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

  def expand(self):
    """Expand the body programmatically (search navigation).

    A collapsed body's inner views were auto-height-sized while hidden
    (``viewport().width() == 0`` -> ~1 line tall), and the manual-click
    path defers the recalc one event-loop tick. A caller about to
    compute geometry (scroll-to-match) needs it NOW: show the body,
    force the layout chain so the views get real widths, then run the
    same refresh the deferred path uses. No-op when already expanded.
    """
    if self.header_button.isChecked() and self.body.isVisible():
      return
    self.header_button.setChecked(True)
    self._sync_body_to_header()
    self.body.layout().activate()
    if self.layout() is not None:
      self.layout().activate()
    self._refresh_inner_heights()

  def ensure_revealed(self):
    """Reveal hidden searchable content (duck-typed protocol).

    ConversationSearch walks a match's ancestors and calls this on any
    that expose it, so the controller needs no knowledge of this class
    or of whether the body is currently collapsed.
    """
    self.expand()

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

  # ---- internals ----------------------------------------------------------

  def _sync_body_to_header(self):
    """Show/hide the body to match the header's checked state."""
    self.body.setVisible(self.header_button.isChecked())
    self._refresh_header()

  def _on_toggled(self):
    self._sync_body_to_header()
    if self.header_button.isChecked():
      # Args/result heights were computed when the body was hidden
      # (viewport().width() == 0), so the inner views cached a tiny
      # single-line height. Force a recalc after Qt has assigned real
      # geometry to the now-visible body. Defer one tick so the
      # layout pass has run before we ask for the viewport width.
      QtCore.QTimer.singleShot(0, self._refresh_inner_heights)

  def _refresh_inner_heights(self):
    for view in (self.args_view, self.result_view):
      refresh = getattr(view, "_auto_height_refresh", None)
      if refresh is None:
        continue
      # The deferred singleShot may fire after the widget's C++ object
      # has been destroyed (e.g. parent QWidget garbage-collected while
      # the callback was still in Qt's event queue). Guard the view
      # access so a stale callback no-ops instead of raising.
      try:
        if not view.isHidden():
          refresh()
      except RuntimeError:
        # 'Internal C++ object already deleted' -- widget is gone.
        return

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

  def _refresh_header(self):
    arrow = "▾" if self.header_button.isChecked() else "▸"
    status, clamped = self._header_status()
    # A clamped status keeps the whole of itself on the tooltip. Otherwise
    # leave the tooltip to the button, which shows the full header text only
    # when it had to elide it -- so an ordinary 'finished' row at a normal
    # width sprouts no redundant hover label.
    tooltip = None
    if clamped:
      tooltip = "%s (%s)" % (self._name, self._status_text())
    self.header_button.set_full_text(
      "%s %s (%s)" % (arrow, self._name, status), tooltip=tooltip)
    self.header_button.setStyleSheet(_COLORS.get(self._color, ""))
