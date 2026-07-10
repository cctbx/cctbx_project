"""Single-line collapsible row used inside MessageBubble for tool calls.

Used for tool calls and tool-approval prompts. Header:
``▸ tool_name (status)`` -- click to expand to ``▾``. Body (hidden by
default): pretty-printed JSON args on top, plain-text result below.
Status text and color are settable so the runner can transition
running → finished/failed/cancelled without restructuring the widget.
"""

import json

from qttbx.qt import QtCore, QtWidgets


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

    self.header_button = QtWidgets.QToolButton(self)
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
        New status text shown in the header parenthetical.
    color : str or None, optional
        One of ``None``, ``'default'``, ``'running'``, ``'error'``, or
        ``'muted'``. ``None`` preserves the current color.
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
        view.
    """
    if not text:
      self.result_view.setPlainText("")
      self.result_view.hide()
      return
    self.result_view.setPlainText(text)
    self.result_view.show()

  def is_running(self):
    """Return True while the call is still in its initial ``running`` state.

    The predicate the turn-cancel sweep uses to find tool cells that never
    reached a terminal state (finished / failed / cancelled) -- their result
    will never arrive, so they would otherwise stay stuck spinning.
    """
    return self._status == "running"

  # ---- internals ----------------------------------------------------------

  def _on_toggled(self):
    self.body.setVisible(self.header_button.isChecked())
    self._refresh_header()
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

  def _refresh_header(self):
    arrow = "▾" if self.header_button.isChecked() else "▸"
    self.header_button.setText("%s %s (%s)" % (arrow, self._name, self._status))
    self.header_button.setStyleSheet(_COLORS.get(self._color, ""))
