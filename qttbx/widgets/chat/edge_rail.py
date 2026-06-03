"""Thin vertical rail with a single toggle button.

ChatWindow puts one EdgeRail on each side. Clicking the button emits
the toggled(bool) signal; ChatWindow connects this to a method that
calls splitter.setSizes(...) to expand or collapse the adjacent panel.
The rail itself owns no panel state -- ChatWindow is the source of
truth and tells the rail what arrow to draw via set_expanded()."""

from qttbx.qt import QtCore, QtWidgets


# Arrows per side and state. 'collapsed' = panel is hidden, button
# points INTO the window (the direction the user would 'open' to).
# 'expanded' = panel is visible, button points OUT (the direction
# the user would 'close' to).
_ARROWS = {
  "left":  {"collapsed": "▶", "expanded": "◀"},
  "right": {"collapsed": "◀", "expanded": "▶"},
}


class EdgeRail(QtWidgets.QFrame):
  """Thin vertical rail with a single panel-toggle button.

  Parameters
  ----------
  side : str
      Which edge the rail sits on; ``"left"`` or ``"right"``.
  tooltip_show : str
      Button tooltip shown while the adjacent panel is collapsed.
  tooltip_hide : str
      Button tooltip shown while the adjacent panel is expanded.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  toggled = QtCore.Signal(bool)  # True = user wants panel expanded.

  def __init__(self, side, tooltip_show, tooltip_hide, parent=None):
    super().__init__(parent)
    assert side in ("left", "right"), side
    self._side = side
    self._tooltip_show = tooltip_show
    self._tooltip_hide = tooltip_hide
    self._expanded = False

    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self.setFixedWidth(18)

    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(0, 6, 0, 0)
    layout.setSpacing(0)

    self.button = QtWidgets.QToolButton(self)
    self.button.setAutoRaise(True)
    self.button.setText(_ARROWS[side]["collapsed"])
    self.button.setToolTip(tooltip_show)
    self.button.setFixedSize(16, 20)
    self.button.clicked.connect(self._on_clicked)

    layout.addWidget(self.button, alignment=QtCore.Qt.AlignHCenter)
    layout.addStretch(1)

  def set_expanded(self, expanded):
    """Update the arrow/tooltip to reflect the panel's expanded state.

    Source-of-truth flip: ChatWindow calls this whenever the adjacent
    panel's actual width crosses the zero/non-zero boundary, keeping
    the arrow consistent with reality regardless of who triggered the
    change (rail / menu / shortcut / splitter drag).

    Does NOT emit ``toggled()`` -- that signal is for user-initiated
    rail clicks only; a feedback loop would otherwise form when
    ChatWindow calls ``set_expanded()`` in response to handling
    ``toggled()``.

    Parameters
    ----------
    expanded : bool
        Whether the adjacent panel is currently expanded.
    """
    self._expanded = bool(expanded)
    state = "expanded" if self._expanded else "collapsed"
    self.button.setText(_ARROWS[self._side][state])
    self.button.setToolTip(
      self._tooltip_hide if self._expanded else self._tooltip_show)

  def _on_clicked(self):
    # Flip the local view; ChatWindow's slot will call set_expanded()
    # back with the authoritative value once it adjusts the splitter.
    new_state = not self._expanded
    self.set_expanded(new_state)
    self.toggled.emit(new_state)
