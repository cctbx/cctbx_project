"""Thin horizontal strip above the conversation.

Shows the title / model / debug log path. No interactive elements --
ChatWindow updates the labels as the conversation, profile, or
debug-log state changes.
"""

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.eliding import ElidingLabel

# Widest the model / debug slots may get. They are secondary chrome, and their
# natural widths are large (a proxy model id runs ~410 px, a debug path
# ~500 px) -- uncapped they would claim most of a wide bar and leave the title
# elided next to them.
_SIDE_SLOT_MAX_WIDTH = 200

# Share of the bar each slot gets when there is not enough for all three.
# Nothing here floors any more, so QHBoxLayout hands the whole shortfall to
# the stretch item: with the title alone stretching, it collapsed to zero
# width at any bar under ~900 px while the debug path kept 500 px of it.
# Giving every slot a stretch spreads the shortfall in these proportions, so
# the title -- the bar's actual content -- keeps roughly two thirds and the
# side slots elide first.
_TITLE_STRETCH = 4
_SIDE_STRETCH = 1


class ChatTopBar(QtWidgets.QFrame):
  """Non-interactive header strip showing title, model, and debug log."""

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self.setFixedHeight(28)
    layout = QtWidgets.QHBoxLayout(self)
    layout.setContentsMargins(8, 2, 8, 2)
    layout.setSpacing(12)

    # Every slot here carries unbounded text: the title is
    # '<profile> / <first 60 chars of first user text>', the model id comes
    # straight through from --model, and the debug slot holds a filesystem
    # path. A plain QLabel reports its full text width as its
    # minimumSizeHint, which floors the centre column -- and with it the
    # whole window -- at a width the user then cannot resize below.
    # ElidingLabel drops that floor to zero and keeps each slot readable by
    # eliding it to whatever width the bar actually gives it.
    self.title_label = ElidingLabel(self)
    self.title_label.setStyleSheet("font-weight: bold;")
    self.model_label = ElidingLabel(self)
    self.model_label.setStyleSheet("color: palette(mid);")
    # Elide from the left so the basename of the debug log stays visible.
    self.debug_label = ElidingLabel(self, mode=QtCore.Qt.ElideLeft)
    self.debug_label.setStyleSheet(
      "color: palette(mid); font-style: italic;")
    self.debug_label.setTextInteractionFlags(
      QtCore.Qt.TextSelectableByMouse)
    self.model_label.setMaximumWidth(_SIDE_SLOT_MAX_WIDTH)
    self.debug_label.setMaximumWidth(_SIDE_SLOT_MAX_WIDTH)

    layout.addWidget(self.title_label, stretch=_TITLE_STRETCH)
    layout.addWidget(self.model_label, stretch=_SIDE_STRETCH)
    layout.addWidget(self.debug_label, stretch=_SIDE_STRETCH)
    self.debug_label.hide()

  def set_title(self, text):
    """Set the title label text, elided to the width the bar gives it."""
    self.title_label.set_full_text(text)

  def set_model(self, text):
    """Set the model label text, elided to the width the bar gives it."""
    self.model_label.set_full_text(text)

  def set_debug_log_path(self, path):
    """Show the debug log path, or hide the slot when there is none.

    Parameters
    ----------
    path : str or pathlib.Path or None
        The debug log path to display. ``None`` hides the slot.
    """
    if not path:
      self.debug_label.set_full_text("")
      self.debug_label.hide()
      return
    # The whole path goes in: _SIDE_SLOT_MAX_WIDTH bounds the slot, the left
    # elide keeps the informative tail, and the label hangs the untruncated
    # path on its own tooltip.
    self.debug_label.set_full_text("debug: %s" % path)
    self.debug_label.show()
