"""Thin horizontal strip above the conversation: title / model / debug
log path. No interactive elements -- ChatWindow updates the labels as
the conversation, profile, or debug-log state changes."""

from qttbx.qt import QtCore, QtWidgets


class ChatTopBar(QtWidgets.QFrame):

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self.setFixedHeight(28)
    layout = QtWidgets.QHBoxLayout(self)
    layout.setContentsMargins(8, 2, 8, 2)
    layout.setSpacing(12)

    self.title_label = QtWidgets.QLabel(self)
    self.title_label.setStyleSheet("font-weight: bold;")
    self.model_label = QtWidgets.QLabel(self)
    self.model_label.setStyleSheet("color: palette(mid);")
    self.debug_label = QtWidgets.QLabel(self)
    self.debug_label.setStyleSheet(
      "color: palette(mid); font-style: italic;")
    # Elide from the left so the basename of the debug log stays visible.
    self.debug_label.setTextInteractionFlags(
      QtCore.Qt.TextSelectableByMouse)

    layout.addWidget(self.title_label, stretch=1)
    layout.addWidget(self.model_label, stretch=0)
    layout.addWidget(self.debug_label, stretch=0)
    self.debug_label.hide()

  def set_title(self, text):
    self.title_label.setText(text or "")

  def set_model(self, text):
    self.model_label.setText(text or "")

  def set_debug_log_path(self, path):
    """path may be a str, pathlib.Path, or None. None hides the slot."""
    if not path:
      self.debug_label.setText("")
      self.debug_label.hide()
      return
    self.debug_label.setText("debug: …%s" % str(path)[-50:])
    self.debug_label.setToolTip(str(path))
    self.debug_label.show()
