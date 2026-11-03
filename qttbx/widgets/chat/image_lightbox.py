"""Modal full-size image viewer. Borderless, click-anywhere or Esc to
close. Opens from MessageBubble image cells."""

from qttbx.qt import QtCore, QtWidgets


class ImageLightbox(QtWidgets.QDialog):
  """Modal full-size image viewer; click anywhere or press Esc to close.

  Parameters
  ----------
  pixmap : QtGui.QPixmap
      The image to display, scaled to fit the screen.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  def __init__(self, pixmap, parent=None):
    super().__init__(parent)
    self.pixmap = pixmap
    self.setModal(True)
    self.setWindowFlags(self.windowFlags() | QtCore.Qt.FramelessWindowHint)
    self.setAttribute(QtCore.Qt.WA_TranslucentBackground, False)
    # Size: 90% of the available screen, preserving aspect.
    screen = QtWidgets.QApplication.primaryScreen()
    target = screen.availableSize() if screen else QtCore.QSize(800, 600)
    target = QtCore.QSize(int(target.width() * 0.9),
                          int(target.height() * 0.9))
    self.resize(target)

    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    self._label = QtWidgets.QLabel(self)
    self._label.setAlignment(QtCore.Qt.AlignCenter)
    self._label.setPixmap(self._scaled_pixmap(pixmap, target))
    self._label.setStyleSheet("background-color: black;")
    layout.addWidget(self._label)

    # Esc closes for free (default QDialog behavior).

  def mousePressEvent(self, event):
    """Close the dialog on any mouse press."""
    # Click anywhere closes. Don't pre-filter by button -- any press
    # ends the modal.
    self._on_clicked()

  def _on_clicked(self):
    self.accept()

  @staticmethod
  def _scaled_pixmap(pm, target_size):
    return pm.scaled(
      target_size,
      QtCore.Qt.KeepAspectRatio,
      QtCore.Qt.SmoothTransformation)
