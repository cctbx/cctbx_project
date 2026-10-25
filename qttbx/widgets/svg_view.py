"""SVG viewer: one SVG file shown whole, read by zooming and scrolling.

The drawing is Qt's own (QtSvg), so there is no SVG code here. The picture
opens fitted to the window and stays fitted as the window is resized until
the user zooms; it is never split into pieces.

    from qttbx.widgets.svg_view import SvgViewer
    viewer = SvgViewer()
    viewer.load("graph_report.svg")
    viewer.show()

or as a window of its own:

    libtbx.python -m qttbx.widgets.svg_view graph_report.svg

QtSvg implements SVG Tiny 1.2. Shapes and text are drawn; markers (the usual
way to put a head on an arrow) are not, so a file meant for this viewer
should draw its arrow heads as shapes.
"""

import os
import sys

from qttbx.qt import QtCore, QtGui, QtWidgets, QT_API

# Qt6 moved the SVG widgets out of QtSvg. Branch on QT_API as qttbx.qt does.
if QT_API == "PySide6":
  from PySide6.QtSvgWidgets import QGraphicsSvgItem
else:
  from PySide2.QtSvg import QGraphicsSvgItem


class SvgView(QtWidgets.QGraphicsView):
  """The picture itself: fitted until zoomed, then scrolled or dragged."""

  ZOOM_STEP = 1.25
  MINIMUM_SCALE = 0.02
  MAXIMUM_SCALE = 16.0

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setScene(QtWidgets.QGraphicsScene(self))
    self._item = None
    self._fitted = True
    self.setBackgroundBrush(QtGui.QBrush(QtCore.Qt.white))
    self.setRenderHints(QtGui.QPainter.Antialiasing
                        | QtGui.QPainter.TextAntialiasing
                        | QtGui.QPainter.SmoothPixmapTransform)
    self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
    self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)

  def load(self, file_name):
    """Show file_name; False, with the view left empty, if it is not an
    SVG that Qt can read."""
    self.scene().clear()
    self._item = None
    if not os.path.isfile(file_name):
      return False
    item = QGraphicsSvgItem(file_name)
    if not item.renderer().isValid():
      return False
    self._item = item
    self.scene().addItem(item)
    self.scene().setSceneRect(item.boundingRect())
    self.fit()
    return True

  def scale_factor(self):
    return self.transform().m11()

  def set_scale_factor(self, factor):
    factor = max(self.MINIMUM_SCALE, min(self.MAXIMUM_SCALE, factor))
    self.setTransform(QtGui.QTransform.fromScale(factor, factor))

  def fit(self):
    self._fitted = True
    if self._item is not None:
      self.fitInView(self.sceneRect(), QtCore.Qt.KeepAspectRatio)

  def zoom_in(self):
    self._fitted = False
    self.set_scale_factor(self.scale_factor() * self.ZOOM_STEP)

  def zoom_out(self):
    self._fitted = False
    self.set_scale_factor(self.scale_factor() / self.ZOOM_STEP)

  def full_size(self):
    self._fitted = False
    self.set_scale_factor(1.0)

  def resizeEvent(self, event):
    super().resizeEvent(event)
    if self._fitted:
      self.fit()

  def wheelEvent(self, event):
    # Ctrl and the wheel zoom; the wheel alone scrolls, as in a browser.
    if event.modifiers() & QtCore.Qt.ControlModifier:
      if event.angleDelta().y() > 0:
        self.zoom_in()
      else:
        self.zoom_out()
      event.accept()
    else:
      super().wheelEvent(event)


class SvgViewer(QtWidgets.QWidget):
  """SvgView under a row of zoom buttons."""

  def __init__(self, parent=None):
    super().__init__(parent)
    self.view = SvgView(self)
    buttons = QtWidgets.QHBoxLayout()
    for label, slot in (("Fit", self.view.fit),
                        ("Zoom in", self.view.zoom_in),
                        ("Zoom out", self.view.zoom_out),
                        ("Full size", self.view.full_size)):
      button = QtWidgets.QPushButton(label, self)
      button.clicked.connect(slot)
      buttons.addWidget(button)
    buttons.addWidget(QtWidgets.QLabel(
      "Drag or scroll to move; the buttons, or Ctrl and the wheel, zoom.",
      self))
    buttons.addStretch(1)
    layout = QtWidgets.QVBoxLayout(self)
    layout.addLayout(buttons)
    layout.addWidget(self.view, 1)

  def load(self, file_name):
    loaded = self.view.load(file_name)
    if loaded:
      self.setWindowTitle(os.path.basename(file_name))
    return loaded


def run(args):
  from qttbx.widgets.font_init import init_default_app_font
  if len(args) != 1:
    print("usage: libtbx.python -m qttbx.widgets.svg_view <file.svg>")
    return 1
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  viewer = SvgViewer()
  if not viewer.load(args[0]):
    print("Not an SVG file that Qt can read: %s" % args[0])
    return 1
  viewer.resize(700, 800)
  viewer.show()
  if QT_API == "PySide6":
    return app.exec()
  return app.exec_()


if __name__ == "__main__":
  sys.exit(run(sys.argv[1:]))
