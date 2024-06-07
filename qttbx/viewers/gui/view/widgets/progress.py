"""
An attempt to make a circular progress indicator. Does not work yet.
"""
from pathlib import Path

from PySide2.QtCore import QTimer
from PySide2.QtSvg import QGraphicsSvgItem
from PySide2.QtWidgets import (
    QDialog,
    QGraphicsScene,
    QGraphicsView,
    QVBoxLayout
)

class CircularProgressIndicator(QDialog):
  def __init__(self):
    super().__init__()
    self.setModal(True)
    self.setWindowTitle("Please wait...")
    self.setFixedSize(200, 200)

    layout = QVBoxLayout()

    # Create QGraphicsView and QGraphicsScene
    self.view = QGraphicsView(self)
    self.scene = QGraphicsScene(self)
    self.view.setScene(self.scene)

    # Add SVG item
    icon_path = Path(__file__).parent / '../assets/icons/material/progress.svg'
    self.svg_item = QGraphicsSvgItem(str(icon_path))
    self.svg_item.setTransformOriginPoint(self.svg_item.boundingRect().width()/2, self.svg_item.boundingRect().height()/2)
    self.scene.addItem(self.svg_item)

    layout.addWidget(self.view)

    self.setLayout(layout)

    # Timer for rotation
    self.timer = QTimer(self)
    self.angle = 0
    self.timer.timeout.connect(self.rotate_svg)
    self.timer.start(50)

  def rotate_svg(self):
    self.angle += 5
    if self.angle == 360:
      self.angle = 0
    self.svg_item.setRotation(self.angle)
