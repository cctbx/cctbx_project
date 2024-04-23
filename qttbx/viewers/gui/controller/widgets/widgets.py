from PySide2.QtWidgets import QMessageBox
import pandas as pd
import numpy as np
from PySide2.QtGui import QColor

from PySide2 import QtCore
from PySide2.QtCore import QObject, QAbstractTableModel,  Qt, QTimer, QPoint, Signal
from PySide2.QtWidgets import QWidget, QHeaderView, QListView,QTableView, QDialog, QLabel, QVBoxLayout, QHBoxLayout,QWidget, QComboBox, QStyle, QStyleOptionComboBox,  QTextEdit
from PySide2.QtWidgets import  QVBoxLayout, QWidget,  QSlider
from PySide2.QtGui import QMouseEvent, QPainter
from PySide2.QtGui import QFontMetrics
from PySide2.QtCore import QModelIndex
from PySide2.QtGui import QCursor

from PySide2.QtWidgets import QApplication, QTableView, QMenu, QAction
from PySide2.QtWidgets import QTableView, QMenu, QAction, QHeaderView
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure



_active_toasts = []

from PySide2.QtWidgets import QApplication, QDialog, QLineEdit, QPushButton, QVBoxLayout, QLabel

from ..controller import Controller
from dataclasses import replace

class InputDialog(QDialog):
  def __init__(self, parent=None,defaults_dict=None):
    super().__init__(parent)
    self.defaults_dict = defaults_dict
    if self.defaults_dict:
       self.defaults_dict = {key:str(value) for key,value in defaults_dict.items()}
    self.setWindowTitle('Create Geometry Edit')

    # Main layout
    mainLayout = QVBoxLayout()

    # First input field with label
    self.firstInputLayout = QVBoxLayout()
    self.firstInputLabel = QLabel("Ideal")
    self.firstInput = QLineEdit(self)
    default = ""
    if self.defaults_dict:
      if "ideal" in self.defaults_dict:
         default = self.defaults_dict["ideal"]
    self.firstInput.setText(default)  # Set default value
    self.firstInputLayout.addWidget(self.firstInputLabel)
    self.firstInputLayout.addWidget(self.firstInput)

    # Second input field with label
    self.secondInputLayout = QVBoxLayout()
    self.secondInputLabel = QLabel("Sigma")
    self.secondInput = QLineEdit(self)
    default = ""
    if self.defaults_dict:
      if "sigma" in self.defaults_dict:
         default = self.defaults_dict["sigma"]
    self.secondInput.setText(default)  # Set default value
    self.secondInputLayout.addWidget(self.secondInputLabel)
    self.secondInputLayout.addWidget(self.secondInput)

    # Third input field with label
    self.thirdInputLayout = QVBoxLayout()
    self.thirdInputLabel = QLabel("Third Input")
    self.thirdInput = QLineEdit(self)
    self.thirdInput.setText("Default 3")  # Set default value
    self.thirdInputLayout.addWidget(self.thirdInputLabel)
    self.thirdInputLayout.addWidget(self.thirdInput)

    # Combine input fields horizontally
    self.inputsLayout = QHBoxLayout()
    self.inputsLayout.addLayout(self.firstInputLayout)
    self.inputsLayout.addLayout(self.secondInputLayout)
    self.inputsLayout.addLayout(self.thirdInputLayout)

    mainLayout.addLayout(self.inputsLayout)

    # Buttons
    self.buttonsLayout = QHBoxLayout()
    self.cancelButton = QPushButton('Cancel', self)
    self.cancelButton.clicked.connect(self.reject)
    self.buttonsLayout.addWidget(self.cancelButton)

    self.acceptButton = QPushButton('Accept', self)
    self.acceptButton.clicked.connect(self.accept)
    self.buttonsLayout.addWidget(self.acceptButton)

    self.anotherButton = QPushButton('Delete', self)
    # Connect this button to another slot if needed
    self.buttonsLayout.addWidget(self.anotherButton)

    mainLayout.addLayout(self.buttonsLayout)

    self.setLayout(mainLayout)

  def otherAction(self):
    # Placeholder for another action
    print("Another Action triggered")

  def exec_(self):
    # Move the dialog to the cursor's current position
    cursorPos = QCursor.pos()
    self.move(cursorPos.x(), cursorPos.y())  # Move the dialog to cursor position

    # Call the base class exec_ method to show the dialog
    return super().exec_()

class ISOWidgetController(Controller):
  def __init__(self,parent=None,view=None,map_ref=None):
    super().__init__(parent=parent,view=view)
    self.map_ref = map_ref
    self.view.map_ref = map_ref

    # Connect
    self.view.slider.valueChanged.connect(self.iso_update)


  def iso_update(self,value):
    #print("slider update: active_map_ref",self.map_ref.id)
    real_value = self.view.slider.slider_iso_to_real(value,self.view.slider.slider_iso_granularity)
    self.view.parent().iso_label.setText(f"ISO: {round(real_value,2)}")
    style = replace(self.map_ref.style,iso=real_value)
    self.map_ref.style = style




class OpacityWidgetController(Controller):
  def __init__(self,parent=None,view=None,map_ref=None):
    super().__init__(parent=parent,view=view)
    self.map_ref = map_ref
    self.view.slider.valueChanged.connect(self.opacity_update)


  def opacity_update(self,value):
    ref_id = self.map_ref.id
    QMessageBox.information(self.view, 'Notice', "Opacity slider not yet connected")
    real_value = value/100
    self.view.label.setText(f"Opacity: {round(real_value,2)}")
