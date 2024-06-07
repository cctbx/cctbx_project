from dataclasses import replace

from PySide2.QtGui import QCursor
from PySide2.QtWidgets import (
    QApplication,
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout
)

from ..controller import Controller

_active_toasts = []


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

  def exec_(self):
    # Get the cursor's current position
    cursorPos = QCursor.pos()

    # Get the dimensions of the dialog
    dialog_width = self.frameSize().width()
    dialog_height = self.frameSize().height()

    # Calculate the new top-left position to center the dialog on the cursor
    new_x = cursorPos.x() - dialog_width // 2
    new_y = cursorPos.y() - dialog_height // 2

    # Get the screen geometry based on the cursor's position
    screen = QApplication.desktop().screenNumber(cursorPos)
    screen_geom = QApplication.desktop().screenGeometry(screen)

    # Adjust the dialog's position to ensure it stays within the screen boundaries
    # Check right boundary
    if new_x + dialog_width > screen_geom.right():
        new_x = screen_geom.right() - dialog_width

    # Check bottom boundary
    if new_y + dialog_height > screen_geom.bottom():
        new_y = screen_geom.bottom() - dialog_height

    # Check left boundary
    if new_x < screen_geom.left():
        new_x = screen_geom.left()

    # Check top boundary
    if new_y < screen_geom.top():
        new_y = screen_geom.top()

    # Move dialog to the adjusted position
    self.move(new_x, new_y)

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
    self.map_ref.id
    QMessageBox.information(self.view, 'Notice', "Opacity slider not yet connected")
    real_value = value/100
    self.view.label.setText(f"Opacity: {round(real_value,2)}")
