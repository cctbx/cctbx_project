from PySide2.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget
)


class EditsEditDialog(QDialog):
  input_names = ["ideal", "sigma"]

  def __init__(self, parent=None, input_names=None, defaults_dict=None):
    super().__init__(parent)
   
    self.setWindowTitle(f'Selection Dialog')
    mainLayout = QVBoxLayout()
    
    
    self.setLayout(mainLayout)

    self.buttonsLayout = QHBoxLayout()
    self.cancelButton = QPushButton('Cancel', self)
    self.cancelButton.clicked.connect(self.reject)
    self.buttonsLayout.addWidget(self.cancelButton)

    self.acceptButton = QPushButton('Accept', self)
    self.acceptButton.clicked.connect(self.accept)
    self.buttonsLayout.addWidget(self.acceptButton)

    mainLayout.addLayout(self.buttonsLayout)

  def collectInputValues(self):
    # Collect values from QLineEdit widgets
    values = {name: field.text() for name, field in self.inputFields.items()}
    values["action"] = self.action
    return values

  def exec_(self):
    # Get the cursor's current position
    cursorPos = QCursor.pos()

    # Move the dialog to the cursor's position initially
    self.move(cursorPos.x(), cursorPos.y())

    # Adjust position to ensure the dialog is fully visible on the screen
    screen = QApplication.desktop().screenNumber(cursorPos)
    screen_geom = QApplication.desktop().screenGeometry(screen)

    # Calculate the dialog's geometry after the initial move
    dialog_geom = self.geometry()
    dialog_x = cursorPos.x()
    dialog_y = cursorPos.y()
    dialog_width = dialog_geom.width()
    dialog_height = dialog_geom.height()

    # Check right boundary
    if dialog_x + dialog_width > screen_geom.right():
        dialog_x = screen_geom.right() - dialog_width

    # Check bottom boundary
    if dialog_y + dialog_height > screen_geom.bottom():
        dialog_y = screen_geom.bottom() - dialog_height

    # Check left boundary
    if dialog_x < screen_geom.left():
        dialog_x = screen_geom.left()

    # Check top boundary
    if dialog_y < screen_geom.top():
        dialog_y = screen_geom.top()

    # Move dialog to the adjusted position
    self.move(dialog_x, dialog_y)
    
    result = super().exec_()
    # Collect values if dialog was accepted
    if result == QDialog.Accepted:
      return self.collectInputValues()
    else:
      return None

class TableFilter(QWidget):
  """
  A filter panel to make a table less overwhelming
  """
  def __init__(self, parent=None):
    super().__init__(parent=parent)
    
    # Main vertical layout
    self.layout = QHBoxLayout()

    spacer = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
    self.layout.addSpacerItem(spacer)

    # Dropdown for filter selections
    generic_filter_layout = QVBoxLayout()
    self.combobox_filter = QComboBox()
    self.combobox_filter.setMaximumWidth(200)
    generic_filter_layout.addWidget(QLabel("Filters"))
    generic_filter_layout.addWidget(self.combobox_filter)
    self.layout.addLayout(generic_filter_layout)


    # Dropdown for component selections
    component_filter_layout = QVBoxLayout()
    self.combobox_comp = QComboBox()
    self.combobox_comp.setMaximumWidth(80)
    component_filter_layout.addWidget(QLabel("Components"))
    component_filter_layout.addWidget(self.combobox_comp)
    self.layout.addLayout(component_filter_layout)

    # Reset
    reset_filter_layout = QVBoxLayout()
    self.reset_button = QPushButton("Reset")
    self.reset_button.setMaximumWidth(80)
    reset_filter_layout.addWidget(QLabel(""))
    reset_filter_layout.addWidget(self.reset_button)
    self.layout.addLayout(reset_filter_layout)

    # Set the layout to the QWidget
    self.setLayout(self.layout)