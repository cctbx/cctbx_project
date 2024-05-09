from pathlib import Path

import pandas as pd
from PySide2.QtWidgets import  QFrame, QDialog, QLabel,QPushButton, QHBoxLayout,QVBoxLayout, QApplication, QMainWindow, QTabWidget, QTableView, QWidget, QComboBox, QSpacerItem, QSizePolicy
from PySide2.QtGui import QStandardItemModel, QStandardItem, QIcon


from ..table import PandasTableView
from ...state.table import PandasTableModel
from ..widgets.tab import GUITab,GUITabWidget



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
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    
    # Main horizontal layout
    self.layout = QHBoxLayout()

    # Dropdown by component
    self.component_layout = QVBoxLayout()
    component_label = QLabel("Component:")
    self.component_layout.addWidget(component_label)
    self.combobox_comp = QComboBox()
    self.combobox_comp.addItem("All")
    self.combobox_comp.setMaximumWidth(80)
    self.component_layout.addWidget(self.combobox_comp)
    self.layout.addLayout(self.component_layout)
    
    # Vertical line as a separator
    line = QFrame()
    line.setFrameShape(QFrame.VLine)
    line.setFrameShadow(QFrame.Sunken)
    self.layout.addWidget(line)

    # Dropdown by percentile
    self.percentile_hlayout = QHBoxLayout()

    # Percentile selection with label on top
    percentile_vlayout = QVBoxLayout()
    percentile_label = QLabel("Percentile:")
    percentile_vlayout.addWidget(percentile_label)
    self.combobox_percentile = QComboBox()
    for percentile in [10, 8, 6, 4, 2]:
        self.combobox_percentile.addItem(str(percentile * 10))
    self.combobox_percentile.setMaximumWidth(80)
    percentile_vlayout.addWidget(self.combobox_percentile)
    self.percentile_hlayout.addLayout(percentile_vlayout)

    # Direction selection with label on top
    direction_vlayout = QVBoxLayout()
    direction_label = QLabel("residual")
    direction_vlayout.addWidget(direction_label)
    self.combobox_percentile_direction = QComboBox()
    self.combobox_percentile_direction.setMaximumWidth(100)
    self.combobox_percentile_direction.addItem("Ascending")
    self.combobox_percentile_direction.addItem("Descending")
    direction_vlayout.addWidget(self.combobox_percentile_direction)
    self.percentile_hlayout.addLayout(direction_vlayout)

    # Add the horizontal layout to the main layout
    self.layout.addLayout(self.percentile_hlayout)

    # Set the layout to the QWidget
    self.setLayout(self.layout)


class HierarchicalFilter(QWidget):
  """
  Filter by hierarchy (chain, resseq, atom)
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.layout = QVBoxLayout()
    self.combo_button_layout = self.build_comboboxes()
    self.layout.addItem(self.combo_button_layout)
    self.setLayout(self.layout)


  def build_comboboxes(self):
    # Horizontal layout for comboboxes and save button
    combo_button_layout = QHBoxLayout()

    # combobox for files
    self.combobox_chain = QComboBox()
    combo_button_layout.addWidget(self.combobox_chain,4)

    # Create a combobox for data level keys
    self.combobox_resseq = QComboBox()
    combo_button_layout.addWidget(self.combobox_resseq,4)

    # Create a combobox for block level keys
    self.combobox_atom = QComboBox()
    combo_button_layout.addWidget(self.combobox_atom,4)

    # #Create the vertical separator
    # separator = QFrame()
    # separator.setFrameShape(QFrame.VLine)
    # separator.setFrameShadow(QFrame.Sunken)
    # combo_button_layout.addWidget(separator,1)

    # # Save button
    # self.save_button = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/save.svg'
    # save_icon = QIcon(str(icon_path))
    # self.save_button.setIcon(save_icon)
    # # self.save_button.setMaximumSize(32, 48)
    # # combo_button_layout.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
    # #self.layout.setAlignment(Qt.AlignVCenter)
    # combo_button_layout.addWidget(self.save_button,1)

    return combo_button_layout

  def setComboBoxToValue(self,comboBox, value):
    index = comboBox.findText(value)
    if index >= 0:  # The value was found
      comboBox.setCurrentIndex(index)
    else:
      print(f"Value '{value}' not found in ComboBox.")
