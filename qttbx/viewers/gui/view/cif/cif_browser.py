from pathlib import Path

from PySide2.QtWidgets import QComboBox, QFrame, QHBoxLayout, QPushButton, QVBoxLayout
from PySide2.QtGui import QIcon

from ..table import PandasTableView
from ..widgets.tab import GUITab


class CifTableView(PandasTableView):

  def __init__(self, parent=None):
    super().__init__(parent=parent)


class CifBrowserTabView(GUITab):
  """
  View cif structure
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.layout = QVBoxLayout()
    self.combo_button_layout = self.build_comboboxes()
    self.layout.addItem(self.combo_button_layout)
    self.table_view = CifTableView()
    self.layout.addWidget(self.table_view)
    self.setLayout(self.layout)


  def build_comboboxes(self):
    # Horizontal layout for comboboxes and save button
    combo_button_layout = QHBoxLayout()

    # combobox for files
    self.combobox_files = QComboBox()
    combo_button_layout.addWidget(self.combobox_files,4)

    # Create a combobox for data level keys
    self.combobox_data = QComboBox()
    combo_button_layout.addWidget(self.combobox_data,4)

    # Create a combobox for block level keys
    self.combobox_block = QComboBox()
    combo_button_layout.addWidget(self.combobox_block,4)

    #Create the vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)
    combo_button_layout.addWidget(separator,1)

    # Save button
    self.save_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/save.svg'
    save_icon = QIcon(str(icon_path))
    self.save_button.setIcon(save_icon)
    # self.save_button.setMaximumSize(32, 48)
    # combo_button_layout.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
    #self.layout.setAlignment(Qt.AlignVCenter)
    combo_button_layout.addWidget(self.save_button,1)

    return combo_button_layout

  def setComboBoxToValue(self,comboBox, value):
    index = comboBox.findText(value)
    if index >= 0:  # The value was found
      comboBox.setCurrentIndex(index)
    else:
      print(f"Value '{value}' not found in ComboBox.")


# class CifBrowserTabView(GUITab):
#   """
#   View cif structure
#   """
#   def __init__(self,parent=None):
#     super().__init__(parent=parent)
#     layout = QVBoxLayout()
#     self.layout = layout





#     # add empty dataframe
#     self.setLayout(layout)
#     self.table_view = CifTableView()
#     layout.addWidget(self.table_view)



