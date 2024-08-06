from pathlib import Path

import qtawesome as qta
from PySide2.QtCore import Qt
from PySide2.QtGui import QIcon
from PySide2.QtCore import QObject, Signal
from PySide2.QtWidgets import (
    QComboBox,
    QFrame,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
    QLabel
)

from ..table import PandasTableView
from ..widgets.tab import GUITab



class CifTableView(PandasTableView):

  def __init__(self, parent=None):
    super().__init__(parent=parent)


class CifBrowserTabView(GUITab):
  """
  View cif structure
  """

  #new_block = Signal(str)
  #was_modified = Signal(bool)

  def __init__(self,parent=None,load_button=True,next_button=False):
    super().__init__(parent=parent)
    self.layout = QVBoxLayout()
    

    # combo boxes
    self.combo_button_layout = self.build_comboboxes(load_button=load_button,next_button=next_button)
    self.layout.addItem(self.combo_button_layout)


    # Cif table
    self.table_view = CifTableView()
    self.layout.addWidget(self.table_view)
    self.setLayout(self.layout)



  def build_comboboxes(self,load_button=True,next_button=False):
    # Horizontal layout for comboboxes and save button
    combo_button_layout = QHBoxLayout()

    # add load button if requested
    if load_button:
      load_button_layout = QVBoxLayout()
      self.load_button = QPushButton()
      load_icon = qta.icon("mdi.plus")
      self.load_button.setIcon(load_icon)
      # self.load_button.setMaximumSize(50, 50)
      # self.load_button.setContentsMargins(10, 10, 0, 0)
      load_button_layout.addWidget(self.load_button)
      #combo_button_layout.addWidget(self.load_button,1)
      combo_button_layout.addItem(load_button_layout)



    # combobox for files
    self.combobox_files = QComboBox()
    #combo_button_layout.addWidget(self.combobox_files,4)

    # add label
    files_layout = QVBoxLayout()
    #files_layout.addWidget(QLabel("  File"))
    files_layout.addWidget(self.combobox_files)
    combo_button_layout.addLayout(files_layout,4)



    # Create a combobox for data level keys
    self.combobox_data_block = QComboBox()
    #combo_button_layout.addWidget(self.combobox_data,4)

    # add label
    self.data_block_layout = QVBoxLayout()
    #data_block_layout.addWidget(QLabel("  Data"))
    self.data_block_layout.addWidget(self.combobox_data_block)
    combo_button_layout.addLayout(self.data_block_layout,4)

    # Create a combobox for block level keys
    self.combobox_data_item = QComboBox()
    #combo_button_layout.addWidget(self.combobox_block,4)

    # add label
    data_item_layout = QVBoxLayout()
    #data_item_layout.addWidget(QLabel("  Item"))
    data_item_layout.addWidget(self.combobox_data_item)
    combo_button_layout.addLayout(data_item_layout,4)

    # #Create the vertical separator
    # separator = QFrame()
    # separator.setFrameShape(QFrame.VLine)
    # separator.setFrameShadow(QFrame.Sunken)
    # combo_button_layout.addWidget(separator,1)

    if next_button:
      # Next button
      next_button_layout = QVBoxLayout()
      self.next_button = QPushButton()
      self.next_button.setToolTip("Move to 'next' instance")
      next_icon = qta.icon("mdi.arrow-right")
      self.next_button.setIcon(next_icon)
      next_button_layout.addWidget(self.next_button)
      combo_button_layout.addLayout(next_button_layout)

    # Edit button
    edit_button_layout = QVBoxLayout()
    self.edit_button = QPushButton()
    edit_icon = qta.icon("mdi.pencil")
    self.edit_button.setIcon(edit_icon)
    edit_button_layout.addWidget(self.edit_button)
    combo_button_layout.addLayout(edit_button_layout)

    # Modified notification label
    modified_label_layout = QVBoxLayout()
    self.modified_label = QLabel("")
    modified_label_layout.addWidget(self.modified_label)
    #combo_button_layout.addWidget(self.modified_label,1)
    combo_button_layout.addLayout(modified_label_layout)
    self.combo_button_layout = combo_button_layout


    # Save button
    save_button_layout = QVBoxLayout()
    self.save_button = QPushButton()
    save_icon = qta.icon("mdi.floppy")
    self.save_button.setIcon(save_icon)
    # self.save_button.setMaximumSize(32, 48)
    # self.save_button.setMinimumSize(32, 32)
    #combo_button_layout.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
    #self.layout.setAlignment(Qt.AlignVCenter)
    #combo_button_layout.addWidget(self.save_button,1)
    save_button_layout.addWidget(self.save_button)
    combo_button_layout.addLayout(save_button_layout)

    return combo_button_layout

  def set_was_modified(self,was_modified):
    if was_modified:
      self.modified_label.setText(" ** modified ** ")
    else:
      self.modified_label.setText("")

  def setComboBoxToValue(self,comboBox, value):
    index = comboBox.findText(value)
    if index >= 0:  # The value was found
      comboBox.setCurrentIndex(index)
    else:
      self.log(f"Value '{value}' not found in ComboBox.")

