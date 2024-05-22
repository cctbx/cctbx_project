
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog

from PySide2 import QtCore
from ..controller import Controller
from .tables import (
  BondTableController,
  AngleTableController,
  DihedralTableController,
  ChiralTableController,
  PlanarityTableController,
  NonbondedTableController
)
from pathlib import Path
import webbrowser
import os
import pandas as pd



class EditsTableTopTabController(Controller):
  """
  Simple controller to manage sub-tabs
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.bonds = BondTableController(parent=self,view=view.bonds)
    self.angles = AngleTableController(parent=self,view=view.angles)
    self.dihedrals = DihedralTableController(parent=self,view=view.dihedrals)
    self.chirals = ChiralTableController(parent=self,view=view.chirals)
    self.planes = PlanarityTableController(parent=self,view=view.planes)
    #self.nonbonded = NonbondedTableController(parent=self,view=view.nonbonded)

    # Start hidden
    self.view.toggle_tab_visible("Bonds",show=False)
    self.view.toggle_tab_visible("Angles",show=False)
    self.view.toggle_tab_visible("Dihedrals",show=False)
    self.view.toggle_tab_visible("Chirals",show=False)
    self.view.toggle_tab_visible("Planes",show=False)


  def write_all_edits(self):
    edits = []

    tables = [self.bonds,self.angles,self.dihedrals]
    for table in tables:
      if table.edits_ref:
        for obj in table.edits_ref.data:
          s = obj.to_edits_string()
          edits.append(s)
    edit_string = "\n".join(edits)
    

    # Suggest a default filename and directory
    home_dir = Path.home()
    defaultFileName = str(Path(home_dir,"restraint_edits.params"))  # Change to a relevant path
    # Open the save file dialog
    fileName, _ = QFileDialog.getSaveFileName(self.view, 'Save File', defaultFileName,
                                              'All Files (*);;Text Files (*.txt)')
    if fileName:
      print(f'File selected for saving: {fileName}')
      with open(fileName,"w") as fh:
        fh.write(edit_string)

      # Open the file in the default application
      file_path = fileName
      file_path = os.path.realpath(file_path)
      url = 'file://' + file_path
      webbrowser.open(url)



