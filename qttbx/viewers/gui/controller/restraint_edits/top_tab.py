from pathlib import Path
import webbrowser
import os

from PySide2.QtWidgets import QFileDialog
import iotbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str

from ...state.base import ObjectFrame
from ..controller import Controller
from .tables import (
  BondTableController,
  AngleTableController,
  DihedralTableController,
  ChiralTableController,
  PlanarityTableController
)
from ...state.edits import *


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
    
    gm_phil = iotbx.phil.parse(
      input_string     = grand_master_phil_str,
      process_includes = True)
    edit_phil = gm_phil.get("geometry_restraints")
    output_edits = edit_phil.extract()
    #edits = self.params.geometry_restraints.edits
    name_classes = {
      "bond":BondEdit,
      "angle":AngleEdit,
      "dihedral":DihedralEdit,
    }
    for name,class_name in name_classes.items():
      table = getattr(self,name+"s")
      phil_objs = []
      if table.df is not None:
        objs = ObjectFrame.from_df(table.df,row_class=class_name)
        for obj in objs:
          phil_obj = obj.as_phil_obj()
          phil_objs.append(phil_obj)
        setattr(output_edits.geometry_restraints.edits,name,phil_objs)
    
    output_phil = edit_phil.format(python_object=output_edits)
    output_phil = edit_phil.fetch_diff(output_phil)
    edit_string = output_phil.as_str()
    print(edit_string)

    # Suggest a default filename and directory
    home_dir = Path.home()
    defaultFileName = str(Path(home_dir,"restraint_edits.params"))  # Change to a relevant path
    # Open the save file dialog
    fileName, _ = QFileDialog.getSaveFileName(self.view, 'Save File', defaultFileName,
                                              'All Files (*);;Text Files (*.txt)')
    if fileName:
      self.log(f'File selected for saving: {fileName}')
      with open(fileName,"w") as fh:
        fh.write(edit_string)

      # # Open the file in the default application
      # file_path = fileName
      # file_path = os.path.realpath(file_path)
      # url = 'file://' + file_path
      # webbrowser.open(url)



