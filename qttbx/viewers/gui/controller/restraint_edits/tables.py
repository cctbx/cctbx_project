
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2 import QtCore
from ..controller import Controller
from ..table import TableController

from ...state.geometry import Geometry
from ...state.ref import GeometryRef
from ...state.table import  PandasTableModel
from ...state.ref import SelectionRef, EditsRef
from ...state.edits import (
  EditData,
  BondEdit,
  AngleEdit,
  DihedralEdit,
)
from pathlib import Path
import pandas as pd



class EditsTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  row_class = EditData # Generic, use subclasses
  restraint_name =  None # Subclass and fill this in (bond,angle,etc)
  #columns_to_include= ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  supress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old"]
  rename_columns = {"ideal_new":"Ideal","sigma_new":"Sigma","action":"Action"}
  column_prefixes_to_include = ["atom_id"]
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.edits_ref = None

    # Signals
    self.state.signals.edits_change.connect(self.update)
    self.view.table_view.removeEdit.connect(self.remove_edit)
    self.view.write_button.clicked.connect(self.write_edits)

  def update(self,edits_ref):
    if edits_ref.data.row_class==self.row_class:
      self.edits_ref = edits_ref
      self.dataframe = edits_ref.data.df
      if len(self.dataframe)>0:
        self.dataframe = self.dataframe.rename(columns=self.rename_columns)
        self.table_model = PandasTableModel(self.dataframe,suppress_columns=self.supress_columns)
        self.parent.view.toggle_tab_visible(self.title,show=True)
      else:
        self.state.signals.remove_ref.emit(edits_ref)
        self.parent.view.toggle_tab_visible(self.title,show=False)

  def remove_edit(self,row_dict):
    if self.edits_ref is not None:
      edit = self.row_class.from_dict(row_dict)
      self.edits_ref.data.remove(edit)
      self.update(self.edits_ref)

  def write_edits(self):
    self.parent.write_all_edits()



class BondTableController(EditsTableTabController):
  title = "Bonds"
  restraint_name = "bond"
  row_class = BondEdit
  supress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old"]
  rename_columns = {"ideal_new":"Ideal","sigma_new":"Sigma","action":"Action"}

class AngleTableController(EditsTableTabController):
  title = "Angles"
  restraint_name = "angle"
  row_class = AngleEdit
  supress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old"]
  rename_columns = {"ideal_new":"Ideal","sigma_new":"Sigma","action":"Action"}

class DihedralTableController(EditsTableTabController):
  title = "Dihedrals"
  restraint_name = "dihedral"
  row_class = DihedralEdit

class ChiralTableController(EditsTableTabController):
  title = "Chirals"
  restraint_name = "chirality"

class PlanarityTableController(EditsTableTabController):
  title = "Planes"
  restraint_name = "plane"

class NonbondedTableController(EditsTableTabController):
  title = "Nonbonded"
  restraint_name = "nonbonded"



