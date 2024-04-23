
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2 import QtCore
from PySide2.QtCore import QObject, QAbstractTableModel,  Qt, QTimer, QPoint, Signal

from ..controller import Controller
from ...state.geometry import Geometry
from ...state.ref import GeometryRef
from ..table import TableController
from ...state.table import PandasTableModel
from ...state.ref import SelectionRef
from ...state.base import ObjectFrame
from ...view.widgets import (
  BondEditDialog,
  AngleEditDialog,
  DihedralEditDialog,
  ChiralEditDialog,
  PlaneEditDialog
)
from ...state.edits import BondEdit
from ...state.ref import EditsRef

from pathlib import Path
import pandas as pd

"""
Controller for a Pandas Table in the restraints tab
"""

class RestraintsTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  restraint_name = None # Subclass and fill this in (bond,angle,etc)
  columns_to_include= ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  column_prefixes_to_include = ["atom_id"]
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.restraints_ref = None


    # Signals
    self.view.table_view.addEdit.connect(self.addEdit)
    self.state.signals.restraints_change.connect(self.update)

  def addEdit(self,row_dict):
    raise NotImplementedError


  def update(self,restraint_ref):
    self.restraints_ref = restraint_ref
    print("updating restraints table with ref:",restraint_ref.id)
    if hasattr(restraint_ref.data,self.restraint_name):
      df = getattr(restraint_ref.data,self.restraint_name)
      cols = ["Chain","Res","Seq"]
      # add prefix columns
      for col in df.columns:
        for prefix in self.column_prefixes_to_include:
          if col.startswith(prefix):
            cols.append(col)
      # add exact match columns
      cols+= [col for col in self.columns_to_include if col in df.columns]
      # add i_seqs (list of all i_seqs)
      cols+= ["i_seqs"]
      # TODO: reorder cols without copy? How.
      self.table_model = PandasTableModel(df[cols].copy(),suppress_columns=["i_seqs"])
      # self.table = table
      # self.table._parent_explicit = self
    print(self.parent.view.hiddenTabs)
    if restraint_ref is not None:
      self.parent.view.toggle_tab_visible(self.title,show=True)


class BondTableController(RestraintsTableTabController):
  title="Bonds"
  restraint_name = "bond"

  def addEdit(self,row_dict):
    dialog = BondEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      print("Making edit")

      row = BondEdit(
        i_seqs = row_dict["i_seqs"],
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,EditsRef):
          if ref.restraints_ref == self.restraints_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = EditsRef(data=objframe,restraints_ref=self.restraints_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class AngleTableController(RestraintsTableTabController):
  title = "Angles"
  restraint_name = "angle"

class DihedralTableController(RestraintsTableTabController):
  title = "Dihedrals"
  restraint_name = "dihedral"

class ChiralTableController(RestraintsTableTabController):
  title = "Chirals"
  restraint_name = "chirality"

class PlanarityTableController(RestraintsTableTabController):
  title = "Planes"
  restraint_name = "plane"

class NonbondedTableController(RestraintsTableTabController):
  title = "Non-bonded"
  restraint_name = "nonbonded"



