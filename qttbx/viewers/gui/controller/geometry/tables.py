
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2 import QtCore
from PySide2.QtCore import QObject, QAbstractTableModel,  Qt, QTimer, QPoint, Signal

from ..controller import Controller
from ..filter import TableFilterController
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
from ...state.edits import (
  BondEdit,
  AngleEdit,
  DihedralEdit
)
from ...state.ref import (
  BondEditsRef,
  AngleEditsRef,
  DihedralEditsRef

)

from pathlib import Path
import pandas as pd

"""
Controller for a Pandas Table in geometry
"""


class GeometryTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  restraint_name = None # Subclass and fill this in (bond,angle,etc)
  column_order =  ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  suppress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old"]
  rename_columns = {"ideal_new":"Ideal","sigma_new":"Sigma","action":"Action"}
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.restraints_ref = None
    self.col_names = None

    # Filter panel
    self.filter = TableFilterController(parent=self,view=self.view.filter)

    # Signals
    self.view.table_view.addEdit.connect(self.addEdit)
    self.state.signals.restraints_change.connect(self.update)
    self.state.signals.filter_update.connect(self.update_from_filter)

  def addEdit(self,row_dict):
    raise NotImplementedError

  def update_from_filter(self,filter):
    
    df = self.dataframe
    comp = filter.current_comp
    direction = filter.current_percentile_direction
    metric = "residual"
    
    # Determine the quantile based on the direction
    if direction == "Ascending":
        quantile = (float(filter.current_percentile) / 100)
    else:  # Descending
        quantile = 1 - (float(filter.current_percentile) / 100)

    # Filter by component
    if comp == "All":
        comp_test = pd.Series([True] * len(df))
    else:
        comp_test = df["Res"] == comp

    # Apply the quantile filtering
    if metric in df.columns:
        percentile_cutoff = df[metric].quantile(quantile)
        if direction == "Ascending":
            percentile_test = df[metric] <= percentile_cutoff
        else:  # Descending
            percentile_test = df[metric] >= percentile_cutoff
    else:
        percentile_test = pd.Series([True] * len(df))

    # Apply both conditions to filter the DataFrame
    df_filtered = df[comp_test & percentile_test]
    df_filtered = df_filtered.rename(columns=self.rename_columns)
    suppress_cols = [c.lower() for c in self.suppress_columns]
    suppress_cols+= [c.capitalize() for c in self.suppress_columns]
    self.table_model = PandasTableModel(df_filtered,suppress_columns=suppress_cols)

  def reform_columns(self,df):
    cols = ["Chain","Res","Seq"]
    # add atom_id prefix columns
    for col in df.columns:
      if col.startswith("atom_id"):
        cols.append(col)

    cols+= [col for col in self.column_order if col in df.columns]
    df = df[cols].copy()
    # add atom_id prefix columns
    for col in cols:
      if col.startswith("atom_id"):
        number = col.split("atom_id")[-1]
        number = number.strip("_")
        col_name = f"Atom {number}"
        self.rename_columns[col] = col_name

    df.rename(columns=self.rename_columns,inplace=True)
    

    # capitalize
    capitalize = {c:c.capitalize() for c in df.columns}
    df.rename(columns=capitalize,inplace=True)

    self.col_names = list(df.columns)
    return df

  def update(self,restraint_ref):
    self.restraints_ref = restraint_ref
    #print("updating restraints table with ref:",restraint_ref.id)
    if hasattr(restraint_ref.data,self.restraint_name):
      df = getattr(restraint_ref.data,self.restraint_name)
      self.dataframe = self.reform_columns(df)
      self.table_model = PandasTableModel(self.dataframe,suppress_columns=[])
      # self.table = table
      # self.table._parent_explicit = self
    if restraint_ref is not None:
      self.parent.view.toggle_tab_visible(self.title,show=True)
    
    # update filter
    self.filter.update_quiet()



  # A generic pandas helper function to match a row_dict with an actual row in a dataframe
  @staticmethod 
  def find_matching_row(df, dict_subset):
    # Convert the dictionary to a pandas Series for comparison
    series = pd.Series(dict_subset)
    
    # Check if each column in the DataFrame contains the corresponding series value
    mask = df[list(dict_subset.keys())].eq(series)
    
    # Find rows where all specified dictionary key-value pairs match
    matching_rows = df[mask.all(axis=1)]
    
    # Return the first matching row or None if no match is found
    return matching_rows.iloc[0] if not matching_rows.empty else None

  @staticmethod
  def get_sel_strings_from_iseqs(mol,i_seqs):
    # return a list of string selections for an iterable of i_seqs
    sel_strings = []
    for i in i_seqs:
      s = mol.sites.select_query_from_i_seqs([i]).phenix_string
      sel_strings.append(s)
    return sel_strings

class BondTableController(GeometryTableTabController):
  title="Bonds"
  restraint_name = "bond"

  def addEdit(self,row_dict):
    df = getattr(self.restraints_ref.data,self.restraint_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = BondEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      print("Making edit")

      row = BondEdit(
        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        action="mod"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,BondEditsRef):
          if ref.restraints_ref == self.restraints_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = BondEditsRef(data=objframe,restraints_ref=self.restraints_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class AngleTableController(GeometryTableTabController):
  title = "Angles"
  restraint_name = "angle"

  def addEdit(self,row_dict):
    df = getattr(self.restraints_ref.data,self.restraint_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = AngleEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      row = AngleEdit(
        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        action="mod"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,AngleEditsRef):
          if ref.restraints_ref == self.restraints_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = AngleEditsRef(data=objframe,restraints_ref=self.restraints_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class DihedralTableController(GeometryTableTabController):
  title = "Dihedrals"
  restraint_name = "dihedral"

  def addEdit(self,row_dict):
    df = getattr(self.restraints_ref.data,self.restraint_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = DihedralEditDialog(defaults_dict=row_dict)
    
    if dialog.exec_():
      row = DihedralEdit(

        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        weight_old = dialog.collectInputValues()["weight"],
        weight_new=row_dict["weight"],
        harmonic_old = dialog.collectInputValues()["harmonic"],
        harmonic_new=row_dict["harmonic"],
        action="mod"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,DihedralEditsRef):
          if ref.restraints_ref == self.restraints_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = DihedralEditsRef(data=objframe,restraints_ref=self.restraints_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class ChiralTableController(GeometryTableTabController):
  title = "Chirals"
  restraint_name = "chirality"

class PlanarityTableController(GeometryTableTabController):
  title = "Planes"
  restraint_name = "plane"

class NonbondedTableController(GeometryTableTabController):
  title = "Non-bonded"
  restraint_name = "nonbonded"



