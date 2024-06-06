

from ..filter import TableFilterController
from ..table import TableController
from ...state.table import PandasTableModel
from ....core.selection import Selection
from ...state.base import ObjectFrame
from ...view.widgets import (
  BondEditDialog,
  AngleEditDialog,
  DihedralEditDialog
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

import pandas as pd

"""
Controller for a Pandas Table in geometry
"""


class GeometryTableTabController(TableController):
  """
  Base class which manipulates the Dataframe
  """
  geometry_name = None # Subclass and fill this in (bond,angle,etc)
  column_order =  ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  suppress_columns = ["i_seqs","sel_strings","ideal_old","sigma_old"]
  rename_columns = {"ideal_new":"Ideal","sigma_new":"Sigma","action":"Action"}
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.geometry_ref = None
    self.col_names = None
    self.debug_flag = False

    # Filter panel
    self.filter = TableFilterController(parent=self,view=self.view.filter)

    # Signals
    self.view.table_view.addEdit.connect(self.addEdit)
    self.state.signals.geometry_change.connect(self.update)
    self.state.signals.filter_update.connect(self.update_from_filter)

  def addEdit(self,row_dict):
    raise NotImplementedError

  def update_from_filter(self,filter_obj,debug_flag):
    df = self.dataframe

    #comp = filter_obj.comp_id
    
    # # Determine the quantile based on the direction
    # if direction == "Ascending":
    #     quantile = (float(filter.current_percentile) / 100)
    # else:  # Descending
    #     quantile = 1 - (float(filter.current_percentile) / 100)

    # Filter by component
    # if comp == "All":
    #     comp_test = pd.Series([True] * len(df))
    # else:
    #     comp_test = df["Res"] == comp

    # # Apply the quantile filtering
    # if metric in df.columns:
    #     percentile_cutoff = df[metric].quantile(quantile)
    #     if direction == "Ascending":
    #         percentile_test = df[metric] <= percentile_cutoff
    #     else:  # Descending
    #         percentile_test = df[metric] >= percentile_cutoff
    # else:
    #     percentile_test = pd.Series([True] * len(df))

    # Apply both conditions to filter the DataFrame
    df_filtered = filter_obj.filter_df(df)
    df_filtered = df_filtered.rename(columns=self.rename_columns)
    suppress_cols = [c.lower() for c in self.suppress_columns]
    suppress_cols+= [c.capitalize() for c in self.suppress_columns]
    table_model = PandasTableModel(df_filtered,suppress_columns=suppress_cols)
    self.table_model = table_model

  def reform_columns(self,df):
    cols = ["i_seqs","Chain","Res","Seq"]
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
        col_name = f"Atom{number}"
        self.rename_columns[col] = col_name

    df.rename(columns=self.rename_columns,inplace=True)
    

    # capitalize
    capitalize = {c:c.capitalize() for c in df.columns if c not in ["i_seqs"]}
    df.rename(columns=capitalize,inplace=True)

    self.col_names = list(df.columns)
    return df

  def update(self,geometry_ref):
    if hasattr(geometry_ref.data,self.geometry_name):
      df = getattr(geometry_ref.data,self.geometry_name)
      self.dataframe = self.reform_columns(df)
      self.table_model = PandasTableModel(self.dataframe,suppress_columns=["i_seqs"])
      # self.table = table
      # self.table._parent_explicit = self
    if geometry_ref is not None:
      self.parent.view.toggle_tab_visible(self.title,show=True)
    
    # update filter
    self.filter.update_quiet()



  # A generic pandas helper function to match a row_dict with an actual row in a dataframe
  @staticmethod 
  def find_matching_row(df, dict_subset):
    # hack to match keys
    new_d = {}
    for key,value in dict_subset.items():
      if "Atom" in key:
        i = key[-1]
        new_d[f"atom_id_{i}"] = value
      elif key.lower() in df.columns:
        new_d[key.lower()] = value
      elif key in df.columns:
        new_d[key] = value
      else:
        assert False, "Unable to match key between named tuple and actual df"

    dict_subset = new_d
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
      selection = Selection.from_i_seqs(mol.sites,[i])
      sel_strings.append(selection.phenix_string)
    return sel_strings


  @staticmethod
  def get_atom_records_from_iseqs(sites,i_seqs):
    # return a list of dictionaries (atom record)
    atom_records= []
    for i in i_seqs:
      record = sites.select_from_i_seqs([i]).to_records_compositional()[0]
      atom_records.append(record)
    return atom_records

  @staticmethod
  def get_labels_compositional_from_iseqs(sites,i_seqs):
    # return a list of string labels which concat items in atom record
    labels= []
    for i in i_seqs:
      label = sites.select_from_i_seqs([i]).to_labels_compositional()[0]
      labels.append(label)
    return labels

class BondTableController(GeometryTableTabController):
  title="Bonds"
  geometry_name = "bond"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = BondEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      sites = self.state.mol.sites
      i_seqs = row_dict["i_seqs"]
      labels_compositional = self.get_labels_compositional_from_iseqs(sites,i_seqs)
      row = BondEdit(
        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        labels_compositional = labels_compositional,
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        action="mod"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,BondEditsRef):
          if ref.geometry_ref == self.geometry_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = BondEditsRef(data=objframe,geometry_ref=self.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class AngleTableController(GeometryTableTabController):
  title = "Angles"
  geometry_name = "angle"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = AngleEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      sites = self.state.mol.sites
      i_seqs = row_dict["i_seqs"]
      labels_compositional = self.get_labels_compositional_from_iseqs(sites,i_seqs)
      row = AngleEdit(
        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        labels_compositional = labels_compositional,

        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        action="mod"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,AngleEditsRef):
          if ref.geometry_ref == self.geometry_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = AngleEditsRef(data=objframe,geometry_ref=self.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class DihedralTableController(GeometryTableTabController):
  title = "Dihedrals"
  geometry_name = "dihedral"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = DihedralEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      sites = self.state.mol.sites
      i_seqs = row_dict["i_seqs"]
      labels_compositional = self.get_labels_compositional_from_iseqs(sites,i_seqs)

      row = DihedralEdit(

        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        labels_compositional = labels_compositional,

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
          if ref.geometry_ref == self.geometry_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = DihedralEditsRef(data=objframe,geometry_ref=self.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

class ChiralTableController(GeometryTableTabController):
  title = "Chirals"
  geometry_name = "chirality"

class PlanarityTableController(GeometryTableTabController):
  title = "Planes"
  geometry_name = "plane"

class NonbondedTableController(GeometryTableTabController):
  title = "Non-bonded"
  geometry_name = "nonbonded"
  
  def addEdit(self,row_dict):
    # Create edit considering it as a bond
    df = getattr(self.geometry_ref.data,self.geometry_name)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    # calc distance
    # TODO: for now just standard
    row_dict["ideal"] = 3
    row_dict["sigma"] = 0.01
    dialog = BondEditDialog(defaults_dict=row_dict)
    if dialog.exec_():
      sites = self.state.mol.sites
      i_seqs = row_dict["i_seqs"]
      labels_compositional = self.get_labels_compositional_from_iseqs(sites,i_seqs)
      row = BondEdit(
        i_seqs = row_dict["i_seqs"],
        sel_strings = self.get_sel_strings_from_iseqs(self.state.mol,row_dict["i_seqs"]),
        labels_compositional = labels_compositional,
        ideal_old=row_dict["ideal"],
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=row_dict["sigma"],
        sigma_new = dialog.collectInputValues()["sigma"],
        action="add"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,BondEditsRef):
          if ref.geometry_ref == self.geometry_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = BondEditsRef(data=objframe,geometry_ref=self.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")
