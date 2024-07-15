
import pandas as pd
from PySide2.QtCore import Qt

from ..filter import TableFilterController
from ..table import TableController
from ...state.table import PandasTableModel
from ....core.selection import Selection
from ...state.base import ObjectFrame
from ...view.widgets import (
  BondEditDialog,
  AngleEditDialog,
  DihedralEditDialog,
  ChiralEditDialog
)
from ...state.edits import (
  BondEdit,
  AngleEdit,
  DihedralEdit,
  ChiralEdit
)
from ...state.ref import (
  BondEditsRef,
  AngleEditsRef,
  DihedralEditsRef
)

from ..filter import CompositeFilter


class GeometryTableTabController(TableController):
  """
  Controller for a Pandas Table in geometry tab
  """
  geometry_type = None # Subclass and fill this in (bond,angle,etc)

  # Define display cols (or col prefixes) and specify order
  display_columns = ['Label','ideal','model',"sigma","delta",'residual'] 

  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.geometry_ref = None
    self.col_names = None
    self.debug_flag = False

    # Filter panel
    self.filter = TableFilterController(parent=self,view=self.view.filter)

    # Signals
    self.view.table_view.addEdit.connect(self.addEdit)
    self.view.edit_controls.edit_button.clicked.connect(self.on_edit_pressed)
    self.state.signals.geometry_change.connect(self.update)
    self.view.table_view.horizontalHeader().sortIndicatorChanged.connect(self.on_sort_indicator_changed)

      
    #self.state.signals.filter_update.connect(self.update_from_filter)

  # def get_suppress_columns(self,df):
  #   suppress_columns = []
  #   for col in self.dataframe.columns:
  #     for name in self.column_prefixes:
  #       added = False
  #       if not added:
  #         if self.split_cif_key_prefix(col) == name:
  #           suppress_columns.append(col)
  #           added = True

      
  #   suppress_columns = ["i_seqs"] + suppress_columns
  #   return suppress_columns
  @property
  def selected(self):
    return self.view.table_view.selected_rows()

  def on_edit_pressed(self):
    sel = self.selected
    if len(sel)==1:
      row_dict = sel.to_dict("records")[0]
      self.addEdit(row_dict)

  def addEdit(self,row_dict):
    raise NotImplementedError

  # def update_from_filter(self,filter_obj,debug_flag):
  #   df = self.dataframe
  #   if df is None:
  #     return
  #   #comp = filter_obj.comp_id
    
  #   # # Determine the quantile based on the direction
  #   # if direction == "Ascending":
  #   #     quantile = (float(filter.current_percentile) / 100)
  #   # else:  # Descending
  #   #     quantile = 1 - (float(filter.current_percentile) / 100)

  #   # Filter by component
  #   # if comp == "All":
  #   #     comp_test = pd.Series([True] * len(df))
  #   # else:
  #   #     comp_test = df["Res"] == comp

  #   # # Apply the quantile filtering
  #   # if metric in df.columns:
  #   #     percentile_cutoff = df[metric].quantile(quantile)
  #   #     if direction == "Ascending":
  #   #         percentile_test = df[metric] <= percentile_cutoff
  #   #     else:  # Descending
  #   #         percentile_test = df[metric] >= percentile_cutoff
  #   # else:
  #   #     percentile_test = pd.Series([True] * len(df))

  #   # Apply both conditions to filter the DataFrame

  #   df_filtered = filter_obj.filter_df(df)
  #   df_filtered = df_filtered.rename(columns=self.rename_columns)
  #   suppress_cols = [c.lower() for c in self.suppress_columns]
  #   suppress_cols+= [c.capitalize() for c in self.suppress_columns]
  #   table_model = PandasTableModel(df_filtered,suppress_columns=suppress_cols)
  #   self.table_model = table_model
  #   #self.filter_obj = filter_obj

  # @property
  # def column_prefixes(self):
  #   # Columns that are maintained in the geometry dataframes
  #   return list(self.state.params.core_map_to_mmcif.keys()) + ["residue_class"]

  # def reform_columns(self,df):
  #   cols = ["i_seqs"]
  #   #prefixes =  ["Label"]

  #   for prefix in self.column_prefixes:
  #     for col in df.columns:
  #       if col.startswith(prefix):
  #         if not col.startswith("id_str"):
  #           cols.append(col)

  #   cols+= [col for col in self.column_order if col in df.columns]
  #   df = df[cols].copy()
    # # add prefix columns
    # for prefix in prefixes:
    #   for col in cols:
    #     if col.startswith(prefix):
    #       number = col.split(prefix)[-1]
    #       number = number.strip("_")
    #       col_name = f"Atom{number}"
    #       #self.rename_columns[col] = col_name

    # df.rename(columns=self.rename_columns,inplace=True)
    

    # # capitalize
    # capitalize = {c:c.capitalize() for c in df.columns if c not in ["i_seqs"]}
    #df.rename(columns=self.rename_columns,inplace=True)

    self.col_names = list(df.columns)
    return df

  def update(self,geometry_ref):
    self.geometry_ref = geometry_ref
    if hasattr(geometry_ref.data,self.geometry_type):
      df = getattr(geometry_ref.data,self.geometry_type)
      self.dataframe = df#self.reform_columns(df)
      #suppress_columns = self.get_suppress_columns(df)
      self.table_model = PandasTableModel(self.dataframe,display_columns=self.display_columns,capitalize=True,remove_underscores=True)
      # self.table = table
      # self.table._parent_explicit = self
    if geometry_ref is not None:
      self.parent.view.toggle_tab_visible(self.title,show=True)
    
    # update filter
    self.filter.update_quiet()
    self.initialize_sort()

  def initialize_sort(self):
    # Initialize sorting by residual
    if self.table_model:
      col_name = None
      if self.column_name_exists("Residual"):
        col_name = "Residual"
      elif self.column_name_exists("Residual 1"):
        col_name = "Residual 1"
      if col_name:
        self.sort_table_by_column_name(col_name,Qt.DescendingOrder)

  # Sorting machineray. Candidate to move up to superclass
  def on_sort_indicator_changed(self,logicalIndex, order):
    #print(f"Sorting changed. Column: {logicalIndex}, Order: {'Ascending' if order == Qt.AscendingOrder else 'Descending'}")
    sorted_by = self.table_model.headerData(logicalIndex, Qt.Horizontal, Qt.DisplayRole)
    order = 'Ascending' if order == Qt.AscendingOrder else 'Descending'
    self.view.sort_label.setText(f"Sorted by: {sorted_by}, Order: {order}")

  def sort_table_by_column_name(self,column_name, order):
    model = self.table_model
    view = self.view.table_view
    logical_index = None
    for i in range(model.columnCount(None)):
      if model.headerData(i, Qt.Horizontal, Qt.DisplayRole) == column_name:
        logical_index = i
        break
    if logical_index is not None:
      view.sortByColumn(logical_index, order)
    else:
      print(f"Column '{column_name}' not found in the model")

  # A generic pandas helper function to match a row_dict with an actual row in a dataframe
  @staticmethod
  def find_matching_row(df, dict_subset):
    # Start with a mask of all True values
    mask = pd.Series([True] * len(df))

    for key, value in dict_subset.items():
      if value is None:
        # Check for NaN values in the DataFrame
        mask &= df[key].isnull()
      else:
        # Check for matching values in the DataFrame
        mask &= df[key] == value

    # Find rows where all specified dictionary key-value pairs match
    matching_rows = df[mask]

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
  geometry_type = "bond"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_type)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = BondEditDialog(defaults_dict=row_dict,action="mod")
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
  geometry_type = "angle"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_type)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = AngleEditDialog(defaults_dict=row_dict,action="mod")
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
  geometry_type = "dihedral"

  def addEdit(self,row_dict):
    df = getattr(self.geometry_ref.data,self.geometry_type)
    #import pdb
    #pdb.set_trace()
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    dialog = DihedralEditDialog(defaults_dict=row_dict,action="mod")
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
  geometry_type = "chirality"

class PlanarityTableController(GeometryTableTabController):
  title = "Planes"
  geometry_type = "plane"

class NonbondedTableController(GeometryTableTabController):
  title = "Non-bonded"
  geometry_type = "nonbonded"
  
  def addEdit(self,row_dict):
    # Create edit considering it as a bond
    df = getattr(self.geometry_ref.data,self.geometry_type)
    row_dict = self.find_matching_row(df,row_dict).to_dict() # get full row from geometry
    # calc distance
    # TODO: for now just standard
    row_dict["ideal"] = 3
    row_dict["sigma"] = 0.01
    dialog = BondEditDialog(defaults_dict=row_dict,action="add")
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
