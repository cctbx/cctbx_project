import pandas as pd
from PySide2 import QtCore
from PySide2.QtCore import QAbstractTableModel, Qt
from PySide2.QtGui import QBrush, QColor
from PySide2.QtWidgets import  QStyledItemDelegate, QAbstractItemView

from .controller import Controller
from ..view.table import CustomDelegate
from ..state.ref import SelectionRef
from ..state.base import ObjectFrame
from ..state.color import Color
from ...core.selection import Selection

"""
Controller for Pandas table generally
"""

class TableController(Controller):
  title = "Table"
  row_class = None
  editable_columns = [] # if not empty, enforced as exclusive
  non_editable_columns = [] # excluded regardless of presence in editable list


  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._dataframe = None
    self._table_model = None
    self._objectframe = None
    self.delegate = None # For changing cosmetics of cells
    self.was_modified = False

    # Signals
    self.view.table_view.mouseReleased.connect(self.on_mouse_released)
    self.view.table_view.editRequested.connect(self.on_edit_requested)



  @property
  def dataframe(self):
    if self._dataframe is None:
      if self._objectframe is not None:
        self._dataframe = self.objectframe.to_df()
    return self._dataframe

  @dataframe.setter
  def dataframe(self,value):
    self._dataframe = value

    # adjust column widths. TODO: Move elsewhere... view?
    self.view.table_view.adjustColumnWidths()
    print("Showing tab due to dataframe setter: ",self.title)
    self.parent.view.toggle_tab_visible(self.title,show=True)

  @property
  def objectframe(self):
    if self._objectframe is None:
      if self.dataframe is not None:
        self._objectframe = ObjectFrame.from_df(self.dataframe,self.row_class)
    return self._objectframe

  @objectframe.setter
  def objectframe(self,value):
    self._objectframe = value

  @property
  def table_model(self):
    return self._table_model

  @table_model.setter
  def table_model(self,value):
    self._table_model = value
    self.view.table_view.setModel(value)
 
    self.reset_delegate()
    self.update_delegate()

  def column_name_exists(self, column_name):
    for i in range(self.table_model.columnCount(None)):
      if self.table_model.headerData(i, Qt.Horizontal, Qt.DisplayRole) == column_name:
        return True
    return False
  

  def on_edit_requested(self,index):
    df = self.table_model.df
    idx = index.row()
    specific_row = next(df.iloc[[idx]].itertuples(index=False), None)
    row_dict = specific_row._asdict()
    dialog = self.view.table_view.edit_dialog(
      row_dict=row_dict,
      edit_allowed_items=self.editable_columns,
      edit_denied_items=self.non_editable_columns)

    if dialog.exec_():
      # Dialog exited, update the dataframe
      new_data = dialog.collectInputValues()
      # Check if a change will occur
      modified_fields = dialog.getModifiedFields()
      print("modified fields:",modified_fields)

      if len(modified_fields)>0:

        # get modified cols
        cols = list(df.columns)
        col_idxs = [cols.index(col) for col in modified_fields.keys()]

        # Get initial
        new_data = df.loc[idx].to_dict()


        # Set new data
        new_data.update(modified_fields)
        df.loc[idx] = new_data

        # Alert that data has changed
        for col_idx in col_idxs:
          color = Color.from_string("powderblue")
          self.table_model.color_map[(idx,col_idx)] = color
          self.table_model.font_map[(idx,col_idx)] = {"bold":True,"italic":True}
          #self.delegate.update_color_map(idx,col_idx,color)
          #self.delegate.update_font_map(idx,col_idx,bold=True,italic=True)
        self.view.set_was_modified(True)
        self.table_model.was_modified = True
        self.update_delegate()

  def reset_delegate(self):
    self.view.table_view.setItemDelegate(QStyledItemDelegate(self.view.table_view))
    self.delegate = self.view.table_view.itemDelegate()

  def update_delegate(self):
    view = self.view.table_view
    delegate = CustomDelegate(parent=view)
    # TODO: Why manual reset necessary?
    delegate.color_map = {}
    delegate.font_map = {}
    delegate.dtype_map = CustomDelegate.get_dtype_map(self.table_model.df)
    view.setItemDelegate(delegate)
    view.setSelectionBehavior(QAbstractItemView.SelectRows)
    self.delegate = delegate

    # update colors
    for (col,row),color in self.table_model.color_map.items():
      self.delegate.update_color_map(col,row,color)

    # update font
    for (col,row),fonts in self.table_model.font_map.items():
      self.delegate.update_font_map(col,row,**fonts)



  def select_row(self, row_index):
    selection_model = self.view.table_view.selectionModel()
    index = self.table_model.index(row_index, 0)  # Select the entire row
    selection_model.select(index, selection_model.Select | selection_model.Rows)


  def on_mouse_released(self):
    selected = self.view.table_view.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)


  def on_selection_changed(self, selected, deselected):
    # Send a selection signal out to focus in graphics
    df_sel = self.view.table_view.selected_rows()
    if len(df_sel)==1:
      if df_sel is not None:
        # switch to atom picking level
        self.state.signals.picking_level.emit("atom")
        flattened_i_seqs = [item for sublist in df_sel['i_seqs'] for item in sublist]
        flattened_i_seqs = [int(e) for e in flattened_i_seqs if pd.notna(e)]
        selection = Selection.from_i_seqs(self.state.mol.sites,flattened_i_seqs)
        ref = SelectionRef(data=selection,model_ref=self.state.active_model_ref,show=False)
        self.state.add_ref(ref)
        self.state.active_selection_ref = ref
      else:
        print("no atoms returned as query")