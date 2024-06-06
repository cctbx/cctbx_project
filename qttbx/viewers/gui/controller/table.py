import pandas as pd
from PySide2 import QtCore

from .controller import Controller
from ..state.ref import SelectionRef
from ..state.table import PandasTableModel
from ..state.base import ObjectFrame
from ...core.selection import Selection

"""
Controller for Pandas table generally
"""

class TableController(Controller):
  title = "Table"
  row_class = None
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._dataframe = None
    self._table_model = None
    self._objectframe = None

    # Signals
    self.view.table_view.mouseReleased.connect(self.on_mouse_released)


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