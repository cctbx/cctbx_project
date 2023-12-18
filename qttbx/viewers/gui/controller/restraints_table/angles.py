from PySide2 import QtCore

from ...view.widgets import  FastTableView, PandasTableModel
from ...state.ref import SelectionRef
from ..controller import Controller



class AnglesTableTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Signals
    self.view.table.mouseReleased.connect(self.on_mouse_released)
    self.state.signals.restraints_change.connect(self.update)



  def update(self,restraint_ref):
    model = PandasTableModel(restraint_ref.dfs["angle"])
    self.view.table.setModel(model)


  def on_selection_changed(self, selected, deselected):
    df_sel = self.view.table.selected_rows()
    if len(df_sel)==1:
      if df_sel is not None:
        if "i_seqs" in df_sel.columns:
          # switch to atom picking level
          self.state.signals.picking_level.emit(1) # element

          flattened_i_seqs = [item for sublist in df_sel['i_seqs'] for item in sublist]
          flattened_i_seqs = [int(e) for e in flattened_i_seqs]
          # TODO: verify that .loc[i_seqs] will in fact work
          sel_sites = self.state.mol.sites.__class__(self.state.mol.sites.loc[flattened_i_seqs],attrs_map=self.state.mol.sites.attrs_map)
          query = self.state.mol.sites._convert_sites_to_query(sel_sites)
          ref = SelectionRef(data=query,model_ref=self.state.active_model_ref,show=False)
          self.state.add_ref(ref)
          self.state.active_selection_ref = ref
      else:
        print("no atoms returned as query")


  def on_mouse_released(self):
    selected = self.view.table.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)
