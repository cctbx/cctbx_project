from PySide2 import QtCore

from ..state.table import PandasTableModel as PandasTable
from .controller import Controller
from ...core.atom_sites import AtomSites

class SitesTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Signals
    self.view.table_view.mouseReleased.connect(self.on_mouse_released)
    self.state.signals.model_change.connect(self.update)


  def update(self,*args):
    if self.state.mol is not None:
      model = PandasTable(self.state.mol.sites)
      self.view.table_view.setModel(model)


  def on_selection_changed(self, selected, deselected):
    df_sel = self.view.table_view.selected_rows()
    if df_sel is not None:
      # switch to atom picking level
      self.view.parent_explicit.viewer_tab_view.selection_controls.combo_box.setCurrentIndex(1)
      sel = df_sel.index
      sel_sites = AtomSites(self.state.mol.sites.loc[sel])
      query = self.state.mol.sites._convert_sites_to_query(sel_sites)
      ref = SelectionRef(data=query,model_ref=self.state.active_model_ref,show=False)
      self.state.add_ref(ref)
      self.state.active_selection_ref = ref

    else:
      print("no atoms returned as query")


  def on_mouse_released(self):
    selected = self.view.table_view.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)
