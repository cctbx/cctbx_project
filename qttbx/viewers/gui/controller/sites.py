from PySide2 import QtCore

from ..view.widgets import  PandasTableModel
from ..state.ref import SelectionRef
from .controller import Controller


class SitesTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Signals
    self.view.table.mouseReleased.connect(self.on_mouse_released)
    self.state.signals.model_change.connect(self.update)

  
  def update(self,*args):
    if self.state.mol is not None:
      model = PandasTableModel(self.state.mol.sites)
      self.view.table.setModel(model)


  def on_selection_changed(self, selected, deselected):
    print("on_selection_changed_atoms()")
    df_sel = self.view.table.selected_rows()
    if df_sel is not None:
      # switch to atom picking level
      self.view.parent_explicit.viewer_tab_view.selection_controls.combo_box.setCurrentIndex(1)
      sel = df_sel.index
      sel_sites = self.state.mol.sites.__class__(self.state.mol.sites.loc[sel],attrs_map=self.state.mol.sites.attrs_map)
      query = self.state.mol.sites._convert_sites_to_query(sel_sites)
      ref = SelectionRef(data=query,model_ref=self.state.active_model_ref,show=False)
      self.state.add_ref(ref)
      self.state.active_selection_ref = ref
      self.state.signals.select.emit(ref)
    
    else:
      print("no atoms returned as query")
   
    
  def on_mouse_released(self):
    print("on_mouse_released_atoms()")
    selected = self.view.table.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)
      