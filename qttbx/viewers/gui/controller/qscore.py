import pandas as pd

from PySide2.QtWidgets import QApplication
from PySide2 import QtCore

from .controller import Controller
from ..state.ref import  QscoreRef
from ..state.results import QscoreResult
from ..view.widgets import   PandasTableModel
from ..state.ref import SelectionRef

class QscoreTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    
    # Signals
    self.view.table.mouseReleased.connect(self.on_mouse_released)
    self.view.process_button.clicked.connect(self.calculate_qscore)
    self.state.signals.results_change.connect(self.update)
    self.state.signals.model_change.connect(self.handle_model_change)
    #self.view.hist_widget.emitter.histogram_click_value.connect(self.handle_hist_click)

    # Flags
    self.header_hidden = False

  def handle_hist_click(self):
    pass

  def handle_model_change(self,model_ref):
    if model_ref is not None:
      if "qscore" in model_ref.results:
        self.update(model_ref)
      else:
        self._show_header()
  
  def update(self):
    model_ref = self.state.active_model_ref
    if "qscore" in model_ref.results:
      results_ref = model_ref.results["qscore"]
      df = pd.DataFrame({"Qscore":results_ref.data.qscore_per_atom})
      df = df.round(3)
      df = pd.concat([df,model_ref.mol.sites],axis=1)
      df.sort_values(by="Qscore",inplace=True,ascending=False)
      
      model = PandasTableModel(df)
      self.view.table.setModel(model)
      #self.view.add_histogram()


  def calculate_qscore(self):
    # NOTE: This is temporary, all results should come through programs
    from ....qscore_standalone.qscore2 import qscore_np
    from iotbx.map_model_manager import map_model_manager
    if self.state.active_model_ref is not None:
      model = self.state.active_model_ref.model
      map_manager = self.state.active_map_ref.map_manager
      mmm = map_model_manager(model=model,map_manager=map_manager)
      q_np_v2 = qscore_np(mmm,n_probes=8,selection=None,version=2,nproc=8)

    # Hide header
    self._hide_header()

    # build result and ref
    qdata =  QscoreResult(program_name="qscore",qscore_per_atom=q_np_v2)
    result_ref = QscoreRef(data=qdata,model_ref = self.state.active_model_ref,selection_ref=None)
    self.state.add_ref(result_ref)

    
        

  def _hide_header(self):
    if not self.header_hidden:
      self._hide_child_layout(self.view.layout, self.view.header_layout)
      QApplication.processEvents()  # Update the UI
      self.header_hidden = True

  def _show_header(self):
    if self.header_hidden:
      self._show_child_layout(self.view.layout, self.view.header_layout)
      QApplication.processEvents()  # Update the UI
      self.header_hidden = False

  def _hide_child_layout(self, layout_parent, layout_child):
    layout_parent.removeItem(layout_child)
    for i in range(layout_child.count()):
      widget = layout_child.itemAt(i).widget()
      if widget is not None:
        widget.hide()

  def _show_child_layout(self, layout_parent, layout_child):
    for i in range(layout_child.count()):
      widget = layout_child.itemAt(i).widget()
      if widget is not None:
        widget.show()
    # Assuming the header layout should be inserted at position 0
    layout_parent.insertLayout(0, layout_child)



  def on_selection_changed(self, selected, deselected):
    df_sel = self.view.table.selected_rows()
    if df_sel is not None:
      df_sel = df_sel.drop(columns=["Qscore"])
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
    selected = self.view.table.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)