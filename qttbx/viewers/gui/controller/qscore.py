import pandas as pd

from PySide2.QtWidgets import QApplication
from PySide2 import QtCore

from cctbx.maptbx.qscore import calc_qscore
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from cctbx.programs.qscore import Program as QscoreProgram
from iotbx.cli_parser import run_program
from libtbx import phil
import numpy as np

from .controller import Controller
from ..state.results import  Result, ResultsRef
from ..state.table import PandasTableModel as PandasTable
from ..state.ref import SelectionRef

class QscoreTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    # Signals
    self.view.table_view.mouseReleased.connect(self.on_mouse_released)
    self.view.process_button.clicked.connect(self.calculate_qscore)
    self.state.signals.results_change.connect(self.update)
    self.state.signals.model_change.connect(self.handle_model_change)
    self.view.combobox.currentIndexChanged.connect(self.update)
    #self.view.hist_widget.emitter.histogram_click_value.connect(self.handle_hist_click)

    # Flags
    self.header_hidden = False

    # Qscore params
    working_phil = phil.parse(QscoreProgram.master_phil_str,process_includes=True)
    params = working_phil.extract()
    self.params = params

  def handle_hist_click(self):
    pass


  def handle_model_change(self,model_ref):
    if model_ref is not None:
      if "qscore" in model_ref.results:
        self.update(model_ref)
      else:
        self._show_header()

  def update(self):
    print("update")
    model_ref = self.state.active_model_ref
    if "qscore" in model_ref.results:
      results_ref = model_ref.results["qscore"]
      df = results_ref.data.results.qscore_dataframe

      # TODO: Fix below. id and positional index not the same
      # TODO: Also, syntax should be more like:
      #        sites["Qscore"] = values
      #        df = sites.per_atom(keep=["Q-score"],sort=[],ascending=False)
      #        # which is...
      #        df = df.sort_values(by="Q-score",inplace=False,ascending=ascending)
      #        df = df[keep + sites.per_atom_cols]
      #        # or
      #        df = sites.per_residue(keep=["Q-score"],agg="first",sort=,ascending=)
      #        # which is...
      #        # df.groupby(sites.residue_cols).agg('first').reset_index())
      #        df = df.sort_values(by="Q-score",inplace=False,ascending=ascending)
      #        df = df[keep + sites.per_residue_cols]
      df["id"] = df.index.values+1
      df = pd.merge(df, model_ref.mol.sites, on='id', how='left')
      df = df.round(3)
      df.sort_values(by="Q-score",inplace=True,ascending=False)
      df = df[["Q-score","Q-Residue","chain_id","resseq","resname","name"]]

      # Check granularity
      index = self.view.combobox.currentIndex()
      label = self.view.combobox.itemText(index)
      print("label:",label)
      if label == "Q-Residues":
        df = df.groupby(["chain_id","resseq","resname"]).agg('first').reset_index()
        df = df[["Q-Residue","chain_id","resseq","resname"]]

      model = PandasTable(df)
      self.view.table_view.setModel(model)
      #self.view.add_histogram()


  def calculate_qscore(self):

    if self.state.active_model_ref is not None:

      dm = DataManager() #TODO: use state dm, not a temporary one
      dm.add_model(
        self.state.active_model_ref.id,
        self.state.active_model_ref.model)
      dm.add_real_map(
        self.state.active_map_ref.id,
        self.state.active_map_ref.map_manager)
      task = QscoreProgram(dm,self.params)
      task.run()
      results = task.get_results()

      # build result and ref
      qresult =  Result(program_name="qscore",
                        model_refs = [self.state.active_model_ref],
                        map_refs = [self.state.active_map_ref],
                        results=results,
                        params=self.params,
                        )
      result_ref = ResultsRef(data=qresult)
      self.state.add_ref(result_ref)


      # Hide header
      self._hide_header()






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
    df_sel = self.view.table_view.selected_rows()
    if df_sel is not None:
      df_sel = df_sel.drop(columns=["Q-score"])
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
    selected = self.view.table_view.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)
