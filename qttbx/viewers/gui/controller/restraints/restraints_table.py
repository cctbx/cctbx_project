
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2 import QtCore
from ..controller import Controller
from ...state.restraints import Restraints
from ...state.ref import RestraintsRef
from ...view.widgets import  PandasTable
from ...state.ref import SelectionRef
from .restraints_staging import RestraintStagingListController

from pathlib import Path
import pandas as pd



class RestraintTableTabController(Controller):
  restraint_name = None # Subclass and fill this in (bond,angle,etc)
  columns_to_include= ['ideal','model',"sigma","delta",'residual',"vdw","action"]
  column_prefixes_to_include = ["atom_id"]
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    # Signals
    self.view.table.mouseReleased.connect(self.on_mouse_released)
    self.state.signals.restraints_change.connect(self.update)

  def update(self,restraint_ref):
    print("updating restraints table with ref:",restraint_ref.id)
    if hasattr(restraint_ref.data,self.restraint_name):
      df = getattr(restraint_ref.data,self.restraint_name)
      cols = ["Chain","Res","Seq"]
      # add prefix columns
      for col in df.columns:
        for prefix in self.column_prefixes_to_include:
          if col.startswith(prefix):
            cols.append(col)
      # add exact match columns
      cols+= [col for col in self.columns_to_include if col in df.columns]
      # add i_seqs (list of all i_seqs)
      cols+= ["i_seqs"]
      # TODO: reorder cols without copy? How.
      table = PandasTable(df[cols].copy(),suppress_columns=["i_seqs"])
      self.table = table
      self.table._parent_explicit = self


  def on_selection_changed(self, selected, deselected):
    df_sel = self.view.table.selected_rows()
    if len(df_sel)==1:
      if df_sel is not None:
        # switch to atom picking level
        self.state.signals.picking_level.emit(1) # element
        flattened_i_seqs = [item for sublist in df_sel['i_seqs'] for item in sublist]
        flattened_i_seqs = [int(e) for e in flattened_i_seqs if pd.notna(e)]
        sel_sites = self.state.mol.sites.iloc[flattened_i_seqs]
        query = self.state.mol.sites._convert_sites_to_query(sel_sites)
        ref = SelectionRef(data=query,model_ref=self.state.active_model_ref,show=False)
        self.state.add_ref(ref)
        self.state.active_selection_ref = ref
      else:
        print("no atoms returned as query")

  @property
  def table(self):
    return self.view.table.model()

  @table.setter
  def table(self,value):
    self.view.table.setModel(value)
    # adjust column widths. TODO: Move to table classes
    self.view.table.adjustColumnWidths()

  def on_mouse_released(self):
    selected = self.view.table.selectionModel().selection()
    deselected = QtCore.QItemSelection()
    self.on_selection_changed(selected, deselected)


class BondTableController(RestraintTableTabController):
  restraint_name = "bond"

class AngleTableController(RestraintTableTabController):
  restraint_name = "angle"

class DihedralTableController(RestraintTableTabController):
  restraint_name = "dihedral"

class ChiralTableController(RestraintTableTabController):
  restraint_name = "chirality"

class PlanarityTableController(RestraintTableTabController):
  restraint_name = "plane"

class NonbondedTableController(RestraintTableTabController):
  restraint_name = "nonbonded"

class RestraintsTableTopTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    #self.staging = RestraintStagingListController(parent=self,view=self.view.staging.list_view)
    self.bonds = BondTableController(parent=self,view=view.bonds)
    self.angles = AngleTableController(parent=self,view=view.angles)
    self.dihedrals = DihedralTableController(parent=self,view=view.dihedrals)
    self.chirals = ChiralTableController(parent=self,view=view.chirals)
    self.planes = PlanarityTableController(parent=self,view=view.planes)
    self.nonbonded = NonbondedTableController(parent=self,view=view.nonbonded)

    # for widget in self.view.widgets:
    #   widget.process_button.clicked.connect(self.process_model)

  #   # Signals
  #   self.state.signals.restraints_change.connect(self.handle_restraints_change)

  #   # Flags
  #   self.header_hidden = True
  #   self._hide_header()

  # def handle_restraints_change(self):
  #   if self.state.active_model_ref.has_restraints:
  #     if not self.header_hidden:
  #       self._hide_header()
  #   else:
  #     if self.header_hidden:
  #       self._show_header()


  # def process_model(self):
  #   if self.state.active_model_ref is None:
  #     msg = QMessageBox()
  #     msg.setWindowTitle("Notification")
  #     msg.setText("Select an active model before generating restraints.")
  #     msg.setIcon(QMessageBox.Information)
  #     msg.setStandardButtons(QMessageBox.Ok)
  #     msg.exec_()
  #   else:
  #     if not self.state.active_model_ref.has_restraints:
  #       ref = self.state.active_model_ref
  #       model = ref.model
  #       model.add_crystal_symmetry_if_necessary()
  #       model.process(make_restraints=True)
  #       if model.get_restraints_manager() is not None:
  #         grm = model.get_restraints_manager().geometry
  #       else:
  #         msg = QMessageBox()
  #         msg.setWindowTitle("Notification")
  #         msg.setText("Failed to process model and make restraints...")
  #         msg.setIcon(QMessageBox.Information)
  #         msg.setStandardButtons(QMessageBox.Ok)
  #         msg.exec_()
  #         return


  #       # write out geo file
  #       geo_file_dir = Path.cwd()
  #       if ref.data.filepath is not None:
  #         geo_file_dir = Path(ref.data.filepath).parent
  #       geo_file_path = geo_file_dir / Path(ref.label+".geo")
  #       grm.write_geo_file(sites_cart=model.get_sites_cart(),
  #                          site_labels = model.get_xray_structure().scatterers().extract_labels(),
  #                          file_name=str(geo_file_path))

  #       # notify
  #       msg = QMessageBox()
  #       msg.setWindowTitle("Notification")
  #       msg.setText("Restraints written to .geo file:\n"+str(geo_file_path))
  #       msg.setIcon(QMessageBox.Information)
  #       msg.setStandardButtons(QMessageBox.Ok)
  #       msg.exec_()

  # def load_from_geo_file(self,geo_file_path):
  #   restraints = Restraints.from_geo_file(geo_file_path)
  #   ref = RestraintsRef(restraints,self.state.active_model_ref)


  # def _hide_header(self):
  #   for widget in self.view.widgets:
  #     self._hide_child_layout(widget.layout, widget.header_layout)
  #     QApplication.processEvents()  # Update the UI
  #     self.header_hidden = True

  # def _show_header(self):
  #   for widget in self.view.widgets:
  #     self._show_child_layout(widget.layout, widget.header_layout)
  #     QApplication.processEvents()  # Update the UI
  #     self.header_hidden = False

  # def _hide_child_layout(self, layout_parent, layout_child):
  #   layout_parent.removeItem(layout_child)
  #   for i in range(layout_child.count()):
  #     widget = layout_child.itemAt(i).widget()
  #     if widget is not None:
  #       widget.hide()

  # def _show_child_layout(self, layout_parent, layout_child):
  #   for i in range(layout_child.count()):
  #     widget = layout_child.itemAt(i).widget()
  #     if widget is not None:
  #       widget.show()
    # Assuming the header layout should be inserted at position 0
    # layout_parent.insertLayout(0, layout_child)


