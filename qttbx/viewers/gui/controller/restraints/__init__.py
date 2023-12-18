
from PySide2.QtWidgets import QApplication, QMessageBox

from PySide2.QtCore import QUrl, QThread, Signal, Slot, QObject, QThreadPool, QRunnable
from ..scroll_entry import ScrollEntryController
from ...view.widgets.scroll_entry import ScrollEntryView
from ..scroll_list import ScrollableListController
from ..controller import Controller
from ...state.restraints import Restraints
from ...state.ref import RestraintsRef, RestraintRef
from ..restraints.bonds import BondTabController

import time

class RestraintsTopTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.bonds = BondTabController(parent=self,view=view.bonds)

    for widget in self.view.widgets:
      widget.process_button.clicked.connect(self.process_model)

    # Signals for the presence of restraints
    self.state.signals.restraints_change.connect(self.handle_restraints_change)

    # Flags
    self.header_hidden = False

  def handle_restraints_change(self):
    if self.state.active_model_ref.has_restraints:
      if not self.header_hidden:
        self._hide_header()
    else:
      if self.header_hidden:
        self._show_header()

  def process_model(self):
    if self.state.active_model_ref is None:
      msg = QMessageBox()
      msg.setWindowTitle("Notification")
      msg.setText("Select an active model before generating restraints.")
      msg.setIcon(QMessageBox.Information)
      msg.setStandardButtons(QMessageBox.Ok)
      msg.exec_()
    else:
      if not self.state.active_model_ref.has_restraints:
        ref = self.state.active_model_ref
        model = ref.model
        model.process(make_restraints=True)
        restraints_ref = RestraintsRef.from_model_ref(self.state.active_model_ref)
        ref.restraints = restraints_ref
        self.state.signals.restraints_change.emit(restraints_ref)



  def _hide_header(self):
    for widget in self.view.widgets:
      self._hide_child_layout(widget.layout, widget.header_layout)
      QApplication.processEvents()  # Update the UI
      self.header_hidden = True

  def _show_header(self):
    for widget in self.view.widgets:
      self._show_child_layout(widget.layout, widget.header_layout)
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
