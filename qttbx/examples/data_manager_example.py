"""qttbx DataManagerWidget -- end-to-end runnable example.

Run:
  /path/to/build/bin/libtbx.python qttbx/examples/data_manager_example.py

A window opens with:

  * a form of basic (non-file) PHIL parameters on the left,
  * the DataManagerWidget on the right -- drag files into it from the
    file manager, or use the "Add files..." button.  Files are added
    to the DataManager and can be bound to compatible PHIL path
    parameters via the per-row "+ add" popup.
  * a "Live PHIL extract" view at the bottom -- click Refresh to see
    the current PHIL state after edits and bindings.

The DataManager itself is exposed as ``window.widget.data_manager``
so the example can be driven programmatically as well.
"""
import sys

from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
  QApplication, QFormLayout, QLabel, QMainWindow,
  QPushButton, QSplitter, QTextEdit, QVBoxLayout, QWidget,
)

import iotbx.phil

from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
from qttbx.widgets.data_manager import DataManagerWidget, DataManagerTableModel
from qttbx.examples._helpers import ProportionalColumns


MASTER_PHIL = """
refinement {
  cycles = 5
    .type = int(value_min=1, value_max=20)
    .short_caption = "Refinement cycles"
  weight = 1.0
    .type = float(value_min=0.0, value_max=10.0)
    .short_caption = "Weight"
  use_ncs = False
    .type = bool
    .short_caption = "Use NCS"
  mode = fast *thorough exhaustive
    .type = choice
    .short_caption = "Refinement mode"
}

data {
  input_model = None
    .type = path
    .style = file_type:pdb
    .short_caption = "Input model"
  reference_model = None
    .type = path
    .style = file_type:pdb
    .short_caption = "Reference model"
  input_data = None
    .type = path
    .style = file_type:hkl
    .short_caption = "Input data"
  output_map = None
    .multiple = True
    .type = path
    .style = file_type:ccp4_map
    .short_caption = "Output map"
}
"""


def make_basic_form(model):
  """Build a form of basic (non-file) PHIL parameters."""
  form = QFormLayout()
  for caption, path in (
      ("Cycles", "refinement.cycles"),
      ("Weight", "refinement.weight"),
      ("Use NCS", "refinement.use_ncs"),
      ("Mode", "refinement.mode"),
  ):
    form.addRow(caption, PhilField(model, path))
  container = QWidget()
  container.setLayout(form)
  return container


def make_extract_view(model):
  """Build a read-only PHIL-extract pane with a Refresh button."""
  container = QWidget()
  layout = QVBoxLayout(container)
  text = QTextEdit()
  text.setReadOnly(True)
  text.setLineWrapMode(QTextEdit.NoWrap)
  refresh_btn = QPushButton("Refresh")
  layout.addWidget(QLabel("Live PHIL extract:"))
  layout.addWidget(text, stretch=1)
  layout.addWidget(refresh_btn)

  def _refresh():
    working = model.get_working_phil(diff_only=False)
    text.setPlainText(working.as_str())

  refresh_btn.clicked.connect(_refresh)
  _refresh()
  return container, _refresh


class ExampleWindow(QMainWindow):
  """Main window: basic-form + DataManagerWidget + live extract."""

  def __init__(self, parent=None):
    QMainWindow.__init__(self, parent)
    self.setWindowTitle("qttbx DataManagerWidget example")
    self.resize(1100, 700)

    master = iotbx.phil.parse(MASTER_PHIL)
    self._model = PhilModel()
    self._model.initialize_model(master)

    self.widget = DataManagerWidget(phil_model=self._model)

    # Column proportions: Filename 60%, Type 20%, Used for 20%; Delete
    # keeps its widget-default fixed width. The helper installs an
    # event filter on the table that reapplies the proportions on
    # every viewport-resize.
    self._column_proportions = ProportionalColumns(
      self.widget._table,
      [
        (DataManagerTableModel.COL_FILENAME, 0.60),
        (DataManagerTableModel.COL_TYPE,     0.20),
        (DataManagerTableModel.COL_USED_FOR, 0.20),
      ])

    extract_view, refresh = make_extract_view(self._model)
    self._refresh_extract = refresh

    top = QSplitter(Qt.Horizontal)
    top.addWidget(make_basic_form(self._model))
    top.addWidget(self.widget)
    top.setSizes([350, 750])

    root = QSplitter(Qt.Vertical)
    root.addWidget(top)
    root.addWidget(extract_view)
    root.setSizes([500, 200])

    self.setCentralWidget(root)

    # Refresh the extract view whenever bindings change so the user
    # can watch PHIL state update without clicking Refresh.
    self.widget.bindingChanged.connect(lambda *_a: self._refresh_extract())
    self.widget.fileAdded.connect(lambda *_a: self._refresh_extract())
    self.widget.fileRemoved.connect(lambda *_a: self._refresh_extract())

def main():
  # Opt into Qt's high-DPI scaling before constructing QApplication so
  # the widget's hardcoded pixel constants (chip paddings, column
  # widths, border widths) get scaled by the system's device pixel
  # ratio. Font sizes already scale via QFontMetrics. These attributes
  # are application-global and must be set before any QApplication.
  if QApplication.instance() is None:
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)
  app = QApplication.instance() or QApplication(sys.argv)
  win = ExampleWindow()
  win.show()
  sys.exit(app.exec_())


if __name__ == "__main__":
  main()
