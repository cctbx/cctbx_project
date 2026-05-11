"""qttbx PHIL widgets -- end-to-end runnable example.

Run:
  /path/to/build/bin/libtbx.python qttbx/examples/phil_widgets_example.py

A window opens with a tree view of an example PHIL scope on the left and a
tab strip of form-style demos on the right. Both panes share the same
PhilModel; edits made in either propagate to the other automatically.

Every implemented widget is exercised at least once:

  Scalar         IntWidget, BoolWidget, FloatWidget, KeyWidget, ChoiceWidget,
                 ChoiceMultiWidget
  Text-shaped    StrWidget (+ StrTextWidget), QstrWidget, PathWidget
  Lists          StringsWidget (+ StringsTextWidget), WordsWidget
                 (+ WordsTextWidget), IntsWidget (+ IntsTextWidget),
                 FloatsWidget (+ FloatsTextWidget)
  Domain         SpaceGroupWidget, UnitCellWidget, AtomSelectionWidget
                 (+ AtomSelectionTextWidget)
  Wrappers       MultipleWidget (definition .multiple = True),
                 RepeatableScopeWidget (scope .multiple = True)

The "Live extract" tab shows the current PHIL settings serialized; click
Refresh to see edits land in the model.
"""
import sys

from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
  QApplication, QFormLayout, QLabel, QMainWindow, QPushButton, QSplitter,
  QTabWidget, QTextEdit, QTreeView, QVBoxLayout, QWidget,
)

import iotbx.phil

from qttbx.phil import PhilModel
from qttbx.widgets.phil import PhilField
from qttbx.widgets.phil.delegate import PhilItemDelegate
from qttbx.widgets.phil.multiple import RepeatableScopeWidget

# Long-text variants -- opt in via widget=.
from qttbx.widgets.phil.str_widget import StrTextWidget
from qttbx.widgets.phil.strings_widget import StringsTextWidget
from qttbx.widgets.phil.words_widget import WordsTextWidget
from qttbx.widgets.phil.ints_widget import IntsTextWidget
from qttbx.widgets.phil.floats_widget import FloatsTextWidget
from qttbx.widgets.phil.atom_selection import AtomSelectionTextWidget

# Domain widgets self-register on import.
import qttbx.widgets.phil.space_group     # noqa: F401  (registers "space_group")
import qttbx.widgets.phil.unit_cell       # noqa: F401  (registers "unit_cell")
import qttbx.widgets.phil.atom_selection  # noqa: F401  (registers "atom_selection")


MASTER_PHIL = """
basic {
  count = 5
    .type = int(value_min=1, value_max=20)
    .short_caption = "Refinement cycles"
  weight = 1.0
    .type = float(value_min=0.0, value_max=10.0)
  use_ncs = True
    .type = bool
    .short_caption = "Use NCS"
  output_label = result
    .type = key
  mode = fast *thorough exhaustive
    .type = choice
    .short_caption = "Refinement mode"
  methods = *fast slow exhaustive
    .type = choice(multi=True)
    .short_caption = "Active methods"
}

text {
  title = "Untitled"
    .type = str
    .short_caption = "Project title"
  tag = "default-tag"
    .type = qstr
  output_dir = None
    .type = path
    .short_caption = "Output directory"
  notes = "Edit me in the long-text tab for ample room."
    .type = str
    .short_caption = "Notes"
}

lists {
  labels = "F SIGF"
    .type = strings
    .short_caption = "MTZ labels"
  ncs_chains = "A B C"
    .type = words
    .short_caption = "NCS chains"
  cell_dims = 10 20 30
    .type = ints(size_min=3, size_max=3)
    .short_caption = "Cell dimensions"
  resolutions = 2.5 1.8 1.0
    .type = floats(value_min=0.5)
    .short_caption = "Resolution shells"
}

crystallography {
  space_group = "P 21 21 21"
    .type = space_group
    .short_caption = "Space group"
  unit_cell = 78.5 78.5 36.7 90 90 90
    .type = unit_cell
    .short_caption = "Unit cell"
  selection = "chain A and resseq 10:50"
    .type = atom_selection
    .short_caption = "Atom selection"
}

data_manager_demo {
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

aux_file = None
  .type = path

ncs_group {
  selection = "chain A"
    .type = str
    .short_caption = "Chain selection"
  weight = 1.0
    .type = float(value_min=0.0, value_max=10.0)
    .short_caption = "Weight"
}
"""


def build_master_phil():
  """Parse the example master scope and mark multi-instance objects."""
  master = iotbx.phil.parse(MASTER_PHIL)
  for obj in master.objects:
    if obj.name in ("aux_file", "ncs_group"):
      obj.multiple = True
  return master


def make_form(model, paths, widgets=None):
  """Build a QWidget hosting a QFormLayout of PhilFields for each path.

  Parameters
  ----------
    model: PhilModel
    paths: list of str
      Dotted PHIL paths to bind PhilFields to.
    widgets: dict of str -> type, optional
      Mapping from path to a non-default widget class (e.g., StrTextWidget).
  """
  widgets = widgets or {}
  panel = QWidget()
  layout = QFormLayout(panel)
  for path in paths:
    kwargs = {"widget": widgets[path]} if path in widgets else {}
    layout.addRow(PhilField(model, path, **kwargs))
  return panel


def make_data_manager_tab(model):
  """Build the 'Data Manager' tab.

  Parameters
  ----------
  model : qttbx.phil.PhilModel

  Returns
  -------
  QWidget
  """
  from qttbx.widgets.data_manager import (
    DataManagerWidget, DataManagerTableModel)
  from qttbx.examples._helpers import ProportionalColumns
  widget = DataManagerWidget(phil_model=model)
  # Match the 60/20/20 layout from data_manager_example.py. Keep a
  # reference on the widget so the helper isn't garbage-collected when
  # this function returns; QObject parenting alone would suffice but
  # an explicit attribute is easier to spot for a reader.
  widget._column_proportions = ProportionalColumns(
    widget._table,
    [
      (DataManagerTableModel.COL_FILENAME, 0.60),
      (DataManagerTableModel.COL_TYPE,     0.20),
      (DataManagerTableModel.COL_USED_FOR, 0.20),
    ])
  return widget


class ExampleWindow(QMainWindow):
  """Top-level showcase window.

  Splits a tree view of the entire PHIL scope from a tabbed form-view
  panel. Both views share the same PhilModel; the bottom-right tab shows
  the live extract for verification.
  """

  def __init__(self):
    QMainWindow.__init__(self)
    self.setWindowTitle("qttbx PHIL widgets -- example")
    self.resize(1200, 800)

    self._model = PhilModel()
    self._model.initialize_model(build_master_phil())

    central = QSplitter(Qt.Horizontal, self)
    self.setCentralWidget(central)

    # Left: tree view (every leaf editable via PhilItemDelegate).
    tree = QTreeView()
    tree.setModel(self._model)
    tree.setItemDelegate(PhilItemDelegate())
    tree.expandAll()
    tree.setColumnWidth(0, 230)
    tree.setAlternatingRowColors(True)
    tree.setEditTriggers(
      QTreeView.DoubleClicked | QTreeView.EditKeyPressed)
    central.addWidget(tree)

    # Right: tabbed forms.
    tabs = QTabWidget()
    central.addWidget(tabs)
    central.setStretchFactor(0, 1)
    central.setStretchFactor(1, 1)

    tabs.addTab(make_form(self._model, [
      "basic.count",
      "basic.weight",
      "basic.use_ncs",
      "basic.output_label",
      "basic.mode",
      "basic.methods",
    ]), "Scalars")

    tabs.addTab(make_form(self._model, [
      "text.title",
      "text.tag",
      "text.output_dir",
      "text.notes",
    ]), "Text")

    tabs.addTab(make_form(self._model, [
      "lists.labels",
      "lists.ncs_chains",
      "lists.cell_dims",
      "lists.resolutions",
    ]), "Lists")

    tabs.addTab(make_form(self._model, [
      "text.notes",
      "lists.labels",
      "lists.ncs_chains",
      "lists.cell_dims",
      "lists.resolutions",
      "crystallography.selection",
    ], widgets={
      "text.notes": StrTextWidget,
      "lists.labels": StringsTextWidget,
      "lists.ncs_chains": WordsTextWidget,
      "lists.cell_dims": IntsTextWidget,
      "lists.resolutions": FloatsTextWidget,
      "crystallography.selection": AtomSelectionTextWidget,
    }), "Long-text variants")

    tabs.addTab(make_form(self._model, [
      "crystallography.space_group",
      "crystallography.unit_cell",
      "crystallography.selection",
    ]), "Domain widgets")

    multi_panel = QWidget()
    mp_layout = QVBoxLayout(multi_panel)
    mp_layout.addWidget(QLabel(
      "<b>Definition multiple</b> -- aux_file (path .multiple = True), "
      "rendered with MultipleWidget:"))
    mp_layout.addWidget(PhilField(self._model, "aux_file"))
    mp_layout.addSpacing(16)
    mp_layout.addWidget(QLabel(
      "<b>Scope multiple</b> -- ncs_group { ... } .multiple = True, "
      "rendered with RepeatableScopeWidget:"))
    mp_layout.addWidget(RepeatableScopeWidget(self._model, "ncs_group"))
    mp_layout.addStretch(1)
    tabs.addTab(multi_panel, "Multiples")

    tabs.addTab(make_data_manager_tab(self._model), "Data Manager")

    self._extract_view = QTextEdit()
    self._extract_view.setReadOnly(True)
    self._extract_view.setStyleSheet("font-family: monospace;")
    refresh_btn = QPushButton("Refresh")
    refresh_btn.clicked.connect(self._refresh_extract)
    extract_panel = QWidget()
    ep_layout = QVBoxLayout(extract_panel)
    ep_layout.addWidget(QLabel(
      "Working PHIL serialized from the model. "
      "Click Refresh after editing to see updated values."))
    ep_layout.addWidget(refresh_btn)
    ep_layout.addWidget(self._extract_view, stretch=1)
    tabs.addTab(extract_panel, "Live extract")
    self._refresh_extract()

  def _refresh_extract(self):
    working = self._model.get_working_phil(diff_only=False)
    self._extract_view.setPlainText(working.as_str())


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
