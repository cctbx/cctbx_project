from __future__ import absolute_import, division, print_function

# Case rows: (label, factory, requires_data, skip_reason)
#   label: human-readable id, "<module>.<ClassName>"
#   factory: callable(parent) -> widget; pass None for requires_data=True
#   requires_data: bool; True means the class needs data we don't synthesize
#   skip_reason: str; only consulted when requires_data is True
#
# Add cases incrementally as the per-module phase progresses.  The harness
# treats this list as authoritative; the spec's tier lists are representative.

CASES = []


def _make_metallicbutton(parent):
  from wxtbx.metallicbutton import MetallicButton
  return MetallicButton(parent, -1, label="x")

CASES.append(("wxtbx.metallicbutton.MetallicButton", _make_metallicbutton, False, ""))


def _add_phil_controls_cases():
  from wxtbx.phil_controls.intctrl import IntCtrl
  from wxtbx.phil_controls.floatctrl import FloatCtrl
  from wxtbx.phil_controls.strctrl import StrCtrl
  from wxtbx.phil_controls.boolctrl import BoolCtrl
  from wxtbx.phil_controls.choice import ChoiceCtrl
  from wxtbx.phil_controls.path import PathCtrl
  from wxtbx.phil_controls.tribool import TriBoolCtrl
  from wxtbx.phil_controls.unit_cell import UnitCellCtrl
  from wxtbx.phil_controls.space_group import SpaceGroupCtrl
  from wxtbx.phil_controls.symop import SymopCtrl
  from wxtbx.phil_controls.array_label import ArrayLabelCtrl
  from wxtbx.phil_controls.tree import PhilTreeCtrl
  CASES.extend([
    ("wxtbx.phil_controls.intctrl.IntCtrl",        lambda p: IntCtrl(p),        False, ""),
    ("wxtbx.phil_controls.floatctrl.FloatCtrl",    lambda p: FloatCtrl(p),      False, ""),
    ("wxtbx.phil_controls.strctrl.StrCtrl",        lambda p: StrCtrl(p),        False, ""),
    ("wxtbx.phil_controls.boolctrl.BoolCtrl",      lambda p: BoolCtrl(p),       False, ""),
    ("wxtbx.phil_controls.choice.ChoiceCtrl",      lambda p: ChoiceCtrl(p, choices=["a", "b"]), False, ""),
    ("wxtbx.phil_controls.path.PathCtrl",          lambda p: PathCtrl(p),       False, ""),
    ("wxtbx.phil_controls.tribool.TriBoolCtrl",    lambda p: TriBoolCtrl(p),    False, ""),
    ("wxtbx.phil_controls.unit_cell.UnitCellCtrl", lambda p: UnitCellCtrl(p),   False, ""),
    ("wxtbx.phil_controls.space_group.SpaceGroupCtrl", lambda p: SpaceGroupCtrl(p), False, ""),
    ("wxtbx.phil_controls.symop.SymopCtrl",        lambda p: SymopCtrl(p),      False, ""),
    ("wxtbx.phil_controls.array_label.ArrayLabelCtrl", lambda p: ArrayLabelCtrl(p), False, ""),
    ("wxtbx.phil_controls.tree.PhilTreeCtrl",      lambda p: PhilTreeCtrl(p),   False, ""),
  ])


def _add_segmented_and_checklist_cases():
  from wxtbx.segmentedctrl import (
    SegmentedButtonControl, SegmentedRadioControl, SegmentedToggleControl,
  )
  from wxtbx.checklistctrl import CheckListCtrl
  CASES.extend([
    ("wxtbx.segmentedctrl.SegmentedButtonControl", lambda p: SegmentedButtonControl(p), False, ""),
    ("wxtbx.segmentedctrl.SegmentedRadioControl",  lambda p: SegmentedRadioControl(p),  False, ""),
    ("wxtbx.segmentedctrl.SegmentedToggleControl", lambda p: SegmentedToggleControl(p), False, ""),
    ("wxtbx.checklistctrl.CheckListCtrl",          lambda p: CheckListCtrl(p),          False, ""),
  ])


def _add_tables_and_listeditor_cases():
  from wxtbx.tables import TableView
  from wxtbx.listeditor import ListEditor
  CASES.extend([
    ("wxtbx.tables.TableView",      lambda p: TableView(p),  False, ""),
    ("wxtbx.listeditor.ListEditor", lambda p: ListEditor(p), False, ""),
  ])


def _add_dialog_cases():
  from wxtbx.misc_dialogs import ChoiceDialog
  from wxtbx.path_dialogs import ConfirmCleanupDialog, DirectoryCleanupDialog
  from wxtbx.symmetry_dialog import SymmetryChoiceDialog, SymmetryDialog
  from wxtbx.heavy_atoms import OccupancyDialog
  from wxtbx.pdb_editor import AddResiduesDialog, SelectChainDialog

  class _CleanupStub(object):
    """Minimal stand-in for libtbx.path.clean_out_directory.

    ConfirmCleanupDialog reads n_files, n_dirs, get_freed_space(), dir_paths
    and file_paths; supply empty values so the dialog renders without
    walking the filesystem.
    """
    n_files = 0
    n_dirs = 0
    dir_paths = ()
    file_paths = ()
    def get_freed_space(self):
      return "0.0 KB"

  def _one_chain():
    import iotbx.pdb
    pdb_in = iotbx.pdb.input(source_info=None, lines=(
      "ATOM      1  N   ALA A   1      11.000  12.000  13.000  "
      "1.00 20.00           N\nEND\n"))
    return list(pdb_in.construct_hierarchy().chains())

  CASES.extend([
    ("wxtbx.misc_dialogs.ChoiceDialog",
        lambda p: ChoiceDialog(title="t", message="m",
                               choices=["a", "b"], parent=p),
        False, ""),
    ("wxtbx.path_dialogs.ConfirmCleanupDialog",
        lambda p: ConfirmCleanupDialog(p, cleanup_obj=_CleanupStub()),
        False, ""),
    ("wxtbx.path_dialogs.DirectoryCleanupDialog",
        lambda p: DirectoryCleanupDialog(p),
        False, ""),
    ("wxtbx.symmetry_dialog.SymmetryChoiceDialog",
        lambda p: SymmetryChoiceDialog(p, title="title"),
        False, ""),
    ("wxtbx.symmetry_dialog.SymmetryDialog",
        lambda p: SymmetryDialog(p),
        False, ""),
    ("wxtbx.heavy_atoms.OccupancyDialog",
        lambda p: OccupancyDialog(p),
        False, ""),
    ("wxtbx.pdb_editor.AddResiduesDialog",
        lambda p: AddResiduesDialog(p),
        False, ""),
    ("wxtbx.pdb_editor.SelectChainDialog",
        lambda p: SelectChainDialog(p, title="t", message="m",
                                    chains=_one_chain()),
        False, ""),
  ])


def _add_plot_cases():
  from wxtbx.plots import (
    plot_container, histogram, image_plot, iotbx_data_plot_base,
    small_plot, plot_frame, loggraph,
  )
  from wxtbx.b_plot import b_plot_panel
  from wxtbx.polygon import ColorBox

  def _stub_table():
    """Build a minimal iotbx.data_plots.table_data with one graph."""
    from iotbx import data_plots
    return data_plots.table_data(
      title="stub",
      column_labels=["x", "y"],
      graph_names=["plot1"],
      graph_columns=[[0, 1]],
      data=[[0.0, 1.0, 2.0], [0.5, 1.5, 2.5]],
    )

  CASES.extend([
    ("wxtbx.plots.plot_container",
        lambda p: plot_container(parent=p),
        False, ""),
    ("wxtbx.plots.histogram",
        lambda p: histogram(parent=p),
        False, ""),
    ("wxtbx.plots.image_plot",
        lambda p: image_plot(parent=p),
        False, ""),
    ("wxtbx.plots.iotbx_data_plot_base",
        lambda p: iotbx_data_plot_base(parent=p, tables=[_stub_table()]),
        False, ""),
    ("wxtbx.plots.small_plot",
        lambda p: small_plot(parent=p, table=_stub_table()),
        False, ""),
    ("wxtbx.plots.plot_frame",
        lambda p: plot_frame(parent=p, title="t"),
        False, ""),
    ("wxtbx.plots.loggraph",
        lambda p: loggraph(parent=p, title="t", tables=[_stub_table()]),
        False, ""),
    ("wxtbx.b_plot.b_plot_panel",
        lambda p: b_plot_panel(parent=p),
        False, ""),
    ("wxtbx.polygon.ColorBox",
        lambda p: ColorBox(parent=p),
        False, ""),
  ])

  # PolygonPanel requires an mmtbx.polygon.output.renderer (which itself
  # needs histogram_data + structure_stats with strict shape constraints).
  # Building a real one is the panel's responsibility, not the harness's.
  CASES.append((
    "wxtbx.polygon.PolygonPanel", None, True,
    "needs mmtbx.polygon.output.renderer with histogram_data and stats"))


def _add_tier4_skips():
  CASES.extend([
    ("wxtbx.xtriage.XtriageFrame",                          None, True, "needs xtriage analysis result"),
    ("wxtbx.pdb_editor.PDBTreeFrame",                       None, True, "needs a model"),
    ("wxtbx.adp_statistics.ADPStatisticsFrame",             None, True, "needs model + ADP statistics"),
    ("wxtbx.anomalous_scattering.AnomPlotFrame",            None, True, "needs miller arrays"),
    ("wxtbx.sequence_view.msa_frame",                       None, True, "needs sequences and structure"),
    ("wxtbx.sequence_view.sequence_frame",                  None, True, "needs sequences"),
    ("wxtbx.ensemble_validation.ensemble_validation_panel", None, True, "needs ensemble results"),
    ("wxtbx.custom_restraints.CustomRestraintsPanel",       None, True, "needs restraints data"),
    ("wxtbx.info_panels.PDBFileInfo",                       None, True, "needs a file"),
    ("wxtbx.info_panels.MapCoeffsInfo",                     None, True, "needs a file"),
    ("wxtbx.info_panels.ReflectionFileInfo",                None, True, "needs a file"),
    ("wxtbx.info_panels.ImageFileInfo",                     None, True, "needs a file"),
    ("wxtbx.polygon.PolygonFrame",                          None, True, "needs structure stats"),
    ("wxtbx.plots.molprobity.multi_criterion_frame",        None, True, "needs MolProbity results"),
    ("wxtbx.plots.molprobity.ramalyze_frame",               None, True, "needs Ramachandran data"),
    ("wxtbx.plots.molprobity.rotalyze_frame",               None, True, "needs rotamer data"),
    ("wxtbx.plots.emringer.emringer_plot_frame",            None, True, "needs EMRinger data"),
    ("wxtbx.polygon_db_viewer",                             None, True, "frame needs database data"),
    ("wxtbx.b_plot.BPlotFrame",                             None, True, "needs B-factor histogram data"),
  ])


_add_phil_controls_cases()
_add_segmented_and_checklist_cases()
_add_tables_and_listeditor_cases()
_add_dialog_cases()
_add_plot_cases()
_add_tier4_skips()
