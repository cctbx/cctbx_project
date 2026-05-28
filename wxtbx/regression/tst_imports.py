from __future__ import absolute_import, division, print_function
import importlib

from wxtbx.regression._warnings import install_wx_deprecation_filters

# Install once at module load so that even our own importlib.import_module
# calls below are covered before the loop starts.  Re-installed per-iteration
# inside exercise() to defeat transitive filter mutations.
install_wx_deprecation_filters()

# Every in-scope wxtbx submodule, listed explicitly so the test is the
# authoritative inventory.  The two mmtbx files are also included.
MODULES = [
  "wxtbx",
  "wxtbx.app",
  "wxtbx.bitmaps",
  "wxtbx.icons",
  "wxtbx.utils",
  "wxtbx.layout",
  "wxtbx.errors",
  "wxtbx.windows",
  "wxtbx.metallicbutton",
  "wxtbx.segmentedctrl",
  "wxtbx.checklistctrl",
  "wxtbx.path_dialogs",
  "wxtbx.misc_dialogs",
  "wxtbx.info_panels",
  "wxtbx.listeditor",
  "wxtbx.tables",
  "wxtbx.symmetry_dialog",
  "wxtbx.browser",
  "wxtbx.heavy_atoms",
  "wxtbx.adp_statistics",
  "wxtbx.anomalous_scattering",
  "wxtbx.b_plot",
  "wxtbx.polygon",
  "wxtbx.polygon_db_viewer",
  "wxtbx.xtriage",
  "wxtbx.pdb_editor",
  "wxtbx.sequence_view",
  "wxtbx.ensemble_validation",
  "wxtbx.custom_restraints",
  "wxtbx.mtz_dump",
  "wxtbx.process_control",
  "wxtbx.monitor",
  "wxtbx.qstat_view",
  "wxtbx.phil_controls",
  "wxtbx.phil_controls.intctrl",
  "wxtbx.phil_controls.floatctrl",
  "wxtbx.phil_controls.strctrl",
  "wxtbx.phil_controls.boolctrl",
  "wxtbx.phil_controls.choice",
  "wxtbx.phil_controls.choice_multi",
  "wxtbx.phil_controls.text_base",
  "wxtbx.phil_controls.numbers",
  "wxtbx.phil_controls.ints",
  "wxtbx.phil_controls.floats",
  "wxtbx.phil_controls.path",
  "wxtbx.phil_controls.tree",
  "wxtbx.phil_controls.simple_dialogs",
  "wxtbx.phil_controls.misc_ctrls",
  "wxtbx.phil_controls.space_group",
  "wxtbx.phil_controls.symop",
  "wxtbx.phil_controls.unit_cell",
  "wxtbx.phil_controls.array_label",
  "wxtbx.phil_controls.tribool",
  "wxtbx.phil_controls.strings",
  "wxtbx.plots",
  "wxtbx.plots.molprobity",
  "wxtbx.plots.emringer",
  "wxtbx.reports.blast",
  "wxtbx.reports.pdb_symmetry",
  "gltbx.wx_viewer",
  "gltbx.wx_viewer_leapmotion",
  "mmtbx.command_line.map_box",
]

def exercise():
  failed = []
  for name in MODULES:
    # Re-prepend our filters: a previous import may have transitively pulled
    # in phaser.deprecated_keywords_bpl, which prepends a broad "default"
    # filter that would mask ours.
    install_wx_deprecation_filters()
    try:
      importlib.import_module(name)
      print("OK   {}".format(name))
    except Exception as e:
      failed.append((name, e))
      print("FAIL {}: {}: {}".format(name, type(e).__name__, e))
  if failed:
    print("{} module(s) failed to import cleanly".format(len(failed)))
    raise SystemExit(1)
  print("OK")

if __name__ == "__main__":
  exercise()
