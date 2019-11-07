from __future__ import absolute_import, division, print_function
import iotbx.pdb
import sys
import mmtbx.utils
import mmtbx.building.loop_closure.utils
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran


def run(args):
  print(args)
  print(args[0])
  print(args[1])
  print(args[2])
  pdb_inp = iotbx.pdb.input(source_info=None, file_name=args[0])
  pdb_h = pdb_inp.construct_hierarchy()
  ref_h = iotbx.pdb.input(source_info=None, file_name=args[1]).construct_hierarchy()
  # print dir(pdb_inp)
  xrs = pdb_h.extract_xray_structure()
  outlier_selection_txt = mmtbx.building.loop_closure.utils. \
    rama_score_selection(ref_h, None, "outlier", 1)
  negate_selection = "all"
  # asc = ref_h.atom_selection_cache()
  if outlier_selection_txt != "" and outlier_selection_txt is not None:
    negate_selection = "not (%s)" % outlier_selection_txt
  # sel = asc.selection(negate_selection)

  minimize_wrapper_for_ramachandran(
      pdb_h,
      xrs,
      ref_h,
      excl_string_selection=negate_selection,
      log=None,
      run_first_minimization_without_reference=True,
      oldfield_weight_scale=1,
      oldfield_plot_cutoff=0.027,
      nonbonded_weight=500,
      reference_sigma=0.5)
  pdb_h.write_pdb_file(file_name=args[2])


if (__name__ == "__main__"):
  print("__name__", __name__)
  run(sys.argv[1:])
