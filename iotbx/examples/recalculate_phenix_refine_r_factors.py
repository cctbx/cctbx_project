"""
Read in an MTZ file produced by phenix.refine, extract the
F-obs-filtered, F-model, and R-free-flags arrays, and calculate R-factors both
for the entire dataset and for resolution shells.  This serves as an example
both for processing MTZ files, and for cctbx.miller functionality.
"""

from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_utils import get_r_free_flags_scores
from iotbx.file_reader import any_file
import sys

def compute_r_factors(fobs, fmodel, flags):
  fmodel, fobs = fmodel.common_sets(other=fobs)
  fmodel, flags = fmodel.common_sets(other=flags)
  fc_work = fmodel.select(~(flags.data()))
  fo_work = fobs.select(~(flags.data()))
  fc_test = fmodel.select(flags.data())
  fo_test = fobs.select(flags.data())
  r_work = fo_work.r1_factor(fc_work)
  r_free = fo_test.r1_factor(fc_test)
  print("r_work = %.4f" % r_work)
  print("r_free = %.4f" % r_free)
  print("")
  flags.setup_binner(n_bins=20)
  fo_work.use_binning_of(flags)
  fc_work.use_binner_of(fo_work)
  fo_test.use_binning_of(fo_work)
  fc_test.use_binning_of(fo_work)
  for i_bin in fo_work.binner().range_all():
    sel_work = fo_work.binner().selection(i_bin)
    sel_test = fo_test.binner().selection(i_bin)
    fo_work_bin = fo_work.select(sel_work)
    fc_work_bin = fc_work.select(sel_work)
    fo_test_bin = fo_test.select(sel_test)
    fc_test_bin = fc_test.select(sel_test)
    if fc_test_bin.size() == 0 : continue
    r_work_bin = fo_work_bin.r1_factor(other=fc_work_bin,
      assume_index_matching=True)
    r_free_bin = fo_test_bin.r1_factor(other=fc_test_bin,
      assume_index_matching=True)
    cc_work_bin = fo_work_bin.correlation(fc_work_bin).coefficient()
    cc_free_bin = fo_test_bin.correlation(fc_test_bin).coefficient()
    legend = flags.binner().bin_legend(i_bin, show_counts=False)
    print("%s  %8d %8d  %.4f %.4f  %.3f %.3f" % (legend, fo_work_bin.size(),
      fo_test_bin.size(), r_work_bin, r_free_bin, cc_work_bin, cc_free_bin))

def run(args):
  mtz_in = any_file(args[0])
  ma = mtz_in.file_server.miller_arrays
  flags = fmodel = fobs = None
  # select the output arrays from phenix.refine.  This could easily be modified
  # to handle MTZ files from other programs.
  for array in ma :
    labels = array.info().label_string()
    if labels.startswith("R-free-flags"):
      flags = array
    elif labels.startswith("F-model"):
      fmodel = abs(array)
    elif labels.startswith("F-obs-filtered"):
      fobs = array
  if (None in [flags, fobs, fmodel]):
    raise RuntimeError("Not a valid phenix.refine output file")
  scores = get_r_free_flags_scores([flags], None)
  test_flag_value = scores.test_flag_values[0]
  flags = flags.customized_copy(data=flags.data()==test_flag_value)
  compute_r_factors(fobs, fmodel, flags)

if (__name__ == "__main__"):
  run(sys.argv[1:])
