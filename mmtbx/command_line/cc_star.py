
from __future__ import division
from libtbx.utils import Sorry, Usage
import libtbx.phil
import sys

def merging_and_model_statistics (fmodel, unmerged_i_obs, n_bins=20) :
  from iotbx import merging_statistics
  # very important: must use original intensities for i_obs, not squared f_obs
  # from fmodel, because French-Wilson treatment is one-way
  i_obs = unmerged_i_obs.merge_equivalents(use_internal_variance=False).array()
  if (i_obs.anomalous_flag()) :
    i_obs = i_obs.average_bijvoet_mates()
  f_model = fmodel.f_model()
  work_sel = f_model.customized_copy(data=fmodel.arrays.work_sel)
  free_sel = f_model.customized_copy(data=fmodel.arrays.free_sel)
  if (f_model.anomalous_flag()) :
    f_model = f_model.average_bijvoet_mates()
    work_sel = work_sel.average_bijvoet_mates()
    free_sel = free_sel.average_bijvoet_mates()
  i_obs, work_sel = i_obs.common_sets(other=work_sel)
  i_obs, free_sel = i_obs.common_sets(other=free_sel)
  i_obs, f_model = i_obs.common_sets(other=f_model)
  i_calc = abs(f_model).f_as_f_sq()
  d_max, d_min = i_calc.d_max_min()
  model_arrays = merging_statistics.model_based_arrays(
    i_obs=i_obs,
    i_calc=i_calc,
    work_sel=work_sel,
    free_sel=free_sel)
  return merging_statistics.dataset_statistics(
    i_obs=unmerged_i_obs,
    crystal_symmetry=i_calc,
    d_min=d_min,
    d_max=d_max,
    n_bins=n_bins,
    model_arrays=model_arrays)

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
unmerged_data = None
  .type = path
unmerged_labels = None
  .type = str
n_bins = 20
  .type = int(value_min=5)
""", process_includes=True)

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""\
mmtbx.cc_star model.pdb data.mtz unmerged_data=data.hkl [n_bins=X] [options]

Implementation of method described in:
Karplus PA & Diederichs K (2012) Science 336:1030-3.

Full parameters:
%s
""" % master_phil.as_str(prefix=" "))
  import mmtbx.utils
  from iotbx import merging_statistics
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    scattering_table="n_gaussian",
    out=out)
  params = cmdline.params
  fmodel = cmdline.fmodel
  if (params.unmerged_data is None) :
    raise Sorry("Please specify unmerged_data file")
  unmerged_i_obs = merging_statistics.select_data(
    file_name=params.unmerged_data,
    data_labels=params.unmerged_labels,
    log=out)
  stats = merging_and_model_statistics(
    fmodel=fmodel,
    unmerged_i_obs=unmerged_i_obs,
    n_bins=params.n_bins)
  stats.show_cc_star(out=out)
  print >> out, ""
  print >> out, "Reference:"
  print >> out, "  Karplus PA & Diederichs K (2012) Science 336:1030-3."
  print >> out, ""
  return stats

if (__name__ == "__main__") :
  run(sys.argv[1:])
