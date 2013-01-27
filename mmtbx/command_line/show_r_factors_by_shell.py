
from __future__ import division
from libtbx.utils import Sorry, Usage
import re
import sys

master_phil = """
mtz_file = None
  .type = path
f_obs_label = F-obs-filtered
  .type = str
f_calc_label = F-model
  .type = str
r_free_label = R-free-flags
  .type = str
n_bins = 10
  .type = int
"""

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""mmtbx.show_r_factors_by_shell DATA_FILE [OPTIONS]

example:
  mmtbx.show_r_factors_by_shell refine_001.mtz 10

Shows a table of R-factors by resolution shell, starting from the reflections
in a typical phenix.refine output MTZ file.

Full parameters:
%s
""" % (master_phil))
  from iotbx import file_reader
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil,
    reflection_file_def="mtz_file",
    integer_def="n_bins")
  params = cmdline.work.extract()
  if (params.mtz_file is None) :
    raise Sorry("No MTZ file supplied!")
  mtz_in = file_reader.any_file(params.mtz_file,
    force_type="hkl",
    raise_sorry_if_errors=True)
  f_calc = f_obs = r_free_flags = None
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().labels
    first_label_non_anom = re.sub("\(.*", "", labels[0])
    if ((labels[0] == params.f_obs_label) or
        (first_label_non_anom == params.f_obs_label)) :
      f_obs = array
    elif ((labels[0] == params.f_calc_label) or
          (first_label_non_anom == params.f_calc_label)) :
      f_calc = abs(array)
    elif ((labels[0] == params.r_free_label) or
          (first_label_non_anom == params.r_free_label)) :
      r_free_flags = array.customized_copy(data=array.data()==1)
  assert (not None in [f_obs, f_calc, r_free_flags])
  assert (f_obs.data().size() == f_calc.data().size())
  r_free_flags = r_free_flags.common_set(other=f_obs)
  f_obs.setup_binner(n_bins=params.n_bins)
  print >> out, ""
  for bin in f_obs.binner().range_used() :
    sele_bin = f_obs.binner().selection(bin)
    f_obs_bin = f_obs.select(sele_bin)
    f_calc_bin = f_calc.select(sele_bin)
    r_free_flags_bin = r_free_flags.select(sele_bin)
    f_obs_work = f_obs_bin.select(~(r_free_flags_bin.data()))
    f_calc_work = f_calc_bin.select(~(r_free_flags_bin.data()))
    f_obs_free = f_obs_bin.select(r_free_flags_bin.data())
    f_calc_free = f_calc_bin.select(r_free_flags_bin.data())
    r_work_bin = f_obs_work.r1_factor(f_calc_work)
    r_free_bin = f_obs_free.r1_factor(f_calc_free)
    print >> out, "  %6.3f - %6.3f  %6.4f  %6.4f" % (f_obs_work.d_max_min()[0],
      f_obs_work.d_min(), r_work_bin, r_free_bin)
  f_obs_work = f_obs.select(~(r_free_flags.data()))
  f_calc_work = f_calc.select(~(r_free_flags.data()))
  f_obs_free = f_obs.select(r_free_flags.data())
  f_calc_free = f_calc.select(r_free_flags.data())
  r_work = f_obs_work.r1_factor(f_calc_work)
  r_free = f_obs_free.r1_factor(f_calc_free)
  print >> out, ""
  print >> out, "  %6.3f - %6.3f  %6.4f  %6.4f" % (f_obs.d_max_min()[0],
      f_obs.d_min(), r_work, r_free)
  print >> out, ""

if (__name__ == "__main__") :
  run(sys.argv[1:])
