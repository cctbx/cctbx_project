
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage
from libtbx import slots_getstate_setstate
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

class r_factor_shell(object):
  __slots__ = ["d_min", "d_max", "r_work", "r_free"]

  def __init__(self, d_min, d_max, r_work, r_free):
    self.d_min = d_min
    self.d_max = d_max
    self.r_work = r_work
    self.r_free = r_free

  def show(self, out, prefix="  "):
    print("%s%6.3f - %6.3f  %6.4f  %6.4f" % (prefix, self.d_max,
      self.d_min, self.r_work, self.r_free), file=out)

class r_factor_shells(slots_getstate_setstate):
  __slots__ = ["shells", "overall"]

  def __init__(self, f_obs, f_calc, r_free_flags, n_bins=10):
    assert (not None in [f_obs, f_calc, r_free_flags])
    assert (f_obs.data().size() == f_calc.data().size())
    r_free_flags = r_free_flags.common_set(other=f_obs)
    f_obs.setup_binner(n_bins=n_bins)
    self.shells = []
    self.overall = None
    for bin in f_obs.binner().range_used():
      sele_bin = f_obs.binner().selection(bin)
      f_obs_bin = f_obs.select(sele_bin)
      f_calc_bin = f_calc.select(sele_bin)
      r_free_flags_bin = r_free_flags.select(sele_bin)
      f_obs_work = f_obs_bin.select(~(r_free_flags_bin.data()))
      f_calc_work = f_calc_bin.select(~(r_free_flags_bin.data()))
      f_obs_free = f_obs_bin.select(r_free_flags_bin.data())
      f_calc_free = f_calc_bin.select(r_free_flags_bin.data())
      shell = r_factor_shell(
        d_min=f_obs_work.d_min(),
        d_max=f_obs_work.d_max_min()[0],
        r_work=f_obs_work.r1_factor(f_calc_work),
        r_free=f_obs_free.r1_factor(f_calc_free))
      self.shells.append(shell)
    f_obs_work = f_obs.select(~(r_free_flags.data()))
    f_calc_work = f_calc.select(~(r_free_flags.data()))
    f_obs_free = f_obs.select(r_free_flags.data())
    f_calc_free = f_calc.select(r_free_flags.data())
    self.overall = r_factor_shell(
      d_min=f_obs.d_min(),
      d_max=f_obs.d_max_min()[0],
      r_work = f_obs_work.r1_factor(f_calc_work),
      r_free = f_obs_free.r1_factor(f_calc_free))

  def show(self, out=sys.stdout, prefix="  "):
    print("", file=out)
    for shell in self.shells :
      shell.show(out=out, prefix=prefix)
    print("", file=out)
    self.overall.show(out=out, prefix=prefix)
    print("", file=out)

def run(args, out=sys.stdout):
  if (len(args) == 0) or ("--help" in args):
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
  if (params.mtz_file is None):
    raise Sorry("No MTZ file supplied!")
  mtz_in = file_reader.any_file(params.mtz_file,
    force_type="hkl",
    raise_sorry_if_errors=True)
  f_calc = f_obs = r_free_flags = None
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().labels
    first_label_non_anom = re.sub(r"\(.*", "", labels[0])
    if ((labels[0] == params.f_obs_label) or
        (first_label_non_anom == params.f_obs_label)):
      f_obs = array
    elif ((labels[0] == params.f_calc_label) or
          (first_label_non_anom == params.f_calc_label)):
      f_calc = abs(array)
    elif ((labels[0] == params.r_free_label) or
          (first_label_non_anom == params.r_free_label)):
      r_free_flags = array.customized_copy(data=array.data()==1)
  shells = r_factor_shells(
    f_obs=f_obs,
    f_calc=f_calc,
    r_free_flags=r_free_flags,
    n_bins=params.n_bins)
  shells.show(out=out)
  return shells

if (__name__ == "__main__"):
  run(sys.argv[1:])
