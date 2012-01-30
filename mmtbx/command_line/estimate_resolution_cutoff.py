
import libtbx.phil
from libtbx.utils import Usage
import sys

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
number_of_shells = 30
  .type = int
criteria = sigmaa cc *both
  .type = choice(multi=False)
min_cc = 0.3
  .type = float
min_sigmaa = 0.5
  .type = float
""", process_includes=True)

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  if (len(args) == 0) :
    raise Usage("""\
mmtbx.estimate_resolution_cutoff model.pdb data.mtz [options]

Attempts to determine a suitable resolution cutoff based on SigmaA and
correlation between F-obs and F-model.  This assumes that the reflections file
extends to a significantly higher resolution than was used in refinement of
the input model.
""")
  import mmtbx.utils
  from mmtbx.scaling.sigmaa_estimation import sigmaa_estimator
  import cctbx.miller
  from scitbx.array_family import flex
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False)
  f_obs = cmdline.fmodel.f_obs()
  f_model = cmdline.fmodel.f_model()
  r_free_flags = cmdline.fmodel.r_free_flags()
  params = cmdline.params
  sigmaa_obj = sigmaa_estimator(
    miller_obs=f_obs,
    miller_calc=f_model,
    r_free_flags=r_free_flags,
    kernel_width_free_reflections=100,
    n_sampling_points=params.number_of_shells)
  sigmaa_bins = sigmaa_obj.show_short(out=out)
  prev_d_min = f_obs.d_max_min()[0]
  sigmaa_cutoff = f_obs.d_min()
  for d_min, sigmaa in zip(sigmaa_bins.resolution, sigmaa_bins.sigmaa) :
    if (sigmaa < params.min_sigmaa) :
      sigmaa_cutoff = prev_d_min
      break
    prev_d_min = d_min
  print >> out, "Estimated cutoff based on SigmaA: %g" % sigmaa_cutoff
  f_obs.setup_binner(n_bins=params.number_of_shells)
  binner = f_obs.binner()
  cc_bins = []
  for i_bin in binner.range_all() :
    selection = binner.selection(i_bin)
    f_obs_bin = f_obs.select(selection)
    f_model_bin = f_model.select(selection)
    f_obs_data = f_obs_bin.data()
    f_model_data = flex.abs(f_model_bin.data())
    cc = flex.linear_correlation(f_obs_data, f_model_data).coefficient()
    cc_bins.append(cc)
  print >> out, "CC(F_obs, F_model) by resolution:"
  cc_binned = cctbx.miller.binned_data(
    binner=binner,
    data=cc_bins,
    data_fmt="%5.3f")
  cc_binned.show(show_unused=False, f=out, prefix="  ")
  cc_cutoff = f_obs.d_min()
  for i_bin in binner.range_used() :
    cc = cc_bins[i_bin]
    bin_range = binner.bin_d_range(i_bin)
    if (cc < params.min_cc) :
      cc_cutoff = bin_range[0]
      break
  print >> out, "Estimated cutoff based on CC: %g" % cc_cutoff

if (__name__ == "__main__") :
  run(sys.argv[1:])
