# LIBTBX_SET_DISPATCHER_NAME phenix.development.fem

from __future__ import division
import mmtbx.f_model
import mmtbx.utils
import mmtbx.maps
import iotbx.phil
import iotbx.pdb
from cctbx import adptbx
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import runtime_utils
import os.path
import time
import sys
import random

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
random_seed=2679941
  .type = int
output {
  file_name = fem.mtz
    .type = path
    .style = hidden
  column_root_label = FEM
    .type = str
    .input_size = 120
    .short_caption = Base MTZ column label
  gui_output_dir = None
    .type = path
    .short_caption = Output directory
    .style = bold output_dir
}
include scope libtbx.phil.interface.tracking_params
""", process_includes=True)
master_params = master_phil # for phenix GUI

def manage_random_seed(random_seed):
  if(random_seed is None):
    random_seed = flex.get_random_seed()
  random.seed(random_seed)
  flex.set_random_seed(random_seed)

def run(args, command_name = "phenix.development.fem", log = sys.stdout):
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=log,
    process_pdb_file=False,
    force_non_anomalous=True,
    create_fmodel=False,
    usage_string="""
  phenix.development.fem model.pdb data.mtz
  phenix.development.fem data.mtz model.pdb f_obs_label=F

Calculate a "feature-enhanced" 2mFo-DFc map.
""")
  params = cmdline.params
  xray_structure = cmdline.xray_structure
  f_obs = cmdline.f_obs
  r_free_flags = cmdline.r_free_flags
  manage_random_seed(random_seed=params.random_seed)
  cs=f_obs.crystal_symmetry()
  crystal_gridding = f_obs.crystal_gridding(
    d_min                   = f_obs.d_min(),
    resolution_factor       = 0.25,
    grid_step               = None,
    symmetry_flags          = None,
    mandatory_factors       = None,
    max_prime               = 5,
    assert_shannon_sampling = True)
  fmodel = mmtbx.f_model.manager(
    f_obs = f_obs,
    r_free_flags = r_free_flags,
    xray_structure = xray_structure)
  fmodel.update_all_scales(update_f_part1_for=None)
  fmodel.show(show_approx=False)
  ### BEGIN: compute THE SIMPLEST possible 2mFo-DFc (the most original one)
  mc_orig = fmodel.electron_density_map(
    update_f_part1 = False).map_coefficients(
      map_type     = "2mFo-DFc",
      isotropize   = True,
      fill_missing = False)
  ### END: compute THE SIMPLEST possible 2mFo-DFc (the most original one)
  print >> log, "r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(), fmodel.r_free())
  ### b-factor sharpen
  xrs = fmodel.xray_structure
  b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
  print >> log, "Max B subtracted from atoms and used to sharpen map:", b_iso_min
  xrs.shift_us(b_shift=-b_iso_min)
  b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
  assert approx_equal(b_iso_min, 0, 1.e-3)
  fmodel.update_xray_structure(
    xray_structure = xrs,
    update_f_calc = True)
  #
  fmodel.update_all_scales(update_f_part1_for="refinement")
  #
  ko = mmtbx.maps.kick(
    fmodel           = fmodel.deep_copy(),
    crystal_gridding = crystal_gridding)
  #### Compute FEM start
  fem = mmtbx.maps.fem(
    ko=ko,
    crystal_gridding=crystal_gridding,
    fmodel=fmodel)
  #### Compute FEM end
  mtz_dataset = mc_orig.as_mtz_dataset(column_root_label="2mFoDFc")
  mtz_dataset.add_miller_array(
    miller_array=ko.complete_set,
    column_root_label="2mFoDFc_FilSharp")
  mtz_dataset.add_miller_array(
    miller_array=ko.map_coefficients,
    column_root_label="KICK")
  mtz_dataset.add_miller_array(
    miller_array=fem,
    column_root_label=params.output.column_root_label)
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = params.output.file_name)
  return os.path.abspath(params.output.file_name)

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=self.args, log=sys.stdout)

def validate_params (params) :
  return mmtbx.utils.validate_input_params(params)

if(__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time: %6.4f"%(time.time()-t0)
