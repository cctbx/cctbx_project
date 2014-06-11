# LIBTBX_SET_DISPATCHER_NAME phenix.fem
# LIBTBX_SET_DISPATCHER_NAME phenix.feature_enhanced_map

from __future__ import division
import mmtbx.command_line
import mmtbx.f_model
import mmtbx.maps
import iotbx.phil
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import runtime_utils
import os.path
import time
import sys
import random
import mmtbx.maps.fem
import mmtbx.masks

def get_master_phil () :
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="""
random_seed=2679941
  .type = int
omit = False
  .type = bool
  .help = Use composite OMIT protocol
sharp=True
  .type=bool
signal_threshold = 0.5
  .type = float
ignore_zero_occupancy_atoms = True
  .type=bool
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
""",
    enable_twin_law=False)

master_phil = get_master_phil()
master_params = master_phil # for phenix GUI

def manage_random_seed(random_seed):
  if(random_seed is None):
    random_seed = flex.get_random_seed()
  random.seed(random_seed)
  flex.set_random_seed(random_seed)

def run(args, command_name = "phenix.development.fem", log = sys.stdout):
  cmdline = mmtbx.command_line.load_model_and_data(
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
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.ignore_zero_occupancy_atoms = params.ignore_zero_occupancy_atoms
  #
  sel = f_obs.data()>0
  f_obs = f_obs.select(sel)
  r_free_flags = r_free_flags.select(sel)
  #
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xray_structure,
    mask_params    = mask_params)
  fmodel.update_all_scales(update_f_part1_for=None)
  fmodel.show(show_approx=False)
  print >> log, "r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(), fmodel.r_free())
  fem = mmtbx.maps.fem.run(fmodel=fmodel, use_omit=params.omit,
    sharp=params.sharp, signal_threshold=params.signal_threshold)
  # output
  mtz_dataset = fem.mc.as_mtz_dataset(column_root_label="2mFoDFc")
  if(fem.mc_resolve is not None):
    mtz_dataset.add_miller_array(
      miller_array=fem.mc_resolve,
      column_root_label="ResolveDM")
  mtz_dataset.add_miller_array(
    miller_array=fem.mc_result,
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
  return mmtbx.command_line.validate_input_params(params)

if(__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time: %6.4f"%(time.time()-t0)
