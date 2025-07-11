"""Calculate FEM maps"""
# LIBTBX_SET_DISPATCHER_NAME phenix.fem
# LIBTBX_SET_DISPATCHER_NAME phenix.feature_enhanced_map

from __future__ import absolute_import, division, print_function
import mmtbx.command_line
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

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="""
random_seed=2679941
  .type = int
use_omit = True
  .type = bool
  .help = Use composite OMIT protocol
sharp=True
  .type=bool
use_unsharp_masking = True
  .type = bool
resolution_factor = 1./4
  .type = float
signal_threshold = 0.5
  .type = float
use_resolve = Auto
  .type = bool
use_max_map = True
  .type = bool
ignore_zero_occupancy_atoms = True
  .type=bool
output {
  file_name_prefix = fem
    .type = str
    .input_size = 400
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
  if not r_free_flags:
    from libtbx.utils import Sorry
    raise Sorry("Please supply a file with r_free flags for FEM")
  manage_random_seed(random_seed=params.random_seed)
  cs=f_obs.crystal_symmetry()
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.ignore_zero_occupancy_atoms = params.ignore_zero_occupancy_atoms
  #
  fem = mmtbx.maps.fem.run(
    f_obs               = f_obs,
    r_free_flags        = r_free_flags,
    xray_structure      = xray_structure,
    use_resolve         = params.use_resolve,
    use_omit            = params.use_omit,
    use_max_map         = params.use_max_map,
    sharp               = params.sharp,
    use_unsharp_masking = params.use_unsharp_masking,
    resolution_factor   = params.resolution_factor,
    log                 = log)
  mtz_file_name = "%s.mtz"%params.output.file_name_prefix
  ccp4_map_file_name = "%s.ccp4"%params.output.file_name_prefix
  fem.write_output_files(
    mtz_file_name      = mtz_file_name,
    ccp4_map_file_name = ccp4_map_file_name,
    fem_label          = params.output.column_root_label,
    orig_label         = "2mFo-DFc")
  return os.path.abspath(mtz_file_name)

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=self.args, log=sys.stdout)

def validate_params(params):
  return mmtbx.command_line.validate_input_params(params)

if(__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print("Time: %6.4f"%(time.time()-t0))

