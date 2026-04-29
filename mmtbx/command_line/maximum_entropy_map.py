"""Compute maximum entropy map from map coefficients"""
# LIBTBX_SET_DISPATCHER_NAME phenix.maximum_entropy_map

from __future__ import absolute_import, division, print_function
import mmtbx.utils
import iotbx.phil
from iotbx import reflection_file_utils
import cctbx.maptbx.mem as mem
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx import runtime_utils
import os.path
import sys

master_params_str="""\
hkl_file_name = None
  .type = path
  .short_caption = Map coefficients
  .style = bold file_type:mtz process_hkl child:map_arrays:label input_file
label = None
  .type = str
  .input_size = 160
  .short_caption = Column labels
  .style = bold renderer:draw_map_arrays_widget noauto
pdb_file_name = None
  .type = path
  .short_caption = Model file
  .style = file_type:pdb
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
solvent_fraction = None
  .type = float
f_000 = None
  .type = float
  .short_caption = F(0,0,0)
lam = 0.05
  .type = float
  .short_caption = Lambda
lambda_increment_factor = 1.02
  .type = float
output_file_name = None
  .type = path
  .short_caption = Output file
  .style = bold output_file file_type:mtz
column_root_label = MEM
  .type = str
  .input_size = 100
  .short_caption = Output label base
output_high_resolution = None
  .type = float
mean_solvent_density = 0.35
  .type = float
max_iterations = 1000
  .type = float
resolution_factor = 0.3
  .type = float
  .short_caption = Grid resolution factor
beta = 0.7
  .type = float
convergence_at_r_factor = 0.05
  .type = float
convergence_r_threshold = 0.1
  .type = float
gui_output_dir = None
  .type = path
  .short_caption = Output directory
  .style = output_dir
include scope libtbx.phil.interface.tracking_params
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

def format_usage_message(log):
  print("-"*79, file=log)
  msg = """\
phenix.max_entropy_map: map modification using Maximum Entropy Method (MEM)

Usage examples:
  phenix.max_entropy_map map_coeffs.mtz
  phenix.max_entropy_map model.pdb map_coeffs.mtz
  phenix.max_entropy_map model.pdb map.mtz label="2FOFCWT,PH2FOFCWT"

Feedback:
  PAfonine@lbl.gov or phenixbb@phenix-online.org"""
  print(msg, file=log)
  print("-"*79, file=log)

def run(args, log):
  timer = user_plus_sys_time()
  format_usage_message(log = log)
  if(len(args)==0): return
  parsed = master_params()
  inputs = mmtbx.utils.process_command_line_args(
    args=args, master_params=parsed, log=log)
  params = inputs.params.extract()
  broadcast(m="Input parameters", log = log)
  inputs.params.show(prefix="  ")
  ###
  xray_structure = None
  if(len(inputs.pdb_file_names)>0):
    broadcast(m="Input model", log = log)
    assert len(inputs.pdb_file_names) == 1
    print("  file name:", inputs.pdb_file_names[0], file=log)
    xray_structure = iotbx.pdb.input(
      file_name = inputs.pdb_file_names[0]).xray_structure_simple()
    assert xray_structure is not None
    xray_structure.show_summary(prefix="  ", f=log)
    mmtbx.utils.setup_scattering_dictionaries(
      scattering_table = params.scattering_table,
      xray_structure   = xray_structure,
      d_min            = 0.25)
    xray_structure.scattering_type_registry().show(prefix="  ", out = log)
  ###
  broadcast(m="Input reflection data", log = log)
  reff = inputs.reflection_file_names
  if(len(reff) > 1):
    raise Sorry("One reflection file should be provided.")
  elif(len(reff) == 0):
    if(params.hkl_file_name is None):
      raise Sorry("No reflection file provided.")
    else: reff = [params.hkl_file_name]
  map_coeffs = reflection_file_utils.extract_miller_array_from_file(
    file_name = reff[0],
    label     = params.label,
    type      = "complex",
    log       = log)
  assert map_coeffs is not None
  map_coeffs.show_comprehensive_summary(prefix="  ", f=log)
  ###
  broadcast(m="MEM calculations begin", log = log)
  f_000 = params.f_000
  solvent_fraction = params.solvent_fraction
  if(f_000 is None):
    f_000_obj = mmtbx.utils.f_000(
      xray_structure       = xray_structure,
      unit_cell_volume     = map_coeffs.unit_cell().volume(),
      solvent_fraction     = params.solvent_fraction,
      mean_solvent_density = params.mean_solvent_density)
    f_000 = f_000_obj.f_000
    solvent_fraction = f_000_obj.solvent_fraction
  print("F(0,0,0): %12.6f"%f_000, file=log)
  if(solvent_fraction is not None):
    print("solvent_fraction: %6.4f" % solvent_fraction, file=log)
  result = mem.run(
    f                       = map_coeffs,
    f_000                   = f_000,
    lam                     = params.lam,
    lambda_increment_factor = params.lambda_increment_factor,
    resolution_factor       = params.resolution_factor,
    verbose                 = True,
    start_map               = "min_shifted",
    max_iterations          = params.max_iterations,
    use_modification        = True,
    beta                    = params.beta,
    convergence_at_r_factor = params.convergence_at_r_factor,
    xray_structure          = xray_structure,
    convergence_r_threshold = params.convergence_r_threshold,
    log                     = log)
  ###
  broadcast(m="Output MEM map coefficients", log = log)
  ind = max(0,reff[0].rfind("."))
  ofn = params.output_file_name
  if (ofn is None):
    ofn = reff[0]+"_mem.mtz" if ind==0 else reff[0][:ind]+"_mem.mtz"
  print("  Output file name:", ofn, file=log)
  result.write_mtz_file(file_name = ofn,
    column_root_label=params.column_root_label,
    d_min=params.output_high_resolution)
  broadcast(m="All done", log=log)
  return os.path.abspath(ofn)

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args), log=sys.stdout)

def validate_params(params):
  if (params.hkl_file_name is None):
    raise Sorry("Please specify an MTZ file containing map coefficients.")
  elif (params.label is None):
    raise Sorry("No column labels specified.")
  if (params.resolution_factor > 0.5):
    raise Sorry("The grid resolution factor must be a decimal number <= 0.5.")
  return True

def finish_job(result):
  output_files, stats = [], []
  if (result is not None) and (os.path.isfile(result)):
    output_files.append((result, "Maximum entropy map"))
  return output_files, stats

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  run(sys.argv[1:], log=log)
  print("Total time: %-8.3f" % timer.elapsed(), file=log)

