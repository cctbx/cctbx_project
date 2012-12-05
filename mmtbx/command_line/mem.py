from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.max_entropy_map

import sys
import mmtbx.utils
from libtbx.utils import user_plus_sys_time
import iotbx.phil
from iotbx import reflection_file_utils
from libtbx.utils import Sorry
import cctbx.maptbx.mem as mem

master_params_str="""\
hkl_file_name = None
  .type = str
label = None
  .type = str
pdb_file_name = None
  .type = str
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
solvent_fraction = None
  .type = float
f_000 = None
  .type = float
lam = 0.05
  .type = float
lambda_increment_factor = 1.02
  .type = float
output_file_name = None
  .type = str
mean_solvent_density = 0.35
  .type = float
max_iterations = 1000
  .type = float
resolution_factor = 0.3
  .type = float
beta = 0.7
  .type = float
convergence_at_r_factor = 0.05
  .type = float
"""

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def format_usage_message(log):
  print >> log, "-"*79
  msg = """\
phenix.max_entropy_map: map modification using Maximum Entropy Method (MEM)

Usage examples:
  phenix.max_entropy_map map_coeffs.mtz
  phenix.max_entropy_map model.pdb map_coeffs.mtz
  phenix.max_entropy_map model.pdb map.mtz label=['2FOFCWT', 'PH2FOFCWT']

Feedback:
  PAfonine@lbl.gov or phenixbb@phenix-online.org"""
  print >> log, msg
  print >> log, "-"*79

def run(args, log):
  timer = user_plus_sys_time()
  format_usage_message(log = log)
  if(len(args)==0): return
  parsed = iotbx.phil.parse(master_params_str)
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
    print >> log, "  file name:", inputs.pdb_file_names[0]
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
  map_coeffs = reflection_file_utils.extract_complex_miller_array_from_file(
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
  print >> log, "F(0,0,0): %12.6f"%f_000
  if(solvent_fraction is not None):
    print >> log, "solvent_fraction: %6.4f" % solvent_fraction
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
    log                     = log)
  ###
  broadcast(m="Output MEM map coefficients", log = log)
  ind = max(0,reff[0].rfind("."))
  ofn = reff[0]+"_mem.mtz" if ind==0 else reff[0][:ind]+"_mem.mtz"
  print >> log, "  Output file name:", ofn
  result.write_mtz_file(file_name = ofn)
  broadcast(m="All done", log=log)

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  run(sys.argv[1:], log=log)
  print >> log, "Total time: %-8.3f" % timer.elapsed()
