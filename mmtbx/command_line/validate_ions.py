
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header
from libtbx import Auto
import sys

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    enable_pdb_interpretation_params=True,
    enable_stop_for_unknowns=False,
    phil_string="""
include scope mmtbx.ions.identify.ion_master_phil
debug = True
  .type = bool
nproc = Auto
  .type = int
""")

def run(args, out=sys.stdout):
  usage_string="""\
mmtbx.validate_ions model.pdb data.mtz [options ...]

Utility to validate ions that have been built into a model, based on local
environment, electron density maps, and atomic properties.
"""
  import mmtbx.ions.identify
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=True,
    create_fmodel=True,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    prefer_anomalous=True)
  fmodel = cmdline.fmodel
  xray_structure = cmdline.xray_structure
  params = cmdline.params
  pdb_hierarchy = cmdline.pdb_hierarchy
  geometry = cmdline.geometry
  make_header("Inspecting ions", out=out)
  manager = mmtbx.ions.identify.create_manager(
    pdb_hierarchy = pdb_hierarchy,
    fmodel=fmodel,
    geometry_restraints_manager=geometry,
    wavelength=params.input.wavelength,
    params=params,
    verbose = params.debug,
    nproc = params.nproc,
    log=out)
  manager.show_current_scattering_statistics(out=out)
  return manager.validate_ions(out = out, debug = params.debug)

if (__name__ == "__main__"):
  run(sys.argv[1:])
