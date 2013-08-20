
from __future__ import division
from libtbx.str_utils import make_header
from libtbx.utils import Usage
from libtbx import Auto
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
include scope mmtbx.ions.ion_master_phil
debug = True
  .type = bool
wavelength = None
  .type = float
nproc = Auto
  .type = int
pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
  stop_for_unknowns = False
    .type = bool
}
""", process_includes=True)

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""\
mmtbx.validate_ions model.pdb data.mtz [options ...]

Utility to validate ions that have been built into a model, based on local
environment, electron density maps, and atomic properties.

Full parameters:
%s
""" % master_phil.as_str(prefix=" "))
  from mmtbx import ions
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    update_f_part1_for="map",
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=True,
    create_fmodel=True,
    prefer_anomalous=True)
  fmodel = cmdline.fmodel
  xray_structure = cmdline.xray_structure
  params = cmdline.params
  if (params.wavelength is None) :
    from iotbx.file_reader import any_file
    pdb_in = any_file(params.input.pdb.file_name[0],
      force_type="pdb")
    wavelength = pdb_in.file_object.extract_wavelength()
    if (wavelength is not None) :
      print >> out, ""
      print >> out, "Using wavelength = %g from PDB header" % wavelength
      params.wavelength = wavelength
  if (params.wavelength is not None) :
    xray_structure.set_inelastic_form_factors(
      photon=params.wavelength,
      table="sasaki")
    fmodel.update_xray_structure(xray_structure, update_f_calc=True)
  pdb_hierarchy = cmdline.pdb_hierarchy
  geometry = cmdline.geometry
  make_header("Inspecting ions", out=out)
  manager = ions.create_manager(
    pdb_hierarchy = pdb_hierarchy,
    fmodel=fmodel,
    geometry_restraints_manager=geometry,
    wavelength=params.wavelength,
    params=params,
    verbose = params.debug,
    nproc = params.nproc,
    log=out)
  manager.show_current_scattering_statistics(out=out)
  return manager.validate_ions(out = out, debug = params.debug)

if (__name__ == "__main__") :
  run(sys.argv[1:])
