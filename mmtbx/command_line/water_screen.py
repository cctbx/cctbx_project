
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header
from libtbx.utils import Sorry
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
include scope mmtbx.ions.svm.svm_phil_str
debug = True
  .type = bool
elements = Auto
  .type = str
use_svm = False
  .type = bool
nproc = Auto
  .type = int
""")

def run(args, out=sys.stdout):
  usage_string = """
mmtbx.water_screen model.pdb data.mtz [options ...]

Utility to flag waters that may actually be elemental ions, based on local
environment, electron density maps, and atomic properties.
"""
  import mmtbx.ions.identify
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=True,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    create_fmodel=True,
    prefer_anomalous=True)
  fmodel = cmdline.fmodel
  xray_structure = cmdline.xray_structure
  params = cmdline.params
  if (params.use_svm):
    if (params.elements is Auto):
      raise Sorry("You must specify elements to consider when using the SVM "+
        "prediction method.")
  pdb_hierarchy = cmdline.pdb_hierarchy
  geometry = cmdline.geometry
  make_header("Inspecting water molecules", out=out)
  manager_class = None
  if (params.use_svm):
    manager_class = mmtbx.ions.svm.manager
  manager = mmtbx.ions.identify.create_manager(
    pdb_hierarchy = pdb_hierarchy,
    fmodel = fmodel,
    geometry_restraints_manager = geometry,
    wavelength = params.input.wavelength,
    params = params,
    verbose = params.debug,
    nproc = params.nproc,
    log = out,
    manager_class = manager_class)
  manager.show_current_scattering_statistics(out=out)
  candidates = Auto
  if (params.elements is not Auto) and (params.elements is not None):
    from cctbx.eltbx import chemical_elements
    lu = chemical_elements.proper_upper_list()
    elements = params.elements.replace(",", " ")
    candidates = elements.split()
    for elem in candidates :
      if (elem.upper() not in lu):
        raise Sorry("Unrecognized element '%s'" % elem)
  results = manager.analyze_waters(
    out = out,
    debug = params.debug,
    candidates = candidates)
  return results, pdb_hierarchy

if (__name__ == "__main__"):
  run(sys.argv[1:])
