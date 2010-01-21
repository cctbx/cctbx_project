# LIBTBX_SET_DISPATCHER_NAME phenix.fmodel

import sys, os
import mmtbx.utils
import iotbx.phil
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import Sorry
import cctbx.xray

legend = """
phenix.fmodel: a tool to compute structure factors, Fmodel:

  Fmodel = scale * exp(AnisoScale) * (Fcalc + k_sol * exp(-b_sol*s^2/4) * Fmask)

  where:

  - Fmodel - total model structure factor (complex value)
  - AnisoScale = -ht*A(-1)*b_cart*A(-1)th/4
  - h - column vector with Miller indices
  - A - orthogonalization matrix
  - b_cart - anisotropic scale matrix
  - t and (-1) denotes transposition and inversion operations
  - scale - overall scale factor
  - Fcalc - structure factors calculated from atomic model
  - k_sol and b_sol - Flat Bulk solvent model parameters
  - Fmask - structure factors calculated from bulk solvent mask

Usage examples:

  1) phenix.fmodel model.pdb high_resolution=1.5

     will result in a file containing complete set of Fmodel = Fcalc computed
     from atomic model up to 1.5A resolution.

  2) phenix.fmodel model.pdb scale=2 k_sol=0.35 b_sol=50 b_cart="1 2 3 0 4 7" high_res=1.5 low_res=10

     will result in a file containing complete set of Fmodel computed using the
     above formula in resolution range 1.5-20.0A.

  3) phenix.fmodel model.pdb high_resolution=1.5 algorithm=direct

     is similar to "1)" but the Fcalc are computed using direct summation algorithm.

  4) phenix.fmodel model.pdb high_res=1.5 format=cns label=FOBS type=real r_free=0.1

     will result in CNS formatted file containing complete set of amplitudes of
     Fmodel = Fcalc computed up to 1.5A resolution, labelled as FOBS, and free-R
     flags with 10% of test reflections. This is a typical command to simulate Fobs.

  5) phenix.fmodel model.pdb high_res=1.5 scattering_table=neutron

     will result in a file containing complete set of Fmodel = Fcalc computed
     from atomic model up to 1.5A resolution using neutron scattering table.

  6) phenix.fmodel model.pdb parameters.txt

     will result in a structure factor file, where Fmodel were computed using
     parameters defined in parameters.txt file. The parameters.txt file can
     contain all or any subset of parameters listed below. Note, that each {
     must have a matching one }.

See below for complete list of available parameters.
"""

fmodel_from_xray_structure_params_str = """\
fmodel
{
  k_sol = 0.0
    .type = float
    .help = Bulk solvent k_sol values
    .short_caption=Bulk solvent K_sol value
    .expert_level=1
  b_sol = 0.0
    .type = float
    .help = Bulk solvent b_sol values
    .short_caption=Bulk solvent B_sol value
    .expert_level=1
  b_cart = 0 0 0 0 0 0
    .type = floats(6)
    .help = Anisotropic scale matrix
    .short_caption = Anisotropic scale matrix
    .expert_level=1
  scale = 1.0
    .type = float
    .help = Overall scale factor
    .expert_level=1
}
structure_factors_accuracy
  .short_caption = Structure factors accuracy
  .style = menu_item auto_align noauto parent_submenu:advanced
{
  include scope mmtbx.f_model.sf_and_grads_accuracy_master_params
}
mask
  .short_caption = Bulk solvent mask
  .style = menu_item auto_align noauto parent_submenu:advanced
{
  include scope mmtbx.masks.mask_master_params
}
"""
fmodel_from_xray_structure_params = iotbx.phil.parse(
  fmodel_from_xray_structure_params_str, process_includes=True)

fmodel_from_xray_structure_master_params_str = """\
high_resolution = None
  .type = float
  .expert_level=1
low_resolution = None
  .type = float
  .expert_level=1
r_free_flags_fraction = None
  .type = float
  .expert_level=1
add_sigmas = False
  .type = bool
  .expert_level=1
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
  .help = Choices of scattering table for structure factors calculations
  .expert_level=1
%s
output
  .short_caption = Reflection output
  .expert_level=0
{
  format = *mtz cns
    .type = choice
    .short_caption = File format
  label = FMODEL
    .type = str
    .short_caption = Data label
    .expert_level=1
  type = real *complex
    .type = choice
    .short_caption = Data type
    .expert_level=1
  file_name = None
    .type = path
    .expert_level=1
    .short_caption = Output reflections file
    .style = bold noauto new_file
}

anomalous_scatterers
  .short_caption = Anomalous sites
  .style = menu_item scrolled parent_submenu:Atom_selections
{
  group
    .optional = True
    .multiple = True
    .short_caption = Anomalous scatterer group
    .style = auto_align
  {
    selection = None
      .type = str
      .short_caption = Atom selection
      .input_size = 400
      .style = selection
    f_prime = 0
      .type = float
      .short_caption = f'
    f_double_prime = 0
      .type = float
      .short_caption = f''
  }
}
"""%fmodel_from_xray_structure_params_str

fmodel_from_xray_structure_master_params = iotbx.phil.parse(
  fmodel_from_xray_structure_master_params_str, process_includes=True)

def set_fp_fdp_for_anomalous_scatterers(pdb_hierarchy, xray_structure,
  anomalous_scatterer_groups):
  scatterers = xray_structure.scatterers()
  for group in anomalous_scatterer_groups:
    iselection = pdb_hierarchy.atom_selection_cache().selection(
      string = group.selection).iselection()
    if(iselection.size() == 0):
      raise Sorry(
        "Empty selection: selection string '%s' does not select any atom."%
        group.selection)
    for i_seq in iselection:
      scatterers[i_seq].fp = group.f_prime
      scatterers[i_seq].fdp = group.f_double_prime

def run(args, log = sys.stdout):
  print >> log, legend
  processed_args = mmtbx.utils.process_command_line_args(args = args, log = log,
    master_params = fmodel_from_xray_structure_master_params)
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(
    file_names = processed_args.pdb_file_names)
  pdb_combined.report_non_unique(out = log)
  print >> log, "-"*79
  print >> log, "\nParameters to compute Fmodel::\n"
  processed_args.params.show(out = log, prefix=" ")
  params = processed_args.params.extract()
  if(len(pdb_combined.unique_file_names) == 0): return
  print >> log, "-"*79
  print >> log, "\nInput PDB file(s):", " ".join(processed_args.pdb_file_names)
  pdb_inp = iotbx.pdb.input(source_info = None,
    lines = flex.std_string(pdb_combined.raw_records))
  cryst1 = pdb_inp.crystal_symmetry_from_cryst1()
  if(cryst1 is None):
    raise Sorry("CRYST1 record in input PDB file is incomplete or missing.")
  else:
    if([cryst1.unit_cell(), cryst1.space_group_info()].count(None) != 0):
      raise Sorry("CRYST1 record in input PDB file is incomplete or missing.")
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.show_summary(f = log)
  if(len(params.anomalous_scatterers.group) != 0):
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    set_fp_fdp_for_anomalous_scatterers(
      pdb_hierarchy              = pdb_hierarchy,
      xray_structure             = xray_structure,
      anomalous_scatterer_groups = params.anomalous_scatterers.group)
  print >> log, "-"*79
  print >> log, "Computing model structure factors, Fmodel:"
  if(params.output.format == "cns"): extension = ".hkl"
  elif(params.output.format == "mtz"): extension = ".mtz"
  ofn = params.output.file_name
  if(ofn is None):
    ofn = os.path.basename(processed_args.pdb_file_names[0])
    if(len(processed_args.pdb_file_names)==1): ofn = ofn + extension
    else: ofn = ofn + "_et_al" + extension
  mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xray_structure,
    add_sigmas     = params.add_sigmas,
    params         = params).write_to_file(file_name = ofn)
  print >> log, "Output file name:", ofn
  print >> log, "All done."
  print >> log, "-"*79


if (__name__ == "__main__"):
  run(args=sys.argv[1:])
