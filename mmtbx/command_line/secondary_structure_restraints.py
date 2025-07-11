"""generate pseudo H-bond restraints for alpha helices, beta sheets, and nucleic acid base pairs"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_restraints

from mmtbx.secondary_structure import sec_str_master_phil_str, \
  sec_str_master_phil, manager
import iotbx.pdb
import mmtbx.model
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
import sys, os

master_phil_str = """
show_all_params = False
  .type = bool
  .style = hidden
filter_outliers = True
  .type = bool
  .style = hidden
format = *phenix phenix_refine phenix_bonds pymol pdb refmac kinemage csv
  .type = choice
  .style = hidden
quiet = False
  .type = bool
  .style = hidden
verbose = -1
  .type = int
  .style = hidden
output_prefix = None
  .type = str
ignore_annotation_in_file = False
  .type = bool
  .help = ignore annotation that is present in input file
file_name = None
  .type = path
  .multiple = True
  .optional = True
  .style = hidden
  %s
""" % sec_str_master_phil_str

master_phil = iotbx.phil.parse(master_phil_str, process_includes=True)

def show_usage():
  help_msg = """\
phenix.secondary_structure_restraints: tool for manipulating with secondary
  structure restraints. It can search and output them in various formats.

Please note, that if HELIX/SHEET records are present in supplied .pdb file,
automatic search will not be executed. These records will be outputted in
choosen format instead.

Usage examples:
  phenix.secondary_structure_restraints model.pdb
  phenix.secondary_structure_restraints model.pdb format=phenix_refine
  phenix.secondary_structure_restraints model.pdb ignore_annotation_in_file=True
  phenix.secondary_structure_restraints model.pdb search_method=from_ca

Full scope of parameters:
  """
  print(help_msg)
  master_phil.show()

def run(args, params=None, out=sys.stdout, log=sys.stderr):
  # params keyword is for running program from GUI dialog
  if ( ((len(args) == 0) and (params is None)) or
       ((len(args) > 0) and ((args[0] == "-h") or (args[0] == "--help"))) ):
    show_usage()
    return

  # parse command-line arguments
  if (params is None):
    pcl = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil_str,
      pdb_file_def="file_name")
    work_params = pcl.work.extract()
  # or use parameters defined by GUI
  else:
    work_params = params
  pdb_files = work_params.file_name

  work_params.secondary_structure.enabled=True
  assert work_params.format in ["phenix", "phenix_refine", "phenix_bonds",
      "pymol", "refmac", "kinemage", "pdb", 'csv']
  if work_params.quiet :
    out = StringIO()

  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_files)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  cs = pdb_structure.crystal_symmetry()

  corrupted_cs = False
  if cs is not None:
    if [cs.unit_cell(), cs.space_group()].count(None) > 0:
      corrupted_cs = True
      cs = None
    elif cs.unit_cell().volume() < 10:
      corrupted_cs = True
      cs = None

  if cs is None:
    if corrupted_cs:
      print("Symmetry information is corrupted, ", file=out)
    else:
      print("Symmetry information was not found, ", file=out)
    print("putting molecule in P1 box.", file=out)
    from cctbx import uctbx
    atoms = pdb_structure.atoms()
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart=atoms.extract_xyz(),
      buffer_layer=3)
    atoms.set_xyz(new_xyz=box.sites_cart)
    cs = box.crystal_symmetry()

  defpars = mmtbx.model.manager.get_default_pdb_interpretation_params()
  defpars.pdb_interpretation.automatic_linking.link_carbohydrates=False
  defpars.pdb_interpretation.c_beta_restraints=False
  defpars.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  defpars.pdb_interpretation.allow_polymer_cross_special_position=True
  # defpars.pdb_interpretation.flip_symmetric_amino_acids=False
  # defpars.pdb_interpretation.sort_atoms=False
  model = mmtbx.model.manager(
      model_input=pdb_structure,
      crystal_symmetry=cs,
      stop_for_unknowns=False)
  model.process(pdb_interpretation_params=defpars)
  pdb_hierarchy = model.get_hierarchy()
  geometry = None
  if pdb_hierarchy.contains_nucleic_acid():
    model.process(make_restraints = True)
    geometry = model.get_restraints_manager().geometry
  if len(pdb_hierarchy.models()) != 1 :
    raise Sorry("Multiple models not supported.")
  ss_from_file = None
  if (hasattr(pdb_structure, "extract_secondary_structure") and
      not work_params.ignore_annotation_in_file):
    ss_from_file = pdb_structure.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    geometry_restraints_manager=geometry,
    sec_str_from_pdb_file=ss_from_file,
    params=work_params.secondary_structure,
    verbose=work_params.verbose)

  # bp_p = nucleic_acids.get_basepair_plane_proxies(
  #     pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.base_pair,
  #     geometry)
  # st_p = nucleic_acids.get_stacking_proxies(
  #     pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.stacking_pair,
  #     geometry)
  # hb_b, hb_a = nucleic_acids.get_basepair_hbond_proxies(pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.base_pair)
  result_out = StringIO()
  # prefix_scope="refinement.pdb_interpretation"
  # prefix_scope=""
  prefix_scope=""
  if work_params.format == "phenix_refine":
    prefix_scope = "refinement.pdb_interpretation"
  elif work_params.format == "phenix":
    prefix_scope = "pdb_interpretation"
  ss_phil = None
  working_phil = m.as_phil_str(master_phil=sec_str_master_phil)
  phil_diff = sec_str_master_phil.fetch_diff(source=working_phil)

  if work_params.format in ["phenix", "phenix_refine"]:
    comment = "\n".join([
      "# These parameters are suitable for use in e.g. phenix.real_space_refine",
      "# or geometry_minimization. To use them in phenix.refine add ",
      "# 'refinement.' if front of pdb_interpretation."])
    if work_params.format == "phenix_refine":
      comment = "\n".join([
      "# These parameters are suitable for use in phenix.refine only.",
      "# To use them in other Phenix tools remove ",
      "# 'refinement.' if front of pdb_interpretation."])
    print(comment, file=result_out)
    if (prefix_scope != ""):
      print("%s {" % prefix_scope, file=result_out)
    if work_params.show_all_params :
      working_phil.show(prefix="  ", out=result_out)
    else :
      phil_diff.show(prefix="  ", out=result_out)
    if (prefix_scope != ""):
      print("}", file=result_out)
  elif work_params.format == "pdb":
    if m.actual_sec_str.fits_in_pdb_format():
      print(m.actual_sec_str.as_pdb_str(), file=result_out)
    else:
      raise Sorry("Annotations could not fit in PDB format.")
  elif work_params.format == "phenix_bonds" :
    raise Sorry("Not yet implemented.")
  elif work_params.format in ["pymol", "refmac", "kinemage", 'csv'] :
    m.show_summary(log=out)
    model.process(make_restraints=True)
    (hb_proxies, hb_angle_proxies, planarity_proxies,
        parallelity_proxies) = m.create_all_new_restraints(
        pdb_hierarchy=pdb_hierarchy,
        grm=model.get_restraints_manager().geometry,
        log=out)
    if hb_proxies.size() > 0:
      if work_params.format == "pymol" :
        file_load_add = "load %s" % work_params.file_name[0]
        # surprisingly, pymol handles filenames with whitespaces without quotes...
        print(file_load_add, file=result_out)
        bonds_in_format = hb_proxies.as_pymol_dashes(
            pdb_hierarchy=pdb_hierarchy)
      elif work_params.format == "kinemage" :
        bonds_in_format = hb_proxies.as_kinemage(
            pdb_hierarchy=pdb_hierarchy)
      elif work_params.format == "csv" :
        bonds_in_format = hb_proxies.as_csv(
            pdb_hierarchy=pdb_hierarchy)
      else :
        bonds_in_format = hb_proxies.as_refmac_restraints(
            pdb_hierarchy=pdb_hierarchy)
      print(bonds_in_format, file=result_out)
    if hb_angle_proxies.size() > 0:
      if work_params.format == "pymol":
        angles_in_format = hb_angle_proxies.as_pymol_dashes(
            pdb_hierarchy=pdb_hierarchy)
        print(angles_in_format, file=result_out)
  result = result_out.getvalue()
  out_prefix = os.path.basename(work_params.file_name[0])
  if work_params.output_prefix is not None:
    out_prefix = work_params.output_prefix
  filename = "%s_ss.eff" % out_prefix
  if work_params.format == "pymol":
    filename = "%s_ss.pml" % out_prefix
  outf = open(filename, "w")
  outf.write(result)
  outf.close()
  print(result, file=out)

  return os.path.abspath(filename)

# =============================================================================
# GUI-specific functions
def validate_params(params):
  pass

if __name__ == "__main__" :
  run(sys.argv[1:])

