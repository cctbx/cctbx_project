from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_restraints

from mmtbx.secondary_structure import sec_str_master_phil_str, manager
import iotbx.pdb
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import Sorry
import cStringIO
import sys, os

def run (args, out=sys.stdout, log=sys.stderr) :
  pdb_files = []
  sources = []
  force_new_annotation = False
  master_phil_str = """
show_all_params = False
  .type = bool
filter_outliers = True
  .type = bool
format = *phenix phenix_refine phenix_bonds pymol refmac kinemage
  .type = choice
quiet = False
  .type = bool
verbose = -1
  .type = int
file_name = None
  .type = path
  .multiple = True
  .optional = True
  .style = hidden
  %s
""" % sec_str_master_phil_str
  pcl = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="file_name")

  from mmtbx.secondary_structure import sec_str_master_phil
  work_params = pcl.work.extract()
  work_params.secondary_structure.enabled=True
  assert work_params.format in ["phenix", "phenix_refine", "phenix_bonds",
      "pymol", "refmac", "kinemage"]
  if work_params.quiet :
    out = cStringIO.StringIO()
  pdb_files = work_params.file_name
  if len(pdb_files) > 0 :
    work_params.file_name.extend(pdb_files)
  if len(pdb_files) == 0 :
    raise Sorry("No PDB files specified.")
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_files)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  cs = pdb_structure.crystal_symmetry()

  corrupted_cs = False
  if cs is not None:
    if [cs.unit_cell(), cs.space_group()].count(None) > 0:
      corrupted_cs = True
      cs = None

  if cs is None:
    if corrupted_cs:
      print >> out, "Symmetry information is corrupted, "
    else:
      print >> out, "Symmetry information was not found, "
    print >> out, "putting molecule in P1 box."
    from cctbx import uctbx
    atoms = pdb_structure.atoms()
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart=atoms.extract_xyz(),
      buffer_layer=3)
    atoms.set_xyz(new_xyz=box.sites_cart)
    cs = box.crystal_symmetry()
  from mmtbx.monomer_library import pdb_interpretation, server
  import mmtbx
  import mmtbx.command_line.geometry_minimization

  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  defpars = mmtbx.command_line.geometry_minimization.master_params().extract()
  defpars.pdb_interpretation.automatic_linking.link_carbohydrates=False
  defpars.pdb_interpretation.c_beta_restraints=False
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    pdb_inp        = pdb_structure,
    crystal_symmetry = cs,
    params         = defpars.pdb_interpretation,
    force_symmetry = True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  geometry = processed_pdb_file.geometry_restraints_manager()
  geometry.pair_proxies(processed_pdb_file.xray_structure().sites_cart())
  pdb_hierarchy.atoms().reset_i_seq()
  if len(pdb_hierarchy.models()) != 1 :
    raise Sorry("Multiple models not supported.")
  m = manager(pdb_hierarchy=pdb_hierarchy,
    geometry_restraints_manager=geometry,
    sec_str_from_pdb_file=pdb_structure.extract_secondary_structure(),
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
  result_out = cStringIO.StringIO()
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
      "# or geometry_minimization. To use theim in phenix.refine add ",
      "# 'refinement.' if front of pdb_interpretation."])
    if work_params.format == "phenix_refine":
      comment = "\n".join([
      "# These parameters are suitable for use in phenix.refine only.",
      "# To use theim in other Phenix tools remove ",
      "# 'refinement.' if front of pdb_interpretation."])
    print >> result_out, comment
    if (prefix_scope != "") :
      print >> result_out, "%s {" % prefix_scope
    if work_params.show_all_params :
      working_phil.show(prefix="  ", out=result_out)
    else :
      phil_diff.show(prefix="  ", out=result_out)
    if (prefix_scope != "") :
      print >> result_out, "}"

  elif work_params.format == "phenix_bonds" :
    raise Sorry("Not yet implemented.")
  elif work_params.format in ["pymol", "refmac", "kinemage"] :
    m.show_summary(log=out)
    (hb_proxies, hb_angle_proxies, planarity_proxies,
        parallelity_proxies) = m.create_all_new_restraints(
        pdb_hierarchy=pdb_hierarchy,
        grm=geometry,
        log=out)
    if hb_proxies.size() > 0:
      if work_params.format == "pymol" :
        bonds_in_format = hb_proxies.as_pymol_dashes(
            pdb_hierarchy=pdb_hierarchy)
      elif work_params.format == "kinemage" :
        bonds_in_format = hb_proxies.as_kinemage(
            pdb_hierarchy=pdb_hierarchy)
      else :
        bonds_in_format = hb_proxies.as_refmac_restraints(
            pdb_hierarchy=pdb_hierarchy)
      print >> result_out, bonds_in_format
  result = result_out.getvalue()
  outf = open("%s_ss.eff" %os.path.basename(work_params.file_name[0]), "w")
  outf.write(result)
  outf.close()
  print >> out, result

if __name__ == "__main__" :
  run(sys.argv[1:])
