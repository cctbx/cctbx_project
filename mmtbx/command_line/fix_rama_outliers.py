from __future__ import division
import iotbx.pdb
import sys
from scitbx.array_family import flex
from mmtbx.rotamer import ramachandran_eval
import mmtbx.utils
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx.monomer_library import idealized_aa
from libtbx.utils import Sorry, null_out
from mmtbx.command_line.geometry_minimization import get_geometry_restraints_manager
from mmtbx.validation.ramalyze import RAMALYZE_OUTLIER
from mmtbx.validation.ramalyze import RAMALYZE_ALLOWED # import dependency
from mmtbx.validation.ramalyze import RAMALYZE_OUTLIER # import dependency
from mmtbx.building.loop_closure.ccd import ccd_python
from mmtbx.building.loop_closure import utils, starting_conformations
from mmtbx.pdbtools import truncate_to_poly_gly
from mmtbx.secondary_structure.build import side_chain_placement, \
    set_xyz_smart

loop_idealization_aster_phil_str = """
loop_idealization
{
  enabled = True
    .type = bool
  change_non_rama_outliers = True
    .type = bool
    .help = Allow changing non-outlier ramachandran residues angles
  output_prefix = bchain_2
    .type = str
}
"""

master_phil = iotbx.phil.parse(loop_idealization_aster_phil_str)

class loop_idealization():
  def __init__(self, pdb_hierarchy, params=None, log=null_out(), verbose=True):
    self.original_pdb_h = pdb_hierarchy
    self.params = self.process_params(params)
    self.log = log
    self.verbose = verbose

  def process_params(self, params):
    if params is None:
      params = master_phil.fetch().extract()
      params.loop_idealization.enabled = True
    if hasattr(params, "model_idealization"):
      p_pars = params.model_idealization
    else:
      assert hasattr(params, "enabled") and hasattr(params, "change_non_rama_outliers"), \
          "Something wrong with parameters passed to model_idealization"
      p_pars = params

    assert isinstance(p_pars.enabled, bool)
    assert isinstance(p_pars.change_non_rama_outliers, bool)
    return p_pars


def get_main_chain_rmsd_range(
    hierarchy, original_h, all_atoms=False, start_res_num=None, end_res_num=None):
  rmsd = 0
  mc_atoms = None
  if all_atoms:
    mc_atoms = ["N", "CA", "C", "O"]
  else:
    mc_atoms = ["N", "CA", "C"]
  for m_atom, ref_atom in zip(hierarchy.atoms(), original_h.atoms()):
    if m_atom.name.strip() in mc_atoms:
      if (start_res_num is None or
          end_res_num is None or
          m_atom.parent().parent().resseq_as_int() in range(start_res_num, end_res_num+1)):
        rmsd += m_atom.distance(ref_atom)**2
  return rmsd**0.5


def minimize_hierarchy(hierarchy, xrs, original_pdb_h,
    list_of_exclusion_resnums, log=None):
  from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
  from mmtbx.refinement.geometry_minimization import run2
  from mmtbx.geometry_restraints import reference

  if log is None:
    log = null_out()
  params_line = grand_master_phil_str
  params = iotbx.phil.parse(
      input_string=params_line, process_includes=True).extract()
  params.pdb_interpretation.peptide_link.ramachandran_restraints = True
  params.pdb_interpretation.c_beta_restraints=True

  processed_pdb_files_srv = mmtbx.utils.\
      process_pdb_file_srv(
          crystal_symmetry= xrs.crystal_symmetry(),
          pdb_interpretation_params = params.pdb_interpretation,
          log=log,
          cif_objects=None)
  processed_pdb_file, junk = processed_pdb_files_srv.\
      process_pdb_files(raw_records=flex.split_lines(hierarchy.as_pdb_string()))
  grm = get_geometry_restraints_manager(
      processed_pdb_file, xrs)

  selection = ""
  for resnum in list_of_exclusion_resnums:
    selection += " and not resseq %d" % resnum
  selection = "name N or name CA or name C or name O " + selection
  print "selection for reference restraints:", selection
  asc = original_pdb_h.atom_selection_cache()
  sel = asc.selection(selection)


  grm.geometry.append_reference_coordinate_restraints_in_place(
      reference.add_coordinate_restraints(
          sites_cart = original_pdb_h.atoms().extract_xyz().select(sel),
          selection  = sel,
          sigma      = 0.5))
  obj = run2(
      restraints_manager       = grm,
      pdb_hierarchy            = hierarchy,
      correct_special_position_tolerance = 1.0,
      max_number_of_iterations = 300,
      number_of_macro_cycles   = 5,
      bond                     = True,
      nonbonded                = True,
      angle                    = True,
      dihedral                 = True,
      chirality                = True,
      planarity                = True,
      fix_rotamer_outliers     = True,
      log                      = log)


def place_side_chains(hierarchy, original_h,
    rotamer_manager, start_res_num, end_res_num):
  ideal_res_dict = idealized_aa.residue_dict()
  asc = original_h.atom_selection_cache()
  gly_atom_names = set([" N  ", " CA ", " C  ", " O  "])
  for rg in hierarchy.residue_groups():
    if rg.resseq_as_int() in range(start_res_num, end_res_num+1):
      # cut extra atoms
      ag = rg.only_atom_group()
      for atom in ag.atoms():
        if (atom.name not in gly_atom_names):
          ag.remove_atom(atom=atom)
      # get ag from original hierarchy

      orig_ag = original_h.select(asc.selection("resseq %d" % rg.resseq_as_int())
          ).models()[0].chains()[0].residue_groups()[0].atom_groups()[0]
      # get ideal
      ideal_ag = ideal_res_dict[ag.resname.lower()].models()[0].chains()[0].\
        residue_groups()[0].atom_groups()[0]
      # print "got to placement"
      side_chain_placement(ag, orig_ag, rotamer_manager)

def get_fixed_moving_parts(pdb_hierarchy, start_res_num, end_res_num):
  original_pdb_h = pdb_hierarchy.deep_copy()
  xrs = original_pdb_h.extract_xray_structure()
  truncate_to_poly_gly(pdb_hierarchy, start_res_num, end_res_num)
  cache = pdb_hierarchy.atom_selection_cache()
  # print "selectioin:", "resid %d through %d" % (start_res_num, end_res_num)
  m_selection = cache.selection("resid %d through %d" % (start_res_num, end_res_num))
  moving_h = pdb_hierarchy.select(m_selection)
  moving_h.reset_atom_i_seqs()
  # print dir(moving_h)
  # STOP()
  m_cache = moving_h.atom_selection_cache()
  # print "len inp h atoms", pdb_hierarchy.atoms().size()
  # print "len moving_h atoms", moving_h.atoms().size()
  moving_ref_atoms_iseqs = []
  # here we need N, CA, C atoms from the end_res_num residue
  sel = m_cache.selection("resid %d and name N" % end_res_num)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_N = a.detached_copy()

  sel = m_cache.selection("resid %d and name CA" % end_res_num)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_CA = a.detached_copy()

  sel = m_cache.selection("resid %d and name C" % end_res_num)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_C = a.detached_copy()

  fixed_ref_atoms = [fixed_N, fixed_CA, fixed_C]

  return moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms

def fix_rama_outlier(
    pdb_hierarchy, out_res_num, prefix="", minimize=True):
  original_pdb_h = pdb_hierarchy.deep_copy()
  rotamer_manager = RotamerEval()
  for ccd_radius, change_all, change_radius in [
      (1, False, 0),
      (2, False, 0),
      # (3, False, 0),
      (2, True, 1),
      # (3, True, 1),
      ]:
  # while ccd_radius <= 3:
    print "  Starting optimization with radius, change_all:", ccd_radius, change_all
    moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms = get_fixed_moving_parts(
        pdb_hierarchy=pdb_hierarchy,
        start_res_num=out_res_num-ccd_radius,
        end_res_num=out_res_num+ccd_radius)
    moving_h_set = None
    if change_all:
      moving_h_set = starting_conformations.get_all_starting_conformations(moving_h, change_radius)
    else:
      moving_h_set = starting_conformations.get_starting_conformations(moving_h)

    if len(moving_h_set) == 0:
      # outlier was fixed before somehow...
      return original_pdb_h

    rotamer_manager = RotamerEval()
    for i, h in enumerate(moving_h_set):
      ccd_obj = ccd_python(fixed_ref_atoms, h, moving_ref_atoms_iseqs)
      ccd_obj.run()
      resulting_rmsd = ccd_obj.resulting_rmsd
      states = ccd_obj.states
      n_iter = ccd_obj.n_iter

      # resulting_rmsd, states, n_iter = ccd(
      #     fixed_ref_atoms, h, moving_ref_atoms_iseqs, moving_h)

      mc_rmsd = get_main_chain_rmsd_range(moving_h, h, all_atoms=True)
      print "Resulting anchor and backbone RMSDs, n_iter for model %d:" % i,
      print resulting_rmsd, ",", mc_rmsd, ",", n_iter
      #
      # setting new coordinates
      #
      moved_with_side_chains_h = pdb_hierarchy.deep_copy()
      set_xyz_smart(moved_with_side_chains_h, h)
      #
      # placing side-chains
      #
      # moved_with_side_chains_h.write_pdb_file(
      #     file_name="%s_before_sc_placement_%d.pdb" % (prefix, i))
      place_side_chains(moved_with_side_chains_h, original_pdb_h,
          rotamer_manager, out_res_num-ccd_radius, out_res_num+ccd_radius)
      # moved_with_side_chains_h.write_pdb_file(
      #     file_name="%s_after_sc_placement_%d.pdb" % (prefix, i))


      #
      # finalizing with geometry_minimization
      #
      if mc_rmsd < 3:
        if minimize:
          print "minimizing..."
          moved_with_side_chains_h.write_pdb_file(
              file_name="%s_result_before_min_%d.pdb" % (prefix, i))
          minimize_hierarchy(moved_with_side_chains_h, xrs, original_pdb_h)
        moved_with_side_chains_h.write_pdb_file(
            file_name="%s_result_minimized_%d.pdb" % (prefix, i))
        final_rmsd = get_main_chain_rmsd_range(moved_with_side_chains_h,
            original_pdb_h, out_res_num-ccd_radius, out_res_num+ccd_radius)
        print "FINAL RMSD after minimization:", final_rmsd
        return moved_with_side_chains_h
    ccd_radius += 1


  print "Epic FAIL: failed to fix rama outlier"
  return original_pdb_h






def get_resnums_of_chain_rama_outliers(pdb_hierarchy, r):
  phi_psi_atoms = utils.get_phi_psi_atoms(pdb_hierarchy)
  result = []
  rama_results = []
  ranges_for_idealization = []
  print "rama outliers for input hierarchy:"
  list_of_reference_exclusion = []
  utils.list_rama_outliers_h(pdb_hierarchy, r)
  for phi_psi_pair, rama_key in phi_psi_atoms:
    ev = utils.rama_evaluate(phi_psi_pair, r, rama_key)
    rama_results.append(ev)
    if ev == RAMALYZE_OUTLIER:
      resnum = phi_psi_pair[0][2].parent().parent().resseq_as_int()
      result.append(resnum)
  return result

def get_list_of_exclusions_for_reference(resnums_of_rama_outliers):
  result = []
  for resnum in resnums_of_rama_outliers:
    result += [resnum-1, resnum, resnum+1]
  return result

def idealize_chain(pdb_hierarchy, prefix, minimize_whole=True):
  occ_groups = pdb_hierarchy.occupancy_groups_simple()
  if len(occ_groups) > 0:
    raise Sorry("Alternative conformations are not supported.")
  xrs = pdb_hierarchy.extract_xray_structure()
  working_h = pdb_hierarchy.deep_copy()
  working_h.reset_atom_i_seqs()
  r = ramachandran_eval.RamachandranEval()
  # phi_psi_atoms = get_phi_psi_atoms(working_h)
  rama_results = []
  ranges_for_idealization = []
  print "rama outliers for input hierarchy:"
  rama_out_resnums = get_resnums_of_chain_rama_outliers(working_h, r)
  list_of_reference_exclusion = get_list_of_exclusions_for_reference(
      rama_out_resnums)

  # list_of_reference_exclusion = []
  # list_rama_outliers_h(working_h, r)
  # for phi_psi_pair, rama_key in phi_psi_atoms:
  #   ev = rama_evaluate(phi_psi_pair, r, rama_key)
  #   rama_results.append(ev)
  #   if ev == "OUTLIER":
  #     resnum = phi_psi_pair[0][2].parent().parent().resseq_as_int()
  #     ranges_for_idealization.append((resnum-2, resnum+2))
  #     list_of_reference_exclusion += [resnum-1, resnum, resnum+1]
  # outlier_indices = [i for i,x in enumerate(rama_results) if x == "OUTLIER"]
  # print outlier_indices
  # n_outliers = len(outlier_indices)
  # print "Ranges for idealization:", ranges_for_idealization
  out_i = 0
  for rama_out_resnum in rama_out_resnums:
    print
    print "Fixing outlier:", rama_out_resnum
    # fn = "%s_before_ccd_%d.pdb" % (prefix, out_i)
    # print "  writing file %s" % fn
    # working_h.write_pdb_file(file_name="%s_before_ccd_%d.pdb" % (prefix, out_i))
    new_h = fix_rama_outlier(
      pdb_hierarchy=working_h,
      out_res_num=rama_out_resnum,
      prefix=prefix,
      minimize=False)
    print "listing outliers after loop minimization"
    utils.list_rama_outliers_h(new_h, r)
    fn = "%s_after_loop_%d.pdb" % (prefix, out_i)
    print "  writing file %s" % fn
    new_h.write_pdb_file(file_name=fn)
    working_h = new_h
    out_i += 1
  if minimize_whole:
    print "minimizing whole thing..."
    minimize_hierarchy(new_h, xrs, pdb_hierarchy,list_of_reference_exclusion, log=None)
    new_h.write_pdb_file(file_name="%s_all_minized.pdb" % prefix)

def run(args):
  # print "args", args


  log = sys.stdout
  inputs = mmtbx.utils.process_command_line_args(args=args,
      master_params=master_phil)
  work_params = inputs.params.extract()
  inputs.params.show(prefix=" ", out=log)
  pdb_file_names = list(inputs.pdb_file_names)
  if len(pdb_file_names) == 0:
    raise Sorry("No PDB file specified")
  work_params.loop_idealization.enabled=True
  # work_params.model_idealization.file_name_before_regularization="before.pdb"
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  pdb_h = pdb_input.construct_hierarchy()


  loop_ideal = loop_idealization(pdb_h, work_params.loop_idealization, log)

  pdb_inp = iotbx.pdb.input(args[0])
  pdb_h = pdb_inp.construct_hierarchy()
  prefix = ""
  minimize = "minimize" in args
  if len(args) > 3:
    prefix = args[3]
  else:
    prefix = args[0][:-4]
  print "prefix=", prefix
  idealize_chain(pdb_h, prefix)

if (__name__ == "__main__"):
  run(sys.argv[1:])
