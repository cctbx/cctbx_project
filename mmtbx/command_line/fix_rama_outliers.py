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
from mmtbx.validation import ramalyze
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
  output_prefix = rama_fixed
    .type = str
  minimize_whole = True
    .type = bool
}
"""

master_phil = iotbx.phil.parse(loop_idealization_aster_phil_str)

class loop_idealization():
  def __init__(self, pdb_hierarchy, params=None, log=null_out(), verbose=True):
    if len(pdb_hierarchy.models()) > 1:
      raise Sorry("Multi-model files are not supported")
    self.original_pdb_h = pdb_hierarchy
    xrs = pdb_hierarchy.extract_xray_structure()
    asc = pdb_hierarchy.atom_selection_cache()
    self.resulting_pdb_h = pdb_hierarchy.deep_copy()
    self.params = self.process_params(params)
    self.log = log
    self.verbose = verbose
    self.r = ramachandran_eval.RamachandranEval()
    ram = ramalyze.ramalyze(pdb_hierarchy=pdb_hierarchy)
    self.p_initial_rama_outliers = ram.out_percent
    self.p_before_minimization_rama_outliers = None
    self.p_after_minimiaztion_rama_outliers = None
    self.ref_exclusion_selection = ""
    for chain in pdb_hierarchy.only_model().chains():
      print >> self.log, "Idealizing chain %s" % chain.id
      selection = "protein and chain %s and (name N or name CA or name C or name O)" % chain.id
      sel = asc.selection("chain %s" % chain.id)
      chain_h = self.original_pdb_h.select(sel)
      m = chain_h.only_model()
      i = 0
      cutted_chain_h = None
      for c in m.chains():
        if i == 0:
          cutted_chain_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(c)
        else:
          print >> self.log, "WARNING!!! Duplicating chain ids! Only the first chain will be processed."
          print >> self.log, "  Removing chain %s with %d residues" % (c.id, len(c.residues()))
          m.remove_chain(c)
        i += 1
      exclusions, ch_h = self.idealize_chain(hierarchy=(cutted_chain_h if cutted_chain_h else chain_h))
      if ch_h is not None:
        set_xyz_smart(self.resulting_pdb_h, ch_h)
        for resnum in exclusions:
          selection += " and not resseq %s" % resnum
      self.ref_exclusion_selection += "(%s) or " % selection
    if len(self.ref_exclusion_selection) > 0:
      self.ref_exclusion_selection = self.ref_exclusion_selection[:-3]
    self.resulting_pdb_h.write_pdb_file(file_name="%s_before_minization.pdb" % self.params.output_prefix)
    ram = ramalyze.ramalyze(pdb_hierarchy=self.resulting_pdb_h)
    self.p_before_minimization_rama_outliers = ram.out_percent
    if self.params.minimize_whole:
      print >> self.log, "minimizing whole thing..."
      print >> self.log, "self.ref_exclusion_selection", self.ref_exclusion_selection
      minimize_hierarchy(self.resulting_pdb_h, xrs, self.original_pdb_h, self.ref_exclusion_selection, log=None)
      # self.resulting_pdb_h.write_pdb_file(file_name="%s_all_minized.pdb" % self.params.output_prefix)
      ram = ramalyze.ramalyze(pdb_hierarchy=self.resulting_pdb_h)
      self.p_after_minimiaztion_rama_outliers = ram.out_percent
    # return new_h      

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

  def idealize_chain(self, hierarchy):
    # check no ac:
    for c in hierarchy.chains():
      if len(c.conformers()) > 1:
        raise Sorry("Alternative conformations are not supported.")
    working_h = hierarchy.deep_copy()
    working_h.reset_atom_i_seqs()
    rama_results = []
    ranges_for_idealization = []
    print >> self.log, "rama outliers for input hierarchy:"
    rama_out_resnums = self.get_resnums_of_chain_rama_outliers(
        working_h)
    if len(rama_out_resnums) == 0:
      return None, None
    # get list of residue numbers that should be excluded from reference
    list_of_reference_exclusion = []
    for resnum in rama_out_resnums:
      start_rn, prev_rn = get_res_nums_around(hierarchy, resnum, 1, 1)
      list_of_reference_exclusion += [start_rn, resnum, prev_rn]
    out_i = 0
    for rama_out_resnum in rama_out_resnums:
      print >> self.log
      print >> self.log, "Fixing outlier:", rama_out_resnum
      new_h = self.fix_rama_outlier(
        pdb_hierarchy=working_h,
        out_res_num=rama_out_resnum,
        prefix=self.params.output_prefix,
        minimize=False)
      print >> self.log, "listing outliers after loop minimization"
      outp = utils.list_rama_outliers_h(new_h, self.r)
      print >> self.log, outp
      fn = "%s_after_loop_%d.pdb" % (self.params.output_prefix, out_i)
      # print >> self.log, "  writing file %s" % fn
      # new_h.write_pdb_file(file_name=fn)
      working_h = new_h
      out_i += 1
    return list_of_reference_exclusion, new_h

  def ccd_solution_is_ok(self, 
      anchor_rmsd, mc_rmsd, ccd_radius, change_all_angles, change_radius):
    adaptive_mc_rmsd = {1:3.0, 2:3.5, 3:4.0}
    return mc_rmsd < adaptive_mc_rmsd[ccd_radius] and anchor_rmsd < 0.3


  def fix_rama_outlier(self, 
      pdb_hierarchy, out_res_num, prefix="", minimize=True):
    

    original_pdb_h = pdb_hierarchy.deep_copy()
    rotamer_manager = RotamerEval()
    all_results = []
    for ccd_radius, change_all, change_radius in [
        (1, False, 0), 
        (2, False, 0), 
        # (3, False, 0), 
        (2, True, 1),
        # (3, True, 1),
        ]:
    # while ccd_radius <= 3:
      print >> self.log, "  Starting optimization with radius, change_all:", ccd_radius, change_all
      moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms = get_fixed_moving_parts(
          pdb_hierarchy=pdb_hierarchy,
          out_res_num=out_res_num, 
          n_following=ccd_radius, 
          n_previous=ccd_radius)
      moving_h_set = None
      if change_all:
        moving_h_set = starting_conformations.get_all_starting_conformations(moving_h, change_radius, log=self.log)
      else:
        moving_h_set = starting_conformations.get_starting_conformations(moving_h, log=self.log)

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
        print >> self.log, "Resulting anchor and backbone RMSDs, n_iter for model %d:" % i, 
        print >> self.log, resulting_rmsd, ",", mc_rmsd, ",", n_iter
        all_results.append((h.deep_copy(), mc_rmsd, resulting_rmsd, n_iter))
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
        placing_range = get_res_nums_around(moved_with_side_chains_h, 
            center_resnum=out_res_num, 
            n_following=ccd_radius, 
            n_previous=ccd_radius, 
            include_intermediate=True)
        place_side_chains(moved_with_side_chains_h, original_pdb_h, 
            rotamer_manager, placing_range)
        # moved_with_side_chains_h.write_pdb_file(
        #     file_name="%s_after_sc_placement_%d.pdb" % (prefix, i))


        #
        # finalizing with geometry_minimization
        #
        # !!! This is the condition of acceptance of transformation!
        # if mc_rmsd < adaptive_mc_rmsd[ccd_radius]:
        if self.ccd_solution_is_ok(
            anchor_rmsd=resulting_rmsd, 
            mc_rmsd=mc_rmsd,
            ccd_radius=ccd_radius, 
            change_all_angles=change_all,
            change_radius=change_radius):
          if minimize:
            print >> self.log, "minimizing..."
            moved_with_side_chains_h.write_pdb_file(
                file_name="%s_result_before_min_%d.pdb" % (prefix, i))
            minimize_hierarchy(moved_with_side_chains_h, xrs, original_pdb_h, self.log)
          moved_with_side_chains_h.write_pdb_file(
              file_name="%s_result_minimized_%d.pdb" % (prefix, i))
          final_rmsd = get_main_chain_rmsd_range(moved_with_side_chains_h, 
              original_pdb_h, placing_range)
          print >> self.log, "FINAL RMSD after minimization:", final_rmsd
          return moved_with_side_chains_h
      ccd_radius += 1


    print >> self.log, "Epic FAIL: failed to fix rama outlier"
    all_results.sort(key=lambda tup: tup[1])
    print >> self.log, "  Options were: (mc_rmsd, resultign_rmsd, n_iter)"
    for i in all_results:
      print >> self.log, i[1:]
    # STOP()
    return original_pdb_h

  def get_resnums_of_chain_rama_outliers(self, pdb_hierarchy):
    phi_psi_atoms = utils.get_phi_psi_atoms(pdb_hierarchy)
    # print "len phi psi atoms", len(phi_psi_atoms)
    result = []
    rama_results = []
    ranges_for_idealization = []
    # print >> self.log, "rama outliers for input hierarchy:"
    list_of_reference_exclusion = []
    outp = utils.list_rama_outliers_h(pdb_hierarchy, self.r)
    print >> self.log, outp
    for phi_psi_pair, rama_key in phi_psi_atoms:
      # print "resseq:", phi_psi_pair[0][2].parent().parent().resseq
      ev = utils.rama_evaluate(phi_psi_pair, self.r, rama_key)
      # print "  ev", ev
      rama_results.append(ev)
      if ev == ramalyze.RAMALYZE_OUTLIER:
        resnum = phi_psi_pair[0][2].parent().parent().resseq
        result.append(resnum)
    return result


def get_main_chain_rmsd_range(
    hierarchy, original_h, all_atoms=False, placing_range=None):
  rmsd = 0
  mc_atoms = None
  if all_atoms:
    mc_atoms = ["N", "CA", "C", "O"]
  else:
    mc_atoms = ["N", "CA", "C"]
  for m_atom, ref_atom in zip(hierarchy.atoms(), original_h.atoms()):
    if m_atom.name.strip() in mc_atoms:
      if (placing_range is None or 
          m_atom.parent().parent().resseq in placing_range):
        rmsd += m_atom.distance(ref_atom)**2
  return rmsd**0.5


def minimize_hierarchy(hierarchy, xrs, original_pdb_h, 
    excl_string_selection, log=None):
  from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
  from mmtbx.refinement.geometry_minimization import run2
  from mmtbx.geometry_restraints import reference

  if log is None:
    log = null_out()
  params_line = grand_master_phil_str  
  params = iotbx.phil.parse(
      input_string=params_line, process_includes=True).extract()
  params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  params.pdb_interpretation.peptide_link.ramachandran_restraints = True
  params.pdb_interpretation.peptide_link.oldfield.weight_scale=3
  params.pdb_interpretation.peptide_link.oldfield.plot_cutoff=0.03
  params.pdb_interpretation.c_beta_restraints=True

  processed_pdb_files_srv = mmtbx.utils.\
      process_pdb_file_srv(
          crystal_symmetry= xrs.crystal_symmetry(),
          pdb_interpretation_params = params.pdb_interpretation,
          stop_for_unknowns         = False,
          log=log,
          cif_objects=None)
  processed_pdb_file, junk = processed_pdb_files_srv.\
      process_pdb_files(raw_records=flex.split_lines(hierarchy.as_pdb_string()))
  grm = get_geometry_restraints_manager(
      processed_pdb_file, xrs)

  asc = original_pdb_h.atom_selection_cache()
  sel = asc.selection(excl_string_selection)


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
    rotamer_manager, placing_range):
  ideal_res_dict = idealized_aa.residue_dict()
  asc = original_h.atom_selection_cache()
  gly_atom_names = set([" N  ", " CA ", " C  ", " O  "])
  for rg in hierarchy.residue_groups():
    if rg.resseq in placing_range:
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

def get_res_nums_around(pdb_hierarchy, center_resnum, n_following, n_previous,
    include_intermediate=False):
  residue_list = list(
      pdb_hierarchy.only_model().only_chain().only_conformer().residues())
  center_index = None
  for i in range(len(residue_list)):
    if residue_list[i].resseq == center_resnum:
      center_index = i
      break
  # print "start/end resids", residue_list[i-n_previous].resseq, residue_list[i+n_following].resseq
  # print "center i, len", center_index, len(residue_list)
  if not include_intermediate:
    return residue_list[max(0,center_index-n_previous)].resseq, \
        residue_list[min(len(residue_list)-1,center_index+n_following)].resseq
  else:
    res = []
    for i in range(max(0,center_index-n_previous), 
        min(len(residue_list)-1,center_index+n_following+1)):
      res.append(residue_list[i].resseq)
    return res


def get_fixed_moving_parts(pdb_hierarchy, out_res_num, n_following, n_previous):
  # limitation: only one  chain in pdb_hierarchy!!!
  original_pdb_h = pdb_hierarchy.deep_copy()
  start_res_num, end_res_num = get_res_nums_around(
      pdb_hierarchy, out_res_num, n_following, n_previous)

  xrs = original_pdb_h.extract_xray_structure()
  truncate_to_poly_gly(pdb_hierarchy, start_res_num, end_res_num)
  cache = pdb_hierarchy.atom_selection_cache()
  # print "selectioin:", "resid %d through %d" % (start_res_num, end_res_num)
  m_selection = cache.selection("resid %s through %s" % (start_res_num, end_res_num))
  moving_h = pdb_hierarchy.select(m_selection)
  moving_h.reset_atom_i_seqs()
  # print dir(moving_h)
  # STOP()
  m_cache = moving_h.atom_selection_cache()
  # print "len inp h atoms", pdb_hierarchy.atoms().size()
  # print "len moving_h atoms", moving_h.atoms().size()
  moving_ref_atoms_iseqs = []
  # here we need N, CA, C atoms from the end_res_num residue
  eff_end_resnum = end_res_num
  sel = m_cache.selection("resid %s" % end_res_num)
  while len(moving_h.select(sel).atoms()) == 0:
    eff_end_resnum -= 1
    sel = m_cache.selection("resid %s" % eff_end_resnum)

  sel = m_cache.selection("resid %s and name N" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_N = a.detached_copy()

  sel = m_cache.selection("resid %s and name CA" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_CA = a.detached_copy()

  sel = m_cache.selection("resid %s and name C" % eff_end_resnum)
  a = moving_h.select(sel).atoms()[0]
  moving_ref_atoms_iseqs.append(a.i_seq)
  fixed_C = a.detached_copy()

  fixed_ref_atoms = [fixed_N, fixed_CA, fixed_C]

  return moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms


def run(args, log=sys.stdout):
  # print "args", args

  inputs = mmtbx.utils.process_command_line_args(args=args,
      master_params=master_phil)
  work_params = inputs.params.extract()
  inputs.params.show(prefix=" ", out=log)
  pdb_file_names = list(inputs.pdb_file_names)
  if len(pdb_file_names) == 0:
    raise Sorry("No PDB file specified")
  work_params.loop_idealization.enabled=True
  # print work_params.loop_idealization.output_prefix
  # STOP()
  # work_params.model_idealization.file_name_before_regularization="before.pdb"
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  pdb_h = pdb_input.construct_hierarchy()


  loop_ideal = loop_idealization(pdb_h, work_params.loop_idealization, log)
  loop_ideal.resulting_pdb_h.write_pdb_file(
      file_name="%s_very_final.pdb" % work_params.loop_idealization.output_prefix)
  print >> log, "Outlier percentages: initial, after ccd, after minimization:"
  print >> log, loop_ideal.p_initial_rama_outliers, 
  print >> log, loop_ideal.p_before_minimization_rama_outliers,
  print >> log, loop_ideal.p_after_minimiaztion_rama_outliers



if (__name__ == "__main__"):
  run(sys.argv[1:])