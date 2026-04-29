from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.utils
from mmtbx.monomer_library import idealized_aa
from libtbx.utils import Sorry, null_out
from mmtbx.validation import ramalyze
from mmtbx.building.loop_closure.ccd import ccd_cpp
from mmtbx.building.loop_closure import utils, starting_conformations
from mmtbx.secondary_structure.build.ss_idealization import side_chain_placement, \
    set_xyz_smart
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran
from cctbx import maptbx
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from mmtbx.conformation_dependent_library import generate_protein_threes
import math
from libtbx import easy_pickle
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager
import mmtbx.refinement.real_space.fit_residues

#import boost_adaptbx.boost.python as bp
#ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
#from mmtbx_validation_ramachandran_ext import rama_eval

from iotbx.pdb.hybrid_36 import hy36encode, hy36decode


from mmtbx.refinement.real_space.individual_sites import minimize_wrapper_with_map
from six.moves import zip
from six.moves import range

loop_idealization_master_phil_str = """
loop_idealization
{
  enabled = True
    .type = bool
    .expert_level = 0
  output_prefix = None
    .type = str
    .expert_level = 0
  minimize_whole = True
    .type = bool
    .expert_level = 1
  fix_rotamers = True
    .type = bool
    .help = Fix rotamers before minimization
    .expert_level = 1
  force_rama_fixes = False
    .type = bool
    .help = If true, the procedure will pick and apply the best variant even \
      if all of them are above thresholds to be picked straight away. \
      Alternatively, when False, the procedure will accept failure and leave \
      a ramachandran outlier intact.
    .expert_level = 1
  include_allowed = False
    .type = bool
    .help = Try to fix residues in allowed Ramachandran region with CCD
  avoid_allowed_region = False
    .type = bool
    .help = If true, CCD will not let residues go to allowed region on \
      Ramachandran plot
  make_all_trans = True
    .type = bool
    .help = If true, procedure will try to get rid of all twisted and \
      cis-peptides making all of them trans.
    .expert_level = 1
  save_states = False
    .type = bool
    .help = Save states of CCD. Generates a states file for every model. \
      Warning! Significantly slower!
    .expert_level = 2
  number_of_ccd_trials = 3
    .type = int
    .help = How many times we are trying to fix outliers in the same chain
    .expert_level = 2
  variant_search_level = 2
    .type = int
    .help = how thoroughly variants will be explored (1-3)
    .expert_level = 2
  variant_number_cutoff = 50
    .type = int
    .help = how many first variants to take from generated
    .expert_level = 2
  variant_deviation_accept_level = low *med high very_high
    .type = choice(multi=False)
    .help = what deviation from starting model is acceptable when screening \
      variants. med is fine for real-space, low is for reciprocal space.
    expert_level = 2
  debug = False
    .type = bool
    .help = Output of intermediate files
    .expert_level = 3
}
"""

master_phil = iotbx.phil.parse(loop_idealization_master_phil_str)

#
# XXX
# XXX  This needs to be much more general tool to be called loop idealization.
# XXX  Right now it is just Ramachandran idealization. This should be
# XXX  separated. Loop idealization should include Cablam and Rama idealization.
# XXX
#


class loop_idealization():
  def __init__(self,
               model, # changed in place
               params=None,
               reference_map=None,
               log=null_out(),
               verbose=False,
               tried_rama_angles={},
               tried_final_rama_angles={},
               n_run=0):
    if model.get_number_of_models() > 1:
      raise Sorry("Multi-model files are not supported")
    self.model = model
    self.original_pdb_h = self.model.get_hierarchy()
    self.reference_map = reference_map
    self.params = self.process_params(params)
    self.log = log
    self.verbose = verbose
    self.ideal_res_dict = idealized_aa.residue_dict()
    self.n_run = n_run

    iaar_rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load(
      rotamers="favored")
    sin_cos_table = scitbx.math.sin_cos_table(n=10000)

    ram = ramalyze.ramalyze(pdb_hierarchy=self.model.get_hierarchy())
    self.p_initial_rama_outliers = ram.out_percent
    self.p_before_minimization_rama_outliers = None
    self.p_after_minimiaztion_rama_outliers = None

    # here we are recording what CCD solutions were used to fix particular
    # outliers to not use the same in the next CCD try.
    # Nested dict. First level:
    # key: chain id, value: dict
    #   key: resid (string), value: list of tried variants.
    self.tried_rama_angles = tried_rama_angles
    self.tried_final_rama_angles = tried_final_rama_angles

    berkeley_count = utils.list_rama_outliers_h(
        self.model.get_hierarchy(),
        self.model.get_ramachandran_manager()).count("\n")
    self.berkeley_p_before_minimization_rama_outliers = \
        berkeley_count/float(self.model.get_hierarchy().overall_counts().n_residues)*100
    n_bad_omegas = utils.n_bad_omegas(self.model.get_hierarchy())

    self.berkeley_p_after_minimiaztion_rama_outliers = self.berkeley_p_before_minimization_rama_outliers
    self.ref_exclusion_selection = ""
    self.number_of_ccd_trials = 0
    # print "logic expr outcome:", (self.number_of_ccd_trials < 10 and self.berkeley_p_before_minimization_rama_outliers > 0.001)
    # print self.number_of_ccd_trials < 10
    print("Rama outliers before idealization:", berkeley_count, self.berkeley_p_before_minimization_rama_outliers, file=self.log)
    if (self.berkeley_p_before_minimization_rama_outliers <= 0.001 and
        (n_bad_omegas<1 and self.params.make_all_trans)):
      print("No ramachandran outliers, skipping CCD step.", file=self.log)
    if self.verbose:
      print("n_bad_omegas", n_bad_omegas)
      print("self.params.make_all_trans",self.params.make_all_trans)
    if not self.params.enabled:
      print("Loop idealization is not enabled, use 'enabled=True'.", file=self.log)
    while (self.number_of_ccd_trials < self.params.number_of_ccd_trials
        and (self.berkeley_p_after_minimiaztion_rama_outliers > 0.001 or
            (n_bad_omegas>=1 and self.params.make_all_trans))
        and self.params.enabled):
      print("CCD try number, outliers:", self.number_of_ccd_trials, self.berkeley_p_before_minimization_rama_outliers, file=self.log)
      processed_chain_ids = []
      # TODO: verify tried_rama_angles is a dict or not
      for chain in self.model.get_master_hierarchy().only_model().chains():
        if chain.id not in self.tried_rama_angles:
          self.tried_rama_angles[chain.id] = {}
        if chain.id not in self.tried_final_rama_angles:
          self.tried_final_rama_angles[chain.id] = {}
        print("Idealizing chain %s" % chain.id, file=self.log)
        if chain.id not in processed_chain_ids:
          processed_chain_ids.append(chain.id)
        else:
          continue
        selection = "protein and chain %s and (name N or name CA or name C or name O)" % chain.id
        sel = self.model.selection("chain %s" % chain.id)
        chain_h = self.model.get_hierarchy().select(sel)
        m = chain_h.only_model()
        i = 0
        cutted_chain_h = None
        for c in m.chains():
          if i == 0:
            cutted_chain_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(c)
          else:
            print("WARNING!!! Duplicating chain ids! Only the first chain will be processed.", file=self.log)
            print("  Removing chain %s with %d residues" % (c.id, len(c.residues())), file=self.log)
            m.remove_chain(c)
          i += 1
        exclusions, ch_h = self.idealize_chain(
            hierarchy=(cutted_chain_h if cutted_chain_h else chain_h),
            include_allowed=self.params.include_allowed,
            tried_rama_angles_for_chain=self.tried_rama_angles[chain.id],
            tried_final_rama_angles_for_chain=self.tried_final_rama_angles[chain.id])
        if ch_h is not None:
          set_xyz_smart(
              # dest_h=self.model.get_hierarchy(),
              dest_h=chain,
              source_h=ch_h)
          self.model.set_sites_cart_from_hierarchy(multiply_ncs=True)
          for resnum in exclusions:
            selection += " and not resseq %s" % resnum
        self.ref_exclusion_selection += "(%s) or " % selection
        if self.verbose:
          print("self.tried_rama_angles", self.tried_rama_angles)
          print("self.tried_final_rama_angles", self.tried_final_rama_angles)
      #
      # dumping and reloading hierarchy to do proper rounding of coordinates
      resulting_pdb_h = iotbx.pdb.input(
          source_info=None,
          lines=self.model.get_hierarchy().as_pdb_or_mmcif_string()).construct_hierarchy()
      self.model.set_sites_cart(resulting_pdb_h.atoms().extract_xyz())

      berkeley_count = utils.list_rama_outliers_h(resulting_pdb_h).count("\n")
      self.berkeley_p_before_minimization_rama_outliers = \
          berkeley_count/float(resulting_pdb_h.overall_counts().n_residues)*100
      ram = ramalyze.ramalyze(pdb_hierarchy=resulting_pdb_h)
      self.p_before_minimization_rama_outliers = ram.out_percent
      if len(self.ref_exclusion_selection) > 0:
        self.ref_exclusion_selection = self.ref_exclusion_selection[:-3]

      duke_count = ram.get_outliers_count_and_fraction()[0]
      if self.params.debug or berkeley_count != duke_count:
        if berkeley_count != duke_count:
          print("Discrepancy between berkeley and duke after ccd:", berkeley_count, duke_count, file=self.log)
        self.model.pdb_or_mmcif_string_info(
          target_filename="%d%s_discrepancy.pdb" % (self.number_of_ccd_trials, self.params.output_prefix),
          write_file=True)

      if self.verbose:
        print("self.params.minimize_whole", self.params.minimize_whole)
      if self.params.minimize_whole:
        print("minimizing whole chain...", file=self.log)
        if self.verbose:
          print("self.ref_exclusion_selection", self.ref_exclusion_selection, file=self.log)
        # print >> sel
        # XXX but first let's check and fix rotamers...
        if self.params.fix_rotamers:
          print("Fixing/checking rotamers in loop idealization...", file=self.log)
          excl_sel = self.ref_exclusion_selection
          if len(excl_sel) == 0:
            excl_sel = None
          non_outliers_for_check = self.model.selection("(%s)" % self.ref_exclusion_selection)
          result = mmtbx.refinement.real_space.fit_residues.run(
              pdb_hierarchy     = self.model.get_hierarchy(),
              crystal_symmetry  = self.model.crystal_symmetry(),
              map_data          = self.reference_map,
              rotamer_manager   = iaar_rotamer_manager,
              sin_cos_table     = sin_cos_table,
              backbone_sample   = False,
              mon_lib_srv       = self.model.get_mon_lib_srv(),
              log               = self.log)
          model.set_sites_cart(
              sites_cart = result.pdb_hierarchy.atoms().extract_xyz())
        if self.params.debug:
          self.model.pdb_or_mmcif_string_info(
            target_filename="%d%s_all_not_minized.pdb" % (self.number_of_ccd_trials,self.params.output_prefix),
            write_file=True)
        if self.reference_map is None:
          minimize_wrapper_for_ramachandran(
              model = self.model,
              original_pdb_h=self.original_pdb_h,
              excl_string_selection=self.ref_exclusion_selection,
              log=None)
        else:
          mwwm = minimize_wrapper_with_map(
              model = self.model,
              target_map=self.reference_map,
              number_of_cycles=1,
              log=self.log)
      if self.params.debug:
        self.model.pdb_or_mmcif_string_info(
          target_filename="%d%s_all_minized.pdb" % (self.number_of_ccd_trials,self.params.output_prefix),
          write_file=True)
      ram = ramalyze.ramalyze(pdb_hierarchy=self.model.get_hierarchy())
      self.p_after_minimiaztion_rama_outliers = ram.out_percent
      berkeley_count = utils.list_rama_outliers_h(
          self.model.get_hierarchy(),
          self.model.get_ramachandran_manager()).count("\n")
      duke_count = ram.get_outliers_count_and_fraction()[0]
      n_bad_omegas = utils.n_bad_omegas(self.model.get_hierarchy())
      self.berkeley_p_after_minimiaztion_rama_outliers = \
          berkeley_count/float(self.model.get_hierarchy().overall_counts().n_residues)*100
      if berkeley_count != duke_count:
        print("Discrepancy between berkeley and duke after min:", berkeley_count, duke_count, file=self.log)
      else:
        print("Number of Rama outliers after min:", berkeley_count, file=self.log)
      print("Number of bad omegas:", n_bad_omegas, file=self.log)
      self.number_of_ccd_trials += 1
    # return new_h
    # return self.tried_rama_angles, self.tried_final_rama_angles

  def process_params(self, params):
    if params is None:
      params = master_phil.fetch().extract()
      params.loop_idealization.enabled = True
    if hasattr(params, "loop_idealization"):
      p_pars = params.loop_idealization
    else:
      assert hasattr(params, "enabled"), \
          "Something wrong with parameters passed to model_idealization"
      p_pars = params
    if p_pars.output_prefix is None:
      p_pars.output_prefix = "rama_fixed"
    assert isinstance(p_pars.enabled, bool)
    return p_pars

  def idealize_chain(self, hierarchy, include_allowed=False, tried_rama_angles_for_chain={},
      tried_final_rama_angles_for_chain={}):
    # check no ac:
    for c in hierarchy.chains():
      if len(c.conformers()) > 1:
        raise Sorry("Alternative conformations are not supported.")
      if "UNK" in c.get_residue_names_padded():
        pass
        # raise Sorry("UNK residues are not supported.")
    working_h = hierarchy.deep_copy()
    working_h.reset_atom_i_seqs()
    rama_results = []
    ranges_for_idealization = []
    print("rama outliers for input hierarchy:", file=self.log)
    rama_out_resnums = self.get_resnums_of_chain_rama_outliers(
        working_h,
        include_allowed=include_allowed)
    if len(rama_out_resnums) == 0:
      return None, None
    # get list of residue numbers that should be excluded from reference
    list_of_reference_exclusion = []
    for resnum in rama_out_resnums:
      excl_res = get_res_nums_around(
          hierarchy, [resnum], 2, 2, include_intermediate=True)
      list_of_reference_exclusion += excl_res
    out_i = 0
    chain_ss_annot = self.model.get_ss_annotation()
    if chain_ss_annot is not None:
      chain_ss_annot = chain_ss_annot.deep_copy()
      chain_ss_annot.remove_empty_annotations(hierarchy=working_h)
    # combine outliers next to each other
    comb_rama_out_resnums = [[rama_out_resnums[0]]]
    for r_out_resnum in rama_out_resnums[1:]:
      if abs(hy36decode(4, r_out_resnum)-hy36decode(4, comb_rama_out_resnums[-1][-1])) < 2:
        # combine
        comb_rama_out_resnums[-1].append(r_out_resnum)
      else:
        # separate
        comb_rama_out_resnums.append([r_out_resnum])

    print("Combined outliers for fixing:", comb_rama_out_resnums, file=self.log)
    for rama_out_resnum in comb_rama_out_resnums:
      print("Fixing outlier:", rama_out_resnum, file=self.log)
      self.log.flush()

      new_h = self.fix_rama_outlier(
        pdb_hierarchy=working_h,
        out_res_num_list=rama_out_resnum,
        prefix=self.params.output_prefix,
        minimize=False,
        ss_annotation=chain_ss_annot,
        tried_rama_angles_for_chain=tried_rama_angles_for_chain,
        tried_final_rama_angles_for_chain=tried_final_rama_angles_for_chain)
      if self.verbose:
        print("listing outliers after loop minimization", file=self.log)
        outp = utils.list_rama_outliers_h(new_h, self.model.get_ramachandran_manager())
        print(outp, file=self.log)
        utils.list_omega_outliers(new_h, self.log)
      self.log.flush()
      working_h = new_h
      out_i += 1
    return list_of_reference_exclusion, new_h

  def rangle_decart_dist(self, angle1, angle2):
    def normalize_angle(angle):
      result = (angle[0], angle[1])
      if angle[0] < 0:
        result = (angle[0]+180, angle[1])
      if angle[1] < 0:
        result = (angle[0], angle[1]+180)
      return result
    a1 = normalize_angle(angle1)
    a2 = normalize_angle(angle2)
    result = math.sqrt((a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 )
    return result

  def ccd_solution_is_duplicated(self,
      final_angles,
      tried_final_rama_angles_for_chain):
    # first checking for repeated solutions:
    # TODO verify these arguments final_fnales and tried_final* are dicts
    for rn, angles in final_angles:
      if rn in tried_final_rama_angles_for_chain:
        for previous_angles in tried_final_rama_angles_for_chain[rn]:
          if self.rangle_decart_dist(angles, previous_angles) < 20:
            print("Rejecting the same solution:", angles, previous_angles)
            return True
    return False

  def ccd_solution_is_ok(self,
      anchor_rmsd, mc_rmsd, n_outliers, ccd_radius,
      change_all_angles, change_radius,
      contains_ss_element, fixing_omega):

    # then checking rmsd
    adaptive_mc_rmsd = {1:3.0, 2:3.5, 3:4.0, 4:4.5, 5:5.5, 6:7.0, 7:8.5, 8:10.0, 9:12.0}
    num_of_run = max(self.number_of_ccd_trials, self.n_run)
    accept_level_coeff = 1
    if self.params.variant_deviation_accept_level == "low":
      accept_level_coeff = 0.5
    elif self.params.variant_deviation_accept_level == "high":
      accept_level_coeff = 1.5
    elif self.params.variant_deviation_accept_level == "very_high":
      accept_level_coeff = 2.5
    for k in adaptive_mc_rmsd:
      adaptive_mc_rmsd[k] = adaptive_mc_rmsd[k] * (1 + 0.3*num_of_run) * accept_level_coeff
    if fixing_omega:
      for k in adaptive_mc_rmsd:
        adaptive_mc_rmsd[k] += 1
    # print "adaptive_mc_rmsd", adaptive_mc_rmsd
    # adaptive_mc_rmsd = {1:2.5, 2:3.0, 3:3.5}
    if anchor_rmsd is None:
      anchor_rmsd = 0
    ss_multiplier = 1
    if contains_ss_element:
      ss_multiplier = 0.4
    # print "checking is_ok, n_out, target rmsd:", n_outliers, adaptive_mc_rmsd[ccd_radius+n_outliers-1]
    target_rmsd = adaptive_mc_rmsd.get(ccd_radius+n_outliers-1, 14.0)*ss_multiplier
    if (mc_rmsd < target_rmsd  and anchor_rmsd < 0.3):
      return True, target_rmsd
    elif ccd_radius == 3 and change_all_angles and change_radius == 2:
      # we are desperate and trying the most extensive search,
      # this deserves relaxed criteria...
      # print "desperate checking", adaptive_mc_rmsd[ccd_radius+n_outliers+1]*ss_multiplier
      target_rmsd = adaptive_mc_rmsd.get(ccd_radius+n_outliers+1, 14.0)*ss_multiplier
      return mc_rmsd < target_rmsd and anchor_rmsd < 0.4, target_rmsd
    return False, target_rmsd

  def fix_rama_outlier(self,
      pdb_hierarchy, out_res_num_list, prefix="", minimize=True,
      ss_annotation=None,
      tried_rama_angles_for_chain={},
      tried_final_rama_angles_for_chain={}):

    def comb_pair_in_bad_pairs(comb_pair, bad_pairs):
      if None in comb_pair:
        return False
      all_combs = [comb_pair]
      all_combs.append((comb_pair[0]-20, comb_pair[1]))
      all_combs.append((comb_pair[0]+20, comb_pair[1]))
      all_combs.append((comb_pair[0], comb_pair[1]-20))
      all_combs.append((comb_pair[0], comb_pair[1]+20))
      all_c_adj = []
      for p in all_combs:
        new_p = p
        if p[0] > 180:
          new_p = (p[0]-180, p[1])
        if p[0] < -180:
          new_p = (p[0]+180, p[1])
        if p[1] > 180:
          new_p = (p[0], p[1]-180)
        if p[0] < -180:
          new_p = (p[0], p[1]+180)
        all_c_adj.append(new_p)
      for p in all_c_adj:
        if p in bad_pairs:
          return True
      return False

    original_pdb_h = pdb_hierarchy.deep_copy()
    original_pdb_h.reset_atom_i_seqs()
    original_pdb_h_asc = original_pdb_h.atom_selection_cache()
    chain_id = original_pdb_h.only_model().only_chain().id
    all_results = []
    # only forward
    # variants_searches = [
    #     #ccd_radius, change_all, change_radius, direction_forward
    #     ((1, False, 0, True ),1),
    #     # ((1, False, 0, False),1),
    #     ((2, False, 0, True ),1),
    #     # ((2, False, 0, False),1),
    #     ((3, False, 0, True ),2),
    #     # ((3, False, 0, False),2),
    #     ((2, True,  1, True ),1),
    #     # ((2, True,  1, False),1),
    #     ((3, True,  1, True ),2),
    #     # ((3, True,  1, False),2),
    #     ((3, True,  2, True ),3),
    #     # ((3, True,  2, False),3),
    # ]
    # only backward
    # variants_searches = [
    #     #ccd_radius, change_all, change_radius, direction_forward
    #     # ((1, False, 0, True ),1),
    #     ((1, False, 0, False),1),
    #     # ((2, False, 0, True ),1),
    #     ((2, False, 0, False),1),
    #     # ((3, False, 0, True ),2),
    #     ((3, False, 0, False),2),
    #     # ((2, True,  1, True ),1),
    #     ((2, True,  1, False),1),
    #     # ((3, True,  1, True ),2),
    #     ((3, True,  1, False),2),
    #     # ((3, True,  2, True ),3),
    #     ((3, True,  2, False),3),
    # ]
    # both
    variants_searches = [
        #ccd_radius, change_all, change_radius, direction_forward
        ((1, False, 0, True ),1),
        ((1, False, 0, False),1),
        ((2, False, 0, True ),1),
        ((2, False, 0, False),1),
        ((3, False, 0, True ),2),
        ((3, False, 0, False),2),
        ((2, True,  1, True ),1),
        ((2, True,  1, False),1),
        ((3, True,  1, True ),2),
        ((3, True,  1, False),2),
        ((3, True,  2, True ),3),
        ((3, True,  2, False),3),
    ]
    decided_variants = []
    for variant, level in variants_searches:
      if level <= self.params.variant_search_level:
        decided_variants.append(variant)

    for ccd_radius, change_all, change_radius, direction_forward in decided_variants:
    # while ccd_radius <= 3:
      fixing_omega = False
      if self.verbose:
        print("  Starting optimization with radius=%d, " % ccd_radius, end=' ', file=self.log)
        print("change_all=%s, change_radius=%d, " % (change_all, change_radius), end=' ', file=self.log)
        print("direction=forward" if direction_forward else "direction=backwards", file=self.log)
        self.log.flush()
      #
      (moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms,
          m_selection, contains_ss_element, anchor_present) = get_fixed_moving_parts(
              pdb_hierarchy=pdb_hierarchy,
              out_res_num_list=out_res_num_list,
              # n_following=1,
              # n_previous=ccd_radius+ccd_radius-1,
              n_following=ccd_radius,
              n_previous=ccd_radius,
              ss_annotation=ss_annotation,
              direction_forward=direction_forward,
              log=self.log if self.verbose else None)
      # print "  moving_ref_atoms_iseqs", moving_ref_atoms_iseqs
      if self.verbose:
        print("  moving_h resseqs:", [x.resseq for x in moving_h.residue_groups()])
      moving_h_set = []
      all_angles_combination_f = starting_conformations.get_all_starting_conformations(
          moving_h,
          change_radius,
          include_allowed=self.params.include_allowed,
          n_outliers=len(out_res_num_list),
          direction_forward=direction_forward,
          cutoff=self.params.variant_number_cutoff,
          change_all=change_all,
          log=self.log if self.verbose else None,
          check_omega=self.params.make_all_trans,
          )

      #

      # print "len(all_angles_combination_f)", len(all_angles_combination_f)
      if len(all_angles_combination_f) == 0:
        print("In starting conformations - outlier was fixed?")
        # return result
      else:
        # here we should filter  first ones that in
        # tried_rama_angles_for_chain
        filter_out = [] # [[tried values],[tried values],...]
        for three in generate_protein_threes(
            hierarchy=moving_h,
            geometry=None):
          if three[1].resseq in list(tried_rama_angles_for_chain.keys()):
            filter_out.append(tried_rama_angles_for_chain[three[1].resseq])
          else:
            filter_out.append((None, None))
        ff_all_angles = []
        if self.verbose:
          print("filter_out", filter_out, file=self.log)
        for comb in all_angles_combination_f:
          good = True
          for comb_pair, bad_pairs in zip(comb, filter_out):
            if bad_pairs == (None, None):
              continue
            # print "comb_pair, bad_pairs", comb_pair, bad_pairs
            # if comb_pair in bad_pairs:
            if comb_pair_in_bad_pairs(comb_pair, bad_pairs):
              good = False
              # print "  Rejecting comb_pair", comb_pair
              break
          if good:
            ff_all_angles.append(comb)
        if self.verbose:
          print("len(all_angles_combination_f)", len(all_angles_combination_f), file=self.log)
          print("len(ff_all_angles)", len(ff_all_angles), file=self.log)
        n_added = 0
        n_all_combination = len(ff_all_angles)
        if n_all_combination == 0:
          print("Strange - got 0 combinations.", file=self.log)
        i_max = min(self.params.variant_number_cutoff, n_all_combination)
        # assert i_max > 0
        step = 0
        if i_max > 1:
          step = float(n_all_combination-1)/float(i_max-1)
        if step < 1:
          step = 1
        for i in range(i_max):
          comb = ff_all_angles[int(round(step*i))]
          setted_h, fixed_omega = starting_conformations.set_rama_angles(
                  moving_h,
                  list(comb),
                  direction_forward=direction_forward,
                  check_omega=self.params.make_all_trans)
          fixing_omega = fixing_omega or fixed_omega
          # print >> self.log, "Model %d, angles:" % i, comb
          if self.params.make_all_trans and utils.n_bad_omegas(setted_h) != 0:
            # Skipping conformation where omega was set incorrectly for
            # some reason.
            print("Model_%d_angles_%s.pdb" % (i, comb), end=' ')
            print("got ", utils.n_bad_omegas(moving_h_set[-1]), "bad omegas")
            moving_h_set[-1].write_pdb_or_mmcif_file(target_filename="Model_%d_angles_%s.pdb" % (i, comb))
            utils.list_omega(moving_h_set[-1], self.log)
            continue
            # assert 0
          moving_h_set.append(setted_h)

      if len(moving_h_set) == 0:
        # outlier was fixed before somehow...
        # or there's a bug in get_starting_conformations
        print("outlier was fixed before somehow", file=self.log)
        return original_pdb_h
      if self.verbose:
        print("self.tried_rama_angles inside", self.tried_rama_angles, file=self.log)
        print("tried_rama_angles_for_chain", tried_rama_angles_for_chain, file=self.log)
        print("checking values", ccd_radius, change_all, change_radius, direction_forward, file=self.log)
      for i, h in enumerate(moving_h_set):
        # if [x in tried_rama_angles_for_chain.keys() for x in out_res_num_list].count(True) > 0:
        #   print >> self.log, "Warning!!! make something here (check angles or so)"
        #   print >> self.log, "Skipping nonstable solution, tried previously:", (ccd_radius, change_all, change_radius, direction_forward, i)
        #   continue
        resulting_rmsd = None
        n_iter = 0
        if anchor_present:
          fixed_ref_atoms_coors = [x.xyz for x in fixed_ref_atoms]
          # print "params to constructor", fixed_ref_atoms, h, moving_ref_atoms_iseqs
          # easy_pickle.dump(file_name="crash.pkl", obj=[
          #     fixed_ref_atoms_coors,
          #     h,
          #     moving_ref_atoms_iseqs,
          #     direction_forward,
          #     self.params.save_states])
          ccd_obj = ccd_cpp(fixed_ref_atoms_coors, h, moving_ref_atoms_iseqs)
          ccd_obj.run(
              direction_forward=direction_forward,
              save_states=self.params.save_states,
              avoid_allowed_region=self.params.avoid_allowed_region)
          resulting_rmsd = ccd_obj.resulting_rmsd
          n_iter = ccd_obj.n_iter

          if self.params.save_states:
            states = ccd_obj.states
            states.write(file_name="%s%s_%d_%s_%d_%i_states.pdb" % (chain_id, out_res_num_list[0], ccd_radius, change_all, change_radius, i)) # PDB OK
        map_target = 0
        if self.reference_map is not None:
          map_target = maptbx.real_space_target_simple(
              unit_cell   = self.model.crystal_symmetry().unit_cell(),
              density_map = self.reference_map,
              sites_cart  = h.atoms().extract_xyz())

        mc_rmsd = get_main_chain_rmsd_range(moving_h, h, all_atoms=True)
        if self.verbose:
          print("Resulting anchor and backbone RMSDs, mapcc, n_iter for model %d, OK_rmsd:" % i, end=' ', file=self.log)
          print(resulting_rmsd, ",", mc_rmsd, ",", map_target, ",", n_iter, end=' ', file=self.log)
          self.log.flush()
        #
        # setting new coordinates
        #
        moved_with_side_chains_h = pdb_hierarchy.deep_copy()

        # setting xyz
        #
        for i_source, i_dest in enumerate(m_selection):
          moved_with_side_chains_h.atoms()[i_dest].set_xyz(h.atoms()[i_source].xyz)

        # set_xyz_smart(
        #     dest_h=moved_with_side_chains_h,
        #     source_h=h)

        #
        # placing side-chains
        #
        # moved_with_side_chains_h.write_pdb_file(
        #     file_name="%s_before_sc_placement_%d.pdb" % (prefix, i))
        placing_range = get_res_nums_around(moved_with_side_chains_h,
            center_resnum_list=out_res_num_list,
            n_following=ccd_radius,
            n_previous=ccd_radius,
            include_intermediate=True,
            avoid_ss_annot=ss_annotation)
        place_side_chains(moved_with_side_chains_h, original_pdb_h, original_pdb_h_asc,
            self.model.get_rotamer_manager(), placing_range, self.ideal_res_dict)
        # moved_with_side_chains_h.write_pdb_file(
        #     file_name="%s_after_sc_placement_%d.pdb" % (prefix, i))


        #
        # finalizing with geometry_minimization
        #

        # determining angles of interest
        # print "Recording picked angle for outliers"
        threes = generate_protein_threes(
          # hierarchy=moving_h,
          hierarchy=h,
          geometry=None)
        start_angles = []
        final_angles = []
        for angle_pair, three in zip(ff_all_angles[int(round(step*i))], threes):
          # print "three[1].resseq in out_res_num_list, angle_pair", three[1].resseq, out_res_num_list, angle_pair
          if three[1].resseq in out_res_num_list:
            # if three[1].resseq not in tried_rama_angles_for_chain.keys():
            #   tried_rama_angles_for_chain[three[1].resseq] = []
            start_angles.append((three[1].resseq, angle_pair))
            ps_angles = three.get_phi_psi_angles()
            final_angles.append((three[1].resseq, tuple(ps_angles)))
            # tried_rama_angles_for_chain[three[1].resseq].append(angle_pair)
            # print >> self.log, "Ended up with", three[1].resseq, "%.1f %.1f" % (ps_angles[0], ps_angles[1])
        # print "Updated tried_rama_angles_for_chain:", tried_rama_angles_for_chain
        if (not self.ccd_solution_is_duplicated(
            final_angles=final_angles,
            tried_final_rama_angles_for_chain=tried_final_rama_angles_for_chain)):
          all_results.append((moved_with_side_chains_h.deep_copy(), mc_rmsd, resulting_rmsd, map_target, n_iter))
        else:
          if self.verbose:
            print("Duplicate.", file=self.log)
          continue
        (ccd_ok, target_rmsd) = self.ccd_solution_is_ok(
            anchor_rmsd=resulting_rmsd,
            mc_rmsd=mc_rmsd,
            n_outliers=len(out_res_num_list),
            ccd_radius=ccd_radius,
            change_all_angles=change_all,
            change_radius=change_radius,
            contains_ss_element=contains_ss_element,
            fixing_omega=fixing_omega)
        if self.verbose:
          print(", %.2f" % target_rmsd, file=self.log)
        if ccd_ok:
          if self.verbose:
            print("Choosen result (mc_rmsd, anchor_rmsd, map_target, n_iter):", mc_rmsd, resulting_rmsd, map_target, n_iter, file=self.log)
          # Save to tried_ccds
          for rn, angles in start_angles:
            if rn not in list(tried_rama_angles_for_chain.keys()):
              tried_rama_angles_for_chain[rn] = []
            tried_rama_angles_for_chain[rn].append(angles)
          # Save final angles
          for rn, angles in final_angles:
            if rn not in list(tried_final_rama_angles_for_chain.keys()):
              tried_final_rama_angles_for_chain[rn] = []
            tried_final_rama_angles_for_chain[rn].append(angles)
          if self.verbose:
            print("Ended up with", final_angles, file=self.log)
            print("Updated tried_rama_angles_for_chain:", tried_rama_angles_for_chain, file=self.log)
            print("Updated tried_final_rama_angles_for_chain:", tried_final_rama_angles_for_chain, file=self.log)

          self.log.flush()
          assert not minimize
          # This is not working....
          if minimize:
            print("minimizing...", file=self.log)
            # moved_with_side_chains_h.write_pdb_file(
            #     file_name="%s_result_before_min_%d.pdb" % (prefix, i))
            if self.reference_map is None:
              minimize_wrapper_for_ramachandran(
                  hierarchy=moved_with_side_chains_h,
                  xrs=xrs,
                  original_pdb_h=original_pdb_h,
                  log=self.log,
                  grm=self.grm,
                  ss_annotation=self.secondary_structure_annotation)
            else:
              mwwm = minimize_wrapper_with_map(
                  pdb_h=moved_with_side_chains_h,
                  xrs=xrs,
                  target_map=self.reference_map,
                  grm=self.grm,
                  ss_annotation=self.secondary_structure_annotation,
                  log=self.log)
          # moved_with_side_chains_h.write_pdb_file(
          #     file_name="%s_result_minimized_%d.pdb" % (prefix, i))
          final_rmsd = get_main_chain_rmsd_range(moved_with_side_chains_h,
              original_pdb_h, placing_range)
          if self.verbose:
            print("FINAL RMSD after minimization:", final_rmsd, file=self.log)
          return moved_with_side_chains_h


    all_results.sort(key=lambda tup: tup[1])
    if self.verbose:
      print("ALL RESULTS:", file=self.log)
      i = 0
      for ar in all_results:
        print(ar[1:], end=' ', file=self.log)
        if ar[2] < 0.4:
          # fn = "variant_%d.pdb" % i
          # ar[0].write_pdb_file(file_name=fn)
          # print fn
          i += 1
        else:
          print("  no output", file=self.log)
    if self.params.force_rama_fixes:
      # find and apply the best varian from all_results. This would be the one
      # with the smallest rmsd given satisfactory closure
      print("Applying the best found variant:", end=' ', file=self.log)
      i = 0
      while i < len(all_results) and all_results[i][2] > 1.5:
        i += 1
      # apply
      # === duplication!!!!
      if i < len(all_results):
        print(all_results[i][1:], file=self.log)
        if minimize:
          print("minimizing...", file=self.log)
          # all_results[i][0].write_pdb_file(
          #     file_name="%s_result_before_min_%d.pdb" % (prefix, i))
          if self.reference_map is None:
            minimize_wrapper_for_ramachandran(
                hierarchy=all_results[i][0],
                xrs=xrs,
                original_pdb_h=original_pdb_h,
                log=self.log,
                grm=self.grm,
                ss_annotation=self.secondary_structure_annotation)
          else:
            mwwm = minimize_wrapper_with_map(
                pdb_h=all_results[i][0],
                xrs=xrs,
                target_map=self.reference_map,
                grm=self.grm,
                ss_annotation=self.secondary_structure_annotation,
                log=self.log)
        # all_results[i][0].write_pdb_file(
        #     file_name="%s_result_minimized_%d.pdb" % (prefix, i))
        final_rmsd = get_main_chain_rmsd_range(all_results[i][0],
            original_pdb_h, placing_range)
        if self.verbose:
          print("FINAL RMSD after minimization:", final_rmsd, file=self.log)
        return all_results[i][0]
      else:
        if self.verbose:
          print(" NOT FOUND!", file=self.log)
          for i in all_results:
            print(i[1:], file=self.log)
      # === end of duplication!!!!

    else:
      if self.verbose:
        print("Epic FAIL: failed to fix rama outlier:", out_res_num_list, file=self.log)
        print("  Options were: (mc_rmsd, resultign_rmsd, n_iter)", file=self.log)
        for i in all_results:
          print(i[1:], file=self.log)
    return original_pdb_h

  def get_resnums_of_chain_rama_outliers(self, pdb_hierarchy, include_allowed=False):
    phi_psi_atoms = utils.get_phi_psi_atoms(pdb_hierarchy, omega=True)
    # pdb_hierarchy.write_pdb_file(file_name="phi_psi_atoms.pdb")
    # print "len phi psi atoms", len(phi_psi_atoms)
    result = []
    rama_results = []
    ranges_for_idealization = []
    # print >> self.log, "rama outliers for input hierarchy:"
    list_of_reference_exclusion = []
    outp = utils.list_rama_outliers_h(pdb_hierarchy, self.model.get_ramachandran_manager(), include_allowed=include_allowed)
    print(outp, file=self.log)
    for phi_psi_pair, rama_key, omega in phi_psi_atoms:
      if phi_psi_pair[0] is not None and phi_psi_pair[1] is not None:
        # print "resseq:", phi_psi_pair[0][2].parent().parent().resseq
        ev = utils.rama_evaluate(phi_psi_pair, self.model.get_ramachandran_manager(), rama_key)
        # print "  ev", ev
        rama_results.append(ev)
        # print phi_psi_pair[0][0].id_str()
        resnum = phi_psi_pair[0][2].parent().parent().resseq
        if ev == ramalyze.RAMALYZE_OUTLIER:
          result.append(resnum)
        if include_allowed and ev == ramalyze.RAMALYZE_ALLOWED:
          result.append(resnum)
        if omega is not None and abs(abs(omega)-180) > 30:
          print("Spotted twisted/cis peptide:", resnum, omega, file=self.log)
          result.append(resnum)
    # STOP()
    return result


def place_side_chains(hierarchy, original_h, original_h_asc,
    rotamer_manager, placing_range, ideal_res_dict):
  # ideal_res_dict = idealized_aa.residue_dict()
  # asc = original_h.atom_selection_cache()
  gly_atom_names = set([" N  ", " CA ", " C  ", " O  "])
  n_ca_c_present = [False, False, False]
  for rg in hierarchy.residue_groups():
    if rg.resseq in placing_range:
      # cut extra atoms
      ag = rg.only_atom_group()
      for atom in ag.atoms():
        if (atom.name not in gly_atom_names):
          ag.remove_atom(atom=atom)
        if atom.name == " N  ":
          n_ca_c_present[0] = True
        elif atom.name == " CA ":
          n_ca_c_present[1] = True
        elif atom.name == " C  ":
          n_ca_c_present[2] = True
      if n_ca_c_present.count(True) < 3:
        # Not all essential atoms present for placing side-chain.
        continue
      # get ag from original hierarchy
      orig_ag = original_h.select(original_h_asc.selection("resseq %s" % rg.resseq)
          ).models()[0].chains()[0].residue_groups()[0].atom_groups()[0]
      # get ideal
      # ideal_ag = ideal_res_dict[ag.resname.lower()].models()[0].chains()[0].\
      #   residue_groups()[0].atom_groups()[0]
      # print "got to placement"
      side_chain_placement(ag, orig_ag, rotamer_manager)

def get_loop_borders(pdb_hierarchy, center_resnum_list, ss_annot):
  """ get loop resum beginning and end around center_resnum"""
  f_start_res_num =-9999
  f_end_res_num = 9999999
  if ss_annot is not None:
    for elem in ss_annot.simple_elements():
      if f_start_res_num < elem.get_end_resseq_as_int() <= hy36decode(4, center_resnum_list[0]):
        # print "  cutting..."
        f_start_res_num = elem.get_end_resseq_as_int()
      if hy36decode(4, center_resnum_list[-1]) <= elem.get_start_resseq_as_int() < f_end_res_num:
        # print "  cutting..."
        f_end_res_num = elem.get_start_resseq_as_int()
  loop_length = f_end_res_num - f_start_res_num
  return f_start_res_num, f_end_res_num


def get_res_nums_around(pdb_hierarchy, center_resnum_list, n_following, n_previous,
    include_intermediate=False, avoid_ss_annot=None, log=None):
  """
  Warning, this function most likely won't work properly with insertion codes
  """
  if log is None:
    log = StringIO()
  working_ss_annot = None
  if avoid_ss_annot is not None:
    working_ss_annot = avoid_ss_annot.deep_copy()
    working_ss_annot.remove_empty_annotations(
        hierarchy=pdb_hierarchy)
  residue_list = list(
      pdb_hierarchy.only_model().only_chain().only_conformer().residues())
  center_index = []
  for i in range(len(residue_list)):
    if residue_list[i].resseq in center_resnum_list:
      center_index.append(i)
      # break
  if not include_intermediate:
    # return residue_list[max(0,center_index-n_previous)].resseq, \
    #     residue_list[min(len(residue_list)-1,center_index+n_following)].resseq
    print("center_index, resnum list", center_index, center_resnum_list, file=log)
    # assert len(center_index) == len(center_resnum_list)
    start_res_num = residue_list[max(0,center_index[0]-n_previous)].resseq_as_int()
    end_res_num = residue_list[min(len(residue_list)-1,center_index[-1]+n_following)].resseq_as_int()
    srn, ern = get_loop_borders(pdb_hierarchy, center_resnum_list, working_ss_annot)
    print("start_res_num, end_res_num", start_res_num, end_res_num, file=log)
    print("srn, ern", srn, ern, file=log)
    # srn, ern = -9999, 9999999
    # So now we have borders of the loop: srn, ern, center_resnum,
    # n_following, n_previous.
    # We combine the above knowledge to adjust the borders keeping the same
    # loop size
    # adjst beginning
    if srn > start_res_num:
      end_res_num += srn - start_res_num
      start_res_num = srn
    # adjust end
    if ern < end_res_num:
      start_res_num = max(start_res_num - (end_res_num-ern), srn)
      end_res_num = ern

    f_start_res_num = start_res_num
    f_end_res_num = end_res_num
    print("srn, ern", srn, ern, file=log)
    print("f_start_res_num, f_end_res_num", f_start_res_num, f_end_res_num, file=log)
    if f_end_res_num == hy36decode(4, center_resnum_list[-1]):
      f_end_res_num += 1
    if f_start_res_num == hy36decode(4,center_resnum_list[0]):
      f_start_res_num -= 1
    print("after f_start_res_num, f_end_res_num", f_start_res_num, f_end_res_num, file=log)
    return hy36encode(4, f_start_res_num), hy36encode(4, f_end_res_num)
  else:
    res = []
    for i in range(max(0,center_index[0]-n_previous),
        min(len(residue_list)-1,center_index[-1]+n_following+1)):
      res.append(residue_list[i].resseq)
    return res

def get_fixed_moving_parts(pdb_hierarchy, out_res_num_list, n_following, n_previous,
    ss_annotation=None, direction_forward=True, log=None):
  # limitation: only one  chain in pdb_hierarchy!!!
  if log is None:
    log = StringIO()
  original_pdb_h = pdb_hierarchy.deep_copy()
  # print >> log, "  out_res_num, n_following, n_previous", out_res_num_list, n_following, n_previous
  start_res_num, end_res_num = get_res_nums_around(pdb_hierarchy, out_res_num_list,
      n_following, n_previous, include_intermediate=False, avoid_ss_annot=ss_annotation, log=log)
  print("  start_res_num, end_res_num", start_res_num, end_res_num, file=log)
  xrs = original_pdb_h.extract_xray_structure()
  pdb_hierarchy.truncate_to_poly_gly()
  cache = pdb_hierarchy.atom_selection_cache()
  # print "POSSIBLE ERROR:", "selectioin:", "(name N or name CA or name C or name O) and resid %s through %s" % (
  #         start_res_num, end_res_num)
  m_selection = cache.iselection(
      "(name N or name CA or name C or name O) and resid %s:%s" % (
          start_res_num, end_res_num))
  # Somewhere here would be the place to tweak n_following, n_previous to
  # exclude SS parts. It would be nice to increase n_prev in case
  # we need to cut on n_following etc.
  # If no ss_annotation is provided, don't filter.
  contains_ss_element = False
  if ss_annotation is not None:
    ss_selection_str = ss_annotation.overall_selection()
    ss_selection = cache.iselection(ss_selection_str)
    intersect = flex.size_t(sorted(list(set(ss_selection) & set(m_selection))))
    if intersect.size() > 0:
      intersect_h = pdb_hierarchy.select(intersect)
      print("Hitting SS element", file=log)
      print(intersect_h.as_pdb_or_mmcif_string(), file=log)
      contains_ss_element = False
      assert intersect_h.atoms_size() > 0, "Wrong atom count in SS intersection"
      # assert 0, "hitting SS element!"

  # print "m_selection=", list(m_selection)
  moving_h = pdb_hierarchy.select(m_selection)
  moving_h.reset_atom_i_seqs()
  # print >> log, "Moving h:"
  # print >> log, moving_h.as_pdb_string()
  # for rg in moving_h.only_chain().residue_groups():
  #   print >> log, "rg:", rg.resseq
  # print dir(moving_h)
  # STOP()
  m_cache = moving_h.atom_selection_cache()
  # print "len inp h atoms", pdb_hierarchy.atoms_size()
  # print "len moving_h atoms", moving_h.atoms_size()
  # here we need N, CA, C atoms from the end_res_num residue
  eff_end_resnum = end_res_num
  if not direction_forward:
    eff_end_resnum = start_res_num
  sel = m_cache.selection("resid %s" % end_res_num)
  int_eff_resnum = hy36decode(4,eff_end_resnum)
  while len(moving_h.select(sel).atoms()) == 0:
    if direction_forward:
      int_eff_resnum -= 1
    else:
      int_eff_resnum += 1
    # print >> log, "resid %d" % int_eff_resnum
    sel = m_cache.selection("resid %d" % int_eff_resnum)
    # print >> log, "sel:", list(sel)
    if int_eff_resnum < -10:
      assert 0
  eff_end_resnum = hy36encode(4, int_eff_resnum)

  anchor_present = True
  moving_ref_atoms_iseqs = []
  fixed_ref_atoms = []
  # print "fixed_ref_atoms:"
  base_sel = "resid %s and name " % eff_end_resnum
  for ssel in ["N", "CA", "C"]:
    sel = m_cache.selection(base_sel+ssel)
    atoms = moving_h.select(sel).atoms()
    if atoms.size() > 0:
      moving_ref_atoms_iseqs.append(atoms[0].i_seq)
      fixed_ref_atoms.append(atoms[0].detached_copy())
    else:
      anchor_present = False
    # print "  ", atoms[0].id_str()
  print("anchor_present", anchor_present, file=log)
  return (moving_h, moving_ref_atoms_iseqs, fixed_ref_atoms, m_selection,
      contains_ss_element, anchor_present)

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
