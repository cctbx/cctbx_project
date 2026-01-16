from __future__ import absolute_import, division, print_function

from libtbx import group_args
import libtbx.phil

from mmtbx.validation.cablam import cablamalyze, fetch_peptide_expectations, \
    fetch_ca_expectations, fetch_motif_contours
from libtbx.utils import Sorry, null_out

from cctbx import geometry_restraints

import itertools

from scitbx.matrix import rotate_point_around_axis

from mmtbx.refinement.geometry_minimization import run2
from mmtbx.building.loop_closure.utils import list_rama_outliers_h
from mmtbx.secondary_structure import manager as ss_manager_class
from mmtbx.geometry_restraints.ramachandran import master_phil as rama_master_phil
import six
from six.moves import range


master_phil_str = '''
cablam_idealization {
  enabled = True
    .type = bool
  nproc = 1
    .type = int
    .help = Parallelization is not implemented
  rotation_angle = 30
    .type = int
    .help = angle to rotate the oxygen while searching for better conformation
  require_h_bond = True
    .type = bool
    .help = require presence of (new) h-bond after rotation
  do_gm = False
    .type = bool
    .help = Run geometry minimization after rotation
  find_ss_after_fixes = True
    .type = bool
    .help = re-evaluate SS after fixing Cablam outliers. May be helpful to \
      identify new or extend previous SS elements
  save_intermediates = False
    .type = bool
    .help = Save all cablam rotation for particular residue in separate file
}
'''

# This is needed to import scope
master_phil = libtbx.phil.parse(master_phil_str)

class cablam_idealization(object):
  def __init__(self, model, params, log):
    """
    model is changed in place
    params - those in master_phil_str without scope name
    """
    self.model = model
    self.params = params
    self.log = log
    self.outliers = None
    self.cablam_contours = fetch_peptide_expectations()
    self.ca_contours = fetch_ca_expectations()
    self.motif_contours = fetch_motif_contours()
    self.n_tried_residues = 0
    self.n_rotated_residues = 0
    self.cablam_fixed_minimized = None
    self.cablam_results = {}

    if not self.params.enabled:
      return

    self.model.process(make_restraints=True)

    print("CaBLAM idealization", file=self.log)

    if self.model.get_hierarchy().models_size() > 1:
      raise Sorry("Multi-model files are not supported")

    self.model.search_for_ncs()
    ncs_obj = self.model.get_ncs_obj()
    nrgl = ncs_obj.get_ncs_restraints_group_list()
    if nrgl.get_n_groups() > 0:
      print(self.model.get_ncs_obj().show_phil_format(), file=self.log)

    self.outliers_by_chain = self.identify_outliers()

    # idealization
    # TODO: verify if outliers by chain is dict
    for chain, outliers in six.iteritems(self.outliers_by_chain):
      b_selection = self.model.selection("chain %s" % chain)
      self.atoms_around = self.model.get_xray_structure().selection_within(7, b_selection).iselection()
      self.cablam_results[chain] = []

      for outlier in outliers:
        rotated_angle = self.fix_cablam_outlier(chain, outlier)
        if rotated_angle is not None:
          self.cablam_results[chain].append((outlier[-1], rotated_angle))
    print("*"*80, file=self.log)

    if self.params.find_ss_after_fixes:
      ss_manager = ss_manager_class(
          pdb_hierarchy=self.model.get_hierarchy(),
          geometry_restraints_manager=self.model.get_restraints_manager().geometry,
          sec_str_from_pdb_file=None,
          params=None,
          mon_lib_srv=self.model.get_mon_lib_srv(),
          verbose=-1,
          log=self.log)
      self.model.get_restraints_manager().geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.model.get_hierarchy(),
          log=self.log)
      self.model.set_ss_annotation(ann=ss_manager.actual_sec_str)

    if params.do_gm:
      self.cablam_fixed_minimized = self._minimize()

  def _get_ca_atom(self, chainid, resid):
    for chain in self.model.get_master_hierarchy().only_model().chains():
      if chain.id.strip() == chainid.strip():
        for rg in chain.residue_groups():
          if rg.resid() == resid:
            for a in rg.atoms():
              if a.name.strip() == "CA":
                return a
    raise Sorry("Something went wrong. Cannot find CA atom.")
    return None

  @staticmethod
  def decode_feedback(feedback):
    if feedback.alpha:
      true_label = "alpha"
    elif feedback.beta:
      true_label = "beta"
    elif feedback.threeten:
      true_label = "threeten"
    else:
      true_label = "loop"
    return true_label

  def fix_cablam_outlier(self, chain, outlier):
    scores = []
    if len(outlier) < 3:
      curresid = outlier[-1].residue.resid()
      prevresid = outlier[-1].prevres.residue.resid()
      curresseq_int = outlier[-1].residue.resseq_as_int()
      prevresseq_int = outlier[-1].prevres.residue.resseq_as_int()
    else:
      print("Don't know how to deal with more than 2 outliers in a row yet. Skipping.", file=self.log)
      return
    # h =  self.model.get_hierarchy()
    # s =  self.model.selection("chain %s and name CA and resid %s" % (chain, prevresid))
    # a1 = self.model.select(s).get_hierarchy().atoms()[0]
    # s =  self.model.selection("chain %s and name CA and resid %s" % (chain, curresid))
    # a2 = self.model.select(s).get_hierarchy().atoms()[0]
    # This is slightly faster, but poorer code. We'll see if it is needed.
    a1 = self._get_ca_atom(chain, prevresid)
    a2 = self._get_ca_atom(chain, curresid)

    print("*"*80, file=self.log)
    print("Atoms for rotation:", chain, prevresid, curresid, file=self.log)
    print(f"Cablam evaluation: {self.decode_feedback(outlier[-1].feedback)}", file=self.log)
    print("*"*80, file=self.log)

    around_str_sel = "chain %s and resid %d:%d" % (chain, prevresseq_int-2, curresseq_int+2)
    chain_around = self.model.select(self.model.selection(around_str_sel))
    assert chain_around.get_number_of_atoms() > 0
    self.atoms_around_cutted = self.atoms_around.deep_copy()
    # angle = 30
    angle = self.params.rotation_angle

    for i in range(int(360/angle)):
      # rotation
      O_atom, N_atom, C_atom = self._rotate_cablam(self.model, chain,
          prevresid, curresid, a1, a2, angle=angle)
      if [O_atom, N_atom, C_atom].count(None) > 0:
        print("Residues are missing essential atom: O, N or C. Skipping.", file=self.log)
        return
      self._rotate_cablam(chain_around, chain,
          prevresid, curresid, a1, a2, angle=angle)
      if self.params.save_intermediates:
        self.model.pdb_or_mmcif_string_info(
            target_filename="out_%s_%d.pdb" % (curresid.strip(), i),
            write_file=True)
      # Score contains tuple of :
      # angle, N Rama outliers, N Cablam outliers, hbonds
      scores.append(self._score_conformation(O_atom, C_atom, N_atom, chain_around, angle*(i+1)))
    print("angle, rama outliers, cablam outliers, hbonds (type, length, angle)", file=self.log)
    for s in scores:
      print(s[0], s[1], s[2], end=' ', file=self.log)
      if len(s[3]) > 0:
        for e in s[3]:
          print("| %s -- > %s, %.2f, %.2f|" % (e[0], e[3].id_str(), e[1], e[2]), end=' ', file=self.log)
      print(file=self.log)
    rot_angle = self._pick_rotation_angle(scores, self.params.require_h_bond)
    # rotate
    self.n_tried_residues += 1
    if rot_angle != 360:
      self.n_rotated_residues += 1
      print("ROTATING by", rot_angle, file=self.log)
      self._rotate_cablam(self.model, chain,
          prevresid, curresid, a1, a2, angle=rot_angle)
    return rot_angle

  def _rotate_cablam(self, model, chain, prevresid, curresid, a1, a2, angle):
    inside = False
    O_atom = None
    N_atom = None
    C_atom = None
    for c in model.get_master_hierarchy().only_model().chains():
      if c.id.strip() == chain.strip():
        for atom in c.atoms():
          if atom.name.strip() == "CA" and atom.parent().parent().resid() == prevresid:
            inside = True
          if atom.name.strip() == "CA" and atom.parent().parent().resid() == curresid:
            inside = False
          if inside and atom.name.strip() in ["N", "CA", "C", "O"]:
            new_xyz = rotate_point_around_axis(
                axis_point_1=a1.xyz,
                axis_point_2=a2.xyz,
                point=atom.xyz,
                angle=angle,
                deg=True)
            atom.set_xyz(new_xyz)
            if atom.name.strip() == "O":
              O_atom = atom
            elif atom.name.strip() == "N":
              N_atom = atom
            elif atom.name.strip() == "C":
              C_atom = atom

        model.set_sites_cart_from_hierarchy(multiply_ncs=True)

        return O_atom, N_atom, C_atom


  def _pick_rotation_angle(self, scores, require_h_bond=True):
    # I want to pick the rotation with H-bond, less Rama outliers and less
    # cablam outliers.
    best = scores[-1] # last, no rotation
    for s in scores[:-1]:
      # [angle, rama, cablam, hbond]
      if require_h_bond and len(s[3]) == 0:
        continue
      if require_h_bond and len(s[3]) > 0:
        # rama better or same, and cablam outliers decreased:
        if s[1] <= best[1] and s[2] < best[2]:
          best = s
      else:
        # same as above, but preserve or increase number of hbonds
        if len(s[3]) >= len(best[3]) and s[1] <= best[1] and s[2] < best[2]:
          best = s
    # check for the claster and try to return the middle one if the case.
    best_position = scores.index(best)
    i = best_position
    while (i < len(scores) and
           scores[i][1] == scores[best_position][1] and
           scores[i][2] == scores[best_position][2] and
           len(scores[i][3]) == len(scores[best_position][3])):
      i += 1
    # If we went past the end, step back one position.
    # if i == len(scores):
    # Backtrack because the condition was broken:
    i -= 1
    middle_pos = (best_position + i) // 2
    # print('lenlenlen', len(scores))
    print('cluster positions:', best_position, i, middle_pos, file=self.log)
    print('cluster angles:', scores[best_position][0], scores[i][0], scores[middle_pos][0], file=self.log)
    return scores[middle_pos][0]

  def _minimize(self):
    m1 = self.model.deep_copy()
    rama_params = rama_master_phil.fetch().extract().ramachandran_plot_restraints
    rama_params.favored = 'oldfield'
    rama_params.allowed = 'oldfield'
    rama_params.outlier = 'oldfield'
    m1.set_ramachandran_plot_restraints(rama_params=rama_params)
    run2(
        restraints_manager=m1.get_restraints_manager(),
        pdb_hierarchy=m1.get_hierarchy(),
        correct_special_position_tolerance=1.0,
        riding_h_manager               = None,
        ncs_restraints_group_list      = [], # These are actually for NCS CONSTRAINTS!
        max_number_of_iterations       = 500,
        number_of_macro_cycles         = 5,
        selection                      = None,
        bond                           = True,
        nonbonded                      = True,
        angle                          = True,
        dihedral                       = True,
        chirality                      = True,
        planarity                      = True,
        parallelity                    = True,
        log = null_out())
    m1.set_sites_cart_from_hierarchy(multiply_ncs=True)
    return m1

  def _score_conformation(self, O_atom, C_atom, N_atom, chain_around, angle):
    # gs = self.model.geometry_statistics()
    # gs.show()
    # print "MOLPROBITY Score:", gs.result().molprobity_score
    # print "Cablam outliers:", gs.result().cablam.outliers
    # print "Clashscore: ", gs.result().clash.score
    hbonds = self._search_hbond(O_atom, C_atom, N_atom, chain_around)
    ro = list_rama_outliers_h(chain_around.get_hierarchy())
    cab_results = cablamalyze(
        pdb_hierarchy=chain_around.get_hierarchy(),
        outliers_only=True,
        out=null_out(),
        quiet=True,
        cablam_contours = self.cablam_contours,
        ca_contours = self.ca_contours,
        motif_contours = self.motif_contours,
        )
    outliers_only = [x for x in cab_results.results if x.feedback.cablam_outlier]
    return (angle, len(ro.split("\n"))-1, len(outliers_only), hbonds)

  def _search_hbond(self, O_atom, C_atom, N_atom, chain_around):
    def good_hbond(angle, distance):
      return angle > 140 and distance < 3.8
    results = []
    atoms = self.model.get_atoms()
    filtered_atoms_around_cutted = []
    for atom in [atoms[i_seq] for i_seq in self.atoms_around_cutted]:
      if atom.distance(O_atom) > 10:
        continue
      # no need to check the same residue, looking for N atom for bonding
      filtered_atoms_around_cutted.append(atom.i_seq)
      if atom.parent() == O_atom.parent() or atom.parent() == N_atom.parent():
        # print "skipping same residue ", atom.id_str()
        continue
      if atom.name.strip() == 'N':
        angle = geometry_restraints.angle(
            sites=[C_atom.xyz, O_atom.xyz, atom.xyz],
            angle_ideal=0,
            weight=1).angle_model
        if good_hbond(angle, atom.distance(O_atom)):
          # print "Potential bond:", atom.id_str(), atom.distance(O_atom), angle
          results.append((O_atom.id_str(), atom.distance(O_atom), angle, atom))
      if atom.name.strip() == 'O':
        # now we want to find attached N atom (another one)
        another_C_atom = atom.parent().get_atom("C")
        if another_C_atom is not None:
          angle = geometry_restraints.angle(
              sites=[another_C_atom.xyz, atom.xyz, N_atom.xyz],
              angle_ideal=0,
              weight=1).angle_model
          if good_hbond(angle, atom.distance(N_atom)):
            # print "Potential backwards bond:", atom.id_str(), atom.distance(N_atom), angle
            results.append((N_atom.id_str(), atom.distance(N_atom), angle, atom))
    self.atoms_around_cutted = filtered_atoms_around_cutted
    return results

  def identify_outliers(self):
    """
    Returns a dictionary of identified Cablam outliers, organized by protein chain.

    This function first identifies individual Cablam outliers and then groups
    consecutive outliers within the same protein chain into segments.

    :returns: outliers_by_chain (dict)
        A dictionary where:
        *   Keys: chain id (strings)
        *   Values: Are lists of lists. Each inner list (`comb`)
            represents a consecutive segment of Cablam outliers found
            within that specific protein chain.
            *   Elements within Inner Lists: Each element is a
                Cablam outlier objec:
                mmtbx.validation.cablam.cablam_result
                These objects contain detailed
                information about the individual residue identified as an
                outlier, such as its residue sequence number, insertion code,
                and references to its `residue` and `prevres` objects in the PDB hierarchy.

    Purpose:
    The structure of `outliers_by_chain` is designed to:
    *   Group neighboring Cablam outliers: This facilitates subsequent
        processing, such as applying idealizations or fixes to continuous
        problematic regions, rather than treating each outlier in isolation.

    Note:
    *   The function explicitly skips any Cablam outliers that belong to
        alternative conformations (i.e., those with a non-empty `altloc`
        attribute). Therefore, `outliers_by_chain` will only contain
        information about outliers from the primary (empty `altloc`) conformation.
    """
    self.starting_cablam_evaluation = cablamalyze(
        pdb_hierarchy=self.model.get_master_hierarchy(),
        outliers_only=True,
        out=null_out(),
        quiet=True,
        cablam_contours = self.cablam_contours,
        ca_contours = self.ca_contours,
        motif_contours = self.motif_contours)
    outliers_only = [x for x in self.starting_cablam_evaluation.results if x.feedback.cablam_outlier]# and x.feedback.c_alpha_geom_outlier]
    # Stores objects of
    outliers_by_chain = {}
    for k, g in itertools.groupby(outliers_only, key=lambda x: x.residue_id()[:2]):
      outliers_by_chain[k] = []
      comb = []
      for i in g:
        if i.altloc.strip() != '':
          print("  ", i, "<--- SKIPPING, alternative conformations.", file=self.log)
          continue
        if len(comb) == 0:
          comb = [i]
        else:
          if (i.resseq_as_int() - comb[-1].resseq_as_int() == 1 or
              (i.resseq_as_int() == comb[-1].resseq_as_int() and i.icode != comb[-1].icode)):
            comb.append(i)
          else:
            outliers_by_chain[k].append(comb)
            comb = [i]
        print("  ", i, file=self.log)
      # comb here is combination of adjacent outliers. If outlier is isolated - it
      # is included as a single one.
      outliers_by_chain[k].append(comb)
    return outliers_by_chain

  def get_results(self):
    flat_cablam_results = [item for sublist in self.cablam_results.values() for item in sublist]
    n_fixed_for_loop = sum(1 for ci, angle in flat_cablam_results if angle != 360 and not (ci.feedback.alpha or ci.feedback.beta or ci.feedback.threeten))
    n_fixed_for_alpha = sum(1 for ci, angle in flat_cablam_results if angle != 360 and ci.feedback.alpha)
    n_fixed_for_beta = sum(1 for ci, angle in flat_cablam_results if angle != 360 and ci.feedback.beta)
    n_fixed_for_threeten = sum(1 for ci, angle in flat_cablam_results if angle != 360 and ci.feedback.threeten)

    n_not_fixed_for_loop = sum(1 for ci, angle in flat_cablam_results if angle == 360 and not (ci.feedback.alpha or ci.feedback.beta or ci.feedback.threeten))
    n_not_fixed_for_alpha = sum(1 for ci, angle in flat_cablam_results if angle == 360 and ci.feedback.alpha)
    n_not_fixed_for_beta = sum(1 for ci, angle in flat_cablam_results if angle == 360 and ci.feedback.beta)
    n_not_fixed_for_threeten = sum(1 for ci, angle in flat_cablam_results if angle == 360 and ci.feedback.threeten)

    # print
    # assert n_fixed_for_loop + n_fixed_for_alpha + n_fixed_for_beta + n_fixed_for_threeten == self.n_rotated_residues

    return group_args(
      model = self.model,
      starting_cablam_evaluation = self.starting_cablam_evaluation,
      cablam_results = self.cablam_results,
      model_minimized = self.cablam_fixed_minimized,
      n_initial_cablam_outliers = self.starting_cablam_evaluation.summary_stats['cablam_outliers'],
      n_tried_residues = self.n_tried_residues,
      n_rotated_residues = self.n_rotated_residues,
      n_fixed_for_loop = n_fixed_for_loop,
      n_fixed_for_alpha = n_fixed_for_alpha,
      n_fixed_for_beta = n_fixed_for_beta,
      n_fixed_for_threeten = n_fixed_for_threeten,

      n_not_fixed_for_loop = n_not_fixed_for_loop,
      n_not_fixed_for_alpha = n_not_fixed_for_alpha,
      n_not_fixed_for_beta = n_not_fixed_for_beta,
      n_not_fixed_for_threeten = n_not_fixed_for_threeten,
      )
