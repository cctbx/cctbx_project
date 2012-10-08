from __future__ import division
import cctbx.geometry_restraints
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.utils import rotatable_bonds
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.refinement import fit_rotamers
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
import sys, math
from mmtbx.ncs import restraints
from libtbx.utils import Sorry
from mmtbx.torsion_restraints import utils
from mmtbx import ncs
import mmtbx.utils
from libtbx import Auto

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 2.5
   .type = float
   .short_caption = Restraint sigma (degrees)
 limit = 15.0
   .type = float
   .short_caption = Restraint limit (degrees)
 slack = 0.0
   .type = float
   .short_caption = Restraint slack (degrees)
 similarity = .80
   .type = float
   .short_caption = Sequence similarity cutoff
 fix_outliers = Auto
   .type = bool
   .short_caption = Fix rotamer outliers first
 check_rotamer_consistency = Auto
   .type = bool
   .short_caption = Check for consistency between NCS-related sidechains
   .help = Check for rotamer differences between NCS matched \
     sidechains and search for best fit amongst candidate rotamers
 target_damping = False
   .type = bool
   .expert_level = 1
 damping_limit = 10.0
   .type = float
   .expert_level = 1
 verbose = True
   .type = bool
 use_cc_for_target_angles = False
   .type = bool
   .expert_level = 4
 filter_phi_psi_outliers = True
   .type = bool
   .expert_level = 4
 silence_warnings = False
   .type = bool
   .expert_level = 4
 restraint_group
  .multiple=True
  .optional=True
  .caption = These atom selections define groups of residues whose dihedral \
    angles will be restrained to be similar.  This is normally done \
    automatically, and the restraints are designed to release dihedral angles \
    which are genuinely different.  You do not have to enter groups now \
    unless you wish to view and/or edit them prior to running phenix.refine.
  .short_caption=Torsion NCS restraint group
  .style = noauto box caption_img:icons/custom/ncs_tb.png
 {
  selection=None
    .type=atom_selection
    .short_caption=Restrained selection
    .multiple=True
    .input_size = 540
    .style = use_list
  b_factor_weight=10
    .type=float
    .short_caption = B factor weight
  coordinate_sigma=0.5
      .type = float
 }
""")

class torsion_ncs(object):
  def __init__(self,
               pdb_hierarchy,
               fmodel,
               geometry,
               sites_cart,
               params,
               b_factor_weight=None,
               coordinate_sigma=None,
               log=None):
    if(log is None): log = sys.stdout
    #parameter initialization
    if params.sigma is None or params.sigma < 0:
      raise Sorry("torsion NCS sigma parameter must be >= 0.0")
    self.sigma = params.sigma
    if params.limit is None or params.limit < 0:
      raise Sorry("torsion NCS limit parameter must be >= 0.0")
    self.limit = params.limit
    self.slack = params.slack
    self.use_cc_for_target_angles = params.use_cc_for_target_angles
    self.filter_phi_psi_outliers = params.filter_phi_psi_outliers
    self.b_factor_weight = b_factor_weight
    self.coordinate_sigma = coordinate_sigma
    self.fmodel = fmodel
    self.ncs_groups = []
    self.dp_ncs = []
    self.phi_list = []
    self.psi_list = []
    self.ncs_dihedral_proxies = None
    self.name_hash = utils.build_name_hash(pdb_hierarchy)
    self.segid_hash = utils.build_segid_hash(pdb_hierarchy)
    self.sym_atom_hash = utils.build_sym_atom_hash(pdb_hierarchy)
    self.params = params
    self.found_ncs = None
    self.log = log
    self.njump = 1
    self.min_length = 10
    self.sa = SidechainAngles(False)
    self.sidechain_angle_hash = self.build_sidechain_angle_hash()
    self.c_alpha_hinges = None
    self.r = rotalyze()
    print >> self.log, "Determining NCS matches..."
    dp_hash = {}
    self.use_segid = False
    i_seq_hash = utils.build_i_seq_hash(pdb_hierarchy)
    chain_hash = utils.build_chain_hash(pdb_hierarchy)
    name_hash = utils.build_name_hash(pdb_hierarchy)
    element_hash = utils.build_element_hash(pdb_hierarchy)
    chains = pdb_hierarchy.models()[0].chains()
    sel_cache = pdb_hierarchy.atom_selection_cache()
    n_ncs_groups = 0
    for i_seq, group in enumerate(self.params.restraint_group):
      n_selections = 0
      for selection in group.selection:
        if(selection is not None):
          n_selections += 1
      if n_selections == 1:
        raise Sorry(
          "Torsion NCS restraint_groups require at least 2 selections")
      elif n_selections > 1:
        n_ncs_groups += 1
    if n_ncs_groups > 0:
      sequences = {}
      padded_sequences = {}
      structures = {}
      alignments = {}
      for restraint_group in params.restraint_group:
        for selection_i in restraint_group.selection:
          sel_atoms_i = (utils.phil_atom_selections_as_i_seqs_multiple(
                           cache=sel_cache,
                           string_list=[selection_i]))
          sel_seq, sel_seq_padded, sel_structures = \
            utils.extract_sequence_and_sites(
            pdb_hierarchy=pdb_hierarchy,
            selection=sel_atoms_i)
          sequences[selection_i] = sel_seq
          padded_sequences[selection_i] = sel_seq_padded
          structures[selection_i] = sel_structures
      for restraint_group in params.restraint_group:
        ncs_set = []
        for selection_i in restraint_group.selection:
          ncs_set.append(selection_i)
          for selection_j in restraint_group.selection:
            if selection_i == selection_j:
              continue
            seq_pair = (sequences[selection_i],
                        sequences[selection_j])
            seq_pair_padded = (padded_sequences[selection_i],
                               padded_sequences[selection_j])
            struct_pair = (structures[selection_i],
                           structures[selection_j])
            residue_match_map = \
              utils._alignment(
                params=params,
                sequences=seq_pair,
                padded_sequences=seq_pair_padded,
                structures=struct_pair,
                log=log)
            key = (selection_i, selection_j)
            alignments[key] = residue_match_map
        self.ncs_groups.append(ncs_set)
      self.alignments = alignments
    else:
      atom_labels = list(pdb_hierarchy.atoms_with_labels())
      segids = flex.std_string([ a.segid for a in atom_labels ])
      self.use_segid = not segids.all_eq('    ')
      ncs_groups_manager = get_ncs_groups(
          pdb_hierarchy=pdb_hierarchy,
          use_segid=self.use_segid,
          params=self.params,
          log=self.log)
      self.ncs_groups = ncs_groups_manager.ncs_groups
      self.alignments = ncs_groups_manager.alignments
      new_ncs_groups = None
      #sort NCS groups alphabetically
      def selection_sort(match_list):
        match_list.sort()
        return match_list[0]
      self.ncs_groups.sort(key=selection_sort)
      if len(self.ncs_groups) > 0:
        new_ncs_groups = "refinement {\n ncs {\n  torsion {\n"
        for ncs_set in self.ncs_groups:
          new_ncs_groups += "   restraint_group {\n"
          for chain in ncs_set:
            new_ncs_groups += "    selection = %s\n" % chain
          if self.b_factor_weight is not None:
            new_ncs_groups += \
              "    b_factor_weight = %f\n" % self.b_factor_weight
          if self.coordinate_sigma is not None:
            new_ncs_groups += \
              "    coordinate_sigma = %f\n" % self.coordinate_sigma
          new_ncs_groups += "   }\n"
        new_ncs_groups += "  }\n }\n}"
      self.found_ncs = new_ncs_groups

    if self.found_ncs is not None:
      for dp in geometry.dihedral_proxies:
        h_atom = False
        for i_seq in dp.i_seqs:
          if element_hash[i_seq] == " H":
            h_atom = True
        if not h_atom:
          dp_hash[dp.i_seqs] = dp

      cbetadev_hash = utils.build_cbetadev_hash(
                        pdb_hierarchy=pdb_hierarchy)
      self.cbeta_proxies = []
      for cp in geometry.chirality_proxies:
        key = ""
        CAsite = None
        Csite = None
        Nsite = None
        CBsite = None
        CAkey = None
        Ckey = None
        Nkey = None
        CBkey = None
        cbeta = True
        wrong_atoms_for_c_beta = False
        used_atoms = []
        for i_seq in cp.i_seqs:
          if self.name_hash[i_seq][0:4] not in \
            [' CA ', ' N  ', ' C  ', ' CB ']:
            cbeta = False
            wrong_atoms_for_c_beta = True
          elif cbeta and self.name_hash[i_seq][0:4] == ' CA ':
            CAkey = self.name_hash[i_seq]
            CAsite = i_seq
          elif cbeta and self.name_hash[i_seq][0:4] == ' CB ':
            CBkey = self.name_hash[i_seq]
            CBsite = i_seq
            try:
              if float(cbetadev_hash[name_hash[i_seq][4:14]]) >= 0.25:
                cbeta = False
            except Exception:
                cbeta = False
          elif cbeta and self.name_hash[i_seq][0:4] == ' C  ':
            Ckey = self.name_hash[i_seq]
            Csite = i_seq
          elif cbeta and self.name_hash[i_seq][0:4] == ' N  ':
            Nkey = self.name_hash[i_seq]
            Nsite = i_seq
        if not cbeta and not wrong_atoms_for_c_beta:
          print >> self.log, "skipping C-beta restraint for %s" % \
            name_hash[CBsite][4:14]
        if cbeta:
          i_seqs = [Csite, Nsite, CAsite, CBsite]
          dp = cctbx.geometry_restraints.dihedral_proxy(
                 i_seqs=i_seqs,
                 angle_ideal=0.0,
                 weight=1.0)
          dp_hash[dp.i_seqs] = dp
          self.cbeta_proxies.append(dp)
          i_seqs = [Nsite, Csite, CAsite, CBsite]
          dp = cctbx.geometry_restraints.dihedral_proxy(
                 i_seqs=i_seqs,
                 angle_ideal=0.0,
                 weight=1.0)
          dp_hash[dp.i_seqs] = dp
          self.cbeta_proxies.append(dp)

      super_hash = {}
      res_match_master = {}
      res_to_selection_hash = {}
      for i, group in enumerate(self.ncs_groups):
        for chain_i in group:
          selection = utils.selection(
                       string=chain_i,
                       cache=sel_cache)
          c_atoms = pdb_hierarchy.select(selection).atoms()
          for atom in c_atoms:
            for chain_j in group:
              if chain_i == chain_j:
                continue
              res_key = self.name_hash[atom.i_seq][4:]
              atom_key = self.name_hash[atom.i_seq][0:4]
              j_match = None
              key = (chain_i, chain_j)
              cur_align = self.alignments.get(key)
              if cur_align is not None:
                j_match = cur_align.get(res_key)
              if j_match is not None:
                j_i_seq = i_seq_hash.get(atom_key+j_match)
                if j_i_seq is None:
                  continue
                if super_hash.get(atom.i_seq) is None:
                  super_hash[atom.i_seq] = dict()
                if super_hash.get(j_i_seq) is None:
                  super_hash[j_i_seq] = dict()
                super_hash[atom.i_seq][chain_j] = j_i_seq
                super_hash[j_i_seq][chain_i] = atom.i_seq
                if res_match_master.get(res_key) is None:
                  res_match_master[res_key] = []
                if res_match_master.get(j_match) is None:
                  res_match_master[j_match] = []
                if j_match not in res_match_master[res_key]:
                  res_match_master[res_key].append(j_match)
                if res_key not in res_match_master[j_match]:
                  res_match_master[j_match].append(res_key)
                res_to_selection_hash[res_key] = chain_i
                res_to_selection_hash[j_match] = chain_j
      self.res_match_master = res_match_master

      for dp in geometry.dihedral_proxies:
        temp = dict()
        for i_seq in dp.i_seqs:
          cur_matches = super_hash.get(i_seq)
          if cur_matches is None:
            continue
          for key in cur_matches.keys():
            try:
              temp[key].append(cur_matches[key])
            except Exception:
              temp[key] = []
              temp[key].append(cur_matches[key])
        dp_match = []
        dp_match.append(dp)
        for key in temp.keys():
          cur_dp_hash = dp_hash.get(tuple(temp[key]))
          if cur_dp_hash is not None:
            dp_match.append(cur_dp_hash)
            dp_hash[tuple(temp[key])] = None
        dp_hash[dp.i_seqs] = None
        if len(dp_match) > 1:
          self.dp_ncs.append(dp_match)

      for cb in self.cbeta_proxies:
        temp = dict()
        for i_seq in cb.i_seqs:
          cur_matches = super_hash.get(i_seq)
          if cur_matches is None:
            continue
          for key in cur_matches.keys():
            try:
              temp[key].append(cur_matches[key])
            except Exception:
              temp[key] = []
              temp[key].append(cur_matches[key])
        dp_match = []
        dp_match.append(cb)
        for key in temp.keys():
          cur_dp_hash = dp_hash.get(tuple(temp[key]))
          if cur_dp_hash is not None:
            dp_match.append(cur_dp_hash)
            dp_hash[tuple(temp[key])] = None
        dp_hash[cb.i_seqs] = None
        if len(dp_match) > 1:
          self.dp_ncs.append(dp_match)

      #initialize tracking hashes
      for dp_set in self.dp_ncs:
        for dp in dp_set:
          angle_atoms = self.get_torsion_atoms(dp)
          angle_resname = self.get_torsion_resname(dp)
          angle_id = self.get_torsion_id(dp)
          #phi
          if angle_atoms == ' C  '+' N  '+' CA '+' C  ':
            phi_id = self.get_torsion_id(dp, phi_psi=True)
            self.phi_list.append(dp.i_seqs)
          #psi
          elif angle_atoms == ' N  '+' CA '+' C  '+' N  ':
            psi_id = self.get_torsion_id(dp, phi_psi=True)
            self.psi_list.append(dp.i_seqs)

      match_counter = {}
      inclusive_range = {}
      for group in self.ncs_groups:
        cur_len = len(group)
        for chain in group:
          match_counter[chain] = cur_len
          inclusive_range[chain] = []

      matched = []
      ncs_match_hash = {}
      for dp_set in self.dp_ncs:
        key_set = []
        for dp in dp_set:
          if len(dp_set) < 2:
            continue
          cur_key = ""
          for i_seq in dp.i_seqs:
            cur_key += self.name_hash[i_seq]
          if cur_key[4:19] == cur_key[23:38] and \
             cur_key[4:19] == cur_key[42:57]:
            key_set.append(cur_key[4:19])
        if len(dp_set) == len(key_set):
          key_set.sort()
          master_key = None
          skip = False
          for i, key in enumerate(key_set):
            if i == 0:
              master_key = key
              if master_key in matched:
                skip = True
              elif ncs_match_hash.get(key) is None:
                ncs_match_hash[key] = []
              elif len(key_set) <= len(ncs_match_hash[key]):
                skip = True
              else:
                ncs_match_hash[key] = []
            elif not skip:
              ncs_match_hash[master_key].append(key)
              matched.append(key)
      self.ncs_match_hash = ncs_match_hash
      self.reduce_redundancies()

      for res in self.ncs_match_hash.keys():
        resnum = res[6:10]
        hash_key = res_to_selection_hash[res]
        cur_len = match_counter[hash_key]
        if len(self.ncs_match_hash[res]) == (cur_len - 1):
          inclusive_range[hash_key].append(int(resnum))
          for res2 in self.ncs_match_hash[res]:
            resnum2 = res2[6:10]
            hash_key = res_to_selection_hash[res2]
            inclusive_range[hash_key].append(int(resnum2))

      #determine ranges
      self.master_ranges = {}
      for key in inclusive_range.keys():
        current = None
        previous = None
        start = None
        stop = None
        self.master_ranges[key] = []
        inclusive_range[key].sort()
        for num in inclusive_range[key]:
          if previous == None:
            start = num
            previous = num
          elif num > (previous + 1):
            finish = previous
            self.master_ranges[key].append( (start, finish) )
            start = num
            finish = None
            previous = num
          else:
            previous = num
        if previous != None:
          finish = previous
          self.master_ranges[key].append( (start, finish) )

      if self.params.verbose:
        self.show_ncs_summary(log=log)
      print >> self.log, "Initializing torsion NCS restraints..."
      self.rama = ramalyze()
      self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          pdb_hierarchy=pdb_hierarchy,
                                          log=log)
    elif(not self.params.silence_warnings):
      print >> log, \
        "** WARNING: No torsion NCS found!!" + \
        "  Please check parameters. **"

  def show_ncs_summary(self, log=None):
    if(log is None): log = sys.stdout
    def get_key_chain_num(res):
      return res[4:]
    sorted_keys = sorted(self.ncs_match_hash, key=get_key_chain_num)
    print >> log, "--------------------------------------------------------"
    print >> log, "Torsion NCS Matching Summary:"
    for key in sorted_keys:
      if key.endswith("    "):
        print_line = key[:-4]
      else:
        print_line = key
      for match in self.ncs_match_hash[key]:
        if match.endswith("    "):
          print_line += " <=> %s" % (match[:-4])
        else:
          print_line += " <=> %s" % (match)
      print >> log, print_line
    print >> log, "--------------------------------------------------------"

  def reduce_redundancies(self):
    #clear out redundancies
    for key in self.ncs_match_hash.keys():
      for key2 in self.ncs_match_hash.keys():
        if key == key2:
          continue
        if key in self.ncs_match_hash[key2]:
          del self.ncs_match_hash[key]

  def get_torsion_atoms(self, dp):
    atoms = ''
    for i_seq in dp.i_seqs:
      atom_name = self.name_hash[i_seq][0:4]
      atoms += atom_name
    return atoms

  def get_torsion_resname(self, dp):
    resname = None
    for i_seq in dp.i_seqs:
      cur_resname = self.name_hash[i_seq][5:8]
      if resname == None:
        resname = cur_resname
      elif cur_resname != resname:
        return None
    return resname

  def get_torsion_id(self, dp, phi_psi=False, chi_only=False):
    id = None
    chi_atoms = False
    atom_list = []
    if phi_psi:
      return self.name_hash[dp.i_seqs[1]][4:]
    for i_seq in dp.i_seqs:
      cur_id = self.name_hash[i_seq][4:]
      atom = self.name_hash[i_seq][:4]
      atom_list.append(atom)
      if id == None:
        id = cur_id
      elif cur_id != id:
        return None
      if chi_only:
        if atom not in [' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OXT']:
          chi_atoms = True
    if chi_only and not chi_atoms:
      return None
    return id

  def generate_dihedral_ncs_restraints(
        self,
        sites_cart,
        pdb_hierarchy,
        log):
    self.ncs_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    target_map_data = None
    if self.fmodel is not None and self.use_cc_for_target_angles:
      target_map_data, residual_map_data = self.prepare_map(
                                             fmodel=self.fmodel)
    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)
    self.rama_outliers = None
    rama_outlier_list = []
    if self.filter_phi_psi_outliers:
      self.rama_outliers = \
        self.get_ramachandran_outliers(pdb_hierarchy)
      for outlier in self.rama_outliers.splitlines():
        temp = outlier.split(':')
        rama_outlier_list.append(temp[0])
    for dp_set in self.dp_ncs:
      if len(dp_set) < 2:
        continue
      angles = []
      cc_s = []
      is_rama_outlier = []
      rotamer_state = []
      wrap_hash = {}
      for i, dp in enumerate(dp_set):
        di = cctbx.geometry_restraints.dihedral(
               sites_cart=sites_cart, proxy=dp)
        angle = di.angle_model
        wrap_chis = self.is_symmetric_torsion(dp)
        if wrap_chis:
          if angle > 90.0 or angle < -90.0:
            sym_i_seq = dp.i_seqs[3] #4th atom
            swap_i_seq = self.sym_atom_hash.get(sym_i_seq)
            if swap_i_seq is not None:
              swap_i_seqs = (dp.i_seqs[0],
                             dp.i_seqs[1],
                             dp.i_seqs[2],
                             swap_i_seq)
              dp_temp = cctbx.geometry_restraints.dihedral_proxy(
                i_seqs=swap_i_seqs,
                angle_ideal=0.0,
                weight=1/self.sigma**2,
                limit=self.limit,
                top_out=TOP_OUT_FLAG,
                slack=self.slack)
              wrap_hash[i] = dp_temp
              di = cctbx.geometry_restraints.dihedral(
                     sites_cart=sites_cart, proxy=dp_temp)
              angle = di.angle_model
            else:
              angle = None
        angles.append(angle)
        rama_out = False
        if (dp.i_seqs in self.phi_list) or (dp.i_seqs in self.psi_list):
          angle_id = self.get_torsion_id(dp, phi_psi=True)
          key = angle_id[4:6].strip()+angle_id[6:10]+' '+angle_id[0:4]
          if key in rama_outlier_list:
            rama_out = True
        is_rama_outlier.append(rama_out)

        if target_map_data is not None:
          tor_iselection = flex.size_t()
          for i_seq in dp.i_seqs:
            tor_iselection.append(i_seq)
          tor_sites_cart = \
            sites_cart.select(tor_iselection)
          di_cc = self.get_sites_cc(sites_cart=tor_sites_cart,
                                    target_map_data=target_map_data)
          cc_s.append(di_cc)

        angle_id = self.get_torsion_id(dp, chi_only=True)
        if angle_id is not None:
          key = angle_id[4:6].strip()+angle_id[6:10]+' '+angle_id[0:4]
          rot_id = model_hash.get(key)
          rotamer_state.append(rot_id)

      target_angles = self.get_target_angles(
                        angles=angles,
                        cc_s=cc_s,
                        is_rama_outlier=is_rama_outlier,
                        rotamer_state=rotamer_state)
      for i, dp in enumerate(dp_set):
        target_angle = target_angles[i]
        angle_atoms = self.get_torsion_atoms(dp)
        angle_resname = self.get_torsion_resname(dp)
        angle_id = self.get_torsion_id(dp)
        cur_dict = self.sidechain_angle_hash.get(angle_resname)
        angle_name = None
        if cur_dict != None:
          angle_name = \
            cur_dict.get(angle_atoms)
        if target_angle is not None:
          angle_atoms = self.get_torsion_atoms(dp)
          angle_resname = self.get_torsion_resname(dp)
          angle_id = self.get_torsion_id(dp)
          cur_dict = self.sidechain_angle_hash.get(angle_resname)
          angle_name = None
          dp_sym = wrap_hash.get(i)
          if dp_sym is not None:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=dp_sym.i_seqs,
              angle_ideal=target_angle,
              weight=1/self.sigma**2,
              limit=self.limit,
              top_out=TOP_OUT_FLAG,
              slack=self.slack)
          else:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=dp.i_seqs,
              angle_ideal=target_angle,
              weight=1/self.sigma**2,
              limit=self.limit,
              top_out=TOP_OUT_FLAG,
              slack=self.slack)
          self.ncs_dihedral_proxies.append(dp_add)
    if len(self.ncs_dihedral_proxies) == 0:
      if (not self.params.silence_warnings) :
        print >> log, \
          "** WARNING: No torsion NCS found!!" + \
          "  Please check parameters. **"
    else:
      print >> log, \
        "Number of torsion NCS restraints: %d" \
          % len(self.ncs_dihedral_proxies)

  def update_dihedral_ncs_restraints(self,
                                     geometry,
                                     sites_cart,
                                     pdb_hierarchy,
                                     log=None):
    if log is None:
      log = sys.stdout
    make_sub_header(
      "Updating torsion NCS restraints",
      out=log)
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          pdb_hierarchy=pdb_hierarchy,
                                          log=log)
    self.add_ncs_dihedral_proxies(geometry=geometry)

  def is_symmetric_torsion(self, dp):
    i_seqs = dp.i_seqs
    resname = self.name_hash[i_seqs[0]][5:8].upper()
    if resname not in \
      ['ASP', 'GLU', 'PHE', 'TYR']:
      return False
    torsion_atoms = []
    for i_seq in i_seqs:
      name = self.name_hash[i_seq]
      atom = name[0:4]
      torsion_atoms.append(atom)
    if resname == 'ASP':
      if torsion_atoms == [' CA ', ' CB ', ' CG ', ' OD1'] or \
         torsion_atoms == [' CA ', ' CB ', ' CG ', ' OD2']:
        return True
    elif resname == 'GLU':
      if torsion_atoms == [' CB ', ' CG ', ' CD ', ' OE1'] or \
         torsion_atoms == [' CB ', ' CG ', ' CD ', ' OE2']:
        return True
    elif resname == 'PHE' or resname == 'TYR':
      if torsion_atoms == [' CA ', ' CB ',' CG ',' CD1'] or \
         torsion_atoms == [' CA ', ' CB ',' CG ',' CD2']:
        return True
    return False

  def get_target_angles(self,
                        angles,
                        cc_s,
                        is_rama_outlier,
                        rotamer_state):
    assert (len(rotamer_state) == len(angles)) or \
           (len(rotamer_state) == 0)
    clusters = {}
    used = []
    target_angles = angles[:]
    if is_rama_outlier.count(False) == 0:
      for i, target in enumerate(target_angles):
        target_angles[i] = None
      return target_angles
    max_i = None
    for i, cc in enumerate(cc_s):
      if is_rama_outlier[i]:
        continue
      if max_i is None:
        max_i = i
      elif max < cc:
        max_i = i
    for i, ang_i in enumerate(angles):
      if i in used:
        continue
      if ang_i is None:
        continue
      for j, ang_j in enumerate(angles):
        if i == j:
          continue
        elif j in used:
          continue
        elif ang_j is None:
          continue
        else:
          if len(rotamer_state) > 0:
            #require rotamers to match to restrain chi angles
            #WHAT ABOUT OUTLIERS?
            if (rotamer_state[i] != rotamer_state[j]) or \
               (rotamer_state[i] == "OUTLIER" and
                rotamer_state[j] == "OUTLIER"):
              continue
          if i not in used:
            clusters[i] = []
            clusters[i].append(i)
            clusters[i].append(j)
            used.append(i)
            used.append(j)
          else:
            clusters[i].append(j)
            used.append(j)
      if i not in used:
        clusters[i] = None
    for key in clusters.keys():
      cluster = clusters[key]
      if cluster is None:
        target_angles[key] = None
      else:
        cluster_angles = []
        for i in cluster:
          if is_rama_outlier[i]:
            cluster_angles.append(None)
          else:
            cluster_angles.append(angles[i])
        if max_i is not None:
          target_angle = angles[max_i]
        else:
          target_angle = utils.get_angle_average(cluster_angles)
        if self.params.target_damping:
          for c in cluster:
            c_dist = utils.angle_distance(angles[c], target_angle)
            if c_dist > self.params.damping_limit:
              d_target = \
                utils.get_angle_average([angles[c], target_angle])
              target_angles[c] = d_target
            else:
              target_angles[c] = target_angle
        else:
          for c in cluster:
            target_angles[c] = target_angle
    return target_angles

  def add_ncs_dihedral_proxies(self, geometry):
    geometry.ncs_dihedral_proxies= \
      self.ncs_dihedral_proxies

  def rotamer_correction_init(self, xray_structure, pdb_hierarchy):
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.unit_cell = xray_structure.unit_cell()
    self.exclude_free_r_reflections = False
    self.torsion_params = \
      fit_rotamers.torsion_search_params().extract().torsion_search
    self.torsion_params.range_start = -10.0
    self.torsion_params.range_stop = 10.0
    self.torsion_params.step = 1.0
    self.torsion_params.min_angle_between_solutions = 0.5
    self.c_alpha_hinges = get_c_alpha_hinges(xray_structure=xray_structure,
                                        pdb_hierarchy=pdb_hierarchy)

  def prepare_map(self, fmodel):
    map_obj = fmodel.electron_density_map()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "2mFo-DFc", use_all_data=(
        not self.exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "mFo-DFc", use_all_data=(
        not self.exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    residual_map_data = fft_map.real_map_unpadded()
    return target_map_data, residual_map_data

  def get_ramachandran_outliers(self, pdb_hierarchy):
    rama_outliers, output_list = \
      self.rama.analyze_pdb(hierarchy=pdb_hierarchy,
                            outliers_only=True)
    return rama_outliers

  def get_rotamer_data(self, pdb_hierarchy):
    rot_list_model, coot_model = \
      self.r.analyze_pdb(hierarchy=pdb_hierarchy)
    model_hash = {}
    model_score = {}
    all_rotamers = {}
    model_chis = {}
    for line in rot_list_model.splitlines():
      res, occ, rotamericity, chi1, chi2, chi3, chi4, name = line.split(':')
      model_hash[res]=name
      model_score[res]=rotamericity
    for key in self.res_match_master.keys():
      res_key = key[5:10]+' '+key[0:4]
      all_rotamers[res_key] = []
      model_rot = model_hash.get(res_key)
      if model_rot is not None and model_rot != "OUTLIER":
        all_rotamers[res_key].append(model_rot)
      for match_res in self.res_match_master[key]:
        j_key = match_res[5:10]+' '+match_res[0:4]
        j_rot = model_hash.get(j_key)
        if j_rot is not None and j_rot != "OUTLIER":
          if j_rot not in all_rotamers[res_key]:
            all_rotamers[res_key].append(j_rot)

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = \
              self.r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = \
                  self.r.sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
                  model_chis[key] = chis
              except Exception:
                print >> self.log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
    return model_hash, model_score, all_rotamers, model_chis

  def fix_rotamer_outliers(self,
                           xray_structure,
                           geometry_restraints_manager,
                           pdb_hierarchy,
                           outliers_only=False,
                           log=None,
                           quiet=False):
    self.last_round_outlier_fixes = 0
    if self.c_alpha_hinges is None:
      self.rotamer_correction_init(xray_structure, pdb_hierarchy)
    sites_cart = xray_structure.sites_cart()
    for atom in pdb_hierarchy.atoms():
      i_seq = atom.i_seq
      atom.xyz = sites_cart[i_seq]
    selection_radius = 5
    fmodel = self.fmodel
    if(log is None): log = self.log
    make_sub_header(
      "Correcting NCS rotamer outliers",
      out=log)

    target_map_data, residual_map_data = self.prepare_map(fmodel=fmodel)

    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)

    fix_list = {}
    rotamer_targets = {}

    for key in self.res_match_master.keys():
      res_key = key[5:10]+' '+key[0:4]
      model_rot = model_hash.get(res_key)
      if model_rot == "OUTLIER":
        rotamer = None
        score = 0.0
        for match_res in self.res_match_master[key]:
          j_key = match_res[5:10]+' '+match_res[0:4]
          j_rot = model_hash.get(j_key)
          j_score = model_score.get(j_key)
          if j_rot is not None and j_score is not None:
            if j_rot != "OUTLIER":
              if rotamer == None:
                rotamer = j_key
                score = j_score
                target = j_rot
              else:
                if j_score > score:
                  rotamer = j_key
                  score = j_score
                  target = j_rot
        if rotamer != None:
          fix_list[res_key] = rotamer
          rotamer_targets[res_key] = target

    sites_cart_moving = xray_structure.sites_cart()
    sites_cart_backup = sites_cart_moving.deep_copy()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          all_dict = \
            self.r.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY", "HOH"]:
              continue
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            if key in fix_list.keys():
              model_rot, m_chis, value = self.r.evaluate_rotamer(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  sites_cart=sites_cart_moving)
              residue_name = key[-3:]
              cur_rotamer = rotamer_targets[key]
              r_chis = self.r.sa.get_rotamer_angles(
                             residue_name=residue_name,
                             rotamer_name=cur_rotamer)
              if m_chis is not None and r_chis is not None:
                status = self.rotamer_search(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  m_chis=m_chis,
                  r_chis=r_chis,
                  rotamer=cur_rotamer,
                  sites_cart_moving=sites_cart_moving,
                  xray_structure=xray_structure,
                  target_map_data=target_map_data,
                  residual_map_data=residual_map_data,
                  key=key,
                  log=log)
                if status:
                  print >> log, "Set %s to %s rotamer" % \
                    (key, cur_rotamer)
                  self.last_round_outlier_fixes += 1

  def get_sites_cc(self,
                   sites_cart,
                   target_map_data):
    t = fit_rotamers.target(sites_cart,
                            self.unit_cell,
                            target_map_data)
    return t


  def get_sidechain_map_correlation(self,
                                    xray_structure,
                                    pdb_hierarchy):
    map_cc_hash = {}
    sigma_cutoff_hash = {}
    fmodel = self.fmodel
    target_map_data, residual_map_data = self.prepare_map(fmodel=fmodel)
    sites_cart_moving = xray_structure.sites_cart()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          all_dict = self.r.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY", "HOH"]:
              continue
            key = atom_group.atoms()[0].pdb_label_columns()[4:]+\
                  atom_group.atoms()[0].segid
            residue_iselection = atom_group.atoms().extract_i_seq()
            sidechain_only_iselection = flex.size_t()
            for i_seq in residue_iselection:
              atom_name = self.name_hash[i_seq][0:4]
              if atom_name not in [' N  ', ' CA ', ' C  ', ' O  ']:
                sidechain_only_iselection.append(i_seq)
            sites_cart_residue = \
              sites_cart_moving.select(sidechain_only_iselection)
            t_test = self.get_sites_cc(sites_cart_residue,
                                       target_map_data)
            map_cc_hash[key] = t_test
            sigma_state = fit_rotamers.all_sites_above_sigma_cutoff(
                            sites_cart_residue,
                            self.unit_cell,
                            target_map_data,
                            1.0)
            sigma_cutoff_hash[key] = sigma_state
    return map_cc_hash, sigma_cutoff_hash

  def fix_rotamer_consistency(self,
                              xray_structure,
                              geometry_restraints_manager,
                              pdb_hierarchy,
                              log=None,
                              quiet=False):
    self.last_round_rotamer_changes = 0
    if self.c_alpha_hinges is None:
      self.rotamer_correction_init(xray_structure, pdb_hierarchy)
    sites_cart = xray_structure.sites_cart()
    for atom in pdb_hierarchy.atoms():
      i_seq = atom.i_seq
      atom.xyz = sites_cart[i_seq]
    fmodel = self.fmodel
    if(log is None): log = self.log
    make_sub_header(
      "Checking NCS rotamer consistency",
      out=log)
    exclude_free_r_reflections = False
    rot_list_model, coot_model = \
      self.r.analyze_pdb(hierarchy=pdb_hierarchy)

    target_map_data, residual_map_data = self.prepare_map(fmodel=fmodel)

    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)

    sites_cart_moving = xray_structure.sites_cart()
    map_cc_hash, sigma_cutoff_hash = \
      self.get_sidechain_map_correlation(xray_structure, pdb_hierarchy)
    cc_candidate_list = []
    for key in self.ncs_match_hash.keys():
      whole_set = []
      value = map_cc_hash.get(key)
      max = None
      max_key = None
      if value is not None:
        whole_set.append( (key, value) )
        max = value
        max_key = key
      for member in self.ncs_match_hash[key]:
        value = map_cc_hash.get(member)
        if value is not None:
          whole_set.append( (member, value) )
          if max is None:
            max = value
            max_key = member
          else:
            if value > max:
              max = value
              max_key = member
      if max is None or max <= 0.0:
        continue
      for set in whole_set:
        cur_key = set[0]
        cur_value  = set[1]
        #fudge factor to account for zero cur_value
        if cur_value <= 0.0:
          cur_value = 0.0001
        percentage = (max - cur_value) / cur_value
        if percentage > 0.2:
          cc_candidate_list.append(cur_key)

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          all_dict = self.r.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY", "HOH"]:
              continue
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            if key in all_rotamers:
              if (len(all_rotamers[key]) >= 2):
                cc_key = atom_group.atoms()[0].pdb_label_columns()[4:]+\
                  atom_group.atoms()[0].segid
                if cc_key not in cc_candidate_list:
                  continue
                model_rot, m_chis, value = self.r.evaluate_rotamer(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  sites_cart=sites_cart_moving)
                residue_name = key[-3:]
                # why do I not try to fix outliers here?
                if model_rot == "OUTLIER":
                  continue
                current_best = model_rot
                #C-alpha prep
                cur_ca = None
                for atom in atom_group.atoms():
                  if atom.name == " CA ":
                    cur_ca = atom.i_seq
                for cur_rotamer in all_rotamers.get(key):
                  if cur_rotamer == model_rot:
                    continue
                  r_chis = self.r.sa.get_rotamer_angles(
                             residue_name=residue_name,
                             rotamer_name=cur_rotamer)
                  if m_chis is not None and r_chis is not None:
                    status = self.rotamer_search(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      m_chis=m_chis,
                      r_chis=r_chis,
                      rotamer=cur_rotamer,
                      sites_cart_moving=sites_cart_moving,
                      xray_structure=xray_structure,
                      target_map_data=target_map_data,
                      residual_map_data=residual_map_data,
                      key=key,
                      log=log)
                    if status:
                      current_best = cur_rotamer
                      atom_dict = all_dict.get(atom_group.altloc)
                      m_chis = \
                        self.r.sa.measureChiAngles(atom_group,
                                                   atom_dict,
                                                   sites_cart_moving)
                if current_best != model_rot:
                  print >> self.log, "Set %s to %s rotamer" % \
                    (key,
                     current_best)
                  self.last_round_rotamer_changes += 1
                else:
                  rotamer, chis, value = self.r.evaluate_rotamer(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      sites_cart=sites_cart_moving)
                  assert rotamer == model_rot

  def rotamer_search(
        self,
        atom_group,
        all_dict,
        m_chis,
        r_chis,
        rotamer,
        sites_cart_moving,
        xray_structure,
        target_map_data,
        residual_map_data,
        key,
        log):
    include_ca_hinge = False
    axis_and_atoms_to_rotate, tardy_labels= \
      rotatable_bonds.axes_and_atoms_aa_specific(
          residue=atom_group,
          mon_lib_srv=self.mon_lib_srv,
          remove_clusters_with_all_h=True,
          include_labels=True,
          log=None)
    if (axis_and_atoms_to_rotate is None) :
      print >> log, "Skipped %s rotamer (TARDY error)" % key
      return False
    assert len(m_chis) == len(r_chis)
    #exclude H-only clusters if necessary
    while len(axis_and_atoms_to_rotate) > len(m_chis):
      axis_and_atoms_to_rotate = \
        axis_and_atoms_to_rotate[:-1]
    assert len(m_chis) == len(axis_and_atoms_to_rotate)
    counter = 0
    residue_iselection = atom_group.atoms().extract_i_seq()
    cur_ca = None
    ca_add = None
    ca_axes = []
    for atom in atom_group.atoms():
      if atom.name == " CA ":
        cur_ca = atom.i_seq
    if cur_ca is not None:
      cur_c_alpha_hinges = self.c_alpha_hinges.get(cur_ca)
      if cur_c_alpha_hinges is not None:
        residue_length = len(tardy_labels)
        for ca_pt in cur_c_alpha_hinges[0]:
          residue_iselection.append(ca_pt)
          tardy_labels.append(self.name_hash[ca_pt][0:4])
        for bb_pt in cur_c_alpha_hinges[1]:
          residue_iselection.append(bb_pt)
          tardy_labels.append(self.name_hash[bb_pt][0:4])
        end_pts = (residue_length, residue_length+1)
        group = []
        for i, value in enumerate(tardy_labels):
          if i not in end_pts:
            group.append(i)
        ca_add = [end_pts, group]
        ca_axes.append(ca_add)
        for ax in axis_and_atoms_to_rotate:
          ca_axes.append(ax)
    sites_cart_residue = \
      sites_cart_moving.select(residue_iselection)
    sites_cart_residue_start = sites_cart_residue.deep_copy()
    selection = flex.bool(
                  len(sites_cart_moving),
                  residue_iselection)
    rev_first_atoms = []

    rev_start = fit_rotamers.rotamer_evaluator(
      sites_cart_start = sites_cart_residue_start,
      unit_cell        = self.unit_cell,
      two_mfo_dfc_map  = target_map_data,
      mfo_dfc_map      = residual_map_data)

    for aa in axis_and_atoms_to_rotate:
      axis = aa[0]
      atoms = aa[1]
      new_xyz = flex.vec3_double()
      angle_deg = r_chis[counter] - m_chis[counter]
      if angle_deg < 0:
        angle_deg += 360.0
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
                    axis_point_1=sites_cart_residue[axis[0]],
                    axis_point_2=sites_cart_residue[axis[1]],
                    point=sites_cart_residue[atom],
                    angle=angle_deg, deg=True)
        sites_cart_residue[atom] = new_xyz
      counter += 1

    #***** TEST *****
    sites_cart_moving.set_selected(
      residue_iselection, sites_cart_residue)
    cur_rotamer, cur_chis, cur_value = self.r.evaluate_rotamer(
      atom_group=atom_group,
      all_dict=all_dict,
      sites_cart=sites_cart_moving)
    assert rotamer == cur_rotamer
    #****************

    if len(ca_axes) == 0:
      eval_axes = axis_and_atoms_to_rotate
    else:
      eval_axes = ca_axes
      include_ca_hinge = True
    for i_aa, aa in enumerate(eval_axes):
      if(i_aa == len(eval_axes)-1):
        sites_aa = flex.vec3_double()
        for aa_ in aa[1]:
          sites_aa.append(sites_cart_residue[aa_])
      elif i_aa == 0 and include_ca_hinge:
        sites_aa = flex.vec3_double()
        for aa_ in aa[1]:
          sites_aa.append(sites_cart_residue[aa_])
      else:
        sites_aa = flex.vec3_double([sites_cart_residue[aa[1][0]]])
      rev_i = fit_rotamers.rotamer_evaluator(
        sites_cart_start = sites_aa,
        unit_cell        = self.unit_cell,
        two_mfo_dfc_map  = target_map_data,
        mfo_dfc_map      = residual_map_data)
      rev_first_atoms.append(rev_i)

    rev = fit_rotamers.rotamer_evaluator(
      sites_cart_start = sites_cart_residue,
      unit_cell        = self.unit_cell,
      two_mfo_dfc_map  = target_map_data,
      mfo_dfc_map      = residual_map_data)

    residue_sites_best = sites_cart_residue.deep_copy()
    residue_sites_best, rotamer_id_best = \
      fit_rotamers.torsion_search(
        residue_evaluator=rev,
        cluster_evaluators=rev_first_atoms,
        axes_and_atoms_to_rotate=eval_axes,
        rotamer_sites_cart=sites_cart_residue,
        rotamer_id_best=rotamer,
        residue_sites_best=residue_sites_best,
        params = self.torsion_params,
        rotamer_id = rotamer,
        include_ca_hinge = include_ca_hinge)
    sites_cart_moving.set_selected(
        residue_iselection, residue_sites_best)
    xray_structure.set_sites_cart(sites_cart_moving)
    cur_rotamer, cur_chis, cur_value = self.r.evaluate_rotamer(
      atom_group=atom_group,
      all_dict=all_dict,
      sites_cart=sites_cart_moving)
    rotamer_match = (cur_rotamer == rotamer)
    if rev_start.is_better(sites_cart=residue_sites_best,
                           percent_cutoff=0.15, verbose=True) and \
       rotamer_match:
      sidechain_only_iselection = flex.size_t()
      for i_seq in residue_iselection:
        atom_name = self.name_hash[i_seq][0:4]
        if atom_name not in [' N  ', ' CA ', ' C  ', ' O  ',
                             ' OXT', ' H  ', ' HA ']:
          sidechain_only_iselection.append(i_seq)
      selection = flex.bool(
                    len(sites_cart_moving),
                    sidechain_only_iselection)
      selection_within = xray_structure.selection_within(
      radius    = 1.0,
      selection = selection)
      #check for bad steric clashes
      created_clash = False
      for i, state in enumerate(selection_within):
        if state:
          if i not in sidechain_only_iselection:
            #print >> self.log, "atom clash: ", self.name_hash[i]
            created_clash = True
      if created_clash:
        sites_cart_moving.set_selected(
          residue_iselection, sites_cart_residue_start)
        xray_structure.set_sites_cart(sites_cart_moving)
        return False
      return True
    else:
      sites_cart_moving.set_selected(
        residue_iselection, sites_cart_residue_start)
      xray_structure.set_sites_cart(sites_cart_moving)
      return False

  def process_ncs_restraint_groups(self, model, processed_pdb_file):
    log = self.log
    ncs_groups = ncs.restraints.groups()
    sites_cart = None

    for param_group in self.params.restraint_group:
      master = param_group.selection[0]
      selection_strings = []
      found_offset = False
      range_text = ""
      for range in self.master_ranges[master]:
        if range_text == "":
          range_text = "(resseq %d:%d" % (range[0], range[1])
        else:
          range_text += " or resseq %d:%d" % (range[0], range[1])
      range_text += ")"
      master = master + " and " + range_text
      for selection in param_group.selection[1:]:
        range_text = ""
        for range in self.master_ranges[selection]:
          if range_text == "":
            range_text = "(resseq %d:%d" % (range[0], range[1])
          else:
            range_text += " or resseq %d:%d" % (range[0], range[1])
        range_text += ")"
        temp_selection = selection + " and " + range_text
        selection_strings.append(temp_selection)
      group = ncs.restraints.group.from_atom_selections(
        processed_pdb              = processed_pdb_file,
        reference_selection_string = master,
        selection_strings          = selection_strings,
        coordinate_sigma           = param_group.coordinate_sigma,
        b_factor_weight            = param_group.b_factor_weight,
        special_position_warnings_only
          = False,
        log = log)
      ncs_groups.members.append(group)
      print >> log
    if (len(ncs_groups.members) == 0):
      print >> log, "No NCS restraint groups specified."
      print >> log
    else:
      model.restraints_manager.torsion_ncs_groups = ncs_groups

  def build_sidechain_angle_hash(self):
    sidechain_angle_hash = {}
    for key in self.sa.atomsForAngle.keys():
      resname = key[0:3].upper()
      if sidechain_angle_hash.get(resname) is None:
        sidechain_angle_hash[resname] = {}
      new_key = ''
      for atom in self.sa.atomsForAngle[key]:
        new_key += atom
      new_value = key[4:]
      sidechain_angle_hash[resname][new_key] = new_value
    #modifications
    sidechain_angle_hash['ILE'][' N   CA  CB  CG2'] = 'chi1'
    sidechain_angle_hash['THR'][' N   CA  CB  CG2'] = 'chi1'
    sidechain_angle_hash['VAL'][' N   CA  CB  CG2'] = 'chi1'
    return sidechain_angle_hash

  def get_number_of_restraints_per_group(self, pdb_hierarchy):
    torsion_counts = {}
    sel_cache = pdb_hierarchy.atom_selection_cache()
    for group in self.ncs_groups:
      for selection in group:
        sel_atoms_i = (utils.phil_atom_selections_as_i_seqs_multiple(
                         cache=sel_cache,
                         string_list=[selection]))
        torsion_counts[selection] = 0
        for dp in self.ncs_dihedral_proxies:
          if dp.i_seqs[0] in sel_atoms_i:
            torsion_counts[selection] += 1
    return torsion_counts

  def get_torsion_rmsd(self, sites_cart):
    self.histogram_under_limit = None
    self.histogram_over_limit = None
    self.torsion_rmsd = None
    self.all_torsion_rmsd = None
    dp_proxies_under_limit = cctbx.geometry_restraints.shared_dihedral_proxy()
    dp_proxies_over_limit = cctbx.geometry_restraints.shared_dihedral_proxy()
    for dp in self.ncs_dihedral_proxies:
      di = cctbx.geometry_restraints.dihedral(
             sites_cart=sites_cart, proxy=dp)
      delta = abs(di.delta)
      if delta <= self.limit:
        dp_proxies_under_limit.append(dp)
      else:
        dp_proxies_over_limit.append(dp)
    torsion_deltas_under_limit = cctbx.geometry_restraints.dihedral_deltas(
                       sites_cart = sites_cart,
                       proxies = dp_proxies_under_limit)
    torsion_deltas_over_limit = cctbx.geometry_restraints.dihedral_deltas(
                       sites_cart = sites_cart,
                       proxies = dp_proxies_over_limit)
    torsion_deltas_all = cctbx.geometry_restraints.dihedral_deltas(
                       sites_cart = sites_cart,
                       proxies = self.ncs_dihedral_proxies)
    if len(torsion_deltas_under_limit) > 0:
      self.histogram_under_limit = \
        flex.histogram(
          data=flex.abs(torsion_deltas_under_limit),
          data_min=0.0,
          data_max=self.limit,
          n_slots=10)
      self.torsion_rmsd = self.calculate_torsion_rmsd(
                            deltas=torsion_deltas_under_limit)
    if ( (len(torsion_deltas_over_limit) > 0) and
         (limit < 180.0) ):
      self.histogram_over_limit = \
        flex.histogram(
          data=flex.abs(torsion_deltas_over_limit),
          data_min=self.limit,
          data_max=math.ceil(
          max(flex.abs(torsion_deltas_over_limit))),
          n_slots=10)
    if len(torsion_deltas_all) > 0:
      self.all_torsion_rmsd = self.calculate_torsion_rmsd(
                                deltas=torsion_deltas_all)

  def calculate_torsion_rmsd(self, deltas):
    assert len(deltas) > 0
    delta_sq_sum = 0.0
    for delta in deltas:
      delta_sq_sum += ( abs(delta)**2 )
    return math.sqrt(delta_sq_sum / len(deltas))

#split out functions
class get_ncs_groups(object):
  def __init__(self,
               pdb_hierarchy,
               use_segid,
               params,
               log):
    ncs_groups = []
    alignments = {}
    used_chains = []
    pair_hash = {}
    chains = pdb_hierarchy.models()[0].chains()
    am = utils.alignment_manager(pdb_hierarchy, use_segid, log)

    for i, chain_i in enumerate(chains):
      found_conformer = False
      for conformer in chain_i.conformers():
        if not conformer.is_protein() and not conformer.is_na():
          continue
        else:
          found_conformer = True
      if not found_conformer:
        continue
      segid_i = utils.get_unique_segid(chain_i)
      if segid_i == None:
        continue
      if (use_segid) :
        chain_i_str = "chain '%s' and segid '%s'" % \
          (chain_i.id, segid_i)
      else :
        chain_i_str = "chain '%s'" % chain_i.id
      for chain_j in chains[i+1:]:
        found_conformer = False
        for conformer in chain_j.conformers():
          if not conformer.is_protein() and not conformer.is_na():
            continue
          else:
            found_conformer = True
        if not found_conformer:
          continue
        segid_j = utils.get_unique_segid(chain_j)
        if segid_j == None:
          continue
        if (use_segid) :
          chain_j_str = "chain '%s' and segid '%s'" % (chain_j.id, segid_j)
        else :
          chain_j_str = "chain '%s'" % chain_j.id
        seq_pair = (am.sequences[chain_i_str],
                    am.sequences[chain_j_str])
        seq_pair_padded = (am.padded_sequences[chain_i_str],
                           am.padded_sequences[chain_j_str])
        struct_pair = (am.structures[chain_i_str],
                       am.structures[chain_j_str])
        residue_match_map = \
          utils._alignment(
            params=params,
            sequences=seq_pair,
            padded_sequences=seq_pair_padded,
            structures=struct_pair,
            log=log)
        key = (chain_i_str, chain_j_str)
        alignments[key] = residue_match_map
        if ( min(len(residue_match_map),
                 chain_i.residue_groups_size(),
                 chain_j.residue_groups_size()) \
             / max(len(residue_match_map),
                   chain_i.residue_groups_size(),
                   chain_j.residue_groups_size()) \
             >= params.similarity ):
          pair_key = (chain_i.id, segid_i)
          match_key = (chain_j.id, segid_j)
          if used_chains is not None:
            if match_key in used_chains:
              continue
          assign_key = None
          if pair_key in used_chains:
            for group_key in pair_hash.keys():
              if pair_key in pair_hash[group_key]:
                assign_key = group_key
          if assign_key is None:
            assign_key = pair_key
          if (not assign_key in pair_hash) :
            pair_hash[assign_key] = []
          pair_hash[assign_key].append(match_key)
          used_chains.append(match_key)

    for key in pair_hash.keys():
      ncs_set = []
      if (use_segid) :
        chain_str = "chain '%s' and segid '%s'" % (key[0], key[1])
      else :
        chain_str = "chain '%s'" % (key[0])
      ncs_set.append(chain_str)
      for add_chain in pair_hash[key]:
        if (use_segid) :
          chain_str = "chain '%s' and segid '%s'" % \
            (add_chain[0], add_chain[1])
        else :
          chain_str = "chain '%s'" % (add_chain[0])
        ncs_set.append(chain_str)
      ncs_groups.append(ncs_set)

    self.alignments = alignments
    self.ncs_groups = ncs_groups

def determine_ncs_groups(pdb_hierarchy,
                         params=None,
                         log=None):
  pdb_hierarchy.reset_i_seq_if_necessary()
  if params is None:
    params = torsion_ncs_params.extract()
  if log is None:
    log = sys.stdout
  atom_labels = list(pdb_hierarchy.atoms_with_labels())
  segids = flex.std_string([ a.segid for a in atom_labels ])
  use_segid = not segids.all_eq('    ')
  ncs_groups_manager = get_ncs_groups(
                         pdb_hierarchy=pdb_hierarchy,
                         use_segid=use_segid,
                         params=params,
                         log=log)
  return ncs_groups_manager.ncs_groups

def check_residues_are_connected (ca_1, ca_2, max_sep=4.0, min_sep=0.) :
  from scitbx import matrix
  ca_1_mat = matrix.col(ca_1.xyz)
  ca_2_mat = matrix.col(ca_2.xyz)
  dist = (ca_1_mat - ca_2_mat).length()
  if (dist > max_sep) or (dist < min_sep) :
    return False
  return True

def get_c_alpha_hinges(xray_structure,
                       pdb_hierarchy):
  c_alphas = []
  c_alpha_hinges = {}
  sites_cart = xray_structure.sites_cart()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          cur_ca = None
          cur_c = None
          cur_o = None
          cur_n = None
          cur_h = None
          for atom in atom_group.atoms():
            if atom.name == " CA ":
              cur_ca = atom
            elif atom.name == " C  ":
              cur_c = atom
            elif atom.name == " N  ":
              cur_n = atom
            elif atom.name == " O  ":
              cur_o = atom
            elif atom.name == " H  ":
              cur_h = atom
          if cur_ca is not None and cur_c is not None and \
             cur_n is not None and cur_o is not None:
            moving_tpl = (cur_n, cur_c, cur_o)
            if cur_h is not None:
              moving_tpl += tuple([cur_h])
            c_alphas.append( (cur_ca, moving_tpl) )
  for i, ca in enumerate(c_alphas):
    if i < 1 or i == (len(c_alphas)-1):
      continue
    current = ca
    previous = c_alphas[i-1]
    next = c_alphas[i+1]
    prev_connected = check_residues_are_connected(previous[0], current[0])
    next_connected = check_residues_are_connected(current[0], next[0])
    if prev_connected and next_connected:
      nodes = (previous[0].i_seq, next[0].i_seq)
      moving = (previous[1][1].i_seq, previous[1][2].i_seq, next[1][0].i_seq)
      if len(next[1]) > 3:
        moving += tuple([next[1][3].i_seq])
      c_alpha_hinges[current[0].i_seq] = [nodes, moving]
  return c_alpha_hinges

# XXX wrapper for running in Phenix GUI
class _run_determine_ncs_groups (object) :
  def __init__ (self, params, pdb_hierarchy) :
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy

  def __call__ (self, *args, **kwds) :
    return determine_ncs_groups(
      params=self.params,
      pdb_hierarchy=self.pdb_hierarchy)
