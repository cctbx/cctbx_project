import cctbx.geometry_restraints
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.refinement import fit_rotamers, real_space
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
import sys
from mmtbx.ncs import restraints
from libtbx.utils import Sorry
from mmtbx.torsion_restraints import utils
from mmtbx import ncs
import mmtbx.utils
from libtbx import Auto

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 3.0
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
   .short_caption = Check for rotamer differences between NCS matched \
     sidechains and search for best fit amongst candidate rotamers
 target_damping = False
   .type = bool
   .expert_level = 1
 damping_limit = 10.0
   .type = float
   .expert_level = 1
 verbose = True
   .type = bool
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
    self.sigma = params.sigma
    self.limit = params.limit
    self.slack = params.slack
    self.b_factor_weight = b_factor_weight
    self.coordinate_sigma = coordinate_sigma
    self.pdb_hierarchy = pdb_hierarchy
    self.fmodel = fmodel
    self.ncs_groups = []
    self.dp_ncs = []
    self.chi_tracker = {}
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
      atom_labels = list(self.pdb_hierarchy.atoms_with_labels())
      segids = flex.std_string([ a.segid for a in atom_labels ])
      self.use_segid = not segids.all_eq('    ')
      ncs_groups_manager = get_ncs_groups(
          pdb_hierarchy=self.pdb_hierarchy,
          use_segid=self.use_segid,
          params=self.params,
          log=self.log)
      self.ncs_groups = ncs_groups_manager.ncs_groups
      self.alignments = ncs_groups_manager.alignments
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

    for dp in geometry.dihedral_proxies:
      h_atom = False
      for i_seq in dp.i_seqs:
        if element_hash[i_seq] == " H":
          h_atom = True
      if not h_atom:
        dp_hash[dp.i_seqs] = dp

    cbetadev_hash = utils.build_cbetadev_hash(
                      pdb_hierarchy=self.pdb_hierarchy)
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
      for i_seq in cp.i_seqs:
        if self.name_hash[i_seq][0:4] not in \
          [' CA ', ' N  ', ' C  ', ' CB ']:
          cbeta = False
        if self.name_hash[i_seq][0:4] == ' CA ':
          CAkey = self.name_hash[i_seq]
          CAsite = i_seq
        elif self.name_hash[i_seq][0:4] == ' CB ':
          CBkey = self.name_hash[i_seq]
          CBsite = i_seq
          try:
            if float(cbetadev_hash[name_hash[i_seq][4:14]]) >= 0.25:
              c_beta = False
              print >> self.log, "skipping C-beta restraint for %s" % \
                name_hash[i_seq][4:14]
          except Exception:
              c_beta = False
        elif self.name_hash[i_seq][0:4] == ' C  ':
          Ckey = self.name_hash[i_seq]
          Csite = i_seq
        elif self.name_hash[i_seq][0:4] == ' N  ':
          Nkey = self.name_hash[i_seq]
          Nsite = i_seq
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
        c_atoms = self.pdb_hierarchy.select(selection).atoms()
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

    for dp_set in self.dp_ncs:
      for dp in dp_set:
        angle_atoms = self.get_torsion_atoms(dp)
        angle_resname = self.get_torsion_resname(dp)
        angle_id = self.get_torsion_id(dp)
        cur_dict = self.sidechain_angle_hash.get(angle_resname)
        if cur_dict != None:
          angle_name = \
            cur_dict.get(angle_atoms)
          if angle_name != None:
            if self.chi_tracker.get(angle_id) is None:
              self.chi_tracker[angle_id] = {}
            self.chi_tracker[angle_id][angle_name] = False

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
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)

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

  def get_torsion_id(self, dp):
    id = None
    for i_seq in dp.i_seqs:
      cur_id = self.name_hash[i_seq][4:]
      if id == None:
        id = cur_id
      elif cur_id != id:
        return None
    return id

  def reset_chi_tracker(self):
    for key1 in self.chi_tracker.keys():
      for key2 in self.chi_tracker[key1].keys():
        self.chi_tracker[key1][key2] = False

  def generate_dihedral_ncs_restraints(self, sites_cart, log):
    self.ncs_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    #print self.sym_atom_hash
    for dp_set in self.dp_ncs:
      if len(dp_set) < 2:
        continue
      angles = []
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
      target_angles = self.get_target_angles(
                        angles=angles)
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
          if cur_dict is not None:
            angle_name = \
              cur_dict.get(angle_atoms)
          if angle_name is not None:
            if angle_name[-1:] == '1':
              self.chi_tracker[angle_id][angle_name] = True
            elif angle_name[-1:].isdigit() :
              current_chi_number = int(angle_name[-1:])
              previous_chi_number = current_chi_number - 1
              previous_chi_id = "chi%d" % previous_chi_number
              previous_chi_state = \
                self.chi_tracker[angle_id].get(previous_chi_id)
              if previous_chi_state == False:
                #if angle_id.endswith("    "):
                #  print >> self.log, \
                #    "skipping %s, chi%d" % \
                #    (angle_id[:-4], current_chi_number)
                #else:
                #  print >> self.log, \
                #    "skipping %s, chi%d" % \
                #    (angle_id, current_chi_number)
                continue
              else:
                self.chi_tracker[angle_id][angle_name] = True
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
                                     log=None):
    if log is None:
      log = sys.stdout
    make_sub_header(
      "Updating torsion NCS restraints",
      out=log)
    self.reset_chi_tracker()
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
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

  def get_target_angles(self, angles):
    clusters = {}
    used = []
    target_angles = angles[:]
    wrapped = [False]*len(angles)
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
          #if utils.angle_distance(ang_i, ang_j) <= self.params.cutoff:
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
          cluster_angles.append(angles[i])
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

  def fix_rotamer_outliers(self,
                           xray_structure,
                           geometry_restraints_manager,
                           outliers_only=False,
                           log=None,
                           quiet=False):
    self.last_round_outlier_fixes = 0
    pdb_hierarchy=self.pdb_hierarchy
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
    r = rotalyze()
    sa = SidechainAngles(False)
    exclude_free_r_reflections = False
    unit_cell = xray_structure.unit_cell()
    mon_lib_srv = mmtbx.monomer_library.server.server()
    rot_list_model, coot_model = r.analyze_pdb(hierarchy=pdb_hierarchy)

    map_obj = self.fmodel.electron_density_map()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "2mFo-DFc", use_all_data=(not exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "mFo-DFc", use_all_data=(not exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    residual_map_data = fft_map.real_map_unpadded()

    brm = real_space.box_refinement_manager(
            xray_structure = xray_structure,
            pdb_hierarchy = pdb_hierarchy,
            target_map = target_map_data,
            geometry_restraints_manager = geometry_restraints_manager)

    model_hash = {}
    model_score = {}
    model_chis = {}
    fix_list = {}
    rotamer_targets = {}
    all_rotamers = {}
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

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = r.sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
                  model_chis[key] = chis
              except Exception:
                print >> log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)

    sites_cart_moving = xray_structure.sites_cart()
    sites_cart_backup = sites_cart_moving.deep_copy()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          all_dict = r.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname == "PRO":
              continue
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            if key in fix_list.keys():
              m_chis = model_chis.get(key)
              r_chis = model_chis.get(fix_list[key])
              if m_chis is not None and r_chis is not None:
                axis_and_atoms_to_rotate= \
                  fit_rotamers.axes_and_atoms_aa_specific(
                      residue=atom_group,
                      mon_lib_srv=mon_lib_srv,
                      remove_clusters_with_all_h=True,
                      log=None)
                if (axis_and_atoms_to_rotate is None) :
                  print >> log, "Skipped %s rotamer (TARDY error)" % key
                  continue
                assert len(m_chis) == len(r_chis)
                #exclude H-only clusters if necessary
                while len(axis_and_atoms_to_rotate) > len(m_chis):
                  axis_and_atoms_to_rotate = \
                    axis_and_atoms_to_rotate[:-1]
                assert len(m_chis) == len(axis_and_atoms_to_rotate)
                counter = 0
                residue_iselection = atom_group.atoms().extract_i_seq()
                sites_cart_residue = \
                  sites_cart_moving.select(residue_iselection)
                sites_cart_residue_start = sites_cart_residue.deep_copy()
                selection = flex.bool(
                              len(sites_cart_moving),
                              residue_iselection)
                rev = fit_rotamers.rotamer_evaluator(
                  sites_cart_start = sites_cart_residue,
                  unit_cell        = unit_cell,
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

                sites_cart_moving.set_selected(
                    residue_iselection, sites_cart_residue)
                xray_structure.set_sites_cart(sites_cart_moving)
                rotamer, value = r.evaluate_rotamer(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  sites_cart=sites_cart_moving)
                brm.refine(selection=selection,
                           monitor_clashscore=False)
                sites_cart_refined_residue = \
                  brm.sites_cart.select(residue_iselection)
                rotamer, value = r.evaluate_rotamer(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  sites_cart=brm.sites_cart)
                xray_structure.set_sites_cart(brm.sites_cart)
                rotamer_match = (rotamer == rotamer_targets[key])
                if rev.is_better(sites_cart_refined_residue) and \
                   rotamer_match:
                  sites_cart_moving.set_selected(
                    residue_iselection, sites_cart_refined_residue)
                  xray_structure.set_sites_cart(sites_cart_moving)
                  mmtbx.utils.assert_xray_structures_equal(
                    x1 = xray_structure,
                    x2 = brm.xray_structure)
                  print >> log, "Set %s to %s rotamer" % \
                    (key, rotamer)
                  self.last_round_outlier_fixes += 1
                else:
                  sites_cart_moving.set_selected(
                    residue_iselection, sites_cart_residue_start)
                  xray_structure.set_sites_cart(sites_cart_moving)

  def fix_rotamer_consistency(self,
                              xray_structure,
                              geometry_restraints_manager,
                              log=None,
                              quiet=False):
    self.last_round_rotamer_changes = 0
    pdb_hierarchy=self.pdb_hierarchy
    sites_cart = xray_structure.sites_cart()
    for atom in pdb_hierarchy.atoms():
      i_seq = atom.i_seq
      atom.xyz = sites_cart[i_seq]
    selection_radius = 5
    fmodel = self.fmodel
    if(log is None): log = self.log
    make_sub_header(
      "Checking NCS rotamer consistency",
      out=log)
    r = rotalyze()
    sa = SidechainAngles(False)
    exclude_free_r_reflections = False
    unit_cell = xray_structure.unit_cell()
    mon_lib_srv = mmtbx.monomer_library.server.server()
    rot_list_model, coot_model = r.analyze_pdb(hierarchy=pdb_hierarchy)

    map_obj = self.fmodel.electron_density_map()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "2mFo-DFc", use_all_data=(not exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
    fft_map = map_obj.fft_map(resolution_factor = 1./4,
      map_type = "mFo-DFc", use_all_data=(not exclude_free_r_reflections))
    fft_map.apply_sigma_scaling()
    residual_map_data = fft_map.real_map_unpadded()

    brm = real_space.box_refinement_manager(
            xray_structure = xray_structure,
            pdb_hierarchy = pdb_hierarchy,
            target_map = target_map_data,
            geometry_restraints_manager = geometry_restraints_manager)

    model_hash = {}
    model_score = {}
    model_chis = {}
    fix_list = {}
    rotamer_targets = {}
    all_rotamers = {}
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
            all_dict = r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = r.sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
                  model_chis[key] = chis
              except Exception:
                print >> log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)

    sites_cart_moving = xray_structure.sites_cart()

    #fix rotamer consistency
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          all_dict = r.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname == "PRO":
              continue
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            if key in all_rotamers:
              if (len(all_rotamers[key]) >= 2):
                m_chis = model_chis.get(key)
                residue_name = key[-3:]
                model_rot = model_hash.get(key)
                current_best = model_rot
                for cur_rotamer in all_rotamers.get(key):
                  if cur_rotamer == model_rot:
                    continue
                  r_chis = sa.get_rotamer_angles(
                             residue_name=residue_name,
                             rotamer_name=cur_rotamer)
                  if m_chis is not None and r_chis is not None:
                    axis_and_atoms_to_rotate= \
                      fit_rotamers.axes_and_atoms_aa_specific(
                          residue=atom_group,
                          mon_lib_srv=mon_lib_srv,
                          remove_clusters_with_all_h=True,
                          log=None)
                    if (axis_and_atoms_to_rotate is None) :
                      print >> log, "Skipped %s rotamer (TARDY error)" % key
                      continue
                    assert len(m_chis) == len(r_chis)
                    #exclude H-only clusters if necessary
                    while len(axis_and_atoms_to_rotate) > len(m_chis):
                      axis_and_atoms_to_rotate = \
                        axis_and_atoms_to_rotate[:-1]
                    assert len(m_chis) == len(axis_and_atoms_to_rotate)
                    counter = 0
                    residue_iselection = atom_group.atoms().extract_i_seq()
                    sites_cart_residue = \
                      sites_cart_moving.select(residue_iselection)
                    sites_cart_residue_start = sites_cart_residue.deep_copy()
                    selection = flex.bool(
                                  len(sites_cart_moving),
                                  residue_iselection)
                    rev = fit_rotamers.rotamer_evaluator(
                      sites_cart_start = sites_cart_residue,
                      unit_cell        = unit_cell,
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

                    sites_cart_moving.set_selected(
                        residue_iselection, sites_cart_residue)
                    xray_structure.set_sites_cart(sites_cart_moving)
                    rotamer, value = r.evaluate_rotamer(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      sites_cart=sites_cart_moving)
                    brm.refine(selection=selection,
                               monitor_clashscore=False)
                    sites_cart_refined_residue = \
                      brm.sites_cart.select(residue_iselection)
                    rotamer, value = r.evaluate_rotamer(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      sites_cart=brm.sites_cart)
                    xray_structure.set_sites_cart(brm.sites_cart)
                    rotamer_match = (rotamer == cur_rotamer)
                    if rev.is_better(sites_cart_refined_residue) and \
                      rotamer_match:
                      sites_cart_moving.set_selected(
                        residue_iselection, sites_cart_refined_residue)
                      xray_structure.set_sites_cart(sites_cart_moving)
                      mmtbx.utils.assert_xray_structures_equal(
                        x1 = xray_structure,
                        x2 = brm.xray_structure)
                      current_best = cur_rotamer
                    else:
                      sites_cart_moving.set_selected(
                        residue_iselection, sites_cart_residue_start)
                      xray_structure.set_sites_cart(sites_cart_moving)
                if current_best != model_rot:
                  print >> self.log, "Set %s to %s rotamer" % \
                    (key,
                     current_best)
                  self.last_round_rotamer_changes += 1
                else:
                  rotamer, value = r.evaluate_rotamer(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      sites_cart=sites_cart_moving)
                  assert rotamer == model_rot

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
        #print >> log, \
        #  "chain %s has conflicting segid values - skipping" % chain_i.id
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

# XXX wrapper for running in Phenix GUI
class _run_determine_ncs_groups (object) :
  def __init__ (self, params, pdb_hierarchy) :
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy

  def __call__ (self, *args, **kwds) :
    return determine_ncs_groups(
      params=self.params,
      pdb_hierarchy=self.pdb_hierarchy)
