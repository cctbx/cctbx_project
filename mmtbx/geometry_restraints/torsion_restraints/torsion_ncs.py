from __future__ import absolute_import, division, print_function
import cctbx.geometry_restraints
from mmtbx.validation import rotalyze
from mmtbx.validation import ramalyze
from mmtbx.validation import analyze_peptides
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from cctbx.array_family import flex
import iotbx.phil
from libtbx.str_utils import make_sub_header
import sys, math
from libtbx.utils import Sorry
from mmtbx.geometry_restraints.torsion_restraints import utils, rotamer_search
from libtbx import Auto
from libtbx.str_utils import line_breaker
from six.moves import zip

# Refactoring notes 11 May 2015
# Torsion NCS restraints should use the same procedure to find NCS groups
# as other NCS restraints. Therefore all search code will be avoided and
# NCS groups should be provided via ncs_obj which is instance of
# cctbx_project.iotbx.ncs.ncs_preprocess.ncs_group_object

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 2.5
   .type = float
   .short_caption = Restraint sigma (degrees)
 limit = 15.0
   .type = float
   .short_caption = Restraint limit (degrees)
 fix_outliers = False
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
 filter_phi_psi_outliers = True
   .type = bool
   .expert_level = 4
 restrain_to_master_chain = False
   .type = bool
   .expert_level = 4
 silence_warnings = False
   .type = bool
   .expert_level = 4
""")

def target(sites_cart_residue, unit_cell, m):
  sites_frac_residue = unit_cell.fractionalize(sites_cart_residue)
  result = 0
  for rsf in sites_frac_residue:
    result += m.eight_point_interpolation(rsf)
  return result

def all_sites_above_sigma_cutoff(sites_cart_residue,
                                 unit_cell,
                                 m,
                                 sigma_cutoff):
  sites_frac_residue = unit_cell.fractionalize(sites_cart_residue)
  for rsf in sites_frac_residue:
    if m.eight_point_interpolation(rsf) < sigma_cutoff:
      return False
  return True

class torsion_ncs(object):
  def __init__(self,
               model,
               fmodel=None,
               params=None,
               selection=None,
               ncs_groups=None,
               alignments=None,
               ncs_dihedral_proxies=None,
               log=None):
    assert model is not None
    if(log is None): log = sys.stdout
    if params is None:
      params = torsion_ncs_params.extract()
    #parameter initialization
    if params.sigma is None or params.sigma < 0:
      raise Sorry("torsion NCS sigma parameter must be >= 0.0")
    self.sigma = params.sigma
    if params.limit is None or params.limit < 0:
      raise Sorry("torsion NCS limit parameter must be >= 0.0")
    # assert ncs_obj is not None
    self.limit = params.limit
    self.selection = selection
    self.model = model
    self.pdb_hierarchy = None
    self.pdb_hierarchy = self.model.get_hierarchy()
    self.ncs_obj = self.model.get_ncs_obj()
    # self.cache = self.model.get_atom_selection_cache()
    assert self.ncs_obj is not None

    #slack is not a user parameter for now
    self.slack = 0.0
    self.filter_phi_psi_outliers = params.filter_phi_psi_outliers
    self.restrain_to_master_chain = params.restrain_to_master_chain
    self.fmodel = fmodel
    self.ncs_restraints_group_list = self.ncs_obj.get_ncs_restraints_group_list()
    self.ncs_groups_selection_string_list = self.ncs_restraints_group_list.get_array_of_str_selections()
    self.log = log
    self.params = params
    self.dp_ncs = None
    self.rotamer_search_manager = None
    self.ncs_dihedral_proxies = ncs_dihedral_proxies
    self.alignments = alignments
    self.sa = SidechainAngles(False)
    #sanity check
    # if self.pdb_hierarchy is not None:
    #   self.pdb_hierarchy.reset_i_seq_if_necessary()
    #   self.cache = self.pdb_hierarchy.atom_selection_cache()
    if self.alignments is None:
      self.find_ncs_groups(pdb_hierarchy=self.pdb_hierarchy)
    if self.pdb_hierarchy is not None:
      self.find_ncs_matches_from_hierarchy(model=self.model)

  def get_alignments(self):
    # This function should use something like common_res_dict already available
    # in ncs_obj or something similar instead of aligning residues itself.
    def get_key(rg):
      resname = rg.atoms()[0].pdb_label_columns()[5:8]
      if resname.upper() == "MSE":
        resname = "MET"
      updated_resname = utils.modernize_rna_resname(resname)
      return rg.atoms()[0].pdb_label_columns()[4:5]+\
          updated_resname+rg.atoms()[0].pdb_label_columns()[8:]+\
          rg.atoms()[0].segid
    alignments = {}
    # for group in self.ncs_groups:
    for i_group, group in enumerate(self.ncs_restraints_group_list):
      for i, isel in enumerate(group.get_iselections_list()):
        for j, jsel in enumerate(group.get_iselections_list()):
          if i <= j:
            continue
          # isel_plus = isel
          # jsel_plus = jsel
          # if self.ncs_obj.exclude_selection is not None:
          #   isel_plus += " and not (%s)" % self.ncs_obj.exclude_selection
          #   jsel_plus += " and not (%s)" % self.ncs_obj.exclude_selection
          # sel_i = self.cache.selection(isel_plus)
          # sel_j = self.cache.selection(jsel_plus)
          h_i = self.pdb_hierarchy.select(isel)
          h_j = self.pdb_hierarchy.select(jsel)
          # chain matching procedure
          # matching_chain_numbers = []
          matching_chains = []
          for ii, i_chain in enumerate(h_i.chains()):
            for jj, j_chain in enumerate(h_j.chains()):
              # print "Checking chains:", i_chain.id, j_chain.id
              # print "  ", i_chain.atoms_size(), j_chain.atoms_size()
              # print "  ", i_chain.is_similar_hierarchy(j_chain)
              # print "  ", i_chain.as_sequence(), j_chain.as_sequence(), i_chain.as_sequence() == j_chain.as_sequence()
              if (i_chain.atoms_size() == j_chain.atoms_size()
                  # and i_chain.is_similar_hierarchy(j_chain)
                  and i_chain.as_sequence() == j_chain.as_sequence()):
                # matching_chain_numbers.append((ii, jj))
                matching_chains.append((i_chain, j_chain))
                break
          # residue matching
          i_chains = h_i.chains()
          j_chains = h_j.chains()
          residue_match_map1 = {}
          residue_match_map2 = {}
          if len(matching_chains) == 0:
            msg = "Failed to find matching chains. "
            msg += "No NCS restraints will be applied.\n"
            msg += "Try to use phenix.simple_ncs_from_pdb or leave ncs_groups "
            msg += "blank to allow \nautomatic search for NCS copies."
            print(msg, file=self.log)
          for ic, jc in matching_chains:
            for rg1, rg2 in zip(ic.residue_groups(), jc.residue_groups()):
              resname1 = rg1.atom_groups()[0].resname
              resname2 = rg2.atom_groups()[0].resname
              if (resname1 != resname2 and
                  # This excludes the case when in one NCS copy residue is MET
                  # and in another - MSE. They will be excluded without
                  # raising Sorry. They could matched, but it is difficult
                  # to figure out in this code how to make it happen.
                  not (resname1 in ["MET", "MSE"] and resname2 in ["MET", "MSE"])):
                msg = "Error in matching procedure: matching "
                msg += "'%s %s' and '%s %s'.\n" % (
                    resname1, rg1.id_str(), resname2, rg2.id_str())
                msg += "Please tell developers about this error"
                msg += " (with files to reproduce it)! "
                msg += "If you need to proceed try to disable NCS."
                raise Sorry(msg)
              residue_match_map1[get_key(rg1)] = get_key(rg2)
              residue_match_map2[get_key(rg2)] = get_key(rg1)
          # This duplication won't be necessary when all the rest code would
          # be able to work without it. 2x less memory...
          ikey = self.ncs_groups_selection_string_list[i_group][i]
          jkey = self.ncs_groups_selection_string_list[i_group][j]
          key1 = (ikey, jkey)
          key2 = (jkey, ikey)
          alignments[key1] = residue_match_map1
          alignments[key2] = residue_match_map2
    return alignments

  def find_ncs_groups(self, pdb_hierarchy):
    print("Determining NCS matches...", file=self.log)
    self.use_segid = False
    chains = pdb_hierarchy.models()[0].chains()
    n_ncs_groups = self.ncs_obj.number_of_ncs_groups
    if n_ncs_groups > 0:
      new_alignments = self.get_alignments()
      self.alignments = new_alignments
      return
    else:
      # ================================================
      # Should never be executed because I'm disabling search functionality
      # in this module
      # ================================================
      assert 0

  def as_cif_block(self, cif_block, hierarchy, scattering_type):
    self.ncs_restraints_group_list.as_cif_block(
        cif_block=cif_block,
        hierarchy=hierarchy,
        scattering_type=scattering_type,
        ncs_type='Torsion NCS')
    return cif_block

  def as_pdb(self, out):
    assert out is not None
    torsion_counts=self.get_number_of_restraints_per_group()
    sites_cart = self.model.get_sites_cart()
    self.get_torsion_rmsd(sites_cart=sites_cart)
    pr = "REMARK   3  "
    print(pr+"TORSION NCS DETAILS.", file=out)
    print(pr+" NUMBER OF NCS GROUPS : %-6d"%len(self.ncs_groups_selection_string_list), file=out)
    for i_group, ncs_group in enumerate(self.ncs_groups_selection_string_list):
      count = 0
      print(pr+" NCS GROUP : %-6d"%(i_group+1), file=out)
      selection_strings = ncs_group
      for selection in selection_strings:
        lines = line_breaker(selection, width=34)
        for i_line, line in enumerate(lines):
          if (i_line == 0):
            print(pr+"   SELECTION          : %s"%line, file=out)
          else:
            print(pr+"                      : %s"%line, file=out)
        count += torsion_counts[selection]
      print(pr+"   RESTRAINED TORSIONS: %-d" % count, file=out)
      if self.torsion_rmsd is not None:
        print(pr+"   BELOW LIMIT RMSD   : %-10.3f" % \
          self.torsion_rmsd, file=out)
      if self.all_torsion_rmsd is not None:
        print(pr+"   ALL RESTRAINT RMSD : %-10.3f" % \
          self.all_torsion_rmsd, file=out)
    if self.histogram_under_limit is not None:
      print(pr + "  Histogram of differences under limit:", file=out)
      self.histogram_under_limit.show(
        f=out,
        prefix=pr+"  ",
        format_cutoffs="%8.3f")
    if self.histogram_over_limit is not None:
      print(pr + "  Histogram of differences over limit:", file=out)
      self.histogram_over_limit.show(
        f=out,
        prefix=pr+"  ",
        format_cutoffs="%8.3f")

  def find_ncs_matches_from_hierarchy(self, model):
    # This is list of dp_match, see below
    self.dp_ncs = []
    self.cb_dp_ncs = []
    self.dihedral_proxies_backup = None
    pdb_hierarchy = model.get_hierarchy()
    self.name_hash = utils.build_name_hash(pdb_hierarchy)
    self.segid_hash = utils.build_segid_hash(pdb_hierarchy)
    self.sym_atom_hash = utils.build_sym_atom_hash(pdb_hierarchy)
    self.sa = SidechainAngles(False)
    self.sidechain_angle_hash = self.build_sidechain_angle_hash()
    self.r = rotalyze.rotalyze(pdb_hierarchy=pdb_hierarchy)
    self.unit_cell = None
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    if self.selection is None:
      self.selection = flex.bool(len(sites_cart), True)

    # XXX Update. Due to possible changes in number of atoms in hierarchy
    # during refinement we have to construct again all internal stuff
    # on every macro-cycle.
    complete_dihedral_proxies = utils.get_complete_dihedral_proxies_2(
        model = self.model)
    # print "Number of complete_dihedral_proxies in find_ncs_matches_from_hierarchy", complete_dihedral_proxies.size()

    if len(self.ncs_restraints_group_list) > 0:
      element_hash = utils.build_element_hash(pdb_hierarchy)
      i_seq_hash = utils.build_i_seq_hash(pdb_hierarchy)
      dp_hash = {}
      for dp in complete_dihedral_proxies:
        h_atom = False
        for i_seq in dp.i_seqs:
          if element_hash[i_seq] == " H":
            h_atom = True
        if not h_atom:
          complete = True
          for i_seq in dp.i_seqs:
            if not self.selection[i_seq]:
              complete = False
          if complete:
            dp_hash[dp.i_seqs] = dp

      super_hash = {}
      res_match_master = {}
      res_to_selection_hash = {}
      # This is for cache/hash generation....
      # print pdb_hierarchy.as_pdb_string()
      for i, group in enumerate(self.ncs_groups_selection_string_list):
        for isel, chain_i in zip(self.ncs_restraints_group_list[i].get_iselections_list(), group):
          # print i, "isel,", list(isel)
          c_atoms = pdb_hierarchy.select(isel).atoms()
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
      resname = None
      atoms_key = None
      for dp in complete_dihedral_proxies:
        temp = dict()
        #filter out unwanted torsions
        atoms = []
        for i_seq in dp.i_seqs:
          atom = self.name_hash[i_seq][:4]
          atoms.append(atom)
          atoms_key = ",".join(atoms)
          resname = self.get_torsion_resname(dp)
        if resname is not None:
          if ( (resname.lower() == 'arg') and
               (atoms_key == ' CD , NE , CZ , NH2') ):
            continue
          elif ( (resname.lower() == 'tyr') and
               (atoms_key == ' CD1, CE1, CZ , OH ' or
                atoms_key == ' CE1, CZ , OH , HH ') ):
            continue
          elif ( (resname.lower() == 'ser') and
               (atoms_key == ' CA , CB , OG , HG ') ):
            continue
          elif ( (resname.lower() == 'thr') and
               (atoms_key == ' CA , CB , OG1, HG1') ):
            continue
          elif ( (resname.lower() == 'cys') and
               (atoms_key == ' CA , CB , SG , HG ') ):
            continue
          elif ( (resname.lower() == 'met') and
               (atoms_key == ' CG , SD , CE ,1HE ' or
                atoms_key == ' CG , SD , CE , HE1') ):
            continue
        ################
        for i_seq in dp.i_seqs:
          cur_matches = super_hash.get(i_seq)
          if cur_matches is None:
            continue
          for key in list(cur_matches.keys()):
            try:
              temp[key].append(cur_matches[key])
            except Exception:
              temp[key] = []
              temp[key].append(cur_matches[key])
        # This is [[dp,bool,bool,bool],[dp, bool,bool,bool],...]
        # where dp - dihedral proxy from complete_dihedral_proxies
        # bool, bool, bool - is phi/psi/omega angle
        dp_match = []
          # cctbx.geometry_restraints.shared_dihedral_proxy()
        dp_match.append([dp, False, False, False])
        for key in list(temp.keys()):
          cur_dp_hash = dp_hash.get(tuple(temp[key]))
          if cur_dp_hash is not None:
            dp_match.append([cur_dp_hash, False, False, False])
            dp_hash[tuple(temp[key])] = None
        dp_hash[dp.i_seqs] = None
        if len(dp_match) > 1:
          self.dp_ncs.append(dp_match)
      #initialize tracking hashes
      for dp_set in self.dp_ncs:
        for dp in dp_set:
          angle_atoms = self.get_torsion_atoms(dp[0])
          angle_resname = self.get_torsion_resname(dp[0])
          angle_id = utils.get_torsion_id(dp=dp[0], name_hash=self.name_hash)
          #phi
          if angle_atoms == ' C  '+' N  '+' CA '+' C  ':
            dp[1] = True
          #psi
          elif angle_atoms == ' N  '+' CA '+' C  '+' N  ':
            dp[2] = True
          #omega
          elif angle_atoms == ' CA '+' C  '+' N  '+' CA ':
            dp[3] = True

      match_counter = {}
      inclusive_range = {}
      for group in self.ncs_groups_selection_string_list:
        cur_len = len(group)
        for chain in group:
          match_counter[chain] = cur_len
          inclusive_range[chain] = []

      matched = []
      ncs_match_hash = {}
      for dp_set in self.dp_ncs:
        if len(dp_set) > 1:
          key_set = []
          for dp in dp_set:
            cur_key = ""
            for i_seq in dp[0].i_seqs:
              cur_key += self.name_hash[i_seq]
            if cur_key[4:19] == cur_key[23:38] and \
               cur_key[4:19] == cur_key[42:57]:
              key_set.append(cur_key[4:19])
          if len(dp_set) == len(key_set): # probably always True
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

      for res in list(self.ncs_match_hash.keys()):
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
      for key in list(inclusive_range.keys()):
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

      if self.ncs_dihedral_proxies is None:
        self.show_ncs_summary(log=self.log)
      if self.ncs_dihedral_proxies is None: #first time run
        print("Initializing torsion NCS restraints...", file=self.log)
      else:
        print("Verifying torsion NCS restraints...", file=self.log)
      self.rama = ramalyze.ramalyze(pdb_hierarchy=pdb_hierarchy)
      self.generate_dihedral_ncs_restraints(
        sites_cart=sites_cart,
        pdb_hierarchy=pdb_hierarchy,
        log=self.log)
    elif(not self.params.silence_warnings):
      # ================================================
      # Should never be executed
      # ================================================
      assert 0
      print("** WARNING: No torsion NCS found!!" + \
        "  Please check parameters. **", file=self.log)

  def show_ncs_summary(self, log=None):
    if(log is None): log = sys.stdout
    def get_key_chain_num(res):
      return res[4:]
    sorted_keys = sorted(self.ncs_match_hash, key=get_key_chain_num)
    print("--------------------------------------------------------", file=log)
    print("Torsion NCS Matching Summary:", file=log)
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
      print(print_line, file=log)
    print("--------------------------------------------------------", file=log)

  def reduce_redundancies(self):
    #clear out redundancies
    for key in list(self.ncs_match_hash.keys()):
      for key2 in list(self.ncs_match_hash.keys()):
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

  def get_chi_id(self, dp):
    atoms = []
    for i_seq in dp.i_seqs:
      atom = self.name_hash[i_seq][:4]
      atoms.append(atom)
    atoms_key = ",".join(atoms)
    resname = self.get_torsion_resname(dp)
    resAtomsToChi = self.sa.resAtomsToChi.get(resname.lower())
    if resAtomsToChi is None:
      chi_id = None
    else:
      chi_id = resAtomsToChi.get(atoms_key)
    return chi_id

  def build_chi_tracker(self, pdb_hierarchy):
    self.current_chi_restraints = {}
    #current_rotamers = self.r.current_rotamers
    current_rotamers = {}
    for rot in self.r.results:
      current_rotamers[rot.id_str()] = rot.rotamer_name
    for key in list(self.ncs_match_hash.keys()):
      search_key = key[4:11]+key[0:4] # SEGIDs?
      rotamer = current_rotamers.get(search_key)
      if rotamer is not None:
        split_rotamer = rotalyze.split_rotamer_names(rotamer=rotamer)
        self.current_chi_restraints[key] = split_rotamer
      key_list = self.ncs_match_hash.get(key)
      for key2 in key_list:
        search_key2 = key2[4:11]+key2[0:4] # SEGIDs?
        rotamer = current_rotamers.get(search_key2)
        if rotamer is not None:
          split_rotamer = rotalyze.split_rotamer_names(rotamer=rotamer)
          self.current_chi_restraints[key2] = split_rotamer

  def generate_dihedral_ncs_restraints(
        self,
        sites_cart,
        pdb_hierarchy,
        log):
    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)
    self.build_chi_tracker(pdb_hierarchy)
    self.ncs_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    target_map_data = None
    #if self.fmodel is not None and self.use_cc_for_target_angles:
    #  target_map_data, residual_map_data = self.prepare_map(
    #                                         fmodel=self.fmodel)
    #rama_outliers = None
    rama_outlier_list = []
    omega_outlier_list = []
    if self.filter_phi_psi_outliers:
      rama_outlier_list = \
        self.get_ramachandran_outliers(pdb_hierarchy)
      #rama_outliers = \
      #  self.get_ramachandran_outliers(pdb_hierarchy)
      #for outlier in rama_outliers.splitlines():
      #  temp = outlier.split(':')
      #  rama_outlier_list.append(temp[0])
      omega_outlier_list = \
        self.get_omega_outliers(pdb_hierarchy)
    torsion_counter = 0

    for dp_set in self.dp_ncs:
      if len(dp_set) < 2:
        continue
      angles = []
      #cc_s = []
      is_rama_outlier = []
      is_omega_outlier = []
      rotamer_state = []
      chi_ids = []
      wrap_hash = {}
      for i, dp in enumerate(dp_set):
        di = cctbx.geometry_restraints.dihedral(
               sites_cart=sites_cart, proxy=dp[0])
        angle = di.angle_model
        wrap_chis = self.is_symmetric_torsion(dp[0])
        if wrap_chis:
          if angle > 90.0 or angle < -90.0:
            sym_i_seq = dp[0].i_seqs[3] #4th atom
            swap_i_seq = self.sym_atom_hash.get(sym_i_seq)
            if swap_i_seq is not None:
              swap_i_seqs = (dp[0].i_seqs[0],
                             dp[0].i_seqs[1],
                             dp[0].i_seqs[2],
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
        if dp[1] or dp[2]:
          angle_id = utils.get_torsion_id(
                       dp=dp[0],
                       name_hash=self.name_hash,
                       phi_psi=True)
          key = angle_id[4:11]+angle_id[0:4] #segid?
          if key in rama_outlier_list:
            rama_out = True
        is_rama_outlier.append(rama_out)

        omega_out = False
        if dp[3]:
          angle_id = utils.get_torsion_id(
                       dp=dp[0],
                       name_hash=self.name_hash,
                       omega=True)
          key1 = \
            angle_id[0][4:6].strip()+angle_id[0][6:10]+' '+angle_id[0][0:4]
          key2 = \
            angle_id[1][4:6].strip()+angle_id[1][6:10]+' '+angle_id[1][0:4]
          if (key1, key2) in omega_outlier_list:
            omega_out = True
        is_omega_outlier.append(omega_out)
        #if target_map_data is not None:
        #  tor_iselection = flex.size_t()
        #  for i_seq in dp.i_seqs:
        #    tor_iselection.append(i_seq)
        #  tor_sites_cart = \
        #    sites_cart.select(tor_iselection)
        #  di_cc = self.get_sites_cc(sites_cart=tor_sites_cart,
        #                            target_map_data=target_map_data)
        #  cc_s.append(di_cc)

        angle_id = utils.get_torsion_id(
                     dp=dp[0],
                     name_hash=self.name_hash,
                     chi_only=True)
        if angle_id is not None:
          split_rotamer_list = self.current_chi_restraints.get(angle_id)
          which_chi = self.get_chi_id(dp[0])
          rotamer_state.append(split_rotamer_list)
          if which_chi is not None:
            chi_ids.append(which_chi)
      target_angles = self.get_target_angles(
                        angles=angles,
                        #cc_s=cc_s,
                        is_rama_outlier=is_rama_outlier,
                        is_omega_outlier=is_omega_outlier,
                        rotamer_state=rotamer_state,
                        chi_ids=chi_ids)
      #if angle_id is not None: # and which_chi is None:
      #print angles
      #print target_angles
      for i, dp in enumerate(dp_set):
        target_angle = target_angles[i]
        angle_atoms = self.get_torsion_atoms(dp[0])
        angle_resname = self.get_torsion_resname(dp[0])
        angle_id = utils.get_torsion_id(dp=dp[0], name_hash=self.name_hash)
        cur_dict = self.sidechain_angle_hash.get(angle_resname)
        angle_name = None
        if cur_dict != None:
          angle_name = \
            cur_dict.get(angle_atoms)
        if target_angle is not None:
          angle_atoms = self.get_torsion_atoms(dp[0])
          angle_resname = self.get_torsion_resname(dp[0])
          angle_id = utils.get_torsion_id(dp=dp[0], name_hash=self.name_hash)
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
              i_seqs=dp[0].i_seqs,
              angle_ideal=target_angle,
              weight=1/self.sigma**2,
              limit=self.limit,
              top_out=TOP_OUT_FLAG,
              slack=self.slack)
          self.ncs_dihedral_proxies.append(dp_add)
          torsion_counter += 1

    if len(self.ncs_dihedral_proxies) == 0:
      if (not self.params.silence_warnings):
        print("** WARNING: No torsion NCS found!!" + \
          "  Please check parameters. **", file=log)
    else:
      print("Number of torsion NCS restraints: %d\n" \
          % len(self.ncs_dihedral_proxies), file=log)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels,
        proxy_label="NCS torsion angle",
        f=sys.stdout):
    self.ncs_dihedral_proxies.show_sorted(
        by_value=by_value,
        sites_cart=sites_cart,
        site_labels=site_labels,
        proxy_label=proxy_label,
        f=f)

  def update_dihedral_ncs_restraints(self,
                                     model,
                                     log=None):
    self.model = model
    if log is None:
      log = sys.stdout
    make_sub_header(
      "Updating torsion NCS restraints",
      out=log)
    self.dp_ncs = None
    if self.dp_ncs is None:
      self.find_ncs_matches_from_hierarchy(model=model)
    else:
      self.generate_dihedral_ncs_restraints(sites_cart=self.model.get_sites_cart(),
                                            pdb_hierarchy=self.model.get_hierarchy(),
                                            log=log)
    # self.add_ncs_dihedral_proxies(geometry=geometry)

  def is_symmetric_torsion(self, dp):
    i_seqs = dp.i_seqs
    resname = self.name_hash[i_seqs[0]][5:8].upper()
    if resname not in \
      ['ASP', 'GLU', 'PHE', 'TYR']: #, 'ASN', 'GLN', 'HIS']:
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
    #elif resname == 'ASN':
    #  if torsion_atoms == [' CA ', ' CB ',' CG ',' OD1'] or \
    #     torsion_atoms == [' CA ', ' CB ',' CG ',' ND2']:
    #    return True
    #elif resname == 'GLN':
    #  if torsion_atoms == [' CB ', ' CG ',' CD ',' OE1'] or \
    #     torsion_atoms == [' CB ', ' CG ',' CD ',' NE2']:
    #    return True
    #elif resname == 'HIS':
    #  if torsion_atoms == [' CA ', ' CB ',' CG ',' ND1'] or \
    #     torsion_atoms == [' CA ', ' CB ',' CG ',' CD2']:
    #    return True
    return False

  def get_target_angles(self,
                        angles,
                        #cc_s,
                        is_rama_outlier,
                        is_omega_outlier,
                        rotamer_state,
                        chi_ids):
    assert (len(rotamer_state) == len(angles)) or \
           (len(rotamer_state) == 0)
    chi_num = None
    # print "chi_ids", chi_ids
    if len(chi_ids) > 0:
      assert len(chi_ids) == chi_ids.count(chi_ids[0])
      if ('oh' not in chi_ids and
          'sh' not in chi_ids and
          'me' not in chi_ids):
        chi_num = chi_ids[0][-1:]
    clusters = {}
    used = []
    target_angles = [None] * len(angles)

    #check for all outliers for current target
    if ( (is_rama_outlier.count(False)  == 0) or
         (is_omega_outlier.count(False) == 0) or
         ( ( (len(rotamer_state)-rotamer_state.count(None)) < 2 and
              len(rotamer_state) > 0)) ):
      for i, target in enumerate(target_angles):
        target_angles[i] = None
      return target_angles
    ###########

    max_i = None
    #for i, cc in enumerate(cc_s):
    #  if is_rama_outlier[i]:
    #    continue
    #  if max_i is None:
    #    max_i = i
    #  elif max < cc:
    #    max_i = i
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
          if len(rotamer_state)> 0:
            if rotamer_state[i] is None:
              continue
            elif rotamer_state[j] is None:
              continue
          nonstandard_chi = False
          is_proline = False
          chi_matching = True
          if len(rotamer_state) > 0:
            if rotamer_state[i][0] in ['UNCLASSIFIED', 'OUTLIER',
                                       'Cg_exo', 'Cg_endo']:
              nonstandard_chi = True
            elif rotamer_state[j][0] in ['UNCLASSIFIED', 'OUTLIER',
                                         'Cg_exo', 'Cg_endo']:
              nonstandard_chi = True
            if rotamer_state[i][0] in ['Cg_exo', 'Cg_endo'] or \
               rotamer_state[j][0] in ['Cg_exo', 'Cg_endo']:
              is_proline = True
          if (chi_num is not None) and not nonstandard_chi:
            chi_counter = int(chi_num)
            while chi_counter > 0:
              if (rotamer_state[i][chi_counter-1] !=
                  rotamer_state[j][chi_counter-1]):
                chi_matching = False
              chi_counter -= 1
          if not chi_matching:
            continue
          if is_proline and \
             rotamer_state[i][0] != rotamer_state[j][0]:
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
    for key in list(clusters.keys()):
      cluster = clusters[key]
      if cluster is None:
        target_angles[key] = None
      else:
        cluster_angles = []
        cluster_outliers = 0
        for i in cluster:
          if is_rama_outlier[i]:
            cluster_angles.append(None)
          elif is_omega_outlier[i]:
            cluster_angles.append(None)
          elif len(rotamer_state) > 0:
            if rotamer_state[i][0] == 'OUTLIER':
              cluster_angles.append(None)
              cluster_outliers += 1
            else:
              cluster_angles.append(angles[i])
          else:
            cluster_angles.append(angles[i])
        if max_i is not None:
          target_angle = angles[max_i]
        else:
          target_angle = utils.get_angle_average(cluster_angles)
        if self.params.target_damping:
          for c in cluster:
            if target_angle is None:
              target_angles[c] = None
            else:
              c_dist = utils.angle_distance(angles[c], target_angle)
              if c_dist > self.params.damping_limit:
                d_target = \
                  utils.get_angle_average([angles[c], target_angle])
                target_angles[c] = d_target
              else:
                target_angles[c] = target_angle
        else:
          if (len(cluster) - cluster_outliers) == 1:
            for c in cluster:
              if rotamer_state[c][0] == 'OUTLIER':
                target_angles[c] = target_angle
              else:
                target_angles[c] = None
          else:
            for c in cluster:
              target_angles[c] = target_angle
        if (self.restrain_to_master_chain):
          if target_angles[cluster[0]] is not None:
            for i,c in enumerate(cluster):
              if i == 0:
                target_angles[c] = None
              else:
                target_angles[c] = angles[cluster[0]]
    return target_angles

  def get_ramachandran_outliers(self, pdb_hierarchy):
    rama_outliers = []
    self.rama = ramalyze.ramalyze(pdb_hierarchy=pdb_hierarchy,
                                  outliers_only=True)
    for r in self.rama.results:
      rama_outliers.append(r.id_str())
    return rama_outliers

  def get_omega_outliers(self, pdb_hierarchy):
    cis_peptides, trans_peptides, omega_outliers = \
      analyze_peptides.analyze(pdb_hierarchy=pdb_hierarchy)
    return omega_outliers

  def get_rotamer_data(self, pdb_hierarchy):
    self.r = rotalyze.rotalyze(pdb_hierarchy=pdb_hierarchy)
    model_hash = {}
    model_score = {}
    all_rotamers = {}
    model_chis = {}
    for rot in self.r.results:
      model_hash[rot.id_str()] = rot.rotamer_name
      model_score[rot.id_str()] = rot.score
    for key in list(self.res_match_master.keys()):
      res_key = key[4:10]+' '+key[0:4]
      all_rotamers[res_key] = []
      model_rot = model_hash.get(res_key)
      if model_rot is not None and model_rot != "OUTLIER":
        all_rotamers[res_key].append(model_rot)
      for match_res in self.res_match_master[key]:
        j_key = match_res[4:10]+' '+match_res[0:4]
        j_rot = model_hash.get(j_key)
        if j_rot is not None and j_rot != "OUTLIER":
          if j_rot not in all_rotamers[res_key]:
            all_rotamers[res_key].append(j_rot)

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = \
              rotalyze.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              #try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = \
                  self.sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = utils.id_str(
                          chain_id=chain.id,
                          resseq=residue_group.resseq,
                          resname=atom_group.resname,
                          icode=residue_group.icode,
                          altloc=atom_group.altloc)
                  #key = '%s%5s %s' % (
                  #    chain.id, residue_group.resid(),
                  #    atom_group.altloc+atom_group.resname)
                  model_chis[key] = chis
    return model_hash, model_score, all_rotamers, model_chis

  def fix_rotamer_outliers(self,
                           xray_structure,
                           geometry_restraints_manager,
                           pdb_hierarchy,
                           outliers_only=False,
                           log=None,
                           quiet=False):
    self.last_round_outlier_fixes = 0
    if self.rotamer_search_manager is None:
      self.rotamer_search_manager = rotamer_search.manager(
                                      pdb_hierarchy=pdb_hierarchy,
                                      xray_structure=xray_structure,
                                      name_hash=self.name_hash,
                                      selection=self.selection,
                                      log=self.log)
    if self.unit_cell is None:
      self.unit_cell = xray_structure.unit_cell()
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

    self.rotamer_search_manager.prepare_map(fmodel=fmodel)

    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)

    fix_list = {}
    rotamer_targets = {}

    for key in list(self.res_match_master.keys()):
      res_key = key[4:11]+key[0:4]
      model_rot = model_hash.get(res_key)
      if model_rot == "OUTLIER":
        rotamer = None
        score = 0.0
        for match_res in self.res_match_master[key]:
          j_key = match_res[4:11]+match_res[0:4]
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
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if not chain.is_protein():
          continue
        for residue_group in chain.residue_groups():
          all_dict = \
            rotalyze.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY"]:
              continue
            key = utils.id_str(
                    chain_id=chain.id,
                    resseq=residue_group.resseq,
                    resname=atom_group.resname,
                    icode=residue_group.icode,
                    altloc=atom_group.altloc)
            #key = '%s%5s %s' % (
            #          chain.id, residue_group.resid(),
            #          atom_group.altloc+atom_group.resname)
            if key in list(fix_list.keys()):
              model_rot, m_chis, value = rotalyze.evaluate_rotamer(
                atom_group=atom_group,
                sidechain_angles=self.sa,
                rotamer_evaluator=self.model.get_rotamer_manager(),
                rotamer_id=self.model.get_rotamer_id(),
                all_dict=all_dict,
                sites_cart=sites_cart_moving)
              residue_name = key[-3:]
              cur_rotamer = rotamer_targets[key]
              r_chis = self.sa.get_rotamer_angles(
                         residue_name=residue_name,
                         rotamer_name=cur_rotamer)
              if m_chis is not None and r_chis is not None:
                status = self.rotamer_search_manager.search(
                  atom_group=atom_group,
                  all_dict=all_dict,
                  m_chis=m_chis,
                  r_chis=r_chis,
                  rotamer=cur_rotamer,
                  sites_cart_moving=sites_cart_moving,
                  xray_structure=xray_structure,
                  key=key)
                if status:
                  print("Set %s to %s rotamer" % \
                    (key, cur_rotamer), file=log)
                  self.last_round_outlier_fixes += 1

  def get_sites_cc(self,
                   sites_cart,
                   target_map_data):
    t = target(sites_cart, self.unit_cell, target_map_data)
    return t


  def get_sidechain_map_correlation(self,
                                    xray_structure,
                                    pdb_hierarchy):
    map_cc_hash = {}
    sigma_cutoff_hash = {}
    fmodel = self.fmodel
    target_map_data, residual_map_data = \
      utils.prepare_map(fmodel=fmodel)
    sites_cart_moving = xray_structure.sites_cart()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        #only works with protein sidechains
        if not chain.is_protein():
          continue
        for residue_group in chain.residue_groups():
          all_dict = rotalyze.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY"]:
              continue
            key = atom_group.atoms()[0].pdb_label_columns()[4:]+\
                  atom_group.atoms()[0].segid
            residue_iselection = atom_group.atoms().extract_i_seq()
            residue_elements = atom_group.atoms().extract_element()
            sidechain_only_iselection = flex.size_t()
            for i, i_seq in enumerate(residue_iselection):
              atom_name = self.name_hash[i_seq][0:4]
              if atom_name not in [' N  ', ' CA ', ' C  ', ' O  '] and \
                 residue_elements[i].strip() not in ['H','D']:
                sidechain_only_iselection.append(i_seq)
            sites_cart_residue = \
              sites_cart_moving.select(sidechain_only_iselection)
            t_test = self.get_sites_cc(sites_cart_residue,
                                       target_map_data)
            map_cc_hash[key] = t_test
            sigma_state = all_sites_above_sigma_cutoff(
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
    # Removed check for existence to prevent a bug. See commit
    # 67af93785e4 for more info.
    self.rotamer_search_manager = rotamer_search.manager(
                                    pdb_hierarchy=pdb_hierarchy,
                                    xray_structure=xray_structure,
                                    name_hash=self.name_hash,
                                    selection=self.selection,
                                    log=self.log)
    if self.unit_cell is None:
      self.unit_cell = xray_structure.unit_cell()
    sites_cart = xray_structure.sites_cart()
    for atom in pdb_hierarchy.atoms():
      i_seq = atom.i_seq
      atom.xyz = sites_cart[i_seq]
    fmodel = self.fmodel
    if(log is None): log = self.log
    make_sub_header(
      "Checking NCS rotamer consistency",
      out=log)

    self.rotamer_search_manager.prepare_map(fmodel=fmodel)

    model_hash, model_score, all_rotamers, model_chis = \
      self.get_rotamer_data(pdb_hierarchy=pdb_hierarchy)

    sites_cart_moving = xray_structure.sites_cart()
    map_cc_hash, sigma_cutoff_hash = \
      self.get_sidechain_map_correlation(xray_structure, pdb_hierarchy)
    cc_candidate_list = []
    for key in list(self.ncs_match_hash.keys()):
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
          if not sigma_cutoff_hash[cur_key]:
            cc_candidate_list.append(cur_key)

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if not chain.is_protein():
          continue
        for residue_group in chain.residue_groups():
          all_dict = rotalyze.construct_complete_sidechain(residue_group)
          for atom_group in residue_group.atom_groups():
            if atom_group.resname in ["PRO", "GLY"]:
              continue
            key = utils.id_str(
                    chain_id=chain.id,
                    resseq=residue_group.resseq,
                    resname=atom_group.resname,
                    icode=residue_group.icode,
                    altloc=atom_group.altloc)
            #key = '%s%5s %s' % (
            #          chain.id, residue_group.resid(),
            #          atom_group.altloc+atom_group.resname)
            if key in all_rotamers:
              if (len(all_rotamers[key]) >= 2):
                cc_key = atom_group.atoms()[0].pdb_label_columns()[4:]+\
                  atom_group.atoms()[0].segid
                if cc_key not in cc_candidate_list:
                  continue
                model_rot, m_chis, value = rotalyze.evaluate_rotamer(
                  atom_group=atom_group,
                  sidechain_angles=self.sa,
                  rotamer_evaluator=self.model.get_rotamer_manager(),
                  rotamer_id=self.model.get_rotamer_id(),
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
                  r_chis = self.sa.get_rotamer_angles(
                             residue_name=residue_name,
                             rotamer_name=cur_rotamer)
                  if m_chis is not None and r_chis is not None:
                    status = self.rotamer_search_manager.search(
                      atom_group=atom_group,
                      all_dict=all_dict,
                      m_chis=m_chis,
                      r_chis=r_chis,
                      rotamer=cur_rotamer,
                      sites_cart_moving=sites_cart_moving,
                      xray_structure=xray_structure,
                      key=key)
                    if status:
                      current_best = cur_rotamer
                      atom_dict = all_dict.get(atom_group.altloc)
                      m_chis = \
                        self.sa.measureChiAngles(atom_group,
                                                 atom_dict,
                                                 sites_cart_moving)
                if current_best != model_rot:
                  print("Set %s to %s rotamer" % \
                    (key,
                     current_best), file=self.log)
                  self.last_round_rotamer_changes += 1
                else:
                  rotamer, chis, value = rotalyze.evaluate_rotamer(
                    atom_group=atom_group,
                    sidechain_angles=self.sa,
                    rotamer_evaluator=self.model.get_rotamer_manager(),
                    rotamer_id=self.model.get_rotamer_id(),
                    all_dict=all_dict,
                    sites_cart=sites_cart_moving)
                  assert rotamer == model_rot

  def build_sidechain_angle_hash(self):
    sidechain_angle_hash = {}
    for key in list(self.sa.atomsForAngle.keys()):
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

  def get_number_of_restraints_per_group(self):
    torsion_counts = {}
    for i_gr, group in enumerate(self.ncs_restraints_group_list):
      for i_sel, selection in enumerate(group.get_iselections_list()):
        key = self.ncs_groups_selection_string_list[i_gr][i_sel]
        torsion_counts[key] = self.ncs_dihedral_proxies.\
            proxy_select(n_seq=self.model.get_number_of_atoms(),
                         iselection=selection).\
                size()
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
         (self.limit < 180.0) ):
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

  def proxy_select(self, nseq, iselection):
    #
    # This is still not proper manager selection. A lot of stuff remains old.
    # It still works because it is being updated every macro-cycle via
    # update_dihedral_ncs_restraints
    #
    import copy
    assert (self.ncs_dihedral_proxies is not None)
    new_ncs_dihedral_proxies = \
        self.ncs_dihedral_proxies.proxy_select(nseq, iselection)
    new_manager = copy.copy(self)
    new_manager.ncs_dihedral_proxies = new_ncs_dihedral_proxies
    new_manager.ncs_restraints_group_list = \
        self.ncs_restraints_group_list.select(flex.bool(nseq, iselection))
    return new_manager

  def remove_reference_dihedrals_in_place(self):
    self.ncs_dihedral_proxies = None

  def get_n_proxies(self):
    if self.ncs_dihedral_proxies is not None:
      return self.ncs_dihedral_proxies.size()
    else:
      return 0

  def size(self):
    return self.get_n_proxies()

  def target_and_gradients(self, sites_cart, unit_cell, gradient_array):
    if unit_cell is None:
      return cctbx.geometry_restraints.dihedral_residual_sum(
        sites_cart=sites_cart,
        proxies=self.ncs_dihedral_proxies,
        gradient_array=gradient_array)
    else:
      return cctbx.geometry_restraints.dihedral_residual_sum(
          unit_cell=unit_cell,
          sites_cart=sites_cart,
          proxies=self.ncs_dihedral_proxies,
          gradient_array=gradient_array)

# XXX wrapper for running in Phenix GUI
class _run_iotbx_ncs_input(object):
  def __init__(self, params, pdb_hierarchy):
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy

  def __call__(self, *args, **kwds):
    return iotbx.ncs.input(hierarchy=self.pdb_hierarchy).\
      print_ncs_phil_param()
