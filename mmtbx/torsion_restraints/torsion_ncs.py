import mmtbx.alignment
from iotbx.pdb import amino_acid_codes
import cctbx.geometry_restraints
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.refinement import fit_rotamers
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
import sys
from libtbx import group_args
from mmtbx.ncs import restraints
from libtbx.utils import Sorry
from mmtbx.torsion_restraints import utils
from mmtbx import ncs

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 1.0
   .type = float
   .short_caption = Restraint sigma (degrees)
 limit = 15.0
   .type = float
   .short_caption = Restraint limit (degrees)
 cutoff = 45.0
   .type = float
   .short_caption = Cutoff (degrees)
 slack = 0.0
   .type = float
   .short_caption = Restraint slack (degrees)
 similarity = .80
   .type = float
   .short_caption = Sequence similarity cutoff
 b_factor_weight = None
   .type=float
   .short_caption = B factor weight
   .expert_level=1
 hydrogens = False
   .type = bool
   .short_caption = Include hydrogens
 main_chain = True
   .type = bool
   .short_caption = Include backbone atoms
 side_chain = True
   .type = bool
   .short_caption = Include sidechain atoms
 fix_outliers = True
   .type = bool
   .short_caption = Fix rotamer outliers first
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
 }
""")

class torsion_ncs(object):
  def __init__(self,
               pdb_hierarchy,
               geometry,
               sites_cart,
               params,
               log=None):
    if(log is None): log = sys.stdout
    self.sigma = params.sigma
    self.limit = params.limit
    self.slack = params.slack
    self.pdb_hierarchy = pdb_hierarchy
    self.ncs_groups = []
    self.dp_ncs = []
    self.sequences = {}
    self.padded_sequences = {}
    self.structures = {}
    self.residue_start = {}
    self.residue_finish = {}
    self.offset_dict = {}
    self.ncs_dihedral_proxies = None
    self.name_hash = utils.build_name_hash(pdb_hierarchy)
    self.segid_hash = utils.build_segid_hash(pdb_hierarchy)
    self.params = params
    self.found_ncs = None
    self.log = log
    self.njump = 1
    self.min_length = 10
    self.min_percent= 95.0
    print >> self.log, "Determining NCS matches..."
    pair_hash = {}
    dp_hash = {}
    used_chains = []
    res_match_hash = {}
    atom_labels = list(self.pdb_hierarchy.atoms_with_labels())
    segids = flex.std_string([ a.segid for a in atom_labels ])
    use_segid = not segids.all_eq('    ')
    i_seq_hash = utils.build_i_seq_hash(pdb_hierarchy)
    chain_hash = utils.build_chain_hash(pdb_hierarchy)
    name_hash = utils.build_name_hash(pdb_hierarchy)
    element_hash = utils.build_element_hash(pdb_hierarchy)
    chains = pdb_hierarchy.models()[0].chains()
    sel_cache = pdb_hierarchy.atom_selection_cache()
    alignments = {}
    n_ncs_groups = 0
    for i_seq, group in enumerate(self.params.restraint_group):
      n_selections = 0
      for selection in group.selection:
        if(selection is not None):
          n_selections += 1
      if n_selections == 1:
        raise Sorry("Torsion NCS restraint_groups require at least 2 selections")
      elif n_selections > 1:
        n_ncs_groups += 1
    if n_ncs_groups > 0:
      for restraint_group in params.restraint_group:
        for selection_i in restraint_group.selection:
          sel_atoms_i = (utils.phil_atom_selections_as_i_seqs_multiple(
                           cache=sel_cache,
                           string_list=[selection_i]))
          sel_seq, sel_seq_padded, sel_structures = \
            self.extract_sequence_and_sites(
            pdb_hierarchy=pdb_hierarchy,
            selection=sel_atoms_i)
          self.sequences[selection_i] = sel_seq
          self.padded_sequences[selection_i] = sel_seq_padded
          self.structures[selection_i] = sel_structures
      for restraint_group in params.restraint_group:
        ncs_set = []
        for selection_i in restraint_group.selection:
          ncs_set.append(selection_i)
          for selection_j in restraint_group.selection:
            if selection_i == selection_j:
              continue
            selections = (selection_i, selection_j)
            residue_match_map = self._alignment(
                                  pdb_hierarchy=pdb_hierarchy,
                                  params=params,
                                  selections=selections,
                                  log=log)
            key = (selection_i, selection_j)
            alignments[key] = residue_match_map
        self.ncs_groups.append(ncs_set)
    else:
      for i, chain_i in enumerate(chains):
        found_conformer = False
        start = 0
        finish = 0
        for conformer in chain_i.conformers():
          if not conformer.is_protein() and not conformer.is_na():
            continue
          else:
            found_conformer = True
            start = conformer.residues()[0].resseq_as_int()
            finish = \
              conformer.residues()[len(conformer.residues())-1].resseq_as_int()
        if not found_conformer:
          continue
        #test for unique segid
        segid = utils.get_unique_segid(chain_i)
        if segid == None:
          print >> log, \
            "chain %s has conflicting segid values - skipping" % chain_i.id
          continue
        if (use_segid) :
          chain_i_str = "chain '%s' and segid '%s'" % \
            (chain_i.id, segid)
        else :
          chain_i_str = "chain '%s'" % chain_i.id

        chain_i_list = [chain_i_str]
        sel_atoms_i = (utils.phil_atom_selections_as_i_seqs_multiple(
                     cache=sel_cache,
                     string_list=chain_i_list))
        chain_seq, chain_seq_padded, chain_structures = \
          self.extract_sequence_and_sites(
          pdb_hierarchy=pdb_hierarchy,
          selection=sel_atoms_i)
        self.sequences[chain_i_str] = chain_seq
        self.padded_sequences[chain_i_str] = chain_seq_padded
        self.structures[chain_i_str] = chain_structures
        self.residue_start[chain_i_str] = start
        self.residue_finish[chain_i_str] = finish

      for i, chain_i in enumerate(chains):
        found_conformer = False
        for conformer in chain_i.conformers():
          if not conformer.is_protein() and not conformer.is_na():
            continue
          else:
            found_conformer = True
        if not found_conformer:
          continue
        #test for unique segid
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
          #test for unique segid
          segid_j = utils.get_unique_segid(chain_j)
          if segid_j == None:
            continue
          if (use_segid) :
            chain_j_str = "chain '%s' and segid '%s'" % (chain_j.id, segid_j)
          else :
            chain_j_str = "chain '%s'" % chain_j.id
          selections = (chain_i_str, chain_j_str)
          residue_match_map = self._alignment(pdb_hierarchy=pdb_hierarchy,
                                params=params,
                                selections=selections,
                                log=log)
          if ( min(len(residue_match_map),
                   chain_i.residue_groups_size(),
                   chain_j.residue_groups_size()) \
               / max(len(residue_match_map),
                     chain_i.residue_groups_size(),
                     chain_j.residue_groups_size()) \
               >= self.params.similarity ):
            key = (chain_i_str, chain_j_str)
            alignments[key] = residue_match_map
            pair_key = (chain_i.id, segid_i)
            match_key = (chain_j.id, segid_j)
            if used_chains is not None:
              if match_key in used_chains:
                continue
            if (not pair_key in pair_hash) :
              pair_hash[pair_key] = []
            pair_hash[pair_key].append(match_key)
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
        self.ncs_groups.append(ncs_set)
      #print self.ncs_groups
      #STOP()

      #calculate sequence offsets
      for ncs_set in self.ncs_groups:
        chain1 = ncs_set[0]
        seq1 = self.sequences[chain1]
        start1 = self.residue_start[chain1]
        for chain2 in ncs_set[1:]:
          seq2 = self.sequences[chain2]
          start2 = self.residue_start[chain2]
          offset = self.find_offset(
                     seq2,
                     start2,
                     seq1,
                     start1)
          self.offset_dict[chain1+chain2] = offset
      #print self.offset_dict
      #STOP()
      new_ncs_groups = "refinement {\n ncs {\n  torsion {\n"
      for ncs_set in self.ncs_groups:
        new_ncs_groups += "   restraint_group {\n"
        for chain in ncs_set:
          new_ncs_groups += "    selection = %s\n" % chain
        if params.b_factor_weight is not None:
          new_ncs_groups += \
            "    b_factor_weight = %f\n" % params.b_factor_weight
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
            cur_align = alignments.get(key)
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
    self.res_match_master = res_match_master

    #symmetric residues - Val, Phe, Leu, Tyr, Asp, Glu

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
      dp_hash[dp.i_seqs] = None
      if len(dp_match) > 1:
        self.dp_ncs.append(dp_match)

    if self.params.verbose:
      self.show_ncs_summary(log=log)
    print >> self.log, "Initializing torsion NCS restraints..."
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)

  def extract_padded_sequence_from_chain(self, chain):
    seq = []
    padded_seq = []
    last_resseq = 0
    is_na = False
    for conformer in chain.conformers():
      if conformer.is_na():
        is_na = True
    for rg in chain.residue_groups():
      resseq = rg.resseq_as_int()
      if (resseq > (last_resseq + 1)) :
        for x in range(resseq - last_resseq - 1) :
          padded_seq.append('X')
      last_resseq = resseq
      resname = rg.unique_resnames()[0]
      if is_na:
        olc = self.get_nucleic_acid_one_letter_code(resname)
      else:
        olc=\
        amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
      padded_seq.append(olc)
    return "".join(padded_seq)

  def extract_sequence_and_sites(self, pdb_hierarchy, selection):
    seq = []
    result = []
    padded_seq = []
    last_resseq = 0
    counter = 0
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        is_na = False
        for conformer in chain.conformers():
          if conformer.is_na():
            is_na = True
        for rg in chain.residue_groups():
          if(len(rg.unique_resnames())==1):
            resname = rg.unique_resnames()[0]
            if is_na:
              olc = self.get_nucleic_acid_one_letter_code(resname)
            else:
              olc= \
              amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
            atoms = rg.atoms()
            i_seqs = utils.get_i_seqs(atoms)
            if(olc!="X") and utils.is_residue_in_selection(i_seqs, selection):
              seq.append(olc)
              resseq = rg.resseq_as_int()
              if (resseq > (last_resseq + 1)) :
                for x in range(resseq - last_resseq - 1) :
                  padded_seq.append('X')
              last_resseq = resseq
              result.append(group_args(i_seq = counter, rg = rg))
              padded_seq.append(olc)
              counter += 1
    return "".join(seq), "".join(padded_seq), result

  def get_nucleic_acid_one_letter_code(self, resname):
    olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
    if olc != "X":
      return "X"
    if resname[0:2] == "  ":
      return resname[2]
    elif resname[0] == " " and (resname[1] == "D" or resname[1] == "d"):
      return resname[2]
    else:
      return resname[0]

  def _alignment(self, pdb_hierarchy,
                    params,
                    selections,
                    log=None):
    if(log is None): log = sys.stdout
    res_match_hash = {}
    model_mseq_res_hash = {}
    ref_mseq_res_hash = {}
    model_seq = self.sequences[selections[0]]
    model_seq_padded = self.padded_sequences[selections[0]]
    model_structures = self.structures[selections[0]]
    ref_seq = self.sequences[selections[1]]
    ref_seq_padded = self.padded_sequences[selections[1]]
    ref_structures = self.structures[selections[1]]
    for struct in model_structures:
      model_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    for struct in ref_structures:
      ref_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    if model_seq == ref_seq:
      pg = mmtbx.alignment.pairwise_global(
             model_seq,
             ref_seq)
    else:
      pg = mmtbx.alignment.pairwise_global(
             model_seq_padded,
             ref_seq_padded)
    offset_i = 0
    offset_j = 0
    i = 0
    j = 0
    seq_j = pg.result2[j]
    for seq_i in pg.result1:
      seq_j = pg.result2[j]
      if seq_i == seq_j and seq_i != 'X' and seq_j != 'X':
        res_match_hash[model_mseq_res_hash[i-offset_i]] = \
          ref_mseq_res_hash[j-offset_j]
        i += 1
        j += 1
      else:
        if seq_i == 'X' and seq_j == 'X':
          i += 1
          j += 1
          offset_i += 1
          offset_j += 1
        elif (seq_i == 'X' and seq_j == '-') or \
             (seq_i == '-' and seq_j == 'X'):
          i += 1
          j += 1
          offset_i += 1
          offset_j += 1
        elif seq_i == 'X':
          i += 1
          j += 1
          offset_i += 1
        elif seq_j == 'X':
          i += 1
          j += 1
          offset_j += 1
        elif seq_i == '-':
          i += 1
          j += 1
          offset_i += 1
        elif seq_j == '-':
          i += 1
          j += 1
          offset_j += 1
        else:
          i += 1
          j += 1
    return res_match_hash

  def show_ncs_summary(self, log=None):
    ncs_match_hash = {}
    matched = []
    if(log is None): log = sys.stdout
    for dp_set in self.dp_ncs:
      key_set = []
      for dp in dp_set:
        if len(dp_set) < 2:
          continue
        cur_key = ""
        for i_seq in dp.i_seqs:
          cur_key += (self.name_hash[i_seq] + self.segid_hash[i_seq])
        if cur_key[5:19] == cur_key[24:38] and \
           cur_key[5:19] == cur_key[43:57]:
          key_set.append(cur_key[5:19])
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
    def get_key_chain_num(res):
      return res[4:]
    sorted_keys = sorted(self.ncs_match_hash, key=get_key_chain_num)
    print >> log, "--------------------------------------------------------"
    print >> log, "Torsion NCS Matching Summary:"
    for key in sorted_keys:
      print_line = key
      for match in self.ncs_match_hash[key]:
        print_line += "  <=====>  %s" % (match)
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

  def generate_dihedral_ncs_restraints(self, sites_cart, log):
    self.ncs_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    for dp_set in self.dp_ncs:
      if len(dp_set) < 2:
        continue
      angles = []
      temp_match = {}
      for dp in dp_set:
        di = cctbx.geometry_restraints.dihedral(sites_cart=sites_cart, proxy=dp)
        angle = di.angle_model
        angles.append(angle)
        temp_match[dp.i_seqs] = angle
      target_angles = self.get_target_angles(angles=angles)
      for dp in dp_set:
        target_angle = target_angles[temp_match[dp.i_seqs]]
        if target_angle is not None:
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
    print >> log, "Updating torsion NCS restraints..."
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)
    self.add_ncs_dihedral_proxies(geometry=geometry)

  def get_target_angles(self, angles):
    clusters = {}
    used = []
    target_angles = {}
    for i in angles:
      if i in used:
        continue
      for j in angles:
        if i == j:
          continue
        elif j in used:
          continue
        else:
          if utils.angle_distance(i, j) <= self.params.cutoff:
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
        target_angle = utils.get_angle_average(cluster)
        if self.params.target_damping:
          for c in cluster:
            c_dist = utils.angle_distance(c, target_angle)
            if c_dist > self.params.damping_limit:
              d_target = utils.get_angle_average([c, target_angle])
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
                           log=None,
                           quiet=False):
    pdb_hierarchy=self.pdb_hierarchy
    if(log is None): log = self.log
    make_sub_header(
      "Correcting NCS rotamer outliers",
      out=log)
    r = rotalyze()
    sa = SidechainAngles(False)
    mon_lib_srv = mmtbx.monomer_library.server.server()
    rot_list_model, coot_model = r.analyze_pdb(hierarchy=pdb_hierarchy)
    model_hash = {}
    model_score = {}
    model_chis = {}
    fix_list = {}
    for line in rot_list_model.splitlines():
      res, rotamericity, chi1, chi2, chi3, chi4, name = line.split(':')
      model_hash[res]=name
      model_score[res]=rotamericity
    for key in self.res_match_master.keys():
      res_key = key[5:6]+key[-5:-1]+' '+key[0:4]
      model_rot = model_hash.get(res_key)
      if model_rot == "OUTLIER":
        rotamer = None
        score = 0.0
        for match_res in self.res_match_master[key]:
          j_key = match_res[5:6]+match_res[-5:-1]+' '+match_res[0:4]
          j_rot = model_hash.get(j_key)
          j_score = model_score.get(j_key)
          if j_rot is not None and j_score is not None:
            if j_rot != "OUTLIER":
              if rotamer == None:
                rotamer = j_key
                score = j_score
              else:
                if j_score > score:
                  rotamer = j_key
                  score = j_score
        if rotamer != None:
          fix_list[res_key] = rotamer

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = sa.measureChiAngles(atom_group, atom_dict)
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

    sites_cart_start = xray_structure.sites_cart()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
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
                assert len(m_chis) == len(axis_and_atoms_to_rotate)
                counter = 0
                residue_iselection = atom_group.atoms().extract_i_seq()
                sites_cart_residue = \
                  xray_structure.sites_cart().select(residue_iselection)
                for aa in axis_and_atoms_to_rotate:
                  axis = aa[0]
                  atoms = aa[1]
                  atom_group.atoms().set_xyz(new_xyz=sites_cart_residue)
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
                  sites_cart_start = sites_cart_start.set_selected(
                        residue_iselection, sites_cart_residue)
                  counter += 1
                xray_structure.set_sites_cart(sites_cart_start)
                print >> log, "Fixed %s rotamer" % key

  def process_ncs_restraint_groups(self, model, processed_pdb_file):
    log = self.log
    ncs_groups = ncs.restraints.groups()
    sites_cart = None

    for param_group in self.params.restraint_group:
      master = param_group.selection[0]
      selection_strings = []
      found_offset = False
      for selection in param_group.selection[1:]:
        offset = self.offset_dict.get(master+selection)
        if offset is None or offset == 0:
          selection_strings.append(selection)
        else:
          start = self.residue_start[selection]
          finish = self.residue_finish[selection]
          temp_selection = selection + " and (resseq %d:%d)" \
                             % (start, finish)
          selection_strings.append(temp_selection)
          found_offset = True
      if found_offset:
        start = self.residue_start[master]
        finish = self.residue_finish[master]
        master = master + " and (resseq %d:%d)" \
                   % (start, finish)

      group = ncs.restraints.group.from_atom_selections(
        processed_pdb              = processed_pdb_file,
        reference_selection_string = master,
        selection_strings          = selection_strings,
        coordinate_sigma           = 0.05,
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

  def find_offset(self,chain1,start1,chain2,start2):
    best_overlap=0
    best_offset=None
    for offset_a in xrange(-len(chain2)+1,len(chain1)+1):
     offset=start1-start2+offset_a*self.njump
     overlap=self.get_overlap(chain1,start1,chain2,start2,offset)
     if overlap and overlap>best_overlap:
       best_overlap=overlap
       best_offset=offset
    return best_offset

  def get_overlap(self,chain1,start1,chain2,start2,offset):
    offset1=0
    offset2=offset
    n_match=0
    mismatches=0
    for i in xrange(len(chain1)):
      resno=i*self.njump+start1+offset1
      res1=chain1[i]
      j=(resno-start2-offset2)//self.njump
      if j>=0 and j<len(chain2):
        res2=chain2[j]
      else:
        res2=None
      if res1!=None and res2!=None:
        if res1==res2:
          n_match+=1
        else:
          mismatches+=1
    if 100.*float(n_match)/float(self.max(1,n_match+mismatches)) < self.min_percent:
      return None
    if n_match >= self.min_length//self.njump:
      return n_match
    else:
      return None

  def max(self,x,y):
    if x>=y: return x
    return y

#split out functions
def determine_ncs_groups(pdb_hierarchy,
                         params=None,
                         log=None):
  pdb_hierarchy.reset_i_seq_if_necessary()
  sel_cache = pdb_hierarchy.atom_selection_cache()
  if params is None:
    params = torsion_ncs_params.extract()
  if log is None:
    log = sys.stdout
  ncs_groups = []
  alignments = {}
  used_chains = []
  sequences = {}
  padded_sequences = {}
  structures = {}
  residue_start = {}
  residue_finish = {}
  pair_hash = {}
  chains = pdb_hierarchy.models()[0].chains()
  for i, chain_i in enumerate(chains):
    found_conformer = False
    start = 0
    finish = 0
    for conformer in chain_i.conformers():
      if not conformer.is_protein() and not conformer.is_na():
        continue
      else:
        found_conformer = True
        start = conformer.residues()[0].resseq_as_int()
        finish = \
          conformer.residues()[len(conformer.residues())-1].resseq_as_int()
    if not found_conformer:
      continue
    chain_i_str = "chain '%s'" % chain_i.id

    chain_i_list = [chain_i_str]
    sel_atoms_i = (utils.phil_atom_selections_as_i_seqs_multiple(
                 cache=sel_cache,
                 string_list=chain_i_list))
    chain_seq, chain_seq_padded, chain_structures = \
      utils.extract_sequence_and_sites(
        pdb_hierarchy=pdb_hierarchy,
        selection=sel_atoms_i)
    sequences[chain_i_str] = chain_seq
    padded_sequences[chain_i_str] = chain_seq_padded
    structures[chain_i_str] = chain_structures
    residue_start[chain_i_str] = start
    residue_finish[chain_i_str] = finish

  for i, chain_i in enumerate(chains):
    found_conformer = False
    for conformer in chain_i.conformers():
      if not conformer.is_protein() and not conformer.is_na():
        continue
      else:
        found_conformer = True
    if not found_conformer:
      continue
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
      chain_j_str = "chain '%s'" % chain_j.id
      selections = (chain_i_str, chain_j_str)
      seq_pair = (sequences[chain_i_str], sequences[chain_j_str])
      seq_pair_padded = (padded_sequences[chain_i_str],
                         padded_sequences[chain_j_str])
      struct_pair = (structures[chain_i_str], structures[chain_j_str])
      residue_match_map = _alignment(pdb_hierarchy=pdb_hierarchy,
                            params=params,
                            selections=selections,
                            sequences=seq_pair,
                            padded_sequences=seq_pair_padded,
                            structures=struct_pair,
                            log=log)
      if ( min(len(residue_match_map),
               chain_i.residue_groups_size(),
               chain_j.residue_groups_size()) \
           / max(len(residue_match_map),
                 chain_i.residue_groups_size(),
                 chain_j.residue_groups_size()) \
           >= params.similarity ):
        key = (chain_i_str, chain_j_str)
        alignments[key] = residue_match_map
        if used_chains is not None:
          if chain_i.id in used_chains:
            continue
        try:
          pair_hash[chain_i.id].append(chain_j.id)
        except Exception:
          pair_hash[chain_i.id] = []
          pair_hash[chain_i.id].append(chain_j.id)
        used_chains.append(chain_j.id)

  for key in pair_hash.keys():
    ncs_set = []
    chain_str = "chain '%s'" % key
    ncs_set.append(chain_str)
    for add_chain in pair_hash[key]:
      chain_str = "chain '%s'" % add_chain
      ncs_set.append(chain_str)
    ncs_groups.append(ncs_set)

  return ncs_groups

def _alignment(pdb_hierarchy,
               params,
               selections,
               sequences,
               padded_sequences,
               structures,
               log=None):
  if(log is None): log = sys.stdout
  res_match_hash = {}
  model_mseq_res_hash = {}
  ref_mseq_res_hash = {}
  model_seq = sequences[0]
  model_seq_padded = padded_sequences[0]
  model_structures = structures[0]
  ref_seq = sequences[1]
  ref_seq_padded = padded_sequences[1]
  ref_structures = structures[1]
  for struct in model_structures:
    model_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
  for struct in ref_structures:
    ref_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
  if model_seq == ref_seq:
    pg = mmtbx.alignment.pairwise_global(
           model_seq,
           ref_seq)
  else:
    pg = mmtbx.alignment.pairwise_global(
           model_seq_padded,
           ref_seq_padded)
  offset_i = 0
  offset_j = 0
  i = 0
  j = 0
  seq_j = pg.result2[j]
  for seq_i in pg.result1:
    seq_j = pg.result2[j]
    if seq_i == seq_j and seq_i != 'X' and seq_j != 'X':
      res_match_hash[model_mseq_res_hash[i-offset_i]] = \
        ref_mseq_res_hash[j-offset_j]
      i += 1
      j += 1
    else:
      if seq_i == 'X' and seq_j == 'X':
        i += 1
        j += 1
        offset_i += 1
        offset_j += 1
      elif (seq_i == 'X' and seq_j == '-') or \
           (seq_i == '-' and seq_j == 'X'):
        i += 1
        j += 1
        offset_i += 1
        offset_j += 1
      elif seq_i == 'X':
        i += 1
        j += 1
        offset_i += 1
      elif seq_j == 'X':
        i += 1
        j += 1
        offset_j += 1
      elif seq_i == '-':
        i += 1
        j += 1
        offset_i += 1
      elif seq_j == '-':
        i += 1
        j += 1
        offset_j += 1
      else:
        i += 1
        j += 1
  return res_match_hash

# XXX wrapper for running in Phenix GUI
class _run_determine_ncs_groups (object) :
  def __init__ (self, params, pdb_hierarchy) :
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy

  def __call__ (self, *args, **kwds) :
    return determine_ncs_groups(
      params=self.params,
      pdb_hierarchy=self.pdb_hierarchy)
