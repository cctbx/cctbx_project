import mmtbx.alignment
from iotbx.pdb import amino_acid_codes
import cctbx.geometry_restraints
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.refinement import fit_rotamers
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
from iotbx.pdb import common_residue_names_get_class
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
import sys, math
from libtbx import Auto
from libtbx import group_args

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 1.0
   .type = float
 limit = 15.0
   .type = float
 cutoff = 45.0
   .type = float
 slack = 0.0
   .type = float
 similarity = .80
   .type = float
 hydrogens = False
   .type = bool
 main_chain = True
   .type = bool
 side_chain = True
   .type = bool
 fix_outliers = True
   .type = bool
 target_damping = False
   .type = bool
 damping_limit = 10.0
   .type = float
 verbose = False
   .type = bool
 edits
   .short_caption = Edit torsion NCS restraints
   .style = menu_item parent_submenu:advanced noauto
   .expert_level = 2
 {
   include scope \
     mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
 }
 remove
   .short_caption = Remove torsion NCS restraints
   .expert_level = 2
   .style = menu_item parent_submenu:advanced noauto
 {
   include scope \
   mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str
 }
 alignment
    .help = Set of parameters for sequence alignment. Defaults are good for most \
            of cases
    .short_caption = Sequence alignment
    .style = box auto_align
 {
  alignment_style =  *local global
    .type = choice
  gap_opening_penalty = 1
    .type = float
  gap_extension_penalty = 1
    .type = float
  similarity_matrix =  blosum50  dayhoff *identity
    .type = choice
 }
 reference_group
  .multiple=True
  .optional=True
  .short_caption=Torsion NCS restraint group
  .style = noauto auto_align menu_item parent_submenu:advanced
 {
  reference=None
    .type=atom_selection
    .short_caption=Reference selection
  selection=None
    .type=atom_selection
    .short_caption=Restrained selection
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
    self.ncs_dihedral_proxies = None
    self.name_hash = self.build_name_hash(pdb_hierarchy)
    self.params = params
    self.log = log
    print >> self.log, "Determining NCS matches..."
    pair_hash = {}
    dp_hash = {}
    used_chains = []
    res_match_hash = {}
    i_seq_hash = self.build_i_seq_hash(pdb_hierarchy)
    chain_hash = self.build_chain_hash(pdb_hierarchy)
    name_hash = self.build_name_hash(pdb_hierarchy)
    chain_id_hash = {}
    chains = pdb_hierarchy.models()[0].chains()
    sel_cache = pdb_hierarchy.atom_selection_cache()
    alignments = {}
    for i, chain_i in enumerate(chains):
      if self.get_chain_type(chain_i) == "HETATM":
        continue
      chain_id_hash[chain_i.id] = chain_i
      chain_i_str = "chain '%s'" % chain_i.id
      chain_i_list = [chain_i_str]
      sel_atoms_i = (self.phil_atom_selections_as_i_seqs_multiple(
                   cache=sel_cache,
                   string_list=chain_i_list))
      for chain_j in chains[i+1:]:
        if self.get_chain_type(chain_j) == "HETATM":
          continue
        chain_j_str = "chain '%s'" % chain_j.id
        chain_j_list = [chain_j_str]
        sel_atoms_j = (self.phil_atom_selections_as_i_seqs_multiple(
                         cache=sel_cache,
                         string_list=chain_j_list))
        selections = (sel_atoms_i, sel_atoms_j)
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
             > self.params.similarity ):
          key = (chain_i_str, chain_j_str)
          alignments[key] = residue_match_map
          if used_chains is not None:
            if chain_i.id in used_chains:
              continue
          try:
            pair_hash[chain_i.id].append(chain_j.id)
          except:
            pair_hash[chain_i.id] = []
            pair_hash[chain_i.id].append(chain_j.id)
          used_chains.append(chain_j.id)

    for key in pair_hash.keys():
      ncs_set = []
      ncs_set.append(key)
      for add_chain in pair_hash[key]:
        ncs_set.append(add_chain)
      self.ncs_groups.append(ncs_set)

    for dp in geometry.dihedral_proxies:
      dp_hash[dp.i_seqs] = dp

    super_hash = {}
    res_match_master = {}
    for i, group in enumerate(self.ncs_groups):
      for chain_i in group:
        working_chain = chain_id_hash[chain_i]
        c_atoms = working_chain.atoms()
        chain_i_str = "chain '%s'" % chain_i
        for atom in c_atoms:
          for chain_j in group:
            if chain_i == chain_j:
              continue
            chain_j_str = "chain '%s'" % chain_j
            res_key = self.name_hash[atom.i_seq][4:]
            atom_key = self.name_hash[atom.i_seq][0:4]
            j_match = None
            key = (chain_i_str, chain_j_str)
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

    for dp in geometry.dihedral_proxies:
      temp = dict()
      for i_seq in dp.i_seqs:
        cur_matches = super_hash.get(i_seq)
        if cur_matches is None:
          continue
        for key in cur_matches.keys():
          try:
            temp[key].append(cur_matches[key])
          except:
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
    if self.params.verbose:
      self.show_ncs_summary(log=log)
    print >> self.log, "Initializing torsion NCS restraints..."
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)

  def get_chain_type(self, chain):
    macro_count = 0
    micro_count = 0
    chain_length = chain.residue_groups_size()
    for conformer in chain.conformers():
      for residue in conformer.residues():
        if (common_residue_names_get_class(residue.resname) == \
            'other' or \
            common_residue_names_get_class(residue.resname) == \
            'common_small_molecule' or \
            common_residue_names_get_class(residue.resname) == \
            'common_element'):
          micro_count += 1
        else:
          macro_count += 1
    if micro_count / chain_length >= .50:
      return "HETATM"
    else:
      return "ATOM"

  def selection(self, string, cache):
    return cache.selection(
      string=string)

  def iselection(self, string, cache=None):
    return self.selection(string=string, cache=cache).iselection()

  def phil_atom_selection_multiple(
        self,
        cache,
        string_list,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    result = []
    for string in string_list:
      if (string is None):
        if (allow_none): return None
        raise Sorry('Atom selection cannot be None:\n  =None')
      elif (string is Auto):
        if (allow_auto): return Auto
        raise Sorry('Atom selection cannot be Auto:\n  %s=Auto')
      try:
          result.append(self.selection(string=string, cache=cache).iselection())
      except KeyboardInterrupt: raise
      except Exception, e: # keep e alive to avoid traceback
        fe = format_exception()
        raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
          'reference_group', string, fe))
      if (raise_if_empty_selection and result.count(True) == 0):
        raise Sorry('Empty atom selection:\n  %s=%s' % (
          'reference_group', string))
    return result

  def phil_atom_selections_as_i_seqs_multiple(self,
                                              cache,
                                              string_list):
    result = []
    iselection = self.phil_atom_selection_multiple(
          cache=cache,
          string_list=string_list,
          raise_if_empty_selection=False)
    for i in iselection:
      if (i.size() == 0):
        raise Sorry("No atom selected")
      for atom in i:
        result.append(atom)
    return result

  def is_residue_in_selection(self, i_seqs, selection):
    for i_seq in i_seqs:
      if i_seq not in selection:
        return False
    return True

  def get_i_seqs(self, atoms):
    i_seqs = []
    for atom in atoms:
      i_seqs.append(atom.i_seq)
    return i_seqs

  def extract_sequence_and_sites(self, pdb_hierarchy, selection):
    seq = []
    result = []
    counter = 0
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          if(len(rg.unique_resnames())==1):
            resname = rg.unique_resnames()[0]
            olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
            atoms = rg.atoms()
            i_seqs = self.get_i_seqs(atoms)
            if(olc!="X") and self.is_residue_in_selection(i_seqs, selection):
              seq.append(olc)
              result.append(group_args(i_seq = counter, rg = rg))
              counter += 1
    return "".join(seq), result

  def _alignment(self, pdb_hierarchy,
                    params,
                    selections,
                    log=None):
    if(log is None): log = sys.stdout
    res_match_hash = {}
    model_mseq_res_hash = {}
    model_seq, model_structures = self.extract_sequence_and_sites(
      pdb_hierarchy=pdb_hierarchy,
      selection=selections[0])
    ref_mseq_res_hash = {}
    ref_seq, ref_structures = self.extract_sequence_and_sites(
      pdb_hierarchy = pdb_hierarchy,
      selection=selections[1])
    for struct in model_structures:
      model_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    for struct in ref_structures:
      ref_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    align_obj = mmtbx.alignment.align(
      seq_a                 = model_seq,
      seq_b                 = ref_seq,
      gap_opening_penalty   = params.alignment.gap_opening_penalty,
      gap_extension_penalty = params.alignment.gap_extension_penalty,
      similarity_function   = params.alignment.similarity_matrix,
      style                 = params.alignment.alignment_style)
    alignment = align_obj.extract_alignment()
    matches = alignment.matches()
    exact_match_selections = alignment.exact_match_selections()
    exact_a = tuple(exact_match_selections[0])
    exact_b = tuple(exact_match_selections[1])
    for i, i_seq in enumerate(alignment.i_seqs_a):
      if i_seq != None:
        if alignment.i_seqs_b[i] != None and matches[i] in ['*','|']:
          res_match_hash[model_mseq_res_hash[i_seq]] = \
            ref_mseq_res_hash[alignment.i_seqs_b[i]]
    #print >> log, "  --> aligning model sequence to reference sequence"
    #alignment.pretty_print(block_size  = 50,
    #                       n_block     = 1,
    #                       top_name    = "model",
    #                       bottom_name = "ref",
    #                       out         = log)
    return res_match_hash

  def show_ncs_summary(self, log=None):
    if(log is None): log = sys.stdout
    print >> log, "--------------------------------------------------------"
    print >> log, "Torsion NCS Matching Summary:"
    for dp_set in self.dp_ncs:
      if len(dp_set) < 2:
        continue
      dp_text = "NCS dihedral:\n"
      for dp in dp_set:
        for i_seq in dp.i_seqs:
          dp_text += self.name_hash[i_seq]
        dp_text += "\n"
      print >> log, dp_text
    print >> log, "--------------------------------------------------------"

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
        "** WARNING: No dihedral NCS restraints found!!" + \
        "  Please check parameters. **"
    else:
      print >> log, \
        "Number of dihedral NCS restraints: %d" \
          % len(self.ncs_dihedral_proxies)


  def update_dihedral_ncs_restraints(self,
                                     geometry,
                                     sites_cart,
                                     log=None):
    if log is None:
      log = sys.stdout
    print >> log, "Updating dihedral NCS restraints..."
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
          if self.angle_distance(i, j) <= self.params.cutoff:
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
        target_angle = self.get_angle_average(cluster)
        if self.params.target_damping:
          for c in cluster:
            c_dist = self.angle_distance(c, target_angle)
            if c_dist > self.params.damping_limit:
              d_target = self.get_angle_average([c, target_angle])
              target_angles[c] = d_target
            else:
              target_angles[c] = target_angle
        else:
          for c in cluster:
            target_angles[c] = target_angle
    return target_angles

  def angle_distance(self, angle1, angle2):
    distance = math.fabs(angle1 - angle2)
    if distance > 180.0:
      distance -= 360.0
    return math.fabs(distance)

  def get_angle_average(self, angles):
    n_angles = len(angles)
    sum = 0.0
    a1 = angles[0]
    if a1 > 180.0:
      a1 -= 360.0
    elif a1 < -180.0:
      a1 += 360.0
    sum += a1
    for angle in angles[1:]:
      a2 = angle
      if (a1 - a2) > 180.0:
        a2 += 360.0
      elif (a1 - a2) < -180.0:
        a2 -= 360.0
      sum += a2
    average = sum / n_angles
    return average

  def build_name_hash(self, pdb_hierarchy):
    i_seq_name_hash = dict()
    for atom in pdb_hierarchy.atoms():
      atom_name = atom.pdb_label_columns()[0:4]
      resname = atom.pdb_label_columns()[5:8]
      updated_resname = self.modernize_rna_resname(resname)
      if common_residue_names_get_class(updated_resname) == "common_rna_dna":
        updated_atom = self.modernize_rna_atom_name(atom=atom_name)
      else:
        updated_atom = atom_name
      key = updated_atom+atom.pdb_label_columns()[4:5]+\
            updated_resname+atom.pdb_label_columns()[8:]
      i_seq_name_hash[atom.i_seq]=key
    return i_seq_name_hash

  def build_i_seq_hash(self, pdb_hierarchy):
    name_i_seq_hash = dict()
    for atom in pdb_hierarchy.atoms():
      atom_name = atom.pdb_label_columns()[0:4]
      resname = atom.pdb_label_columns()[5:8]
      updated_resname = self.modernize_rna_resname(resname)
      if common_residue_names_get_class(updated_resname) == "common_rna_dna":
        updated_atom = self.modernize_rna_atom_name(atom=atom_name)
      else:
        updated_atom = atom_name
      key = updated_atom+atom.pdb_label_columns()[4:5]+\
            updated_resname+atom.pdb_label_columns()[8:]
      name_i_seq_hash[key]=atom.i_seq
    return name_i_seq_hash

  def build_chain_hash(self, pdb_hierarchy):
    chain_hash = dict()
    for chain in pdb_hierarchy.chains():
      for atom in chain.atoms():
        chain_hash[atom.i_seq] = chain.id
    return chain_hash

  def modernize_rna_resname(self, resname):
    if common_residue_names_get_class(resname,
         consider_ccp4_mon_lib_rna_dna=True) == "common_rna_dna" or \
       common_residue_names_get_class(resname,
         consider_ccp4_mon_lib_rna_dna=True) == "ccp4_mon_lib_rna_dna":
      tmp_resname = resname.strip()
      if len(tmp_resname) == 1:
        return "  "+tmp_resname
      elif len(tmp_resname) == 2:
        if tmp_resname[0:1].upper() == 'D':
          return " "+tmp_resname.upper()
        elif tmp_resname[1:].upper() == 'D':
          return " D"+tmp_resname[0:1].upper()
        elif tmp_resname[1:].upper() == 'R':
          return "  "+tmp_resname[0:1].upper()
    return resname

  def extract_sequence_from_chain(self, chain):
    seq = []
    res_seq = []
    for rg in chain.residue_groups():
      if(len(rg.unique_resnames())==1):
        resname = rg.unique_resnames()[0]
        olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
        if(olc!="X"):
          seq.append(olc)
          res_seq.append(rg.resid())
    return "".join(seq), res_seq

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
              except:
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
                      remove_clusters_with_all_h=False,
                      log=None)
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
