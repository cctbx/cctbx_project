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
from libtbx.str_utils import make_sub_header
import sys, math

TOP_OUT_FLAG = True

torsion_ncs_params = iotbx.phil.parse("""
 sigma = 1.0
   .type = float
 limit = 7.5
   .type = float
 cutoff = 20.0
   .type = float
 slack = 0.0
   .type = float
 hydrogens = False
   .type = bool
 main_chain = True
   .type = bool
 side_chain = True
   .type = bool
 fix_outliers = True
   .type = bool
 verbose = False
   .type = bool
 edits
   .short_caption = Edit reference model restraints
   .style = menu_item parent_submenu:reference_model auto_align noauto
   .expert_level = 2
 {
   include scope \
     mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
 }
 remove
   .short_caption = Remove geometry restraints from reference model
   .expert_level = 2
   .style = menu_item parent_submenu:reference_model auto_align noauto
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
  alignment_style =  local *global
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
  .short_caption=Reference group
  .style = noauto auto_align menu_item parent_submenu:reference_model
 {
  reference=None
    .type=str
    .short_caption=Reference selection
    .style = selection
  selection=None
    .type=str
    .short_caption=Restrained selection
    .style = selection
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
    super_hash = {}
    pair_hash = {}
    dp_hash = {}
    match_master = []
    used_chains = []
    i_seq_hash = self.build_i_seq_hash(pdb_hierarchy)
    chain_hash = self.build_chain_hash(pdb_hierarchy)
    chains = pdb_hierarchy.models()[0].chains()
    for i, chain_i in enumerate(chains):
      seq_i = self.extract_sequence_from_chain(chain=chain_i)
      if len(seq_i) == 0:
        continue
      chain_hash[chain_i.id] = chain_i
      for chain_j in chains[i+1:]:
        seq_j = self.extract_sequence_from_chain(chain=chain_j)
        if len(seq_j) == 0:
          continue
        align_obj = mmtbx.alignment.align(
          seq_a                 = seq_i,
          seq_b                 = seq_j,
          gap_opening_penalty   = params.alignment.gap_opening_penalty,
          gap_extension_penalty = params.alignment.gap_extension_penalty,
          similarity_function   = params.alignment.similarity_matrix,
          style                 = params.alignment.alignment_style)
        if (align_obj.score()/len(seq_i)) >= .80:
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
    for i, group in enumerate(self.ncs_groups):
      match_master.append(dict())
      for chain_id in group:
        working_chain = chain_hash[chain_id]
        c_atoms = working_chain.atoms()
        for atom in c_atoms:
          key = self.name_hash[atom.i_seq][0:8]+self.name_hash[atom.i_seq][-5:]
          try:
            match_master[i][key].append(atom.i_seq)
          except:
            match_master[i][key] = []
            match_master[i][key].append(atom.i_seq)
    for group in match_master:
      for key in group.keys():
        for i in group[key]:
          super_hash[i] = dict()
          for j in group[key]:
            if i == j:
              continue
            else:
              super_hash[i][chain_hash[j]] = j
    for dp in geometry.dihedral_proxies:
      temp = dict()
      for i_seq in dp.i_seqs:
        try:
          cur_matches = super_hash[i_seq]
        except:
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
        try:
          if dp_hash[tuple(temp[key])] is not None:
            dp_match.append(dp_hash[tuple(temp[key])])
            dp_hash[tuple(temp[key])] = None
        except:
          continue
      dp_hash[dp.i_seqs] = None
      if len(dp_match) > 1:
        self.dp_ncs.append(dp_match)
    if self.params.verbose:
      self.show_ncs_summary(log=log)
    print >> self.log, "Initializing torsion NCS restraints..."
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)
    res_match_master = []
    for i, group in enumerate(match_master):
      res_match_master.append(dict())
      for atom_key in group.keys():
        if atom_key[0:4] == " CA ":
          for i_seq in group[atom_key]:
            chain = chain_hash[i_seq]
            try:
              if chain not in res_match_master[i][atom_key[4:]]:
                res_match_master[i][atom_key[4:]].append(chain)
            except:
              res_match_master[i][atom_key[4:]] = []
              res_match_master[i][atom_key[4:]].append(chain)
    self.res_match_master = res_match_master
    #print res_match_master
    #STOP()
    #print >> log, "NCS summary: %d, %d" % (
    #  len(geometry.dihedral_proxies),
    #  len(self.ncs_dihedral_proxies))

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
      #dp_text = "NCS dihedral:\n"
      angles = []
      temp_match = {}
      for dp in dp_set:
        di = cctbx.geometry_restraints.dihedral(sites_cart=sites_cart, proxy=dp)
        angle = di.angle_model
        angles.append(angle)
        temp_match[dp.i_seqs] = angle
        #for i_seq in dp.i_seqs:
        #  dp_text += self.name_hash[i_seq]
        #dp_text += "\n"
      #print >> log, dp_text
      target_angles = self.get_target_angles(angles=angles)
      #print >> log, target_angles
      for dp in dp_set:
        target_angle = target_angles[temp_match[dp.i_seqs]]
        if target_angle is not None:
          dp_add = cctbx.geometry_restraints.dihedral_proxy(
            i_seqs=dp.i_seqs,
            angle_ideal=target_angle,
            weight=1/self.sigma**2,
            limit=self.limit,
            top_out=TOP_OUT_FLAG)
          self.ncs_dihedral_proxies.append(dp_add)

  def update_dihedral_ncs_restraints(self,
                                     geometry,
                                     sites_cart,
                                     log=None):
    if log is None:
      log = sys.stdout
    self.generate_dihedral_ncs_restraints(sites_cart=sites_cart,
                                          log=log)
    self.add_ncs_dihedral_proxies(geometry=geometry)
    print >> log, "Updating dihedral NCS restraints..."

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
    for rg in chain.residue_groups():
      if(len(rg.unique_resnames())==1):
        resname = rg.unique_resnames()[0]
        olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
        if(olc!="X"):
          seq.append(olc)
    return "".join(seq)

  def add_ncs_dihedral_proxies(self, geometry):
    geometry.reference_dihedral_proxies= \
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
    for group in self.res_match_master:
      for key in group.keys():
        for chain in group[key]:
          try:
            res_key = chain+key[-5:-1]+' '+key[0:4]
            if model_hash[res_key] == "OUTLIER":
              rotamer = None
              score = 0.0
              for chain_j in group[key]:
                if chain_j == chain:
                  continue
                j_key = chain_j+key[-5:-1]+' '+key[0:4]
                try:
                  if model_hash[j_key] != "OUTLIER":
                    if rotamer == None:
                      rotamer = j_key
                      score = model_score[j_key]
                    else:
                      if model_score[j_key] > score:
                        rotamer = j_key
                        score = model_score[j_key]
                except:
                  pass
              if rotamer != None:
                fix_list[res_key] = rotamer
          except:
            pass

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
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            try:
              if key in fix_list.keys():
                axis_and_atoms_to_rotate=fit_rotamers.axes_and_atoms_aa_specific(
                      residue=atom_group,
                      mon_lib_srv=mon_lib_srv,
                      remove_clusters_with_all_h=False,
                      log=None)
                m_chis = model_chis[key]
                r_chis = model_chis[fix_list[key]]
                assert len(m_chis) == len(r_chis)
                assert len(m_chis) == len(axis_and_atoms_to_rotate)
                counter = 0
                residue_iselection = atom_group.atoms().extract_i_seq()
                sites_cart_residue = xray_structure.sites_cart().select(residue_iselection)
                for aa in axis_and_atoms_to_rotate:
                  axis = aa[0]
                  atoms = aa[1]
                  atom_group.atoms().set_xyz(new_xyz=sites_cart_residue)
                  new_xyz = flex.vec3_double()
                  angle_deg = r_chis[counter] - m_chis[counter]
                  if angle_deg < 0:
                    angle_deg += 360.0
                  for atom in atoms:
                    new_xyz = fit_rotamers.rotate_point_around_axis(
                                axis_point_1=sites_cart_residue[axis[0]],
                                axis_point_2=sites_cart_residue[axis[1]],
                                point=sites_cart_residue[atom],
                                angle_deg=angle_deg)
                    sites_cart_residue[atom] = new_xyz
                  sites_cart_start = sites_cart_start.set_selected(
                        residue_iselection, sites_cart_residue)
                  counter += 1
                xray_structure.set_sites_cart(sites_cart_start)
                print >> log, "Fixed %s rotamer" % key
            except:
              pass
