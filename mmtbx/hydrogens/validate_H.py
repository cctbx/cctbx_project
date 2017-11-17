from __future__ import division
from iotbx.pdb import common_residue_names_get_class
from libtbx import group_args
#from libtbx.utils import Sorry
from mmtbx.rotamer import rotamer_eval

def is_hydrogen(atom):
  answer = False
  if (atom.element.strip().upper() == 'H'):
    answer = True
  return answer

def is_deuterium(atom):
  answer = False
  if (atom.element.strip().upper() == 'D'):
    answer = True
  return answer

class validate_H():
  """ This class is for the validation of H and D atoms, especially for models
  obtained by neutron diffraction."""
  def __init__(self, model):
    self.model = model
    self.pdb_hierarchy = self.model.get_hierarchy()
  # results
    self.overall_counts_hd = None
    self.hd_exchanged_sites = None
    self.hd_sites_analysis = None
    self.renamed = None
    self.missing_HD_atoms = None

  def validate_inputs(self):
    if not self.model.has_hd:
      #raise Sorry("There are no H or D atoms in the model.")
      return 0

  def get_missing_h_in_residue(self, residue, mon_lib_srv):
    missing = []
    ca_xyz, xyz = None, None
    mlq = rotamer_eval.mon_lib_query(residue.resname.strip().upper(), mon_lib_srv)
    if mlq is not None:
      #hd_aliases = mlq.hydrogen_deuterium_aliases()
      atom_name_list = []
      for atom in residue.atoms():
        if atom.name == " CA ":
          ca_xyz = atom.xyz
        if (not atom.element_is_hydrogen()):
          xyz = atom.xyz
        atom_name = atom.name.strip().upper()
        atom_name_list.append(atom_name)
        if is_deuterium(atom):
          atom_name_list.append(atom_name.replace('D','H',1))
        #if atom_name in hd_aliases:
        #  atom_name_list.append(hd_aliases[atom_name])
      if not ca_xyz:
        ca_xyz = xyz
      reference_H = []
      atom_dict = mlq.atom_dict()
      for at in atom_dict:
        reference_H.append(at)
      reference_non_H = []
      for non in mlq.non_hydrogen_atoms():
        reference_non_H.append(non.atom_id.strip().upper())
      reference_list = [x for x in reference_H if x not in reference_non_H]
      alternative_names = [
        ('HA1', 'HA2', 'HA3'),
        ('HB1', 'HB2', 'HB3'),
        ('HG1', 'HG2', 'HG3'),
        ('HD1', 'HD2', 'HD3'),
        ('HE1', 'HE2', 'HE3'),
        ('HG11', 'HG12', 'HG13')
        ]
      for alts in alternative_names:
        if (alts[0] in reference_list and alts[1] in reference_list):
          if (atom_dict[alts[0]].type_energy == 'HCH2' and
              atom_dict[alts[1]].type_energy == 'HCH2'):
            reference_list.append(alts[2])
            reference_list.remove(alts[0])
            if alts[2].replace('H','D',1) in atom_name_list:
              atom_name_list.append(alts[2])
      for atom_name in reference_list:
        if atom_name not in atom_name_list:
          if (atom_name == 'H' and
            ('H1' in atom_name_list and 'H2' in atom_name_list and 'H3' in atom_name_list)):
            continue
          atom_temp = atom_name.replace("*", "'")
          if atom_name.upper() == "O1P":
            atom_temp = "OP1"
          elif atom_name.upper() == "O2P":
            atom_temp = "OP2"
          if atom_temp not in atom_name_list:
            missing.append(atom_name)
    return xyz, missing

  def missing_H(self):
    missing_HD_atoms = []
    from mmtbx.rotamer import rotamer_eval
    #protein = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
    #  "modified_rna_dna", "ccp4_mon_lib_rna_dna"]
    get_class = common_residue_names_get_class
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          missing = []
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            if (get_class(name=residue.resname ) == 'common_water'): continue
            #if (get_class(name=residue.resname)) not in protein: continue
            residue_id = residue.id_str()
            xyz, missing_list = self.get_missing_h_in_residue(
              residue=residue,
              mon_lib_srv=self.model.get_mon_lib_srv())
            if missing_list:
              for atom_name in missing_list:
                if atom_name not in missing:
                  missing.append(atom_name)
          if missing:
            missing_HD_atoms.append((residue_id, missing, xyz))
    self.missing_HD_atoms = missing_HD_atoms

  def get_atom(self, residue_group, name):
    for atom_group in residue_group.atom_groups():
      atom = atom_group.get_atom(name.strip())
      if atom: return atom
    return None

  def find_closest_hd(self, atom, residue_group, neighbors, minimum_distance):
    exchange_atom = None
    hd_without_altloc = []
    if (is_hydrogen(atom)):
      element = 'D'
    else:
      element = 'H'
    for _atom in residue_group.atoms():
      if (_atom.parent().altloc == ''):
        hd_without_altloc.append(_atom.name)
      if _atom.element.strip().upper()!= element: continue
      if neighbors:
        if _atom.i_seq not in neighbors: continue
      current_distance = _atom.distance(atom)
      if current_distance <= minimum_distance:
        minimum_distance = current_distance
        exchange_atom = _atom
    return exchange_atom, hd_without_altloc

  def find_exchanged_pair(self, atom_H, atom_D, residue_group, fsc1):
    minimum_distance = atom_H.distance(atom_D)
    if (atom_H.parent().altloc == '' and atom_D.parent().altloc != ''):
      atom = atom_D
    #elif (atom_D.parent().altloc == '' and atom_H.parent().altloc == ''):
    #  STOP()
    else:
      atom = atom_H
    exchange_atom, hd_without_altloc = self.find_closest_hd(
      atom = atom,
      residue_group = residue_group,
      neighbors = list(fsc1[atom.i_seq]),
      minimum_distance = minimum_distance)
    if exchange_atom is not None:
      if (is_hydrogen(atom)):
        atom_H_new, atom_D_new = atom, exchange_atom
      else:
        atom_H_new, atom_D_new = exchange_atom, atom
      self.rename_hd(
        atom_H = atom_H_new,
        atom_D = atom_D_new,
        hd_without_altloc = hd_without_altloc,
        residue_group = residue_group)
    else:
      return None, None
    return atom_H_new, atom_D_new

  def rename_hd(self, atom_H, atom_D, hd_without_altloc, residue_group):
    atom_D_old_name = atom_D.name
    atom_D_new_name = atom_H.name.replace('H','D',1)
    if atom_D_old_name == atom_D_new_name:
      return
    if atom_D_new_name not in hd_without_altloc:
      atom_with_same_target_name = self.get_atom(
        residue_group = residue_group,
        name = atom_D_new_name)
      if atom_with_same_target_name is not None:
        atom_with_same_target_name.set_name(atom_D_old_name)
        self.renamed.append({
          'atom':atom_with_same_target_name, 'oldname': atom_D_new_name})
      atom_D.set_name(atom_D_new_name)
      self.renamed.append({'atom':atom_D, 'oldname': atom_D_old_name})
    else:
      atom_H_old_name = atom_H.name
      atom_H_new_name = atom_D.name.replace('D','H',1)
      if atom_H_old_name == atom_H_new_name:
        return
      atom_with_same_target_name = self.get_atom(
        residue_group = residue_group,
        name = atom_H_new_name)
      if atom_with_same_target_name is not None:
        atom_with_same_target_name.set_name(atom_H_old_name)
        self.renamed.append({
          'atom':atom_with_same_target_name, 'oldname': atom_H_new_name})
      atom_H.set_name(atom_H_new_name)
      self.renamed.append({'atom':atom_H, 'oldname': atom_H_old_name})

  def get_exchanged_sites_and_curate_swapped(self, pdb_hierarchy):
    self.renamed = []
    geometry_restraints = self.model.get_restraints_manager().geometry
    fsc1 = geometry_restraints.shell_sym_tables[1].full_simple_connectivity()
    protein = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
      "modified_rna_dna", "ccp4_mon_lib_rna_dna"]
    get_class = common_residue_names_get_class
    hd_exchanged_sites = {}
    for residue_group in pdb_hierarchy.residue_groups():
      for atom in residue_group.atoms():
        resname = atom.parent().resname
        if (get_class(name=resname) in protein):
          if (is_hydrogen(atom)):
            atom_H = atom
            name_atom_H = atom_H.name
            name_atom_D = name_atom_H.replace('H','D',1)
            atom_D = self.get_atom(
              residue_group = residue_group,
              name = name_atom_D)
            if atom_D is None: continue
            if atom_H.xyz == atom_D.xyz:
              hd_exchanged_sites[atom_H.i_seq] = [atom_H, atom_D]
            else:
              atom_H_new, atom_D_new = self.find_exchanged_pair(
                atom_H        = atom_H,
                atom_D        = atom_D,
                residue_group = residue_group,
                fsc1          = fsc1)
              if atom_H_new is not None:
                hd_exchanged_sites[atom_H_new.i_seq] = [atom_H_new, atom_D_new]
    self.hd_exchanged_sites = hd_exchanged_sites

  def analyze_hd_sites(self):
    sites_different_xyz = []
    sites_different_b = []
    sites_sum_occ_not_1 = []
    sites_occ_sum_no_scattering = []
    rotatable_hd_selection = self.model.rotatable_hd_selection()
    eps_xyz = 0.001
    eps_b = 0.01
    delta_occ_sum = 0.001
    occ_h_zero_scattering = 0.64
    eps_occ_zero_scatt = 0.05
    max_distance_between_rotatable_H = 1.0
    for iseq in self.hd_exchanged_sites:
      atom_H = self.hd_exchanged_sites[iseq][0]
      atom_D = self.hd_exchanged_sites[iseq][1]
      # H/D at different positions
      delta_xyz = atom_H.distance(atom_D)
      if (delta_xyz >= eps_xyz):
        sites_different_xyz.append(
          (atom_H.id_str(),atom_D.id_str(),delta_xyz,atom_H.xyz,atom_D.xyz))
      # H/D with different B
      delta_b = abs(atom_H.b - atom_D.b)
      if (delta_b >= eps_b):
        delta_b = atom_H.b - atom_D.b
        sites_different_b.append(
          (atom_H.id_str(), atom_D.id_str(), delta_b, atom_H.xyz, atom_D.xyz))
      # H/D with sum of occupancies lt or gt 1
      occupancy_sum = atom_H.occ + atom_D.occ
      if (abs(1-occupancy_sum) >= delta_occ_sum):
        sites_sum_occ_not_1.append(
          (atom_H.id_str(),atom_D.id_str(),occupancy_sum,atom_H.xyz,atom_D.xyz))
      # rotatable H/D with zero scattering sum
      if ((atom_H.i_seq in rotatable_hd_selection) and
          (atom_D.i_seq in rotatable_hd_selection)):
        if (atom_H.distance(atom_D) < max_distance_between_rotatable_H):
          if ((abs(atom_H.occ-occ_h_zero_scattering) <= eps_occ_zero_scatt)
              and (abs(atom_D.occ-(1-occ_h_zero_scattering)) <= eps_occ_zero_scatt)):
            sites_occ_sum_no_scattering.append(
              (atom_H.id_str(),atom_D.id_str(),atom_H.occ,atom_D.occ,atom_H.xyz))

    self.hd_sites_analysis = group_args(
      sites_different_xyz = sites_different_xyz,
      sites_different_b   = sites_different_b,
      sites_sum_occ_not_1 = sites_sum_occ_not_1,
      sites_occ_sum_no_scattering = sites_occ_sum_no_scattering)

  def get_overall_counts(self):
    xray_structure = self.model.get_xray_structure()
    selection=xray_structure.scatterers().extract_scattering_types()!="unknown"
    xrs_sel = xray_structure.select(selection)
    sct = xrs_sel.scatterers().extract_scattering_types()
    count_h = (sct=="H").count(True)
    count_d = (sct=="D").count(True)
    count_water = 0
    resname_classes = self.pdb_hierarchy.overall_counts().resname_classes
    if 'common_water' in resname_classes.keys():
      count_water = resname_classes['common_water']
    return count_h, count_d, count_water

  def count_hd_atoms(self):
    count_h, count_d, n_water = self.get_overall_counts()
    protein = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
      "modified_rna_dna", "ccp4_mon_lib_rna_dna"]
    get_class = common_residue_names_get_class
    count_hd_atoms_protein = 0
    count_h_protein, count_d_protein = 0, 0
    count_h_water, count_d_water = 0, 0
    count_water = 0
    count_water_0h, count_water_1h, count_water_2h = 0, 0, 0
    count_water_altconf = 0
    count_water_no_oxygen = 0
    hd_atoms_with_occ_0 = []
    for residue_group in self.pdb_hierarchy.residue_groups():
      for resname in residue_group.unique_resnames():
        if (get_class(name=resname) == 'common_water'):
          count_water +=1
          count_hd_in_rg, count_o_in_rg = 0, 0
          for atom in residue_group.atoms():
            is_alt_conf = False
            # XXX No break down for water in alt conf for now
            if (atom.parent().altloc != ''):
              count_water_altconf +=1
              is_alt_conf = True
              break
            else:
              if (atom.element_is_hydrogen()):
                count_hd_in_rg += 1
              elif (atom.element.strip().upper() == 'O'):
                count_o_in_rg += 1
          if not is_alt_conf:
            if   count_hd_in_rg == 1 and count_o_in_rg == 1:
              count_water_1h += 1
            elif count_hd_in_rg == 2 and count_o_in_rg == 1:
              count_water_2h += 1
            elif count_hd_in_rg == 0 and count_o_in_rg == 1:
              count_water_0h += 1
            elif count_o_in_rg == 0:
              count_water_no_oxygen += 1
      for atom in residue_group.atoms():
        resname = atom.parent().resname
        if (get_class(name=resname) in protein):
          if (not atom.element_is_hydrogen()):
            continue
          count_hd_atoms_protein += 1
          if (atom.occ == 0):
            hd_atoms_with_occ_0.append((atom.id_str(), atom.xyz))
          if (is_hydrogen(atom)):
            count_h_protein += 1
          elif (is_deuterium(atom)):
            count_d_protein += 1
        elif (get_class(name=resname) == 'common_water'):
          #if (not atom.element_is_hydrogen()): continue
          if (is_hydrogen(atom)):
            count_h_water += 1
          elif (is_deuterium(atom)):
            count_d_water += 1
    assert (count_hd_atoms_protein == count_h_protein + count_d_protein)
    assert (count_water_1h + count_water_2h + count_water_0h + \
      count_water_altconf + count_water_no_oxygen == count_water)
    assert (count_water == n_water)

    self.overall_counts_hd = group_args(
      count_h               = count_h,
      count_d               = count_d,
      count_h_protein       = count_h_protein,
      count_d_protein       = count_d_protein,
      count_h_water         = count_h_water,
      count_d_water         = count_d_water,
      count_water           = count_water,
      count_water_0h        = count_water_0h,
      count_water_1h        = count_water_1h,
      count_water_2h        = count_water_2h,
      count_water_altconf   = count_water_altconf,
      count_water_no_oxygen = count_water_no_oxygen,
      hd_atoms_with_occ_0   = hd_atoms_with_occ_0
      )

  def run(self):
    self.get_exchanged_sites_and_curate_swapped(
      pdb_hierarchy = self.pdb_hierarchy)
    self.count_hd_atoms()
    if self.hd_exchanged_sites:
      self.analyze_hd_sites()
    self.missing_H()

  def get_results(self):
    return group_args(
        overall_counts_hd     = self.overall_counts_hd,
        hd_exchanged_sites    = self.hd_exchanged_sites,
        hd_sites_analysis     = self.hd_sites_analysis,
        pdb_hierarchy_curated = self.pdb_hierarchy,
        renamed               = self.renamed,
        missing_HD_atoms      = self.missing_HD_atoms
        )
