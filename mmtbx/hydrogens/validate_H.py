from __future__ import absolute_import, division, print_function
from iotbx.pdb import common_residue_names_get_class
from libtbx import group_args
from libtbx.str_utils import make_sub_header
from mmtbx.validation import restraints
from mmtbx.validation.molprobity import mp_geo
#from libtbx.utils import Sorry
from mmtbx.rotamer import rotamer_eval
from six.moves import zip

protein = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
  "modified_rna_dna", "ccp4_mon_lib_rna_dna"]

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

def get_atom_info_if_hd(atoms_info):
  # returns atom_info for H atom, always last H atom in list of atoms_info items
  atom_info_hd = None
  for atom_info in atoms_info:
    element = atom_info.element.strip().upper()
    if element == 'H' or element == 'D':
      atom_info_hd = atom_info
  return atom_info_hd

class validate_H(object):
  """ This class is for the validation of H and D atoms, especially for models
  obtained by neutron diffraction."""
  def __init__(self, model, use_neutron_distances):
  # input
    self.model = model
    self.use_neutron_distances = use_neutron_distances
  # derived
    self.pdb_hierarchy = self.model.get_hierarchy()
  # results
    self.overall_counts_hd = None
    self.hd_exchanged_sites = None
    self.hd_sites_analysis = None
    self.renamed = None
    self.missing_HD_atoms = None
    self.outliers_bonds = None
    self.outliers_angles = None
    self.bond_results = None
    self.curated_hierarchy = None

  def validate_inputs(self):
    if not self.model.has_hd():
      #raise Sorry("There are no H or D atoms in the model.")
      return 0
    # ensure that grm exists
    self.model.get_restraints_manager()

  def get_missing_h_in_residue(self, residue, mon_lib_srv):
    missing = []
    ca_xyz, xyz, xyzh = None, None, None
    mlq = rotamer_eval.mon_lib_query(residue.resname.strip().upper(), mon_lib_srv)
    if mlq is not None:
      #hd_aliases = mlq.hydrogen_deuterium_aliases()
      atom_name_list = []
      for atom in residue.atoms():
        atom_name = atom.name.strip().upper()
        if atom_name == "CA":
          ca_xyz = atom.xyz
        if (not atom.element_is_hydrogen()):
          xyz = atom.xyz
        else:
          atom_name_list.append(atom_name)
          xyzh = atom.xyz
        if is_deuterium(atom):
          atom_name_list.append(atom_name.replace('D','H',1))
        #if atom_name in hd_aliases:
        #  atom_name_list.append(hd_aliases[atom_name])
      if not ca_xyz:
        ca_xyz = xyz
      if not ca_xyz:
        ca_xyz = xyzh
      # Step 1: Get list of expected H and non-H atoms
      reference = []
      atom_dict = mlq.atom_dict()
      for at in atom_dict:
        reference.append(at)
      # Step 2: Get list of expected non-H atoms
      reference_non_H = []
      for non in mlq.non_hydrogen_atoms():
        reference_non_H.append(non.atom_id.strip().upper())
      # Step 3: Get list of expected H atoms only
      reference_hydrogens = [x for x in reference if x not in reference_non_H]
      # Step 4: There can be naming differences: HB1+HB2 or HB2+HB3
      alternative_names = [
        ('HA1', 'HA2', 'HA3'),
        ('HB1', 'HB2', 'HB3'),
        ('HG1', 'HG2', 'HG3'),
        ('HD1', 'HD2', 'HD3'),
        ('HE1', 'HE2', 'HE3'),
        ('HG11', 'HG12', 'HG13')
        ]
      for alts in alternative_names:
        if (alts[0] in reference_hydrogens and alts[1] in reference_hydrogens):
          if (atom_dict[alts[0]].type_energy == 'HCH2' and
              atom_dict[alts[1]].type_energy == 'HCH2'):
            reference_hydrogens.append(alts[2])
            reference_hydrogens.remove(alts[0])
#            if alts[2].replace('H','D',1) in atom_name_list:
#              atom_name_list.append(alts[2])
      for atom_name in reference_hydrogens:
        if atom_name not in atom_name_list:
          if (atom_name == 'H' and
            ('H1' in atom_name_list and 'H2' in atom_name_list and 'H3' in atom_name_list)):
            continue
          if (atom_name.endswith('*')):
            atom_temp = atom_name.replace("*", "'")
          elif (atom_name.endswith('*1')):
            atom_temp = atom_name.replace("*1", "'")
          elif (atom_name.endswith('*2')):
            atom_temp = atom_name.replace("*2", "''")
          else:
            atom_temp = atom_name

#          if atom_name.upper() == "O1P":
#            atom_temp = "OP1"
#          elif atom_name.upper() == "O2P":
#            atom_temp = "OP2"
          if atom_temp not in atom_name_list:
            if atom_temp.endswith("'"):
              missing.append(atom_temp)
            else:
              missing.append(atom_name)
    return ca_xyz, missing

  def missing_hydrogens(self):
    missing_HD_atoms = []
    #from mmtbx.rotamer import rotamer_eval
    get_class = common_residue_names_get_class
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          id_strings, missing, conformers, xyzs = [], [], [], []
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            if (get_class(name=residue.resname ) == 'common_water'): continue
            xyz, missing_list = self.get_missing_h_in_residue(
              residue     = residue,
              mon_lib_srv = self.model.get_mon_lib_srv())
            if missing_list:
              missing.append(missing_list)
              conformers.append(conformer.altloc)
              xyzs.append(xyz)
              id_strings.append(residue.id_str())
          if missing:
            # if all conformers lack the same H atoms, only add once
            if len([list(tupl) for tupl in {tuple(item) for item in missing }]) == 1:
                missing_HD_atoms.append((id_strings[0],
                                         missing[0],
                                         ", ".join(conformers),
                                         xyzs[0]))
            # otherwise, add for each conformer
            else:
              for missing_list, conformer, id_str, xyz in zip(
                                        missing, conformers, id_strings, xyzs):
                missing_HD_atoms.append((id_str,
                                         missing_list,
                                         conformer,
                                         xyz))
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
        self.renamed.append(
                            (atom_with_same_target_name.id_str(), # id_str
                             atom_with_same_target_name.name,     # new name
                             atom_D_new_name,                     # old name
                             atom_with_same_target_name.xyz))     # xyz
      atom_D.set_name(atom_D_new_name)
      self.renamed.append(
                          (atom_D.id_str(),
                           atom_D.name,
                           atom_D_old_name,
                           atom_D.xyz))
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
        self.renamed.append(
                            (atom_with_same_target_name.id_str(),
                             atom_with_same_target_name.name,
                             atom_H_new_name,
                             atom_with_same_target_name.xyz))
      atom_H.set_name(atom_H_new_name)
      self.renamed.append(
                          (atom_H.id_str(),
                           atom_H.name,
                           atom_H_old_name,
                           atom_H.xyz))

  def get_exchanged_sites_and_curate_swapped(self, pdb_hierarchy):
    self.renamed = []
    geometry_restraints = self.model.get_restraints_manager().geometry
    fsc1 = geometry_restraints.shell_sym_tables[1].full_simple_connectivity()
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
    if self.renamed is not None:
      self.curated_hierarchy = pdb_hierarchy

  def analyze_hd_sites(self):
    sites_different_xyz = []
    sites_different_b = []
    sites_sum_occ_not_1 = []
    sites_occ_sum_no_scattering = []
    rotatable_hd_selection = self.model.rotatable_hd_selection()
    eps_xyz = 0.001
    eps_b = 0.01
    delta_occ_sum = 0.001
    occ_h_zero_scattering = 0.64 # value for which sum occ H and D is zero
    eps_occ_zero_scatt = 0.05
    # For rotatable H, H and D may be at different positions
    # However, when they are close to each other, cancellation may occur
    # Introduce max distance, corresponds to approx. 45 deg between O-D and O-H
    max_distance_between_rotatable_H = 0.8
    for iseq in self.hd_exchanged_sites:
      atom_H = self.hd_exchanged_sites[iseq][0]
      atom_D = self.hd_exchanged_sites[iseq][1]
      # H/D at different positions
      delta_xyz = atom_H.distance(atom_D)
      if (delta_xyz >= eps_xyz):
        sites_different_xyz.append(
                                  (atom_H.id_str(),
                                   atom_D.id_str(),
                                   delta_xyz,
                                   atom_H.xyz,
                                   atom_D.xyz))
      # H/D with different B
      delta_b = abs(atom_H.b - atom_D.b)
      if (delta_b >= eps_b):
        delta_b = atom_H.b - atom_D.b
        sites_different_b.append(
                                (atom_H.id_str(),
                                 atom_D.id_str(),
                                 delta_b,
                                 atom_H.xyz,
                                 atom_D.xyz))
      # H/D with sum of occupancies lt or gt 1
      occupancy_sum = atom_H.occ + atom_D.occ
      if (abs(1-occupancy_sum) >= delta_occ_sum):
        sites_sum_occ_not_1.append(
                                  (atom_H.id_str(),
                                   atom_D.id_str(),
                                   occupancy_sum,
                                   atom_H.xyz,
                                   atom_D.xyz))
      # rotatable H/D with zero scattering sum, if closer than cut off apart
      if ((atom_H.i_seq in rotatable_hd_selection) and
          (atom_D.i_seq in rotatable_hd_selection)):
        if (atom_H.distance(atom_D) < max_distance_between_rotatable_H):
          if ((abs(atom_H.occ-occ_h_zero_scattering) <= eps_occ_zero_scatt)
              and
              (abs(atom_D.occ-(1-occ_h_zero_scattering))<= eps_occ_zero_scatt)):
            sites_occ_sum_no_scattering.append(
                                              (atom_H.id_str(),
                                               atom_D.id_str(),
                                               atom_H.occ,
                                               atom_D.occ,
                                               atom_H.xyz,
                                               atom_D.xyz))

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
    get_class = common_residue_names_get_class
    count_hd_atoms_protein = 0
    count_h_protein, count_d_protein = 0, 0
    count_h_water, count_d_water = 0, 0
    count_water = 0
    count_water_0h, count_water_1h, count_water_2h = 0, 0, 0
    count_water_more_h = 0
    count_water_altconf = 0
    count_water_no_oxygen = 0
    hd_atoms_with_occ_0 = []
    single_hd_atoms_occ_lt_1 = []
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
            elif count_hd_in_rg > 2:
              count_water_more_h += 1
      for atom in residue_group.atoms():
        resname = atom.parent().resname
        if (get_class(name=resname) in protein):
          if (not atom.element_is_hydrogen()):
            continue
          count_hd_atoms_protein += 1
          if (atom.occ == 0):
            hd_atoms_with_occ_0.append((atom.id_str(), atom.xyz))
          if (atom.occ <1 and atom.occ > 0 and atom.parent().altloc == ''):
            single_hd_atoms_occ_lt_1.append(
                                            (atom.id_str(),
                                             atom.occ,
                                             atom.xyz))
          if (is_hydrogen(atom)):    count_h_protein += 1
          elif (is_deuterium(atom)): count_d_protein += 1
        elif (get_class(name=resname) == 'common_water'):
          if (is_hydrogen(atom)):    count_h_water += 1
          elif (is_deuterium(atom)): count_d_water += 1
    assert (count_hd_atoms_protein == count_h_protein + count_d_protein)
    assert (count_water_1h + count_water_2h + count_water_0h + \
      count_water_altconf + count_water_no_oxygen + count_water_more_h == count_water)
    assert (count_water == n_water)

    count_h_other = count_h - count_h_protein - count_h_water
    count_d_other = count_d - count_d_protein - count_d_water

    self.overall_counts_hd = group_args(
      count_h                  = count_h,
      count_d                  = count_d,
      count_h_protein          = count_h_protein,
      count_d_protein          = count_d_protein,
      count_h_water            = count_h_water,
      count_d_water            = count_d_water,
      count_h_other            = count_h_other,
      count_d_other            = count_d_other,
      count_water              = count_water,
      count_water_0h           = count_water_0h,
      count_water_1h           = count_water_1h,
      count_water_2h           = count_water_2h,
      count_water_altconf      = count_water_altconf,
      count_water_no_oxygen    = count_water_no_oxygen,
      hd_atoms_with_occ_0      = hd_atoms_with_occ_0,
      single_hd_atoms_occ_lt_1 = single_hd_atoms_occ_lt_1
      )

  def get_hd_state(self):
    hd_state = None
    sel_h = self.pdb_hierarchy.atom_selection_cache().selection(
      string = 'element H')
    sel_d = self.pdb_hierarchy.atom_selection_cache().selection(
      string = 'element D')
    if sel_h.count(True) != 0 and sel_d.count(True) == 0:
      hd_state = 'all_h'
    elif sel_h.count(True) == 0 and sel_d.count(True) != 0:
      hd_state = 'all_d'
    else:
      hd_state = 'h_and_d'
    return hd_state

  def bond_angle_outliers(self):
    get_class = common_residue_names_get_class
    rc = restraints.combined(
           pdb_hierarchy  = self.pdb_hierarchy,
           xray_structure = self.model.get_xray_structure(),
           geometry_restraints_manager = self.model.get_restraints_manager().\
                                         geometry,
           ignore_hd      = False, # important
           outliers_only  = False,
           use_segids_in_place_of_chainids = False)

    bond_mean_delta, n_bonds, bond_mean = 0, 0, 0

    # bond outliers involving hydrogens
    outliers_bonds = []
    for result in rc.bonds.results:
      atom_info_hd = get_atom_info_if_hd(atoms_info = result.atoms_info)
      # Consider only H/D atoms
      if atom_info_hd is not None:
        # Calculate mean bond length and delta for non-water
        # --> used to get rough idea if H are at X-ray or neutron bond lengths.
        if (get_class(name=atom_info_hd.resname) != 'common_water'):
          bond_mean_delta = bond_mean_delta + result.delta
          bond_mean = bond_mean + result.model
          n_bonds += 1
        if result.is_outlier():
          atoms_str = mp_geo.get_atoms_str(atoms_info=result.atoms_info)
          outliers_bonds.append(
                            (atom_info_hd.id_str(),
                             atoms_str,
                             result.model,
                             result.delta,
                             result.target,
                             atom_info_hd.xyz) )
    self.outliers_bonds = outliers_bonds

    if n_bonds:
      bond_mean_delta = bond_mean_delta/n_bonds
      bond_mean = bond_mean/n_bonds

    xray_distances_used = False
    # value 0.08 was obtained by checking all 123 neutron models deposited
    # until Sep 2017 and by analysing delta
    if (bond_mean_delta >= 0.08 and self.use_neutron_distances):
      xray_distances_used = True

    self.bond_results = group_args(
      bond_mean_delta = bond_mean_delta,
      bond_mean = bond_mean,
      xray_distances_used = xray_distances_used
      )

    # angle outliers involving hydrogens
    outliers_angles = []
    for result in rc.angles.results:
      atom_info_hd = get_atom_info_if_hd(atoms_info = result.atoms_info)
      # Consider only H/D atoms
      if atom_info_hd is not None:
        if result.is_outlier():
          atoms_str = mp_geo.get_atoms_str(atoms_info=result.atoms_info)
          outliers_angles.append(
                            [atom_info_hd.id_str(),
                             atoms_str,
                             result.model,
                             result.delta,
                             result.target,
                             atom_info_hd.xyz] )
    self.outliers_angles = outliers_angles

  def run(self):
    # Get overall counts of H and D
    self.count_hd_atoms()
    # Find bond and angle outliers
    self.bond_angle_outliers()
    # Find missing H atoms
    self.missing_hydrogens()
    # if H and D are both present, analyse and curate potential H/D states
    if self.get_hd_state() == 'h_and_d':
      self.get_exchanged_sites_and_curate_swapped(
        pdb_hierarchy = self.pdb_hierarchy.deep_copy())
    # If H/D sites are present, analyze for mismatches
    if self.hd_exchanged_sites:
      self.analyze_hd_sites()

  def get_results(self):
    count_exchanged_sites = None
    if self.hd_exchanged_sites:
      count_exchanged_sites = len(self.hd_exchanged_sites)
    return group_args(
        overall_counts_hd     = self.overall_counts_hd,
        count_exchanged_sites = count_exchanged_sites,
        hd_sites_analysis     = self.hd_sites_analysis,
        renamed               = self.renamed,
        missing_HD_atoms      = self.missing_HD_atoms,
        outliers_bonds        = self.outliers_bonds,
        outliers_angles       = self.outliers_angles,
        bond_results          = self.bond_results
        )

  def get_curated_hierarchy(self):
    return self.curated_hierarchy

class validate_H_results(object):
  '''
  Controller class for displaying results from validate_H
  '''
  def __init__(self, results, log=None):
    self.results = results
    self.log = log

  def formatted_print(self, prefix, values, log):
    maxlen = max([ len(value[0]) for value in values ])
    for value in values:
      base_format = '%s%%-%ds: %%i' % (prefix, maxlen)
      print(base_format % value, file=log)

  def print_overall_results(self, overall_counts_hd, prefix='', log=None):
    if (log is None):
      log = self.log

    oc = overall_counts_hd

    make_sub_header('H/D atoms in the input model', out=log)
    self.hd_overall_values = [
      ('Total number of hydrogen atoms' , oc.count_h),
      ('Total number of deuterium atoms' , oc.count_d),
      ('Number of H atoms (protein)' , oc.count_h_protein),
      ('Number of D atoms (protein)' , oc.count_d_protein),
      ('Number of H atoms (water)' , oc.count_h_water),
      ('Number of D atoms (water)' , oc.count_d_water),
      ('Number of H atoms (other)' , oc.count_h_other),
      ('Number of D atoms (other)' , oc.count_d_other),
    ]
    self.formatted_print(prefix, self.hd_overall_values, log)

    make_sub_header('Water molecules', out=log)
    self.hd_water_values = [
      ('Number of water', oc.count_water),
      ('Number of water with 0 H (or D)', oc.count_water_0h),
      ('Number of water with 1 H (or D)', oc.count_water_1h),
      ('Number of water with 2 H (or D)', oc.count_water_2h),
      ('Number of water in alternative conformation', oc.count_water_altconf),
      ('Number of water without oxygen atom', oc.count_water_no_oxygen)
    ]
    self.formatted_print(prefix, self.hd_water_values, log)

  def print_renamed(self, renamed, prefix='', log=None):
    if (log is None):
      log = self.log

    make_sub_header('The following atoms were renamed:', out=log)
    for entry in renamed:
      id_str = entry[0]
      oldname = entry[2]
      newname = entry[1]
      print('%s%s atom %s --> %s' % (prefix, id_str, oldname, newname),
            file=log)

  def export_renamed_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.renamed):
      for entry in self.results.renamed:
        table.append([entry[0], entry[2], entry[1], None, entry[-1]])
    return table

  def print_atoms_occ_lt_1(self, hd_atoms_with_occ_0, single_hd_atoms_occ_lt_1,
                           prefix='', log=None):
    if (log is None):
      log = self.log

    if hd_atoms_with_occ_0:
      make_sub_header('H (or D) atoms with zero occupancy', out=log)
      for item in hd_atoms_with_occ_0:
        print('%s%s' % (prefix, item[0]), file=log)
    if single_hd_atoms_occ_lt_1:
      make_sub_header('H (or D) atoms with occupancy < 1', out=log)
      for item in single_hd_atoms_occ_lt_1:
        print('%s%s with occupancy %s' % (prefix, item[0], item[1]),
              file=log)

  def export_occupancies_0_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.overall_counts_hd.hd_atoms_with_occ_0):
      for item in self.results.overall_counts_hd.hd_atoms_with_occ_0:
        table.append([item[0].split('"')[1], None, item[-1]])
    return table

  def export_occupancies_lt_1_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.overall_counts_hd.single_hd_atoms_occ_lt_1):
      for item in self.results.overall_counts_hd.single_hd_atoms_occ_lt_1:
        table.append([item[0].split('"')[1], item[1], None, item[-1]])
    return table

  def print_results_hd_sites(
      self, count_exchanged_sites, hd_sites_analysis, overall_counts_hd,
      prefix='', log=None):
    if (log is None):
      log = self.log

    sites_different_xyz = hd_sites_analysis.sites_different_xyz
    sites_different_b   = hd_sites_analysis.sites_different_b
    sites_sum_occ_not_1 = hd_sites_analysis.sites_sum_occ_not_1
    sites_occ_sum_no_scattering = hd_sites_analysis.sites_occ_sum_no_scattering

    make_sub_header('H/D EXCHANGED SITES', out=log)
    self.hd_exchange_values = [
      ('Number of H/D exchanged sites', count_exchanged_sites),
      ('Number of atoms modelled only as H',
       overall_counts_hd.count_h_protein - count_exchanged_sites),
      ('Number of atoms modelled only as D',
       overall_counts_hd.count_d_protein - count_exchanged_sites)
    ]
    self.formatted_print(prefix, self.hd_exchange_values, log)

    if sites_different_xyz:
      print('\n%sH/D pairs not at identical positions:' % prefix, file=log)
      for item in sites_different_xyz:
        print('%s  %s and  %s at distance %.3f' % \
          (prefix, item[0][5:-1], item[1][5:-1], item[2]), file=log)

    if sites_different_b:
      print('\n%sH/D pairs without identical ADPs:' % prefix, file=log)
      for item in sites_different_b:
        print('%s  %s and %s ' % (prefix, item[0][5:-1], item[1][5:-1]),
              file=log)

    if sites_sum_occ_not_1:
      print('\n%sH/D pairs with occupancy sum != 1:' % prefix, file=log)
      for item in sites_sum_occ_not_1:
        print('%s  %s  and %s with occupancy sum %s' %
              (prefix, item[0][5:-1], item[1][5:-1], item[2]), file=log)

    if sites_occ_sum_no_scattering:
      print('\n%sRotatable H/D pairs with zero scattering occupancy sum:' %
            prefix, file=log)
      for item in sites_occ_sum_no_scattering:
        print('%s  %s with occ %s and  %s with occ %s' %
              (prefix, item[0][5:-1], item[2], item[1][5:-1], item[3]),
              file=log)

  def export_sites_different_xyz_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    result = None
    if (self.results.hd_sites_analysis):
      result = self.results.hd_sites_analysis.sites_different_xyz
    if (result):
      for item in result:
        table.append([item[0][5:-1], item[1][5:-1], item[2], None, item[-1]])
    return table

  def export_sites_different_b_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    result = None
    if (self.results.hd_sites_analysis):
      result = self.results.hd_sites_analysis.sites_different_b
    if (result):
      for item in result:
        table.append([item[0][5:-1], item[1][5:-1], None, item[-1]])
    return table

  def export_sites_sum_occ_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    result = None
    if (self.results.hd_sites_analysis):
      result = self.results.hd_sites_analysis.sites_sum_occ_not_1
    if (result):
      for item in result:
        table.append([item[0][5:-1], item[1][5:-1], item[2], None, item[-1]])
    return table

  def export_sites_occ_scattering_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    result = None
    if (self.results.hd_sites_analysis):
      result = self.results.hd_sites_analysis.sites_occ_sum_no_scattering
    if (result):
      for item in result:
        table.append([item[0][5:-1], item[2], item[1][5:-1], item[3],
                      None, item[-1]])
    return table

  def print_missing_HD_atoms(self, missing_HD_atoms, prefix, log=None):
    if (log is None):
      log = self.log

    make_sub_header('MISSING H or D atoms', out=log)
    for item in missing_HD_atoms:
      print('%s%s conformer %s : %s ' % (prefix, item[0][8:-1], item[2], ", ".join(item[1])),
            file=log)

  def export_missing_HD_atoms_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.missing_HD_atoms):
      for item in self.results.missing_HD_atoms:
        table.append([item[0][8:-1], item[2], ', '.join(item[1]),
                      None, item[-1]])
    return table

  def print_outliers_bonds_angles(self, outliers_bonds, outliers_angles,
                                  prefix='', log=None):
    if (log is None):
      log = self.log

    if outliers_bonds:
      make_sub_header('Bond outliers', out=log)
      for item in outliers_bonds:
        print('%s%s, Bond %s, observed: %.3f, delta from target: %.3f' % \
          (prefix, item[0], item[1], item[2], item[3]), file=log)
    if outliers_angles:
      make_sub_header('Angle outliers', out=log)
      for item in outliers_angles:
        print('%s%s, Angle %s, observed: %.3f, delta from target: %.3f' % \
          (prefix, item[0], item[1], item[2], item[3]), file=self.log)

  def export_outliers_bonds_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.outliers_bonds):
      for item in self.results.outliers_bonds:
        table.append([item[0], item[1], item[2], item[3], None, item[-1]])
    return table

  def export_outliers_angles_for_wxGUI(self):
    # last element should be xyz for residue/atom
    table = list()
    if (self.results.outliers_angles):
      for item in self.results.outliers_angles:
        table.append([item[0], item[1], item[2], item[3], None, item[-1]])
    return table

  def print_xray_distance_warning(self):
    print('*'*79, file=self.log)
    print(
      'WARNING: Model has a majority of X-H bonds with X-ray bond lengths.\n \
      Input was to use neutron distances. Please check your model carefully.',
      file=self.log)
    print('*'*79, file=self.log)

  def print_results(self, results=None, prefix= '', log=None):
    if (results is None):
      results = self.results
    if (log is not None):
      self.log = log
    assert (results is not None)
    assert (self.log is not None)

    overall_counts_hd  = results.overall_counts_hd
    count_exchanged_sites = results.count_exchanged_sites
    renamed            = results.renamed
    hd_sites_analysis  = results.hd_sites_analysis
    missing_HD_atoms   = results.missing_HD_atoms
    hd_atoms_with_occ_0 = overall_counts_hd.hd_atoms_with_occ_0
    single_hd_atoms_occ_lt_1 = overall_counts_hd.single_hd_atoms_occ_lt_1
    outliers_bonds     = results.outliers_bonds
    outliers_angles    = results.outliers_angles
    bond_results       = results.bond_results

    if overall_counts_hd:
      self.print_overall_results(overall_counts_hd, prefix=prefix, log=log)
    if renamed:
      self.print_renamed(renamed, prefix=prefix, log=log)
    if hd_atoms_with_occ_0 or single_hd_atoms_occ_lt_1:
      self.print_atoms_occ_lt_1(hd_atoms_with_occ_0, single_hd_atoms_occ_lt_1,
                                prefix=prefix, log=log)
    if count_exchanged_sites is not None:
      self.print_results_hd_sites(
        count_exchanged_sites, hd_sites_analysis, overall_counts_hd,
        prefix=prefix, log=log)
    if missing_HD_atoms:
      self.print_missing_HD_atoms(missing_HD_atoms, prefix=prefix, log=log)
    if outliers_bonds or outliers_angles:
      self.print_outliers_bonds_angles(outliers_bonds, outliers_angles,
                                       prefix=prefix, log=log)
    if bond_results.xray_distances_used:
      self.print_xray_distance_warning()
