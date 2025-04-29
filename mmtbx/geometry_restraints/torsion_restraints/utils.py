from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, null_out
from iotbx.pdb import common_residue_names_get_class
from iotbx.pdb import input
from cctbx.array_family import flex
import math
import sys
from six.moves import zip

def process_reference_files(
      reference_file_list,
      log=None):
  if log is None:
    log = sys.stdout
  reference_hierarchy_list = []
  for file in reference_file_list:
    pdb_io = input(file)
    cur_hierarchy = pdb_io.construct_hierarchy()
    cur_hierarchy.reset_i_seq_if_necessary()
    ter_indices = pdb_io.ter_indices()
    if ter_indices is not None:
      check_for_internal_chain_ter_records(
        pdb_hierarchy=cur_hierarchy,
        ter_indices=ter_indices)
    reference_hierarchy_list.append(cur_hierarchy)
  return reference_hierarchy_list

def get_reference_dihedral_proxies(
      reference_hierarchy_list,
      reference_file_list,
      mon_lib_srv=None,
      ener_lib=None,
      crystal_symmetry=None,
      restraint_objects=None,
      monomer_parameters=None,
      log=None):
  from mmtbx.monomer_library import server
  if log is None:
    log = sys.stdout
  if mon_lib_srv is None:
    mon_lib_srv = server.server()
  if ener_lib is None:
    ener_lib = server.ener_lib()
  reference_dihedral_proxies = {}
  for file_name, pdb_hierarchy in zip(reference_file_list,
                                      reference_hierarchy_list):
    dihedral_proxies = get_complete_dihedral_proxies(
                         pdb_hierarchy=pdb_hierarchy,
                         mon_lib_srv=mon_lib_srv,
                         ener_lib=ener_lib,
                         crystal_symmetry=crystal_symmetry,
                         restraint_objects=restraint_objects,
                         monomer_parameters=monomer_parameters,
                         log=log)
    reference_dihedral_proxies[file_name]=dihedral_proxies
  return reference_dihedral_proxies

def get_complete_dihedral_proxies_2(
    model,
    log=None):
  from six.moves import cStringIO as StringIO
  from cctbx.geometry_restraints import dihedral_proxy_registry
  from mmtbx.conformation_dependent_library import generate_protein_threes
  if log is None:
    log = StringIO
  dihedral_registry = dihedral_proxy_registry(
      strict_conflict_handling=True)
  dihedral_registry.initialize_table()
  grm = model.get_restraints_manager().geometry
  dihedral_proxies = grm.get_dihedral_proxies().deep_copy()
  for p in dihedral_proxies:
    dihedral_registry.add_if_not_duplicated(p)
  for three in generate_protein_threes(
      hierarchy=model.get_hierarchy(),
      geometry=None):
    proxies = three.get_dummy_dihedral_proxies(
        only_psi_phi_pairs=False)
    for p in proxies:
      dihedral_registry.add_if_not_duplicated(p)
  return dihedral_registry.proxies

def get_complete_dihedral_proxies(
      pdb_hierarchy,
      mon_lib_srv=None,
      ener_lib=None,
      crystal_symmetry=None,
      restraint_objects=None,
      monomer_parameters=None,
      log=None):
  #
  # This function is called only for reference files, that were not processed
  # yet. For the main file only get_dihedrals_and_phi_psi below is called.
  # Still used for reference model torsion restraints
  #
  import mmtbx.model
  from mmtbx.monomer_library import server
  assert pdb_hierarchy is not None
  if log is None:
    log = sys.stdout
  if mon_lib_srv is None:
    mon_lib_srv = server.server()
  if ener_lib is None:
    ener_lib = server.ener_lib()
  if pdb_hierarchy is not None:
    if pdb_hierarchy.fits_in_pdb_format():
      raw_records = pdb_hierarchy.as_pdb_string()
    else:
      raw_records = pdb_hierarchy.as_mmcif_string()
  if raw_records is not None:
    if (isinstance(raw_records, str)):
      raw_records = flex.split_lines(raw_records)
  work_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  work_params.pdb_interpretation.c_beta_restraints=False
  work_params.pdb_interpretation.automatic_linking.link_none=True
  work_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
  pdb_inp = input(lines=raw_records, source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      restraint_objects=restraint_objects,
      monomer_parameters=monomer_parameters,
      stop_for_unknowns=False,
      log=null_out())
  model.process(pdb_interpretation_params=work_params,
    make_restraints=True)
  return get_dihedrals_and_phi_psi(model)

def get_dihedrals_and_phi_psi(model):
  from cctbx.geometry_restraints import dihedral_proxy_registry
  dihedral_registry = dihedral_proxy_registry(
      strict_conflict_handling=True)
  dihedral_registry.initialize_table()
  from mmtbx.conformation_dependent_library import generate_protein_threes
  grm = None
  try:
    grm = model._processed_pdb_file.geometry_restraints_manager()
  except AttributeError:
    grm = model.get_restraints_manager().geometry
  assert grm is not None
  dihedral_proxies = grm.get_dihedral_proxies().deep_copy()
  for p in dihedral_proxies:
    dihedral_registry.add_if_not_duplicated(p)
  for three in generate_protein_threes(
      hierarchy=model.get_hierarchy(),
      geometry=None):
    proxies = three.get_dummy_dihedral_proxies(
        only_psi_phi_pairs=False)
    for p in proxies:
      dihedral_registry.add_if_not_duplicated(p)
  return dihedral_registry.proxies

def modernize_rna_resname(resname):
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
    elif tmp_resname in ["ADE", "CYT", "GUA", "URI"]:
      return "  "+tmp_resname[0:1].upper()
  return resname

def modernize_rna_atom_name(atom):
   new_atom = atom.replace('*',"'")
   if new_atom == " O1P":
     new_atom = " OP1"
   elif new_atom == " O2P":
     new_atom = " OP2"
   return new_atom

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    atom_name = atom.pdb_label_columns()[0:4]
    resname = atom.pdb_label_columns()[5:8]
    updated_resname = modernize_rna_resname(resname)
    if common_residue_names_get_class(updated_resname) == "common_rna_dna":
      updated_atom = modernize_rna_atom_name(atom=atom_name)
    else:
      updated_atom = atom_name
    key = updated_atom+atom.pdb_label_columns()[4:5]+\
          updated_resname+atom.pdb_label_columns()[8:]+\
          atom.segid
    i_seq_name_hash[atom.i_seq]=key
  return i_seq_name_hash

def build_i_seq_hash(pdb_hierarchy):
  name_i_seq_hash = dict()
  for atom in pdb_hierarchy.atoms():
    atom_name = atom.pdb_label_columns()[0:4]
    resname = atom.pdb_label_columns()[5:8]
    updated_resname = modernize_rna_resname(resname)
    if common_residue_names_get_class(updated_resname) == "common_rna_dna":
      updated_atom = modernize_rna_atom_name(atom=atom_name)
    else:
      updated_atom = atom_name
    key = updated_atom+atom.pdb_label_columns()[4:5]+\
          updated_resname+atom.pdb_label_columns()[8:]+\
          atom.segid
    name_i_seq_hash[key]=atom.i_seq
  return name_i_seq_hash

def id_str(chain_id,
            resseq,
            resname,
            icode,
            altloc,
            segid=None,
            ignore_altloc=False):
  base = "%2s%4s%1s" % (chain_id, resseq, icode)
  if (not ignore_altloc):
    base += "%1s" % altloc
  else :
    base += " "
  base += "%3s" % resname
  if (segid is not None):
    base += " segid='%4s'" % segid
  return base


def build_xyz_hash(pdb_hierarchy):
  name_xyz_hash = dict()
  for atom in pdb_hierarchy.atoms():
    name_xyz_hash[atom.pdb_label_columns()]=atom.xyz
  return name_xyz_hash

def build_resid_hash(pdb_hierarchy):
  resid_hash = dict()
  for rg in pdb_hierarchy.residue_groups():
    resid = rg.resseq_as_int()
    for atom in rg.atoms():
      resid_hash[atom.i_seq]=resid
  return resid_hash

def build_i_seq_xyz_hash(pdb_hierarchy):
  i_seq_xyz_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_xyz_hash[atom.i_seq] = atom.xyz
  return i_seq_xyz_hash

def build_element_hash(pdb_hierarchy):
  i_seq_element_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_element_hash[atom.i_seq]=atom.element
  return i_seq_element_hash

def build_chain_hash(pdb_hierarchy):
  chain_hash = dict()
  for chain in pdb_hierarchy.chains():
    for atom in chain.atoms():
      chain_hash[atom.i_seq] = chain.id
  return chain_hash

def build_segid_hash(pdb_hierarchy):
  segid_hash = dict()
  for atom in pdb_hierarchy.atoms():
    segid_hash[atom.i_seq] = atom.segid
  return segid_hash

def build_sym_atom_hash(pdb_hierarchy):
  sym_atom_hash = dict()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          if residue.resname.upper() in ['ASP', 'GLU', 'PHE', 'TYR']:
            if residue.resname.upper() == 'ASP':
              atom1 = ' OD1'
              atom2 = ' OD2'
            elif residue.resname.upper() == 'GLU':
              atom1 = ' OE1'
              atom2 = ' OE2'
            elif residue.resname.upper() in ['PHE', 'TYR']:
              atom1 = ' CD1'
              atom2 = ' CD2'
            atom1_i_seq = None
            atom2_i_seq = None
            for atom in residue.atoms():
              if atom.name == atom1:
                atom1_i_seq = atom.i_seq
              elif atom.name == atom2:
                atom2_i_seq = atom.i_seq
            if atom1_i_seq != None and atom2_i_seq != None:
              sym_atom_hash[atom1_i_seq] = atom2_i_seq
              sym_atom_hash[atom2_i_seq] = atom1_i_seq
  return sym_atom_hash

def angle_distance(angle1, angle2):
  distance = math.fabs(angle1 - angle2)
  if distance > 180.0:
    distance -= 360.0
  return math.fabs(distance)

def get_angle_average(angles):
  local_angles = []
  for angle in angles:
    if angle is not None:
      local_angles.append(angle)
  n_angles = len(local_angles)
  if n_angles < 1:
    return None
  sum = 0.0
  a1 = local_angles[0]
  if a1 > 180.0:
    a1 -= 360.0
  elif a1 < -180.0:
    a1 += 360.0
  sum += a1
  for angle in local_angles[1:]:
    a2 = angle
    if (a1 - a2) > 180.0:
      a2 += 360.0
    elif (a1 - a2) < -180.0:
      a2 -= 360.0
    sum += a2
  average = sum / n_angles
  return average

def check_for_internal_chain_ter_records(
      pdb_hierarchy,
      ter_indices):
  if ter_indices.size() == 0: return
  chains = pdb_hierarchy.chains()
  atoms = pdb_hierarchy.atoms()
  chain_ter_matches = {}
  chain_ranges = {}
  for chain in chains:
    if not chain.is_protein() and not chain.is_na():
      continue
    min = None
    max = None
    for atom in chain.atoms():
      if min is not None:
        if atom.i_seq < min:
          min = atom.i_seq
      else:
        min = atom.i_seq
      if max is not None:
        if atom.i_seq > max:
          max = atom.i_seq
      else:
        max = atom.i_seq
    if chain_ranges.get(chain.id) is None:
      chain_ranges[chain.id] = []
    chain_ranges[chain.id].append( (min, max) )

  #find min/max for all chains with same id
  reduced_chain_ranges = {}
  for key in chain_ranges.keys():
    min_all = None
    max_all = None
    ranges = chain_ranges[key]
    for min, max in ranges:
      if min_all is not None:
        if min < min_all:
          min_all = min
      else:
        min_all = min
      if max_all is not None:
        if max > max_all:
          max_all = max
      else:
        max_all = max
    reduced_chain_ranges[key] = (min_all, max_all)
  for ter_id in ter_indices:
    for key in reduced_chain_ranges.keys():
      min, max = reduced_chain_ranges[key]
      if ter_id > min and ter_id < max:
        raise Sorry("chain '%s' contains one or more "%(key)+
                    "errant TER cards.\nPlease remove and try again.")

def get_torsion_id(dp,
                   name_hash,
                   phi_psi=False,
                   chi_only=False,
                   omega=False):
  #
  # used in torsion_ncs
  id = None
  chi_atoms = False
  atom_list = []
  altloc = None
  if phi_psi:
    return name_hash[dp.i_seqs[1]][4:]
  elif omega:
    #LIMITATION: doesn't work with segIDs currently
    return name_hash[dp.i_seqs[0]][4:], \
           name_hash[dp.i_seqs[3]][4:]
  for i_seq in dp.i_seqs:
    cur_id = name_hash[i_seq][4:]
    atom = name_hash[i_seq][:4]
    atom_list.append(atom)
    cur_altloc = name_hash[i_seq][4:5]
    if id == None:
      id = cur_id
    if cur_altloc != " " and altloc:
      altloc = cur_altloc
    elif cur_id != id:
      return None
    resname = cur_id[1:4]
    if common_residue_names_get_class(resname,
         consider_ccp4_mon_lib_rna_dna=True) != "common_amino_acid":
      return None
    if chi_only:
      if atom not in [' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OXT']:
        chi_atoms = True
  if chi_only and not chi_atoms:
    return None
  return id

def get_c_alpha_hinges(pdb_hierarchy,
                       xray_structure=None,
                       selection=None):
  #
  # used in rotamer_search.py
  c_alphas = []
  c_alpha_hinges = {}
  if xray_structure is not None:
    sites_cart = xray_structure.sites_cart()
  else:
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
  if selection is None:
    selection = flex.bool(len(sites_cart), True)
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
            if( (not selection[cur_ca.i_seq]) or
                (not selection[cur_c.i_seq])  or
                (not selection[cur_n.i_seq])  or
                (not selection[cur_o.i_seq]) ):
              continue
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

def check_residues_are_connected(ca_1, ca_2, max_sep=4.0, min_sep=0.):
  from scitbx import matrix
  ca_1_mat = matrix.col(ca_1.xyz)
  ca_2_mat = matrix.col(ca_2.xyz)
  dist = (ca_1_mat - ca_2_mat).length()
  if (dist > max_sep) or (dist < min_sep):
    return False
  return True

def prepare_map(
      fmodel,
      exclude_free_r_reflections=False):
  map_obj = fmodel.electron_density_map()
  fft_map = map_obj.fft_map(resolution_factor = 1./4,
    map_type = "2mFo-DFc", use_all_data=(
      not exclude_free_r_reflections))
  fft_map.apply_sigma_scaling()
  target_map_data = fft_map.real_map_unpadded()
  fft_map = map_obj.fft_map(resolution_factor = 1./4,
    map_type = "mFo-DFc", use_all_data=(
      not exclude_free_r_reflections))
  fft_map.apply_sigma_scaling()
  residual_map_data = fft_map.real_map_unpadded()
  return target_map_data, residual_map_data
