from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdFMCS
from collections import defaultdict
from rdkit.Chem import Lipinski
from scitbx.array_family import flex

"""
Utility functions to work with rdkit

Functions:
  convert_model_to_rdkit: Convert cctbx model to rdkit mol
  convert_elbow_to_rdkit: Convert elbow molecule to rdkit mol
  mol_to_3d: Generate 3D conformer of an rdkit mol
  mol_to_2d: Generate 2D conformer of an rdkit mol
  mol_from_smiles: Generate rdkit mol from smiles string
  match_mol_indices: Match atom indices of different mols
"""

#def get_rdkit_bond_type(cif_order):
#  """
#  Maps CCTBX CIF bond orders to RDKit BondType enum.
#  """
#  if not cif_order: return Chem.BondType.SINGLE
#  order = cif_order.lower()
#  if 'sing' in order: return Chem.BondType.SINGLE
#  if 'doub' in order: return Chem.BondType.DOUBLE
#  if 'trip' in order: return Chem.BondType.TRIPLE
#  if 'deloc' in order or 'arom' in order: return Chem.BondType.AROMATIC
#  return Chem.BondType.SINGLE


def get_rdkit_bond_type(cif_order, elements=None):
  """
  Maps CCTBX CIF bond orders to RDKit BondType enum.
  elements: A set/list of the two element symbols involved, e.g. {'C', 'O'}
  """
  if not cif_order: return Chem.BondType.SINGLE
  order = cif_order.lower()

  if 'sing' in order: return Chem.BondType.SINGLE
  if 'doub' in order: return Chem.BondType.DOUBLE
  if 'trip' in order: return Chem.BondType.TRIPLE
  if 'arom' in order: return Chem.BondType.AROMATIC

  if 'deloc' in order:
    # Special handling for P and S to avoid valence errors (e.g. Valence 6 on P)
    if elements and any(e in elements for e in ['P', 'S', 'CL']):
        return Chem.BondType.SINGLE

    # For Carbon (Carboxylates), 1.5 is safe and correct
    return Chem.BondType.ONEANDAHALF

  return Chem.BondType.SINGLE

def is_amide_bond(mol, bond):
  """
  Checks if a bond is a C-N amide bond to prevent cutting peptides.
  """
  a1 = bond.GetBeginAtom()
  a2 = bond.GetEndAtom()

  c_atom, n_atom = None, None
  if a1.GetSymbol() == 'C' and a2.GetSymbol() == 'N':
    c_atom, n_atom = a1, a2
  elif a2.GetSymbol() == 'C' and a1.GetSymbol() == 'N':
    c_atom, n_atom = a2, a1

  if c_atom is None: return False

  for nbr in c_atom.GetNeighbors():
    if nbr.GetIdx() == n_atom.GetIdx(): continue
    if nbr.GetSymbol() == 'O':
      bond_to_o = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
      if bond_to_o is not None and bond_to_o.GetBondType() == Chem.BondType.DOUBLE:
        return True
  return False

def get_rdkit_mol_from_atom_group_and_cif_obj(atom_group, cif_object):
  atoms_ligand = atom_group.atoms()

  # Mappings
  cctbx_to_rdkit = {}
  rdkit_to_cctbx = {}
  name_to_rdkit = {}

  # NEW: Map atom names to element strings (e.g., "CA" -> "C")
  name_to_element = {}

  mol = Chem.RWMol()

  # --- Build RDKit Nodes (Atoms) ---
  for cctbx_atom in atoms_ligand:
    atom_idx = cctbx_atom.i_seq

    # Get Element string (e.g., "C", "N", "P")
    element_str = cctbx_atom.element.strip().capitalize()
    if not element_str:
        # Fallback if element is empty
        element_str = cctbx_atom.name.strip()[0:1].capitalize()

    rd_atom = Chem.Atom(element_str)
    rd_atom.SetProp("_Name", cctbx_atom.name)

    rd_idx = mol.AddAtom(rd_atom)

    cctbx_to_rdkit[atom_idx] = rd_idx
    rdkit_to_cctbx[rd_idx] = atom_idx

    # Store mappings
    clean_name = cctbx_atom.name.strip()
    name_to_rdkit[clean_name] = rd_idx
    name_to_element[clean_name] = element_str.upper() # Store as "C", "P", etc.

  # --- Build Bonds from CIF ---
  if cif_object and hasattr(cif_object, "bond_list"):
    for bond in cif_object.bond_list:
      # These are STRINGS (names)
      atom_name_1 = bond.atom_id_1.strip()
      atom_name_2 = bond.atom_id_2.strip()

      # Use .type if available, otherwise assume single
      order_str = getattr(bond, 'type', 'sing')

      if atom_name_1 in name_to_rdkit and atom_name_2 in name_to_rdkit:
        rd_i = name_to_rdkit[atom_name_1]
        rd_j = name_to_rdkit[atom_name_2]

        if mol.GetBondBetweenAtoms(rd_i, rd_j) is None:
          # LOOKUP ELEMENTS from the map we built
          el1 = name_to_element.get(atom_name_1, 'C')
          el2 = name_to_element.get(atom_name_2, 'C')

          # Pass element set to the helper
          r_type = get_rdkit_bond_type(order_str, elements={el1, el2})

          mol.AddBond(rd_i, rd_j, r_type)

  # Apply charge fixes
  mol = fix_charges(mol)

  try:
    Chem.SanitizeMol(mol)
  except ValueError:
    print("Warning: Sanitization failed. Proceeding with unsanitized molecule.")
    mol.UpdatePropertyCache(strict=False)

  return mol, rdkit_to_cctbx

def fix_charges(mol):
  mol.UpdatePropertyCache(strict=False)
  for atom in mol.GetAtoms():
    anum = atom.GetAtomicNum()
    # GetValence(Explicit) replaces GetExplicitValence()
    val = atom.GetValence(Chem.ValenceType.EXPLICIT)
    charge = atom.GetFormalCharge()

    # Fix Nitrogen: 4 bonds, neutral -> +1
    if anum == 7 and val == 4 and charge == 0:
      atom.SetFormalCharge(1)

    # Fix Phosphorus/Sulfur:
    # Since we mapped 'deloc' -> SINGLE, we might have P with 4 single bonds.
    # Neutral P cannot have 4 bonds. P+ can.
    if anum == 15 and val == 4 and charge == 0:
      atom.SetFormalCharge(1)

  return mol

def get_cctbx_isel_for_rigid_components(atom_group,
                                        cif_object,
                                        filter_lone_linkers=True):
  mol, rdkit_to_cctbx = get_rdkit_mol_from_atom_group_and_cif_obj(
    atom_group = atom_group,
    cif_object = cif_object)
  cctbx_rigid_components = get_rigid_components(
    mol, rdkit_to_cctbx, filter_lone_linkers)

  return cctbx_rigid_components

def get_rigid_components(mol,
                         rdkit_to_cctbx,
                         filter_lone_linkers=True):

  # Identify Rotatable Bonds
  rotatable_pattern = Lipinski.RotatableBondSmarts
  matches = mol.GetSubstructMatches(rotatable_pattern)

  candidate_cut_bonds = []
  min_heavy_atoms = 2

  # Map: (bond_index, atom_index) -> Size of the fragment this atom ends up in
  # if bond is cut
  fragment_size_map = {}

  for u, v in matches:
    bond = mol.GetBondBetweenAtoms(u, v)
    if bond is None: continue

    if is_amide_bond(mol, bond): continue

    bidx = bond.GetIdx()

    # Cut this bond alone
    test_mol = Chem.FragmentOnBonds(mol, [bidx], addDummies=False)

    # Get indices of atoms in fragments
    frag_indices_tuples = Chem.GetMolFrags(test_mol, asMols=False)

    heavy_counts = []

    # Calculate heavy atom counts for this specific cut
    current_cut_sizes = {} # map atom_idx -> size

    for frag_tuple in frag_indices_tuples:
      # Count heavy atoms in this fragment
      h_count = 0
      for atom_idx in frag_tuple:
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() > 1:
          h_count += 1

      heavy_counts.append(h_count)

      # Store the size for every atom in this fragment
      for atom_idx in frag_tuple:
        current_cut_sizes[atom_idx] = h_count

    # Check validity
    if all(hc >= min_heavy_atoms for hc in heavy_counts):
      candidate_cut_bonds.append(bidx)
      # Save the size map for this bond index
      for atom_idx, size in current_cut_sizes.items():
        fragment_size_map[(bidx, atom_idx)] = size

  bonds_to_skip = set()
  if filter_lone_linkers:
    bonds_to_skip = _filter_isolated_linkers(candidate_cut_bonds, mol, fragment_size_map)

  # Finalize the list
  final_bonds_to_cut = [b for b in candidate_cut_bonds if b not in bonds_to_skip]

  # Fragment
  if not final_bonds_to_cut:
    return [flex.size_t(list(rdkit_to_cctbx.values()))]

  fragmented_mol = Chem.FragmentOnBonds(mol, final_bonds_to_cut, addDummies=False)
  raw_fragments = Chem.GetMolFrags(fragmented_mol, asMols=False)

  # Convert to CCTBX format
  cctbx_rigid_components = []

  merged_fragments = raw_fragments
  for frag in merged_fragments:
    component_indices = flex.size_t()
    for rd_idx in frag:
      global_idx = rdkit_to_cctbx[rd_idx]
      component_indices.append(global_idx)
    cctbx_rigid_components.append(component_indices)

  return cctbx_rigid_components

def _filter_isolated_linkers(candidate_cut_bonds, mol, fragment_size_map):
  # 1. Map bonds to atoms
  cuts_per_atom = defaultdict(list)
  for bidx in candidate_cut_bonds:
    bond = mol.GetBondWithIdx(bidx)
    cuts_per_atom[bond.GetBeginAtomIdx()].append(bidx)
    cuts_per_atom[bond.GetEndAtomIdx()].append(bidx)

  # 2. Identify Safe vs Unsafe atoms
  safe_atoms = set()
  unsafe_atoms = []

  for atom in mol.GetAtoms():
    if atom.GetAtomicNum() <= 1: continue

    idx = atom.GetIdx()
    heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    heavy_degree = len(heavy_neighbors)
    num_cuts = len(cuts_per_atom[idx])

    if num_cuts < heavy_degree:
      safe_atoms.add(idx)
    else:
      unsafe_atoms.append(idx)

  # --- SORTING ORDER ---
  # Process atoms adjacent to rings FIRST.
  # Otherwise, a chain neighbor might "steal" the linker before it attaches to
  # the ring.
  def priority_sort_key(atom_idx):
      atom = mol.GetAtomWithIdx(atom_idx)
      has_ring_neighbor = any(n.IsInRing() for n in atom.GetNeighbors())
      # Python sorts False (0) before True (1). We want True first, so we invert.
      return (not has_ring_neighbor, atom_idx)

  unsafe_atoms.sort(key=priority_sort_key)
  # ------------------------------

  bonds_to_skip = set()

  # 3. Rescue Unsafe Atoms
  for atom_idx in unsafe_atoms:
    if atom_idx in safe_atoms: continue

    atom = mol.GetAtomWithIdx(atom_idx)
    possible_bonds = cuts_per_atom[atom_idx]
    if not possible_bonds: continue

    best_bond_to_save = -1
    best_score = -float('inf')

    for bidx in possible_bonds:
      bond = mol.GetBondWithIdx(bidx)
      neighbor = bond.GetOtherAtom(atom)
      n_idx = neighbor.GetIdx()

      score = 0

      # CRITERIA 1: Ring Priority
      if neighbor.IsInRing(): score += 1000

      # CRITERIA 2: Pair Unsafe Atoms
      if n_idx not in safe_atoms:
        score += 50
      else:
        score -= 10

      # CRITERIA 3: Fragment Size
      if (bidx, n_idx) in fragment_size_map:
        score -= fragment_size_map[(bidx, n_idx)]

      # CRITERIA 4: Atomic Num (Tie-breaker)
      score += (0.1 * neighbor.GetAtomicNum())

      if score > best_score:
        best_score = score
        best_bond_to_save = bidx

    if best_bond_to_save != -1:
      bonds_to_skip.add(best_bond_to_save)
      safe_atoms.add(atom_idx)
      # Also mark neighbor as safe (connection established)
      bond = mol.GetBondWithIdx(best_bond_to_save)
      safe_atoms.add(bond.GetOtherAtom(atom).GetIdx())

  return bonds_to_skip

def get_prop_safe(rd_obj, prop):
  if prop not in rd_obj.GetPropNames(): return False
  return rd_obj.GetProp(prop)

def get_cc_cartesian_coordinates(cc_cif, label='pdbx_model_Cartn_x_ideal', ignore_question_mark=False):
  rc = []
  for i, (code, monomer) in enumerate(cc_cif.items()):
    atom = monomer.get_loop_or_row('_chem_comp_atom')
    # if atom is None: return rc
    for j, tmp in enumerate(atom.iterrows()):
      if label=='pdbx_model_Cartn_x_ideal':
        xyz = (tmp.get('_chem_comp_atom.pdbx_model_Cartn_x_ideal'),
               tmp.get('_chem_comp_atom.pdbx_model_Cartn_y_ideal'),
               tmp.get('_chem_comp_atom.pdbx_model_Cartn_z_ideal'),
               )
      elif label=='model_Cartn_x':
        xyz = (tmp.get('_chem_comp_atom.model_Cartn_x'),
               tmp.get('_chem_comp_atom.model_Cartn_y'),
               tmp.get('_chem_comp_atom.model_Cartn_z'),
               )
      rc.append(xyz)
      if not ignore_question_mark and '?' in xyz[-1]: return None
  # print(rc)
  return rc

def read_chemical_component_filename(filename):
  from iotbx import cif
  bond_order_ccd = {
    1.5:Chem.rdchem.BondType.AROMATIC,
    'SING': Chem.rdchem.BondType.SINGLE,
    'DOUB': Chem.rdchem.BondType.DOUBLE,
    'TRIP': Chem.rdchem.BondType.TRIPLE,
  }
  bond_order_rdkitkey = {value:key for key,value in bond_order_ccd.items()}
  ccd = cif.reader(filename).model()
  lookup={}
  def is_coordinates(x):
    return x!=('?', '?', '?')
  xyzs = get_cc_cartesian_coordinates(ccd, ignore_question_mark=True)
  xyzs = list(filter(is_coordinates, xyzs))
  if xyzs is None or len(xyzs)==0:
    xyzs = get_cc_cartesian_coordinates(ccd, label='model_Cartn_x', ignore_question_mark=True)
    xyzs = list(filter(is_coordinates, xyzs))
  if xyzs is None or len(xyzs)==0:
    for code, monomer in ccd.items():
      break
    raise Sorry('''
  Generating H restraints from Chemical Components for %s failed. Please supply
  restraints.
  ''' % code)
  for i, (code, monomer) in enumerate(ccd.items()):
    molecule = Chem.Mol()
    desc = monomer.get_loop_or_row('_chem_comp')
    rwmol = Chem.RWMol(molecule)
    atom = monomer.get_loop_or_row('_chem_comp_atom')
    # if atom is None: continue
    conformer = Chem.Conformer(atom.n_rows())
    for j, tmp in enumerate(atom.iterrows()):
      new = Chem.Atom(tmp.get('_chem_comp_atom.type_symbol').capitalize())
      new.SetFormalCharge(int(tmp.get('_chem_comp_atom.charge')))
      for prop in ['atom_id', 'type_symbol']:
        new.SetProp(prop, tmp.get('_chem_comp_atom.%s' % prop, '?'))
      rdatom = rwmol.AddAtom(new)
      if xyzs[j][0] in ['?']:
        pass
      else:
        xyz = (float(xyzs[j][0]), float(xyzs[j][1]), float(xyzs[j][2]))
        conformer.SetAtomPosition(rdatom, xyz)
      lookup[tmp.get('_chem_comp_atom.atom_id')]=j
    bond = monomer.get_loop_or_row('_chem_comp_bond')
    if bond:
      for tmp in bond.iterrows():
        atom1 = tmp.get('_chem_comp_bond.atom_id_1')
        atom2 = tmp.get('_chem_comp_bond.atom_id_2')
        atom1 = lookup.get(atom1)
        atom2 = lookup.get(atom2)
        order = tmp.get('_chem_comp_bond.value_order')
        order = bond_order_ccd[order]
        rwmol.AddBond(atom1, atom2, order)
  rwmol.AddConformer(conformer)
  # Chem.SanitizeMol(rwmol)
  # from rdkit.Chem.PropertyMol import PropertyMol
  molecule = rwmol.GetMol()
  # molecule = PropertyMol(molecule)
#  print(dir(molecule))
#  print(desc)
#  print(dir(desc))
#  for key, item in desc.items():
#    key = key.split('.')[1]
#    print(key,list(item))
#    molecule.SetProp(key,item[0])
#    print(molecule.HasProp(key))
#    print(molecule.GetProp(key))
#  print(dir(molecule.GetPropNames()))
#  print(molecule.GetPropsAsDict())
  # print(dir(rwmol))
  return molecule

def get_molecule_from_resname(resname):
  import os
  from mmtbx.chemical_components import get_cif_filename
  filename = get_cif_filename(resname)
  if not os.path.exists(filename): return None
  try:
    molecule = read_chemical_component_filename(filename)
  except Exception as e:
    print(e)
    return None
  return molecule

def mol_from_chemical_component(code):
  from mmtbx.chemical_components import get_cif_filename
  rc = get_cif_filename(code)
  molecule = read_chemical_component_filename(rc)
  return molecule

def convert_model_to_rdkit(cctbx_model):
  """
  Convert a cctbx model molecule object to an
  rdkit molecule object

  TODO: Bond type is always unspecified
  """
  assert cctbx_model.restraints_manager is not None, "Restraints manager must be set"

  mol = Chem.Mol()
  rwmol = Chem.RWMol(mol)
  conformer = Chem.Conformer(cctbx_model.get_number_of_atoms())

  for i,atom in enumerate(cctbx_model.get_atoms()):
    element = atom.element.strip().upper()
    if element =="D":
      element = "H"
    else:
      element = element
    atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(element)
    rdatom = Chem.Atom(atomic_number)
    rdatom.SetFormalCharge(atom.charge_as_int())
    rdatom_idx = rwmol.AddAtom(rdatom)
    conformer.SetAtomPosition(rdatom_idx,atom.xyz)

  rm = cctbx_model.restraints_manager
  grm = rm.geometry
  bonds_simple, bonds_asu = grm.get_all_bond_proxies()
  bond_proxies = bonds_simple.get_proxies_with_origin_id()
  for bond_proxy in bond_proxies:
    begin, end = bond_proxy.i_seqs
    order = Chem.rdchem.BondType.UNSPECIFIED
    rwmol.AddBond(int(begin),int(end),order)

  rwmol.AddConformer(conformer)
  mol = rwmol.GetMol()
  return mol

def convert_elbow_to_rdkit(elbow_mol):
  """
  Convert elbow molecule object to an
  rdkit molecule object

  TODO: Charge
  """
  # elbow bond order to rdkit bond orders
  bond_order_elbowkey = {
    1.5:Chem.rdchem.BondType.AROMATIC,
    1: Chem.rdchem.BondType.SINGLE,
    2: Chem.rdchem.BondType.DOUBLE,
    3: Chem.rdchem.BondType.TRIPLE,
  }
  bond_order_rdkitkey = {value:key for key,value in bond_order_elbowkey.items()}
  atoms = list(elbow_mol)
  mol = Chem.Mol()
  rwmol = Chem.RWMol(mol)
  conformer = Chem.Conformer(len(atoms))

  for i,atom in enumerate(atoms):
    xyz = atom.xyz
    atomic_number = atom.number
    rdatom = rwmol.AddAtom(Chem.Atom(int(atomic_number)))
    conformer.SetAtomPosition(rdatom,xyz)

  for i,bond in enumerate(elbow_mol.bonds):
    bond_atoms = list(bond)
    start,end = atoms.index(bond_atoms[0]), atoms.index(bond_atoms[1])
    order = bond_order_elbowkey[bond.order]
    rwmol.AddBond(int(start),int(end),order)

  rwmol.AddConformer(conformer)
  mol = rwmol.GetMol()
  return mol

def convert_rdkit_to_elbow(rwmol):
  from elbow.chemistry.SimpleMoleculeClass import SimpleMoleculeClass
  from elbow.chemistry.xyzClass import xyzClass
  positions = molecule.GetConformer().GetPositions()
  smc = SimpleMoleculeClass()
  for i, atom in enumerate(smc):
    atom.xyz = xyzClass(positions[i])
    atom.record_name = 'LIG'
    atom.chainID = 'A'
    atom.segID = ''
  smc.SetOriginalFormat('PDB')
  assert 0

def enumerate_bonds(mol):
  idx_set_bonds = {frozenset((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())) for bond in mol.GetBonds()}
  # check that the above approach matches the more exhaustive approach used for angles/torsion
  idx_set = set()
  for atom in mol.GetAtoms():
    for neigh1 in atom.GetNeighbors():
      idx0,idx1 = atom.GetIdx(), neigh1.GetIdx()
      s = frozenset([idx0,idx1])
      if len(s)==2:
        if idx0>idx1:
            idx0,idx1 = idx1,idx0
            idx_set.add(s)
  assert idx_set == idx_set_bonds
  return idx_set_bonds

def enumerate_angles(mol):
  idx_set = set()
  for atom in mol.GetAtoms():
    for neigh1 in atom.GetNeighbors():
      for neigh2 in neigh1.GetNeighbors():
        idx0,idx1,idx2 = atom.GetIdx(), neigh1.GetIdx(),neigh2.GetIdx()
        s = (idx0,idx1,idx2)
        if len(set(s))==3:
          if idx0>idx2:
            idx0,idx2 = idx2,idx0
          idx_set.add((idx0,idx1,idx2))
  return idx_set

def enumerate_torsions(mol):
  idx_set = set()
  for atom0 in mol.GetAtoms():
    idx0 = atom0.GetIdx()
    for atom1 in atom0.GetNeighbors():
      idx1 = atom1.GetIdx()
      for atom2 in atom1.GetNeighbors():
        idx2 = atom2.GetIdx()
        if idx2==idx0:
          continue
        for atom3 in atom2.GetNeighbors():
          idx3 = atom3.GetIdx()
          if idx3 == idx1 or idx3 == idx0:
            continue
          s = (idx0,idx1,idx2,idx3)
          if len(set(s))==4:
            if idx0<idx3:
              idx_set.add((idx0,idx1,idx2,idx3))
            else:
              idx_set.add((idx3,idx2,idx1,idx0))
  return idx_set

def mol_to_3d(mol):
  """
  Convert and rdkit mol to 3D coordinates
  """
  assert len(mol.GetConformers())==0, "mol already has conformer"
  param = rdDistGeom.ETKDGv3()
  conf_id = rdDistGeom.EmbedMolecule(mol,clearConfs=True)
  return mol

def mol_to_2d(mol):
  """
  Convert and rdkit mol to 2D coordinates
  """
  mol = Chem.Mol(mol) # copy to preserve original coords
  ret = Chem.rdDepictor.Compute2DCoords(mol)
  return mol

def mol_from_smiles(smiles, embed3d=False, addHs=True, removeHs=False, verbose=False):
  """
  Convert a smiles string to rdkit mol
  """
  ps = Chem.SmilesParserParams()
  ps.removeHs=removeHs
  rdmol = Chem.MolFromSmiles(smiles, ps)
  if verbose: print('rdmol',rdmol)
  if rdmol is None: return rdmol
  if verbose: print('rdmol',rdmol.Debug())
  if addHs: rdmol = Chem.AddHs(rdmol)
  if verbose: print('rdmol',rdmol.Debug())
  if embed3d: rdmol = mol_to_3d(rdmol)
  if verbose: print('rdmol',rdmol.Debug())
  if removeHs: rdmol = Chem.RemoveHs(rdmol)
  if verbose: print('rdmol',rdmol.Debug())
  Chem.SetHybridization(rdmol)
  if verbose: print('rdmol',rdmol.Debug())
  rdmol.UpdatePropertyCache()
  if verbose: print('rdmol',rdmol.Debug())
  return rdmol

def match_mol_indices(mol_list):
  """
  Match atom indices of molecules.

  Args:
      mol_list (list): a list of rdkit mols

  Returns:
      match_list: (list): a list of tuples
                          Each entry is a match beween in mols
                          Each value is the atom index for each mol
  """
  mol_list = [Chem.Mol(mol) for mol in mol_list]
  mcs_SMARTS = rdFMCS.FindMCS(mol_list)
  smarts_mol = Chem.MolFromSmarts(mcs_SMARTS.smartsString)
  match_list = [x.GetSubstructMatch(smarts_mol) for x in mol_list]
  return list(zip(*match_list))

def is_amino_acid(molecule):
  atom_names = ['N', 'CA', 'C', 'O']
  bond_names = [['CA', 'N'],
                ['C', 'O'],
                ['C', 'CA'],
               ]
  acount=0
  for atom in molecule.GetAtoms():
    if get_prop_safe(atom, 'atom_id') in atom_names:
      acount+=1
  bcount=0
  for bond in molecule.GetBonds():
    names = [get_prop_safe(bond.GetBeginAtom(), 'atom_id'),
             get_prop_safe(bond.GetEndAtom(), 'atom_id'),
            ]
    names.sort()
    if names in bond_names:
      bcount+=1
  if acount==4 and bcount==3:
    return True
  return False

def is_nucleic_acid(molecule):
  atom_names = ['P', "'O5'", 'OP1', 'OP2',
                "O3'", "C3'"]
  bond_names = [["O5'", 'P'],
                ['OP1', 'P'],
                ['OP2', 'P'],
                ["C3'", "O3'"],
               ]
  acount=0
  for atom in molecule.GetAtoms():
    if get_prop_safe(atom, 'atom_id') in atom_names:
      acount+=1
  bcount=0
  for bond in molecule.GetBonds():
    names = [get_prop_safe(bond.GetBeginAtom(), 'atom_id'),
             get_prop_safe(bond.GetEndAtom(), 'atom_id'),
            ]
    names.sort()
    if names in bond_names:
      bcount+=1
  if acount==5 and bcount==4:
    return True
  return False

if __name__ == '__main__':
  import sys
  if sys.argv[1:]:
    read_chemical_component_filename(sys.argv[1])
  else:
    for smiles_string in ['CD',
                          'Cd',
                          '[Cd]',
                          'CC',
                          'c1ccc1',
      ]:
      mol = mol_from_smiles(smiles_string, verbose=True)
      try:
        print(mol.Debug())
      except Exception: pass
    assert 0
