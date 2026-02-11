from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdFMCS
from collections import defaultdict

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

def get_rdkit_bond_type(cif_order):
  """
  Maps CCTBX CIF bond orders to RDKit BondType enum.
  """
  if not cif_order: return Chem.BondType.SINGLE
  order = cif_order.lower()
  if 'sing' in order: return Chem.BondType.SINGLE
  if 'doub' in order: return Chem.BondType.DOUBLE
  if 'trip' in order: return Chem.BondType.TRIPLE
  if 'deloc' in order or 'arom' in order: return Chem.BondType.AROMATIC
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

def _filter_isolated_linkers(candidate_cut_bonds, mol, fragment_size_map):
  """
  Helper function for ligand fragmentation
  If bonds are cut when there is a rotational degree of freedom, the number of
  logical fragments may be too large. An example is HEX, which has 5 such bonds.
  However, one can consider HEX as two fragments - with 3 carbons each - fused
  together. Accordingly, this filter will cut HEX into two fragments.
  The filter is optional in case one really wants to split at all applicable
  bonds.
  """
  cuts_per_atom = defaultdict(list)
  for bidx in candidate_cut_bonds:
    bond = mol.GetBondWithIdx(bidx)
    cuts_per_atom[bond.GetBeginAtomIdx()].append(bidx)
    cuts_per_atom[bond.GetEndAtomIdx()].append(bidx)

  bonds_to_skip = set()

  for atom_idx, bond_indices in cuts_per_atom.items():
    # Only process atoms being cut from multiple sides
    if len(bond_indices) < 2: continue

    atom = mol.GetAtomWithIdx(atom_idx)

    heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    heavy_degree = len(heavy_neighbors)

    # If we are about to completely isolate this linker atom
    if heavy_degree >= 2 and len(bond_indices) == heavy_degree:

      best_bond_to_save = -1
      # Initialize with a very low number
      best_neighbor_score = -999999

      for bidx in bond_indices:
        bond = mol.GetBondWithIdx(bidx)
        neighbor = bond.GetOtherAtom(atom)
        neighbor_idx = neighbor.GetIdx()

        score = 0

        # --- CRITERIA 1: RINGS (Highest Priority) ---
        # Always stick to the ring if possible.
        if neighbor.IsInRing():
          score += 10000

        # --- CRITERIA 2: FRAGMENT SIZE (Symmetry Breaker) ---
        # Look up the size of the fragment the neighbor would belong to
        # if we proceeded with this cut.
        # We prefer saving bonds to SMALLER fragments (ends of chains).
        # We subtract the size, so smaller size = higher score.
        if (bidx, neighbor_idx) in fragment_size_map:
          frag_size = fragment_size_map[(bidx, neighbor_idx)]
          score -= (frag_size * 100)

        # --- CRITERIA 3: ATOM PROPERTIES (Tie Breaker) ---
        score += neighbor.GetAtomicNum()

        if score > best_neighbor_score:
          best_neighbor_score = score
          best_bond_to_save = bidx

      if best_bond_to_save != -1:
        bonds_to_skip.add(best_bond_to_save)

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
