from __future__ import absolute_import, division, print_function
import time
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
from rdkit import Chem
from rdkit.Chem import Lipinski, rdmolops
from scitbx.array_family import flex
import mmtbx.ligands.rdkit_utils as rdkit_utils

# ------------------------------------------------------------------------------

def to_py(obj):
  return [list(x) for x in obj]

def normalize(x):
  return sorted([tuple(sorted(v)) for v in x])

# ------------------------------------------------------------------------------

def run():
  run_test_01()
  run_test_02()

# ------------------------------------------------------------------------------

def run_test_01():
  pdb_inp = iotbx.pdb.input(lines=pdb_str_01.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  ligand_isel = model.iselection('resname EPE')
  rigid_comps_isels = run_fragmentation(ligand_isel, model)

  expected = [
    [0, 30, 31],
    [1, 2, 3, 4, 8, 9, 15, 16, 17, 18, 19, 20, 21, 22],
    [5, 23, 24],
    [6, 13, 25, 26, 29],
    [7, 27, 28],
    [10, 11, 12, 14]
    ]

  assert normalize(to_py(rigid_comps_isels)) == normalize(expected)

# ------------------------------------------------------------------------------

def run_test_02():
  pdb_inp = iotbx.pdb.input(lines=pdb_str_02.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  ligand_isel = model.iselection('resname NAG')
  rigid_comps_isels = run_fragmentation(ligand_isel, model)
  assert len(rigid_comps_isels) == 3

#  ph =  model.get_hierarchy()
#  atoms = ph.atoms()
#  for rigid_comp in rigid_comps_isels:
#    print('fragment')
#    print(list(rigid_comp))
#    for idx in rigid_comp:
#      print(atoms[idx].name)

# ------------------------------------------------------------------------------

def run_fragmentation(ligand_isel, model):
  ph =  model.get_hierarchy()
  atoms = ph.atoms()
  ph_ligand = ph.select(ligand_isel)
  atoms_ligand = ph_ligand.atoms()
  resname = atoms_ligand[0].parent().resname

  # Mappings
  # cctbx_idx -> rdkit_idx
  # rdkit_idx -> cctbx_idx
  # atom_name -> rdkit_idx (For CIF lookup)
  cctbx_to_rdkit = {}
  rdkit_to_cctbx = {}
  name_to_rdkit = {}

  mol = Chem.RWMol()

  # To find the CIF, we need the residue name.
  # We assume the selection belongs to one residue type.

  # --- 1. Build RDKit Nodes (Atoms) ---
  for i, atom_idx in enumerate(ligand_isel):
    cctbx_atom = atoms[atom_idx]

    # Handle Element
    element = cctbx_atom.element.strip().capitalize()
    #if not element: element = cctbx_atom.name.strip()[0:1]

    rd_atom = Chem.Atom(element)
    # Store name for debugging and potential mapping checks
    rd_atom.SetProp("_Name", cctbx_atom.name)

    rd_idx = mol.AddAtom(rd_atom)

    cctbx_to_rdkit[atom_idx] = rd_idx
    rdkit_to_cctbx[rd_idx] = atom_idx

    # Clean name for dictionary lookup (e.g. " C1 " -> "C1")
    name_to_rdkit[cctbx_atom.name.strip()] = rd_idx

  # 2. Build Bonds from CIF (Primary source)
  mon_lib_srv = model.get_mon_lib_srv()
  cif_object, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
    residue_name=resname, atom_names=atoms_ligand.extract_name())

  if cif_object and hasattr(cif_object, "bond_list"):
    for bond in cif_object.bond_list:
      atom_1 = bond.atom_id_1
      atom_2 = bond.atom_id_2
      order_str = bond.type

      if atom_1 in name_to_rdkit and atom_2 in name_to_rdkit:
        rd_i = name_to_rdkit[atom_1]
        rd_j = name_to_rdkit[atom_2]
        # Check existing to prevent duplicates
        if mol.GetBondBetweenAtoms(rd_i, rd_j) is None:
          mol.AddBond(rd_i, rd_j, rdkit_utils.get_rdkit_bond_type(order_str))

  # 4. Sanitize (Important for implicit valence calculation)
  try:
    Chem.SanitizeMol(mol)
  except ValueError:
    mol.UpdatePropertyCache(strict=False)

  # 5. Identify Rotatable Bonds
  rotatable_pattern = Lipinski.RotatableBondSmarts
  matches = mol.GetSubstructMatches(rotatable_pattern)

  bonds_to_cut = []
  rotatable_bonds = []
  min_heavy_atoms = 2
  for u, v in matches:
    bond = mol.GetBondBetweenAtoms(u, v)
    if bond is None: continue

    # Filter 1: Protect Amides
    if rdkit_utils.is_amide_bond(mol, bond): continue

    bidx = bond.GetIdx()

    # Cut this bond alone to test fragment sizes
    test_mol = rdmolops.FragmentOnBonds(mol, [bidx], addDummies=False)
    frags = Chem.GetMolFrags(test_mol, asMols=True, sanitizeFrags=False)

    # Count heavy atoms in each fragment
    heavy_counts = [
      sum(1 for a in frag.GetAtoms() if a.GetAtomicNum() > 1)
      for frag in frags
    ]

    # Keep bond only if BOTH sides are "real" fragments
    if all(hc >= min_heavy_atoms for hc in heavy_counts):
      bonds_to_cut.append(bidx)
      rotatable_bonds.append((u, v, bidx))

  # 6. Fragment
  if not bonds_to_cut:
    return [flex.size_t(iselection)]

  fragmented_mol = Chem.FragmentOnBonds(mol, bonds_to_cut, addDummies=False)
  raw_fragments = Chem.GetMolFrags(fragmented_mol, asMols=False) #maybe: sanitizeFrags=False?
  #frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=False)


  # 8. Convert to CCTBX format
  cctbx_rigid_components = []

  merged_fragments = raw_fragments
  for frag in merged_fragments:
    component_indices = flex.size_t()
    for rd_idx in frag:
      global_idx = rdkit_to_cctbx[rd_idx]
      component_indices.append(global_idx)
    cctbx_rigid_components.append(component_indices)

  return cctbx_rigid_components

# ------------------------------------------------------------------------------

pdb_str_01 = '''
CRYST1   14.227   18.379   15.820  90.00  90.00  90.00 P 1
HETATM    1  C10 EPE A 501       6.676   6.409   9.262  1.00 37.96           C
HETATM    2  C2  EPE A 501       7.070   8.694   7.451  1.00 35.02           C
HETATM    3  C3  EPE A 501       6.896  10.111   6.903  1.00 38.00           C
HETATM    4  C5  EPE A 501       5.000   9.315   5.530  1.00 26.48           C
HETATM    5  C6  EPE A 501       5.266   7.934   6.099  1.00 27.99           C
HETATM    6  C7  EPE A 501       6.110  11.478   5.000  1.00 29.40           C
HETATM    7  C8  EPE A 501       7.475  12.145   5.088  1.00 33.85           C
HETATM    8  C9  EPE A 501       6.026   6.573   7.884  1.00 35.27           C
HETATM    9  N1  EPE A 501       5.792   7.971   7.465  1.00 32.60           N
HETATM   10  N4  EPE A 501       6.147  10.182   5.652  1.00 31.52           N
HETATM   11  O1S EPE A 501       8.783   5.973  10.820  1.00 58.15           O
HETATM   12  O2S EPE A 501       9.038   5.000   8.637  1.00 44.87           O
HETATM   13  O3S EPE A 501       9.227   7.370   8.911  1.00 50.91           O
HETATM   14  O8  EPE A 501       7.411  13.379   5.779  1.00 38.77           O
HETATM   15  S   EPE A 501       8.500   6.171   9.381  1.00 61.74           S
HETATM   16  H21 EPE A 501       7.440   8.751   8.346  1.00 35.02           H
HETATM   17  H22 EPE A 501       7.721   8.227   6.904  1.00 35.02           H
HETATM   18  H31 EPE A 501       7.781  10.490   6.784  1.00 38.00           H
HETATM   19  H32 EPE A 501       6.453  10.637   7.587  1.00 38.00           H
HETATM   20  H51 EPE A 501       4.235   9.689   5.995  1.00 26.48           H
HETATM   21  H52 EPE A 501       4.744   9.212   4.600  1.00 26.48           H
HETATM   22  H61 EPE A 501       5.890   7.484   5.507  1.00 27.99           H
HETATM   23  H62 EPE A 501       4.434   7.436   6.069  1.00 27.99           H
HETATM   24  H71 EPE A 501       5.435  12.039   5.413  1.00 29.40           H
HETATM   25  H72 EPE A 501       5.846  11.376   4.072  1.00 29.40           H
HETATM   26  H81 EPE A 501       8.073  11.514   5.518  1.00 33.85           H
HETATM   27  H82 EPE A 501       7.800  12.249   4.180  1.00 33.85           H
HETATM   28  H91 EPE A 501       6.583   6.149   7.213  1.00 35.27           H
HETATM   29  H92 EPE A 501       5.173   6.111   7.873  1.00 35.27           H
HETATM   30  HO8 EPE A 501       7.233  13.215   6.594  1.00 38.77           H
HETATM   31 H101 EPE A 501       6.504   7.190   9.811  1.00 37.96           H
HETATM   32 H102 EPE A 501       6.308   5.637   9.719  1.00 37.96           H
END
'''

pdb_str_02 = '''
CRYST1   13.832   18.833   15.877  90.00  90.00  90.00 P 1
HETATM    1  C1  NAG A   1       7.363   8.528   6.805  1.00 20.00           C
HETATM    2  C2  NAG A   1       6.856   9.616   7.755  1.00 20.00           C
HETATM    3  C3  NAG A   1       7.487   9.443   9.135  1.00 20.00           C
HETATM    4  C4  NAG A   1       7.302   8.019   9.648  1.00 20.00           C
HETATM    5  C5  NAG A   1       7.837   7.027   8.616  1.00 20.00           C
HETATM    6  C6  NAG A   1       7.615   5.580   9.000  1.00 20.00           C
HETATM    7  C7  NAG A   1       6.197  11.818   6.804  1.00 20.00           C
HETATM    8  C8  NAG A   1       6.694  13.155   6.337  1.00 20.00           C
HETATM    9  N2  NAG A   1       7.121  10.952   7.253  1.00 20.00           N
HETATM   10  O1  NAG A   1       6.686   8.627   5.590  1.00 20.00           O
HETATM   11  O3  NAG A   1       6.915  10.382  10.035  1.00 20.00           O
HETATM   12  O4  NAG A   1       8.008   7.862  10.877  1.00 20.00           O
HETATM   13  O5  NAG A   1       7.166   7.233   7.355  1.00 20.00           O
HETATM   14  O6  NAG A   1       6.228   5.274   9.138  1.00 20.00           O
HETATM   15  O7  NAG A   1       5.000  11.544   6.771  1.00 20.00           O
HETATM   16  H1  NAG A   1       8.325   8.669   6.646  1.00 20.00           H
HETATM   17  H2  NAG A   1       5.884   9.484   7.863  1.00 20.00           H
HETATM   18  H3  NAG A   1       8.458   9.627   9.053  1.00 20.00           H
HETATM   19  H4  NAG A   1       6.339   7.846   9.801  1.00 20.00           H
HETATM   20  H5  NAG A   1       8.805   7.184   8.494  1.00 20.00           H
HETATM   21  H61 NAG A   1       8.072   5.394   9.847  1.00 20.00           H
HETATM   22  H62 NAG A   1       8.005   5.000   8.313  1.00 20.00           H
HETATM   23  H81 NAG A   1       7.658  13.134   6.224  1.00 20.00           H
HETATM   24  H82 NAG A   1       6.462  13.833   6.992  1.00 20.00           H
HETATM   25  H83 NAG A   1       6.278  13.375   5.488  1.00 20.00           H
HETATM   26  HN2 NAG A   1       7.957  11.213   7.247  1.00 20.00           H
HETATM   27  HO1 NAG A   1       7.112   8.194   5.000  1.00 20.00           H
HETATM   28  HO3 NAG A   1       7.377  10.385  10.737  1.00 20.00           H
HETATM   29  HO4 NAG A   1       8.832   7.992  10.763  1.00 20.00           H
HETATM   30  HO6 NAG A   1       5.839   5.379   8.392  1.00 20.00           H
END
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
