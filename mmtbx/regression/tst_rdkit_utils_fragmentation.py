from __future__ import absolute_import, division, print_function
import time
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
from rdkit import Chem
from rdkit.Chem import Lipinski
from scitbx.array_family import flex
import mmtbx.ligands.rdkit_utils as rdkit_utils
from collections import defaultdict

# ------------------------------------------------------------------------------
# Helper functions

def to_py(obj):
  return [list(x) for x in obj]

def normalize(x):
  return sorted([tuple(sorted(v)) for v in x])

# ------------------------------------------------------------------------------
# Tests to be executed

def run():
  run_test_01()
  run_test_02()
  run_test_03()
  run_test_04()
  run_test_05()
  run_test_06()

# ------------------------------------------------------------------------------

def run_test_01():
  print('test01...')
  rigid_comps_isels = compute_fragments(pdb_str_01, 'resname EPE', 3)

  expected = [
    [0, 10, 11, 12, 14, 30, 31],
    [1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 27, 28],
    [6, 13, 25, 26, 29]
    ]

  assert normalize(to_py(rigid_comps_isels)) == normalize(expected)

  # again without filtering for lone linkers

  rigid_comps_isels = compute_fragments(pdb_str_01, 'resname EPE', 6, filter_lone_linkers=False)

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
  print('test02...')
  rigid_comps_isels = compute_fragments(pdb_str_02, 'resname NAG', 3)

# ------------------------------------------------------------------------------

def run_test_03():
  print('test03...')
  rigid_comps_isels = compute_fragments(pdb_str_03, 'resname BNL', 2)

# ------------------------------------------------------------------------------

def run_test_04():
  print('test04...')
  rigid_comps_isels = compute_fragments(pdb_str_04, 'resname STI', 6)

# ------------------------------------------------------------------------------

def run_test_05():
  print('test05...')
  rigid_comps_isels = compute_fragments(pdb_str_05, 'resname HEX', 2)

  expected = [
    [0, 1, 2, 6, 7, 8, 9, 10, 11, 12],
    [3, 4, 5, 13, 14, 15, 16, 17, 18, 19]
    ]

  assert normalize(to_py(rigid_comps_isels)) == normalize(expected)

# ------------------------------------------------------------------------------

def run_test_06():
  print('test06...')
  pdb_inp = iotbx.pdb.input(lines=cif_str_06.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = cif_restraints_06).model()
  # bla.cif does not exist, but cif_objects needs a filename in first position
  # of the tuple
  cif_objects = [('bla.cif', cif_object)]
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    restraint_objects = cif_objects,
    log = null_out())
  model.process(make_restraints=True)
  ligand_isel = model.iselection('resname A1EF7')
  rigid_comps_isels = run_fragmentation(ligand_isel, model)
  #print(len(rigid_comps_isels))

  assert len(rigid_comps_isels) == 4

# ------------------------------------------------------------------------------

def compute_fragments(pdb_str, sel_str, expected, filter_lone_linkers=True):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    stop_for_unknowns = False,
    log         = null_out())
  model.process(make_restraints=True)
  ligand_isel = model.iselection(sel_str)
  rigid_comps_isels = run_fragmentation(ligand_isel, model, filter_lone_linkers)
  print(len(rigid_comps_isels))
  assert len(rigid_comps_isels) == expected

#  ph =  model.get_hierarchy()
#  atoms = ph.atoms()
#  for rigid_comp in rigid_comps_isels:
#    print('fragment')
#    print(list(rigid_comp))
#    for idx in rigid_comp:
#      print(atoms[idx].name)

  return rigid_comps_isels

# ------------------------------------------------------------------------------

def get_mol_from_iselection(iselection, mon_lib_srv, atoms):

  atoms_ligand = atoms.select(iselection)
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

  # --- Build RDKit Nodes (Atoms) ---
  for i, atom_idx in enumerate(iselection):
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

  # Build Bonds from CIF (Primary source)
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

  # Sanitize (Important for implicit valence calculation)
  try:
    Chem.SanitizeMol(mol)
  except ValueError:
    mol.UpdatePropertyCache(strict=False)

  return mol, rdkit_to_cctbx

# ------------------------------------------------------------------------------

def run_fragmentation(ligand_isel, model, filter_lone_linkers=True):
  ph =  model.get_hierarchy()
  atoms = ph.atoms()
  #ph_ligand = ph.select(ligand_isel)
  #atoms_ligand = ph_ligand.atoms()
  #resname = atoms_ligand[0].parent().resname

  mon_lib_srv = model.get_mon_lib_srv()
  mol, rdkit_to_cctbx = get_mol_from_iselection(ligand_isel, mon_lib_srv, atoms)

  # Identify Rotatable Bonds
  rotatable_pattern = Lipinski.RotatableBondSmarts
  matches = mol.GetSubstructMatches(rotatable_pattern)

  #bonds_to_cut = []
  #rotatable_bonds = []
  candidate_cut_bonds = []
  min_heavy_atoms = 2

 # Map: (bond_index, atom_index) -> Size of the fragment this atom ends up in if bond is cut
  fragment_size_map = {}

  for u, v in matches:
    bond = mol.GetBondBetweenAtoms(u, v)
    if bond is None: continue

    if rdkit_utils.is_amide_bond(mol, bond): continue

    bidx = bond.GetIdx()

    # Cut this bond alone
    test_mol = Chem.FragmentOnBonds(mol, [bidx], addDummies=False)

    # Get indices of atoms in fragments (asMols=False returns tuple of tuples of indices)
    frag_indices_tuples = Chem.GetMolFrags(test_mol, asMols=False)

    heavy_counts = []
    is_valid_cut = True

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
    bonds_to_skip = _linker_isolation_filter(candidate_cut_bonds, mol, fragment_size_map)

  # Finalize the list
  final_bonds_to_cut = [b for b in candidate_cut_bonds if b not in bonds_to_skip]

  # Fragment
  if not final_bonds_to_cut:
    return [flex.size_t(iselection)]

  fragmented_mol = Chem.FragmentOnBonds(mol, final_bonds_to_cut, addDummies=False)
  raw_fragments = Chem.GetMolFrags(fragmented_mol, asMols=False) #maybe: sanitizeFrags=False?
  #frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=False)

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


def _linker_isolation_filter(candidate_cut_bonds, mol, fragment_size_map):
  # 2. Filter Candidates to Prevent Linker Isolation
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

pdb_str_03 = '''
REMARK BIPHENYL BNL
REMARK This compound has one rotatable bond
CRYST1   15.143   14.093   17.907  90.00  90.00  90.00 P 1
HETATM    1  C1  BNL A   1       8.751   6.093   8.064  1.00 20.00           C
HETATM    2  C12 BNL A   1       5.584   8.353  10.905  1.00 20.00           C
HETATM    3  C13 BNL A   1       5.991   7.741  12.073  1.00 20.00           C
HETATM    4  C14 BNL A   1       6.847   6.659  12.019  1.00 20.00           C
HETATM    5  C15 BNL A   1       7.297   6.188  10.791  1.00 20.00           C
HETATM    6  C16 BNL A   1       6.895   6.790   9.597  1.00 20.00           C
HETATM    7  C17 BNL A   1       6.034   7.887   9.676  1.00 20.00           C
HETATM    8  C2  BNL A   1       7.386   6.284   8.287  1.00 20.00           C
HETATM    9  C3  BNL A   1       6.506   5.976   7.246  1.00 20.00           C
HETATM   10  C4  BNL A   1       6.977   5.502   6.028  1.00 20.00           C
HETATM   11  C5  BNL A   1       8.330   5.323   5.829  1.00 20.00           C
HETATM   12  C6  BNL A   1       9.217   5.617   6.844  1.00 20.00           C
HETATM   13  H1  BNL A   1       9.366   6.292   8.753  1.00 20.00           H
HETATM   14  H3  BNL A   1       5.577   6.095   7.372  1.00 20.00           H
HETATM   15  H12 BNL A   1       5.000   9.093  10.940  1.00 20.00           H
HETATM   16  H13 BNL A   1       5.685   8.060  12.907  1.00 20.00           H
HETATM   17  H14 BNL A   1       7.127   6.238  12.815  1.00 20.00           H
HETATM   18  H15 BNL A   1       7.883   5.447  10.762  1.00 20.00           H
HETATM   19  H17 BNL A   1       5.748   8.312   8.881  1.00 20.00           H
HETATM   20  H4  BNL A   1       6.368   5.302   5.336  1.00 20.00           H
HETATM   21  H5  BNL A   1       8.648   5.000   5.000  1.00 20.00           H
HETATM   22  H6  BNL A   1      10.143   5.496   6.711  1.00 20.00           H
END
'''

pdb_str_04 = '''
REMARK keep amide intact, break into rings, but prevent fragmenting into
REMARK too many single HA components
CRYST1   19.627   25.695   26.928  90.00  90.00  90.00 P 1
HETATM    1  C1  STI A   1      13.077  18.057   5.775  1.00 20.00           C
HETATM    2  C11 STI A   1       8.283  17.460   8.920  1.00 20.00           C
HETATM    3  C12 STI A   1       9.302  17.708   8.018  1.00 20.00           C
HETATM    4  C14 STI A   1       8.844  18.064  13.083  1.00 20.00           C
HETATM    5  C15 STI A   1       9.057  16.981  13.935  1.00 20.00           C
HETATM    6  C16 STI A   1       8.171  16.749  14.983  1.00 20.00           C
HETATM    7  C17 STI A   1       7.079  17.594  15.171  1.00 20.00           C
HETATM    8  C18 STI A   1       6.877  18.669  14.319  1.00 20.00           C
HETATM    9  C19 STI A   1       7.752  18.928  13.261  1.00 20.00           C
HETATM   10  C2  STI A   1      13.874  19.113   6.168  1.00 20.00           C
HETATM   11  C20 STI A   1       7.510  20.100  12.352  1.00 20.00           C
HETATM   12  C22 STI A   1       8.333  14.338  15.503  1.00 20.00           C
HETATM   13  C23 STI A   1       8.502  13.348  16.609  1.00 20.00           C
HETATM   14  C25 STI A   1       9.522  12.407  16.540  1.00 20.00           C
HETATM   15  C26 STI A   1       9.689  11.472  17.551  1.00 20.00           C
HETATM   16  C27 STI A   1       8.841  11.460  18.656  1.00 20.00           C
HETATM   17  C28 STI A   1       7.811  12.395  18.713  1.00 20.00           C
HETATM   18  C29 STI A   1       7.644  13.330  17.702  1.00 20.00           C
HETATM   19  C4  STI A   1      12.556  19.520   7.987  1.00 20.00           C
HETATM   20  C46 STI A   1       9.008  10.432  19.749  1.00 20.00           C
HETATM   21  C49 STI A   1       7.849   8.735  21.058  1.00 20.00           C
HETATM   22  C5  STI A   1      11.689  18.473   7.683  1.00 20.00           C
HETATM   23  C50 STI A   1       6.666   7.795  21.080  1.00 20.00           C
HETATM   24  C52 STI A   1       6.817   7.554  18.703  1.00 20.00           C
HETATM   25  C53 STI A   1       8.008   8.483  18.687  1.00 20.00           C
HETATM   26  C54 STI A   1       5.590   5.911  19.994  1.00 20.00           C
HETATM   27  C6  STI A   1      11.976  17.737   6.540  1.00 20.00           C
HETATM   28  C7  STI A   1      10.505  18.168   8.539  1.00 20.00           C
HETATM   29  C9  STI A   1       9.607  18.090  10.671  1.00 20.00           C
HETATM   30  N10 STI A   1       8.401  17.640  10.236  1.00 20.00           N
HETATM   31  N13 STI A   1       9.754  18.292  12.020  1.00 20.00           N
HETATM   32  N21 STI A   1       8.380  15.650  15.844  1.00 20.00           N
HETATM   33  N3  STI A   1      13.631  19.849   7.262  1.00 20.00           N
HETATM   34  N48 STI A   1       7.912   9.454  19.782  1.00 20.00           N
HETATM   35  N51 STI A   1       6.726   6.840  19.974  1.00 20.00           N
HETATM   36  N8  STI A   1      10.657  18.355   9.864  1.00 20.00           N
HETATM   37  O29 STI A   1       8.158  13.952  14.352  1.00 20.00           O
HETATM   38  H11 STI A   1      13.282  17.563   5.000  1.00 20.00           H
HETATM   39  H21 STI A   1      14.627  19.330   5.645  1.00 20.00           H
HETATM   40  H41 STI A   1      12.379  20.033   8.754  1.00 20.00           H
HETATM   41  H61 STI A   1      11.419  17.017   6.287  1.00 20.00           H
HETATM   42 H111 STI A   1       7.452  17.145   8.583  1.00 20.00           H
HETATM   43 H121 STI A   1       9.186  17.570   7.100  1.00 20.00           H
HETATM   44 H131 STI A   1      10.534  18.611  12.259  1.00 20.00           H
HETATM   45 H151 STI A   1       9.802  16.409  13.797  1.00 20.00           H
HETATM   46 H171 STI A   1       6.479  17.435  15.880  1.00 20.00           H
HETATM   47 H181 STI A   1       6.136  19.236  14.458  1.00 20.00           H
HETATM   48 H201 STI A   1       6.854  20.695  12.751  1.00 20.00           H
HETATM   49 H202 STI A   1       7.178  19.785  11.496  1.00 20.00           H
HETATM   50 H203 STI A   1       8.341  20.585  12.217  1.00 20.00           H
HETATM   51 H211 STI A   1       8.562  15.848  16.679  1.00 20.00           H
HETATM   52 H251 STI A   1      10.108  12.405  15.804  1.00 20.00           H
HETATM   53 H261 STI A   1      10.389  10.842  17.492  1.00 20.00           H
HETATM   54 H281 STI A   1       7.223  12.397  19.451  1.00 20.00           H
HETATM   55 H291 STI A   1       6.947  13.959  17.760  1.00 20.00           H
HETATM   56 H461 STI A   1       9.053  10.897  20.613  1.00 20.00           H
HETATM   57 H462 STI A   1       9.858   9.958  19.622  1.00 20.00           H
HETATM   58 H491 STI A   1       7.769   9.375  21.793  1.00 20.00           H
HETATM   59 H492 STI A   1       8.674   8.225  21.186  1.00 20.00           H
HETATM   60 H501 STI A   1       6.654   7.310  21.928  1.00 20.00           H
HETATM   61 H502 STI A   1       5.839   8.315  21.016  1.00 20.00           H
HETATM   62 H521 STI A   1       6.900   6.909  17.972  1.00 20.00           H
HETATM   63 H522 STI A   1       6.000   8.072  18.560  1.00 20.00           H
HETATM   64 H531 STI A   1       8.041   8.950  17.829  1.00 20.00           H
HETATM   65 H532 STI A   1       8.830   7.962  18.783  1.00 20.00           H
HETATM   66 H541 STI A   1       5.000   6.116  20.744  1.00 20.00           H
HETATM   67 H542 STI A   1       5.923   5.000  20.091  1.00 20.00           H
HETATM   68 H543 STI A   1       5.084   5.979  19.163  1.00 20.00           H
END
'''

pdb_str_05 = '''
CRYST1   12.668   16.763   14.843  90.00  90.00  90.00 P 1
HETATM    1  C1  HEX A   1       6.725  11.189   5.766  1.00 20.00           C
HETATM    2  C2  HEX A   1       6.047   9.866   6.032  1.00 20.00           C
HETATM    3  C3  HEX A   1       6.621   9.117   7.207  1.00 20.00           C
HETATM    4  C4  HEX A   1       5.954   7.793   7.486  1.00 20.00           C
HETATM    5  C5  HEX A   1       6.529   7.043   8.662  1.00 20.00           C
HETATM    6  C6  HEX A   1       5.863   5.713   8.946  1.00 20.00           C
HETATM    7  H11 HEX A   1       6.309  11.620   5.000  1.00 20.00           H
HETATM    8  H12 HEX A   1       7.668  11.039   5.580  1.00 20.00           H
HETATM    9  H13 HEX A   1       6.637  11.763   6.547  1.00 20.00           H
HETATM   10  H21 HEX A   1       6.123   9.303   5.231  1.00 20.00           H
HETATM   11  H22 HEX A   1       5.092  10.026   6.196  1.00 20.00           H
HETATM   12  H31 HEX A   1       6.546   9.682   8.006  1.00 20.00           H
HETATM   13  H32 HEX A   1       7.575   8.960   7.042  1.00 20.00           H
HETATM   14  H41 HEX A   1       6.028   7.229   6.687  1.00 20.00           H
HETATM   15  H42 HEX A   1       5.000   7.950   7.651  1.00 20.00           H
HETATM   16  H51 HEX A   1       6.455   7.609   9.461  1.00 20.00           H
HETATM   17  H52 HEX A   1       7.484   6.887   8.497  1.00 20.00           H
HETATM   18  H61 HEX A   1       6.522   5.000   8.883  1.00 20.00           H
HETATM   19  H62 HEX A   1       5.154   5.555   8.300  1.00 20.00           H
HETATM   20  H63 HEX A   1       5.486   5.724   9.843  1.00 20.00           H
END
'''

cif_str_06 = '''
data_default
_cell.length_a                    17.376
_cell.length_b                    14.809
_cell.length_c                    15.646
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      4026.047
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A

loop_
  _chem_comp.id
   A1EF7

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   HETATM 1 C02 . A1EF7 A 502 ? 5.86700 7.04300 5.78300 1.000 31.68000 C ? A ? . C02 1
   HETATM 2 C03 . A1EF7 A 502 ? 5.02400 7.82100 5.00000 1.000 29.68000 C ? A ? . C03 1
   HETATM 3 C04 . A1EF7 A 502 ? 5.00000 9.21300 5.19400 1.000 28.51000 C ? A ? . C04 1
   HETATM 4 C05 . A1EF7 A 502 ? 5.80000 9.80900 6.15600 1.000 26.18000 C ? A ? . C05 1
   HETATM 5 C06 . A1EF7 A 502 ? 6.62800 9.00600 6.92600 1.000 28.42000 C ? A ? . C06 1
   HETATM 6 C07 . A1EF7 A 502 ? 6.67800 7.61900 6.74600 1.000 32.03000 C ? A ? . C07 1
   HETATM 7 C08 . A1EF7 A 502 ? 7.59800 6.73000 7.62300 1.000 29.38000 C ? A ? . C08 1
   HETATM 8 C10 . A1EF7 A 502 ? 9.68800 6.64100 9.01800 1.000 29.91000 C ? A ? . C10 1
   HETATM 9 C11 . A1EF7 A 502 ? 9.35800 5.53500 9.79800 1.000 29.51000 C ? A ? . C11 1
   HETATM 10 C12 . A1EF7 A 502 ? 10.33000 5.00000 10.60900 1.000 28.64000 C ? A ? . C12 1
   HETATM 11 C13 . A1EF7 A 502 ? 11.57900 5.60100 10.64600 1.000 29.13000 C ? A ? . C13 1
   HETATM 12 C14 . A1EF7 A 502 ? 11.89800 6.71100 9.86900 1.000 28.71000 C ? A ? . C14 1
   HETATM 13 C15 . A1EF7 A 502 ? 10.92500 7.23100 9.05100 1.000 26.38000 C ? A ? . C15 1
   HETATM 14 C17 . A1EF7 A 502 ? 12.37600 8.61600 7.78700 1.000 25.22000 C ? A ? . C17 1
   HETATM 15 N09 . A1EF7 A 502 ? 8.76500 7.33200 8.15300 1.000 30.55000 N ? A ? . N09 1
   HETATM 16 O16 . A1EF7 A 502 ? 11.09000 8.34900 8.23400 1.000 25.01000 O ? A ? . O16 1
   HETATM 17 O18 . A1EF7 A 502 ? 7.42900 5.59000 7.91800 1.000 33.25000 O ? A ? . O18 1
   HETATM 18 SE01 . A1EF7 A 502 ? 6.00000 5.11100 5.62900 1.000 40.80000 Se ? A ? . SE01 1
'''

cif_restraints_06 = '''
# electronic Ligand Builder and Optimisation Workbench (eLBOW)
#   - a module of PHENIX version dev-svn-
#   - file written: Wed Feb 11 14:16:34 2026
#
#   Ligand name: "~{N}-(2-methoxyphenyl)-2-selanyl-benzamide"
#   Random seed: 3628800
#   SMILES string: COc1ccccc1NC(=O)c2ccccc2[SeH]
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
A1EF7      A1EF7 '"~{N}-(2-methoxyphenyl)-2-selanyl-benzamide"' non-polymer 31 18 .
#
data_comp_A1EF7
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
A1EF7       C10    C   CR6    0    .       1.1993    0.3700    1.0690
A1EF7       C17    C   CH3    0    .       3.9714   -1.4646    0.0561
A1EF7       C15    C   CR6    0    .       2.2637   -0.4871    1.3902
A1EF7       C14    C   CR16   0    .       2.9495   -0.3230    2.5784
A1EF7       C02    C   CR6    0    .      -3.2012   -0.6796   -1.2568
A1EF7       C03    C   CR16   0    .      -4.0825   -0.5466   -2.3282
A1EF7       C04    C   CR16   0    .      -3.7019    0.1796   -3.4632
A1EF7       C05    C   CR16   0    .      -2.4518    0.7664   -3.5246
A1EF7       C06    C   CR16   0    .      -1.5684    0.6338   -2.4521
A1EF7       C07    C   CR6    0    .      -1.9526   -0.0968   -1.3112
A1EF7       C08    C   C      0    .      -0.9634   -0.2536   -0.1016
A1EF7       C11    C   CR16   0    .       0.8297    1.3927    1.9505
A1EF7       C12    C   CR16   0    .       1.5129    1.5589    3.1383
A1EF7       C13    C   CR16   0    .       2.5823    0.6992    3.4630
A1EF7       N09    N   NH1    0    .       0.3122    0.4369   -0.1117
A1EF7       O16    O   O2     0    .       2.6360   -1.5322    0.4906
A1EF7       O18    O   O      0    .      -1.2950   -0.8728    0.8567
A1EF7      SE01    SE  SE     0    .      -3.7737   -1.5718    0.3437
A1EF7       H1     H   HCH3   0    .       4.0359   -1.8092   -0.9348
A1EF7       H2     H   HCH3   0    .       4.5913   -2.0845    0.6945
A1EF7       H3     H   HCH3   0    .       4.3028   -0.4735    0.1013
A1EF7       H4     H   HCR6   0    .       3.7812   -0.9926    2.8278
A1EF7       H5     H   HCR6   0    .      -5.0710   -1.0089   -2.2828
A1EF7       H6     H   HCR6   0    .      -4.3769    0.2811   -4.2799
A1EF7       H7     H   HCR6   0    .      -2.1475    1.3526   -4.4468
A1EF7       H8     H   HCR6   0    .      -0.6420    1.0676   -2.4933
A1EF7       H9     H   HCR6   0    .      -0.0039    2.0631    1.6952
A1EF7       H10    H   HCR6   0    .       1.2234    2.3605    3.8301
A1EF7       H11    H   HCR6   0    .       3.1129    0.8269    4.3842
A1EF7       H12    H   HNH1   0    .       0.6057    0.9443   -0.9350
A1EF7       H13    H   H      0    .      -4.6784   -0.7366    1.0523
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
A1EF7   C03     C04   aromatic      1.400 0.020
A1EF7   C03     C02   aromatic      1.394 0.020
A1EF7   C04     C05   aromatic      1.382 0.020
A1EF7  SE01     C02   single        1.920 0.020
A1EF7   C02     C07   aromatic      1.379 0.020
A1EF7   C05     C06   aromatic      1.396 0.020
A1EF7   C07     C06   aromatic      1.408 0.020
A1EF7   C07     C08   single        1.570 0.020
A1EF7   C08     O18   double        1.188 0.020
A1EF7   C08     N09   single        1.451 0.020
A1EF7   C17     O16   single        1.406 0.020
A1EF7   N09     C10   single        1.478 0.020
A1EF7   O16     C15   single        1.428 0.020
A1EF7   C10     C15   aromatic      1.404 0.020
A1EF7   C10     C11   aromatic      1.400 0.020
A1EF7   C15     C14   aromatic      1.382 0.020
A1EF7   C11     C12   aromatic      1.380 0.020
A1EF7   C14     C13   aromatic      1.401 0.020
A1EF7   C12     C13   aromatic      1.410 0.020
A1EF7   C17     H1    single        1.051 0.020
A1EF7   C17     H2    single        1.084 0.020
A1EF7   C17     H3    single        1.046 0.020
A1EF7   C14     H4    single        1.097 0.020
A1EF7   C03     H5    single        1.092 0.020
A1EF7   C04     H6    single        1.064 0.020
A1EF7   C05     H7    single        1.134 0.020
A1EF7   C06     H8    single        1.024 0.020
A1EF7   C11     H9    single        1.100 0.020
A1EF7   C12     H10   single        1.098 0.020
A1EF7   C13     H11   single        1.071 0.020
A1EF7   N09     H12   single        1.011 0.020
A1EF7  SE01     H13   single        1.421 0.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
A1EF7   N09     C10     C11         108.15 3.000
A1EF7   C11     C10     C15         120.14 3.000
A1EF7   N09     C10     C15         131.71 3.000
A1EF7   H3      C17     H2          109.58 3.000
A1EF7   H3      C17     H1          109.38 3.000
A1EF7   H2      C17     H1          109.42 3.000
A1EF7   H3      C17     O16         109.46 3.000
A1EF7   H2      C17     O16         109.48 3.000
A1EF7   H1      C17     O16         109.50 3.000
A1EF7   O16     C15     C14         119.94 3.000
A1EF7   C14     C15     C10         120.04 3.000
A1EF7   O16     C15     C10         120.02 3.000
A1EF7   H4      C14     C13         120.06 3.000
A1EF7   H4      C14     C15         119.97 3.000
A1EF7   C13     C14     C15         119.97 3.000
A1EF7  SE01     C02     C07         119.95 3.000
A1EF7   C07     C02     C03         120.13 3.000
A1EF7  SE01     C02     C03         119.78 3.000
A1EF7   H5      C03     C04         119.95 3.000
A1EF7   H5      C03     C02         120.00 3.000
A1EF7   C04     C03     C02         120.05 3.000
A1EF7   H6      C04     C05         119.94 3.000
A1EF7   H6      C04     C03         119.94 3.000
A1EF7   C05     C04     C03         120.13 3.000
A1EF7   H7      C05     C06         120.26 3.000
A1EF7   H7      C05     C04         119.87 3.000
A1EF7   C06     C05     C04         119.87 3.000
A1EF7   H8      C06     C07         119.96 3.000
A1EF7   H8      C06     C05         120.10 3.000
A1EF7   C07     C06     C05         119.94 3.000
A1EF7   C08     C07     C06         120.26 3.000
A1EF7   C08     C07     C02         119.85 3.000
A1EF7   C06     C07     C02         119.89 3.000
A1EF7   O18     C08     N09         119.95 3.000
A1EF7   N09     C08     C07         120.07 3.000
A1EF7   O18     C08     C07         119.83 3.000
A1EF7   H9      C11     C12         120.10 3.000
A1EF7   H9      C11     C10         119.95 3.000
A1EF7   C12     C11     C10         119.95 3.000
A1EF7   H10     C12     C13         120.01 3.000
A1EF7   H10     C12     C11         119.99 3.000
A1EF7   C13     C12     C11         120.01 3.000
A1EF7   H11     C13     C12         120.08 3.000
A1EF7   H11     C13     C14         120.03 3.000
A1EF7   C12     C13     C14         119.89 3.000
A1EF7   H12     N09     C08         120.00 3.000
A1EF7   H12     N09     C10         119.94 3.000
A1EF7   C08     N09     C10         120.04 3.000
A1EF7   C15     O16     C17         114.02 3.000
A1EF7   H13    SE01     C02         109.42 3.000
#
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
A1EF7 CONST_01       C13     C14     C15     C10           0.00   0.0 0
A1EF7 CONST_02       C13     C12     C11     C10           0.01   0.0 0
A1EF7 CONST_03       C12     C11     C10     C15          -0.01   0.0 0
A1EF7 CONST_04       C12     C13     C14     C15          -0.00   0.0 0
A1EF7 CONST_05       C11     C10     C15     C14           0.01   0.0 0
A1EF7 CONST_06       C11     C12     C13     C14          -0.00   0.0 0
A1EF7 CONST_07       C05     C04     C03     C02           0.00   0.0 0
A1EF7 CONST_08       C05     C06     C07     C02           0.01   0.0 0
A1EF7 CONST_09       C06     C05     C04     C03          -0.00   0.0 0
A1EF7 CONST_10       C06     C07     C02     C03          -0.01   0.0 0
A1EF7 CONST_11       C07     C02     C03     C04           0.00   0.0 0
A1EF7 CONST_12       C07     C06     C05     C04          -0.00   0.0 0
A1EF7 CONST_13       N09     C10     C15     C14         179.38   0.0 0
A1EF7 CONST_14       C08     C07     C02     C03        -179.86   0.0 0
A1EF7 CONST_15      SE01     C02     C03     C04         175.76   0.0 0
A1EF7 CONST_16       C08     C07     C06     C05         179.86   0.0 0
A1EF7 CONST_17      SE01     C02     C07     C06        -175.75   0.0 0
A1EF7 CONST_18       O16     C15     C10     C11         179.81   0.0 0
A1EF7 CONST_19       N09     C10     C11     C12        -179.52   0.0 0
A1EF7 CONST_20       O16     C15     C14     C13        -179.81   0.0 0
A1EF7 CONST_21       H4      C14     C15     C10        -179.96   0.0 0
A1EF7 CONST_22       H10     C12     C11     C10        -179.98   0.0 0
A1EF7 CONST_23       H9      C11     C10     C15         179.93   0.0 0
A1EF7 CONST_24       H11     C13     C14     C15         179.97   0.0 0
A1EF7 CONST_25       H6      C04     C03     C02        -179.92   0.0 0
A1EF7 CONST_26       H8      C06     C07     C02         179.95   0.0 0
A1EF7 CONST_27       H7      C05     C04     C03         179.92   0.0 0
A1EF7 CONST_28       H5      C03     C04     C05         179.99   0.0 0
A1EF7 CONST_29       C07     C08     N09     C10         175.47   0.0 0
A1EF7 CONST_30       C08     N09     C10     C15          93.68   0.0 0
A1EF7 CONST_31       C11     C10     N09     C08         -86.89   0.0 0
A1EF7 CONST_32       H12     N09     C10     C15         -84.73   0.0 0
A1EF7 CONST_33       O18     C08     N09     C10           0.03   0.0 0
A1EF7 Var_01         C14     C15     O16     C17         -58.41  30.0 2
A1EF7 Var_02         O18     C08     C07     C02          -0.08  30.0 2
A1EF7 Var_03         H1      C17     O16     C15        -145.99  30.0 3
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
A1EF7 plan-1    C10 0.020
A1EF7 plan-1    C15 0.020
A1EF7 plan-1    C14 0.020
A1EF7 plan-1    C11 0.020
A1EF7 plan-1    C12 0.020
A1EF7 plan-1    C13 0.020
A1EF7 plan-1    N09 0.020
A1EF7 plan-1    O16 0.020
A1EF7 plan-1     H4 0.020
A1EF7 plan-1     H9 0.020
A1EF7 plan-1    H10 0.020
A1EF7 plan-1    H11 0.020
A1EF7 plan-2    C02 0.020
A1EF7 plan-2    C03 0.020
A1EF7 plan-2    C04 0.020
A1EF7 plan-2    C05 0.020
A1EF7 plan-2    C06 0.020
A1EF7 plan-2    C07 0.020
A1EF7 plan-2    C08 0.020
A1EF7 plan-2     H5 0.020
A1EF7 plan-2     H6 0.020
A1EF7 plan-2     H7 0.020
A1EF7 plan-2     H8 0.020
A1EF7 plan-2   SE01 0.020
A1EF7 plan-3    C10 0.020
A1EF7 plan-3    C07 0.020
A1EF7 plan-3    C08 0.020
A1EF7 plan-3    N09 0.020
A1EF7 plan-3    O18 0.020
A1EF7 plan-3    H12 0.020
'''
# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
