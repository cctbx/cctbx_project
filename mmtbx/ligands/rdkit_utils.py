from __future__ import absolute_import, division, print_function

import copy

from libtbx.utils import Sorry

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdFMCS
from collections import defaultdict
from rdkit.Chem import Lipinski
from scitbx.array_family import flex

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL) # Only show critical errors

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
# ------------------------------------------------------------------------------
def get_atom_element(atom):
  pse = Chem.GetPeriodicTable()
  return pse.GetElementSymbol(atom.GetAtomicNum())

def is_hydrogen(molecule, i):
  atom=molecule.GetAtomWithIdx(i)
  return get_atom_element(atom) in ['H', 'D', 'T']

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

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

def approximate_residue_molecule(model, residue_group, altloc=""):
  """
  Fallback when residue_molecule fails: the residue conformer's atoms (blank plus
  altloc), bonds from the restraint file (by atom name), with its explicit
  single/double/triple orders (others single; to a metal: dative, from the other
  atom), no charges, ring perception only (not sanitized). Returns (mol,
  rdkit_to_iseq).
  """
  from scitbx.array_family import flex as _flex
  ags = [ag for ag in residue_group.atom_groups() if ag.altloc.strip() in ("", altloc)]
  resname = ([ag for ag in ags if ag.altloc.strip()] or ags)[0].resname.strip()
  conf = [a for ag in ags for a in ag.atoms()]
  comp, ani = model.get_mon_lib_srv().get_comp_comp_id_and_atom_name_interpretation(
    residue_name=resname, atom_names=_flex.std_string([a.name for a in conf]))
  names = [a.name.strip() for a in conf]
  if ani is not None:
    names = [(m or n).strip() for m, n in zip(ani.mon_lib_names(), names)]
  mol = Chem.RWMol()
  idx, rdkit_to_iseq, is_metal = {}, {}, set()
  for a, n in zip(conf, names):
    if a.element.strip().upper() in _metal_elements:
      is_metal.add(n)
    e = a.element.strip().capitalize() or a.name.strip()[:1]
    atom = Chem.Atom("H" if e == "D" else e)
    atom.SetNoImplicit(True)
    atom.SetProp("_Name", a.name.strip())
    k = mol.AddAtom(atom)
    idx[n] = k
    rdkit_to_iseq[k] = a.i_seq
  orders = {"sing": Chem.BondType.SINGLE, "doub": Chem.BondType.DOUBLE,
    "trip": Chem.BondType.TRIPLE}
  if comp is not None:
    for b in comp.bond_list:
      a1, a2 = b.atom_id_1.strip(), b.atom_id_2.strip()
      if a1 in idx and a2 in idx and mol.GetBondBetweenAtoms(idx[a1], idx[a2]) is None:
        if a1 in is_metal and a2 not in is_metal:
          a1, a2 = a2, a1
        if a2 in is_metal and a1 not in is_metal:
          mol.AddBond(idx[a1], idx[a2], Chem.BondType.DATIVE)
        else:
          mol.AddBond(idx[a1], idx[a2],
            orders.get((b.type or "").strip().lower()[:4], Chem.BondType.SINGLE))
  mol = mol.GetMol()
  mol.UpdatePropertyCache(strict=False)
  Chem.FastFindRings(mol)
  return mol, rdkit_to_iseq

def residue_rigid_components(model, residue_group, altloc="", filter_lone_linkers=True,
                             filename=None):
  """
  Rigid components of one residue conformer: from residue_molecule's fragment_mol
  (caps as implicit H; the residue's metals; added H have no i_seq and are left out
  of the components),
  or, if it fails, from approximate_residue_molecule. Returns group_args:
  components (flex.size_t of model i_seqs), mol, frags (rdkit indices, for
  drawing), approximate ("approximate: <reason>" or None), molecule (the
  residue_molecule result).
  """
  from libtbx import group_args
  try:
    r = residue_molecule(model, residue_group, altloc=altloc)
  except Exception as e:
    r = group_args(ok=False, reason="%s: %s" % (type(e).__name__, e))
  if r.ok:
    mol, rdkit_to_iseq, approximate = r.fragment_mol, r.fragment_to_iseq, None
  else:
    mol, rdkit_to_iseq = approximate_residue_molecule(model, residue_group, altloc)
    approximate = "approximate: %s" % r.reason
  components, mol, frags = get_rigid_components(mol, rdkit_to_iseq,
    filter_lone_linkers, filename)
  return group_args(components=components, mol=mol, frags=frags,
    approximate=approximate, molecule=r)

def get_cctbx_isel_for_rigid_components(model,
                                        residue_group,
                                        altloc="",
                                        filter_lone_linkers=True,
                                        filename=None):
  return residue_rigid_components(model, residue_group, altloc,
    filter_lone_linkers, filename).components

# ------------------------------------------------------------------------------

# Common non-metal atomic numbers (metalloids Si/As/Te included so they
# aren't treated as metals). Any atom outside this set is treated as a
# metal for fragmentation purposes: coordination-style bonds to metals
# aren't rotatable in the crystallographic sense, so we don't cut across
# them (e.g. Fe-O in an oFo ligand must stay a single rigid fragment).
_NON_METAL_Z = frozenset([
  1, 2,               # H, He
  5, 6, 7, 8, 9, 10,  # B, C, N, O, F, Ne
  14, 15, 16, 17, 18, # Si, P, S, Cl, Ar
  33, 34, 35, 36,     # As, Se, Br, Kr
  52, 53, 54,         # Te, I, Xe
  85, 86,             # At, Rn
])

def _bond_touches_metal(bond):
  return (bond.GetBeginAtom().GetAtomicNum() not in _NON_METAL_Z
       or bond.GetEndAtom().GetAtomicNum()   not in _NON_METAL_Z)

def get_rigid_components(mol,
                         rdkit_to_cctbx,
                         filter_lone_linkers=True,
                         filename=None):

  # Identify Rotatable Bonds
  rotatable_pattern = Lipinski.RotatableBondSmarts
  try:
    matches = mol.GetSubstructMatches(rotatable_pattern)
  except Exception as e:
    print('Failed to fragment the molecule.')
    frags = Chem.GetMolFrags(mol, asMols=False)
    return [flex.size_t(sorted(rdkit_to_cctbx.values()))], mol, list(frags)

  candidate_cut_bonds = []
  min_heavy_atoms = 2

  # Map: (bond_index, atom_index) -> Size of the fragment this atom ends up in
  # if bond is cut
  fragment_size_map = {}

  for u, v in matches:
    bond = mol.GetBondBetweenAtoms(u, v)
    if bond is None: continue

    if is_amide_bond(mol, bond): continue
    if _bond_touches_metal(bond): continue

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
    rdkit_frags = list(Chem.GetMolFrags(mol, asMols=False))
    draw_colored_fragments(mol, rdkit_frags, filename=filename)
    return [flex.size_t(sorted(rdkit_to_cctbx.values()))], mol, rdkit_frags

  fragmented_mol = Chem.FragmentOnBonds(mol, final_bonds_to_cut, addDummies=False)
  raw_fragments = list(Chem.GetMolFrags(fragmented_mol, asMols=False))

  # Convert to CCTBX format
  cctbx_rigid_components = []

  for frag in raw_fragments:
    component_indices = flex.size_t()
    for rd_idx in frag:
      # atoms without a model i_seq (H added from the restraint file) are left out
      if rd_idx in rdkit_to_cctbx:
        component_indices.append(rdkit_to_cctbx[rd_idx])
    if component_indices.size():
      cctbx_rigid_components.append(component_indices)

  draw_colored_fragments(mol, raw_fragments, filename=filename)

  return cctbx_rigid_components, mol, raw_fragments

# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------

def build_drawing_mol_with_missing(frag_mol, cif_object, missing_names):
  '''Return (draw_mol, missing_idxs): a copy of frag_mol with the named missing
  heavy atoms appended (element + atomLabel from cif_object) and their cif bonds
  added where both endpoints exist. Appending preserves existing atom indices.'''
  if not missing_names or cif_object is None or not hasattr(cif_object, 'atom_list'):
    return Chem.Mol(frag_mol), set()
  rw = Chem.RWMol(frag_mol)
  name_to_idx = {}
  for i in range(rw.GetNumAtoms()):
    a = rw.GetAtomWithIdx(i)
    if a.HasProp('_Name'):
      name_to_idx[a.GetProp('_Name').strip()] = i
  name_to_element = {}
  for ca in cif_object.atom_list:
    name_to_element[ca.atom_id.strip()] = ca.type_symbol.strip()
  missing_idxs = set()
  for name in missing_names:
    if name in name_to_idx:
      continue
    el = name_to_element.get(name)
    if not el:
      continue
    try:
      atom = Chem.Atom(el.capitalize())
    except (ValueError, RuntimeError):
      continue
    atom.SetProp('_Name', name)
    atom.SetProp('atomLabel', name)
    idx = rw.AddAtom(atom)
    name_to_idx[name] = idx
    missing_idxs.add(idx)
  if hasattr(cif_object, 'bond_list'):
    for bond in cif_object.bond_list:
      n1 = bond.atom_id_1.strip()
      n2 = bond.atom_id_2.strip()
      if n1 in name_to_idx and n2 in name_to_idx:
        if rw.GetBondBetweenAtoms(name_to_idx[n1], name_to_idx[n2]) is None:
          r_type = get_rdkit_bond_type(
            getattr(bond, 'type', 'sing'),
            elements={name_to_element.get(n1, 'C'), name_to_element.get(n2, 'C')})
          rw.AddBond(name_to_idx[n1], name_to_idx[n2], r_type)
  draw_mol = rw.GetMol()
  for atom in draw_mol.GetAtoms():
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(0)
    atom.UpdatePropertyCache(strict=False)
  return draw_mol, missing_idxs

# ------------------------------------------------------------------------------

_FIG_SUPERSAMPLE = 2
_FIG_BOND_PX = 60   # target bond length (logical px); sets the fixed figure scale
_FIG_MAX_PX = 800   # cap on either canvas dimension (logical px) before supersample

def draw_colored_fragments(mol, rdkit_frags, filename, use_atom_names=False,
                           frag_ccs=None, missing_atom_idxs=None, note=None):
  """
  1. Removes all H atoms.
  2. Strips all charges and implicit H properties (forces clean drawing).
  3. Maps colors from the original fragmented indices to the new clean molecule.
  frag_ccs: optional list of CC floats, one per fragment in rdkit_frags order.
            When provided, each CC value is drawn at the centroid of its fragment.
  note: optional text written below the figure (e.g. "approximate: <reason>").
  """
  if filename is None: return

  # 1. Work on a copy
  mol_viz = Chem.Mol(mol)

  # 2. Tag atoms with their original index so we can trace them back to rdkit_frags
  for atom in mol_viz.GetAtoms():
    atom.SetIntProp("orig_idx", atom.GetIdx())

  # 3. Remove Hydrogens
  # implicitOnly=False ensures we remove the explicit H atoms
  mol_viz = Chem.RemoveHs(mol_viz, implicitOnly=False, sanitize=False)

  # 4. Cleanup Loop
  # This prevents the "OH" or "S+" labels. Force RDKit to draw atoms as-is.
  for atom in mol_viz.GetAtoms():
    # Force RDKit to believe there are no implicit Hydrogens
    atom.SetNoImplicit(True)
    # Force explicit H count to 0
    atom.SetNumExplicitHs(0)
    # Optional: Remove formal charges (e.g. make S+ look like S)
    atom.SetFormalCharge(0)
    atom.SetNumRadicalElectrons(0)
    # Update property cache to accept these "weird" valences
    atom.UpdatePropertyCache(strict=False)

    # --- Labeling Logic ---
    if use_atom_names and atom.HasProp("_Name"):
      # RDKit looks for "atomLabel" to override the symbol
      name = atom.GetProp("_Name").strip()
      atom.SetProp("atomLabel", name)

  # 5. Coordinate Generation
  AllChem.Compute2DCoords(mol_viz)

  # 6. Color Mapping Logic
  # Neutral single-hue (brown/tan) shades used when no RSCC data is available.
  # Deliberately avoids green/amber/red (reserved for the RSCC traffic-light)
  # and grey (reserved for missing atoms), so fragment colouring can't be
  # mistaken for a density verdict. Ordered light/dark-alternating so adjacent
  # fragments contrast (between-fragment bonds are also left uncoloured).
  palette = [
    (0.878, 0.788, 0.651), (0.824, 0.706, 0.549), (0.788, 0.659, 0.463),
    (0.902, 0.835, 0.722), (0.749, 0.627, 0.439), (0.847, 0.761, 0.604),
  ]

  # Per-fragment color: traffic-light by CC when provided, else neutral shades.
  frag_colors = {}
  for i in range(len(rdkit_frags)):
    if frag_ccs is not None and i < len(frag_ccs):
      cc = frag_ccs[i]
      if   cc >= 0.8: frag_colors[i] = (0.6, 1.0, 0.6)
      elif cc >= 0.6: frag_colors[i] = (1.0, 1.0, 0.6)
      else:           frag_colors[i] = (1.0, 0.6, 0.6)
    else:
      frag_colors[i] = palette[i % len(palette)]

  # Map: Original_Index -> Fragment_ID
  old_idx_to_frag_id = {}
  for frag_idx, frag_atoms in enumerate(rdkit_frags):
    for old_idx in frag_atoms:
      old_idx_to_frag_id[old_idx] = frag_idx

  # Map: New_Atom_Index -> Color
  new_atom_highlights = {}
  new_bond_highlights = {}
  new_idx_to_frag_id = {}

  for atom in mol_viz.GetAtoms():
    if atom.HasProp("orig_idx"):
      old_idx = atom.GetIntProp("orig_idx")

      if old_idx in old_idx_to_frag_id:
        frag_id = old_idx_to_frag_id[old_idx]
        color = frag_colors[frag_id]

        new_idx = atom.GetIdx()
        new_atom_highlights[new_idx] = color
        new_idx_to_frag_id[new_idx] = frag_id

  # Grey out missing atoms (in no fragment).
  new_missing_idxs = set()
  if missing_atom_idxs:
    missing_set = set(missing_atom_idxs)
    for atom in mol_viz.GetAtoms():
      if atom.HasProp("orig_idx") and atom.GetIntProp("orig_idx") in missing_set:
        new_idx = atom.GetIdx()
        new_atom_highlights[new_idx] = (0.7, 0.7, 0.7)
        new_missing_idxs.add(new_idx)

  # 7. Highlight Bonds
  for bond in mol_viz.GetBonds():
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()

    if a1 in new_idx_to_frag_id and a2 in new_idx_to_frag_id:
      if new_idx_to_frag_id[a1] == new_idx_to_frag_id[a2]:
        new_bond_highlights[bond.GetIdx()] = new_atom_highlights[a1]
    elif a1 in new_missing_idxs or a2 in new_missing_idxs:
      new_bond_highlights[bond.GetIdx()] = (0.7, 0.7, 0.7)

  # 8. Draw
  # Draw every ligand at a FIXED bond length so atoms and highlight blobs render
  # at the same physical scale regardless of size; the canvas is then fitted to
  # the molecule's 2D bounding box at that scale (+ padding). This makes the
  # figure size a true visual cue of molecule size: a small ligand (e.g. EDO)
  # yields a small figure, a large one a large figure. rdkit is NOT allowed to
  # fit-to-canvas (which previously blew tiny molecules up to fill the frame).
  bond_px = _FIG_BOND_PX * _FIG_SUPERSAMPLE
  conf = mol_viz.GetConformer()
  xs = [conf.GetAtomPosition(i).x for i in range(mol_viz.GetNumAtoms())]
  ys = [conf.GetAtomPosition(i).y for i in range(mol_viz.GetNumAtoms())]
  bb_w = max(xs) - min(xs) if xs else 1.0
  bb_h = max(ys) - min(ys) if ys else 1.0
  # Pixels per molecular unit from the median bond length in the 2D conformer.
  bond_lens = []
  for bond in mol_viz.GetBonds():
    p1 = conf.GetAtomPosition(bond.GetBeginAtomIdx())
    p2 = conf.GetAtomPosition(bond.GetEndAtomIdx())
    bond_lens.append(((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2) ** 0.5)
  bond_lens.sort()
  median_bond = bond_lens[len(bond_lens) // 2] if bond_lens else 1.5
  ppu = bond_px / (median_bond if median_bond > 1e-6 else 1.5)
  pad = int(0.9 * bond_px)          # clears highlight blobs and atom labels
  min_w = int(3.4 * bond_px)        # floor: keep the legend row from clipping
  min_h = int(2.4 * bond_px)        # floor: linear ligands get vertical room
  max_px = _FIG_MAX_PX * _FIG_SUPERSAMPLE
  width  = max(min_w, min(max_px, int(bb_w * ppu) + 2 * pad))
  height = max(min_h, min(max_px, int(bb_h * ppu) + 2 * pad))
  drawer = rdMolDraw2D.MolDraw2DCairo(width, height)

  opts = drawer.drawOptions()
  opts.fillHighlights = True
  opts.fixedBondLength = bond_px  # pin the scale; do not fit-to-canvas
  opts.padding = 0.05  # tight fit; badges + legend are added via PIL below

  drawer.DrawMolecule(
    mol_viz,
    highlightAtoms=list(new_atom_highlights.keys()),
    highlightAtomColors=new_atom_highlights,
    highlightBonds=list(new_bond_highlights.keys()),
    highlightBondColors=new_bond_highlights
  )

  # Collect pixel-space centroids of each fragment so we can place numbered
  # badges. GetDrawCoords must be called after DrawMolecule but before
  # FinishDrawing.
  frag_pixel_centroids = {}
  if frag_ccs is not None:
    from rdkit.Geometry import rdGeometry
    conf = mol_viz.GetConformer()
    frag_mol_pts = {}
    for new_idx, frag_id in new_idx_to_frag_id.items():
      pos = conf.GetAtomPosition(new_idx)
      frag_mol_pts.setdefault(frag_id, []).append((pos.x, pos.y))
    for frag_id, pts in frag_mol_pts.items():
      cx = sum(x for x, y in pts) / len(pts)
      cy = sum(y for x, y in pts) / len(pts)
      canvas_pt = drawer.GetDrawCoords(rdGeometry.Point2D(cx, cy))
      frag_pixel_centroids[frag_id] = (canvas_pt.x, canvas_pt.y)

  drawer.FinishDrawing()
  png_bytes = drawer.GetDrawingText()

  # 9. Composite numbered badges at each fragment centroid and a wrapped
  #    legend below the molecule, so every CC stays linked to its fragment
  #    regardless of molecule shape or fragment count.
  has_missing = bool(missing_atom_idxs)
  if (frag_ccs is not None and frag_pixel_centroids) or has_missing:
    try:
      from PIL import Image, ImageDraw as PILDraw, ImageFont
      import io as _io

      legend_font_size = 20 * _FIG_SUPERSAMPLE
      badge_font_size = 16 * _FIG_SUPERSAMPLE
      badge_r = 13 * _FIG_SUPERSAMPLE
      badge_text_gap = 6 * _FIG_SUPERSAMPLE
      try:
        legend_font = ImageFont.load_default(size=legend_font_size)
      except TypeError:
        legend_font = ImageFont.load_default()
      try:
        badge_font = ImageFont.load_default(size=badge_font_size)
      except TypeError:
        badge_font = ImageFont.load_default()

      def cc_text_rgb(cc):
        if cc >= 0.8: return (40, 140, 40)
        if cc >= 0.6: return (180, 110, 0)
        return (180, 40, 40)

      tmp_draw = PILDraw.Draw(Image.new('RGBA', (1, 1)))
      def text_size(s, font):
        try:
          bb = tmp_draw.textbbox((0, 0), s, font=font)
          return bb[2] - bb[0], bb[3] - bb[1]
        except AttributeError:
          return tmp_draw.textsize(s, font=font)

      def draw_badge(d, cx, cy, n):
        d.ellipse(
          (cx - badge_r, cy - badge_r, cx + badge_r, cy + badge_r),
          fill=(255, 255, 255, 235), outline=(40, 40, 40), width=2 * _FIG_SUPERSAMPLE)
        try:
          d.text((cx, cy), str(n), fill=(0, 0, 0),
                 font=badge_font, anchor='mm')
        except TypeError:
          tw, th = text_size(str(n), badge_font)
          d.text((cx - tw // 2, cy - th // 2), str(n),
                 fill=(0, 0, 0), font=badge_font)

      # Number fragments 1..N left-to-right by centroid x.
      sorted_frag_ids = sorted(
        frag_pixel_centroids,
        key=lambda fid: frag_pixel_centroids[fid][0])
      frag_id_to_n = {fid: i + 1 for i, fid in enumerate(sorted_frag_ids)}

      # Legend entries in the same left-to-right order.
      entries = []
      for fid in sorted_frag_ids:
        if fid >= len(frag_ccs):
          continue
        cc = frag_ccs[fid]
        cc_str = '%.2f' % cc
        cc_w, _ = text_size(cc_str, legend_font)
        entry_w = 2 * badge_r + badge_text_gap + cc_w
        entries.append((frag_id_to_n[fid], cc_str, cc_text_rgb(cc), entry_w))

      prefix = 'RSCC  '
      sep = '   '
      prefix_w, _ = text_size(prefix, legend_font)
      sep_w, _ = text_size(sep, legend_font)
      row_h = max(legend_font_size, 2 * badge_r) + 6 * _FIG_SUPERSAMPLE

      # Wrap entries into rows.
      max_row_w = width - 20 * _FIG_SUPERSAMPLE
      rows = []
      cur_row = []
      cur_w = prefix_w
      for i, (_, _, _, entry_w) in enumerate(entries):
        add_w = entry_w + (sep_w if cur_row else 0)
        if cur_row and cur_w + add_w > max_row_w:
          rows.append(cur_row)
          cur_row = [i]
          cur_w = entry_w
        else:
          cur_row.append(i)
          cur_w += add_w
      if cur_row: rows.append(cur_row)

      legend_pad_top = 8 * _FIG_SUPERSAMPLE
      legend_pad_bot = 6 * _FIG_SUPERSAMPLE
      n_legend_rows = len(rows) + (1 if has_missing else 0)
      legend_h = legend_pad_top + row_h * n_legend_rows + legend_pad_bot
      canvas_h = height + legend_h

      img = Image.open(_io.BytesIO(png_bytes)).convert('RGBA')
      new_img = Image.new('RGBA', (width, canvas_h), (255, 255, 255, 255))
      new_img.paste(img, (0, 0))
      draw = PILDraw.Draw(new_img)

      # Badges on the molecule.
      for fid, (cx, cy) in frag_pixel_centroids.items():
        if fid >= len(frag_ccs):
          continue
        draw_badge(draw, int(cx), int(cy), frag_id_to_n[fid])

      # Legend: same badge graphic + tier-coloured CC value.
      y_row = height + legend_pad_top
      for row_i, row in enumerate(rows):
        x = 10 * _FIG_SUPERSAMPLE
        row_mid = y_row + row_h // 2
        if row_i == 0:
          try:
            draw.text((x, row_mid), prefix, fill=(80, 80, 80),
                      font=legend_font, anchor='lm')
          except TypeError:
            draw.text((x, y_row), prefix, fill=(80, 80, 80), font=legend_font)
          x += prefix_w
        for j, idx in enumerate(row):
          n, cc_str, color, entry_w = entries[idx]
          if j > 0:
            x += sep_w
          draw_badge(draw, x + badge_r, row_mid, n)
          tx = x + 2 * badge_r + badge_text_gap
          try:
            draw.text((tx, row_mid), cc_str, fill=color,
                      font=legend_font, anchor='lm')
          except TypeError:
            draw.text((tx, y_row), cc_str, fill=color, font=legend_font)
          x += entry_w
        y_row += row_h

      if has_missing:
        x = 10 * _FIG_SUPERSAMPLE
        row_mid = y_row + row_h // 2
        draw.ellipse(
          (x, row_mid - badge_r, x + 2 * badge_r, row_mid + badge_r),
          fill=(179, 179, 179), outline=(40, 40, 40),
          width=2 * _FIG_SUPERSAMPLE)
        tx = x + 2 * badge_r + badge_text_gap
        try:
          draw.text((tx, row_mid), '= missing atoms', fill=(80, 80, 80),
                    font=legend_font, anchor='lm')
        except TypeError:
          draw.text((tx, y_row), '= missing atoms', fill=(80, 80, 80),
                    font=legend_font)
        y_row += row_h

      out = _io.BytesIO()
      new_img.save(out, format='PNG')
      png_bytes = out.getvalue()
    except ImportError:
      pass  # PIL not available; save without annotations

  # 10. A note below the figure
  if note:
    try:
      from PIL import Image, ImageDraw as PILDraw, ImageFont
      import io as _io
      import textwrap
      font_size = 16 * _FIG_SUPERSAMPLE
      try:
        font = ImageFont.load_default(size=font_size)
      except TypeError:
        font = ImageFont.load_default()
      img = Image.open(_io.BytesIO(png_bytes)).convert('RGBA')
      w, h = img.size
      lines = textwrap.wrap(note, max(20, int(w / (0.55 * font_size))))
      row = font_size + 6 * _FIG_SUPERSAMPLE
      new_img = Image.new('RGBA', (w, h + row * len(lines) + 8 * _FIG_SUPERSAMPLE),
        (255, 255, 255, 255))
      new_img.paste(img, (0, 0))
      draw = PILDraw.Draw(new_img)
      y = h + 4 * _FIG_SUPERSAMPLE
      for line in lines:
        draw.text((10 * _FIG_SUPERSAMPLE, y), line, fill=(160, 40, 40), font=font)
        y += row
      out = _io.BytesIO()
      new_img.save(out, format='PNG')
      png_bytes = out.getvalue()
    except ImportError:
      pass

  # 11. Save
  with open(filename, 'wb') as f:
    f.write(png_bytes)

# ------------------------------------------------------------------------------

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
  return rc

def overlay_two_molecules(mol1, mol2):
  from rdkit.Chem import rdFMCS

  params = rdFMCS.MCSParameters()
  params.AtomTyper = rdFMCS.AtomCompare.CompareElements
  params.BondTyper = rdFMCS.BondCompare.CompareOrder
  params.BondCompareParameters.RingMatchesRingOnly = True
  params.BondCompareParameters.CompleteRingsOnly = True

  res = rdFMCS.FindMCS([mol1, mol2], params)

  highlightAtoms_mol1 = mol1.GetSubstructMatch(res.queryMol)
  print(highlightAtoms_mol1)
  highlightAtoms_mol2 = mol2.GetSubstructMatch(res.queryMol)
  print(highlightAtoms_mol2)

  return list(zip(highlightAtoms_mol1, highlightAtoms_mol2))


def read_chemical_component_smiles(filename):
  from iotbx import cif
  ccd = cif.reader(filename).model()
  for i, (code, monomer) in enumerate(ccd.items()):
    # molecule = Chem.Mol()
    desc = monomer.get_loop_or_row('_pdbx_chem_comp_descriptor')
    for j, tmp in enumerate(desc.iterrows()):
      smiles_type = tmp.get('_pdbx_chem_comp_descriptor.type')
      if smiles_type=='SMILES_CANONICAL':
        smiles = tmp.get('_pdbx_chem_comp_descriptor.descriptor')
        return smiles

def read_chemical_component_filename_new(filename):
  smiles=read_chemical_component_smiles(filename)
  print('smiles',smiles)
  mol1=mol_from_smiles(smiles)
  print(mol1)
  mol2=read_chemical_component_filename_old(filename)
  print(mol2)
  print(Chem.FindMolChiralCenters(mol1, includeUnassigned=True))
  print(Chem.FindMolChiralCenters(mol2, includeUnassigned=True))

  matches = overlay_two_molecules(mol1, mol2)
  print(matches)

  for k, (i, j) in enumerate(matches):
    # print(dir(mol1.GetAtomWithIdx(i)))
    print(k,i,j,mol1.GetAtomWithIdx(i).GetSymbol(), mol2.GetAtomWithIdx(j).GetSymbol())

  assert 0

def read_chemical_component_filename(filename, verbose=False):
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
  Chem.SanitizeMol(rwmol)
  # from rdkit.Chem.PropertyMol import PropertyMol
  # molecule = PropertyMol(molecule)
  molecule = rwmol.GetMol()
  for key, item in desc.items():
    key = key.split('.')[1]
    molecule.SetProp(key,item[0])
  molecule.SetProp('id', code.upper())
  if verbose:
    for atom in molecule.GetAtoms():
      print(atom, atom.GetPropsAsDict())
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

def molecule_from_chemical_component(code):
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

  for i, atom in enumerate(cctbx_model.get_atoms()):
    element = atom.element.strip().upper()
    if element =="D":
      element = "H"
    else:
      element = element
    atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(element.capitalize())
    rdatom = Chem.Atom(atomic_number)
    rdatom.SetFormalCharge(atom.charge_as_int())
    rdatom_idx = rwmol.AddAtom(rdatom)
    conformer.SetAtomPosition(rdatom_idx,atom.xyz)

  atoms=cctbx_model.get_atoms()
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

def _generate_models_from_residues(cctbx_model):
  for rg in cctbx_model.get_hierarchy().residue_groups():
    assert len(rg.atom_groups())==1
    chain=rg.parent().id
    resseq=rg.resseq.strip()
    s=f'chain {chain} and resseq {resseq}'
    sel = cctbx_model.selection(s)
    tm = cctbx_model.select(sel)
    yield tm

def _generate_models_from_fragments(cctbx_model):
  assert 0
  from mmtbx.conformation_dependent_library import generate_protein_fragments
  pdb_hierarchy=cctbx_model.get_hierarchy()
  geometry_restraints_manager=cctbx_model.get_restraints_manager().geometry
  selections=[[]]
  for j, threes in enumerate(generate_protein_fragments(pdb_hierarchy,
                                                        geometry_restraints_manager,
                                                        length=2,
                                                        include_non_linked=True,
                                                        include_non_standard_peptides=True,
                                                        include_d_amino_acids=True,
                                                        # verbose=1,
                                                        )):
    if threes.are_linked():
      # selections[-1].append(threes[0].resseq.strip())
      chain_resseq=(threes[1].parent().parent().id, threes[1].resseq.strip())
      if chain_resseq not in selections[-1]: selections[-1].append(chain_resseq)
    else:
      chain_resseq=(threes[0].parent().parent().id, threes[0].resseq.strip())
      if chain_resseq not in selections[-1]: selections[-1].append(chain_resseq)
      selections.append([])
      chain_resseq=(threes[1].parent().parent().id, threes[1].resseq.strip())
      if chain_resseq not in selections[-1]: selections[-1].append(chain_resseq)
    # yield threes.are_linked()
  for st in selections:
    if len(st)==1:
      chain, resseq = st[0]
      s=f'chain {chain} and resseq {resseq}'
    else:
      chain, resseq1 = st[0]
      tmp, resseq2 = st[-1]
      s=f'chain {chain} and resseq {resseq1}:{resseq2}'
    sel = cctbx_model.selection(s)
    tm = cctbx_model.select(sel)
    yield tm

def convert_model_to_rdkit_molecules(cctbx_model):
  mols=[]
  mods=[]
  for tm in _generate_models_from_residues(cctbx_model):
    mods.append(tm)
    mol = convert_model_to_rdkit(tm)
    mols.append(mol)
  return mols, mods

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
  Convert a rdkit mol to 3D coordinates
  """
  assert len(mol.GetConformers())==0, "mol already has conformer"
  param = rdDistGeom.ETKDGv3()
  conf_id = rdDistGeom.EmbedMolecule(mol, clearConfs=True)
  return mol

def mol_to_2d(mol):
  """
  Convert a rdkit mol to 2D coordinates
  """
  mol = Chem.Mol(mol) # copy to preserve original coords
  ret = Chem.rdDepictor.Compute2DCoords(mol)
  return mol

def molecule_to_image(molecule, file_name, size=(2000,1500), draw2d=True):
  from rdkit.Chem import AllChem
  from rdkit.Chem import Draw
  molecule=copy.deepcopy(molecule)
  if draw2d: AllChem.Compute2DCoords(molecule)
  image = Draw.MolToImage(molecule, size=size)
  return image

def print_coordinates(mol):
  conformer = mol.GetConformer()
  positions=conformer.GetPositions()
  print("Atomic coordinates (x, y, z): %s" % conformer.Is3D())
  for i, atom in enumerate(mol.GetAtoms()):
    position = conformer.GetAtomPosition(i)
    new_xyz=positions[i]
    print(f"{atom.GetSymbol()} {position.x:.4f} {position.y:.4f} {position.z:.4f}")

def populate_molecule(rdmol, embed3d=True, addHs=True, removeHs=False, verbose=False):
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
  if embed3d:
    conf_id = AllChem.EmbedMolecule(rdmol, AllChem.ETKDG())
    if conf_id == -1:
      raise Sorry("Could not generate 3D coordinates")
      # return None
    AllChem.MMFFOptimizeMolecule(rdmol, confId=conf_id)
    if verbose:
      print_coordinates(rdmol)
  return rdmol

def mol_from_smiles(smiles, embed3d=True, addHs=True, removeHs=False, verbose=False):
  """
  Convert a smiles string to rdkit mol
  """
  ps = Chem.SmilesParserParams()
  ps.removeHs=removeHs
  rdmol = Chem.MolFromSmiles(smiles, ps)
  if verbose: print('rdmol',rdmol)
  if rdmol is None:
    raise Sorry(f'invalid SMILES {smiles}')
    # return rdmol
  return populate_molecule(rdmol,
                           embed3d=embed3d,
                           addHs=addHs,
                           removeHs=removeHs,
                           verbose=verbose)

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

# ------------------------------------------------------------------------------
# RDKit molecule of one residue conformer of a processed model

_metal_elements = frozenset("""LI BE NA MG AL K CA SC TI V CR MN FE CO NI CU ZN GA RB
  SR Y ZR NB MO TC RU RH PD AG CD IN SN CS BA LA CE PR ND PM SM EU GD TB DY HO ER
  TM YB LU HF TA W RE OS IR PT AU HG TL PB BI PO FR RA AC TH PA U NP PU""".split())

# X-H lengths (A) for the H caps and added H; positions do not enter the bond orders
_cap_length = {"C": 1.09, "N": 1.01, "O": 0.96, "S": 1.34, "SE": 1.47, "P": 1.42,
  "B": 1.19, "SI": 1.48}

def _unit(v):
  import math
  n = math.sqrt(sum([c * c for c in v]))
  if n < 1.e-6:
    return None
  return tuple([c / n for c in v])

def _h_site(xyz_x, directions, k, element):
  """An H on the atom at xyz_x along the mean direction, turned a little for the k-th H."""
  d = _unit([sum([v[c] for v in directions]) for c in range(3)]) if directions else None
  if d is None:
    d = (1.0, 0.0, 0.0)
  # spread several H on one atom: tilt by k around a perpendicular axis
  p = _unit((d[1], -d[0], 0.0)) or (0.0, 1.0, 0.0)
  t = 0.6 * k
  d = _unit([d[c] + t * p[c] for c in range(3)])
  r = _cap_length.get(element, 1.0)
  return tuple([xyz_x[c] + r * d[c] for c in range(3)])

def residue_molecule(model, residue_group, altloc="", fsc0=None):
  """
  RDKit molecule of one residue conformer (blank-altloc atoms plus those of altloc)
  of a processed model, with bond orders and formal charges from
  rdkit.Chem.rdDetermineBonds.DetermineBondOrders.

  Atoms: the conformer's atoms (_Name = atom name), metals left out. An H whose
  heavy atom is absent is left out (noted). Bonds: the restraints' bonds within the
  conformer. H caps (along the bond, or away from the present neighbours): one per
  bond to an atom outside the residue (polymer bond, covalent link; links counted
  once per partner atom (chain, residue, name), whatever its altloc), one per
  restraint-file heavy atom missing from the model (on its present neighbour,
  unless that neighbour is linked: the link replaces the leaving atom), one on an
  unlinked, uncharged N or C of a peptide restraint file whose file bonds leave one
  valence open (the absent neighbour residue): explicit orders summing to 2 on N or
  3 on C (N with only CA and H; not an imine N=CA), else two file neighbours (e.g.
  a monomer-library in-chain file without OXT).

  Hydrogens: the restraint file is the reference protonation. Per heavy atom, the
  file's H minus those replaced by link caps (one per cap, less the atom's missing
  restraint-file heavy neighbours; none for a polymer bond of a polymer restraint
  file, e.g. the peptide N-C), compared with the model's H by count:
    - fewer on a metal-bound atom: deprotonation (-1 each);
    - fewer elsewhere: the file's H added (no i_seq), noted in h_differences;
    - more: failure (an H the restraint file does not have).
  Without any H in the conformer: the file's H on every atom (no deprotonation).

  Total charge: the file's formal charges (_chem_comp_atom.charge) over the atoms
  present, minus the deprotonations.
  File first: with formal charges, and an explicit single/double/triple/aromatic
  order in the file for every bond between the heavy atoms present, the molecule
  takes the file's orders and charges (a deprotonated atom's charge lowered by its
  deprotonations); used (source "restraint file") if it sanitizes without radicals
  at that total. Otherwise DetermineBondOrders at that total, in the input atom
  order, then RDKit's canonical order, then 10 seeded random orders (it depends on
  the order); the first that succeeds and agrees with the file is taken
  (search["order"]); none: failure.
  Without formal charges: the totals -4..+4; valid: DetermineBondOrders succeeds and
  the bonds agree with the file; totals whose structure charges a carbon set aside
  when others are valid; one left: taken, charge_certain False; else failure.
  Agreement: the file's explicit single/double/triple bonds between heavy atoms,
  bonds aromatic in RDKit apart; bonds from one atom to neighbours of one element
  (carboxylate, nitro, guanidinium, phosphate) compared as a set of orders.

  Returns group_args: ok, reason (why not ok), mol (None unless ok; never
  unsanitized; metal-free), rdkit_to_iseq {mol index: model i_seq} (caps and added
  H have none), fragment_mol (mol with the caps as implicit H and the residue's
  metals, dative bonds from the residue's atoms), fragment_to_iseq {fragment_mol
  index: model i_seq}, caps [dict(index, on, kind: linked/missing, partner)],
  linked, missing_neighbour, metal_bound (model i_seqs), added_h (mol indices),
  hydrogens ("model", "completed from the restraint file", "restraint file (no H
  in the model)"), h_differences [dict(atom, model, restraint_file, kind:
  added/deprotonation/extra, added)], uncertain_atoms (names of heteroatoms whose
  H come from the restraint file; uncertain_iseqs: their i_seqs), total_charge, total_charge_source (restraint
  file / formal charges / search), charge_certain, search (calls, seconds, valid,
  set_aside, disagree, order), charge_notes, differences (information: bonds against the file's
  explicit bonds one by one, aromatic apart; formal charges against the file's,
  resonance-equivalent atoms apart), seconds, resname.
  """
  import time
  from libtbx import group_args
  from rdkit.Chem import rdDetermineBonds
  from mmtbx.monomer_library import cif_types
  t0 = time.time()
  result = group_args(ok=False, reason=None, mol=None, fragment_mol=None,
    rdkit_to_iseq={}, fragment_to_iseq={}, caps=[], linked=[], missing_neighbour=[],
    metal_bound=[], added_h=[], hydrogens=None, h_differences=[], uncertain_atoms=[],
    uncertain_iseqs=[],
    total_charge=None, total_charge_source=None, charge_certain=None, charge_notes=[],
    differences=dict(bonds=[], charges=[]), seconds=None, resname=None, search=None)
  def fail(reason):
    result.reason = reason
    result.mol = None
    result.fragment_mol = None
    result.seconds = time.time() - t0
    return result
  atoms = model.get_hierarchy().atoms()
  if fsc0 is None:
    fsc0 = model.get_restraints_manager().geometry.shell_sym_tables[0] \
      .full_simple_connectivity()
  def el(i):
    return atoms[i].element.strip().upper()
  def is_h(i):
    return el(i) in ("H", "D")
  def visible(i):
    return atoms[i].parent().altloc.strip() in ("", altloc)
  ags = [ag for ag in residue_group.atom_groups() if ag.altloc.strip() in ("", altloc)]
  if not ags:
    return fail("no atoms for altloc %r" % altloc)
  resname = ([ag for ag in ags if ag.altloc.strip()] or ags)[0].resname.strip()
  result.resname = resname
  conf = [a for ag in ags for a in ag.atoms()]
  conf_seqs = set([a.i_seq for a in conf])
  srv = model.get_mon_lib_srv()
  from scitbx.array_family import flex as _flex
  comp, ani = srv.get_comp_comp_id_and_atom_name_interpretation(residue_name=resname,
    atom_names=_flex.std_string([a.name for a in conf]))
  if comp is None:
    return fail("no restraint file for %s" % resname)
  names = [a.name.strip() for a in conf]
  if ani is not None:
    mon = ani.mon_lib_names()
    names = [(m or n).strip() for m, n in zip(mon, names)]
  dict_name = dict([(a.i_seq, n) for a, n in zip(conf, names)])
  d_atoms = dict([(a.atom_id.strip(), a) for a in comp.atom_list])
  d_el = dict([(n, (a.type_symbol or "").strip().upper()) for n, a in d_atoms.items()])
  d_nb = dict([(n, []) for n in d_atoms])
  for b in comp.bond_list:
    a1, a2 = b.atom_id_1.strip(), b.atom_id_2.strip()
    if a1 in d_nb and a2 in d_nb:
      d_nb[a1].append(a2)
      d_nb[a2].append(a1)
  def d_h(n):
    return [x for x in d_nb.get(n, []) if d_el.get(x) in ("H", "D")]
  # atoms present (metals left out), bonds, links, metal contacts
  metals = [a.i_seq for a in conf if el(a.i_seq) in _metal_elements]
  present = [a.i_seq for a in conf if el(a.i_seq) not in _metal_elements]
  # H whose heavy atom is absent (e.g. HO3 of a GOL without O3): left out
  orphans = [i for i in present if is_h(i) and not [k for k in fsc0[i]
    if visible(k) and k in conf_seqs and not is_h(k)]]
  if orphans:
    result.charge_notes.append("H without its heavy atom, left out: %s" % " ".join(
      [atoms[i].name.strip() for i in orphans]))
  present = [i for i in present if i not in orphans]
  present_set = set(present)
  heavy = [i for i in present if not is_h(i)]
  has_h = len(heavy) < len(present)
  # links: one per partner atom (chain, residue, name), whatever its altloc (a
  # neighbour split into conformers is still one partner)
  rg_key = (residue_group.parent().id, residue_group.resseq, residue_group.icode)
  def residue_key(k):
    rg = atoms[k].parent().parent()
    return (rg.parent().id, rg.resseq, rg.icode)
  bonds, links, metal_bound = set(), [], set()
  for i in present:
    partners = set()
    for k in fsc0[i]:
      if residue_key(k) == rg_key:
        if visible(k) and k in present_set:
          bonds.add(tuple(sorted([i, k])))
        elif visible(k) and el(k) in _metal_elements:
          metal_bound.add(i)
        continue
      if el(k) in _metal_elements:
        metal_bound.add(i)
        continue
      key = residue_key(k) + (atoms[k].name.strip(),)
      if key not in partners:
        partners.add(key)
        links.append((i, k))
  # restraint-file heavy atoms missing from the model: cap their present neighbour
  present_names = dict([(dict_name[i], i) for i in present])
  linked = set([i for i, k in links])
  missing_caps, n_missing_nb = [], {}
  for n, a in sorted(d_atoms.items()):
    if n in present_names or d_el.get(n) in ("H", "D") or d_el.get(n) in _metal_elements:
      continue
    for nb in d_nb[n]:
      i = present_names.get(nb)
      if i is not None and not is_h(i):
        n_missing_nb[i] = n_missing_nb.get(i, 0) + 1
        if i not in linked:
          missing_caps.append((i, n))
  # caps that replace a restraint-file H: links other than polymer bonds of a
  # polymer restraint file
  group = (comp.chem_comp.group or "").strip().upper()
  polymer = group in ("L-PEPTIDE", "D-PEPTIDE", "PEPTIDE", "M-PEPTIDE", "P-PEPTIDE",
    "DNA", "RNA")
  polymer_pairs = (("N", "C"), ("C", "N"), ("P", "O3'"), ("O3'", "P"), ("P", "O3*"),
    ("O3*", "P"))
  # a peptide chain end where the file leaves the link's valence open (unlinked N or
  # C, uncharged, explicit orders summing to 2 on N or 3 on C, e.g. N-CA, N-H of a
  # peptide-linking file; not an imine N=CA): the absent neighbour residue capped
  if polymer and group != "DNA" and group != "RNA":
    kinds = {"sing": 1, "doub": 2, "trip": 3}
    file_orders = {}
    for b in comp.bond_list:
      t = kinds.get((b.type or "").strip().lower()[:4])
      for x in (b.atom_id_1.strip(), b.atom_id_2.strip()):
        file_orders.setdefault(x, []).append(t)
    for n, valence, where in (("N", 3, "preceding"), ("C", 4, "following")):
      i = present_names.get(n)
      o = file_orders.get(n, [])
      q = cif_types.formal_charge_and_problem(getattr(d_atoms.get(n), "charge", None))[0]
      # orders not explicit (e.g. monomer library 'coval'): two file neighbours
      open_ = (sum(o) == valence - 1) if None not in o else len(o) == 2
      if i is not None and i not in linked and o and not q and open_:
        missing_caps.append((i, "(no %s residue)" % where))
  caps_on = {}
  for i, k in links:
    if polymer and (atoms[i].name.strip(), atoms[k].name.strip()) in polymer_pairs:
      continue
    caps_on[i] = caps_on.get(i, 0) + 1
  # hydrogens against the restraint file, per heavy atom
  deprotonations, deprotonated, to_add, extra = 0, {}, {}, []
  for i in heavy:
    n = dict_name[i]
    hs = d_h(n)
    replaced = max(0, caps_on.get(i, 0) - n_missing_nb.get(i, 0))
    want = max(0, len(hs) - replaced)
    mine = [k for k in fsc0[i] if k in present_set and is_h(k)]
    if len(mine) == want:
      continue
    diff = dict(atom=n, model=len(mine), restraint_file=want, kind=None, added=[])
    if len(mine) > want:
      diff["kind"] = "extra"
      extra.append("%s: %d H in the model, %d in the restraint file" % (n, len(mine), want))
    elif has_h and i in metal_bound:
      diff["kind"] = "deprotonation"
      deprotonations += want - len(mine)
      deprotonated[i] = want - len(mine)
      result.charge_notes.append("%s: %d restraint-file H missing on a metal-bound "
        "atom, counted as deprotonation" % (n, want - len(mine)))
    else:
      mine_names = set([dict_name.get(k) for k in mine])
      new = [h for h in hs if h not in mine_names][:want - len(mine)]
      to_add[i] = new
      diff["kind"] = "added"
      diff["added"] = new
      if el(i) != "C":
        result.uncertain_atoms.append(n)
        result.uncertain_iseqs.append(i)
    if has_h:
      result.h_differences.append(diff)
  if not has_h:
    result.hydrogens = "restraint file (no H in the model)"
  elif to_add:
    result.hydrogens = "completed from the restraint file"
    result.charge_notes.append("H added from the restraint file: %s" % " ".join(
      [h for i in heavy for h in to_add.get(i, [])]))
  else:
    result.hydrogens = "model"
  result.linked = sorted(linked)
  result.missing_neighbour = sorted(set([i for i, n in missing_caps]))
  result.metal_bound = sorted(metal_bound)
  if extra:
    return fail("H not in the restraint file: %s" % "; ".join(extra))
  # formal charges
  values = [cif_types.formal_charge_and_problem(getattr(a, "charge", None))
    for a in comp.atom_list]
  have_formal = len([q for q, p in values if q is not None]) > 0
  # molecule
  mol = Chem.RWMol()
  conformer = Chem.Conformer()
  xyz = {}
  idx = {}
  def add_atom(element, name, site):
    e = element.capitalize()
    atom = Chem.Atom("H" if e == "D" else e)
    if e == "D":
      atom.SetIsotope(2)
    atom.SetNoImplicit(True)
    atom.SetProp("_Name", name)
    k = mol.AddAtom(atom)
    xyz[k] = site
    return k
  for i in present:
    idx[i] = add_atom(el(i), atoms[i].name.strip(), atoms[i].xyz)
    result.rdkit_to_iseq[idx[i]] = i
  for i, k in bonds:
    mol.AddBond(idx[i], idx[k], Chem.BondType.SINGLE)
  def directions_away(i):
    x = atoms[i].xyz
    return [tuple([x[c] - atoms[k].xyz[c] for c in range(3)])
      for k in fsc0[i] if k in present_set]
  n_on = {}
  def add_h(i, site_dirs, name):
    n_on[i] = n_on.get(i, 0) + 1
    site = _h_site(atoms[i].xyz, site_dirs, n_on[i] - 1, el(i))
    k = add_atom("H", name, site)
    mol.AddBond(idx[i], k, Chem.BondType.SINGLE)
    return k
  # added H before the caps: fragment_mol drops the caps from the end
  added_charge_names = []
  for i in heavy:
    for h in to_add.get(i, []):
      result.added_h.append(add_h(i, directions_away(i), h))
      added_charge_names.append(h)
  for i, k in links:
    d = tuple([atoms[k].xyz[c] - atoms[i].xyz[c] for c in range(3)])
    c = add_h(i, [d], "CAP")
    mol.GetAtomWithIdx(c).SetProp("cap", "linked")
    result.caps.append(dict(index=c, on=i, kind="linked", partner=k))
  for i, n in missing_caps:
    c = add_h(i, directions_away(i), "CAP")
    mol.GetAtomWithIdx(c).SetProp("cap", "missing")
    result.caps.append(dict(index=c, on=i, kind="missing", partner=n))
  for k, site in xyz.items():
    conformer.SetAtomPosition(k, site)
  mol.AddConformer(conformer, assignId=True)
  # total charge over the atoms present
  charge_names = [dict_name[i] for i in present] + added_charge_names
  unknown = [n for n in charge_names if n not in d_atoms]
  if unknown:
    result.charge_notes.append("not in the restraint file: %s" % " ".join(unknown))
  def first_line(e):
    return str(e).strip().splitlines()[0] if str(e).strip() else type(e).__name__
  formal_total = None
  if have_formal:
    formal_total = 0
    for n in charge_names:
      if n in d_atoms:
        q, problem = cif_types.formal_charge_and_problem(getattr(d_atoms[n], "charge", None))
        if q is None:
          result.charge_notes.append("%s: %s, counted as 0" % (n, problem))
        else:
          formal_total += q
    formal_total -= deprotonations
  else:
    partial = [d_atoms[n].partial_charge for n in charge_names if n in d_atoms]
    if [p for p in partial if p is not None]:
      result.charge_notes.append("sum of partial charges %.3f" % sum(
        [p for p in partial if p is not None]))
  # the file's explicit bond orders between heavy atoms present
  rd_of_name = dict([(dict_name[i], idx[i]) for i in present])
  orders = {"sing": 1, "doub": 2, "trip": 3}
  explicit, file_type = {}, {}
  for b in comp.bond_list:
    a1, a2 = b.atom_id_1.strip(), b.atom_id_2.strip()
    t = (b.type or "").strip().lower()[:4]
    if a1 in rd_of_name and a2 in rd_of_name and \
        d_el.get(a1) not in ("H", "D") and d_el.get(a2) not in ("H", "D"):
      file_type[frozenset([a1, a2])] = t
      if t in orders:
        explicit[(a1, a2)] = orders[t]
  rd_order = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
  def disagreements(m):
    got = {}
    for (a1, a2) in explicit:
      bond = m.GetBondBetweenAtoms(rd_of_name[a1], rd_of_name[a2])
      if bond is not None and not bond.GetIsAromatic():
        got[(a1, a2)] = rd_order.get(bond.GetBondType(), 0)
    bad = set([p for p in got if got[p] != explicit[p]])
    # a file double bond drawn by RDKit as a charge-separated single, e.g. a
    # sulfoxide S=O as [S+]-[O-], is the same structure
    for p in list(bad):
      if explicit[p] == 2 and got[p] == 1:
        q = sorted([m.GetAtomWithIdx(rd_of_name[a]).GetFormalCharge() for a in p])
        if q == [-1, 1]:
          bad.discard(p)
    # bonds from one atom to neighbours of one element: orders compared as a set
    sets = {}
    for (a1, a2) in got:
      sets.setdefault((a1, d_el[a2]), []).append((a1, a2))
      sets.setdefault((a2, d_el[a1]), []).append((a1, a2))
    for members in sets.values():
      if len(members) > 1 and sorted([explicit[p] for p in members]) == \
          sorted([got[p] for p in members]):
        bad -= set(members)
    return ["%s-%s: restraint file %d, RDKit %d" % (p[0], p[1], explicit[p], got[p])
      for p in sorted(bad)]
  def attempt(total, order=None):
    """DetermineBondOrders on mol, the atoms renumbered by order (then back)."""
    m = Chem.RWMol(mol if order is None else Chem.RenumberAtoms(mol, order))
    try:
      rdDetermineBonds.DetermineBondOrders(m, charge=total, allowChargedFragments=True,
        embedChiral=True)
      Chem.SanitizeMol(m)
    except Exception as e:
      return None, first_line(e)
    if order is not None:
      back = [0] * len(order)
      for k, i in enumerate(order):
        back[i] = k
      m = Chem.RWMol(Chem.RenumberAtoms(m, back))
    return m, None
  def file_charge(i):
    return cif_types.formal_charge_and_problem(getattr(d_atoms[dict_name[i]], "charge",
      None))[0] if dict_name[i] in d_atoms else None
  def from_file(total, drop=()):
    """mol with the file's bond orders and formal charges (less the deprotonations and the charges on drop); None unless sanitized without radicals at total."""
    m = Chem.RWMol(mol)
    kinds = {"sing": Chem.BondType.SINGLE, "doub": Chem.BondType.DOUBLE,
      "trip": Chem.BondType.TRIPLE, "arom": Chem.BondType.AROMATIC}
    for i, k in bonds:
      if is_h(i) or is_h(k):
        continue
      bond = m.GetBondBetweenAtoms(idx[i], idx[k])
      t = kinds[file_type[frozenset([dict_name[i], dict_name[k]])]]
      bond.SetBondType(t)
      if t == Chem.BondType.AROMATIC:
        bond.SetIsAromatic(True)
        m.GetAtomWithIdx(idx[i]).SetIsAromatic(True)
        m.GetAtomWithIdx(idx[k]).SetIsAromatic(True)
    for i in present:
      q = None if i in drop else file_charge(i)
      m.GetAtomWithIdx(idx[i]).SetFormalCharge((q or 0) - deprotonated.get(i, 0))
    try:
      Chem.SanitizeMol(m)
      Chem.AssignRadicals(m)
    except Exception:
      return None
    if [a for a in m.GetAtoms() if a.GetNumRadicalElectrons()] or \
        Chem.GetFormalCharge(m) != total:
      return None
    return m
  result.search = dict(calls=0, seconds=0.0, valid=[], set_aside=[], disagree={},
    order=None)
  heavy_bonds = [frozenset([dict_name[i], dict_name[k]]) for i, k in bonds
    if not is_h(i) and not is_h(k)]
  all_explicit = not [p for p in heavy_bonds if file_type.get(p) not in
    ("sing", "doub", "trip", "arom")]
  def from_formal(total, drop=()):
    """(mol, source, None) from the file's charges less drop, or (None, None, reason)."""
    if all_explicit:
      m = from_file(total, drop)
      if m is not None:
        return m, "restraint file", None
    # the input atom order, RDKit's canonical order, 10 seeded random orders
    import random
    probe = Chem.Mol(mol)
    probe.UpdatePropertyCache(strict=False)
    ranks = list(Chem.CanonicalRankAtoms(probe, breakTies=True))
    trials = [("input", None), ("canonical", sorted(range(len(ranks)),
      key=lambda i: ranks[i]))]
    for seed in range(10):
      order = list(range(len(ranks)))
      random.Random(seed).shuffle(order)
      trials.append(("random (seed %d)" % seed, order))
    first = None
    for label, order in trials:
      result.search["calls"] += 1
      m, error = attempt(total, order)
      bad = disagreements(m) if m is not None else None
      if first is None:
        first = (error, bad)
      if m is not None and not bad:
        result.search["order"] = label
        return m, "formal charges", None
    error, bad = first
    tail = " (also in the canonical and 10 random atom orders)"
    if error is not None:
      return None, None, "DetermineBondOrders fails for %s with the formal total %d: " \
        "%s%s" % (resname, total, error, tail)
    return None, None, "bond orders disagree with the restraint file for %s at the " \
      "formal total %d: %s%s" % (resname, total, "; ".join(bad), tail)
  m = None
  if formal_total is not None:
    m, source, formal_failed = from_formal(formal_total)
    used_total = formal_total
    if m is None:
      # a charge on a carbon with four bonds cannot be drawn (e.g. GeoStd's -1 on the
      # CH2 next to a sulfonium): the file's other charges, without it
      drop = [i for i in present if el(i) == "C" and file_charge(i) and
        mol.GetAtomWithIdx(idx[i]).GetDegree() == 4]
      if drop:
        used_total = formal_total - sum([file_charge(i) for i in drop])
        m, source, failed = from_formal(used_total, set(drop))
        if m is not None:
          result.charge_notes.append("formal charge ignored on %s (a carbon with four "
            "bonds): total %d, not %d" % (" ".join(["%s (%+d)" % (dict_name[i],
            file_charge(i)) for i in drop]), used_total, formal_total))
    if m is not None:
      if source == "formal charges" and result.search["order"] != "input":
        result.charge_notes.append("DetermineBondOrders succeeded in the %s atom order" %
          result.search["order"])
      mol, total, certain = m, used_total, True
    else:
      # the restraint file's formal charges cannot be built: search the totals as if
      # there were none
      result.charge_notes.append("restraint file formal charges inconsistent, total "
        "searched: %s" % formal_failed)
  if m is None:
    no_formal = "no formal charges" if formal_total is None else \
      "formal charges inconsistent"
    t_search = time.time()
    valid = {}
    totals = range(-10, 11)
    for t in totals:
      result.search["calls"] += 1
      m, e = attempt(t)
      if m is None:
        continue
      bad = disagreements(m)
      if bad:
        result.search["disagree"][t] = bad
      else:
        valid[t] = m
    result.search["seconds"] = time.time() - t_search
    result.search["valid"] = sorted(valid)
    # set aside totals whose structure charges a carbon (e.g. an acetate at -3 as
    # C[C-]([O-])[O-]) when other totals are valid without one
    def charged_carbon(m):
      return [a for a in m.GetAtoms() if a.GetAtomicNum() == 6 and a.GetFormalCharge()]
    plausible = dict([(t, m) for t, m in valid.items() if not charged_carbon(m)])
    if plausible and len(plausible) < len(valid):
      result.search["set_aside"] = sorted([t for t in valid if t not in plausible])
      result.charge_notes.append("totals %s set aside (charged carbon)" % " ".join(
        ["%+d" % t for t in result.search["set_aside"]]))
      valid = plausible
    if len(valid) > 1:
      return fail("ambiguous total charge for %s (%s): %s" % (resname, no_formal,
        " ".join(["%+d" % t for t in sorted(valid)])))
    if not valid:
      d = result.search["disagree"]
      return fail("no valid structure for %s with total charges %+d..%+d (%s)%s" % (
        resname, totals[0], totals[-1], no_formal, (": bond orders disagree with the restraint file at %s"
        % "; ".join(["%+d (%s)" % (t, ", ".join(d[t])) for t in sorted(d)])) if d else ""))
    total = list(valid)[0]
    mol, source, certain = valid[total], "search", False
  result.total_charge = total
  result.total_charge_source = source
  result.charge_certain = certain
  result.mol = mol.GetMol()
  # differences from the restraint file (information)
  for (a1, a2), o in sorted(explicit.items()):
    bond = result.mol.GetBondBetweenAtoms(rd_of_name[a1], rd_of_name[a2])
    if bond is None or bond.GetIsAromatic():
      continue
    got = rd_order.get(bond.GetBondType(), 0)
    if got != o:
      want = {1: "single", 2: "double", 3: "triple"}
      result.differences["bonds"].append("%s-%s: restraint file %s, RDKit %s" % (
        a1, a2, want[o], str(bond.GetBondType()).lower()))
  if have_formal:
    group = dict([(i, i) for i in heavy])
    def find(i):
      while group[i] != i:
        i = group[i]
      return i
    # resonance: same element sharing a heavy neighbour
    for c in heavy:
      nbs = [k for k in fsc0[c] if k in group]
      for a in nbs:
        for b in nbs:
          if a < b and el(a) == el(b):
            group[find(a)] = find(b)
    sums = {}
    for i in heavy:
      q, problem = cif_types.formal_charge_and_problem(getattr(d_atoms.get(dict_name[i]),
        "charge", None)) if dict_name[i] in d_atoms else (None, None)
      s = sums.setdefault(find(i), [0, 0, []])
      s[0] += q or 0
      s[1] += result.mol.GetAtomWithIdx(idx[i]).GetFormalCharge()
      s[2].append(dict_name[i])
    for g, (qd, qr, members) in sorted(sums.items()):
      if qd != qr:
        result.differences["charges"].append("%s: restraint file %+d, RDKit %+d" % (
          " ".join(members), qd, qr))
  # fragment molecule: caps as implicit H; the residue's metals, dative bonds
  frag = Chem.RWMol(result.mol)
  for cap in sorted(result.caps, key=lambda c: -c["index"]):
    x = frag.GetAtomWithIdx(idx[cap["on"]])
    x.SetNumExplicitHs(x.GetNumExplicitHs() + 1)
    frag.RemoveAtom(cap["index"])
  result.fragment_to_iseq = dict(result.rdkit_to_iseq)
  f_idx = dict([(i, idx[i]) for i in present])
  for i in metals:
    if not visible(i):
      continue
    atom = Chem.Atom(el(i).capitalize())
    atom.SetNoImplicit(True)
    atom.SetProp("_Name", atoms[i].name.strip())
    f_idx[i] = frag.AddAtom(atom)
    result.fragment_to_iseq[f_idx[i]] = i
    frag.GetConformer().SetAtomPosition(f_idx[i], atoms[i].xyz)
  for i in metals:
    for k in fsc0[i]:
      if i in f_idx and k in f_idx and not is_h(k) and \
          frag.GetBondBetweenAtoms(f_idx[i], f_idx[k]) is None:
        if el(k) in _metal_elements:
          frag.AddBond(f_idx[i], f_idx[k], Chem.BondType.SINGLE)
        else:
          frag.AddBond(f_idx[k], f_idx[i], Chem.BondType.DATIVE)
  try:
    Chem.SanitizeMol(frag)
    result.fragment_mol = frag.GetMol()
  except Exception as e:
    return fail("sanitization of the fragment molecule failed: %s" % e)
  result.ok = True
  result.seconds = time.time() - t0
  return result

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
