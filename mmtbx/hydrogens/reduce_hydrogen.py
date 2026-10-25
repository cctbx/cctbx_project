from __future__ import absolute_import, division, print_function
import six
import sys, time
from libtbx.utils import Sorry
import mmtbx.model
import mmtbx.monomer_library.server
import iotbx.pdb
import boost_adaptbx.boost.python as bp
from libtbx.utils import null_out
from libtbx import group_args
from scitbx import matrix
import scitbx.cubicle_neighbors
from cctbx.array_family import flex
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_residue_group
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_atom_group
from cctbx import geometry_restraints
from cctbx.geometry_restraints.linking_class import linking_class
#
from cctbx.maptbx.box import shift_and_box_model
# residue part of hd_exchange.site_key(), for the log of renamed input H
from mmtbx.hydrogens.hd_exchange import residue_part, residue_parts
import math
import itertools
from cctbx import crystal
from cctbx import sgtbx
from cctbx.eltbx import covalent_radii

ext = bp.import_ext("cctbx_geometry_restraints_ext")
get_class = iotbx.pdb.common_residue_names_get_class

def get_ligand_interactions(model, dist_min, cutoff_cno, cutoff_sp):
  """
  Finds distance pairs between ligand atoms and non-ligand/other ligand atoms.

  Args:
  model: mmtbx.model.manager object
  dist_min: Minimum distance threshold
  cutoff_cno: Distance cutoff for pairs exclusively containing C/N/O (
    or halogens)
  cutoff_sp: Distance cutoff for pairs containing S or P

  Returns:
      List of tuples containing (i_seq, j_seq) for interacting atoms.
  """
  pdb_hierarchy = model.get_hierarchy()
  xrs = model.get_xray_structure()
  n_seq = xrs.scatterers().size()
  atom_type = [-1] * n_seq
  is_sp = [False] * n_seq
  ligand_id_counter = 1
  get_class = iotbx.pdb.common_residue_names_get_class
  # 1. Pre-computation: Classify all atoms
  for model_ in pdb_hierarchy.models():
    for chain in model_.chains():
      for residue_group in chain.residue_groups():
        # Check for single-atom ions and filter out hydrogens
        non_h_atoms = [a for a in residue_group.atoms()
                       if a.element.strip().upper() not in ["H", "D"]]
        if len(non_h_atoms) == 0: continue # Ignore entirely
        resname = residue_group.unique_resnames()[0]
        r_class = get_class(resname)
        is_water = (r_class == "common_water")
        is_protein = (r_class in ["common_amino_acid", "modified_amino_acid"])
        is_na = (r_class in ["common_rna_dna", "modified_rna_dna",
                             "ccp4_mon_lib_rna_dna"])
        is_single_atom = (len(non_h_atoms) == 1)
        # Skip conditions
        if is_water or is_single_atom: pass # Left as -1
        # Ligand Assignment
        elif not is_protein and not is_na:
          for a in residue_group.atoms():
            e = a.element.strip().upper()
            if e not in ["H", "D"]:
              atom_type[a.i_seq] = ligand_id_counter
              if e in ["S", "P"]:
                is_sp[a.i_seq] = True
          ligand_id_counter += 1
        # Non-ligand (Protein/Nucleic Acid) Assignment
        else:
          for a in residue_group.atoms():
            e = a.element.strip().upper()
            if e not in ["H", "D"]:
              atom_type[a.i_seq] = 0
              if e in ["S", "P"]:
                is_sp[a.i_seq] = True
  # Setup distance comparisons
  max_cutoff = max(cutoff_cno, cutoff_sp)
  max_cutoff_sq = max_cutoff ** 2
  min_cutoff_sq = dist_min ** 2
  cutoff_cno_sq = cutoff_cno ** 2
  cutoff_sp_sq = cutoff_sp ** 2
  # 2. Spatial Search using CCTBX neighbors_fast_pair_generator
  asu_mappings = xrs.asu_mappings(buffer_thickness=max_cutoff)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=max_cutoff
  )
  pairs = []
  # 3. Fast inner loop
  for pair in pair_generator:
    i_seq = pair.i_seq
    j_seq = pair.j_seq
    ti = atom_type[i_seq]
    tj = atom_type[j_seq]
    # Check a: Exclude ignored entities (water, H, single-atoms)
    if ti == -1 or tj == -1: continue
    # Check b: Must involve at least one ligand
    if ti == 0 and tj == 0: continue
    # Check c: Exclude intra-ligand bonds (same ligand ID and same symmetry
    # operator)
    if ti == tj and pair.j_sym == 0: continue
    dist_sq = pair.dist_sq
    # Check d: Minimum distance threshold
    if dist_sq < min_cutoff_sq: continue
    # Check e: Element-specific cutoff check
    if is_sp[i_seq] or is_sp[j_seq]:
      if dist_sq <= cutoff_sp_sq:
        pairs.append((i_seq, j_seq))
    else:
      if dist_sq <= cutoff_cno_sq:
        pairs.append((i_seq, j_seq))
  return pairs

def get_incorrect_hydrogens_for_bond(model, i_seq_A, i_seq_B):
  """
  Given an mmtbx.model.manager and the sequence indices (i_seqs) of two atoms
  presumed to be covalently bonded, determines whether their current protonation
  states are consistent with the formation of that bond.

  It applies a two-step clearance:
  1. Severe Clash Override: Purges any hydrogen physically occupying the
     incoming bond vector.
  2. Valency Quota: Identifies remaining excess hydrogens and removes the ones
     most sterically hindered by the new bond neighborhood.
  """
  atoms = model.get_atoms()
  sites_cart = model.get_sites_cart()
  atom_A = atoms[i_seq_A]
  atom_B = atoms[i_seq_B]
  #
  # Retrieve Restraints Manager and extract connectivity
  grm = model.get_restraints_manager()
  bond_proxies_simple, _ = grm.geometry.get_all_bond_proxies(
    sites_cart=sites_cart)
  #
  connectivity = {}
  for bp in bond_proxies_simple:
    i, j = bp.i_seqs
    connectivity.setdefault(i, []).append(j)
    connectivity.setdefault(j, []).append(i)
  #
  mon_lib_srv = model.get_mon_lib_srv()
  #
  def get_bond_order(a1, a2):
    if a1.parent().id_str() == a2.parent().id_str():
      resname = a1.parent().resname.strip().upper()
      comp = mon_lib_srv.get_comp_comp_id_direct(resname)
      if comp is not None:
        for bond in comp.bond_list:
          id1, id2 = bond.atom_id_1.strip(), bond.atom_id_2.strip()
          n1, n2 = a1.name.strip(), a2.name.strip()
          if (id1 == n1 and id2 == n2) or (id1 == n2 and id2 == n1):
            btype = bond.type.lower()
            if 'double' in btype: return 2.0
            if 'triple' in btype: return 3.0
            if 'deloc' in btype or 'aromatic' in btype: return 1.5
            return 1.0
    return 1.0
  #
  def get_formal_charge(atom):
    c_str = atom.charge.strip()
    if c_str:
      try:
        rv = int(c_str[-1]+c_str[:-1]) if c_str[-1] in ['+','-'] else int(c_str)
        return rv
      except ValueError:
        pass
    resname = atom.parent().resname.strip().upper()
    comp = mon_lib_srv.get_comp_comp_id_direct(resname)
    if comp is not None:
      for a in comp.atom_list:
        if a.atom_id.strip() == atom.name.strip():
          if hasattr(a, 'charge'):
            c_lib = str(a.charge).strip()
            if c_lib and c_lib != '.':
              try:
                return int(c_lib[-1]+c_lib[:-1]) if c_lib[-1] in ['+','-'] else int(c_lib)
              except ValueError:
                pass
    return 0
  #
  def get_ideal_valence(atom):
    el = atom.element.strip().upper()
    charge = get_formal_charge(atom)
    if el == 'C': return 4
    if el == 'N': return 3 + charge
    if el == 'O': return 2 + charge
    if el == 'S': return 2 + charge
    if el == 'P': return 5
    if el in ['F', 'CL', 'BR', 'I']: return 1
    return 0
  #
  excess_h_iseqs = []
  #
  # Distance threshold for an impossible geometric overlap (e.g. H pointing
  # directly into the incoming heavy atom)
  SEVERE_CLASH_DIST = 1.5
  #
  for target_iseq in [i_seq_A, i_seq_B]:
    target_atom = atoms[target_iseq]
    ideal_val = get_ideal_valence(target_atom)
    if ideal_val == 0: continue
    neighbors = connectivity.get(target_iseq, [])
    heavy_order_sum = 0.0
    h_neighbors = []
    for n_iseq in neighbors:
      n_atom = atoms[n_iseq]
      if n_atom.element.strip().upper() in ['H', 'D', 'T']:
        h_neighbors.append(n_atom)
      else:
        heavy_order_sum += get_bond_order(target_atom, n_atom)
    other_iseq = i_seq_B if target_iseq == i_seq_A else i_seq_A
    other_site = sites_cart[other_iseq]
    if other_iseq not in neighbors: heavy_order_sum += 1.0
    heavy_order_rounded = int(math.floor(heavy_order_sum))
    expected_h = max(0, ideal_val - heavy_order_rounded)
    # ----------------------------------------------------------------------
    # STEP 1: Severe Clash Override
    # If a hydrogen occupies the incoming bond's vector it must be removed,
    # even if removing it drops the atom below its expected valency.
    # ----------------------------------------------------------------------
    surviving_h_neighbors = []
    for h_atom in h_neighbors:
      h_site = sites_cart[h_atom.i_seq]
      dist_to_other = math.sqrt((h_site[0] - other_site[0])**2 +
                                (h_site[1] - other_site[1])**2 +
                                (h_site[2] - other_site[2])**2)
      if dist_to_other < SEVERE_CLASH_DIST:
        # Severe clash detected; mark for removal immediately
        excess_h_iseqs.append(h_atom.i_seq)
      else:
        surviving_h_neighbors.append(h_atom)
    # ----------------------------------------------------------------------
    # STEP 2: Valency Quota Check
    # Now evaluate ONLY the surviving hydrogens against the valency limit.
    # ----------------------------------------------------------------------
    current_h_count = len(surviving_h_neighbors)
    if current_h_count > expected_h:
      num_to_remove = int(current_h_count - expected_h)
      clash_set_iseqs = [other_iseq]
      iters=connectivity.get(target_iseq, []) + connectivity.get(other_iseq, [])
      for n_iseq in iters:
        if n_iseq != target_iseq and n_iseq != other_iseq:
          if atoms[n_iseq].element.strip().upper() not in ['H', 'D', 'T']:
            clash_set_iseqs.append(n_iseq)
      def min_distance_to_clash_set(h_atom):
        h_site = sites_cart[h_atom.i_seq]
        min_dist = float('inf')
        for c_iseq in clash_set_iseqs:
          c_site = sites_cart[c_iseq]
          dist = math.sqrt((h_site[0] - c_site[0])**2 +
                           (h_site[1] - c_site[1])**2 +
                           (h_site[2] - c_site[2])**2)
          if dist < min_dist: min_dist = dist
        return min_dist
      surviving_h_neighbors.sort(
        key=lambda h: (min_distance_to_clash_set(h), h.i_seq))
      excess_h_iseqs.extend(
        [h.i_seq for h in surviving_h_neighbors[:num_to_remove]])
  #
  return excess_h_iseqs

def workaround_003(model):
  pairs = get_ligand_interactions(
    model=model, dist_min=1.1, cutoff_cno=1.6, cutoff_sp=1.9)
  remove_selection = []
  for pair in pairs:
    badH_i_seqs = get_incorrect_hydrogens_for_bond(
      model   = model,
      i_seq_A = pair[0],
      i_seq_B = pair[1])
    if len(badH_i_seqs)>0:
      remove_selection.extend(badH_i_seqs)
  removed = 0
  if len(remove_selection)>0:
    badH_i_seqs = flex.size_t(badH_i_seqs)
    removed = badH_i_seqs.size()
    keep_selection = ~flex.bool(model.size(), badH_i_seqs)
    model = model.select(keep_selection)
  return model, removed

_ener_lib_types = None

def get_ccd_energy_types(cc_cif):
  '''
  Energy types for a CCD entry, which has none: element plus the H the CCD puts
  on the atom (free, neutral molecule). Without a type the probe atom info gets
  no radius and no acceptor (2e4i BGM). A type ener_lib does not know -> None.
  '''
  global _ener_lib_types
  if _ener_lib_types is None:
    from mmtbx.monomer_library import server
    _ener_lib_types = set(server.ener_lib().lib_atom.keys())
  element = dict((a.atom_id, a.type_symbol.strip().upper())
    for a in cc_cif.get('_chem_comp_atom', []))
  neighbors = {}
  for b in cc_cif.get('_chem_comp_bond', []):
    neighbors.setdefault(b.atom_id_1, []).append(b.atom_id_2)
    neighbors.setdefault(b.atom_id_2, []).append(b.atom_id_1)
  result = {}
  for name, e in element.items():
    bonded = [element.get(x) for x in neighbors.get(name, [])]
    n_h = sum(1 for x in bonded if x in ('H', 'D'))
    if e in ('H', 'D'):
      parent = bonded[0] if len(bonded) == 1 else None
      t = {'O': 'HOH1', 'N': 'HNH1', 'S': 'HSH1'}.get(parent, 'HCH')
    elif e == 'O': t = 'OH1' if n_h else 'O'
    elif e == 'N': t = {0: 'N', 1: 'NH1'}.get(n_h, 'NH2')
    elif e == 'S': t = 'SH1' if n_h else 'S'
    else:          t = e
    result[name] = t if t in _ener_lib_types else None
  return result

_library_xh = None

def _get_library_xh():
  '''
  X-H lengths of the monomer library, read from its own entries:
  {'xray': {...}, 'neutron': {...}}, each keyed 'C' (sp3 C-H), 'C_aromatic'
  (aromatic C-H), 'N' (N-H), 'O' (O-H), 'S' (S-H). The library has no X-H
  for other elements.
  '''
  global _library_xh
  if _library_xh is None:
    from mmtbx.monomer_library import server
    srv = server.server()
    sources = {
      'C':          ('ALA', 'CA', 'HA'),
      'C_aromatic': ('PHE', 'CD1', 'HD1'),
      'N':          ('ALA', 'N', 'H'),
      'O':          ('SER', 'OG', 'HG'),
      'S':          ('CYS', 'SG', 'HG')}
    xray, neutron = {}, {}
    for key, (code, x, h) in sources.items():
      for b in srv.get_comp_comp_id_direct(code).bond_list:
        if set((b.atom_id_1, b.atom_id_2)) == set((x, h)):
          xray[key] = float(b.value_dist)
          neutron[key] = float(b.value_dist_neutron)
    assert sorted(xray) == sorted(neutron) == sorted(sources)
    _library_xh = {'xray': xray, 'neutron': neutron}
  return _library_xh

def get_library_xray_xh():
  '''X-ray X-H lengths of the monomer library (see _get_library_xh)'''
  return _get_library_xh()['xray']

def get_library_neutron_xh():
  '''Neutron X-H lengths of the monomer library (see _get_library_xh)'''
  return _get_library_xh()['neutron']

def get_h_restraints(resname, strict=True):
  from mmtbx.monomer_library import cif_types
  from mmtbx.chemical_components import get_cif_dictionary
  from mmtbx.ligands.rdkit_utils import get_molecule_from_resname
  molecule = get_molecule_from_resname(resname)
  if molecule is None: return None
  cc_cif = get_cif_dictionary(resname)
  cc = cc_cif['_chem_comp'][0]
  hs = []
  hsi = []
  chem_comp = cif_types.chem_comp(
    id=cc.id,
    three_letter_code=cc.three_letter_code,
    name=cc.name,
    group=cc.type,
    number_atoms_all=0, #cc.number_atoms_all,
    number_atoms_nh=0, #cc.number_atoms_nh,
    desc_level=".")
  comp_comp_id = cif_types.comp_comp_id(source_info=None, chem_comp=chem_comp)
  lookup = {}
  energy_types = get_ccd_energy_types(cc_cif)
  for i, a in enumerate(cc_cif.get('_chem_comp_atom',[])):
    lookup[a.atom_id]=i
    lookup[i]=a.atom_id
    if a.type_symbol in ['H', 'D']:
      hs.append(a.atom_id)
      hsi.append(i)
    comp_comp_id.atom_list.append(cif_types.chem_comp_atom(
      atom_id=a.atom_id,
      type_symbol=a.type_symbol,
      type_energy=energy_types.get(a.atom_id),
      # partial_charge=a.partial_charge,
      ))
  conf =  molecule.GetConformer()
  from rdkit import Chem # needed import
  from rdkit.Chem import rdMolTransforms
  for b in cc_cif.get('_chem_comp_bond',[]):
    if strict:
      if (b.atom_id_1 not in hs and
          b.atom_id_2 not in hs): continue
    if ( b.atom_id_1 not in lookup or
         b.atom_id_2 not in lookup): continue
    atom_idx1=lookup[b.atom_id_1]
    atom_idx2=lookup[b.atom_id_2]
    bl = rdMolTransforms.GetBondLength(conf, atom_idx1, atom_idx2)
    # X-H: the library's X-ray and neutron lengths for the heavy atom
    # (aromatic from the RDKit molecule); elements the library has no X-H
    # for: 0.9 x the CCD length. Heavy-atom bonds: the CCD length
    # (placement does not move heavy atoms, but the grm stays on the model)
    if b.atom_id_1 in hs or b.atom_id_2 in hs: bl *= .9
    value_dist_neutron = None
    if (b.atom_id_1 in hs) != (b.atom_id_2 in hs):
      x_idx = atom_idx2 if b.atom_id_1 in hs else atom_idx1
      x_atom = molecule.GetAtomWithIdx(x_idx)
      key = x_atom.GetSymbol().upper()
      if key == 'C' and x_atom.GetIsAromatic(): key = 'C_aromatic'
      xray = get_library_xray_xh().get(key)
      if xray is not None: bl = xray
      neutron = get_library_neutron_xh().get(key)
      if neutron is not None: value_dist_neutron = '%0.3f' % neutron
    comp_comp_id.bond_list.append(cif_types.chem_comp_bond(
      atom_id_1=b.atom_id_1,
      atom_id_2=b.atom_id_2,
      type=b.value_order,
      value_dist='%0.3f' % bl,
      value_dist_esd=".1",
      value_dist_neutron=value_dist_neutron))

  from mmtbx.ligands.rdkit_utils import enumerate_angles
  for angle in enumerate_angles(molecule):
    if strict:
      if angle[0] in hsi or angle[2] in hsi:
        av = rdMolTransforms.GetAngleDeg(conf, angle[0], angle[1], angle[2])
      else: continue
    else:
      av = rdMolTransforms.GetAngleDeg(conf, angle[0], angle[1], angle[2])
    if ( angle[0] not in lookup or
         angle[2] not in lookup): continue
    comp_comp_id.angle_list.append(cif_types.chem_comp_angle(
      atom_id_1=lookup[angle[0]],
      atom_id_2=lookup[angle[1]],
      atom_id_3=lookup[angle[2]],
      value_angle='%0.1f' % av,
      value_angle_esd="1"))

  from mmtbx.ligands.rdkit_utils import enumerate_torsions
  for i, angle in enumerate(enumerate_torsions(molecule)):
    if strict:
      if angle[0] in hsi or angle[3] in hsi:
        av = rdMolTransforms.GetDihedralDeg(conf, angle[0], angle[1], angle[2], angle[3])
      else: continue
    else:
      av = rdMolTransforms.GetDihedralDeg(conf, angle[0], angle[1], angle[2], angle[3])
    if ( angle[0] not in lookup or
         angle[3] not in lookup): continue
    comp_comp_id.tor_list.append(cif_types.chem_comp_tor(
      id='Var_%03d' % i,
      atom_id_1=lookup[angle[0]],
      atom_id_2=lookup[angle[1]],
      atom_id_3=lookup[angle[2]],
      atom_id_4=lookup[angle[3]],
      value_angle='%0.1f' % av,
      value_angle_esd='1',
      period='1'))
  return comp_comp_id

def atom_in_restraints(name, cc_cif):
  for a in cc_cif.get('_chem_comp_atom',[]):
    if a.atom_id==name:
      return a
  return None

def chain_start_residue_groups(model, peptide_bond_cutoff=2.0):
  '''
  Residue groups that start a chain or follow a chain break: first residue
  group of each chain, and each residue group whose N is not peptide-bonded
  (C(i-1)-N(i) distance >= peptide_bond_cutoff, same as add_terminal_oxt) to
  C of the preceding residue group. Alternative conformations: bonded if any
  C-N pair is bonded. Cost is linear in number of residue groups.
  '''
  cutoff_sq = peptide_bond_cutoff**2
  result = []
  for chain in model.chains():
    prev_c_sites = None
    for rg in chain.residue_groups():
      n_sites, c_sites = [], []
      for ag in rg.atom_groups():
        n = ag.get_atom('N')
        if n: n_sites.append(n.xyz)
        c = ag.get_atom('C')
        if c: c_sites.append(c.xyz)
      bonded = False
      if prev_c_sites:
        for a in prev_c_sites:
          for b in n_sites:
            if (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 < cutoff_sq:
              bonded = True
              break
          if bonded: break
      if not bonded:
        result.append(rg)
      prev_c_sites = c_sites
  return result

# ==============================================================================

def heavy_neighbors_in_restraints(atom, mon_lib_srv):
  '''
  Names of the non-H atoms bonded to atom in the dictionary placement uses,
  else in the CCD. None if neither knows the residue.
  '''
  name = atom.name.strip()
  mlq = None
  if mon_lib_srv is not None: # None: propeller called before placement
    mlq = mon_lib_query(residue=atom.parent(), mon_lib_srv=mon_lib_srv,
      raise_sorry=False)[0]
  if mlq is not None:
    types = dict((a.atom_id, a.type_symbol) for a in mlq.atom_list)
    bonds = [(b.atom_id_1, b.atom_id_2) for b in mlq.bond_list]
  else:
    from mmtbx.chemical_components import get_cif_dictionary
    cc_cif = get_cif_dictionary(atom.parent().resname)
    if cc_cif is None: return None
    types = {}
    for b in cc_cif.get('_chem_comp_bond',[]):
      for i in (b.atom_id_1, b.atom_id_2):
        a = atom_in_restraints(i, cc_cif)
        if a is not None: types[i] = a.type_symbol
    bonds = [(b.atom_id_1, b.atom_id_2)
             for b in cc_cif.get('_chem_comp_bond',[])]
  rc = []
  for b in bonds:
    for this, other in (b, b[::-1]):
      if (this.strip() == name
          and types.get(other, '').strip().upper() not in ['H', 'D']):
        rc.append(other)
  return rc

# ==============================================================================

def _terminal_h(mlq, atom_dict):
  '''
Identify hydrogens in a peptide-like component that belong only to the
free form: the hydrogen on OXT and all but one hydrogen on N.
Within a chain, N forms a peptide bond with the preceding residue;
at a true N-terminus, place_n_terminal_propeller adds the terminal
hydrogens later.

Identify these hydrogens by bond connectivity rather than atom names,
which vary. For example, 557 geostd entries have two N-bound hydrogens but
no atom named 'H2' (e.g. 9FZ uses H3/H4).
  '''
  on = {'N': [], 'OXT': []}
  for b in mlq.bond_list:
    for x, y in [(b.atom_id_1, b.atom_id_2), (b.atom_id_2, b.atom_id_1)]:
      a = atom_dict.get(y)
      if x in on and a is not None and a.type_symbol in ('H', 'D'):
        on[x].append(y)
  remove = list(on['OXT'])
  if len(on['N']) > 1:
    # keep the amide H: 'H'/'D' if the dictionary has it, else the first one
    order = [a.atom_id for a in mlq.atom_list]
    keep = sorted(on['N'], key=lambda h: (h not in ('H', 'D'), order.index(h)))[0]
    remove.extend(h for h in on['N'] if h != keep)
  return remove

# metals do not count as bonded to N (coordination, not a covalent bond)
nonmetals = ('C', 'N', 'O', 'S', 'P', 'SE', 'B', 'SI', 'F', 'CL', 'BR', 'I')

def _blank_amide_group(rg, prev_rg):
  '''
  The blank atom group of rg if its amide H depends on no split atom: N and
  CA only in the blank group, exactly one C of the previous residue bonded to
  N, and no H on N yet. That H then goes once into the blank group instead of
  into each altloc (6bhe VAL A 274: side chain in A and B, GLN A 273 in A
  only, so H/B had no C). Without a bonded C the H turns with the split side
  chain and stays split. None otherwise.
  '''
  blank = [ag for ag in rg.atom_groups() if ag.altloc == '']
  if not blank or prev_rg is None: return None
  for ag in rg.atom_groups():
    if ag.altloc == '': continue
    if ag.get_atom('N') is not None or ag.get_atom('CA') is not None:
      return None
  for ag in rg.atom_groups():
    if ag.get_atom('H') is not None: return None
  n = blank[0].get_atom('N')
  if n is None or blank[0].get_atom('CA') is None: return None
  cs = [ag.get_atom('C') for ag in prev_rg.atom_groups()]
  cs = [c for c in cs if c is not None and c.distance(n) < 2.0]
  if len(cs) != 1: return None
  return blank[0]

def _split_ca_groups(rg, n):
  '''
  The atom groups with CA if n has no altloc, CA is in two or more altlocs
  and none blank, and the blank group has no N-terminal H. H1/H2/H3 then
  depend on the split CA and go into each of these groups (9grf ALA M 1: in
  the blank group they followed CA A only, and riding H stopped at CA B).
  Empty otherwise.
  '''
  if n.parent().altloc != '': return []
  ags = rg.atom_groups()
  for ag in ags:
    if ag.altloc != '': continue
    if ag.get_atom('CA') is not None: return []
    for name in ('H', 'H1', 'H2', 'H3', 'D', 'D1', 'D2', 'D3'):
      if ag.get_atom(name) is not None: return []
  result = [ag for ag in ags if ag.get_atom('CA') is not None]
  return result if len(result) > 1 else []

def _bonded_to_other_residue(atom, rg, model_atoms, model_sites):
  '''
  True if a non-metal atom of another residue is within bonding distance of
  atom.
  '''
  d2 = (model_sites - atom.xyz).dot()
  for i in (d2 < 1.75**2).iselection():
    a = model_atoms[i]
    if (a.element.strip().upper() in nonmetals
        and a.parent().parent().memory_id() != rg.memory_id()):
      return True
  return False

def _na_terminal_h(ag, rg, expected_h, model_atoms, model_sites):
  '''
  The terminal hydroxyl H of a standard nucleotide that the polymer
  dictionaries pdb_interpretation uses for DNA/RNA chains lack: HO3' when
  O3' is present and bonds to no other residue (a 3' end), HO5' when the
  residue has no P and its O5' bonds to no other residue (a 5'-OH end).
  Neither is a protonation question (pKa above 12), and pdb_interpretation
  restrains both. The atoms are looked up in ag's conformer: the atom groups
  of rg with a blank altloc or ag's (all of them for a blank ag), so a
  backbone in the blank group is seen from a split base. expected_h: the
  dictionary H names, not repeated.
  '''
  altloc = ag.altloc.strip()
  conformer = {}
  for g in rg.atom_groups():
    if altloc and g.altloc.strip() not in ('', altloc): continue
    for a in g.atoms():
      conformer.setdefault(a.name.strip().upper(), a)
  names = []
  for o_name, h_name in (("O3'", "HO3'"), ("O5'", "HO5'")):
    if h_name in expected_h or o_name not in conformer: continue
    if o_name == "O5'" and 'P' in conformer: continue
    if not _bonded_to_other_residue(conformer[o_name], rg, model_atoms,
                                    model_sites):
      names.append(h_name)
  return names

def _dictionary_h_on_n(ag, mon_lib_srv):
  '''
  False if the dictionary bonds "H" to an atom other than N (geostd ASJ: C).
  '''
  if get_class(name=ag.resname) == 'common_amino_acid': return True
  mlq = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv, raise_sorry=False)[0]
  if mlq is None: return True
  for b in mlq.bond_list:
    if 'H' in (b.atom_id_1, b.atom_id_2):
      return 'N' in (b.atom_id_1, b.atom_id_2)
  return True

def mon_lib_query(residue, mon_lib_srv, construct_h_restraints=True, raise_sorry=True):
  # if get_class(residue.resname) in ['common_rna_dna']:
  #   md = get_h_restraints(residue.resname)
  #   return md
  # if print_time: print(residue.resname, get_class(residue.resname))
  if residue.resname == 'UNL':
    return None, None
  md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
    residue_name=residue.resname,
    atom_names=residue.atoms().extract_name())
  cif_object=None
  # if md is None:
  #   md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
  #     residue_name='%s_EL' % residue.resname,
  #     atom_names=residue.atoms().extract_name(),
  #     ad_hoc_single_atom_residues=True)
  if md is None:
    md = get_h_restraints(residue.resname, strict=False)
    if md is None:
      if raise_sorry:
        raise Sorry('Entity "%s" not found in CCD (or GeoStd). Please supply restraints.' % residue.resname)
      else:
        return None, None
    from six.moves import cStringIO as StringIO
    input_string='data_comp_list\n'
    input_string+=str(md.chem_comp.as_cif_loop())
    f=StringIO()
    md.show(f=f)
    # use strip in case 3-letter code has only 2 letters (e.g. DI)
    input_string += '\ndata_comp_%s\n' % residue.resname.strip()
    input_string += '\n%s' % f.getvalue()
    cif_object = iotbx.cif.reader(input_string=input_string).model()
  return md, cif_object

def load_restraint_objects_into_server(model, mon_lib_srv):
  '''
  Load the model's restraint CIF objects (user CIFs, and the CCD objects
  reduce2 attaches) into mon_lib_srv, so that mon_lib_query() finds these
  residues there and does not build a CCD placeholder for them
  (get_h_restraints(): placement-only restraints with library X-H lengths
  and approximate energy types) that would override the user's dictionary.
  '''
  ro = model.get_restraint_objects()
  if ro:
    for fname, cif_obj in ro:
      if cif_obj is not None:
        try:
          mon_lib_srv.process_cif_object(
            cif_object=cif_obj, file_name=fname)
        except Exception:
          pass

def with_auto_restraint_object(ro, resname, cif_object,
                               auto_restraint_names):
  '''
  New list: the restraint objects ro plus the dictionary that mon_lib_query()
  generated from the CCD for resname, as auto_<resname>. ro is not changed:
  select() and deep_copy() pass the same list to the copy.
  '''
  auto_name = 'auto_%s' % resname
  if auto_name not in auto_restraint_names:
    auto_restraint_names.append(auto_name)
  return list(ro or []) + [(auto_name, cif_object)]

def add_auto_restraint_object(model, resname, cif_object,
                              auto_restraint_names):
  '''
  Add the restraint dictionary that mon_lib_query() generated from the CCD
  for resname to the model's restraint objects, as auto_<resname>. The model
  gets a new list (see with_auto_restraint_object).
  '''
  model.set_restraint_objects(with_auto_restraint_object(
    model.get_restraint_objects(), resname, cif_object, auto_restraint_names))

def remove_auto_restraint_objects(model, auto_restraint_names):
  '''
  Remove the auto-generated restraint dictionaries from the model's restraint
  objects without unsetting its restraints manager (see
  place_hydrogens.remove_auto_restraint_objects).
  '''
  if not auto_restraint_names: return
  ro = model.get_restraint_objects()
  if not ro: return
  model._restraint_objects = [
    (name, obj) for name, obj in ro
    if name not in auto_restraint_names]

def add_ccd_restraint_objects(model):
  '''
  For a model with H that is interpreted without placement (reduce2
  approach=optimize): one mon_lib_query() per residue name, as placement does
  per residue. Residues found only in the CCD get placement's CCD-generated
  restraint objects; drop them after the interpretation with
  remove_auto_restraint_objects().
  Returns group_args(no_restraints (residue names without restraints, as
  placement's no_H_placed_mlq), auto_restraint_names).
  '''
  mon_lib_srv = model.get_mon_lib_srv()
  load_restraint_objects_into_server(model, mon_lib_srv)
  no_restraints = set()
  auto_restraint_names = []
  seen = set()
  for ag in model.get_hierarchy().atom_groups():
    if ag.resname in seen: continue
    seen.add(ag.resname)
    mlq, cif_object = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv,
      raise_sorry=False)
    if mlq is None:
      no_restraints.add(ag.resname)
    elif cif_object:
      add_auto_restraint_object(model, ag.resname, cif_object,
        auto_restraint_names)
  return group_args(no_restraints=no_restraints,
    auto_restraint_names=auto_restraint_names)

# ==============================================================================

def residue_label(atom):
  '''
  'resname chain resseq' of an atom's residue, as in the reports.
  '''
  rg = atom.parent().parent()
  return '%s %s %s' % (atom.parent().resname.strip(),
    rg.parent().id.strip(), rg.resseq.strip())

def atoms_not_in_dictionary(model):
  '''
  Heavy atoms whose name is not in their residue's dictionary, in a processed
  model (energy type False; a residue without restraints has ''). Typically a
  dictionary older than the CCD (7a5m 2L5: geostd N4/C12/C11/O4, model
  N/CA/C/O).
  :return: [(residue_label(), [atom names])], in model order.
  '''
  te = getattr(model, '_type_energies', None)
  if te is None: return []
  found, order = {}, []
  for a in model.get_hierarchy().atoms():
    if a.element_is_hydrogen() or te[a.i_seq] not in (False, 'False'):
      continue
    resid = residue_label(a)
    if resid not in found:
      found[resid] = []
      order.append(resid)
    if a.name.strip() not in found[resid]:
      found[resid].append(a.name.strip())
  return [(r, found[r]) for r in order]

def atoms_not_in_dictionary_text(residues, n_h, h_text, h_labels=(),
                                 cap=12):
  '''
  Report of atoms_not_in_dictionary() (placement and reduce2). h_text says
  what happened to the n_h H on them, e.g. 'H atoms on them were not placed'
  (no clause if n_h is 0). At most cap residues are listed, all h_labels.
  '''
  header = '%d residue%s' % (len(residues), 's' if len(residues) > 1 else '')
  if n_h: header += '; %d %s' % (n_h, h_text)
  lines = ['Atom names not in the restraint dictionary of their residue '
    '(dictionary older', 'than the CCD?) - %s:' % header]
  for resid, names in residues[:cap]:
    lines.append(' %s: %s' % (resid, ' '.join(names)))
  if len(residues) > cap:
    lines.append(' (+%d more)' % (len(residues) - cap))
  lines.extend('   ' + label.strip() for label in h_labels)
  resnames = []
  for resid, names in residues:
    if resid.split()[0] not in resnames: resnames.append(resid.split()[0])
  lines.append('Rename these atoms to match the dictionary, or supply a '
    'restraint CIF for %s.' % ', '.join(resnames))
  return '\n'.join(lines)

def get_output_crystal_symmetry(model):
  '''
  Crystal symmetry to write out for a model that is boxed for H placement,
  in the precedence of the model's writer (_figure_out_cs_to_output()): the
  model's unit cell symmetry if set (e.g. the placeholder on a model
  returned by reduce2, also after process()); else its own symmetry; else
  the CRYST1 1 1 1 P 1 placeholder if the input file had a cell record that
  cctbx discards as nonsense (cryo-EM convention); else None. Call before
  process(), which boxes a model that has no symmetry.
  '''
  if model.unit_cell_crystal_symmetry() is not None:
    return model.unit_cell_crystal_symmetry()
  cs = model.crystal_symmetry()
  if (cs is not None) and (cs.unit_cell() is not None):
    return cs
  model_input = model.get_model_input()
  if model_input is None:
    return None
  if hasattr(model_input, 'cif_block'):
    has_cell = '_cell.length_a' in model_input.cif_block
  else:
    has_cell = any(line.startswith('CRYST1')
                   for line in model_input.crystallographic_section())
  if not has_cell:
    return None
  return crystal.symmetry((1, 1, 1, 90, 90, 90), 'P 1')

# ==============================================================================
# H on atoms bound to metals: thiols (A4), backbone amide N and T/U N3 (A5)

# alkali, alkaline earth, transition metals (Zn, Cd, Hg named), lanthanides,
# actinides, Al Ga In Sn Tl Pb Bi; not B Si Ge As Sb Te. Membership makes an
# atom a coordination candidate only. Left out: Rf-Cn and Fm-Lr, which never
# occur in models, and Sg/No would match atom names when the element is
# guessed from the name.
METAL_ELEMENTS = frozenset((
  'LI NA K RB CS FR BE MG CA SR BA RA '
  'SC TI V CR MN FE CO NI CU ZN Y ZR NB MO TC RU RH PD AG CD '
  'HF TA W RE OS IR PT AU HG '
  'LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU '
  'AC TH PA U NP PU AM CM BK CF ES '
  'AL GA IN SN TL PB BI').split())

# N-H next to these: alkali and alkaline-earth contacts are mostly ionic, with
# no evidence for deprotonation (A5)
NH_METAL_ELEMENTS = METAL_ELEMENTS - frozenset(
  'LI NA K RB CS FR BE MG CA SR BA RA'.split())

# images of one metal atom closer than this are one metal site
METAL_SITE_IMAGE_DISTANCE = 2.0

# X-metal window around the sum of the cctbx covalent radii (provisional: it
# covers every S- and N-metal contact of the A4 B0 and A5 models; Ni, Co,
# Cd, Hg untested for S)
THIOL_METAL_WINDOW = (-0.55, 0.25)

THIOL_H_SUPPLIED_NOTE = ('kept as supplied; thiols bound to metals are usually '
  'deprotonated')

METAL_NH_SUPPLIED_NOTE = ('kept as supplied; backbone amide N and '
  'thymine/uracil N3 bound to metals are usually deprotonated')

# thymine/uracil residues whose N3 carries H3
THYMINE_URACIL = frozenset(('DT', 'T', 'U', 'DU'))

def _atom_element(atom):
  e = atom.element.strip().upper()
  if not e:
    e = (atom.determine_chemical_element_simple() or '').strip().upper()
  return e

def _atom_label(atom):
  return atom.id_str().replace('pdb=', '').replace('"', '').strip()

def _covalent_radius(element):
  try:
    return covalent_radii.table(element.title()).radius()
  except (RuntimeError, ValueError):
    return None

def contact_symmetry(crystal_symmetry):
  '''
  The symmetry for contacts across symmetry: None without a cell or with the
  CRYST1 1 1 1 P 1 placeholder.
  '''
  cs = crystal_symmetry
  if (cs is None or cs.unit_cell() is None or cs.space_group_info() is None):
    return None
  if all(abs(x - 1) < 1.e-3 for x in cs.unit_cell().parameters()[:3]):
    return None
  return cs

def _compatible(a, b):
  alt_a, alt_b = a.parent().altloc, b.parent().altloc
  return alt_a == '' or alt_b == '' or alt_a == alt_b

class metal_h_sites(object):
  '''
  X-metal contacts of the H-carrying atoms X that a subclass selects in a
  processed model, and the protonation verdict for each H site on them.

  Candidates (subclass, _candidates()): (i_x, i_h, judged); judged False
  gives contacts only (reported). Heavy neighbours come from the bond
  restraints, not counting metals or metal-coordination proxies.
  Contacts: bond proxies with origin "metal coordination" (simple and asu), and
  atoms of self.metals within THIOL_METAL_WINDOW of the summed cctbx
  covalent radii; same MODEL; across crystal_symmetry (None: no symmetry).
  Verdict per H site, from the metals of a compatible conformer (blank or the
  same altloc) with occupancy > 0. One metal site is a connected component of
  the images of one metal atom (one name in one residue_group), with an edge
  between two images closer than METAL_SITE_IMAGE_DISTANCE, whatever the
  symmetry operator; its coverage is the sum of the occupancies of its
  distinct altlocs (an altloc counted once), and the largest coverage is
  compared with the population the H site represents: X's occupancy if
  the H is blank or X is in the H's altloc; for a blank X and an H in
  altloc X, the occupancy of the heavy atoms of conformer X.
    self.covered coverage >= population (1.e-6)
    unresolved   coverage < population; the H's group has no heavy atoms;
                 the heavy atoms of conformer X differ in occupancy
    not counted  no such metal (incompatible altloc, occupancy 0)

  contacts: group_args(i_x, i_metal, sym_op, image, distance, evidence)
  h_sites:  group_args(i_h, i_x, verdict, reason, population, coverage,
                       i_metal (of the largest site), contacts, proxy,
                       sites: group_args(key, contacts, occ by altloc))
  reported: group_args(i_h, i_x, contacts) of the candidates not judged,
            with counted contacts
  '''
  metals = METAL_ELEMENTS
  covered = None
  state = None

  def _candidates(self, atoms, elements, heavy, h_on):
    raise NotImplementedError

  def __init__(self, model, crystal_symmetry=None):
    self.contacts = []
    self.h_sites = []
    self.reported = []
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    self.atoms = atoms
    elements = [_atom_element(a) for a in atoms]
    model_index = {}
    for i_model, m in enumerate(hierarchy.models()):
      for a in m.atoms(): model_index[a.i_seq] = i_model
    grm = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart()
    simple, asu = grm.get_all_bond_proxies(sites_cart=sites_cart)
    metal_origin = linking_class().get_origin_id('metal coordination')
    asu_mappings = None
    if asu.size() > 0:
      asu_mappings = grm.pair_proxies().bond_proxies.asu_mappings()
    is_h = lambda i: elements[i] in ('H', 'D')
    heavy, h_on, proxied = {}, {}, {}
    def add(i, j, origin_id, rt_mx_ji):
      for x, y, rt in ((i, j, rt_mx_ji), (j, i, rt_mx_ji.inverse())):
        if origin_id == metal_origin:
          proxied.setdefault(x, []).append((y, rt))
        elif is_h(y):
          h_on.setdefault(x, []).append(y)
        elif not is_h(x) and elements[y] not in METAL_ELEMENTS:
          heavy.setdefault(x, []).append(y)
    for p in simple:
      add(p.i_seqs[0], p.i_seqs[1], p.origin_id, sgtbx.rt_mx())
    for p in asu:
      add(p.i_seq, p.j_seq, p.origin_id, asu_mappings.get_rt_mx_ji(pair=p))
    candidates = self._candidates(atoms, elements, heavy, h_on)
    if not candidates: return
    parents = sorted(set(i_x for i_x, i_h, judged in candidates))
    # contacts: proxies, then geometric; one per metal image
    working_cs = model.crystal_symmetry()
    by_x = {}
    def image_of(i_m, rt, cs):
      if rt.is_unit_mx(): return matrix.col(atoms[i_m].xyz)
      uc = cs.unit_cell()
      return matrix.col(uc.orthogonalize(rt * uc.fractionalize(atoms[i_m].xyz)))
    def record(i_x, i_m, rt, cs, evidence):
      image = image_of(i_m, rt, cs)
      for c in by_x.setdefault(i_x, []):
        if c.i_metal == i_m and (c.image - image).length() < 0.01: return
      c = group_args(i_x = i_x, i_metal = i_m, sym_op = str(rt), image = image,
        distance = (image - matrix.col(atoms[i_x].xyz)).length(),
        evidence = evidence)
      by_x[i_x].append(c)
      self.contacts.append(c)
    for i_x in parents:
      for i_m, rt in proxied.get(i_x, []):
        if elements[i_m] not in self.metals: continue
        record(i_x, i_m, rt, working_cs, 'metal coordination proxy')
    metals = [i for i, e in enumerate(elements) if e in self.metals
              and _covalent_radius(e) is not None]
    if metals:
      subset = parents + metals
      max_cut = max(_covalent_radius(elements[i]) for i in parents) + \
        max(_covalent_radius(elements[i]) for i in metals) + \
        THIOL_METAL_WINDOW[1]
      xyz = flex.vec3_double([atoms[i].xyz for i in subset])
      cs = crystal_symmetry
      if cs is None:
        lo, hi = xyz.min(), xyz.max()
        size = max(hi[k] - lo[k] for k in range(3)) + 3 * max_cut + 10
        cs = crystal.symmetry(unit_cell=(size, size, size, 90, 90, 90),
          space_group_symbol='P1')
        xyz = xyz - lo + (max_cut, max_cut, max_cut)
      uc = cs.unit_cell()
      frac = uc.fractionalize(xyz)
      am = cs.special_position_settings().asu_mappings(
        buffer_thickness=max_cut)
      am.process_sites_frac(original_sites=frac)
      pt = crystal.pair_asu_table(asu_mappings=am)
      pt.add_all_pairs(distance_cutoff=max_cut)
      n_parents = len(parents)
      for i, js in enumerate(pt.table()):
        for j, j_syms in js.items():
          for group in j_syms:
            for j_sym in group:
              rt_ji = am.get_rt_mx(i_seq=i, i_sym=0).inverse().multiply(
                am.get_rt_mx(i_seq=j, i_sym=j_sym))
              for x, y, rt in ((i, j, rt_ji), (j, i, rt_ji.inverse())):
                if x >= n_parents or y < n_parents: continue
                i_x, i_m = subset[x], subset[y]
                if model_index[i_x] != model_index[i_m]: continue
                d = uc.distance(frac[x], rt * frac[y])
                r = _covalent_radius(elements[i_x]) + \
                  _covalent_radius(elements[i_m])
                if not (r + THIOL_METAL_WINDOW[0] <= d <= r + THIOL_METAL_WINDOW[1]):
                  continue
                if crystal_symmetry is None: rt = sgtbx.rt_mx()
                record(i_x, i_m, rt, crystal_symmetry, 'geometric')
    # verdicts
    for i_x, i_h, judged in candidates:
      h_atom, x_atom = atoms[i_h], atoms[i_x]
      alt = h_atom.parent().altloc
      population, why = None, ''
      if alt == '' or x_atom.parent().altloc == alt:
        population = x_atom.occ
      else:
        occ = [a.occ for ag in x_atom.parent().parent().atom_groups()
               if ag.altloc == alt for a in ag.atoms() if not is_h(a.i_seq)]
        if not occ:
          why = 'altloc group without heavy atoms'
        elif max(occ) - min(occ) > 1.e-6:
          why = 'heavy atoms of conformer %s differ in occupancy' % alt
        else:
          population = occ[0]
      sites, counted, reasons = [], [], []
      for c in by_x.get(i_x, []):
        m = atoms[c.i_metal]
        if not _compatible(m, h_atom):
          reasons.append('%s: altloc %s, H in %s' % (_atom_label(m),
            m.parent().altloc, alt))
          continue
        if m.occ <= 0:
          reasons.append('%s: occupancy 0' % _atom_label(m))
          continue
        counted.append(c)
      if not judged:
        if counted:
          self.reported.append(group_args(i_h = i_h, i_x = i_x,
            contacts = counted))
        continue
      # metal sites: connected components of the images of one metal atom,
      # an edge between images closer than METAL_SITE_IMAGE_DISTANCE
      # (union-find, so the order of the contacts does not matter)
      key = lambda c: (atoms[c.i_metal].parent().parent().memory_id(),
                       atoms[c.i_metal].name)
      root = list(range(len(counted)))
      def find(k):
        while root[k] != k:
          root[k] = root[root[k]]
          k = root[k]
        return k
      for k, c in enumerate(counted):
        for l in range(k):
          if (key(c) == key(counted[l]) and (c.image - counted[l].image).length()
              < METAL_SITE_IMAGE_DISTANCE):
            root[find(k)] = find(l)
      components = {}
      for k, c in enumerate(counted):
        components.setdefault(find(k), []).append(c)
      for k in sorted(components):
        site = group_args(key = key(components[k][0]),
          contacts = components[k], occ = {})
        for c in components[k]:
          site.occ[atoms[c.i_metal].parent().altloc] = atoms[c.i_metal].occ
        sites.append(site)
      proxy = any(c.evidence != 'geometric' for c in counted)
      i_metal = None
      if not counted:
        verdict, coverage, reason = 'not counted', 0., '; '.join(reasons)
      else:
        coverage, i_metal = max((sum(site.occ.values()),
          min(c.i_metal for c in site.contacts)) for site in sites)
        if population is None:
          verdict, reason = 'unresolved', why
        elif coverage >= population - 1.e-6:
          verdict, reason = self.covered, ''
        else:
          verdict, reason = 'unresolved', 'coverage below population'
      self.h_sites.append(group_args(i_h = i_h, i_x = i_x, verdict = verdict,
        reason = reason, population = population, coverage = coverage,
        i_metal = i_metal, contacts = counted, proxy = proxy, sites = sites))

  def contact_text(self, c):
    a = self.atoms
    return '%s - %s %.2f A%s (%s)' % (_atom_label(a[c.i_x]),
      _atom_label(a[c.i_metal]), c.distance,
      '' if c.sym_op == 'x,y,z' else ' sym ' + c.sym_op, c.evidence)

  def unresolved_message(self, site):
    a = self.atoms
    h, x = a[site.i_h], a[site.i_x]
    if site.population is None:
      return ('%s: the population of %s is unknown (%s). %s was added using '
        'the default %s state; review this site.' % (_atom_label(h),
        h.name.strip(), site.reason, h.name.strip(), self.state))
    alt = h.parent().altloc
    what = x.name.strip() if alt == '' else '%s conformer %s' % (
      x.name.strip(), alt)
    metal = a[site.i_metal]
    return ('%s: %s occupancy %.2f exceeds coordinating %s occupancy %.2f. '
      'Protonation may differ between bound and unbound populations. %s was '
      'added using the default %s state; review this site.' % (
      _atom_label(h), what, site.population, metal.name.strip(),
      site.coverage, h.name.strip(), self.state))

class thiol_metal_sites(metal_h_sites):
  '''
  Thiols bound to metals (A4). Thiol: an S with exactly one non-metal heavy
  neighbour (counted by residue and atom name, so altloc copies count once)
  in the bond restraints of the H site's conformer, carrying that H. Covered:
  thiolate.
  '''
  covered = 'thiolate'
  state = 'thiol'

  def _candidates(self, atoms, elements, heavy, h_on):
    candidates = []
    for i_s, hs in h_on.items():
      if elements[i_s] != 'S': continue
      for i_h in hs:
        nb = set((atoms[j].parent().parent().memory_id(), atoms[j].name)
                 for j in heavy.get(i_s, [])
                 if _compatible(atoms[j], atoms[i_h]))
        if len(nb) == 1: candidates.append((i_s, i_h, True))
    return candidates

def residue_is_amino_acid(resname, mon_lib_srv):
  '''
  The residue-name classes (common or modified amino acid); for other names
  the monomer library's classification from the backbone connectivity
  (is_peptide()), also for user restraints: 9FZ and 9G2 (group "ligand")
  count, the NH2 cap does not.
  '''
  resname = resname.strip()
  if get_class(name=resname) in ('common_amino_acid', 'modified_amino_acid'):
    return True
  if mon_lib_srv is None: return False
  ccid = mon_lib_srv.get_comp_comp_id_direct(resname)
  return ccid is not None and ccid.is_peptide()

class metal_nh_sites(metal_h_sites):
  '''
  N-H bound to metals (A5): backbone amide N with a preceding residue and
  thymine/uracil N3 are judged (covered: deprotonated); other N with H, not
  His ring N, are reported only. Metals: NH_METAL_ELEMENTS. Each H is
  classified with the heavy neighbours of its own conformer.
  '''
  metals = NH_METAL_ELEMENTS
  covered = 'deprotonated'
  state = 'N-H'

  def __init__(self, model, crystal_symmetry=None):
    self._mon_lib_srv = model.get_mon_lib_srv()
    self._amino_acid = {}
    metal_h_sites.__init__(self, model, crystal_symmetry)

  def is_amino_acid(self, resname):
    # residue_is_amino_acid(), cached per name
    resname = resname.strip()
    if resname not in self._amino_acid:
      self._amino_acid[resname] = residue_is_amino_acid(resname,
        self._mon_lib_srv)
    return self._amino_acid[resname]

  def n_class(self, atoms, heavy, i_n, i_h):
    '''
    The class of the N-H site (N i_n, H i_h): 'amide' (a peptide bond: N and
    the other residue amino acids, N bonded to the CA of its own residue and
    to the atom C of the other; no observed O needed), 'N3' (thymine/uracil
    N3), 'his' (His ring N, left to the Optimizer) or 'other'. Neighbours:
    only those compatible with the H's conformer.
    '''
    n, h = atoms[i_n], atoms[i_h]
    ag = n.parent()
    name, resname = n.name.strip(), ag.resname.strip()
    if resname == 'HIS' and name in ('ND1', 'NE2'):
      return 'his'
    if name == 'N3' and resname in THYMINE_URACIL:
      return 'N3'
    if name != 'N' or not self.is_amino_acid(resname):
      return 'other'
    own = ag.parent().memory_id()
    nb = [atoms[j] for j in heavy.get(i_n, []) if _compatible(atoms[j], h)]
    if not any(a.name.strip() == 'CA' and a.parent().parent().memory_id() == own
               for a in nb):
      return 'other'
    for c in nb:
      if (c.name.strip() == 'C' and c.parent().parent().memory_id() != own
          and self.is_amino_acid(c.parent().resname)):
        return 'amide'
    return 'other'

  def _candidates(self, atoms, elements, heavy, h_on):
    candidates = []
    for i_n, hs in h_on.items():
      if elements[i_n] != 'N': continue
      for i_h in hs:
        kind = self.n_class(atoms, heavy, i_n, i_h)
        if kind == 'his': continue
        candidates.append((i_n, i_h, kind in ('amide', 'N3')))
    return candidates

PRIMARY_AMINE = 'primary amine'

def _unsaturated_atoms(resname, mon_lib_srv, cache):
  '''
  Names of the atoms in a double, triple, aromatic or delocalized bond of the
  restraint dictionary, and of the CCD where it describes the residue
  (geostd amino acids give every bond as "coval").
  '''
  resname = resname.strip()
  if resname in cache: return cache[resname]
  names = set()
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  if cc is not None:
    for b in cc.bond_list:
      if str(b.type).strip().lower() in ('double', 'triple', 'aromatic',
                                         'deloc'):
        names.update((b.atom_id_1.strip(), b.atom_id_2.strip()))
  if _ccd_describes(cc, resname):
    from mmtbx.chemical_components import get_cif_dictionary
    try: cc_cif = get_cif_dictionary(resname)
    except Exception: cc_cif = None
    for b in (cc_cif or {}).get('_chem_comp_bond', []):
      if (str(b.value_order).strip().upper() in ('DOUB', 'TRIP', 'AROM') or
          str(getattr(b, 'pdbx_aromatic_flag', 'N')).strip().upper() == 'Y'):
        names.update((b.atom_id_1.strip(), b.atom_id_2.strip()))
  cache[resname] = names
  return names

# a placed NH2 is accepted within these of its targets (N-H in A, angles in
# deg; 6 deg: C-N-M about 95-128 deg). h_h: distinct H; the H-N-H tolerance
# sets their spacing (tetrahedral: 1.45 A at the X-ray N-H of 0.89 A)
NH2_TOLERANCE = group_args(d = 0.02, angle = 6., h_h = 1.0)

def place_metal_nh2(r_n, r_c, r_m, targets, fixed=None):
  """
  The two H of a metal-bound NH2, constructed from the frame C, N and the
  metal (r_m: its Cartesian image, symmetry already applied), and accepted
  only within NH2_TOLERANCE of each H's own targets (targets.d, targets.cnh:
  N-H and C-N-H per H; targets.hnh: H-N-H, also the M-N-H target, the metal
  in the place of the third H). A construction, not an exact solution of
  every dictionary angle: each H lies at its N-H distance, half the H-N-H
  target out of the C-N-M plane, on the in-plane direction that riding H's
  2tetra would give for its C-N-H and M-N-H targets. The side is fixed, so
  the names are: the first H on the side of (C-N) x (M-N), the second on
  the other (convert_metal_bound_amines(): the first kept name in
  dictionary order, e.g. HZ2, H2, takes the first side).
  fixed: the coordinates of a first H that stays (one supplied H): the
  second H then takes the remaining tetrahedral direction away from C, the
  metal and that H.
  Returns ([xyz, xyz], None) in the order of the targets, or (None, reason).
  """
  cnm = None
  try:
    u10, u20 = (r_c - r_n).normalize(), (r_m - r_n).normalize()
    cnm = math.degrees(math.acos(max(-1., min(1., u10.dot(u20)))))
    v0 = u10.cross(u20)
    if v0.length() < 1.e-3:
      return None, 'degenerate frame, C-N-M %.0f deg' % cnm
    v0 = v0.normalize()
    if fixed is not None:
      u_f = (fixed - r_n).normalize()
      w = -(u10 + u20 + u_f)
      if w.length() < 1.e-3:
        return None, 'degenerate frame with the supplied H'
      hs = [fixed, r_n + w.normalize() * targets.d[1]]
    else:
      c0 = u10.dot(u20)
      half = math.radians(targets.hnh) * 0.5
      c2 = math.cos(math.radians(targets.hnh))
      hs = []
      for sign, d, cnh in zip((1, -1), targets.d, targets.cnh):
        c1 = math.cos(math.radians(cnh))
        a = (c1 - c0*c2) / (1 - c0*c0)
        b = (c2 - c0*c1) / (1 - c0*c0)
        d0 = a * u10 + b * u20
        if d0.length() < 1.e-3:
          return None, 'geometry not placeable, C-N-M %.0f deg' % cnm
        d0 = d0.normalize()
        hs.append(r_n + d * (math.cos(half) * d0 + sign * math.sin(half) * v0))
  except (ArithmeticError, ValueError):
    return None, 'geometry not placeable%s' % (
      '' if cnm is None else ', C-N-M %.0f deg' % cnm)
  problem = nh2_placement_problem(r_n, r_c, r_m, hs, targets,
    check_d = (fixed is None, True))
  if problem is not None: return None, problem
  return hs, None

def nh2_placement_problem(r_n, r_c, r_m, hs, targets, check_d=(True, True)):
  """
  None if the two H sites hs make an NH2 within NH2_TOLERANCE of targets,
  else the reason. check_d: whether each H's N-H is checked (not for a
  supplied H).
  """
  cnm = (r_c - r_n).angle(r_m - r_n, deg=True)
  if len(hs) != 2 or not all(math.isfinite(x) for h in hs for x in h):
    return 'geometry not placeable, C-N-M %.0f deg' % cnm
  if (hs[0] - hs[1]).length() < NH2_TOLERANCE.h_h:
    return 'geometry not placeable, C-N-M %.0f deg (H...H %.2f A)' % (cnm,
      (hs[0] - hs[1]).length())
  out = []
  for h, d, cnh, with_d in zip(hs, targets.d, targets.cnh, check_d):
    checks = [('C-N-H', (r_c - r_n).angle(h - r_n, deg=True), cnh,
               NH2_TOLERANCE.angle),
              ('M-N-H', (r_m - r_n).angle(h - r_n, deg=True), targets.hnh,
               NH2_TOLERANCE.angle)]
    if with_d:
      checks.insert(0, ('N-H', (h - r_n).length(), d, NH2_TOLERANCE.d))
    for what, value, target, tol in checks:
      if abs(value - target) > tol:
        out.append('%s %.1f (target %.1f)' % (what, value, target))
  hnh = (hs[0] - r_n).angle(hs[1] - r_n, deg=True)
  if abs(hnh - targets.hnh) > NH2_TOLERANCE.angle:
    out.append('H-N-H %.1f (target %.1f)' % (hnh, targets.hnh))
  if out:
    return 'geometry not placeable, C-N-M %.0f deg: %s' % (cnm, ', '.join(out))
  return None

class metal_amine_sites(metal_h_sites):
  '''
  Primary amines bound to metals (A5 C1b). Primary amine: an N with exactly
  one non-metal heavy neighbour, a C, in the H's conformer, in no planarity
  restraint, with neither the N nor that neighbour in a double, triple,
  aromatic or delocalized bond (amide, amidine, guanidinium, aromatic amine,
  enamine). Each H on it is judged (covered: NH2). amine[i_h]: the branch
  that decided (PRIMARY_AMINE or the reason it is not one) and the heavy
  neighbour i_c (None unless it is one atom). Metals: NH_METAL_ELEMENTS.
  '''
  metals = NH_METAL_ELEMENTS
  covered = 'NH2'
  state = 'NH3+'

  def __init__(self, model, crystal_symmetry=None, mon_lib_srv=None):
    self._mon_lib_srv = mon_lib_srv or model.get_mon_lib_srv()
    self._unsaturated = {}
    grm = model.get_restraints_manager().geometry
    self._planar = set(i for p in (grm.planarity_proxies or [])
                       for i in p.i_seqs)
    self.amine = {}
    metal_h_sites.__init__(self, model, crystal_symmetry)

  def amine_class(self, atoms, heavy, i_n, i_h):
    n, h = atoms[i_n], atoms[i_h]
    nb = sorted(set(j for j in heavy.get(i_n, []) if _compatible(atoms[j], h)))
    names = set((atoms[j].parent().parent().memory_id(), atoms[j].name)
                for j in nb)
    if len(names) != 1:
      return group_args(branch='%d heavy neighbours' % len(names), i_c=None)
    if i_n in self._planar:
      return group_args(branch='in a planarity restraint', i_c=None)
    c = atoms[nb[0]]
    if _atom_element(c) != 'C':
      return group_args(branch='heavy neighbour %s is not C' % c.name.strip(),
        i_c=None)
    unsaturated = lambda a: a.name.strip() in _unsaturated_atoms(
      a.parent().resname, self._mon_lib_srv, self._unsaturated)
    if unsaturated(n):
      branch = 'N in a multiple or aromatic bond'
    elif unsaturated(c):
      branch = '%s in a multiple or aromatic bond' % c.name.strip()
    else:
      branch = PRIMARY_AMINE
    return group_args(branch=branch, i_c=nb[0] if len(nb) == 1 else None)

  def _candidates(self, atoms, elements, heavy, h_on):
    candidates = []
    for i_n, hs in h_on.items():
      if elements[i_n] != 'N': continue
      for i_h in sorted(set(hs)):
        self.amine[i_h] = self.amine_class(atoms, heavy, i_n, i_h)
        if self.amine[i_h].branch == PRIMARY_AMINE:
          candidates.append((i_n, i_h, True))
    return candidates

def _h_name(h):
  # name with its altloc, as reported_metal_h()
  altloc = h.parent().altloc.strip()
  return h.name.strip() + (' altloc %s' % altloc if altloc else '')

def conformer_groups(items, altloc):
  '''
  The H of one N per complete conformer, as _compatible() pairs them:
  [(altloc, items, mixed)], for each altloc of the H the H of that altloc
  plus the blank ones; one blank group if no H has an altloc. mixed: the
  group has blank and alternate H.
  '''
  alts = sorted(set(altloc(x) for x in items) - set(['']))
  if not alts: return [('', list(items), False)]
  result = []
  for alt in alts:
    group = [x for x in items if altloc(x) in ('', alt)]
    result.append((alt, group, any(altloc(x) == '' for x in group)))
  return result

def dictionary_nh3_plus(model, n, i_n, mon_lib_srv):
  '''
  Whether the dictionary makes N n (i_seq i_n) NH3+: energy type NT3 (also
  from mod_NH3) or formal charge +1.
  '''
  types = getattr(model, '_type_energies', None)
  if types is not None and i_n < len(types) and types[i_n].strip() == 'NT3':
    return True
  try:
    cc = mon_lib_srv.get_comp_comp_id_direct(n.parent().resname.strip())
  except Exception:
    cc = None
  if cc is None: return False
  a = cc.atom_dict().get(n.name.strip())
  return a is not None and a.formal_charge() == 1

def protect_supplied_metal_amines(model, crystal_symmetry=None,
                                  use_neutron_distances=None):
  '''
  approach=optimize with metal_bound_amines_as_nh2 (bypasses placement): a
  primary amine that its dictionary makes NH3+ but that carries two H, with
  a metal of NH_METAL_ELEMENTS in the window (any occupancy, alternates,
  symmetry; as convert_metal_bound_amines()), keeps both H where they are:
  they are taken out of the riding H manager before its first idealization
  (all other H keep their parameterization, as the Optimizer would build
  it). No H is added or deleted. H are counted per complete conformer
  (conformer_groups()). Limit: a blank H shared with an NH3 conformer also
  belongs to that conformer and may be moved by its rotator; its input
  position is not guaranteed. The NH2 conformer's own alternate H keep
  their input positions; such a site is listed as partly preserved.
  :return: ([(label, contacts)] of the NH2 kept, [(label, contacts)] of
    those partly preserved)
  '''
  elements = model.get_hierarchy().atoms().extract_element()
  if not any(e.strip().upper() in NH_METAL_ELEMENTS for e in elements):
    return [], []
  sites = metal_amine_sites(model, crystal_symmetry)
  in_window = {}
  for c in sites.contacts:
    in_window.setdefault(c.i_x, []).append(c)
  atoms = sites.atoms
  by_n = {}
  for site in sites.h_sites:
    by_n.setdefault(site.i_x, []).append(site.i_h)
  mon_lib_srv = model.get_mon_lib_srv()
  listed, partial, protected = [], [], []
  for i_n in sorted(by_n):
    if i_n not in in_window: continue
    n = atoms[i_n]
    if not dictionary_nh3_plus(model, n, i_n, mon_lib_srv): continue
    # two H in a complete conformer (blank H counted in each)
    conformers = conformer_groups(by_n[i_n], lambda i: atoms[i].parent().altloc)
    groups = [(alt, hs) for alt, hs, mixed in conformers if len(hs) == 2]
    if not groups: continue
    kept = set(alt for alt, hs in groups)
    others = [alt for alt, hs, mixed in conformers if alt not in kept]
    contact = '; '.join(sites.contact_text(c) for c in in_window[i_n])
    for alt, hs in groups:
      protected.extend(hs)
      shared = [atoms[i].name.strip() for i in hs
                if atoms[i].parent().altloc == '']
      if shared and others:
        own = [atoms[i].name.strip() for i in hs
               if atoms[i].parent().altloc != '']
        partial.append(('%s conformer %s: supplied NH2 partly preserved (%s '
          'kept; shared %s belongs to conformer %s\'s NH3 and may move)' % (
          _atom_label(n), alt, ', '.join(own), ', '.join(shared),
          ', '.join(others)), contact))
      else:
        listed.append(('%s (%s)' % (_atom_label(n), ', '.join(
          _h_name(atoms[i]) for i in hs)), contact))
  if not protected: return listed, partial
  if model.riding_h_manager is None:
    model.setup_riding_h_manager(idealize = False)
    idealize = True
  else:
    idealize = False
  rhm = model.riding_h_manager
  if rhm is None: return listed, partial
  for i in protected:
    rhm.h_parameterization[i] = None
  rhm.parameterization_cpp = rhm.get_parameterization_cpp(
    h_parameterization = rhm.h_parameterization)
  if idealize:
    model.idealize_h_riding(use_neutron_distances = use_neutron_distances)
  return listed, partial

def supplied_metal_h(sites, is_supplied):
  '''
  [(H label, contacts, connectivity)] of the supplied H (is_supplied(atom))
  on judged sites with a counted metal; connectivity: whether the parent
  has a metal-coordination proxy.
  '''
  result = []
  for site in sites.h_sites:
    if site.verdict == 'not counted': continue
    h = sites.atoms[site.i_h]
    if not is_supplied(h): continue
    result.append((_atom_label(h),
      '; '.join(sites.contact_text(c) for c in site.contacts),
      'metal-coordination proxy' if site.proxy else
      'no metal-coordination proxy'))
  return result

def reported_metal_h(sites):
  '''
  [(label, contacts)] of the candidates that are reported only, one entry
  per parent atom: its label with the names of its H (and their altloc).
  '''
  result, by_x = [], {}
  for r in sites.reported:
    if r.i_x not in by_x:
      by_x[r.i_x] = [r, []]
      result.append(by_x[r.i_x])
    h = sites.atoms[r.i_h]
    altloc = h.parent().altloc.strip()
    by_x[r.i_x][1].append(h.name.strip() +
      (' altloc %s' % altloc if altloc else ''))
  return [('%s (%s)' % (_atom_label(sites.atoms[r.i_x]), ', '.join(names)),
    '; '.join(sites.contact_text(c) for c in r.contacts))
    for r, names in result]

def reported_metal_h_rows(sites, key):
  '''
  reported_metal_h() before it is joined, for a list made at the end of
  placement (its atoms may be gone by then): [(parent key, parent label,
  [(H key, H name and altloc)], contacts)], keys from key(atom).
  '''
  rows, by_x = [], {}
  for r in sites.reported:
    x, h = sites.atoms[r.i_x], sites.atoms[r.i_h]
    if r.i_x not in by_x:
      by_x[r.i_x] = (key(x), _atom_label(x), [],
        '; '.join(sites.contact_text(c) for c in r.contacts))
      rows.append(by_x[r.i_x])
    altloc = h.parent().altloc.strip()
    by_x[r.i_x][2].append((key(h), h.name.strip() +
      (' altloc %s' % altloc if altloc else '')))
  return rows

def _metal_h_text(titles, not_placed, unresolved, supplied, hint,
                  reported=()):
  lines = []
  if not_placed:
    lines.append(titles.not_placed)
    for label, contact in not_placed:
      lines.append('  %s  (%s)' % (label, contact))
  if unresolved:
    lines.append(titles.unresolved)
    for msg in unresolved:
      lines.append('  ' + msg)
  if supplied:
    lines.append(titles.supplied)
    for label, contact, connectivity in supplied:
      lines.append('  %s  (%s); %s' % (label, contact, connectivity))
    if hint:
      lines.append('  to drop it: ' + hint)
  if reported:
    lines.append(titles.reported)
    for label, contact in reported:
      lines.append('  %s  (%s)' % (label, contact))
  return '\n'.join(lines)

def metal_thiol_text(not_placed=(), unresolved=(), supplied=(), hint=None):
  '''
  Report of thiols bound to metals: H not placed [(label, contacts)],
  messages of unresolved sites, supplied H kept (supplied_metal_h()), and
  how to drop them (hint); '' if there is nothing to report.
  '''
  return _metal_h_text(group_args(
    not_placed = 'Thiol H not placed: the S is bound to a metal, taken as '
      'thiolate',
    unresolved = 'Thiol H placed at metal sites with unresolved protonation',
    supplied = 'Supplied thiol H at metal sites: %s' % THIOL_H_SUPPLIED_NOTE),
    not_placed, unresolved, supplied, hint)

def metal_nh_text(not_placed=(), unresolved=(), supplied=(), reported=(),
                  hint=None):
  '''
  Report of N-H at metals (metal_nh_sites): H not placed, unresolved sites,
  supplied H kept, and the other N-H next to a metal, left unchanged
  [(label, contacts)]; '' if there is nothing to report.
  '''
  return _metal_h_text(group_args(
    not_placed = 'N-H not placed: backbone amide N or thymine/uracil N3 bound '
      'to a metal, taken as deprotonated',
    unresolved = 'N-H placed at metal sites with unresolved protonation',
    supplied = 'Supplied N-H at metal sites: %s' % METAL_NH_SUPPLIED_NOTE,
    reported = 'N-H next to a metal, not changed'),
    not_placed, unresolved, supplied, hint, reported)

def metal_amine_text(converted=(), unresolved=(), supplied=(),
                     chain_starts=(), supplied_nh2=(), supplied_nh2_partly=()):
  '''
  Report of metal_bound_amines_as_nh2: NH3+ placed as NH2 [(label, contacts)],
  messages of sites left as NH3+, supplied NH3+ kept [(label, contacts)],
  amino-acid chain starts on a metal not built as NH3 [(label, contacts)];
  '' if there is nothing to report.
  '''
  lines = []
  for title, entries in (
      ('Coordinating NH3+ placed as NH2 (the metal takes the place of the H '
       'not placed; the two H are not parameterized for riding)', converted),
      ('Coordinating NH3+ left as NH3+ (unresolved)', unresolved),
      ('Supplied NH3+ at metal sites: kept as supplied; coordinating primary '
       'amines are usually NH2', supplied),
      ('Supplied NH2 at metal sites: kept as supplied, not completed to NH3',
       supplied_nh2),
      ('Supplied NH2 at metal sites, partly preserved (not completed to NH3)',
       supplied_nh2_partly),
      ('Metal-bound chain start, terminus not confirmed (H as '
       'n_terminal_charge places them)', chain_starts)):
    if not entries: continue
    lines.append(title)
    for e in entries:
      lines.append('  ' + (e if isinstance(e, str) else '%s  (%s)' % e))
  return '\n'.join(lines)

# ==============================================================================

def get_reduce_pdb_interpretation_params(use_neutron_distances):
  '''
  Create pdb_interpretation parameter scope.
  Do this in a function so other programs (reduce2) can use the same parameters
  '''
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.restraints_library.cdl=False # XXX this triggers a bug !=360
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
  p.pdb_interpretation.use_neutron_distances = use_neutron_distances
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  #p.pdb_interpretation.automatic_linking.link_metals = True
  p.pdb_interpretation.automatic_linking.link_residues = True
  p.pdb_interpretation.automatic_linking.exclude_hydrogens_from_bonding_decisions = True
  #
  return p

# ==============================================================================

def h1_h2_from_A_X_d_angles(A, X, d, ang1_deg, ang2_deg):
    """
    AI generated code.

    In 3D. Find coordinates of point H1 and H2, if we know both are distanced by
    d from point X, we know the angle (H1,X,A) and (H2,X,A), and we know that H1
    and H2 do not overlap, and all H1,H2,A,X belong to the plane. Coordinates of
    points A,X are known. Need a function in python. Angles are in degrees. Only
    inputs are: coordinates of A and X, distance d and two angles H1XA and H2XA.
    Use pure python. If solution is not unique, just pick any one.

    Pure-python 3D construction of one valid pair (H1, H2), choosing an
    arbitrary plane.

    Inputs:
      A, X: (x,y,z) tuples/lists
      d: distance from X to each Hi
      ang1_deg: angle(H1, X, A) in degrees
      ang2_deg: angle(H2, X, A) in degrees

    Returns:
      (H1, H2) as two (x,y,z) tuples

    Raises:
      ValueError on degenerate/unsatisfiable cases (e.g. A==X, d<0, forced
      overlap).
    """
    # -------- basic vector ops --------
    def v_sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def v_add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def v_mul(s, v): return (s*v[0], s*v[1], s*v[2])
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def norm(v): return math.sqrt(dot(v, v))
    def unit(v, eps=1e-12):
      n = norm(v)
      if n < eps:
        raise ValueError("Cannot normalize near-zero vector.")
      return (v[0]/n, v[1]/n, v[2]/n)
    def proj_perp(ref, u):
      # component of ref perpendicular to u: ref - (ref·u)u
      return v_sub(ref, v_mul(dot(ref, u), u))
    def dist(a, b):
      return norm(v_sub(a, b))
    # -------- input checks --------
    if d < 0:
      raise ValueError("d must be non-negative.")
    AX = v_sub(A, X)
    if norm(AX) < 1e-12:
      raise ValueError("A and X coincide; angle(H, X, A) is undefined.")
    # u = direction from X to A
    u = unit(AX)
    # Choose an arbitrary but deterministic plane through line XA
    # by picking a reference axis not parallel to u.
    ref = (0.0, 0.0, 1.0)
    if abs(dot(u, ref)) > 0.99:
      ref = (0.0, 1.0, 0.0)
    v_raw = proj_perp(ref, u)
    if norm(v_raw) < 1e-12:
      # extremely rare fallback
      ref = (1.0, 0.0, 0.0)
      v_raw = proj_perp(ref, u)
    v = unit(v_raw)
    def point_for(theta_deg, sign):
      th = math.radians(theta_deg)
      c = math.cos(th)
      s = math.sin(th)
      direction = v_add(v_mul(c, u), v_mul(sign * s, v))
      return v_add(X, v_mul(d, direction))
    # Try sign combinations until H1 and H2 don't overlap
    candidates = [
        (+1.0, -1.0),
        (-1.0, +1.0),
        (+1.0, +1.0),
        (-1.0, -1.0),
    ]
    for s1, s2 in candidates:
      H1 = point_for(ang1_deg, s1)
      H2 = point_for(ang2_deg, s2)
      if dist(H1, H2) > 1e-9:
        return H1, H2
    # If we got here, overlap is unavoidable for these inputs
    # (e.g. d=0, angles 0/180, etc.)
    raise ValueError(
    "Cannot produce distinct H1 and H2 for these inputs (overlap unavoidable).")

def find_H1_H2(X, d, angle_deg):
  """
  AI generated code.

  Prompt: In 3D. Find coordinates of point H1 and H2, if we know both are
  distanced by
  d from point X, we know the angle (H1,X,H2). Coordinates of
  point X are known. Need a function in python. Angle is in degrees. Only
  inputs are: coordinates X, distance d and angle H1XH2.
  Use pure python. If solution is not unique, just pick any one.
  Return one pair (H1, H2) in 3D such that:
    |H1 - X| = d
    |H2 - X| = d
    angle(H1, X, H2) = angle_deg  (in degrees)

  Picks a convenient solution in the XY-plane.
  """
  x, y, z = map(float, X)
  d = float(d)
  if d < 0: raise ValueError("d must be non-negative")
  # If d == 0, both points must coincide with X
  if d == 0: return (x, y, z), (x, y, z)
  theta = math.radians(angle_deg)
  # Choose direction for H1 along +x axis from X
  H1 = (x + d, y, z)
  # Place H2 at the requested angle from H1 around X in the XY-plane
  H2 = (x + d * math.cos(theta),
        y + d * math.sin(theta),
        z)
  return H1, H2

def _water_oh_targets(model, use_neutron_distances):
  '''O-H length and H-O-H angle from the HOH restraints; constants if absent'''
  d = 0.98 if use_neutron_distances else 0.85
  angle = 103.91
  comp = model.get_mon_lib_srv().get_comp_comp_id_direct("HOH")
  if comp is not None:
    for bond in comp.bond_list:
      value = bond.value_dist_neutron if use_neutron_distances \
        else bond.value_dist
      if value is not None:
        d = value
        break
    for a in comp.angle_list:
      if a.value_angle is not None:
        angle = a.value_angle
        break
  return d, angle

def workaround_002(model, selection, use_neutron_distances=False):
  d, angle = _water_oh_targets(model, use_neutron_distances)
  h = model.get_hierarchy()
  atoms = h.atoms()
  for m in h.models():
    for c in m.chains():
      for con in c.conformers():
        for r in con.residues():
          if not get_class(name=r.resname) == "common_water": continue
          three = r.atoms().extract_i_seq()
          #assert three in selection
          # only this residue's own O: exactly one O and two H/D
          o = [a for a in r.atoms() if a.element.strip().upper() == "O"]
          ij = [a.i_seq for a in r.atoms()
                if a.element.strip().upper() in ("H", "D")]
          if len(o) != 1 or len(ij) != 2: continue
          p1, p2 = find_H1_H2(X=o[0].xyz, d=d, angle_deg=angle)
          atoms[ij[0]].xyz = p1
          atoms[ij[1]].xyz = p2
  # direct xyz edits: keep the X-ray structure in step
  model.set_sites_cart(model.get_sites_cart())

def _chiral_volume(c, a, b, h):
  c, a, b, h = [matrix.col(x) for x in (c, a, b, h)]
  return (a - c).dot((b - c).cross(h - c))

# CH2 groups whose two H the CCD's ideal coordinates name with the opposite
# hand to the same entry's model coordinates. The model hand is the one reduce
# places and the one deposited models carry, so these read the model set.
# Verified group by group against reduce (1bzs, 2h34). Only where an authority
# exists: the CCD sets disagree in thousands of ligand groups, and nucleotides
# are the other way round (for C5' the ideal set is the right one), so this is
# a list of known-bad entries, not a preference for either set.
_CCD_MODEL_CH2 = {
  'ARG': ('CB', 'CG'),
  'ILE': ('CG1',),
  'LEU': ('CB',),
  'MET': ('CB', 'CG'),
  'MSE': ('CB', 'CG'),
  }

def _h_groups(elements, bond_pairs, sites, n_h, n_heavy):
  '''
  Heavy atoms with exactly n_h H and n_heavy heavy neighbours, all sites known.
  (2, 2) are the CH2 centres, (3, 1) the propellers (CH3, NH3).
  '''
  neighbors = {}
  for i, j in bond_pairs:
    neighbors.setdefault(i, []).append(j)
    neighbors.setdefault(j, []).append(i)
  result = []
  for p, ns in neighbors.items():
    if elements.get(p) in ('H', 'D'): continue
    hs = tuple(sorted(n for n in ns if elements.get(n) in ('H', 'D')))
    hv = tuple(sorted(n for n in ns if elements.get(n) not in ('H', 'D')))
    if len(hs) == n_h and len(hv) == n_heavy and all(n in sites for n in (p,)+hv+hs):
      result.append((p, hv, hs, sites))
  return result

def _dictionary_sites(source_info, resname):
  '''
  Ideal sites of resname from its restraint file, if it has coordinates
  (geostd ligands). One file can define several ligands with the same names.
  '''
  if not source_info or not source_info.startswith("file:"): return None
  import iotbx.cif
  try: cif_model = iotbx.cif.reader(file_path=source_info[5:].strip()).model()
  except Exception: return None
  for block_name, block in cif_model.items():
    if not ("_chem_comp_atom.atom_id" in block and "_chem_comp_atom.x" in block):
      continue
    names = block["_chem_comp_atom.atom_id"]
    if "_chem_comp_atom.comp_id" in block:
      comp_ids = [c.strip() for c in block["_chem_comp_atom.comp_id"]]
    elif block_name == "comp_%s" % resname:
      comp_ids = [resname] * len(names)
    else: continue
    xyz = zip(*[block["_chem_comp_atom.%s" % k] for k in "xyz"])
    sites = dict((name.strip('"'), tuple(float(v) for v in t))
      for comp_id, name, t in zip(comp_ids, names, xyz)
      if comp_id == resname and "?" not in t and "." not in t)
    if sites: return sites
  return None

def _ccd_describes(cc, resname):
  '''
  False if resname's restraints come from a user file: the CCD may use the code
  for another molecule (CCD LIG is C15H11N3). geostd/mon_lib and reduce2's own
  auto_<code> dictionaries are built from the CCD.
  '''
  if cc is None: return True
  parts = str(cc.source_info).replace("file:", "").strip().replace("\\", "/").split("/")
  return ("chem_data" in parts or
          (parts[-1].startswith("auto_") and parts[-1][5:].strip() == resname))

def _h_references(resname, mon_lib_srv, cache, n_h=2, n_heavy=2):
  '''
  CH2 centres (2, 2) or propellers (3, 1) with reference sites. CCD first: it
  defines PDB names, and geostd amino acids carry no coordinates. Then the
  restraint dictionary, for ligands whose H names differ from the CCD (VPH).
  geostd MAN names H61/H62 opposite to the CCD. A user restraint file is used
  alone.

  Sites are the CCD's ideal coordinates, except for the CH2 groups in
  _CCD_MODEL_CH2 (see there). Each set is used whole, never mixed in a group.
  For propellers the two CCD sets agree in 94% of groups, so the ideal set is
  taken as it comes.
  '''
  cache_key = (resname, n_h)
  if cache_key in cache: return cache[cache_key]
  from mmtbx.chemical_components import get_cif_dictionary
  result = []
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  cc_cif = None
  if _ccd_describes(cc, resname):
    try: cc_cif = get_cif_dictionary(resname)
    except Exception: cc_cif = None
  if cc_cif:
    elements = {}
    sites = {'model': {}, 'ideal': {}}
    for a in cc_cif.get('_chem_comp_atom', []):
      elements[a.atom_id] = a.type_symbol.strip().upper()
      for key, fmt in (('model', "model_Cartn_%s"),
                       ('ideal', "pdbx_model_Cartn_%s_ideal")):
        t = [getattr(a, fmt % k, "?") for k in "xyz"]
        if "?" not in t and "." not in t:
          sites[key][a.atom_id] = tuple(float(v) for v in t)
    bond_pairs = [(b.atom_id_1, b.atom_id_2)
      for b in cc_cif.get('_chem_comp_bond', [])]
    # the listed groups first, so name_prochiral_h takes them over the ideal set
    fix = _CCD_MODEL_CH2.get(resname, ()) if n_h == 2 else ()
    if fix:
      result += [g for g in _h_groups(elements, bond_pairs, sites['model'],
                                      n_h, n_heavy) if g[0] in fix]
    result += _h_groups(elements, bond_pairs, sites['ideal'], n_h, n_heavy)
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  if cc is not None:
    sites = _dictionary_sites(cc.source_info, resname)
    if sites:
      elements = dict((a.atom_id, a.type_symbol.strip().upper())
        for a in cc.atom_list)
      result += _h_groups(elements, [(b.atom_id_1, b.atom_id_2)
        for b in cc.bond_list], sites, n_h, n_heavy)
  cache[cache_key] = result
  return result

def name_prochiral_h(hierarchy, mon_lib_srv, kinds=(2, 3)):
  '''
  Riding places the two H of a CH2 in processing order, so about half get each
  other's name (1akg: 16 of 38 vs CCD), and it numbers a propeller (CH3, NH3)
  the other way round than the CCD does, every time: the H are still superposed
  when parameterization.check_propeller_order looks at them, so the order comes
  out of the riding frame instead. Swap names to match the reference chirality.
  Only names change; no atom moves, and the swap survives a later riding
  idealization, which follows the H positions it finds.

  kinds selects the groups: 2 the CH2 centres, 3 the propellers. The CH2 are
  named before exclude_H_on_links, which picks the H a link drops by position,
  the propellers after it: a linked methyl (BGS CS in 2b5z) is a CH2 by then,
  no longer matches a three-H reference and keeps the surviving names it had.

  Returns the swaps made, in order, as (parent, h1, h2) atoms; h1 and h2
  carry their new names.
  '''
  cache, done, swaps = {}, set(), []
  for m in hierarchy.models():
    for c in m.chains():
      for conformer in c.conformers():
        for r in conformer.residues():
          resname = r.resname.strip()
          groups = []
          if 2 in kinds:
            groups += [(2, g) for g in
                        _h_references(resname, mon_lib_srv, cache, 2, 2)]
          if 3 in kinds:
            groups += [(3, g) for g in
                        _h_references(resname, mon_lib_srv, cache, 3, 1)]
          if not groups: continue
          atoms = dict((a.name.strip(), a) for a in r.atoms())
          for n_h, (p, hv, hs, sites) in groups:
            if not all(n in atoms for n in (p,)+hv+hs): continue
            # first source wins. Keyed on the parent and its H atoms: a shared
            # parent can carry its own H in each conformer (partly split GOL)
            key = (atoms[p].i_seq,) + tuple(sorted(atoms[n].i_seq for n in hs))
            if key in done: continue
            # heavy neighbours 2.4 A: C-SE is 1.96 and a modelled one reaches
            # 2.04 (2h34 MSE), S-S 2.03; far enough below any non-bonded
            # contact. H 1.6 A: bonded to the parent, not an H of a neighbour
            # with the same name (HT7: CCD HA on CA, model HA on CB, 2.1 A)
            if any(atoms[p].distance(atoms[n]) > 2.4 for n in hv): continue
            if any(atoms[p].distance(atoms[n]) > 1.6 for n in hs): continue
            done.add(key)
            if n_h == 2:
              # the side of the two heavy neighbours the first H is named on
              ref = (hv[0], hv[1], hs[0])
              swap = (hs[0], hs[1])
            else:
              # the turning sense of the first two H, seen from the neighbour
              ref = (hs[0], hs[1], hv[0])
              swap = (hs[1], hs[2])
            v_ideal = _chiral_volume(sites[p], *[sites[n] for n in ref])
            v_model = _chiral_volume(atoms[p].xyz,
                                     *[atoms[n].xyz for n in ref])
            if abs(v_ideal) < 0.5 or abs(v_model) < 0.5: continue # flat
            if (v_ideal > 0) != (v_model > 0):
              h1, h2 = atoms[swap[0]], atoms[swap[1]]
              h1.name, h2.name = h2.name, h1.name
              swaps.append((atoms[p], h1, h2))
  return swaps

def _bond_orders(resname, mon_lib_srv, cache):
  '''
  Double (2) and triple (3) bonds by atom-name pair; anything else counts 1.
  Restraint dictionary first, CCD overrides (PDB names) unless the dictionary
  is a user file.
  '''
  if resname in cache: return cache[resname]
  orders = {}
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  if cc is not None:
    for b in cc.bond_list:
      o = {"double": 2, "triple": 3}.get(str(b.type).strip().lower())
      if o: orders[frozenset((b.atom_id_1, b.atom_id_2))] = o
  from mmtbx.chemical_components import get_cif_dictionary
  cc_cif = None
  if _ccd_describes(cc, resname):
    try: cc_cif = get_cif_dictionary(resname)
    except Exception: cc_cif = None
  if cc_cif:
    for b in cc_cif.get('_chem_comp_bond', []):
      key = frozenset((b.atom_id_1, b.atom_id_2))
      o = {"DOUB": 2, "TRIP": 3}.get(str(b.value_order).strip().upper())
      if o: orders[key] = o
      else: orders.pop(key, None)
  cache[resname] = orders
  return orders

class place_hydrogens():
  '''
  Add H atoms to a model

  Parameters
  ----------
  use_neutron_distances : bool
    use neutron distances instead of X-ray

  adp_scale : float
    scale factor for isotropic B of H atoms.
    B(H-atom) = adp_scale * B(parent non-H atom)

  keep_existing_H : bool
    keep existing H atoms in model, only place missing H

  na_terminal_h : bool
    place the terminal hydroxyl H of standard nucleic acids that the
    polymer dictionaries (geostd GD, CD ...) lack: HO3' on an O3' bonded to
    no other residue (a 3' end) and HO5' on the O5' of a residue without P
    (a 5'-OH end). pdb_interpretation restrains both. qmready=True uses it.
  n_terminal_charge : str
    where to place N-terminal NH3: 'residue_one' (residue with resseq 1),
    'first_in_chain' (first residue of each chain), 'all_chain_starts'
    (first residue of each chain and each residue after a chain break, i.e.
    N not peptide-bonded to C of preceding residue), 'no_charge' (nowhere)

  mon_lib_srv : monomer library server or None
    server for the dictionary queries; with None, a new one. The model's
    cached server is not used: it can know a CCD-only residue from an earlier
    placement whose placeholder the processed model no longer has.

  metal_bound_amines_as_nh2 : bool
    a coordinating NH3+ becomes NH2 (convert_metal_bound_amines())
  '''

# ------------------------------------------------------------------------------

  def __init__(self,
               model,
               use_neutron_distances = False,
               n_terminal_charge     = 'residue_one',
               adp_scale             = 1,
               exclude_water         = True,
               stop_for_unknowns     = False,
               keep_existing_H       = False,
               validate_e            = False,
               print_time            = False,
               mon_lib_srv           = None,
               na_terminal_h         = False,
               metal_bound_amines_as_nh2 = False):
    self.model                 = model
    self.na_terminal_h         = na_terminal_h
    self.use_neutron_distances = use_neutron_distances
    self.n_terminal_charge     = n_terminal_charge
    self.adp_scale             = adp_scale
    self.exclude_water         = exclude_water
    self.stop_for_unknowns     = stop_for_unknowns
    self.keep_existing_H       = keep_existing_H
    self.validate_e            = validate_e
    self.print_time            = print_time
    self.mon_lib_srv           = mon_lib_srv
    self.metal_bound_amines_as_nh2 = metal_bound_amines_as_nh2
    #
    self.no_H_placed_mlq        = list()
    self.site_labels_disulfides = list()
    self.site_labels_no_para    = list()
    self.site_labels_tertiary_amide = list()
    self.site_labels_missing_neighbor = list()
    self.residues_missing_neighbor    = list()
    self.residues_unknown_atoms       = list()
    self.site_labels_unknown_atoms    = list()
    # Names of restraint dictionaries auto-generated for unknown ligands during
    # placement; these are throwaway (purpose-built for H placement) and are
    # removed from the model before returning so they don't leak into
    # downstream geometry validation.
    self.auto_restraint_names   = list()
    # the received model's restraint objects plus the auto-generated ones
    self.restraint_objects      = None
    # server of the queries in this run (see mon_lib_srv)
    self.query_mon_lib_srv      = None
    # residues that got H but have no restraints in the processed model
    self.residues_no_restraints = list()
    self.residue_keys_new_h     = set()
    # (model, chain, resseq, icode, resname) of no_H_placed_mlq
    self.residue_keys_no_restraints = set()
    #self.charged_atoms          = list()
    self.sl_removed             = list()
    # thiols bound to metals (thiol_metal_sites): H not placed, with the
    # contact; messages of unresolved sites; supplied H kept, with the
    # contact and whether the S has a metal-coordination proxy
    self.site_labels_thiolate   = list()
    self.thiol_h_unresolved     = list()
    self.thiol_h_supplied_kept  = list()
    # the same for backbone amide N and T/U N3 bound to metals
    # (metal_nh_sites), and the other N-H next to a metal, kept
    self.site_labels_metal_nh   = list()
    self.metal_nh_unresolved    = list()
    self.metal_nh_supplied_kept = list()
    self.metal_nh_reported      = list()
    # (target list, reported_metal_h_rows()) for metal_nh_reported,
    # filled in _finish_metal_lists() from the H that are left
    self._metal_nh_reported_src = list()
    # metal_bound_amines_as_nh2: sites placed as NH2 (group_args, see
    # convert_metal_bound_amines()), messages of sites left as NH3+,
    # supplied NH3+ kept and chain starts not built as NH3, (label, contacts)
    self._metal_amine_sites      = None
    self.metal_amine_converted   = list()
    self.metal_amine_converted_text = list()
    self.metal_amine_unresolved  = list()
    self.metal_amine_supplied    = list()
    self.metal_amine_supplied_nh2 = list()
    # H placed or kept as supplied by convert_metal_bound_amines() (keys),
    # and the N it listed
    self._protected_h            = set()
    self._metal_amine_listed_n   = set()
    self.metal_amine_chain_starts = list()
    # H that riding H refused although their parent has all its heavy
    # neighbours
    self.site_labels_riding_refused = list()
    # how to drop the supplied H, for the report (set by reduce2)
    self.supplied_h_hint        = None
    # the input's own symmetry for contacts (contact_symmetry()); a caller
    # that boxed the model sets input_crystal_symmetry (reduce2: the output's,
    # None without a cell); left as from_model, it is read from the model
    self._contact_cs            = None
    self.input_crystal_symmetry = 'from_model'
    # renames of H that were in the input (keep_existing_H), in the order
    # made: (residue part as hd_exchange.site_key(), parent altloc, old name,
    # new name); each exchanges the two names
    self.renamed_existing_h     = list()
    self._existing_h            = set()
    # input positions of those H, by the same identity
    self._existing_h_xyz        = dict()
    self.n_H_initial            = 0
    self.n_H_final              = 0

    if self.print_time:
      self.time_rebox_model        = None
      self.time_remove_element_X   = None
      self.time_add_missing_H      = None
      self.time_terminal_propeller = None
      self.time_make_grm           = None
      self.time_remove_isolated    = None
      self.time_riding_manager     = None
      self.time_remove_H_nopara    = None
      self.time_reset              = None
      self.time_idealize           = None
      self.time_remove_H_on_links  = None

# ------------------------------------------------------------------------------

  def remove_auto_restraint_objects(self):
    '''
    Remove the placeholder restraint dictionaries that were auto-generated for
    unknown ligands during H placement (see add_missing_H_atoms_at_bogus_position).

    The grm has already been built at this point, so the placement restraints
    have served their purpose. We assign self.model._restraint_objects directly
    instead of calling set_restraint_objects(): the latter unsets the restraints
    manager, and we need the already-built grm (and riding-H manager) to remain
    intact for the rest of run() and for callers that reuse it.
    '''
    remove_auto_restraint_objects(self.model, self.auto_restraint_names)

# ------------------------------------------------------------------------------

  def run(self):
    '''
    Function that places H atoms
    '''
    self._run()
    self._finish_metal_lists()

  def _run(self):

    # Create symmetry if necessary
    # ------------------------------
    model_has_bogus_cs = False
    t0 = time.time()
    # before the box: contacts across the input's own symmetry only
    cs = self.input_crystal_symmetry
    if cs == 'from_model': cs = get_output_crystal_symmetry(self.model)
    self._contact_cs = contact_symmetry(cs)
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      # Box in place: callers get their coordinates back unmoved.
      self.model = shift_and_box_model(model = self.model, shift_model = False)
      model_has_bogus_cs = True
      #self.model.add_crystal_symmetry_if_necessary() # this is slower than shift_and_box_model!!!!
    self.time_rebox_model = round(time.time()-t0, 2)


    # Don't stop if model contains element X atoms
    # This needs more discussion before being made final
    # ---------------------------------------------
    t0 = time.time()
    if ' X' in self.model.get_hierarchy().atoms().extract_element():
      self.model = self.model.select(~self.model.selection('element X'))
    self.time_remove_element_X = round(time.time()-t0, 2)


    # Remove existing H if requested
    # ------------------------------
    self.model.get_xray_structure()
    self.n_H_initial = self.model.get_hd_selection().count(True)
    if not self.keep_existing_H:
      self.model = self.model.select(~self.model.get_hd_selection())
    else:
      self._existing_h_xyz = self._h_identities()
      self._existing_h = set(self._existing_h_xyz)

    # Add missing H atoms and place them at bogus position
    # ----------------------------------------------------
    t0 = time.time()
    pdb_hierarchy = self.add_missing_H_atoms_at_bogus_position(
      exclude_water = self.exclude_water)
    self.time_add_missing_H = round(time.time()-t0, 2)

    # Place N-terminal propeller hydrogens
    # TODO double check N-terminal position for PRO residues
    # ------------------------------------
    if self.n_terminal_charge in ['residue_one', 'first_in_chain',
                                  'all_chain_starts']:
      t0 = time.time()
      self.place_n_terminal_propeller(pdb_hierarchy = pdb_hierarchy)
      self.time_terminal_propeller = round(time.time()-t0, 2)

    pdb_hierarchy.sort_atoms_in_place()
    pdb_hierarchy.atoms().reset_serial()
    # Make new model obj and get restraints manager
    # ---------------------------------------------
    p = get_reduce_pdb_interpretation_params(self.use_neutron_distances)
    ro = self.restraint_objects
    t0 = time.time()
    self.model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = pdb_hierarchy,
      stop_for_unknowns = self.stop_for_unknowns,
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = ro,
      log               = null_out())
    # plain_pairs_radius: that table is for ADP restraints only and ignores
    # model boundaries
    self.model.process(pdb_interpretation_params=p,
                       make_restraints=True,
                       plain_pairs_radius=None,
                       # retain_zero_dihedrals=True,
                       )

    self.time_make_grm = round(time.time()-t0, 2)
    self.find_atoms_not_in_dictionary()
    # Drop the throwaway auto-generated ligand restraints now that the grm has
    # been built. They were only needed to place H on unknown ligands; keeping
    # them on the model would leak idealized placement-only geometry (generic
    # X-H lengths, esd=1/period=1 torsions) into any downstream
    # re-interpretation.
    self.remove_auto_restraint_objects()
    self.exclude_H_on_metal_thiolates()
    self.exclude_H_on_metal_bound_N()
    if self.metal_bound_amines_as_nh2:
      self.convert_metal_bound_amines()
    # before riding H: an H on an esterified O has no H-O-X angle to the second
    # heavy neighbour (2e4i BGM 14 -> BGM 15)
    self.exclude_H_on_esterified_O()
    self.exclude_H_without_link_angle()
    # Return if no H have been placed
    sel_h = self.model.get_hd_selection()
    if sel_h.count(True) == 0: return

    # Remove isolated H atoms
    # -----------------------
    # (when heavy atom is missing, H needs not to be placed)
    t0 = time.time()
    sel_isolated = self.model.isolated_atoms_selection()
    sel_lone_H = sel_h & sel_isolated
    self.find_residues_without_restraints(sel_isolated)
    # As h_parameterization will not include these, they can be removed in the
    # next step; for book-keeping it is useful to keep track of lone H as a
    # selection
    #if not sel_lone_H.all_eq(False):
    #  self.model = self.model.select(~sel_lone_H)
    self.time_remove_isolated = round(time.time()-t0, 2)

    sel_h = self.model.get_hd_selection()

    # Restraints a link could not resolve by atom name
    # -----------------------------------------------
    self.add_link_h_restraints()
    self.add_na_terminal_h_restraints()

    # Setup riding H manager
    # ----------------------
    t0 = time.time()
    self.setup_riding_h_manager()
    riding_h_manager = self.model.riding_h_manager
    if riding_h_manager is None:
      return
    self.time_riding_manager = round(time.time()-t0, 2)
    # Remove H that could not be parameterized
    # ----------------------------------------
    t0 = time.time()
    sel_h_in_para = flex.bool(
      [bool(x) for x in riding_h_manager.h_parameterization])
    sel_h_not_in_para = sel_h_in_para.exclusive_or(sel_h)

    water_selection = self.model.solvent_selection()

    if not self.exclude_water and water_selection.count(True)>0:
      workaround_002(
        model     = self.model,
        selection = water_selection.iselection(),
        use_neutron_distances = self.use_neutron_distances)
      water_selection = self.model.solvent_selection()
    # water H are placed by workaround_002, not by riding: keep them and
    # leave them out of the report of unplaced H. Only waters with an O: the
    # H of a water without O are removed like any other unplaced H.
    water_without_o = flex.bool(water_selection.size(), False)
    for rg in self.model.get_hierarchy().residue_groups():
      if get_class(name=rg.atom_groups()[0].resname) != "common_water": continue
      if any(a.element.strip().upper() == "O" for a in rg.atoms()): continue
      water_without_o.set_selected(rg.atoms().extract_i_seq(), True)
    sel_h_not_in_para = sel_h_not_in_para.set_selected(
      water_selection & ~water_without_o, False)
    # keep_existing_H: the H/D of residues without restraints stay as in the
    # input (lone H, not reported)
    sel_kept = self.unrestrained_hd_selection() | \
      self._protected_h_selection()
    sel_h_not_in_para.set_selected(sel_kept, False)
    sel_lone_H = sel_lone_H & ~sel_kept
    # no need to display lone H atoms in the log, so remove from labels
    sel_h_not_in_para_but_not_lone = sel_h_not_in_para.exclusive_or(sel_lone_H)
    # Classify the unplaceable H so the report reflects the real cause:
    #  - on a tertiary-amide backbone N (e.g. the ring N of an internal
    #    proline-type residue such as HYP): omitted because the N is already
    #    fully substituted.
    #  - missing a neighbouring heavy atom: the H's parent atom lacks an expected
    #    heavy neighbour (an incomplete side chain truncated before its next
    #    atom), so e.g. a methylene's rotation is undefined.
    #  - everything else: genuinely could not be parameterized.
    # Only walk the bonds when there is something to classify (usually none).
    self.site_labels_no_para = list()
    self.site_labels_tertiary_amide = list()
    self.site_labels_missing_neighbor = list()
    self.residues_missing_neighbor = list()
    self.site_labels_unknown_atoms = list()
    self.site_labels_riding_refused = list()
    unknown = set(r for r, names in self.residues_unknown_atoms)
    if not sel_h_not_in_para_but_not_lone.all_eq(False):
      grm = self.model.get_restraints_manager().geometry
      bps, asu = grm.get_all_bond_proxies(sites_cart=self.model.get_sites_cart())
      atoms = self.model.get_atoms()
      elements = self.model.get_hierarchy().atoms().extract_element()
      bonds = {}
      for proxy in list(bps) + list(asu):
        if   isinstance(proxy, ext.bond_simple_proxy): i,j = proxy.i_seqs
        elif isinstance(proxy, ext.bond_asu_proxy):    i,j = proxy.i_seq, proxy.j_seq
        else: continue
        bonds.setdefault(i, []).append(j)
        bonds.setdefault(j, []).append(i)
      tertiary = self.h_on_tertiary_amide_n(bonds, atoms, elements)
      def _heavy_neighbors(iseq):
        return [m for m in set(bonds.get(iseq, [])) if elements[m] not in ('H','D')]
      placed = flex.size_t(self.place_anchorless_h(
        sel_h_not_in_para_but_not_lone.iselection(), bonds, _heavy_neighbors))
      sel_h_not_in_para.set_selected(placed, False)
      sel_h_not_in_para_but_not_lone.set_selected(placed, False)
      seen_residues = set()
      for atom in self.model.get_hierarchy().atoms().select(
          sel_h_not_in_para_but_not_lone):
        label = atom.id_str().replace('pdb=','').replace('"','')
        if atom.i_seq in tertiary:
          self.site_labels_tertiary_amide.append(label)
          continue
        rg = atom.parent().parent()
        resid = '%s %s %s' % (atom.parent().resname.strip(),
          rg.parent().id.strip(), rg.resseq.strip())
        # its residue has atoms the dictionary does not know: not incomplete
        if resid in unknown:
          self.site_labels_unknown_atoms.append(label)
          continue
        parents = _heavy_neighbors(atom.i_seq)
        # parent atom missing an expected heavy neighbour -> incomplete residue
        if len(parents) == 1 and len(_heavy_neighbors(parents[0])) < 2:
          # unless it has all its dictionary neighbours: riding H refused the
          # H (an NH2 with an ideal H-N-H of 107-111 deg)
          expected = heavy_neighbors_in_restraints(atoms[parents[0]],
            self.query_mon_lib_srv)
          if (expected is not None and
              len(expected) <= len(_heavy_neighbors(parents[0]))):
            self.site_labels_riding_refused.append(label)
            continue
          self.site_labels_missing_neighbor.append(label)
          if resid not in seen_residues:
            seen_residues.add(resid)
            self.residues_missing_neighbor.append(resid)
        else:
          self.site_labels_no_para.append(label)
    # lone H of a water without O are left out of the labels above (lone H
    # are not displayed); report them as missing their heavy atom
    for atom in self.model.get_hierarchy().atoms().select(
        water_without_o & sel_lone_H & sel_h_not_in_para):
      self.site_labels_missing_neighbor.append(
        atom.id_str().replace('pdb=','').replace('"',''))
      rg = atom.parent().parent()
      resid = '%s %s %s' % (atom.parent().resname.strip(),
        rg.parent().id.strip(), rg.resseq.strip())
      if resid not in self.residues_missing_neighbor:
        self.residues_missing_neighbor.append(resid)
    if not sel_h_not_in_para.all_eq(False):
      self.model = self.model.select(~sel_h_not_in_para)
    self.time_remove_H_nopara = round(time.time()-t0, 2)
    # Reset occupancies, ADPs and idealize H atom positions
    # -----------------------------------------------------
    t0 = time.time()
    sel_kept = self.unrestrained_hd_selection()
    kept_occ = self.model.get_occ().select(sel_kept)
    kept_b = self.model.get_b_iso().select(sel_kept)
    self.model.reset_adp_for_hydrogens(scale = self.adp_scale, keep_aniso=True)
    self.model.reset_occupancy_for_hydrogens_simple()
    if sel_kept.count(True) > 0:
      self.model.set_occupancies(values = kept_occ, selection = sel_kept)
      self.model.set_b_iso(values = kept_b, selection = sel_kept,
        keep_aniso = True)
    self.time_reset = round(time.time()-t0, 2)
    t0 = time.time()
    # None, not False: callers outside reduce2 keep the scattering-table rule
    self.model.idealize_h_riding(
      use_neutron_distances = True if self.use_neutron_distances else None)
    self.keep_h_on_own_sites()
    self.time_idealize = round(time.time()-t0, 2)

    # Remove H atoms that are involved in links (bonds, metal coordination, etc)
    # --------------------------------------------------------------------------
    t0 = time.time()
    # CH2 names must be stereo-correct before a link picks which H to drop
    self._log_swaps(name_prochiral_h(self.model.get_hierarchy(),
      self.model.get_mon_lib_srv(), kinds = (2,)))
    self.exclude_H_on_links()
    # propellers only now: a methyl that lost an H to a link is a CH2, and the
    # names that survived it stay as they are
    self._log_swaps(name_prochiral_h(self.model.get_hierarchy(),
      self.model.get_mon_lib_srv(), kinds = (3,)))
    self.time_remove_H_on_links = round(time.time()-t0, 2)


    # TODO: this should be ideally done *after* reduce optimization
    #if not self.exclude_water:
    #  self.model.add_hydrogens(1., occupancy=0.)

    # Its cases are now handled by linking, bond orders and the HIS exception;
    # on 4us9 it removed the hemiaminal H that must stay (tst_add_hydrogen_12).
    #self.model, _ = workaround_003(model = self.model)

    # List missing H
    mon_lib_srv = self.model.get_mon_lib_srv()
    for m in self.model.get_hierarchy().models():
      for c in m.chains():
        for con in c.conformers():
          for r in con.residues():
            ma = r.missing_atoms(mon_lib_srv = mon_lib_srv, mode="h_only")
            if ma is not None and len(ma)>0:
              msg="chain %s resseq %s resname %s misses:"
              if 0: # Hold off printing untill verbosity is added
                print(msg%(c.id, r.resseq, r.resname), ma)
    self.n_H_final = self.model.get_hd_selection().count(True)

    if self.print_time:
      self.print_times()

  # ----------------------------------------------------------------------------

  def _h_identities(self):
    '''
    {(residue part, altloc, name): xyz} of every H/D; residue part as in
    hd_exchange.site_key()
    '''
    h = self.model.get_hierarchy()
    parts = residue_parts(h)
    return dict((parts[a.parent().memory_id()] + (a.parent().altloc, a.name),
                 a.xyz)
                for a in h.atoms() if a.element.strip().upper() in ('H', 'D'))

  # ----------------------------------------------------------------------------

  def keep_h_on_own_sites(self):
    '''
    keep_existing_H: riding idealization rebuilds the H of a parent from the
    ideal dihedral, so a kept propeller can come back turned by one site, each
    H on another one's site. Per parent and altloc, deal the group's ideal
    sites to its H so that the kept H move least in total (all permutations,
    first found on a tie); new H take the remaining sites.
    '''
    if not self._existing_h_xyz: return
    riding_h_manager = self.model.riding_h_manager
    if riding_h_manager is None: return
    h = self.model.get_hierarchy()
    atoms = h.atoms()
    parts = residue_parts(h)
    groups = {}
    for i_seq, rc in enumerate(riding_h_manager.h_parameterization):
      if rc is None: continue
      groups.setdefault((rc.a0, atoms[i_seq].parent().altloc), []).append(
        i_seq)
    sites_cart = self.model.get_sites_cart()
    changed = False
    for group in groups.values():
      if len(group) < 2: continue
      old = [self._existing_h_xyz.get(parts[atoms[i].parent().memory_id()] +
               (atoms[i].parent().altloc, atoms[i].name)) for i in group]
      if all(x is None for x in old): continue
      sites = [matrix.col(sites_cart[i]) for i in group]
      def moved(order):
        return sum((sites[j] - matrix.col(x)).length()
                   for x, j in zip(old, order) if x is not None)
      best = min(itertools.permutations(range(len(group))), key=moved)
      if list(best) == list(range(len(group))): continue
      for i, j in zip(group, best):
        sites_cart[i] = sites[j].elems
      changed = True
    if changed:
      self.model.set_sites_cart(sites_cart)

  # ----------------------------------------------------------------------------

  def _log_rename(self, parent, atom, old_name, partner=None, part=None):
    '''
    Log the rename of atom (it carries the new name) if it or its swap
    partner (now named old_name) was in the input.
    '''
    if not self._existing_h: return
    if part is None: part = residue_part(atom.parent())
    moved = [(atom, old_name)]
    if partner is not None: moved.append((partner, atom.name))
    old_keys = [part + (a.parent().altloc, n) for a, n in moved]
    existed = [k in self._existing_h for k in old_keys]
    if not any(existed): return
    for k in old_keys: self._existing_h.discard(k)
    for (a, n), e in zip(moved, existed):
      if e: self._existing_h.add(part + (a.parent().altloc, a.name))
    self.renamed_existing_h.append(
      (part, parent.parent().altloc, old_name, atom.name))

  # ----------------------------------------------------------------------------

  def _log_swaps(self, swaps):
    '''
    Log the swaps of name_prochiral_h() that involve input H
    '''
    if not (swaps and self._existing_h): return
    parts = residue_parts(self.model.get_hierarchy())
    for parent, h1, h2 in swaps:
      self._log_rename(parent, h1, h2.name, partner=h2,
        part=parts[h1.parent().memory_id()])

  # ----------------------------------------------------------------------------

  def unrestrained_hd_selection(self):
    '''
    With keep_existing_H, the H/D of the residues without restraints
    (no_H_placed_mlq), by residue_group and resname: the same name can have
    restraints elsewhere. Otherwise nothing.
    '''
    sel = flex.bool(self.model.size(), False)
    if not (self.keep_existing_H and self.residue_keys_no_restraints):
      return sel
    for m in self.model.get_hierarchy().models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            if ((m.id, chain.id, rg.resseq, rg.icode, ag.resname) not in
                self.residue_keys_no_restraints): continue
            for a in ag.atoms():
              if a.element_is_hydrogen(): sel[a.i_seq] = True
    return sel

  # ----------------------------------------------------------------------------

  def find_residues_without_restraints(self, sel_isolated):
    '''
    Residues that got H but have no bond at all in the processed model: their
    H go with the lone H, which are not reported otherwise. A truncated
    residue keeps the bonds of the atoms it has and is not listed.
    '''
    self.residues_no_restraints = list()
    for m in self.model.get_hierarchy().models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          if ((m.id, chain.id, rg.resseq, rg.icode) not in
              self.residue_keys_new_h): continue
          if not sel_isolated.select(rg.atoms().extract_i_seq()).all_eq(True):
            continue
          resid = '%s %s %s' % (rg.atom_groups()[0].resname.strip(),
            chain.id.strip(), rg.resseq.strip())
          if resid not in self.residues_no_restraints:
            self.residues_no_restraints.append(resid)

  # ----------------------------------------------------------------------------

  def add_na_terminal_h_restraints(self):
    '''
    The dihedral riding H needs for the HO5' of a 5'-OH end (na_terminal_h):
    pdb_interpretation gives a 3' end C4'-C3'-O3'-HO3' but a 5' end without
    P no C4'-C5'-O5'-HO5'. Added where missing (standard nucleotides only):
    staggered, ideal 180, esd 30, periodicity 3 (the geostd DNA dictionaries'
    value); atoms taken from the H's conformer.
    '''
    if not self.na_terminal_h: return
    atoms = self.model.get_hierarchy().atoms()
    hs = [a for a in atoms if a.name.strip() == "HO5'"
          and get_class(name=a.parent().resname) == 'common_rna_dna']
    if not hs: return
    grm = self.model.get_restraints_manager().geometry
    in_dihedral = set()
    for proxy in grm.get_dihedral_proxies():
      in_dihedral.update(proxy.i_seqs)
    new = geometry_restraints.shared_dihedral_proxy()
    for h in hs:
      if h.i_seq in in_dihedral: continue
      ag = h.parent()
      altloc = ag.altloc.strip()
      def find(name):
        # the H's conformer: blank or its altloc; any group for a blank H
        for g in ag.parent().atom_groups():
          if altloc and g.altloc.strip() not in ('', altloc): continue
          for a in g.atoms():
            if a.name.strip() == name: return a
        return None
      o5, c5, c4 = find("O5'"), find("C5'"), find("C4'")
      if None in (o5, c5, c4): continue
      new.append(geometry_restraints.dihedral_proxy(
        i_seqs      = [c4.i_seq, c5.i_seq, o5.i_seq, h.i_seq],
        angle_ideal = 180.0,
        weight      = 1./30.**2,
        periodicity = 3))
    # the in_dihedral skip keeps this free of duplicates
    if new.size():
      grm.add_dihedrals_in_place(new, check_for_duplicates=False)

  def add_link_h_restraints(self):
    '''
    The peptide link defines the angle C-N-H (TRANS, 124.3 deg) and the amide
    plane C-N-CA-H, both by atom name.
    A residue in the chain whose amide H is called anything else - 9FZ H3,
    9G2 H6, 1133 of the 2211 peptide-like geostd entries with an N-H - gets no
    such restraint: pdb_interpretation only retries the v2/v3 spelling of the
    name and then drops the definition. Nothing then fixes the H direction,
    riding calls it a rotatable amine H, and it can end up on the preceding
    carbonyl C (5nxq, 17 deg). Add both to this model's restraints. The plane
    is added whole: without the H it holds three atoms, which are coplanar
    whatever they do.

    Only reduce2 places H from these restraints, so this is not the whole fix:
    refinement has the same gap, and the proper place for it is the name
    resolution in pdb_interpretation.
    '''
    srv = self.model.get_mon_lib_srv()
    link = getattr(srv, 'link_link_id_dict', {}).get('TRANS')
    if link is None: return
    defs = [a for a in link.angle_list
            if a.atom_id_2.strip() == 'N' and a.atom_id_3.strip() in ('H', 'D')
            and a.value_angle is not None and a.value_angle_esd]
    if not defs: return
    angle_ideal, esd = defs[0].value_angle, defs[0].value_angle_esd
    # CA-N-H: the dictionary of a free amino acid gives the sp3 amine value
    # (9FZ 108.5), which does not belong to an amide N. In the plane the three
    # angles add up to 360, so the link fixes this one too.
    ca = [a for a in link.angle_list
          if a.atom_id_2.strip() == 'N' and a.atom_id_3.strip() == 'CA'
          and a.value_angle is not None]
    ca_ideal = 360. - angle_ideal - ca[0].value_angle if ca else None
    # the plane that carries the H, if the link has one (TRANS plane2)
    planes = {}
    for p in getattr(link, 'plane_list', []):
      planes.setdefault(p.plane_id, []).append(p)
    plane_esd = None
    for rows in planes.values():
      names = set(r.atom_id.strip() for r in rows)
      if names >= set(['C', 'N', 'CA']) and names & set(['H', 'D']):
        plane_esd = min(r.dist_esd for r in rows if r.dist_esd)
        break
    grm = self.model.get_restraints_manager().geometry
    atoms = self.model.get_atoms()
    elements = self.model.get_hierarchy().atoms().extract_element()
    bps, asu = grm.get_all_bond_proxies(sites_cart = self.model.get_sites_cart())
    bonds = {}
    for proxy in list(bps) + list(asu):
      if   isinstance(proxy, ext.bond_simple_proxy): i,j = proxy.i_seqs
      elif isinstance(proxy, ext.bond_asu_proxy):    i,j = proxy.i_seq, proxy.j_seq
      else: continue
      bonds.setdefault(i, []).append(j)
      bonds.setdefault(j, []).append(i)
    new, planarities, stale = [], [], []
    known = set()
    for proxy in grm.angle_proxies:
      i, j, k = proxy.i_seqs
      known.add((min(i, k), j, max(i, k)))
    # H/A must not get C/B of the previous residue (7a5m ALA A 103: N without
    # altloc, CA and H in A and B)
    altloc = lambda i: atoms[i].parent().altloc.strip()
    same_conf = lambda i, j: '' in (altloc(i), altloc(j)) or altloc(i) == altloc(j)
    for i_h, atom in enumerate(atoms):
      if elements[i_h].strip() not in ('H', 'D'): continue
      heavy = [j for j in set(bonds.get(i_h, []))
               if elements[j].strip() not in ('H', 'D')]
      if len(heavy) != 1: continue
      i_n = heavy[0]
      if atoms[i_n].name.strip() != 'N': continue
      rg_n = atoms[i_n].parent().parent()
      for i_c in set(bonds.get(i_n, [])):
        if atoms[i_c].name.strip() != 'C': continue
        if atoms[i_c].parent().parent().memory_id() == rg_n.memory_id(): continue
        if not same_conf(i_c, i_h): continue
        key = (min(i_c, i_h), i_n, max(i_c, i_h))
        if key in known: continue
        known.add(key)
        new.append(geometry_restraints.angle_proxy(
          i_seqs      = [i_c, i_n, i_h],
          angle_ideal = angle_ideal,
          weight      = 1./esd**2))
        i_ca = [j for j in set(bonds.get(i_n, []))
                if atoms[j].name.strip() == 'CA'
                and atoms[j].parent().parent().memory_id() == rg_n.memory_id()
                and same_conf(j, i_h)]
        # a covalent link (cyclic peptide, 1evd MDH 5 C - ACB 1 N) has the bond
        # but not C-N-CA, which riding H needs for an H with two heavy neighbours
        if ca and ca[0].value_angle_esd and i_ca:
          key_ca = (min(i_c, i_ca[0]), i_n, max(i_c, i_ca[0]))
          if key_ca not in known:
            known.add(key_ca)
            new.append(geometry_restraints.angle_proxy(
              i_seqs      = [i_c, i_n, i_ca[0]],
              angle_ideal = ca[0].value_angle,
              weight      = 1./ca[0].value_angle_esd**2))
        if plane_esd and i_ca:
          i_seqs = [i_c, i_n, i_ca[0], i_h]
          planarities.append(geometry_restraints.planarity_proxy(
            i_seqs  = flex.size_t(i_seqs),
            weights = flex.double([1./plane_esd**2]*len(i_seqs))))
        if ca_ideal is not None and i_ca:
          stale.append((i_ca[0], i_n, i_h))
    # replace the amine CA-N-H of each of those residues by the amide value:
    # proxy_remove drops a proxy only when all of its atoms are selected, so
    # selecting the three atoms takes that angle and nothing else
    for i_ca, i_n, i_h in stale:
      weight = [p.weight for p in grm.angle_proxies
                if tuple(sorted(p.i_seqs)) == tuple(sorted((i_ca, i_n, i_h)))]
      if not weight: continue
      sel = flex.bool(len(atoms), False)
      for i in (i_ca, i_n, i_h): sel[i] = True
      grm.remove_angles_in_place(selection=sel)
      new.append(geometry_restraints.angle_proxy(
        i_seqs      = [i_ca, i_n, i_h],
        angle_ideal = ca_ideal,
        weight      = weight[0]))
    if new: grm.add_angles_in_place(new)
    if planarities: grm.add_planarities_in_place(planarities)

  def exclude_H_on_esterified_O(self):
    '''
    Remove H from an O that has two heavy neighbours: the dictionary of a free
    molecule gives the O a hydroxyl H, but here it is an ester or an ether -
    O3' of a nucleotide inside a chain carries the next phosphate. A polymer
    bond has origin_id 0, so exclude_H_on_links never sees these.
    Also remove H from the terminal O/S of a phosphodiester that bridges to
    another residue (pKa ~1): free-nucleotide dictionaries carry it (CCD BGM
    HOP2 in 2e4i, geostd GS HSP2 in 5zmn).
    Diesters within one residue keep their dictionary H (35G; 58 geostd
    entries have one).
    '''
    grm = self.model.get_restraints_manager().geometry
    bps, asu = grm.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    atoms = self.model.get_hierarchy().atoms()
    elements = atoms.extract_element()
    is_h = lambda i: elements[i].strip() in ('H', 'D')
    # metal coordination is not an ester bond (A4: O-H on a metal keeps its H)
    metal = linking_class().get_origin_id('metal coordination')
    # one pass over the bonds, keeping only what O, S and P carry
    h_on, heavy_on = {}, {}
    for proxy in list(bps) + list(asu):
      if proxy.origin_id == metal: continue
      if   isinstance(proxy, ext.bond_simple_proxy): i,j = proxy.i_seqs
      elif isinstance(proxy, ext.bond_asu_proxy):    i,j = proxy.i_seq, proxy.j_seq
      else: continue
      for i_x, other in ((i, j), (j, i)):
        if elements[i_x].strip() not in ('O', 'S', 'P'): continue
        if is_h(other): h_on.setdefault(i_x, set()).add(other)
        else:           heavy_on.setdefault(i_x, set()).add(other)
    remove = set()
    for i_o, hs in h_on.items():
      if elements[i_o].strip() == 'O' and len(heavy_on.get(i_o, ())) > 1:
        remove.update(hs)
    residue = lambda i: atoms[i].parent().parent().memory_id()
    for i_p, neighbors in heavy_on.items():
      if elements[i_p].strip() != 'P': continue
      os_ = [k for k in neighbors if elements[k].strip() in ('O', 'S')]
      bridged = dict((k, heavy_on.get(k, set()) - {i_p}) for k in os_)
      bridged = dict((k, v) for k, v in bridged.items() if v)
      if len(bridged) < 2: continue
      r_p = residue(i_p)
      if all(residue(m) == r_p for v in bridged.values() for m in v): continue
      for k in os_:
        if k not in bridged: remove.update(h_on.get(k, ()))
    remove = flex.size_t(sorted(remove))
    if remove.size():
      sel = flex.bool(elements.size(), True).set_selected(remove, False)
      self.model = self.model.select(sel)

  def exclude_H_on_metal_thiolates(self):
    '''
    A thiol bound to a metal is taken as a thiolate (A4): an H that reduce2
    built on it is removed when thiol_metal_sites() finds the site covered
    (thiolate), and kept and reported when it is not (unresolved). Supplied
    H (keep_existing_H) are kept, with a note. Geometric contacts and
    metal-coordination proxies decide alike; exclude_H_on_links() leaves
    metal coordination to this step.
    '''
    self._exclude_H_at_metal_sites(thiol_metal_sites, METAL_ELEMENTS,
      self.site_labels_thiolate, self.thiol_h_unresolved,
      self.thiol_h_supplied_kept)

  def exclude_H_on_metal_bound_N(self):
    '''
    A backbone amide N (with a preceding residue) or a thymine/uracil N3 bound
    to a metal of NH_METAL_ELEMENTS is taken as deprotonated (A5): as
    exclude_H_on_metal_thiolates(), with metal_nh_sites(). Other N-H next to
    such a metal (not His ring N) are kept and listed.
    '''
    self._exclude_H_at_metal_sites(metal_nh_sites, NH_METAL_ELEMENTS,
      self.site_labels_metal_nh, self.metal_nh_unresolved,
      self.metal_nh_supplied_kept, self.metal_nh_reported)

  def _atom_key(self, parts, a):
    return parts[a.parent().memory_id()] + (a.parent().altloc, a.name)

  def _nh3_plus(self, n, i_n):
    return dictionary_nh3_plus(self.model, n, i_n, self.query_mon_lib_srv)

  def _h_order(self, n):
    # dictionary order of the H on n; H1-H3 of an N-terminal propeller first
    try:
      cc = self.query_mon_lib_srv.get_comp_comp_id_direct(
        n.parent().resname.strip())
    except Exception:
      cc = None
    names = [a.atom_id.strip() for a in cc.atom_list] if cc is not None else []
    def order(h):
      name = h.name.strip()
      if n.name.strip() == 'N' and name in ('H1', 'H2', 'H3'):
        return (0, 0, name)
      return (1, names.index(name) if name in names else len(names), name)
    return order

  def convert_metal_bound_amines(self):
    """
    metal_bound_amines_as_nh2 (A5 C1b): a primary amine (metal_amine_sites)
    that the dictionary makes NH3+ (three H; energy type NT3 or formal charge
    +1) and whose H sites are covered by one metal position (also through
    symmetry: the contact's image) is placed as NH2. place_metal_nh2()
    constructs and checks the two H first; only an accepted site loses the
    first H in dictionary order (H1 of an N-terminal propeller) and gets the
    coordinates; the two H keep their names and are left out of riding H
    (setup_riding_h_manager()). Not accepted: unresolved, listed, no H
    removed.
    No supplied H (keep_existing_H) is deleted. H are counted per complete
    conformer (conformer_groups()). Two supplied H with any metal in the
    window keep their coordinates and names, are not completed to NH3 and
    are not checked, whatever the coverage or frame (with approach=optimize:
    protect_supplied_metal_amines()), except when the H of a conformer are
    split over blank and alternate groups: such a conformer is unresolved
    and completed as with the option off, the supplied H re-idealized. With
    approach=optimize, a blank H shared with an NH3 conformer also belongs
    to that conformer and may be moved by its rotator; its input position is
    not guaranteed. The NH2 conformer's own alternate H keep their input
    positions. With one supplied
    H the second is constructed around it, the supplied H fixed; not
    accepted: the supplied H is retained and re-idealized, the group
    completed as with the option off. Three supplied H stay and are listed
    (riding H idealizes them as without the option). Left as they are and
    listed: alternate metal positions or alternate positions of the heavy
    neighbour (no unique frame), sites with unresolved protonation; an
    amino-acid N without a preceding residue that does not carry three H
    (chain start, terminus not confirmed). A neutral NH2 keeps both H (step
    B lists it).
    """
    elements = self.model.get_hierarchy().atoms().extract_element()
    if not any(e.strip().upper() in NH_METAL_ELEMENTS for e in elements):
      return
    sites = metal_amine_sites(self.model, self._contact_cs,
      mon_lib_srv = self.query_mon_lib_srv)
    if not sites.h_sites: return
    hierarchy = self.model.get_hierarchy()
    atoms = hierarchy.atoms()
    parts = residue_parts(hierarchy)
    key = lambda a: self._atom_key(parts, a)
    supplied = lambda a: self.keep_existing_H and key(a) in self._existing_h
    in_window = {}
    for c in sites.contacts:
      in_window.setdefault(c.i_x, []).append(c)
    by_n = {}
    for site in sites.h_sites:
      by_n.setdefault(site.i_x, []).append(site)
    grm = self.model.get_restraints_manager().geometry
    bond_ideal = {}
    for p in grm.get_all_bond_proxies(sites_cart = self.model.get_sites_cart())[0]:
      bond_ideal[tuple(sorted(p.i_seqs))] = p.distance_ideal
    angle_ideal = {}
    for p in grm.angle_proxies:
      i, j, k = p.i_seqs
      angle_ideal[(min(i, k), j, max(i, k))] = p.angle_ideal
    remove, install = [], {}
    for i_n in sorted(by_n):
      n = atoms[i_n]
      label = _atom_label(n)
      # per complete conformer (blank H counted in each)
      for alt, group, mixed in conformer_groups(by_n[i_n],
          lambda s: atoms[s.i_h].parent().altloc):
        hs = [atoms[s.i_h] for s in group]
        own = [s for s in group if atoms[s.i_h].parent().altloc == alt] or group
        contacts = list(own[0].contacts)
        window = in_window.get(i_n, [])
        if not window: continue
        contact = '; '.join(sites.contact_text(c) for c in (contacts or window))
        if len(hs) != 3 or not self._nh3_plus(n, i_n):
          if (contacts and len(hs) < 3 and n.name.strip() == 'N' and
              residue_is_amino_acid(n.parent().resname,
                                    self.query_mon_lib_srv)):
            self.metal_amine_chain_starts.append((label, contact))
          continue
        order = self._h_order(n)
        sup = sorted([h for h in hs if supplied(h)], key=order)
        new = sorted([h for h in hs if not supplied(h)], key=order)
        if len(sup) == 3:
          self.metal_amine_supplied.append(('%s (%s)' % (label,
            ', '.join(_h_name(h) for h in sup)), contact))
          continue
        if mixed:
          # one H shared by the conformers: no H removed or placed
          self.metal_amine_unresolved.append('%s: H in blank and alternate '
            'conformers: completed and re-idealized as with the option off'
            % label)
          continue
        if len(sup) == 2:
          # kept as supplied, never completed to NH3; not placed or checked
          remove.extend(h.i_seq for h in new)
          self._protected_h.update(key(h) for h in sup)
          self._metal_amine_listed_n.add(key(n))
          self.metal_amine_supplied_nh2.append(('%s (%s)' % (label,
            ', '.join(_h_name(h) for h in sup)), contact))
          continue
        verdicts = set(s.verdict for s in group)
        if verdicts == set(['not counted']): continue
        reason, message, placed = None, None, None
        if verdicts != set([sites.covered]):
          s = [s for s in group if s.verdict != sites.covered][0]
          # the message names the H: listed as it is
          message = sites.unresolved_message(s)
        elif len(contacts) != 1:
          reason = 'no unique metal position (%s)' % contact
        elif sites.amine[group[0].i_h].i_c is None:
          reason = 'alternate positions of its heavy neighbour; no unique frame'
        else:
          i_c = sites.amine[group[0].i_h].i_c
          # a supplied H first: it stays where it is
          kept = sup + new[1:]
          targets = group_args(
            d = [bond_ideal.get(tuple(sorted((i_n, h.i_seq)))) for h in kept],
            cnh = [angle_ideal.get((min(i_c, h.i_seq), i_n, max(i_c, h.i_seq)))
                   for h in kept],
            hnh = angle_ideal.get((min(kept[0].i_seq, kept[1].i_seq), i_n,
                                   max(kept[0].i_seq, kept[1].i_seq)), 109.47))
          if None in targets.d + targets.cnh:
            reason = 'no ideal N-H or C-N-H for the H'
          else:
            placed, reason = place_metal_nh2(matrix.col(n.xyz),
              matrix.col(atoms[i_c].xyz), matrix.col(contacts[0].image),
              targets, fixed = matrix.col(sup[0].xyz) if sup else None)
        if placed is None:
          text = message or '%s: %s' % (label, reason)
          if sup:
            text += ('; %s retained and re-idealized (completed as with the '
              'option off)' % sup[0].name.strip())
          self.metal_amine_unresolved.append(text)
          continue
        drop = new[:1]
        remove.extend(h.i_seq for h in drop)
        for h, xyz in zip(kept, placed):
          if not supplied(h): install[h.i_seq] = xyz.elems
          self._protected_h.add(key(h))
        self._metal_amine_listed_n.add(key(n))
        text = '%s: %s not placed' % (label, ', '.join(h.name.strip()
          for h in drop))
        if sup:
          text += '; %s kept as supplied, %s placed' % (sup[0].name.strip(),
            kept[1].name.strip())
        self.metal_amine_converted.append(group_args(label = label,
          contact = contact, text = text, n = key(n),
          kept = [key(h) for h in kept], xyz = [x.elems for x in placed]))
    if install:
      sites_cart = self.model.get_sites_cart()
      for i, xyz in install.items():
        sites_cart[i] = xyz
      self.model.set_sites_cart(sites_cart)
    if remove:
      sel = flex.bool(atoms.size(), True).set_selected(
        flex.size_t(sorted(remove)), False)
      self.model = self.model.select(sel)

  def _protected_h_selection(self):
    # the H that convert_metal_bound_amines() placed or kept as supplied:
    # riding H leaves them where they are
    sel = flex.bool(self.model.size(), False)
    if not self._protected_h: return sel
    hierarchy = self.model.get_hierarchy()
    parts = residue_parts(hierarchy)
    for a in hierarchy.atoms():
      if self._atom_key(parts, a) in self._protected_h:
        sel[a.i_seq] = True
    return sel

  def setup_riding_h_manager(self):
    """
    Riding H manager, idealized at once (as model.setup_riding_h_manager()).
    The H of convert_metal_bound_amines() are taken out of it before that
    first idealization; all other H keep their parameterization.
    """
    protected = self._protected_h_selection()
    if protected.count(True) == 0:
      self.model.setup_riding_h_manager(use_ideal_dihedral = True)
      return
    self.model.setup_riding_h_manager(use_ideal_dihedral = True,
      idealize = False)
    rhm = self.model.riding_h_manager
    if rhm is None: return
    for i in protected.iselection():
      rhm.h_parameterization[i] = None
    rhm.parameterization_cpp = rhm.get_parameterization_cpp(
      h_parameterization = rhm.h_parameterization)
    self.model.idealize_h_riding()

  def _finish_metal_lists(self):
    # end of run(): the lists name what is left
    hierarchy = self.model.get_hierarchy()
    parts = residue_parts(hierarchy)
    index = dict((self._atom_key(parts, a), a.i_seq) for a in hierarchy.atoms())
    # by the identities when the list was made (prochiral naming only swaps
    # names within a group)
    for target, rows in self._metal_nh_reported_src:
      for x_key, label, hs, contacts in rows:
        names = [name for k, name in hs if k in index]
        if names and x_key not in self._metal_amine_listed_n:
          target.append(('%s (%s)' % (label, ', '.join(names)), contacts))
    self._metal_nh_reported_src = []
    self.metal_amine_converted_text = [(s.text, s.contact)
      for s in self.metal_amine_converted]

  def _exclude_H_at_metal_sites(self, site_class, metals, not_placed,
                                unresolved, supplied_kept, reported=None):
    # The H that reduce2 built at a covered site of site_class are removed and
    # listed (not_placed: label, contacts); unresolved sites keep the H, with
    # the message; supplied H are kept (supplied_kept); reported: the
    # candidates site_class only reports, listed in _finish_metal_lists().
    elements = self.model.get_hierarchy().atoms().extract_element()
    if not any(e.strip().upper() in metals for e in elements): return
    sites = site_class(self.model, self._contact_cs)
    hierarchy = self.model.get_hierarchy()
    atoms = hierarchy.atoms()
    parts = residue_parts(hierarchy)
    if reported is not None:
      self._metal_nh_reported_src.append((reported, reported_metal_h_rows(
        sites, lambda a: self._atom_key(parts, a))))
    if not sites.h_sites: return
    def supplied(a):
      return (self.keep_existing_H and parts[a.parent().memory_id()] +
        (a.parent().altloc, a.name) in self._existing_h)
    supplied_kept.extend(supplied_metal_h(sites, supplied))
    remove = []
    for site in sites.h_sites:
      if site.verdict == 'not counted': continue
      h = atoms[site.i_h]
      label = _atom_label(h)
      contact = '; '.join(sites.contact_text(c) for c in site.contacts)
      if supplied(h):
        continue
      elif site.verdict == sites.covered:
        remove.append(site.i_h)
        not_placed.append((label, contact))
      else:
        unresolved.append(sites.unresolved_message(site))
    if remove:
      sel = flex.bool(atoms.size(), True).set_selected(
        flex.size_t(sorted(remove)), False)
      self.model = self.model.select(sel)

  def find_atoms_not_in_dictionary(self):
    '''
    For the report: see atoms_not_in_dictionary().
    '''
    self.residues_unknown_atoms = atoms_not_in_dictionary(self.model)

  def exclude_H_without_link_angle(self):
    '''
    A link gives an angle only to the H it keeps (NAG-ASN: HD21-ND2-C1) and
    deletes the other (HD22), which exclude_H_on_links removes only after
    riding H. Riding H stops on that mismatch when the linked atom is itself
    the second neighbour of an H (2fd6: NAG C1 still carries O1 and HO1).
    Remove the H without the angle when another H on the same parent has it.
    Only neighbours in another residue: a gap within one dictionary still
    raises.
    '''
    grm = self.model.get_restraints_manager().geometry
    bps, asu = grm.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    atoms = self.model.get_hierarchy().atoms()
    elements = atoms.extract_element()
    is_h = lambda i: elements[i].strip() in ('H', 'D')
    residue = lambda i: atoms[i].parent().parent().memory_id()
    metal = linking_class().get_origin_id('metal coordination')
    parent, heavy = {}, {}
    for proxy in list(bps) + list(asu):
      # metal coordination takes no H (A4)
      if proxy.origin_id == metal: continue
      if   isinstance(proxy, ext.bond_simple_proxy): i,j = proxy.i_seqs
      elif isinstance(proxy, ext.bond_asu_proxy):    i,j = proxy.i_seq, proxy.j_seq
      else: continue
      if is_h(i) and is_h(j): continue
      if is_h(i):   parent[i] = j
      elif is_h(j): parent[j] = i
      elif residue(i) != residue(j):
        heavy.setdefault(i, set()).add(j)
        heavy.setdefault(j, set()).add(i)
    if not heavy: return
    angle_to = {}
    for ap in grm.get_all_angle_proxies():
      i, j, k = ap.i_seqs
      for h, x in ((i, k), (k, i)):
        if is_h(h) and parent.get(h) == j and not is_h(x):
          angle_to.setdefault(h, set()).add(x)
    hs_on = {}
    for h, p in parent.items():
      if p in heavy: hs_on.setdefault(p, []).append(h)
    altloc = lambda i: atoms[i].parent().altloc.strip()
    remove = set()
    for p, hs in hs_on.items():
      if len(hs) < 2: continue
      for b in heavy[p]:
        # one H per altloc, each with its own neighbour (7a5m ALA A 103)
        hs_b = [h for h in hs
                if '' in (altloc(h), altloc(b)) or altloc(h) == altloc(b)]
        kept = [h for h in hs_b if b in angle_to.get(h, ())]
        if kept: remove.update(h for h in hs_b if h not in kept)
    if remove:
      sel = flex.bool(elements.size(), True).set_selected(
        flex.size_t(sorted(remove)), False)
      self.model = self.model.select(sel)

  def place_anchorless_h(self, iseqs, bonds, heavy_neighbors):
    '''
    Place H that riding H cannot parameterize because the molecule has only two
    heavy atoms (PEO, MOH): no third atom anchors the dihedral. Ideal bond and
    angle; H on the lower-i_seq atom at 0/120/240 deg, on the other at 180 deg
    (staggered/trans); They stay non-riding.
    Returns the i_seqs placed.
    '''
    grm = self.model.get_restraints_manager().geometry
    mon_lib_srv = self.model.get_mon_lib_srv()
    atoms = self.model.get_hierarchy().atoms()
    def _dict_heavy_degree(atom):
      # heavy neighbours of this atom in its dictionary, None if unknown
      try: cc = mon_lib_srv.get_comp_comp_id_direct(atom.parent().resname.strip())
      except Exception: cc = None
      if cc is None: return None
      ad = cc.atom_dict()
      name = atom.name.strip()
      if name not in ad: return None
      n = 0
      for b in cc.bond_list:
        for x, y in [(b.atom_id_1, b.atom_id_2), (b.atom_id_2, b.atom_id_1)]:
          if x == name and y in ad and ad[y].type_symbol not in ('H', 'D'): n += 1
      return n
    by_parent = {}
    for ih in iseqs:
      parents = heavy_neighbors(ih)
      if len(parents) != 1: continue
      i0 = parents[0]
      partners = heavy_neighbors(i0)
      if len(partners) != 1: continue
      i1 = partners[0]
      if heavy_neighbors(i1) != [i0]: continue
      # a truncated residue can also leave two bonded atoms: trust the dictionary
      if _dict_heavy_degree(atoms[i0]) != 1 or _dict_heavy_degree(atoms[i1]) != 1:
        continue
      by_parent.setdefault((i0, i1), []).append(ih)
    if not by_parent: return []
    targets = set(ih for hs in by_parent.values() for ih in hs)
    dist_ideal, angle_ideal = {}, {}
    sites_cart = self.model.get_sites_cart()
    bps, asu = grm.get_all_bond_proxies(sites_cart=sites_cart)
    for p in bps:
      for ih in set(p.i_seqs) & targets:
        dist_ideal[ih] = p.distance_ideal
    for p in grm.angle_proxies:
      if p.i_seqs[0] in targets: angle_ideal[p.i_seqs[0]] = p.angle_ideal
      if p.i_seqs[2] in targets: angle_ideal[p.i_seqs[2]] = p.angle_ideal
    placed = []
    for (i0, i1), hs in by_parent.items():
      r0, r1 = matrix.col(sites_cart[i0]), matrix.col(sites_cart[i1])
      # frame shared by both ends of the fragment
      axis = (r1 - r0) if i0 < i1 else (r0 - r1)
      p = axis.ortho().normalize()
      q = axis.normalize().cross(p)
      u = (r1 - r0).normalize()
      offset = 0 if i0 < i1 else 180
      for k, ih in enumerate(sorted(hs)):
        if ih not in dist_ideal or ih not in angle_ideal: continue
        a = math.radians(angle_ideal[ih])
        phi = math.radians(offset + 120*k)
        d = u*math.cos(a) + (p*math.cos(phi) + q*math.sin(phi))*math.sin(a)
        sites_cart[ih] = r0 + d*dist_ideal[ih]
        placed.append(ih)
    self.model.set_sites_cart(sites_cart)
    return placed

  # ----------------------------------------------------------------------------

  def place_n_terminal_propeller(self, pdb_hierarchy):
    '''
    Place NH3 at residue #1, at first residue in chain, or at all chain starts
    (first residue in chain and first residue after a chain break)
    Changes hierarchy in place
    '''
    # an N bonded to another residue is no terminus (cyclic peptide: 1evd ACB 1
    # N - MDH 5 C)
    for m in pdb_hierarchy.models():
      model_atoms = m.atoms()
      model_sites = model_atoms.extract_xyz()
      if self.n_terminal_charge == 'all_chain_starts':
        rgs_list = chain_start_residue_groups(m)
      else:
        rgs_list = [chain.residue_groups()[0] for chain in m.chains()]
      candidates = []
      for rgs in rgs_list:
        # by default, place NH3 only at residue with resseq 1
        if (self.n_terminal_charge == 'residue_one' and rgs.resseq_as_int() != 1):
          continue
        # SAC in 5xdq, 5zcp. Never needs propeller. Also AYA
        n = None
        for ag in rgs.atom_groups():
          n = ag.get_atom('N') # assumes atom name "N"
          if n: break
        if not n: continue
        bonds = heavy_neighbors_in_restraints(n, self.query_mon_lib_srv)
        if bonds is None: continue
        heavies = 2
        if ag.resname in ['PRO']: # needs a PRO child lookup
          heavies = 3
        if len(bonds) >= heavies: continue
        candidates.append((rgs, ag, n))
      if not candidates: continue
      # Atoms within 1.75 A of candidate N. Few candidates: scan all atoms per
      # N. Many (all_chain_starts, models can have very many chain starts):
      # one spatial lookup for all N, as scan is O(N_atoms) per N. Crossover
      # is at ~15 candidates.
      cutoff_sq = 1.75**2
      neighbors = None
      if len(candidates) > 10:
        n_sites = flex.vec3_double([n.xyz for rgs, ag, n in candidates])
        neighbors = scitbx.cubicle_neighbors.cubicle_neighbors(
          main_sites_cart = model_sites,
          cubicle_edge    = 5).neighbors_of(
            other_sites_cart   = n_sites,
            distance_cutoff_sq = cutoff_sq)
      for i_n, (rgs, ag, n) in enumerate(candidates):
        if neighbors is None:
          nb = ((model_sites - n.xyz).dot() < cutoff_sq).iselection()
        elif i_n in neighbors:
          nb = flex.size_t(list(neighbors[i_n]))
          nb = nb.select((model_sites.select(nb) - n.xyz).dot() < cutoff_sq)
        else:
          nb = []
        if [i for i in nb
            if model_atoms[i].element.strip().upper() in nonmetals
            and model_atoms[i].parent().parent().memory_id() != rgs.memory_id()]:
          continue
        # remove the amide H also for peptide-like residues of class 'other'
        # (UNK, ACB): add_n_terminal_hydrogens keeps an existing H in place of
        # H1, and H, H2, H3 do not match mod_NH3
        if (get_class(name=ag.resname) in
            ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']
            or (ag.get_atom('CA') and ag.get_atom('C'))):
          for ag_h in rgs.atom_groups():
            h = ag_h.get_atom('H')
            if h and _dictionary_h_on_n(ag_h, self.query_mon_lib_srv):
              ag_h.remove_atom(h)
        # TODO make the function below smart, so it
        # 1) knows when to add H1H2H3 or not
        # 2) renames H to H1 (so no need to remove it beforehand)
        split_ca = _split_ca_groups(rgs, n)
        for ag_ca in split_ca:
          c = ag_ca.get_atom('C')
          if c is None:
            c = ([x.get_atom('C') for x in rgs.atom_groups()
                  if x.altloc == '' and x.get_atom('C')] or [None])[0]
          if c is None: continue
          add_n_terminal_hydrogens_to_atom_group(ag_ca,
            n_ca_c=(n, ag_ca.get_atom('CA'), c))
        if not split_ca:
          rc = add_n_terminal_hydrogens_to_residue_group(rgs) # rc is always empty list?

  # ----------------------------------------------------------------------------

  def add_missing_H_atoms_at_bogus_position(self, exclude_water):
    '''
    Add missing H atoms at bogus positions to a copy of the model's
    pdb_hierarchy and return the copy.
    All H atoms are placed at center of coordinates + (0.5, 0.5, 0.5)
    The translation is necessary because sometimes the center of coordinates
    coincides with the position of a heavy atom.

    In one residue/entity, all newly placed H are superposed, they will be
    moved to their expected position later.
    '''
    # TODO temporary fix until v3 names are in mon lib
    alternative_names = [
      ('HA1', 'HA2', 'HA3'),
      ('HB1', 'HB2', 'HB3'),
      ('HG1', 'HG2', 'HG3'),
      ('HD1', 'HD2', 'HD3'),
      ('HE1', 'HE2', 'HE3'),
      ('HG11', 'HG12', 'HG13')
      ]
    # end TODO
    # With keep_existing_H, self.model is still the caller's model: H go into
    # a copy of its hierarchy, auto objects into a new list
    pdb_hierarchy = self.model.get_hierarchy().deep_copy()
    self.restraint_objects = self.model.get_restraint_objects()
    mon_lib_srv = self.mon_lib_srv
    if mon_lib_srv is None:
      mon_lib_srv = mmtbx.monomer_library.server.server()
    self.query_mon_lib_srv = mon_lib_srv
    load_restraint_objects_into_server(self.model, mon_lib_srv)
    #XXX This breaks for 1jxt, residue 2, TYR
    for m in pdb_hierarchy.models():
      model_atoms = m.atoms()
      model_sites = model_atoms.extract_xyz()
      for chain in m.chains():
        chain_rgs = chain.residue_groups()
        for i_rg, rg in enumerate(chain_rgs):
          n_atom_groups = len(rg.atom_groups())
          ags = rg.atom_groups()
          # blank + one altloc, same resname: one conformer. Missing H are
          # added once, to the blank group; per group, each got the other's H
          single_altloc = (n_atom_groups == 2 and
            sorted(ag.altloc == '' for ag in ags) == [False, True] and
            ags[0].resname == ags[1].resname)
          amide_ag = None
          if n_atom_groups > 2:
            amide_ag = _blank_amide_group(rg,
              chain_rgs[i_rg-1] if i_rg > 0 else None)
          add_amide_h = False
          for ag in ags:
            # blank + two or more altlocs: blank-parent H go into each altloc
            # (not the amide H if it depends on no split atom, amide_ag)
            if n_atom_groups > 2 and ag.altloc == '':
              continue
            if single_altloc and ag.altloc != '':
              continue
            #print list(ag.atoms().extract_name())
            if exclude_water:
              if(get_class(name=ag.resname) == "common_water"): continue
            actual = [a.name.strip().upper() for a in ag.atoms()]
            if single_altloc:
              actual += [a.name.strip().upper()
                         for other in ags if other.altloc != ''
                         for a in other.atoms()]
            #
            mlq, cif_object = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv, raise_sorry=False)
            if mlq is None:
              if ag.resname not in self.no_H_placed_mlq:
                self.no_H_placed_mlq.append(ag.resname)
              self.residue_keys_no_restraints.add(
                (m.id, chain.id, rg.resseq, rg.icode, ag.resname))
              continue

            if cif_object:
              self.restraint_objects = with_auto_restraint_object(
                self.restraint_objects, ag.resname, cif_object,
                self.auto_restraint_names)

            expected_h = []
            #expected_ha = []
            atom_dict = mlq.atom_dict()
            def _remove_atoms(atom_dict, names):
              remove=[]
              for k,v in atom_dict.items():
                if k in names:
                  remove.append(k)
              if remove:
                for r in remove:
                  del atom_dict[r]
              return atom_dict
            #
            # don't add polymer H atoms. Terminal H atoms added elsewhere
            #
            # test_for_peptide wants C=O; a reduced C-terminus (ASJ, 4dco) has
            # C-O, and its N is in the chain only if bonded to the previous C
            n = ag.get_atom('N')
            if (mlq.test_for_peptide(atom_dict) or
                (get_class(name=ag.resname) == 'modified_amino_acid'
                 and 'CA' in atom_dict and n is not None
                 and _bonded_to_other_residue(n, rg, model_atoms, model_sites))):
              atom_dict = _remove_atoms(
                atom_dict, _terminal_h(mlq, atom_dict))
            # HO3' used to go here for every RNA/DNA residue, because in a
            # polymer O3' carries the next phosphate. It is placed now and
            # exclude_H_on_esterified_O drops it where the O really is
            # esterified: a free nucleotide (AMP, IMP, U5P) and a 3' end keep
            # their hydroxyl H.
            for k, v in six.iteritems(atom_dict):
              if(v.type_symbol=="H"):
                expected_h.append(k)
              #else:
              #  expected_ha.append(k)
            if (self.na_terminal_h and
                get_class(name=ag.resname) == 'common_rna_dna'):
              expected_h.extend(_na_terminal_h(ag, rg, expected_h,
                                               model_atoms, model_sites))
            #print('expected H', expected_h)
            #
            # TODO start
            # temporary fix until v3 names are in mon lib
            # Not for modified amino acids: their dictionaries follow the CCD
            # names and pdb_interpretation has no v3 mapping for them, so the
            # renamed H got no restraint and was dropped (216 geostd codes,
            # 2hi2 MEA CB, 9eor HT7 CA).
            if (get_class(name=ag.resname) in
                ['common_amino_acid', 'd_amino_acid']):
              for altname in alternative_names:
                if (altname[0] in expected_h and altname[1] in expected_h):
                  if (atom_dict[altname[0]].type_energy == 'HCH2' and
                      atom_dict[altname[1]].type_energy == 'HCH2'):
                    expected_h.append(altname[2])
                    expected_h.remove(altname[0])
                    #print('renamed %s to %s' % (altname[0], altname[2]))
            # TODO end
            missing_h = list(set(expected_h).difference(set(actual)))
            if (amide_ag is not None and 'H' in missing_h and
                any(set([b.atom_id_1, b.atom_id_2]) == set(['N', 'H'])
                    for b in mlq.bond_list)):
              missing_h.remove('H')
              add_amide_h = True
            if 0: print(ag.resname, missing_h)
            #new_xyz = ag.atoms().extract_xyz().mean()
            new_xyz = flex.double(ag.atoms().extract_xyz().mean()) + \
              flex.double([0.5,0.5,0.5])
            new_xyz = tuple(new_xyz)

            hetero = ag.atoms()[0].hetero
            segid = ag.atoms()[0].segid

            if missing_h:
              self.residue_keys_new_h.add(
                (m.id, chain.id, rg.resseq, rg.icode))
            for mh in missing_h:
              if len(mh) < 4: mh = (' ' + mh).ljust(4)
              a = (iotbx.pdb.hierarchy.atom()
                .set_name(new_name=mh)
                .set_element(new_element="H")
                .set_xyz(new_xyz=new_xyz)
                .set_hetero(new_hetero=hetero)
                .set_segid(new_segid=segid))

              ag.append_atom(a)
          if add_amide_h:
            new_xyz = tuple(flex.double(amide_ag.atoms().extract_xyz().mean())
              + flex.double([0.5,0.5,0.5]))
            amide_ag.append_atom(iotbx.pdb.hierarchy.atom()
              .set_name(new_name=' H  ')
              .set_element(new_element="H")
              .set_xyz(new_xyz=new_xyz)
              .set_hetero(new_hetero=amide_ag.atoms()[0].hetero)
              .set_segid(new_segid=amide_ag.atoms()[0].segid))
            self.residue_keys_new_h.add((m.id, chain.id, rg.resseq, rg.icode))
    return pdb_hierarchy

# ------------------------------------------------------------------------------
#
#  def validate_electrons(self):
#    from elbow.quantum import electrons
#    atom_valences = electrons.electron_distribution(
#      self.model.get_hierarchy(), # needs to be altloc free
#      self.model.get_restraints_manager().geometry,
#      verbose=False,
#    )
#    atom_valences.validate(ignore_water=True, raise_if_error=False)
#    self.charged_atoms = atom_valences.get_charged_atoms()

# ------------------------------------------------------------------------------

  def h_on_tertiary_amide_n(self, bonds, atoms, elements):
    """Return {h_iseq: n_iseq} for H/D atoms bonded to an amino-acid backbone N
    that is a tertiary amide. In the H's conformer the N already has three heavy
    neighbours (CA, an N-substituent such as a methyl or the proline-type ring
    carbon, and the peptide C), so it cannot carry an H. Only heavy neighbours
    that share the H's conformer (blank or same altloc) are counted, so a
    backbone split at CA is not mistaken for three substituents. The caller
    passes an already-built {iseq: [neighbor iseqs]} map to avoid re-walking the
    bond proxies.
    """
    result = {}
    for n_iseq, neighbors in bonds.items():
      atom_n = atoms[n_iseq]
      if atom_n.name.strip() != 'N': continue
      if get_class(name=atom_n.parent().resname) not in (
        'common_amino_acid', 'modified_amino_acid', 'd_amino_acid'): continue
      neighbors = set(neighbors)
      for k in neighbors:
        if elements[k] not in ["H", "D"]: continue
        h_altloc = atoms[k].parent().altloc
        n_heavy = sum(1 for m in neighbors if elements[m] not in ["H", "D"]
          and atoms[m].parent().altloc in ('', h_altloc))
        if n_heavy >= 3:
          result[k] = n_iseq
    return result

  def exclude_H_on_links(self, verbose=False):
    """Remove H atoms bound to heavy atoms that form a link

    An exception are HD1 and HE2 of HIS. The mover functionality in reduce will
    take care of those.

    TODO: Could restraints manager have a list of links with relevant information?
          Then we don't have to loop through all proxies here.
    """
    from mmtbx.ligands.chemistry import get_valences
    origin_ids = linking_class()
    grm = self.model.get_restraints_manager()
    bond_proxies_simple, asu = grm.geometry.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    atoms = self.model.get_atoms()
    elements = self.model.get_hierarchy().atoms().extract_element()
    exclusion_iseqs = list()
    exclusion_dict = dict()
    link_partners = dict()
    # Metal coordination is decided elsewhere (A4): thiols by
    # exclude_H_on_metal_thiolates(), other H on metal-bound atoms stay. It
    # neither marks a link nor counts as a bond here.
    metal = origin_ids['metal coordination']
    metal_partners = dict()
    all_proxies = []
    for proxy in list(bond_proxies_simple) + list(asu):
      if proxy.origin_id == metal:
        if isinstance(proxy, ext.bond_simple_proxy): i, j = proxy.i_seqs
        else:                                        i, j = proxy.i_seq, proxy.j_seq
        metal_partners.setdefault(i, set()).add(j)
        metal_partners.setdefault(j, set()).add(i)
        continue
      all_proxies.append(proxy)
    # Loop through bond proxies to find links (i.e. proxies with origin_id != 0)
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      if proxy.origin_id != 0:
        exclusion_iseqs.extend([i,j])
        exclusion_dict[i] = proxy.origin_id
        exclusion_dict[j] = proxy.origin_id
        if isinstance(proxy, ext.bond_simple_proxy): # asu partner needs rt_mx
          link_partners.setdefault(i, []).append(j)
          link_partners.setdefault(j, []).append(i)
    sel_remove = flex.size_t()

    # Find H atoms bound to linked atoms
    removed_dict = {}
    parent_dict = {}
    bonds = {}
    bond_lengths = {}
    ideal_distance = {}
    # Exception for HIS HD1 and HE2, unless the ring N carries a covalent link:
    # then the mover has no protonation choice left (4us9 HIS 752 NE2 - 3PL C).
    # Metal coordination stays with the optimizer's ion lock-down.
    def _his_exception(i_h, i_parent):
      if atoms[i_h].parent().resname != 'HIS': return False
      if atoms[i_h].name.strip() not in ['HD1','DD1', 'HE2', 'DE2']: return False
      return (i_parent not in exclusion_iseqs or
              exclusion_dict.get(i_parent) == origin_ids['metal coordination'])
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      bonds.setdefault(i,[])
      bonds[i].append(j)
      bonds.setdefault(j,[])
      bonds[j].append(i)
      ideal_distance[frozenset((i, j))] = proxy.distance_ideal
      if _his_exception(i, j): continue
      if _his_exception(j, i): continue
      if(elements[i] in ["H","D"] and j in exclusion_iseqs):
        if i not in sel_remove:
          sel_remove.append(i)
          bond_lengths[i] = proxy.distance_ideal
          removed_dict[i] = exclusion_dict[j]
          parent_dict[i]=j
      if(elements[j] in ["H","D"] and i in exclusion_iseqs):
        if j not in sel_remove:
          sel_remove.append(j)
          bond_lengths[j] = proxy.distance_ideal
          removed_dict[j] = exclusion_dict[i]
          parent_dict[j]=i
    # A peptide-bonded N-substituted (e.g. N-methylated) amino acid is a
    # tertiary amide: its backbone N has three heavy neighbours (CA, the
    # N-substituent, the previous residue's C) and cannot carry a backbone H.
    # The peptide bond is plain covalent geometry (origin_id 0), so it is not
    # flagged as a link above; flag the leftover H here so the valence check
    # below removes it. (Reuses the bonds map already built above.)
    for k, n_iseq in self.h_on_tertiary_amide_n(bonds, atoms, elements).items():
      if k not in sel_remove:
        sel_remove.append(k)
        removed_dict[k] = 'a tertiary amide'
        parent_dict[k] = n_iseq
    # remove H atoms NOT to remove - double negative!
    #verbose=True
    if verbose:
      print('removed_dict',removed_dict)
      for i_seq in sel_remove:
        print('remove?',atoms[i_seq].quote())
    remove_from_sel_remove=[]
    # Parent with >= 2 other heavy neighbours: the H on the link partner's site
    # goes first (8oji VPH H4 sits 0.86 A from S1). Terminal parents (NH3,
    # CH3) have conventional names: old order, lowest number survives.
    def _drop_order(i_seq):
      parent = parent_dict[i_seq]
      partners = link_partners.get(parent, [])
      heavy = [n for n in set(bonds.get(parent, []))
               if elements[n].strip() not in ('H', 'D') and n not in partners]
      if not partners or len(heavy) < 2: return float('inf')
      return min(atoms[i_seq].distance(atoms[p]) for p in partners)
    # Count bond orders, not neighbours: an acyl-enzyme ester carbon (CA, =O,
    # link to SER OG) has no room for H (4jxg). Links count 1.
    mon_lib_srv = self.model.get_mon_lib_srv()
    order_cache = {}
    def _residue(atom):
      rg = atom.parent().parent()
      return (rg.parent().id, rg.resid(), atom.parent().resname.strip())
    def _bond_order(i, j):
      if _residue(atoms[i]) != _residue(atoms[j]): return 1
      orders = _bond_orders(atoms[i].parent().resname.strip(), mon_lib_srv,
                            order_cache)
      order = orders.get(frozenset((atoms[i].name.strip(), atoms[j].name.strip())), 1)
      # A multiple bond the link reaction consumed is single in the model (5p9j
      # CAA=CAD at 1.60 A). 0.1 A is half the single-double gap.
      if (order > 1 and atoms[i].distance(atoms[j]) >
          ideal_distance.get(frozenset((i, j)), 0) + 0.1): return 1
      return order
    def _altloc(i_seq):
      return atoms[i_seq].parent().altloc
    # Count one conformer: H copies in A and B are not simultaneous bonds of a
    # shared parent. A blank-altloc H is checked against the fullest conformer.
    def _number_of_bonds(j_seq, i_h):
      nbs = bonds[j_seq]
      alts = [_altloc(i_h)] if _altloc(i_h) else (
        sorted(set(_altloc(n) for n in nbs) - set([''])) or [''])
      return max(sum(_bond_order(j_seq, n) for n in nbs
                     if _altloc(n) in ('', alt)) for alt in alts)
    missing_neighbor = 'a missing heavy neighbour'
    for i_seq in sorted(reversed(list(sel_remove)), key=_drop_order):
      j_seq=parent_dict[i_seq]
      # need to add the use of atomic charge
      valences=get_valences(elements[j_seq])
      number_of_bonds=_number_of_bonds(j_seq, i_seq)
      if number_of_bonds in valences:
        # remove this H from delection
        remove_from_sel_remove.append(i_seq) # ??
        del removed_dict[i_seq]
      else:
        # too few bonds: a heavy neighbour is missing, not the link's fault
        # (3evj NAG A 2 C1 has no C2)
        if number_of_bonds < min(valences):
          removed_dict[i_seq] = missing_neighbor
        bonds[j_seq].remove(i_seq)
        bonds[i_seq].remove(j_seq)

    #print(remove_from_sel_remove)
    fsc0=grm.geometry.shell_sym_tables[0].full_simple_connectivity()
    fsc1=grm.geometry.shell_sym_tables[1].full_simple_connectivity()
    #fsc2=grm.geometry.shell_sym_tables[2].full_simple_connectivity()
    for _i in remove_from_sel_remove:
      parent = fsc0[_i][0]
      first_neighbors = fsc1[_i]
      # same conformer as the H
      fn_filtered = [item for item in first_neighbors if item not in sel_remove
                     and _altloc(item) in ('', _altloc(_i))
                     and item not in metal_partners.get(parent, ())]
      def _sibling(k):
        return fsc0[k][0] == parent and _altloc(k) == _altloc(_i)
      n_kept = sum(1 for k in remove_from_sel_remove if _sibling(k))
      coordp = matrix.col(atoms[parent].xyz)
      # Two H kept on an sp3 atom with two heavy neighbours: a methyl that
      # became a CH2 via a link (e.g. B0I CB1/CB2 in 6b17, linked to a thioether
      # S). Place both H tetrahedrally, symmetric across the plane of the two
      # heavy neighbours, so the real C-C and C-S(link) directions are respected
      # (not eclipsed). The single-H branches below cannot handle this: they
      # would place both H at the same in-plane bisector site.
      if n_kept == 2 and len(fn_filtered) == 2:
        siblings = sorted(k for k in remove_from_sel_remove if _sibling(k))
        u1 = (matrix.col(atoms[fn_filtered[0]].xyz) - coordp).normalize()
        u2 = (matrix.col(atoms[fn_filtered[1]].xyz) - coordp).normalize()
        anti = -(u1 + u2).normalize()       # bisects H-C-H, points away from neighbours
        perp = (u1.cross(u2)).normalize()    # normal to the heavy-atom plane
        beta = math.radians(54.735)          # half of the tetrahedral angle 109.47
        s = 1.0 if siblings.index(_i) == 0 else -1.0
        d = anti*math.cos(beta) + perp*(s*math.sin(beta))
        atoms[_i].xyz = coordp + d.normalize()*bond_lengths[_i]
        continue
      if n_kept > 1:
        continue  # other multi-H cases: leave the riding-idealized positions
      # --- a single H is kept on `parent` ---
      # A lone H kept on a linked backbone N is a peptide-like secondary amide
      # H, not an N-terminus: rename a leftover propeller name (H1/H2/H3,
      # D1/D2/D3) to the plain 'H'/'D' (e.g. cyclic peptide linked through its
      # N-terminal N, as in 3njw).
      parent_atom = atoms[parent]
      if (parent_atom.name.strip() == 'N' and
          get_class(name=parent_atom.parent().resname) in
            ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid'] and
          atoms[_i].name.strip() in ['H1', 'H2', 'H3', 'D1', 'D2', 'D3']):
        old_name = atoms[_i].name
        atoms[_i].name = ' %s  ' % atoms[_i].element.strip()
        self._log_rename(parent_atom, atoms[_i], old_name)
      # improve geometry of the single kept H
      # TODO make sure all atoms are in same conformer
      # TODO check that neighbor atoms are all non H
      if len(fn_filtered) == 3:
        #print('tetrahedral geometry')
        coord1 = matrix.col(atoms[fn_filtered[0]].xyz)
        coord2 = matrix.col(atoms[fn_filtered[1]].xyz)
        coord3 = matrix.col(atoms[fn_filtered[2]].xyz)
        orth = (coord2-coord1).cross(coord3-coord1).normalize()
        vol = orth.dot(coordp-coord1)
        if vol > 0: orth = -orth
        atoms[_i].xyz = coordp - orth*bond_lengths[_i]
      if len(fn_filtered) == 2:
        #print('flat geometry')
        coord1 = matrix.col(atoms[fn_filtered[0]].xyz)
        coord2 = matrix.col(atoms[fn_filtered[1]].xyz)
        half = ((coord1 - coordp).normalize() + (coord2 - coordp).normalize())
        atoms[_i].xyz = coordp-half.normalize()*bond_lengths[_i]

    if remove_from_sel_remove:
      # direct xyz edits above: keep the X-ray structure in step
      self.model.set_sites_cart(self.model.get_sites_cart())
      sel_remove=list(sel_remove)
      for r in remove_from_sel_remove:
        sel_remove.remove(r)
        if verbose: print('keep',atoms[r].quote())
      sel_remove=flex.size_t(sel_remove)
    #
    def _removed_label(v):
      return v if isinstance(v, str) else origin_ids.get_origin_key(v)
    sl_removed = []
    for atom in self.model.get_hierarchy().atoms().select(sel_remove):
      label = atom.id_str().replace('pdb=','').replace('"','')
      if removed_dict[atom.i_seq] == missing_neighbor:
        rg = atom.parent().parent()
        resid = '%s %s %s' % (atom.parent().resname.strip(),
          rg.parent().id.strip(), rg.resseq.strip())
        self.site_labels_missing_neighbor.append(label)
        if resid not in self.residues_missing_neighbor:
          self.residues_missing_neighbor.append(resid)
        continue
      sl_removed.append((label, _removed_label(removed_dict[atom.i_seq])))
    #
    if sel_remove:
      self.model = self.model.select(~flex.bool(self.model.size(), sel_remove))
    self.sl_removed = sl_removed
    self.exclusion_iseqs = exclusion_iseqs


# ------------------------------------------------------------------------------

  def unknown_atoms_text(self, verbose=False):
    '''
    The report of residues_unknown_atoms (show(), reduce2's description).
    '''
    return atoms_not_in_dictionary_text(self.residues_unknown_atoms,
      len(self.site_labels_unknown_atoms), 'H atoms on them were not placed',
      self.site_labels_unknown_atoms if verbose else ())

# ------------------------------------------------------------------------------

  def metal_thiol_text(self):
    '''
    The report of thiols bound to metals (show(), reduce2's description);
    '' if there is none. supplied_h_hint, if set, says how to drop the
    supplied H.
    '''
    return metal_thiol_text(self.site_labels_thiolate, self.thiol_h_unresolved,
      self.thiol_h_supplied_kept, self.supplied_h_hint)

# ------------------------------------------------------------------------------

  def metal_nh_text(self):
    '''
    The report of N-H at metals (show(), reduce2's description); '' if
    there is none.
    '''
    return metal_nh_text(self.site_labels_metal_nh, self.metal_nh_unresolved,
      self.metal_nh_supplied_kept, self.metal_nh_reported, self.supplied_h_hint)

# ------------------------------------------------------------------------------

  def metal_amine_text(self):
    '''
    The report of metal_bound_amines_as_nh2 (show(), reduce2's description);
    '' if there is none.
    '''
    return metal_amine_text(self.metal_amine_converted_text,
      self.metal_amine_unresolved, self.metal_amine_supplied,
      self.metal_amine_chain_starts, self.metal_amine_supplied_nh2)

# ------------------------------------------------------------------------------

  def show(self, log, verbose=False):
    '''
    Informative output. With verbose=True the per-atom lists of H that were not
    placed on incomplete side chains are printed in full; by default only a
    count + residue summary is shown.
    '''
    if log is None: log = sys.stdout
    #
    if (not self.keep_existing_H and self.n_H_initial):
      msg = 'Number of hydrogen atoms trimmed from input model: %s \n'
      print(msg % self.n_H_initial, file=log)
    #
    msg = 'Number of hydrogen atoms added to the input model: %s \n'
    print(msg % self.n_H_final, file=log)
    #
    if self.no_H_placed_mlq:
      msg = '''
No H atoms were placed on the following residues because no restraints
were found:'''
      print(msg, file=log)
      for resname in self.no_H_placed_mlq:
        print(resname, file=log)
    #
    if self.residues_no_restraints:
      msg = ('\nNo H atoms were placed on the following residues because the '
        'processed model has no restraints for them:')
      print(msg, file=log)
      for resid in self.residues_no_restraints:
        print(resid, file=log)
    #
    if self.site_labels_disulfides:
      msg = '''
The following cysteine HG atoms were not placed because the sulfur atom is
involved in a disulfide bond'''
      print(msg, file=log)
      for label in self.site_labels_disulfides:
        print(label, file=log)
    #
    if self.site_labels_tertiary_amide:
      msg = '''
The following backbone H atoms were not placed because the residue's N atom
is a tertiary amide (already bonded to three heavy atoms, e.g. an internal
proline-type or N-substituted residue)'''
      print(msg, file=log)
      for label in self.site_labels_tertiary_amide:
        print(label, file=log)
    #
    if self.residues_unknown_atoms:
      print('', file=log)
      print(self.unknown_atoms_text(verbose), file=log)
    #
    for text in (self.metal_thiol_text(), self.metal_nh_text(),
                 self.metal_amine_text()):
      if text:
        print('', file=log)
        print(text, file=log)
    #
    if self.site_labels_missing_neighbor:
      n_atoms = len(self.site_labels_missing_neighbor)
      residues = self.residues_missing_neighbor
      print('', file=log)
      msg = ('%d H atoms were not placed because a neighbouring heavy atom is '
        'missing\n(incomplete side chains) - %d residues:')
      print(msg % (n_atoms, len(residues)), file=log)
      cap = 12
      tokens = list(residues[:cap])
      more = '(+%d more)' % (len(residues) - cap) if len(residues) > cap else None
      line = ''
      for t in tokens:
        add = (', ' if line else '') + t
        if line and len(' ' + line + add) > 78:
          print(' ' + line, file=log)
          line = t
        else:
          line += add
      if more is not None:
        add = (' ' if line else '') + more
        if line and len(' ' + line + add) > 78:
          print(' ' + line, file=log)
          line = more
        else:
          line += add
      if line:
        print(' ' + line, file=log)
      if verbose:
        for label in self.site_labels_missing_neighbor:
          print(label, file=log)
    #
    if self.site_labels_riding_refused:
      msg = '''
The following H atoms were not placed because riding H could not parameterize
them (their parent atom has all its heavy neighbours)'''
      print(msg, file=log)
      for label in self.site_labels_riding_refused:
        print(label, file=log)
    #
    if self.site_labels_no_para:
      msg = '''
The following H atoms were not placed because they could not be parameterized
(not enough restraints information)'''
      print(msg, file=log)
      for label in self.site_labels_no_para:
        print(label, file=log)
#    if self.charged_atoms:
#      msg = '''
#The following heavy atom have an unusual electron count. This could be because
#heavy atoms or H atoms are missing.'''
#      print(msg, file=log)
#      for item in self.charged_atoms:
#        idstr = item[0].id_str().replace('pdb=','').replace('"','')
#        if 'HOH' in idstr: continue
#        print(idstr, item[1])

    if self.sl_removed:
      print()
      msg = '''Atom %s was not placed because it is involved in %s'''
      for item in self.sl_removed:
        print(msg % (item[0], item[1]), file=log)

# ------------------------------------------------------------------------------

  def get_model(self):
    return self.model

# ------------------------------------------------------------------------------

  def get_counts(self):
    return group_args(
      number_h_final  = self.n_H_final,
      no_H_placed_mlq = self.no_H_placed_mlq,
      residues_no_restraints = self.residues_no_restraints,
      site_labels_disulfides = self.site_labels_disulfides,
      site_labels_no_para = self.site_labels_no_para,
      site_labels_tertiary_amide = self.site_labels_tertiary_amide,
      site_labels_missing_neighbor = self.site_labels_missing_neighbor,
      residues_missing_neighbor = self.residues_missing_neighbor,
      residues_unknown_atoms = self.residues_unknown_atoms,
      site_labels_unknown_atoms = self.site_labels_unknown_atoms,
      site_labels_thiolate = self.site_labels_thiolate,
      thiol_h_unresolved = self.thiol_h_unresolved,
      thiol_h_supplied_kept = self.thiol_h_supplied_kept,
      site_labels_metal_nh = self.site_labels_metal_nh,
      metal_nh_unresolved = self.metal_nh_unresolved,
      metal_nh_supplied_kept = self.metal_nh_supplied_kept,
      metal_nh_reported = self.metal_nh_reported,
      metal_amine_converted = self.metal_amine_converted_text,
      metal_amine_unresolved = self.metal_amine_unresolved,
      metal_amine_supplied = self.metal_amine_supplied,
      metal_amine_supplied_nh2 = self.metal_amine_supplied_nh2,
      metal_amine_chain_starts = self.metal_amine_chain_starts,
      site_labels_riding_refused = self.site_labels_riding_refused)

# ------------------------------------------------------------------------------

  def get_times(self):
    return group_args(
      time_rebox_model        = self.time_rebox_model,
      time_remove_element_X   = self.time_remove_element_X,
      time_add_missing_H      = self.time_add_missing_H,
      time_terminal_propeller = self.time_terminal_propeller,
      time_make_grm           = self.time_make_grm,
      time_remove_isolated    = self.time_remove_isolated,
      time_riding_manager     = self.time_riding_manager,
      time_remove_H_nopara    = self.time_remove_H_nopara,
      time_reset              = self.time_reset,
      time_idealize           = self.time_idealize,
      time_remove_H_on_links  = self.time_remove_H_on_links)

# ------------------------------------------------------------------------------

  def print_times(self):
    print('Detailed timings:')
    print("Rebox model:", self.time_rebox_model)
    print('Remove element X:', self.time_remove_element_X)
    print("Add missing H at bogus position:", self.time_add_missing_H)
    print('Add N-terminal propeller:', self.time_terminal_propeller)
    print("Get new model obj and grm:", self.time_make_grm )
    print("Remove isolated H:", self.time_remove_isolated)
    print("Setup Riding manager:", self.time_riding_manager)
    print("Remove H that were not parameterized:", self.time_remove_H_nopara)
    print("Reset adp, occ:", self.time_reset)
    print("idealize H positions:", self.time_idealize)
    print("Remove H on links:", self.time_remove_H_on_links)
    print()

# ==============================================================================
