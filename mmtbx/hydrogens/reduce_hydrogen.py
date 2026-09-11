from __future__ import absolute_import, division, print_function
import six
import sys, time
from libtbx.utils import Sorry
import mmtbx.model
import iotbx.pdb
import boost_adaptbx.boost.python as bp
from libtbx.utils import null_out
from libtbx import group_args
from scitbx import matrix
from cctbx.array_family import flex
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_residue_group
from cctbx.geometry_restraints.linking_class import linking_class
#
from cctbx.maptbx.box import shift_and_box_model
import math
from cctbx import crystal

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
  for i, a in enumerate(cc_cif.get('_chem_comp_atom',[])):
    lookup[a.atom_id]=i
    lookup[i]=a.atom_id
    if a.type_symbol in ['H', 'D']:
      hs.append(a.atom_id)
      hsi.append(i)
    comp_comp_id.atom_list.append(cif_types.chem_comp_atom(
      atom_id=a.atom_id,
      type_symbol=a.type_symbol,
      # type_energy=a.type_energy,
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
    comp_comp_id.bond_list.append(cif_types.chem_comp_bond(
      atom_id_1=b.atom_id_1,
      atom_id_2=b.atom_id_2,
      type=b.value_order,
      value_dist='%0.3f' % (bl*.9),
      value_dist_esd=".1"))

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

def bonds_in_restraints(atom, exclude_hydrogens=False):
  from mmtbx.chemical_components import get_cif_dictionary
  cc_cif = get_cif_dictionary(atom.parent().resname)
  rc=[]
  for b in cc_cif.get('_chem_comp_bond',[]):
    if b.atom_id_1.strip()==atom.name.strip():
      if exclude_hydrogens:
        a=atom_in_restraints(b.atom_id_2, cc_cif)
        if a.type_symbol in ['H', 'D']: continue
      rc.append(b.atom_id_2)
    if b.atom_id_2.strip()==atom.name.strip():
      if exclude_hydrogens:
        a=atom_in_restraints(b.atom_id_1, cc_cif)
        if a.type_symbol in ['H', 'D']: continue
      rc.append(b.atom_id_1)
  return rc

# ==============================================================================

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

def workaround_002(model, selection):
  h = model.get_hierarchy()
  atoms = h.atoms()
  for m in h.models():
    for c in m.chains():
      for con in c.conformers():
        for r in con.residues():
          if not get_class(name=r.resname) == "common_water": continue
          three = r.atoms().extract_i_seq()
          #assert three in selection
          ij=[]
          for atom in r.atoms():
            e = atom.element.strip().upper()
            if e == "O": X = atom.xyz
            else:        ij.append(atom.i_seq)
          if len(ij)!=2: continue
          p1, p2 = find_H1_H2(X=X, d=0.85, angle_deg=103.91)
          atoms[ij[0]].xyz = p1
          atoms[ij[1]].xyz = p2

def _chiral_volume(c, a, b, h):
  c, a, b, h = [matrix.col(x) for x in (c, a, b, h)]
  return (a - c).dot((b - c).cross(h - c))

def _ch2_groups(elements, bond_pairs, sites):
  '''Heavy atoms with exactly two H and two heavy neighbours, all sites known.'''
  neighbors = {}
  for i, j in bond_pairs:
    neighbors.setdefault(i, []).append(j)
    neighbors.setdefault(j, []).append(i)
  result = []
  for p, ns in neighbors.items():
    if elements.get(p) in ('H', 'D'): continue
    hs = tuple(sorted(n for n in ns if elements.get(n) in ('H', 'D')))
    hv = tuple(sorted(n for n in ns if elements.get(n) not in ('H', 'D')))
    if len(hs) == 2 and len(hv) == 2 and all(n in sites for n in (p,)+hv+hs):
      result.append((p, hv, hs, sites))
  return result

def _dictionary_sites(source_info):
  '''Ideal sites from a restraint file, if it has coordinates (geostd ligands).'''
  if not source_info or not source_info.startswith("file:"): return None
  import iotbx.cif
  try: cif_model = iotbx.cif.reader(file_path=source_info[5:].strip()).model()
  except Exception: return None
  for block in cif_model.values():
    if "_chem_comp_atom.atom_id" in block and "_chem_comp_atom.x" in block:
      xyz = zip(*[block["_chem_comp_atom.%s" % k] for k in "xyz"])
      return dict((name.strip('"'), tuple(float(v) for v in t))
        for name, t in zip(block["_chem_comp_atom.atom_id"], xyz)
        if "?" not in t and "." not in t)
  return None

def _ch2_references(resname, mon_lib_srv, cache):
  '''
  CH2 centres with ideal sites. CCD first: it defines PDB names, and geostd
  amino acids carry no coordinates. Then the restraint dictionary, for ligands
  whose H names differ from the CCD (VPH). geostd MAN names H61/H62 opposite
  to the CCD.
  '''
  if resname in cache: return cache[resname]
  from mmtbx.chemical_components import get_cif_dictionary
  result = []
  try: cc_cif = get_cif_dictionary(resname)
  except Exception: cc_cif = None
  if cc_cif:
    sites, elements = {}, {}
    for a in cc_cif.get('_chem_comp_atom', []):
      elements[a.atom_id] = a.type_symbol.strip().upper()
      t = [getattr(a, "pdbx_model_Cartn_%s_ideal" % k, "?") for k in "xyz"]
      if "?" not in t and "." not in t:
        sites[a.atom_id] = tuple(float(v) for v in t)
    result += _ch2_groups(elements, [(b.atom_id_1, b.atom_id_2)
      for b in cc_cif.get('_chem_comp_bond', [])], sites)
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  if cc is not None:
    sites = _dictionary_sites(cc.source_info)
    if sites:
      elements = dict((a.atom_id, a.type_symbol.strip().upper())
        for a in cc.atom_list)
      result += _ch2_groups(elements, [(b.atom_id_1, b.atom_id_2)
        for b in cc.bond_list], sites)
  cache[resname] = result
  return result

def name_prochiral_h(hierarchy, mon_lib_srv):
  '''
  Riding places the two H of a CH2 in processing order, so about half get each
  other's name (1akg: 16 of 38 vs CCD). Swap names to match ideal chirality.
  '''
  cache, done = {}, set()
  for m in hierarchy.models():
    for c in m.chains():
      for conformer in c.conformers():
        for r in conformer.residues():
          groups = _ch2_references(r.resname.strip(), mon_lib_srv, cache)
          if not groups: continue
          atoms = dict((a.name.strip(), a) for a in r.atoms())
          for p, hv, hs, sites in groups:
            if not all(n in atoms for n in (p,)+hv+hs): continue
            if atoms[p].i_seq in done: continue # first source wins
            if any(atoms[p].distance(atoms[n]) > 2.0 for n in hv+hs): continue
            done.add(atoms[p].i_seq)
            h1, h2 = atoms[hs[0]], atoms[hs[1]]
            v_ideal = _chiral_volume(sites[p], sites[hv[0]], sites[hv[1]],
                                     sites[hs[0]])
            v_model = _chiral_volume(atoms[p].xyz, atoms[hv[0]].xyz,
                                     atoms[hv[1]].xyz, h1.xyz)
            if abs(v_ideal) < 0.5 or abs(v_model) < 0.5: continue # flat
            if (v_ideal > 0) != (v_model > 0):
              h1.name, h2.name = h2.name, h1.name

def _bond_orders(resname, mon_lib_srv, cache):
  '''
  Double (2) and triple (3) bonds by atom-name pair; anything else counts 1.
  Restraint dictionary first, CCD overrides (PDB names).
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
               print_time            = False):
    self.model                 = model
    self.use_neutron_distances = use_neutron_distances
    self.n_terminal_charge     = n_terminal_charge
    self.adp_scale             = adp_scale
    self.exclude_water         = exclude_water
    self.stop_for_unknowns     = stop_for_unknowns
    self.keep_existing_H       = keep_existing_H
    self.validate_e            = validate_e
    self.print_time            = print_time
    #
    self.no_H_placed_mlq        = list()
    self.site_labels_disulfides = list()
    self.site_labels_no_para    = list()
    self.site_labels_tertiary_amide = list()
    self.site_labels_missing_neighbor = list()
    self.residues_missing_neighbor    = list()
    # Names of restraint dictionaries auto-generated for unknown ligands during
    # placement; these are throwaway (purpose-built for H placement) and are
    # removed from the model before returning so they don't leak into
    # downstream geometry validation.
    self.auto_restraint_names   = list()
    #self.charged_atoms          = list()
    self.sl_removed             = list()
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
    if not self.auto_restraint_names: return
    ro = self.model.get_restraint_objects()
    if not ro: return
    self.model._restraint_objects = [
      (name, obj) for name, obj in ro
      if name not in self.auto_restraint_names]

# ------------------------------------------------------------------------------

  def run(self):
    '''
    Function that places H atoms
    '''

    # Create symmetry if necessary
    # ------------------------------
    model_has_bogus_cs = False
    t0 = time.time()
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)
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

    # Add missing H atoms and place them at bogus position
    # ----------------------------------------------------
    t0 = time.time()
    pdb_hierarchy = self.add_missing_H_atoms_at_bogus_position(
      exclude_water = self.exclude_water)
    self.time_add_missing_H = round(time.time()-t0, 2)

    # Place N-terminal propeller hydrogens
    # TODO double check N-terminal position for PRO residues
    # ------------------------------------
    if self.n_terminal_charge in ['residue_one', 'first_in_chain']:
      t0 = time.time()
      self.place_n_terminal_propeller(pdb_hierarchy = pdb_hierarchy)
      self.time_terminal_propeller = round(time.time()-t0, 2)

    pdb_hierarchy.sort_atoms_in_place()
    pdb_hierarchy.atoms().reset_serial()
    # Make new model obj and get restraints manager
    # ---------------------------------------------
    p = get_reduce_pdb_interpretation_params(self.use_neutron_distances)
    ro = self.model.get_restraint_objects()
    t0 = time.time()
    self.model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = pdb_hierarchy,
      stop_for_unknowns = self.stop_for_unknowns,
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = ro,
      log               = null_out())
    self.model.process(pdb_interpretation_params=p,
                       make_restraints=True,
                       # retain_zero_dihedrals=True,
                       )

    self.time_make_grm = round(time.time()-t0, 2)
    # Drop the throwaway auto-generated ligand restraints now that the grm has
    # been built. They were only needed to place H on unknown ligands; keeping
    # them on the model would leak idealized placement-only geometry (0.9x
    # bonds, esd=1/period=1 torsions) into any downstream re-interpretation.
    self.remove_auto_restraint_objects()
    # Return if no H have been placed
    sel_h = self.model.get_hd_selection()
    if sel_h.count(True) == 0: return

    # Remove isolated H atoms
    # -----------------------
    # (when heavy atom is missing, H needs not to be placed)
    t0 = time.time()
    sel_isolated = self.model.isolated_atoms_selection()
    sel_lone_H = sel_h & sel_isolated
    # As h_parameterization will not include these, they can be removed in the
    # next step; for book-keeping it is useful to keep track of lone H as a
    # selection
    #if not sel_lone_H.all_eq(False):
    #  self.model = self.model.select(~sel_lone_H)
    self.time_remove_isolated = round(time.time()-t0, 2)

    sel_h = self.model.get_hd_selection()

    # Setup riding H manager
    # ----------------------
    t0 = time.time()
    self.model.setup_riding_h_manager(use_ideal_dihedral = True)
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
        selection = water_selection.iselection())
      water_selection = self.model.solvent_selection()
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
      seen_residues = set()
      for atom in self.model.get_hierarchy().atoms().select(
          sel_h_not_in_para_but_not_lone):
        label = atom.id_str().replace('pdb=','').replace('"','')
        if atom.i_seq in tertiary:
          self.site_labels_tertiary_amide.append(label)
          continue
        parents = _heavy_neighbors(atom.i_seq)
        # parent atom missing an expected heavy neighbour -> incomplete residue
        if len(parents) == 1 and len(_heavy_neighbors(parents[0])) < 2:
          self.site_labels_missing_neighbor.append(label)
          rg = atom.parent().parent()
          resid = '%s %s %s' % (atom.parent().resname.strip(),
            rg.parent().id.strip(), rg.resseq.strip())
          if resid not in seen_residues:
            seen_residues.add(resid)
            self.residues_missing_neighbor.append(resid)
        else:
          self.site_labels_no_para.append(label)
    if not sel_h_not_in_para.all_eq(False):
      sel_h_not_in_para = sel_h_not_in_para.set_selected(water_selection, False)
      self.model = self.model.select(~sel_h_not_in_para)
    self.time_remove_H_nopara = round(time.time()-t0, 2)
    # Reset occupancies, ADPs and idealize H atom positions
    # -----------------------------------------------------
    t0 = time.time()
    self.model.reset_adp_for_hydrogens(scale = self.adp_scale, keep_aniso=True)
    self.model.reset_occupancy_for_hydrogens_simple()
    self.time_reset = round(time.time()-t0, 2)
    t0 = time.time()
    self.model.idealize_h_riding()
    self.time_idealize = round(time.time()-t0, 2)

    # Remove H atoms that are involved in links (bonds, metal coordination, etc)
    # --------------------------------------------------------------------------
    t0 = time.time()
    # CH2 names must be stereo-correct before a link picks which H to drop
    name_prochiral_h(self.model.get_hierarchy(), self.model.get_mon_lib_srv())
    self.exclude_H_on_links()
    self.time_remove_H_on_links = round(time.time()-t0, 2)


    # TODO: this should be ideally done *after* reduce optimization
    #if not self.exclude_water:
    #  self.model.add_hydrogens(1., occupancy=0.)

    self.model, _ = workaround_003(model = self.model)

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

  def place_n_terminal_propeller(self, pdb_hierarchy):
    '''
    Place NH3 at residue #1 or at first residue in chain
    Changes hierarchy in place
    '''
    for m in pdb_hierarchy.models():
      for chain in m.chains():
        rgs = chain.residue_groups()[0]
        # by default, place NH3 only at residue with resseq 1
        if (self.n_terminal_charge == 'residue_one' and rgs.resseq_as_int() != 1):
          continue
        elif (self.n_terminal_charge == 'first_in_chain'):
          pass
        # SAC in 5xdq, 5zcp. Never needs propeller. Also AYA
        n = None
        for ag in rgs.atom_groups():
          n = ag.get_atom('N') # assumes atom name "N"
          if n: break
        if not n: continue
        bonds = bonds_in_restraints(n, exclude_hydrogens=True)
        heavies = 2
        if ag.resname in ['PRO']: # needs a PRO child lookup
          heavies = 3
        if len(bonds) >= heavies: continue
        if (get_class(name=ag.resname) in
            ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
          for ag_h in rgs.atom_groups():
            h = ag_h.get_atom('H')
            if h:
              ag_h.remove_atom(h)
        # TODO make the function below smart, so it
        # 1) knows when to add H1H2H3 or not
        # 2) renames H to H1 (so no need to remove it beforehand)
        rc = add_n_terminal_hydrogens_to_residue_group(rgs) # rc is always empty list?

  # ----------------------------------------------------------------------------

  def add_missing_H_atoms_at_bogus_position(self, exclude_water):
    '''
    Add missing H atoms at bogus positions to the pdb_hierarchy

    This procedure changes the hierarchy in place.
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
    pdb_hierarchy = self.model.get_hierarchy()
    mon_lib_srv = self.model.get_mon_lib_srv()
    # Load any user-provided restraint CIF objects into the server so that
    # mon_lib_query() finds them and does not generate auto_* CIF objects
    # with empty type_energy values that would overwrite the user entries.
    ro = self.model.get_restraint_objects()
    if ro:
      for fname, cif_obj in ro:
        if cif_obj is not None:
          try:
            mon_lib_srv.process_cif_object(
              cif_object=cif_obj, file_name=fname)
          except Exception:
            pass
    #XXX This breaks for 1jxt, residue 2, TYR
    for m in pdb_hierarchy.models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          n_atom_groups = len(rg.atom_groups())
          for ag in rg.atom_groups():
            if n_atom_groups > 2 and ag.altloc == '':
              continue
            #print list(ag.atoms().extract_name())
            if exclude_water:
              if(get_class(name=ag.resname) == "common_water"): continue
            actual = [a.name.strip().upper() for a in ag.atoms()]
            #
            mlq, cif_object = mon_lib_query(residue=ag, mon_lib_srv=mon_lib_srv, raise_sorry=False)
            if mlq is None:
              self.no_H_placed_mlq.append(ag.resname)
              continue

            if cif_object:
              ro = self.model.get_restraint_objects()
              if ro is None: ro=[]
              auto_name = 'auto_%s' % ag.resname
              ro.append((auto_name, cif_object))
              self.model.set_restraint_objects(ro)
              if auto_name not in self.auto_restraint_names:
                self.auto_restraint_names.append(auto_name)

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
            if mlq.test_for_peptide(atom_dict):
              atom_dict = _remove_atoms(atom_dict, ['H2', 'HXT'])
            elif mlq.test_for_rna_dna(atom_dict):
              atom_dict = _remove_atoms(atom_dict, ["HO3'", 'HO3*'])
            for k, v in six.iteritems(atom_dict):
              if(v.type_symbol=="H"):
                expected_h.append(k)
              #else:
              #  expected_ha.append(k)
            #print('expected H', expected_h)
            #
            # TODO start
            # temporary fix until v3 names are in mon lib
            if (get_class(name=ag.resname) in
                ['common_amino_acid', 'modified_amino_acid', 'd_amino_acid']):
              for altname in alternative_names:
                if (altname[0] in expected_h and altname[1] in expected_h):
                  if (atom_dict[altname[0]].type_energy == 'HCH2' and
                      atom_dict[altname[1]].type_energy == 'HCH2'):
                    expected_h.append(altname[2])
                    expected_h.remove(altname[0])
                    #print('renamed %s to %s' % (altname[0], altname[2]))
            # TODO end
            missing_h = list(set(expected_h).difference(set(actual)))
            if 0: print(ag.resname, missing_h)
            #new_xyz = ag.atoms().extract_xyz().mean()
            new_xyz = flex.double(ag.atoms().extract_xyz().mean()) + \
              flex.double([0.5,0.5,0.5])
            new_xyz = tuple(new_xyz)

            hetero = ag.atoms()[0].hetero
            segid = ag.atoms()[0].segid

            for mh in missing_h:
              if len(mh) < 4: mh = (' ' + mh).ljust(4)
              a = (iotbx.pdb.hierarchy.atom()
                .set_name(new_name=mh)
                .set_element(new_element="H")
                .set_xyz(new_xyz=new_xyz)
                .set_hetero(new_hetero=hetero)
                .set_segid(new_segid=segid))

              ag.append_atom(a)
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
    all_proxies = [p for p in bond_proxies_simple]
    for proxy in asu:
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
    for proxy in all_proxies:
      if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
      elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
      else: assert 0 # never goes here
      bonds.setdefault(i,[])
      bonds[i].append(j)
      bonds.setdefault(j,[])
      bonds[j].append(i)
      # Exception for HIS HD1 and HE2
      if (atoms[i].parent().resname == 'HIS' and
        atoms[i].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
      if (atoms[j].parent().resname == 'HIS' and
        atoms[j].name.strip() in ['HD1','DD1', 'HE2', 'DE2']): continue
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
      return orders.get(frozenset((atoms[i].name.strip(), atoms[j].name.strip())), 1)
    for i_seq in sorted(reversed(list(sel_remove)), key=_drop_order):
      j_seq=parent_dict[i_seq]
      # need to add the use of atomic charge
      valences=get_valences(elements[j_seq])
      number_of_bonds=sum(_bond_order(j_seq, n) for n in bonds[j_seq])
      if number_of_bonds in valences:
        # remove this H from delection
        remove_from_sel_remove.append(i_seq) # ??
        del removed_dict[i_seq]
      else:
        bonds[j_seq].remove(i_seq)
        bonds[i_seq].remove(j_seq)

    #print(remove_from_sel_remove)
    fsc0=grm.geometry.shell_sym_tables[0].full_simple_connectivity()
    fsc1=grm.geometry.shell_sym_tables[1].full_simple_connectivity()
    #fsc2=grm.geometry.shell_sym_tables[2].full_simple_connectivity()
    for _i in remove_from_sel_remove:
      parent = fsc0[_i][0]
      first_neighbors = fsc1[_i]
      fn_filtered = [item for item in first_neighbors if item not in sel_remove]
      n_kept = sum(1 for k in remove_from_sel_remove if fsc0[k][0] == parent)
      coordp = matrix.col(atoms[parent].xyz)
      # Two H kept on an sp3 atom with two heavy neighbours: a methyl that
      # became a CH2 via a link (e.g. B0I CB1/CB2 in 6b17, linked to a thioether
      # S). Place both H tetrahedrally, symmetric across the plane of the two
      # heavy neighbours, so the real C-C and C-S(link) directions are respected
      # (not eclipsed). The single-H branches below cannot handle this: they
      # would place both H at the same in-plane bisector site.
      if n_kept == 2 and len(fn_filtered) == 2:
        siblings = sorted(k for k in remove_from_sel_remove if fsc0[k][0] == parent)
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
        atoms[_i].name = ' %s  ' % atoms[_i].element.strip()
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
      sel_remove=list(sel_remove)
      for r in remove_from_sel_remove:
        sel_remove.remove(r)
        if verbose: print('keep',atoms[r].quote())
      sel_remove=flex.size_t(sel_remove)
    #
    def _removed_label(v):
      return v if isinstance(v, str) else origin_ids.get_origin_key(v)
    sl_removed = [(atom.id_str().replace('pdb=','').replace('"',''),
                   _removed_label(removed_dict[atom.i_seq]))
        for atom in self.model.get_hierarchy().atoms().select(sel_remove)]
    #
    if sel_remove:
      self.model = self.model.select(~flex.bool(self.model.size(), sel_remove))
    self.sl_removed = sl_removed
    self.exclusion_iseqs = exclusion_iseqs


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
      site_labels_disulfides = self.site_labels_disulfides,
      site_labels_no_para = self.site_labels_no_para,
      site_labels_tertiary_amide = self.site_labels_tertiary_amide,
      site_labels_missing_neighbor = self.site_labels_missing_neighbor,
      residues_missing_neighbor = self.residues_missing_neighbor)

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
