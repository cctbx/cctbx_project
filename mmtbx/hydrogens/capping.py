'''
Capping H on dangling bonds (reduce2 capping=True, for QM).

A dangling bond is a valence of a non-H atom left unsatisfied because a
non-H partner is missing from the model: a truncated side chain, a ligand
with missing atoms, a residue so incomplete that qmready's own OXT/NH3+
treatment skipped it. Capping puts H on such atoms, as many as the missing
bonds need. Only H atoms are ever built; a missing non-H atom is never
rebuilt, so a missing internal atom (one with two present dangling
neighbours) is refused (rule 7).

Detection (find_dangling_atoms), on a model with a restraints manager:

1. Bond table from all bond proxies, simple and ASU, without the metal
   coordination origin id: a coordination bond satisfies no valence (same
   exclusion as reduce_hydrogen.exclude_H_without_link_angle).
2. For every non-H, non-metal atom in a residue with a dictionary, a valence
   budget:
   - dictionary partners: the atoms bonded to it in the dictionary, H
     included, metals excluded (a partner whose element is a metal, or a
     bond typed metal), each with its bond order;
   - V = sum of the orders of all partners. The CCD describes the free,
     complete molecule (C carries O and OXT, N carries H and H2, a 5' P
     carries OP3), so V is the atom's full valence; a link partner in
     another residue simply takes the place of a CCD leaving atom;
   - missing: the heavy partners absent from the atom's conformer (its own
     altloc plus blank; for a blank atom any altloc). An atom with no
     missing heavy partner is not examined further: its protonation is
     reduce2's;
   - used = sum of the orders of the present partners + 1 per bonded heavy
     neighbour that is not a dictionary partner of the atom or is in
     another residue (peptide, glycosidic, disulfide, isopeptide, covalent
     ligand) + 1 per bonded H. A neighbour split into altlocs counts once,
     the first altloc seen;
   - n_fill = V - used.
3. n_fill > 0: the atom is dangling. Examples: a backbone C without O
   (V 4: CA, O double, OXT; used CA + N of the next residue = 2) -> 2; ASP
   CG without OD2 but linked to a LYS NZ (V 4; used CB + OD1 double + NZ =
   4) -> 0; an internal nucleotide P (no OP3, O3' of the previous residue
   linked) -> 0; a 5' P without OP3 -> 1; HEM pyrrole N with FE (FE is not
   a partner, the coordination bond is not counted) -> not examined.
4. coordination = bonded non-metal neighbours (H included) + n_fill: the
   number of substituents after filling. 2, 3 and 4 can be placed; anything
   else (a hypervalent S or P missing a double-bonded partner) is skipped
   with a reason: not an error, QM users can decide.
5. A dangling atom with no bonded neighbour at all (a lone CB) cannot be
   oriented: Sorry naming the atom; delete it or complete the residue.
   Water, ions and single-atom residues have no heavy partner to miss.
6. Residues without a dictionary are listed in unknown_residues. With
   guess_for_unknowns=True every non-H atom in them gets the same budget
   with V by element and formal charge (C 4, N 3+q, O 2+q, S 2+q, P 5,
   halogens 1; other elements and single-atom residues are skipped with a
   reason) and used = bonded non-metal neighbours, H included, one each.
   Bonds come from the bond table if the residue has any, else from
   covalent radii (sum + 0.4 A). Without guessing they are only listed,
   for the caller to Sorry.
7. Internal gaps: a missing heavy atom in the missing_names of two or more
   altloc-compatible records of one residue group (for example TYR without
   CE1: CD1 and CZ; DC without C2': C1' and C3'; ARG without CG only: CB
   and CD). Their
   caps would point into the same hole, 0.5-1.3 A apart. Detection still
   returns the records and lists the gap in internal_gaps;
   cap_dangling_atoms refuses it with a Sorry naming the atom and its
   neighbours, before placing anything: rebuild the atom or remove the
   fragment. A truncation (for example ARG at CB, LYS at CD, a C without O)
   leaves only one present neighbour of each missing atom and is no gap.
8. Hypervalent parents: an H on P, or on an S with 3 or more substituents
   (sulfate, sulfonate, sulfone, sulfoxide, sulfonium), changes the
   chemistry (a P-H is a
   phosphonate, not a phosphate) and the charge module cannot count it, so
   such records are listed in hypervalent and cap_dangling_atoms refuses
   them with a Sorry: add the missing oxygen (for example OP3 on a 5' P
   deposited without it, or the O3' lost at an internal chain break) or
   remove the fragment. A thiol or thioether S (2 substituents) is capped.
   Nucleic-acid chain ends are not dangling: qmready places HO3' and HO5'
   there (reduce_hydrogen na_terminal_h); a 5' phosphate keeps its charge.

9. Cut bonds (cut_bonds, {i_seq: n}): bonds of the model's atoms to atoms
   that are not in the model, a fragment cut out of a larger model (QM
   fragments of qrefine; cut_bonds_by_distance finds them in the parent).
   An atom without a record gets one with n_fill n, coordination its bonded
   neighbours + n, and its free dictionary H names (H2 on a backbone N, H
   on a PRO N, HG on a cut SG, HO3' on a cut O3'); an atom with a record
   keeps it, its n_fill already counts the absent partner (a cut backbone
   C: V 4, used CA + O double = 3, fill 1). An atom the dictionary pass
   found saturated (n_fill <= 0), H and metals are ignored. In guess mode
   n is added to the guessed n_fill.

Dictionary and bond orders. Atom and bond lists, and the atom names, are the
CCD's wherever reduce_hydrogen._ccd_describes() holds, else the restraint
dictionary's (a user CIF, or a residue the CCD does not have): geostd amino
acids lack OXT and H2, so only the CCD gives V as above. Bond orders are
reduce_hydrogen._bond_orders(): CCD DOUB/TRIP, or dictionary double/triple,
2 and 3; everything else counts 1.

Aromatic rings and the deloc fallback. The CCD writes aromatic rings Kekule
style, so a ring atom's two neighbours across a missing ring atom get 2 and 1
H (for example TYR without CE1: CD1 one, CZ two) and the fragment stays
closed-shell.
Only a dictionary without CCD backing (a user CIF) can type a bond aromatic or
deloc; such a bond counts as 1 (the deloc fallback), i.e. one H. Where that
leaves the fragment open-shell, reduce2's charge step reports it in the
Capping report as 'open-shell cap on <atom>'.

Missing names. The missing_names of a record (and the 'missing' list of the
report) are the missing dictionary heavy partners except the CCD leaving atoms
that a link replaces: a backbone C bonded to the next residue is not missing
OXT, a 5' P linked to the previous O3' is not missing OP3. A leaving atom is
missing only when nothing takes its place. Only heavy atoms are ever in
missing_names.

Building (cap_dangling_atoms): the geometry follows the steric number,
coordination plus lone pairs (N 1, O 2, S 2 when its coordination is 2 or
less; at most 4): 4 is tetrahedral, 3 is 120 degrees, 2 is linear. Directions
from fill_directions() about the bonded neighbours counted in the record; free
rotations are sampled in 10 degree steps, starting at the nearest off-axis
environment atom, against the atoms of the parent's MODEL within 4 A but the
parent and its bonded neighbours (no symmetry mates), and the one keeping the
H farthest from them wins. A bonded neighbour is at the same distance from
every rotation, so leaving it out lets the 1-3 atoms decide: a CH3 on a CB is
staggered, and the result does not depend on the frame of the input. X-H is
the monomer library's (C with 3 substituents: aromatic C-H), else the sum of
the covalent radii. Each H joins its parent's atom group with the parent's
occupancy, B, segid and hetero flag. Afterwards a cap closer than 1.0 A to
any atom of its MODEL but its parent, symmetry mates included, is a Sorry
(the clash check is symmetry-aware), as in reduce2's qmready clash check.
Several MODELs are capped independently: no atom of one MODEL orients or
clashes with a cap of another. The bond table of the detection result is
reused.

Names of the caps, unique within the conformer and at most 4 characters:
1. dictionary H names that are free (HB2, HB3 on a truncated CB);
2. 'H' + the missing heavy atom's name (HCG, HOP3);
3. HC1, HC2, ... when 2 is too long or taken.

Guess mode (guess_for_unknowns, reduce2's capping_guess_for_unknowns): V by
element and formal charge, bonds from the bond table or covalent radii, as in
rule 6. The bonds found by distance that the bond table lacks are returned as
guessed_bonds and, with the cap H, given to the charge calculation by
extra_sigma_bonds().

Sorry conditions. Here: an internal gap (rule 7; at most 20 listed); a cap
closer than 1.0 A to any altloc-compatible atom of its MODEL other than the
cap's own parent; a dangling atom with no bonded neighbour; no covalent
radius for the parent's element (xh_length); records that do not belong to
the model given to cap_dangling_atoms. In reduce2: capping without qmready;
approach=remove; capping_guess_for_unknowns without capping; a residue
without a dictionary when not guessing.

CCD environment. Without the CCD (a cctbx without chem_data's chemical
components) a residue the restraint dictionary knows falls back to it: no
OXT, H2 or OP3 and no bond orders, so V can be short. find_dangling_atoms
then writes one line per residue name to its log.

Chemistry notes. V comes from the CCD, which fixes a protonation state for
the free molecule; protonation itself is reduce2's. Following the rules
exactly gives some surprising but closed-shell, charge-preserving caps:
- HIS without CE1: ND1 has V 4 (CG, CE1 double, HD1) and gets 2 H, an
  NH3+-like centre on an sp2 CG (and NE2 also misses CE1: an internal gap,
  refused);
- ASP without OD1, the double-bonded O: CG gets 2 H, CH2-O(-) with OD2 bare;
- a C with OXT but no O: 2 H, CH2-OXT(-).

Limitations.
- The capped model cannot be processed by pdb_interpretation (the cap names
  are invented): it is for QM use. Capping is the last step, and the model
  manager reduce2 returns is unprocessed.
- A blank-altloc parent whose dictionary partner exists in only one altloc is
  treated as complete (the partner is present in some conformer).
- reduce2's optimizer rotates H (the terminal HO3' and HO5' included) against
  the asymmetric unit only; the qmready clash check is symmetry-aware, so a
  hydroxyl H at a crystal contact can still end in its Sorry. A QM model is
  usually run without crystal symmetry (no CRYST1), where this cannot occur.
  Documented by mmtbx/regression/tst_reduce2_symmetry_optimizer.py (an
  expected failure until the optimizer is fixed).
- reduce2's total charge is wrong for a model with alternate conformations
  (electrons needs an altloc-free hierarchy and gets the whole model): a
  pre-existing reduce2 issue, independent of capping, documented by
  mmtbx/regression/tst_reduce2_altloc_charge.py (an expected failure).
- A residue reduced to N and CA gets no NH3+ from qmready (the propeller needs
  C) and its N is not examined by capping.
- reduce2 drops an input H on a truncated atom even with keep_existing_H=True;
  capping adds it back under its dictionary name.
- A missing internal atom is not capped (rule 7); the remedy is to rebuild
  it or remove the fragment.
'''
from __future__ import absolute_import, division, print_function
from libtbx import group_args
from libtbx.utils import Sorry
from cctbx.geometry_restraints.linking_class import linking_class
from cctbx import crystal, uctbx
from mmtbx.ligands.electrons import atom_property
from mmtbx.hydrogens import reduce_hydrogen
import iotbx.pdb.hierarchy
import math
from scitbx import matrix
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp

ext = bp.import_ext("cctbx_geometry_restraints_ext")

# valence by element for residues without a dictionary; N, O and S add the
# formal charge
GUESS_VALENCE = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 5,
                 'F': 1, 'CL': 1, 'BR': 1, 'I': 1}
GUESS_CHARGED = ('N', 'O', 'S')

# guessed bond: distance below the sum of the covalent radii plus this
GUESS_BOND_TOLERANCE = 0.4
GUESS_BOND_CUTOFF = 2.5

# steric numbers: substituents plus lone pairs
TETRAHEDRAL, TRIGONAL, LINEAR = 4, 3, 2
LONE_PAIRS = {'N': 1, 'O': 2}
TETRAHEDRAL_ANGLE = 109.5
TRIGONAL_ANGLE = 120.
CAP_BOND_LENGTH = 1.0 # A, H placed this far for scoring a rotation

# building the caps
H_COVALENT_RADIUS = 0.31 # A, for X-H of elements without a library length
ENVIRONMENT_RADIUS = 4.0 # A, sites a free rotation keeps away from
CLASH_DISTANCE = 1.0 # A, as reduce2's qmready clash check

# ==============================================================================

def atom_key(atom):
  '''
  Identity of an atom that survives sorting the hierarchy (i_seq does not):
  (model id, chain id, resseq, icode, resname, altloc, name).
  '''
  ag = atom.parent()
  rg = ag.parent()
  chain = rg.parent()
  return (chain.parent().id.strip(), chain.id.strip(), rg.resseq.strip(),
          rg.icode.strip(), ag.resname.strip(), ag.altloc.strip(),
          atom.name.strip())

def _residue_label(ag):
  '''
  'resname chain resseqicode' of an atom group, as reduce2 prints residues.
  '''
  rg = ag.parent()
  return '%s %s %s%s' % (ag.resname.strip(), rg.parent().id.strip(),
    rg.resseq.strip(), rg.icode.strip())

def _atom_label(atom):
  '''
  Residue label, atom name and altloc, for skipped atoms.
  '''
  label = '%s %s' % (_residue_label(atom.parent()), atom.name.strip())
  altloc = atom.parent().altloc.strip()
  if altloc: label += ' (altloc %s)' % altloc
  return label

def bond_table(model):
  '''
  {i_seq: set(i_seq)} from all bond proxies, simple and ASU, without metal
  coordination (it satisfies no valence).
  '''
  grm = model.get_restraints_manager().geometry
  bps, asu = grm.get_all_bond_proxies(sites_cart = model.get_sites_cart())
  metal = linking_class().get_origin_id('metal coordination')
  table = {}
  for proxy in list(bps) + list(asu):
    if proxy.origin_id == metal: continue
    if   isinstance(proxy, ext.bond_simple_proxy): i,j = proxy.i_seqs
    elif isinstance(proxy, ext.bond_asu_proxy):    i,j = proxy.i_seq, proxy.j_seq
    else: continue
    table.setdefault(i, set()).add(j)
    table.setdefault(j, set()).add(i)
  return table

# ==============================================================================

def _dictionary_atoms(resname, mon_lib_srv, log=None):
  '''
  (atoms, bonds) of resname's dictionary, None if no dictionary knows it:
  atoms [(atom_id, type_symbol, is leaving atom)] in dictionary order, bonds
  [(atom_id_1, atom_id_2, type)]. CCD first where it describes the residue
  (the free molecule, with OXT, H2 and OP3, flagged as leaving atoms), else
  the restraint dictionary (no leaving atoms, bond types). With log, a line
  when the CCD should describe a residue the restraint dictionary knows but
  has no entry for it (the valences may then be incomplete).
  '''
  from mmtbx.chemical_components import get_cif_dictionary
  try: cc = mon_lib_srv.get_comp_comp_id_direct(resname)
  except Exception: cc = None
  cc_cif = None
  if reduce_hydrogen._ccd_describes(cc, resname):
    try: cc_cif = get_cif_dictionary(resname)
    except Exception: cc_cif = None
  if cc_cif and cc_cif.get('_chem_comp_atom'):
    atoms = [(a.atom_id, a.type_symbol,
              getattr(a, 'pdbx_leaving_atom_flag', 'N') == 'Y')
             for a in cc_cif['_chem_comp_atom']]
    bonds = [(b.atom_id_1, b.atom_id_2, '')
             for b in cc_cif.get('_chem_comp_bond', [])]
  elif cc is not None and len(cc.atom_list):
    if log is not None and reduce_hydrogen._ccd_describes(cc, resname):
      print('Capping: no CCD entry for %s; its valences come from the '
        'restraint dictionary and may be incomplete (no OXT, H2, OP3, bond '
        'orders)' % resname, file=log)
    atoms = [(a.atom_id, a.type_symbol, False) for a in cc.atom_list]
    bonds = [(b.atom_id_1, b.atom_id_2, str(b.type).strip().lower())
             for b in cc.bond_list]
  else:
    return None
  return atoms, bonds

def _dictionary_partners(resname, mon_lib_srv, is_metal, cache, log=None):
  '''
  {name: [(partner name, partner is H, partner is a leaving atom)]} of
  resname's dictionary (_dictionary_atoms), metal partners and metal bonds
  left out; None if no dictionary knows resname.
  '''
  if resname in cache: return cache[resname]
  entries = _dictionary_atoms(resname, mon_lib_srv, log=log)
  if entries is None:
    cache[resname] = None
    return None
  atoms, bonds = entries
  elements = dict((n, e.strip().upper()) for n, e, leaving in atoms)
  leaving = set(n for n, e, is_leaving in atoms if is_leaving)
  partners = {}
  for n1, n2, bond_type in bonds:
    if bond_type == 'metal': continue
    e1, e2 = elements.get(n1, ''), elements.get(n2, '')
    if is_metal(e1) or is_metal(e2): continue
    partners.setdefault(n1, []).append((n2, e2 in ('H', 'D'), n2 in leaving))
    partners.setdefault(n2, []).append((n1, e1 in ('H', 'D'), n1 in leaving))
  cache[resname] = partners
  return partners

def _conformer_names(rg):
  '''
  {altloc: set(names)} of a residue group's atom groups.
  '''
  result = {}
  for ag in rg.atom_groups():
    result.setdefault(ag.altloc.strip(), set()).update(
      a.name.strip() for a in ag.atoms())
  return result

def _present_names(atom, names_by_altloc):
  '''
  Names in the atom's conformer: its own altloc plus blank; for a blank
  atom, any altloc.
  '''
  altloc = atom.parent().altloc.strip()
  if altloc:
    return names_by_altloc.get(altloc, set()) | names_by_altloc.get('', set())
  result = set()
  for names in names_by_altloc.values(): result |= names
  return result

def _conformer_neighbours(atom, j_seqs, atoms, is_metal):
  '''
  Bonded non-metal neighbours of atom in one conformer: for an atom in an
  altloc, blank or the same altloc; for a blank atom, blank or the first
  altloc seen (CA A and CA B on a blank CB are one neighbour).
  '''
  altloc = atom.parent().altloc.strip()
  result = []
  for j in sorted(j_seqs):
    other = atoms[j]
    if is_metal(reduce_hydrogen._atom_element(other)): continue
    other_altloc = other.parent().altloc.strip()
    if other_altloc:
      if not altloc: altloc = other_altloc
      if other_altloc != altloc: continue
    result.append(other)
  return result

def _guess_neighbours(atom, rg_has_bonds, j_seqs, atoms, sites, is_metal):
  '''
  Neighbours for guessing: the bond table if the residue has bonds in it,
  else every atom closer than the sum of the covalent radii + 0.4 A.
  '''
  if rg_has_bonds:
    return _conformer_neighbours(atom, j_seqs, atoms, is_metal)
  r_atom = reduce_hydrogen._covalent_radius(reduce_hydrogen._atom_element(atom))
  if r_atom is None: return []
  close = []
  d2 = (sites - atom.xyz).dot()
  for j in (d2 < GUESS_BOND_CUTOFF**2).iselection():
    if j == atom.i_seq: continue
    r = reduce_hydrogen._covalent_radius(
      reduce_hydrogen._atom_element(atoms[j]))
    if r is None: continue
    if d2[j]**0.5 < r_atom + r + GUESS_BOND_TOLERANCE: close.append(j)
  return _conformer_neighbours(atom, close, atoms, is_metal)

def _record(atom, n_fill, neighbours, missing_names, h_names, source,
            records, skipped):
  '''
  Append a dangling_atom record, or a skipped entry if the coordination
  cannot be placed. Sorry if nothing is bonded to orient the H.
  '''
  if not neighbours:
    raise Sorry('Capping: %s needs %d H but has no bonded neighbour to '
      'orient them; delete it or complete the residue.' % (
      _atom_label(atom), n_fill))
  coordination = len(neighbours) + n_fill
  if coordination not in (2, 3, 4):
    skipped.append((_atom_label(atom), 'needs %d H, coordination %d cannot '
      'be placed' % (n_fill, coordination)))
    return
  records.append(group_args(
    key           = atom_key(atom),
    i_seq         = atom.i_seq,
    n_fill        = n_fill,
    coordination  = coordination,
    missing_names = missing_names,
    h_names       = h_names,
    source        = source))

def find_dangling_atoms(model, guess_for_unknowns=False, log=None,
                        cut_bonds=None):
  '''
  Non-H atoms with a dangling bond (see the module docstring), in model
  order. model: mmtbx.model.manager with a restraints manager. log: where to
  say that a residue's CCD entry could not be read (None: silent).
  cut_bonds: {i_seq: n}, bonds of the model's atoms to atoms not in the
  model (rule 9); an atom without a record gets one with n_fill n.
  Returns group_args(records, unknown_residues, skipped, guessed_bonds,
  internal_gaps, bond_table):
    records: group_args(key, i_seq, n_fill, coordination, missing_names,
      h_names, source) per dangling atom; key is atom_key(), i_seq in model;
      missing_names the missing dictionary partners, h_names the dictionary
      H partners absent from the conformer (free to reuse); source
      'dictionary' or 'guess'. missing_names omits the CCD leaving atoms
      that a link replaces (OXT, OP3 ...): a peptide-bonded C is not
      missing OXT
    unknown_residues: labels ('UNL A 1') of residues without a dictionary
    skipped: (atom label, reason) of dangling atoms that get no H
    guessed_bonds: (i_seq, j_seq), i_seq < j_seq, of the bonds guessing
      found by distance that the bond table does not have (see
      extra_sigma_bonds())
    internal_gaps: group_args(residue_key, missing_name, parent_keys) per
      internal gap (_internal_gaps()); cap_dangling_atoms() refuses them
    bond_table: bond_table(model), for cap_dangling_atoms() to reuse
  '''
  mon_lib_srv = model.get_mon_lib_srv()
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  sites = atoms.extract_xyz()
  table = bond_table(model)
  properties = atom_property()
  is_metal = lambda element: properties.is_metal(element.strip().upper())
  partners_cache, orders_cache = {}, {}
  records, unknown_residues, skipped = [], [], []
  guessed_bonds = set()
  saturated = set() # examined by the dictionary pass, n_fill <= 0
  for rg in hierarchy.residue_groups():
    names_by_altloc = None
    for ag in rg.atom_groups():
      resname = ag.resname.strip()
      partners = _dictionary_partners(resname, mon_lib_srv, is_metal,
                                      partners_cache, log=log)
      if partners is None:
        label = _residue_label(ag)
        if label not in unknown_residues: unknown_residues.append(label)
        if guess_for_unknowns:
          _guess_atom_group(ag, table, atoms, sites, is_metal, records,
                            skipped, guessed_bonds, cut_bonds=cut_bonds)
        continue
      orders = reduce_hydrogen._bond_orders(resname, mon_lib_srv, orders_cache)
      order = lambda a, b: orders.get(frozenset((a, b)), 1)
      for atom in ag.atoms():
        if atom.element_is_hydrogen(): continue
        if is_metal(reduce_hydrogen._atom_element(atom)): continue
        name = atom.name.strip()
        own = partners.get(name)
        if not own: continue
        if names_by_altloc is None: names_by_altloc = _conformer_names(rg)
        present = _present_names(atom, names_by_altloc)
        missing = [p for p, is_h, l in own if not is_h and p not in present]
        if not missing: continue
        valence = sum(order(name, p) for p, is_h, l in own)
        heavy_partners = set(p for p, is_h, l in own if not is_h)
        used = sum(order(name, p) for p in heavy_partners if p in present)
        neighbours = _conformer_neighbours(atom, table.get(atom.i_seq, ()),
                                           atoms, is_metal)
        n_link = 0
        for other in neighbours:
          if other.element_is_hydrogen():
            used += 1
          elif (other.name.strip() not in heavy_partners or
                other.parent().parent().memory_id() != rg.memory_id()):
            # link partner, in place of a leaving atom
            n_link += 1
        used += n_link
        n_fill = valence - used
        if n_fill <= 0:
          saturated.add(atom.i_seq)
          continue
        # the leaving atoms the links replace are not missing (C without O,
        # peptide-bonded: O, not OXT)
        for p, is_h, is_leaving in own:
          if n_link and is_leaving and p in missing:
            missing.remove(p)
            n_link -= 1
        h_names = [p for p, is_h, l in own if is_h and p not in present]
        _record(atom, n_fill, neighbours, missing, h_names, 'dictionary',
                records, skipped)
      if cut_bonds:
        # rule 9: atoms that lost a bond to an atom outside the model
        recorded = set(r.i_seq for r in records)
        for atom in ag.atoms():
          n_cut = cut_bonds.get(atom.i_seq, 0)
          if n_cut <= 0: continue
          if atom.element_is_hydrogen(): continue
          if is_metal(reduce_hydrogen._atom_element(atom)): continue
          if atom.i_seq in recorded or atom.i_seq in saturated: continue
          own = partners.get(atom.name.strip(), [])
          if names_by_altloc is None: names_by_altloc = _conformer_names(rg)
          present = _present_names(atom, names_by_altloc)
          h_names = [p for p, is_h, l in own if is_h and p not in present]
          neighbours = _conformer_neighbours(atom, table.get(atom.i_seq, ()),
                                             atoms, is_metal)
          _record(atom, n_cut, neighbours, [], h_names, 'dictionary',
                  records, skipped)
  return group_args(
    records          = records,
    unknown_residues = unknown_residues,
    skipped          = skipped,
    guessed_bonds    = sorted(guessed_bonds),
    internal_gaps    = _internal_gaps(records, atoms),
    hypervalent      = [r for r in records
                        if _is_hypervalent(atoms[r.i_seq], r.coordination)],
    bond_table       = table)

HYPERVALENT_S_COORDINATION = 3
HYDROXYL_LIKE_PARENTS = ('O', 'N', 'S')

def _is_hypervalent(atom, coordination):
  '''
  True for a parent an H cannot cap without changing the chemistry: P (a
  P-H is a phosphonate, not a phosphate) and S with HYPERVALENT_S_COORDINATION
  or more substituents after filling (sulfate, sulfonate, sulfone, sulfoxide,
  sulfonium). A thiol
  or thioether S is not hypervalent.
  '''
  element = atom.element.strip().upper()
  if element == 'P': return True
  return element == 'S' and coordination >= HYPERVALENT_S_COORDINATION

def _hypervalent_line(record, atom):
  '''
  For example: 'Cannot cap DT A 7 P: it is missing O5' and an H on P would be
  a phosphonate, not a phosphate. Add the missing atom or remove ...'
  '''
  element = atom.element.strip().upper()
  if element == 'P': what = 'a phosphonate, not a phosphate'
  else: what = 'a sulfur hydride, not the oxidized or charged S it is'
  if record.missing_names: missing = ' '.join(record.missing_names)
  else: missing = '%d bonded atom(s)' % record.n_fill
  hint = " (OP3 at a 5' end)" if 'OP3' in record.missing_names else ''
  return ('Cannot cap %s: it is missing %s and an H on %s would be %s. Add '
    'the missing atom%s or remove the fragment before running capping.' % (
    _atom_label(atom), missing, element, what, hint))

def _internal_gaps(records, atoms):
  '''
  The missing heavy atoms that two or more altloc-compatible records of one
  residue group miss (for example TYR without CE1: CD1 and CZ; ARG without
  CG only: CB
  and CD): capping both parents would put two H into the same hole.
  [group_args(residue_key, missing_name, parent_keys)] in record order;
  residue_key is atom_key() without altloc and name.
  '''
  by_rg, order = {}, []
  for record in records:
    rg_id = atoms[record.i_seq].parent().parent().memory_id()
    if rg_id not in by_rg:
      by_rg[rg_id] = []
      order.append(rg_id)
    by_rg[rg_id].append(record)
  result = []
  for rg_id in order:
    rg_records = by_rg[rg_id]
    names = []
    for record in rg_records:
      for name in record.missing_names:
        if name not in names: names.append(name)
    for name in names:
      missing_it = [r for r in rg_records if name in r.missing_names]
      if len(missing_it) < 2: continue
      altlocs = sorted(set(r.key[5] for r in missing_it if r.key[5]))
      seen = []
      for altloc in altlocs or ['']:
        parents = [r.key for r in missing_it
                   if not altloc or r.key[5] in ('', altloc)]
        if len(parents) < 2 or parents in seen: continue
        seen.append(parents)
        result.append(group_args(
          residue_key  = parents[0][:5],
          missing_name = name,
          parent_keys  = parents))
  return result

def _internal_gap_line(gap):
  '''
  For example: 'Cannot cap chain A TYR 5: atom CE1 is missing but its
  neighbours CD1 and CZ are present; capping both would overlap. Rebuild ...'
  '''
  model_id, chain_id, resseq, icode, resname = gap.residue_key
  fields = [chain_id, resname, resseq + icode]
  altlocs = sorted(set(k[5] for k in gap.parent_keys if k[5]))
  if altlocs: fields.append(altlocs[0])
  label = 'chain %s' % ' '.join(f for f in fields if f)
  if model_id: label = 'model %s %s' % (model_id, label)
  names = [k[6] for k in gap.parent_keys]
  if len(names) == 2:
    neighbours, which = '%s and %s' % tuple(names), 'both'
  else:
    neighbours = '%s and %s' % (', '.join(names[:-1]), names[-1])
    which = 'all %d' % len(names)
  return ('Cannot cap %s: atom %s is missing but its neighbours %s are '
    'present; capping %s would overlap. Rebuild the missing atom (or remove '
    'the fragment) before running capping.' % (label, gap.missing_name,
    neighbours, which))

def extra_sigma_bonds(result):
  '''
  {i_seq: n} for electrons.run(extra_sigma_bonds=...) from a
  find_dangling_atoms() result: n_fill per record, plus 1 on both atoms of
  each guessed bond (n_fill per record is the number of cap H, one sigma
  bond each). electrons takes its bonds from the restraints manager,
  which has none for a residue without restraints, so every bond found by
  distance has to be given to it as a sigma bond on both atoms, exactly as a
  cap H is (an H bonded by distance too: otherwise it is an isolated H).
  '''
  extra = {}
  for record in result.records:
    extra[record.i_seq] = extra.get(record.i_seq, 0) + record.n_fill
  for i_seq, j_seq in result.guessed_bonds:
    extra[i_seq] = extra.get(i_seq, 0) + 1
    extra[j_seq] = extra.get(j_seq, 0) + 1
  return extra

def cut_bonds_by_distance(hierarchy, selection):
  '''
  {i_seq: n} for find_dangling_atoms(cut_bonds=...) on
  hierarchy.select(selection): the bonds, by covalent distance (rule 6:
  sum of the covalent radii + GUESS_BOND_TOLERANCE), between a selected and
  an unselected heavy non-metal atom, counted per selected atom. The keys
  are the positions of the selected atoms in selection.iselection(), the
  i_seqs of hierarchy.select(selection) after reset_i_seq(). No symmetry:
  the hierarchy is put in a P1 box larger than the cutoff. An unselected atom
  of a residue group that has selected atoms is an alternate conformer of
  the selection (fragments are residue-group complete), not a cut partner.
  '''
  atoms = hierarchy.atoms()
  sites = atoms.extract_xyz()
  selected_groups = set()
  for i in selection.iselection():
    selected_groups.add(atoms[i].parent().parent().memory_id())
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=sites, buffer_layer=5)
  cs = crystal.symmetry(unit_cell=box.unit_cell, space_group_symbol='P1')
  xrs = hierarchy.extract_xray_structure(crystal_symmetry=cs)
  xrs.set_sites_cart(box.sites_cart)
  pst = xrs.pair_asu_table(distance_cutoff=GUESS_BOND_CUTOFF
    ).extract_pair_sym_table()
  properties = atom_property()
  position = {}
  for k, i_seq in enumerate(selection.iselection()): position[i_seq] = k
  result = {}
  for i, pair_sym_dict in enumerate(pst):
    for j in pair_sym_dict.keys():
      if selection[i] == selection[j]: continue
      outside = atoms[j if selection[i] else i]
      if outside.parent().parent().memory_id() in selected_groups: continue
      pair = []
      for a in (atoms[i], atoms[j]):
        if a.element_is_hydrogen(): break
        e = reduce_hydrogen._atom_element(a)
        if properties.is_metal(e): break
        r = reduce_hydrogen._covalent_radius(e)
        if r is None: break
        pair.append(r)
      if len(pair) != 2: continue
      d = (matrix.col(sites[i]) - matrix.col(sites[j])).length()
      if d >= pair[0] + pair[1] + GUESS_BOND_TOLERANCE: continue
      inside = i if selection[i] else j
      k = position[inside]
      result[k] = result.get(k, 0) + 1
  return result

def _has_internal_bonds(rg, table):
  '''
  True if the bond table bonds any two atoms of residue group rg.
  '''
  i_seqs = set(a.i_seq for a in rg.atoms())
  return any(table.get(i, set()) & i_seqs for i in i_seqs)

def _guess_atom_group(ag, table, atoms, sites, is_metal, records, skipped,
                      guessed_bonds, cut_bonds=None):
  '''
  Valence budget by element and formal charge for an atom group without a
  dictionary (rule 6). The bonds found by distance that the bond table does
  not have are added to guessed_bonds as (i_seq, j_seq), i_seq < j_seq.
  '''
  heavy = [a for a in ag.atoms() if not a.element_is_hydrogen()]
  if len(heavy) == 1:
    skipped.append((_atom_label(heavy[0]),
      'single-atom residue without a dictionary, valence not guessed'))
    return
  rg_has_bonds = _has_internal_bonds(ag.parent(), table)
  for atom in heavy:
    element = reduce_hydrogen._atom_element(atom)
    if element not in GUESS_VALENCE:
      skipped.append((_atom_label(atom),
        'no valence to guess for element %s' % element))
      continue
    valence = GUESS_VALENCE[element]
    if element in GUESS_CHARGED: valence += atom.charge_as_int()
    bonded = table.get(atom.i_seq, ())
    neighbours = _guess_neighbours(atom, rg_has_bonds, bonded, atoms, sites,
                                   is_metal)
    for other in neighbours:
      if other.i_seq not in bonded:
        guessed_bonds.add(tuple(sorted((atom.i_seq, other.i_seq))))
    n_fill = valence - len(neighbours)
    if cut_bonds: n_fill += cut_bonds.get(atom.i_seq, 0)
    if n_fill <= 0: continue
    _record(atom, n_fill, neighbours, [], [], 'guess', records, skipped)

# ==============================================================================

def steric_number(element, coordination):
  '''
  Substituents plus lone pairs, at most 4: N has 1 lone pair, O 2, S 2 when
  it has up to 2 substituents (else 0); everything else 0.
  '''
  element = element.strip().upper()
  lone_pairs = LONE_PAIRS.get(element, 0)
  if element == 'S' and coordination <= 2: lone_pairs = 2
  return min(TETRAHEDRAL, coordination + lone_pairs)

def _clearance(directions, parent, environment):
  '''
  Minimum distance from an H at CAP_BOND_LENGTH along any of directions to
  any environment site (infinite for an empty environment).
  '''
  if len(environment) == 0: return float('inf')
  result = float('inf')
  for d in directions:
    h = parent + d * CAP_BOND_LENGTH
    result = min(result, flex.min((environment - h.elems).norms()))
  return result

def _best_candidate(candidates, parent, environment):
  '''
  The candidate (a list of directions) with the largest clearance; ties go to
  the first.
  '''
  best, best_score = None, -1.
  for candidate in candidates:
    score = _clearance(candidate, parent, environment)
    if best is None or score > best_score + 1e-9:
      best, best_score = candidate, score
  return best

def _fixed_candidates(units, steric, n_fill):
  '''
  Candidate direction lists when 2 or 3 neighbours are known.
  '''
  if len(units) == 3:
    return [[-(units[0] + units[1] + units[2]).normalize()]]
  opposite = -(units[0] + units[1]).normalize()
  if steric == TRIGONAL: return [[opposite]]
  # the two remaining tetrahedral positions: the bisector tilted out of the
  # plane of the known bonds, to either side
  normal = units[0].cross(units[1]).normalize()
  tilt = math.radians(TETRAHEDRAL_ANGLE / 2)
  both = [opposite * math.cos(tilt) + normal * math.sin(tilt) * sign
          for sign in (1, -1)]
  if n_fill == 2: return [both]
  return [[d] for d in both]

def _cone_reference(axis, parent, environment):
  '''
  Start of the rotation about axis: the direction, perpendicular to axis, of
  the environment site nearest to parent that is off the axis; any
  perpendicular (axis.ortho()) without one. Taken from the molecule, so the
  sampled rotations, and the rotation chosen, do not depend on the frame.
  '''
  if len(environment) > 0:
    d = (environment - parent.elems).norms()
    for i in flex.sort_permutation(d):
      v = matrix.col(environment[i]) - parent
      perpendicular = v - axis * v.dot(axis)
      if perpendicular.length() > 0.1: return perpendicular.normalize()
  return axis.ortho().normalize()

def _cone_candidates(unit, steric, n_fill, step_deg, reference=None):
  '''
  Candidate direction lists when 1 neighbour is known: n_fill positions on a
  cone around the axis opposite to it, rotated in step_deg steps over the
  period of the arrangement, starting at reference (a unit vector
  perpendicular to the axis; None: axis.ortho()).
  '''
  axis = -unit
  if reference is None: p = axis.ortho().normalize()
  else:                 p = reference
  q = axis.cross(p)
  if steric == TETRAHEDRAL:
    half_angle, spacing = 180. - TETRAHEDRAL_ANGLE, 120.
  else:
    half_angle, spacing = 180. - TRIGONAL_ANGLE, 180.
  # all positions filled: the arrangement repeats every spacing
  period = spacing if n_fill * spacing == 360. else 360.
  cos_half = math.cos(math.radians(half_angle))
  sin_half = math.sin(math.radians(half_angle))
  candidates = []
  for i in range(int(round(period / step_deg))):
    phi = i * step_deg
    positions = []
    for k in range(n_fill):
      t = math.radians(phi + k * spacing)
      positions.append(axis * cos_half +
                       (p * math.cos(t) + q * math.sin(t)) * sin_half)
    candidates.append(positions)
  return candidates

def fill_directions(parent, known, steric, n_fill, environment, step_deg=10):
  '''
  Unit directions (matrix.col) from parent for n_fill new substituents.
  parent: matrix.col; known: matrix.col sites of the bonded neighbours;
  steric: TETRAHEDRAL, TRIGONAL or LINEAR; n_fill <= steric - len(known);
  environment: flex.vec3_double of sites to keep away from (parent and known
  not in it: at a fixed distance from every rotation they only tie the
  score). Where a rotation or a choice is free, the one maximising the
  minimum distance from an H at 1 A to the environment wins, ties to the
  smallest angle (no environment: the first). Rotations about a single known
  neighbour start at the nearest off-axis environment site
  (_cone_reference()), so the result does not depend on the frame.
  '''
  if len(known) == 0:
    raise ValueError('fill_directions needs at least one known neighbour')
  if not 1 <= n_fill <= steric - len(known):
    raise ValueError('cannot add %d substituents to %d known with steric '
                     'number %d' % (n_fill, len(known), steric))
  units = [(k - parent).normalize() for k in known]
  if len(units) == 1:
    if steric == LINEAR: candidates = [[-units[0]]]
    else:
      candidates = _cone_candidates(units[0], steric, n_fill, step_deg,
        _cone_reference(-units[0], parent, environment))
  elif len(units) <= 3:
    candidates = _fixed_candidates(units, steric, n_fill)
  else:
    raise ValueError('too many known neighbours: %d' % len(units))
  return _best_candidate(candidates, parent, environment)

# ==============================================================================

def cap_names(atom_group, h_names, missing_names, n_fill, parent=None):
  '''
  n_fill names for the caps on an atom of atom_group, unique in its
  conformer (atom_group and the atom groups of the residue sharing its
  altloc, or any for a blank one) and at most 4 characters: the free
  dictionary H names (h_names) in order, then, when parent (the atom) is an
  O, N or S, 'H' + the parent's name (the hydroxyl convention: HO5' on O5',
  HO3B on O3B), then 'H' + each missing heavy atom's name (HCG on a CB
  without CG), then HC1, HC2, ...
  '''
  altloc = atom_group.altloc.strip()
  taken = set()
  for ag in atom_group.parent().atom_groups():
    other = ag.altloc.strip()
    if altloc and other and other != altloc: continue
    taken.update(a.name.strip() for a in ag.atoms())
  candidates = list(h_names)
  if (parent is not None and
      parent.element.strip().upper() in HYDROXYL_LIKE_PARENTS):
    candidates.append('H' + parent.name.strip())
  candidates += ['H' + n for n in missing_names]
  result = []
  for name in candidates:
    if len(result) == n_fill: break
    if len(name) <= 4 and name not in taken:
      result.append(name)
      taken.add(name)
  i = 1
  while len(result) < n_fill:
    name = 'HC%d' % i
    if name not in taken:
      result.append(name)
      taken.add(name)
    i += 1
  return result

def xh_length(element, coordination, use_neutron_distances):
  '''
  X-H length of a cap: the monomer library's (C with 3 substituents takes
  the aromatic C-H), else the sum of the covalent radii.
  '''
  element = element.strip().upper()
  if use_neutron_distances: library = reduce_hydrogen.get_library_neutron_xh()
  else:                     library = reduce_hydrogen.get_library_xray_xh()
  key = element
  if element == 'C' and coordination == TRIGONAL: key = 'C_aromatic'
  if key in library: return library[key]
  radius = reduce_hydrogen._covalent_radius(element)
  if radius is None:
    raise Sorry('Capping: no covalent radius for element %s' % element)
  return radius + H_COVALENT_RADIUS

def _missing_heavy_names(ag, mon_lib_srv, missing_leaving, cache):
  '''
  Heavy atoms of ag's dictionary absent from its conformer, in dictionary
  order. Leaving atoms only if a dangling atom of the residue misses them
  (missing_leaving): an OXT replaced by the peptide bond is not missing.
  '''
  resname = ag.resname.strip()
  if resname not in cache:
    cache[resname] = _dictionary_atoms(resname, mon_lib_srv)
  entries = cache[resname]
  if entries is None: return []
  present = _present_names(ag.atoms()[0], _conformer_names(ag.parent()))
  result = []
  for name, element, is_leaving in entries[0]:
    if element.strip().upper() in ('H', 'D') or name in present: continue
    if is_leaving and name not in missing_leaving: continue
    result.append(name)
  return result

def _report_line(parent, names, missing):
  '''
  For example 'chain A ARG 123 B: CB +3 H (HB2 HB3 HCG); missing CG CD NE CZ
  NH1 NH2'; one space between the fields, a blank altloc left out
  ('ARG 1: CB').
  '''
  ag = parent.parent()
  rg = ag.parent()
  fields = [rg.parent().id.strip(), ag.resname.strip(),
            rg.resseq.strip() + rg.icode.strip(), ag.altloc.strip()]
  line = 'chain %s: %s +%d H (%s)' % (' '.join(f for f in fields if f),
    parent.name.strip(), len(names), ' '.join(names))
  if missing: line += '; missing %s' % ' '.join(missing)
  return line

def _new_h(name, xyz, parent):
  '''
  H atom named name at xyz with the occupancy, B, segid and hetero flag of
  parent; the name padded as reduce_hydrogen pads H names.
  '''
  if len(name) < 4: name = (' ' + name).ljust(4)
  atom = iotbx.pdb.hierarchy.atom()
  atom.set_name(new_name=name)
  atom.set_element(new_element='H')
  atom.set_xyz(new_xyz=xyz)
  atom.set_occ(new_occ=parent.occ)
  atom.set_b(new_b=parent.b)
  atom.set_segid(new_segid=parent.segid)
  atom.set_hetero(new_hetero=parent.hetero)
  return atom

def _model_index(hierarchy):
  '''
  flex.size_t: the index of each atom's MODEL in hierarchy.models(), by
  i_seq (i_seq must be current).
  '''
  result = flex.size_t(hierarchy.atoms_size(), 0)
  for i, m in enumerate(hierarchy.models()):
    result.set_selected(m.atoms().extract_i_seq(), i)
  return result

def _check_clashes(model, caps):
  '''
  Sorry if a cap is closer than CLASH_DISTANCE to any atom of its MODEL but
  its parent, symmetry mates included. caps: [(cap, parent)] after
  reset_i_seq.
  '''
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  model_index = _model_index(hierarchy)
  xrs = hierarchy.extract_xray_structure(
    crystal_symmetry=reduce_hydrogen.contact_symmetry(
      model.crystal_symmetry()))
  pst = xrs.pair_asu_table(distance_cutoff=CLASH_DISTANCE
    ).extract_pair_sym_table()
  parent_of = dict((h.i_seq, parent.i_seq) for h, parent in caps)
  uc = xrs.unit_cell()
  frac = xrs.sites_frac()
  worst = None
  for i, pair_sym_dict in enumerate(pst):
    for j, sym_ops in pair_sym_dict.items():
      if i not in parent_of and j not in parent_of: continue
      if parent_of.get(i) == j or parent_of.get(j) == i: continue
      if model_index[i] != model_index[j]: continue
      if not reduce_hydrogen._compatible(atoms[i], atoms[j]): continue
      for rt_mx in sym_ops:
        d = uc.distance(frac[i], rt_mx * frac[j])
        if worst is None or d < worst[0]: worst = (d, i, j, rt_mx)
  if worst is None: return
  d, i, j, rt_mx = worst
  # the cap first; rt_mx maps j next to i, so swapping inverts it
  if i not in parent_of: i, j, rt_mx = j, i, rt_mx.inverse()
  other = atoms[j].id_str()
  if not rt_mx.is_unit_mx(): other += ' (symmetry mate %s)' % str(rt_mx)
  raise Sorry('Capping H %s clashes with %s (%.2f A)' % (
    atoms[i].id_str(), other, d))

def cap_dangling_atoms(model, result, use_neutron_distances=True,
                       mon_lib_srv=None):
  '''
  Add the capping H of find_dangling_atoms() records to the hierarchy of
  model, the model the records come from (their i_seqs are used). Each H
  goes to its parent's atom group (so its altloc) with the parent's
  occupancy, B, segid and hetero flag; i_seq and serial are reset after.
  Sorry, before anything is placed, if result has internal gaps (a missing
  atom two present dangling atoms of the residue miss: their caps would
  overlap) or hypervalent records (P, or S with 3 or more substituents: an
  H there changes the chemistry); at most 20 of each are listed.
  Directions from fill_directions(), about the parent's bonded neighbours in
  its conformer and away from the sites of its MODEL within
  ENVIRONMENT_RADIUS, the parent and those neighbours excluded; length
  xh_length(); names cap_names(). Sorry if a cap ends up closer than
  CLASH_DISTANCE to an atom other than its parent, symmetry mates included.
  Geometry by steric number (steric_number()): tetrahedral, 120 degrees or
  linear; free rotations sampled in 10 degree steps against the atoms of the
  parent's MODEL within ENVIRONMENT_RADIUS but the parent and its bonded
  neighbours (no symmetry mates; the clash check is symmetry-aware).
  The bond table is result.bond_table when there, else built. mon_lib_srv:
  the monomer library server for the dictionary lookups of the report
  (default model.get_mon_lib_srv(); give the server of the processed model
  the records came from when model is an unprocessed copy).
  Names: free dictionary H names, then 'H' + the parent's name for an O, N
  or S parent, then 'H' + missing heavy atom name, then HC1, HC2 ... The
  result cannot be processed by pdb_interpretation (invented
  names): it is for QM use, and capping is the last change to the hierarchy.
  Returns group_args(added, n_added, n_atoms, n_residues, lines, warnings):
    added: [(parent atom_key, [cap names])], in record order
    lines: one report line per parent (_report_line), the missing atoms
      being the residue's missing heavy atoms in dictionary order, without
      the CCD leaving atoms that a link replaces (OXT, OP3 ...)
  '''
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  gaps = getattr(result, 'internal_gaps', None) or []
  text = [_internal_gap_line(gap) for gap in gaps[:20]]
  if len(gaps) > 20: text.append('... and %d more' % (len(gaps) - 20))
  hyper = getattr(result, 'hypervalent', None) or []
  text += [_hypervalent_line(r, atoms[r.i_seq]) for r in hyper[:20]]
  if len(hyper) > 20: text.append('... and %d more' % (len(hyper) - 20))
  if text: raise Sorry('\n'.join(text))
  sites = model.get_sites_cart()
  table = getattr(result, 'bond_table', None)
  if table is None: table = bond_table(model)
  model_index = _model_index(hierarchy)
  if mon_lib_srv is None: mon_lib_srv = model.get_mon_lib_srv()
  properties = atom_property()
  is_metal = lambda element: properties.is_metal(element.strip().upper())
  # leaving atoms missed by a dangling atom, per residue group
  missing_leaving = {}
  for record in result.records:
    rg = atoms[record.i_seq].parent().parent()
    missing_leaving.setdefault(rg.memory_id(), set()).update(
      record.missing_names)
  # the residues' missing atoms, before any cap is appended
  missing_by_record, cache = [], {}
  for record in result.records:
    parent = atoms[record.i_seq]
    if atom_key(parent) != record.key:
      raise Sorry('Capping: the records do not belong to this model '
        '(%s at i_seq %d)' % (_atom_label(parent), record.i_seq))
    if record.source == 'dictionary':
      missing_by_record.append(_missing_heavy_names(parent.parent(),
        mon_lib_srv, missing_leaving[parent.parent().parent().memory_id()],
        cache))
    else:
      missing_by_record.append([])
  added, lines, caps, residues = [], [], [], set()
  for record, missing in zip(result.records, missing_by_record):
    parent = atoms[record.i_seq]
    rg = parent.parent().parent()
    # the neighbours find_dangling_atoms counted
    j_seqs = table.get(record.i_seq, ())
    if record.source == 'guess':
      known = _guess_neighbours(parent, _has_internal_bonds(rg, table),
        j_seqs, atoms, sites, is_metal)
    else:
      known = _conformer_neighbours(parent, j_seqs, atoms, is_metal)
    if len(known) + record.n_fill != record.coordination:
      raise RuntimeError('Capping: %s, key %s: %d neighbours + %d H when '
        'building, coordination %d when detecting' % (_atom_label(parent),
        str(record.key), len(known), record.n_fill, record.coordination))
    center = matrix.col(parent.xyz)
    # the parent's MODEL, without the parent and its bonded neighbours: at a
    # fixed distance from every rotation, they would only tie the score
    d = (sites - parent.xyz).norms()
    near = ((d < ENVIRONMENT_RADIUS) & (d > 0.01) &
            (model_index == model_index[record.i_seq]))
    near[record.i_seq] = False
    for j in j_seqs: near[j] = False
    for a in known: near[a.i_seq] = False
    environment = sites.select(near)
    element = reduce_hydrogen._atom_element(parent)
    directions = fill_directions(center, [matrix.col(a.xyz) for a in known],
      steric_number(element, record.coordination), record.n_fill,
      environment)
    length = xh_length(element, record.coordination, use_neutron_distances)
    names = cap_names(parent.parent(), record.h_names, record.missing_names,
                      record.n_fill, parent=parent)
    for name, direction in zip(names, directions):
      h = _new_h(name, (center + direction * length).elems, parent)
      parent.parent().append_atom(h)
      caps.append((h, parent))
    added.append((record.key, names))
    lines.append(_report_line(parent, names, missing))
    residues.add(rg.memory_id())
  hierarchy.atoms().reset_i_seq()
  hierarchy.atoms().reset_serial()
  if caps: _check_clashes(model, caps)
  return group_args(
    added      = added,
    n_added    = len(caps),
    n_atoms    = len(added),
    n_residues = len(residues),
    lines      = lines,
    warnings   = [])

def cap_model(model, guess_for_unknowns=False, cut_bonds=None,
              use_neutron_distances=True, calculate_charge=True,
              dangling=None, mon_lib_srv=None, log=None):
  '''
  The capping sequence of reduce2 capping=True on model: find_dangling_atoms
  (unless dangling, a result for a model with the same atoms, is given: then
  model may be unprocessed), Sorry on residues without a dictionary unless
  guess_for_unknowns, the total charge of the capped molecule from
  electrons.run(model=, extra_sigma_bonds=) before anything is placed
  (calculate_charge; needs a processed model) with the open-shell warnings,
  then cap_dangling_atoms on the model's hierarchy in place.
  Returns group_args(dangling, capped, total_charge, formal_charges,
  warnings): dangling the find_dangling_atoms result used, capped the
  cap_dangling_atoms result, total_charge and formal_charges (the electrons
  distribution) None when not calculated, warnings the
  'open-shell cap on <atom>' lines.
  '''
  from mmtbx.ligands import electrons
  if dangling is None:
    dangling = find_dangling_atoms(model, guess_for_unknowns=guess_for_unknowns,
                                   log=log, cut_bonds=cut_bonds)
  if dangling.unknown_residues and not guess_for_unknowns:
    raise Sorry('capping=True: no dictionary for %s; set '
      'capping_guess_for_unknowns=True to cap them by element valence'
      % ', '.join(dangling.unknown_residues))
  total_charge, formal_charges, warnings = None, None, []
  if calculate_charge:
    extra = extra_sigma_bonds(dangling)
    formal_charges = electrons.run(
      model                 = model,
      extra_sigma_bonds     = extra,
      return_formal_charges = True)
    total_charge = formal_charges.get_total_charge()
    atoms = model.get_hierarchy().atoms()
    for i_seq in sorted(extra):
      if formal_charges.get(i_seq, 0):
        warnings.append('open-shell cap on %s' % atoms[i_seq].id_str())
  capped = cap_dangling_atoms(model, dangling,
    use_neutron_distances=use_neutron_distances, mon_lib_srv=mon_lib_srv)
  return group_args(
    dangling       = dangling,
    capped         = capped,
    total_charge   = total_charge,
    formal_charges = formal_charges,
    warnings       = warnings)
