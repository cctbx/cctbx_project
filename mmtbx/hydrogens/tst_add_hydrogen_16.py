from __future__ import absolute_import, division, print_function
import iotbx.pdb, iotbx.cif
import mmtbx.model
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.monomer_library import server
from mmtbx.hydrogens import reduce_hydrogen

def run():
  test_001()
  test_002()
  test_003()
  test_004()
  test_005()
  test_006()
  test_007()
  test_008()
  test_009()
  test_010()
  test_011()
  test_012()
  test_013()
  test_014()
  test_015()
  test_016()
  test_017()
  test_018()
  test_019()
  test_020()
  test_021()
  test_022()
  test_023()

# ------------------------------------------------------------------------------

def place_atoms(pdb_str):
  '''Place H (no optimization); return all atoms. Names repeat between
  residues, so only single-residue models may be keyed by name (place()).'''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  return list(obj.get_model().get_hierarchy().atoms())

def place(pdb_str):
  '''Place H in a one-residue model; return the atoms by name.'''
  return {a.name.strip(): a for a in place_atoms(pdb_str)}

def h_names_on(pdb_str, parent_name):
  '''Place H (no optimization); return names of the H bonded to parent_name.'''
  atoms = place_atoms(pdb_str)
  parent = [a for a in atoms if a.name.strip() == parent_name][0]
  return sorted(a.name.strip() for a in atoms
    if a.element.strip() == 'H' and a.distance(parent) < 1.2)

def test_001():
  '''
    A restraint file with several ligands: each ligand's CH2 names are checked
    against its own ideal coordinates. AC3 and AD0 (codes not in the CCD) share
    atom names, with H21/H22 swapped between them. The AD0 H already match AD0,
    so name_prochiral_h must leave them; it used to read AC3, the first block.
  '''
  file_name = "tst_add_hydrogen_16.cif"
  with open(file_name, "w") as f:
    f.write(cif_str_001)
  mon_lib_srv = server.server()
  mon_lib_srv.process_cif_object(
    cif_object = iotbx.cif.reader(input_string=cif_str_001).model(),
    file_name  = file_name)
  h = iotbx.pdb.input(lines=pdb_str_001.split("\n"), source_info=None).construct_hierarchy()
  h.atoms().reset_i_seq()
  names = list(h.atoms().extract_name())
  reduce_hydrogen.name_prochiral_h(h, mon_lib_srv)
  assert list(h.atoms().extract_name()) == names, (
    names, list(h.atoms().extract_name()))

def get_user_lig_srv():
  '''User restraint file for LIG, a code the CCD uses for another molecule.'''
  file_name = "tst_add_hydrogen_16_lig.cif"
  with open(file_name, "w") as f:
    f.write(cif_str_002)
  mon_lib_srv = server.server()
  mon_lib_srv.process_cif_object(
    cif_object = iotbx.cif.reader(input_string=cif_str_002).model(),
    file_name  = file_name)
  return mon_lib_srv

def test_002():
  '''
    A user restraint file wins over an unrelated CCD entry with the same code.
    CCD LIG (C15H11N3) has a CH2 at C17 with the same names as this user LIG
    (propane C5-C17-C20) but the opposite handedness; the H match the user file,
    so name_prochiral_h must leave them (it used CCD geometry and swapped them).
  '''
  mon_lib_srv = get_user_lig_srv()
  h = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None).construct_hierarchy()
  h.atoms().reset_i_seq()
  names = list(h.atoms().extract_name())
  reduce_hydrogen.name_prochiral_h(h, mon_lib_srv)
  assert list(h.atoms().extract_name()) == names, (
    names, list(h.atoms().extract_name()))

def test_003():
  '''
    Bond orders come from the user file only: the user LIG has single bonds, CCD
    LIG's double bonds (C5=C6, C20=C21, ...) must not leak in.
  '''
  mon_lib_srv = get_user_lig_srv()
  orders = reduce_hydrogen._bond_orders("LIG", mon_lib_srv, {})
  assert orders == {}, sorted(tuple(sorted(k)) for k in orders)

def test_004():
  '''
    6EL (5yj1) passes test_for_peptide (N, CA, C, O) but its H2 sits on ring
    CG, not on N. The N-terminal H2 filter removed it by name: CG got 1 H.
  '''
  names = h_names_on(pdb_str_004, 'CG')
  assert names == ['H2', 'H3'], names

def test_005():
  '''
    Control: 0TD H2 is the amine H on N; a residue that is not first in the
    chain keeps only H on N.
  '''
  names = h_names_on(pdb_str_005, 'N')
  assert names == ['H'], names

def test_006():
  '''
    GLZ (aminoacetaldehyde) passes test_for_peptide; its HXT is the aldehyde
    H on C (no OXT). The C-terminal HXT filter removed it by name.
  '''
  names = h_names_on(pdb_str_006, 'C')
  assert names == ['HXT'], names

def check_h(atoms, h, a0, a1, angle_ideal):
  '''H sits on a0 at a plausible X-H length and H-a0-a1 angle.'''
  d = atoms[h].distance(atoms[a0])
  assert 0.8 < d < 1.15, (h, d)
  angle = atoms[a0].angle(atoms[h], atoms[a1], deg=True)
  assert approx_equal(angle, angle_ideal, eps=1.0), (h, angle)

def test_007():
  '''
    PEO (H2O2): nothing beyond the O-O bond anchors the H dihedral, so riding H
    cannot parameterize HO1/HO2 and they were deleted (1ng4).
  '''
  atoms = place(pdb_str_007)
  assert sorted(n for n in atoms if n.startswith('H')) == ['HO1', 'HO2'], list(atoms)
  check_h(atoms, 'HO1', 'O1', 'O2', 100.92)
  check_h(atoms, 'HO2', 'O2', 'O1', 100.89)

def test_008():
  '''
    MOH (methanol): same for the methyl and OH H. EOH (ethanol, anchored by the
    third heavy atom) is the control.
  '''
  atoms = place(pdb_str_008)
  hs = sorted(n for n in atoms if n.startswith('H'))
  assert hs == ['H1', 'H2', 'H3', 'HO'], hs
  for h, angle_ideal in [('H1', 107.79), ('H2', 112.98), ('H3', 112.94)]:
    check_h(atoms, h, 'C', 'O', angle_ideal)
  check_h(atoms, 'HO', 'O', 'C', 108.03)
  for h1, h2 in [('H1', 'H2'), ('H1', 'H3'), ('H2', 'H3')]:
    assert atoms[h1].distance(atoms[h2]) > 1.4, (h1, h2)
  hs = sorted(n for n in place(pdb_str_008_control) if n.startswith('H'))
  assert hs == ['H11', 'H12', 'H21', 'H22', 'H23', 'HO'], hs

def h_on(atoms, resname, resseq, parent_name):
  '''Names of the H bonded to <parent_name> of one residue.'''
  parent = [a for a in atoms
    if a.name.strip() == parent_name
    and a.parent().resname.strip() == resname
    and a.parent().parent().resseq.strip() == resseq]
  assert len(parent) == 1, (resname, resseq, parent_name)
  return sorted(a.name.strip() for a in atoms
    if a.element.strip() == 'H' and a.distance(parent[0]) < 1.3)

def test_009():
  '''
    5nxq: 9FZ and 9G2 are in the polypeptide chain, but their dictionaries are
    the free amino acid (N is NH2: H3/H4 and H6/H7). Both H stayed, so N had
    four bonds. The terminal-H filter only knew the name H2.
  '''
  atoms = place_atoms(pdb_str_009)
  assert len(h_on(atoms, '9FZ', '4', 'N')) == 1, h_on(atoms, '9FZ', '4', 'N')
  assert len(h_on(atoms, '9G2', '10', 'N')) == 1, h_on(atoms, '9G2', '10', 'N')
  # the amino acids in the same chain are unchanged
  assert h_on(atoms, 'ILE', '5', 'N') == ['H'], h_on(atoms, 'ILE', '5', 'N')
  assert h_on(atoms, 'GLU', '11', 'N') == ['H'], h_on(atoms, 'GLU', '11', 'N')
  # the C9-C10 cross-link between the two ligands still removes one methyl H
  assert len(h_on(atoms, '9FZ', '4', 'C9')) == 2, h_on(atoms, '9FZ', '4', 'C9')
  assert len(h_on(atoms, '9G2', '10', 'C10')) == 2, h_on(atoms, '9G2', '10', 'C10')

def atom_of(atoms, resname, resseq, name):
  '''The one atom of that residue.'''
  found = [a for a in atoms if a.name.strip() == name
    and a.parent().resname.strip() == resname
    and a.parent().parent().resseq.strip() == resseq]
  assert len(found) == 1, (resname, resseq, name)
  return found[0]

def test_013():
  '''
    5nxq again: the amide H of a ligand in the chain. The TRANS link defines
    C-N-H by name, so a residue whose H is called anything else (9FZ H3, 9G2
    H6) got no such restraint, riding called it a rotatable amine H and it
    ended up 17 deg from the preceding C. The amide H now sits where ILE 5's
    does (the control): the three angles at N add up to 360 and it is trans to
    the preceding O. CA-N-H comes from the link as well, since the dictionary
    of the free amino acid gives N the sp3 amine value (9FZ 108.5 deg).
  '''
  atoms = place_atoms(pdb_str_013)
  for resname, resseq, h_name, prev_resname, prev_resseq in (
      ('9FZ', '4',  'H3', 'ILE', '3'),
      ('9G2', '10', 'H6', 'LEU', '9'),
      ('ILE', '5',  'H',  '9FZ', '4')):
    n  = atom_of(atoms, resname, resseq, 'N')
    h  = atom_of(atoms, resname, resseq, h_name)
    ca = atom_of(atoms, resname, resseq, 'CA')
    c  = atom_of(atoms, prev_resname, prev_resseq, 'C')
    o  = atom_of(atoms, prev_resname, prev_resseq, 'O')
    why = (resname, h_name)
    assert approx_equal(n.angle(c, h, deg=True), 124.3, eps=1.0), why
    assert approx_equal(n.angle(ca, h, deg=True), 114.0, eps=1.0), why
    assert approx_equal(n.angle(c, h, deg=True) + n.angle(ca, h, deg=True)
                        + n.angle(c, ca, deg=True), 360.0, eps=1.0), why
    d = matrix.col(h.xyz) - matrix.col(n.xyz)
    plane = ((matrix.col(c.xyz) - matrix.col(n.xyz)).cross(
              matrix.col(ca.xyz) - matrix.col(n.xyz))).normalize()
    assert abs(plane.dot(d)) < 0.1, (why, 'H is out of the amide plane')
    assert 0.8 < n.distance(h) < 1.1, (why, n.distance(h))

def test_014():
  '''
    The v2->v3 H rename (X1 -> X3 when the dictionary has X1 and X2) is right
    only where pdb_interpretation maps the v3 name back onto the dictionary:
    standard and D-amino acids. A modified amino acid's dictionary already
    follows the CCD, so the renamed H had no restraint and was dropped
    (2hi2 MEA CB kept HB2 only; 9eor HT7 CA kept HC2 only; its HB1 and HB2
    are not even on the same atom). LEU is the control.
  '''
  atoms = place_atoms(pdb_str_014)
  assert h_on(atoms, 'MEA', '1', 'CB') == ['HB1', 'HB2'], h_on(atoms, 'MEA', '1', 'CB')
  assert h_on(atoms, 'LEU', '3', 'CB') == ['HB2', 'HB3'], h_on(atoms, 'LEU', '3', 'CB')
  atoms = place_atoms(pdb_str_014_ht7)
  assert h_on(atoms, 'HT7', '3', 'CA') == ['HB1', 'HC2'], h_on(atoms, 'HT7', '3', 'CA')
  assert h_on(atoms, 'HT7', '3', 'CB') == ['HA'], h_on(atoms, 'HT7', '3', 'CB')

def test_015():
  '''
  2e4i: BGM 14 -> BGM 15, two nucleotides whose H come from the CCD (no
  geostd entry). HO3' of BGM 14 was removed only after the riding H setup,
  which stopped with a Sorry: no H-O3'-P angle, and P of BGM 15 is the
  neighbour of an H parent (HOP2 on OP2). A BGM next to a standard
  nucleotide passed. The 3' end BGM 15 keeps its hydroxyl H.
  '''
  atoms = place_atoms(pdb_str_015)
  assert h_on(atoms, 'BGM', '14', "O3'") == [], h_on(atoms, 'BGM', '14', "O3'")
  assert h_on(atoms, 'BGM', '15', "O3'") == ["HO3'"], h_on(
    atoms, 'BGM', '15', "O3'")

def test_016():
  '''
  The CCD free nucleotide gives BGM an H on OP2. Inside a chain the phosphate
  is a diester (pKa ~1), so the H goes once P bridges to another residue:
  BGM 15 (O3' of BGM 14 + its own O5'). BGM 14 is the 5' end of the fragment;
  its P bridges only through O5', so it keeps HOP2. Free 35G (cGMP, geostd)
  is a diester within one residue and keeps its HOP2 as well.
  '''
  atoms = place_atoms(pdb_str_015)
  assert h_on(atoms, 'BGM', '15', 'OP1') == []
  assert h_on(atoms, 'BGM', '15', 'OP2') == [], h_on(atoms, 'BGM', '15', 'OP2')
  assert h_on(atoms, 'BGM', '14', 'OP2') == ['HOP2'], h_on(
    atoms, 'BGM', '14', 'OP2')
  atoms = place_atoms(pdb_str_016)
  assert h_on(atoms, '35G', '1', 'O2P') == ['HOP2'], h_on(atoms, '35G', '1', 'O2P')

def test_017():
  '''
  BGM has no geostd entry; its restraints come from the CCD, which has no
  energy types. Without a type the probe atom info had no radius and no
  acceptor: every BGM atom gave "Could not find atom info", the element
  fallback made polar H 1.22 A, and no O was an acceptor. The types now come
  from element + the H the CCD puts on the atom.
  '''
  from libtbx import group_args
  from mmtbx.probe import Helpers
  md = reduce_hydrogen.get_h_restraints('BGM', strict=False)
  types = dict((a.atom_id, a.type_energy) for a in md.atom_list)
  expected = {'O6': 'O', 'OP1': 'O', 'OP2': 'OH1', 'N1': 'NH1', 'N2': 'NH2',
    'N7': 'N', 'H1': 'HNH1', 'HOP2': 'HOH1', "H1'": 'HCH', 'BR': 'BR', 'P': 'P'}
  for name, t in expected.items():
    assert types[name] == t, (name, types[name], t)
  known = server.ener_lib().lib_atom
  assert all(t in known for t in types.values()), types
  # an element ener_lib does not know gets no type, not a bad one
  cc_cif = {
    '_chem_comp_atom': [group_args(atom_id=n, type_symbol=e) for n, e in
      (('X1', 'XX'), ('S1', 'S'), ('HS1', 'H'))],
    '_chem_comp_bond': [group_args(atom_id_1='X1', atom_id_2='S1'),
                        group_args(atom_id_1='S1', atom_id_2='HS1')]}
  assert reduce_hydrogen.get_ccd_energy_types(cc_cif) == {
    'X1': None, 'S1': 'SH1', 'HS1': 'HSH1'}
  # probe atom info on the placed fragment
  pdb_inp = iotbx.pdb.input(lines=pdb_str_015.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  model = obj.get_model()
  bps, asu = model.get_restraints_manager().geometry.get_all_bond_proxies(
    sites_cart = model.get_sites_cart())
  neighbors = Helpers.getBondedNeighborLists(model.get_atoms(), bps, asu)
  from mmtbx.hydrogens import default_probe_phil
  ret = Helpers.getExtraAtomInfo(model=model, bondedNeighborLists=neighbors,
    probePhil=default_probe_phil())
  assert 'Could not find atom info' not in ret.warnings, ret.warnings
  info = dict(((a.parent().parent().resseq.strip(), a.name.strip()),
    ret.extraAtomInfo.getMappingFor(a)) for a in model.get_atoms())
  for name in ('O6', 'OP1', "O4'"):
    assert info[('15', name)].isAcceptor, name
  assert approx_equal(info[('15', 'H1')].vdwRadius, 1.05)
  assert approx_equal(info[('15', 'C2')].vdwRadius, 1.70)

def test_018():
  '''
  7bev (20 NMR models): model.process() builds, by default, a 5 A pair table
  for ADP restraints that ignores model boundaries. Superposed models made it
  6.5 GB; placement built it twice (process, then again after a select).
  Neither place_hydrogens nor reduce2's reinterpretation (approach=optimize,
  flipkins) needs it. Two models of a 7bev fragment.
  '''
  import iotbx.phil
  from mmtbx.programs import reduce2
  pdb_inp = iotbx.pdb.input(lines=pdb_str_018.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  model = obj.get_model()
  grm = model.get_restraints_manager().geometry
  grm.pair_proxies(sites_cart=model.get_sites_cart())
  assert grm.plain_pairs_radius is None
  assert grm.plain_pair_sym_table is None
  n_h = [sum(1 for a in m.atoms() if a.element.strip() == 'H')
    for m in model.get_hierarchy().models()]
  assert len(n_h) == 2 and n_h[0] == n_h[1] > 0, n_h
  # reduce2's reinterpretation with restraints
  params = iotbx.phil.parse(reduce2.master_phil_str,
    process_includes=True).extract()
  task = reduce2.Program(data_manager=None, params=params, logger=null_out())
  task.model = model
  task._ReinterpretModel(make_restraints=True)
  grm = task.model.get_restraints_manager().geometry
  grm.pair_proxies(sites_cart=task.model.get_sites_cart())
  assert grm.plain_pair_sym_table is None

def test_019():
  '''
  2fd6 ASN U 172 - NAG C 1: the NAG-ASN link gives an angle to C1 only for
  HD21 and deletes HD22, which went only after riding H. The deposited NAG
  still has O1, so C1 is the second neighbour of HO1, and riding H stopped:
  "angle restraints involving an H atom are missing". HD22 now goes first.
  1a7s ASN 114 - NAG 401 (no O1) is the control: HD21 only, as before.
  '''
  atoms = place_atoms(pdb_str_019)
  # HO1 of the misplaced O1 (0.79 A from ND2) is close too
  hd = [n for n in h_on(atoms, 'ASN', '172', 'ND2') if n.startswith('HD')]
  assert hd == ['HD21'], hd
  atoms = place_atoms(pdb_str_019_control)
  assert h_on(atoms, 'ASN', '114', 'ND2') == ['HD21'], h_on(
    atoms, 'ASN', '114', 'ND2')

def test_020():
  '''
  7a5m ALA A 103: N, C and O have no altloc, CA and CB are in A and B, SER
  A 102 is all A/B. Each altloc gets its amide H. add_link_h_restraints gave
  H/A the angle C/B(102)-N-H as well (and H/B C/A), and riding H stopped:
  "angle restraints involving an H atom are missing ... N ALA A 103 / CA A".
  Each H must lie in the amide plane of its own altloc.
  '''
  atoms = place_atoms(pdb_str_020)
  a = dict(((x.parent().parent().resseq.strip(), x.parent().altloc,
    x.name.strip()), matrix.col(x.xyz)) for x in atoms)
  n = a[('103', '', 'N')]
  angle = lambda x, z: (x-n).angle(z-n, deg=True)
  for alt in 'AB':
    h, c, ca = a[('103', alt, 'H')], a[('102', alt, 'C')], a[('103', alt, 'CA')]
    assert approx_equal(abs(h-n), 0.86, eps=0.01)
    # planar: the three angles around N add up to 360
    assert approx_equal(angle(c, h) + angle(ca, h) + angle(c, ca), 360,
      eps=1), (alt, angle(c, h), angle(ca, h), angle(c, ca))

def test_021():
  '''
  An atom name the dictionary does not know (here ASP CG misnamed CGX; in
  7a5m the backbone of 2L5/2L6, whose geostd names predate the CCD) is
  reported as such; the H next to it used to be listed under "incomplete
  side chains".
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_021.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  counts = obj.get_counts()
  assert counts.residues_unknown_atoms == [('ASP A 11', ['CGX'])], \
    counts.residues_unknown_atoms
  assert sorted(l.split()[0] for l in counts.site_labels_unknown_atoms) == \
    ['HB2', 'HB3'], counts.site_labels_unknown_atoms
  assert counts.residues_missing_neighbor == [], counts.residues_missing_neighbor

def test_022():
  '''
  1evd is a cyclic peptide: ACB 1 N is bonded to MDH 5 C, so it is no
  N-terminus. It got an NH3 propeller (mod_NH3 then made H2 conflict, and
  after that all three H were dropped). Now it keeps one amide H in the
  C-N-CA plane; the link has no C-N-CA angle, which riding H needs, so
  add_link_h_restraints adds it. 1fb9 (two models) is the control: the
  N-terminus keeps H1, H2, H3 in every model.
  '''
  atoms = place_atoms(pdb_str_022)
  a = dict(((x.parent().resname.strip(), x.name.strip()), x) for x in atoms)
  n = a[('ACB', 'N')]
  hs = [x for x in atoms if x.element.strip() == 'H' and x.distance(n) < 1.2]
  assert [x.name.strip() for x in hs] == ['H'], [x.name.strip() for x in hs]
  nv = matrix.col(n.xyz)
  v = lambda x: matrix.col(x.xyz) - nv
  h, c, ca = v(hs[0]), v(a[('MDH', 'C')]), v(a[('ACB', 'CA')])
  total = h.angle(c, deg=True) + h.angle(ca, deg=True) + c.angle(ca, deg=True)
  assert approx_equal(total, 360, eps=1), total
  atoms = place_atoms(pdb_str_022_control)
  for m in ('   1', '   2'):
    n = [x for x in atoms if x.name.strip() == 'N'
         and x.parent().parent().resseq.strip() == '1'
         and x.parent().parent().parent().parent().id == m][0]
    names = sorted(x.name.strip() for x in atoms if x.element.strip() == 'H'
      and x.distance(n) < 1.2 and x.parent().parent().parent().parent().id == m)
    assert names == ['H1', 'H2', 'H3'], (m, names)

def test_023():
  '''
  4dco VAL C 504 - ASJ C 505: geostd ASJ calls a C H "H" (the amide H on N
  are HN, HNA). The peptide link's C-N-H angle and plane are matched by name,
  so they took the C H, and riding H stopped: "angle and bond restraints are
  conflicting". pdb_interpretation now skips a link restraint on an H that
  the residue's dictionary bonds elsewhere; the amide H gets its angle from
  add_link_h_restraints and lies in the C-N-CA plane.
  ASJ has no C=O, so test_for_peptide missed it and N kept both dictionary H
  (superposed). The control, ASJ without the VAL, keeps both.
  '''
  atoms = place_atoms(pdb_str_023)
  a = dict(((x.parent().resname.strip(), x.name.strip()), x) for x in atoms)
  n, c = a[('ASJ', 'N')], a[('ASJ', 'C')]
  hs = [x for x in atoms if x.element.strip() == 'H' and x.distance(n) < 1.2]
  assert len(hs) == 1, [x.name.strip() for x in hs]
  nv = matrix.col(n.xyz)
  v = lambda x: matrix.col(x.xyz) - nv
  h, cv, ca = v(hs[0]), v(a[('VAL', 'C')]), v(a[('ASJ', 'CA')])
  total = h.angle(cv, deg=True) + h.angle(ca, deg=True) + cv.angle(ca, deg=True)
  assert approx_equal(total, 360, eps=1), total
  assert sorted(x.name.strip() for x in atoms if x.element.strip() == 'H'
    and x.distance(c) < 1.2) == ['H', 'HA']
  free = "\n".join(l for l in pdb_str_023.split("\n") if 'VAL' not in l)
  atoms = place_atoms(free)
  n = [x for x in atoms if x.name.strip() == 'N'][0]
  hs = [x for x in atoms if x.element.strip() == 'H' and x.distance(n) < 1.2]
  assert sorted(x.name.strip() for x in hs) == ['HN', 'HNA'], \
    [x.name.strip() for x in hs]
  assert hs[0].distance(hs[1]) > 1.4

def hand(atoms, resname, resseq, parent, first, second, third):
  '''
  Sign of the chiral volume (first, second, third) at parent. For a CH2 pass
  the two heavy neighbours and the '2' hydrogen: which side that H is named on.
  For a propeller pass two H and the heavy neighbour: the turning sense of the
  names around the axis.
  '''
  a = dict((x.name.strip(), matrix.col(x.xyz)) for x in atoms
    if x.parent().resname.strip() == resname
    and x.parent().parent().resseq.strip() == resseq)
  p = a[parent]
  return 1 if ((a[first]-p).cross(a[second]-p)).dot(a[third]-p) > 0 else -1

def test_010():
  '''
    LEU HB2/HB3 came out swapped (7IT7 fragment). The CH2 reference read the
    CCD's ideal coordinates, and for ARG CB/CG, ILE CG1, LEU CB and MET CB/CG
    those contradict the same entry's model coordinates. The model coordinates
    are what reduce and deposited models follow; all three groups here are -1.
    GLU CB/CG (the CCD agrees with itself) are the control.
  '''
  atoms = place_atoms(pdb_str_010)
  assert hand(atoms, 'LEU', '107', 'CB', 'CA', 'CG', 'HB2') == -1
  assert hand(atoms, 'GLU', '41', 'CB', 'CA', 'CG', 'HB2') == -1
  assert hand(atoms, 'GLU', '41', 'CG', 'CB', 'CD', 'HG2') == -1

def test_011():
  '''
    A propeller came out as the mirror of the CCD every time: the H are still
    superposed when check_propeller_order looks at them, so the order came from
    the riding frame (HB1 got n=2, HB3 n=0). CH3 and NH3 alike; the CCD's two
    coordinate sets agree for 94% of propellers, and here both say -1.
  '''
  atoms = place_atoms(pdb_str_011)
  assert hand(atoms, 'ALA', '42', 'CB', 'HB1', 'HB2', 'CA') == -1
  assert hand(atoms, 'LYS', '87', 'NZ', 'HZ1', 'HZ2', 'CE') == -1
  atoms = place_atoms(pdb_str_010)
  assert hand(atoms, 'LEU', '107', 'CD1', 'HD11', 'HD12', 'CG') == -1
  assert hand(atoms, 'LEU', '107', 'CD2', 'HD21', 'HD22', 'CG') == -1
  # the CH2 of the same residue is unchanged
  assert hand(atoms, 'LEU', '107', 'CB', 'CA', 'CG', 'HB2') == -1

def test_012():
  '''
    AMP (5msd) is a free nucleotide: its O3' is a hydroxyl, but HO3' was
    dropped from every residue that tests as RNA/DNA, where O3' usually
    carries the next phosphate. The H is placed now and removed only where
    the O really is esterified (exclude_H_on_esterified_O).

    The DNA control keeps no HO3' anywhere, including the 3'-terminal DC 3:
    the monomer library entry of a polymer nucleotide has no HO3' to place.
    Giving a 3' end its hydroxyl H is a separate gap, in the library rather
    than here - the RNA entries A/U/G/C have no HO3' either.
  '''
  atoms = place_atoms(pdb_str_012)
  assert h_on(atoms, 'AMP', '1201', "O3'") == ["HO3'"], h_on(
    atoms, 'AMP', '1201', "O3'")
  assert h_on(atoms, 'AMP', '1201', "O2'") == ["HO2'"]
  atoms = place_atoms(pdb_str_012_control)
  for resname, resseq in (('DC', '1'), ('DG', '2'), ('DC', '3')):
    assert h_on(atoms, resname, resseq, "O3'") == [], (resname, resseq)
  # 5GP has the same ligand-style dictionary as AMP, but here its O3' carries
  # the phosphate of C 2 (4TNA), so the H is placed and then removed again
  atoms = place_atoms(pdb_str_012_linked)
  assert h_on(atoms, '5GP', '1', "O3'") == [], h_on(atoms, '5GP', '1', "O3'")
  assert h_on(atoms, '5GP', '1', "O2'") == ["HO2'"]

# ------------------------------------------------------------------------------

cif_str_001 = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
AC3 AC3 'ethanol-like, first ligand' ligand 9 3 .
AD0 AD0 'ethanol-like, second ligand' ligand 9 3 .

data_comp_AC3
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   AC3 C1   C   1.5200   0.0000   0.0000
   AC3 C2   C   0.0000   0.0000   0.0000
   AC3 O3   O  -0.4773   1.3480   0.0000
   AC3 H11  H   1.8641   0.8604   0.2868
   AC3 H12  H   1.8678  -0.5797   0.6956
   AC3 H13  H   1.8641  -0.2868  -0.8604
   AC3 H21  H  -0.3231  -0.4572  -0.7921
   AC3 H22  H  -0.3231  -0.4572   0.7921
   AC3 H3   H  -1.2837   1.6168   0.0000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   AC3 C1   C2   single 1.520 0.020
   AC3 C2   O3   single 1.430 0.020
   AC3 C1   H11  single 0.970 0.020
   AC3 C1   H12  single 0.970 0.020
   AC3 C1   H13  single 0.970 0.020
   AC3 C2   H21  single 0.970 0.020
   AC3 C2   H22  single 0.970 0.020
   AC3 O3   H3   single 0.850 0.020

data_comp_AD0
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   AD0 C1   C   1.5200   0.0000   0.0000
   AD0 C2   C   0.0000   0.0000   0.0000
   AD0 O3   O  -0.4773   1.3480   0.0000
   AD0 H11  H   1.8641   0.8604   0.2868
   AD0 H12  H   1.8678  -0.5797   0.6956
   AD0 H13  H   1.8641  -0.2868  -0.8604
   AD0 H21  H  -0.3231  -0.4572   0.7921
   AD0 H22  H  -0.3231  -0.4572  -0.7921
   AD0 H3   H  -1.2837   1.6168   0.0000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   AD0 C1   C2   single 1.520 0.020
   AD0 C2   O3   single 1.430 0.020
   AD0 C1   H11  single 0.970 0.020
   AD0 C1   H12  single 0.970 0.020
   AD0 C1   H13  single 0.970 0.020
   AD0 C2   H21  single 0.970 0.020
   AD0 C2   H22  single 0.970 0.020
   AD0 O3   H3   single 0.850 0.020
"""

pdb_str_001 = """
HETATM    1  C1  AD0 A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    2  C2  AD0 A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  O3  AD0 A   1      -0.477   1.348   0.000  1.00 20.00           O
HETATM    4  H21 AD0 A   1      -0.323  -0.457   0.792  1.00 20.00           H
HETATM    5  H22 AD0 A   1      -0.323  -0.457  -0.792  1.00 20.00           H
END
"""

cif_str_002 = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG LIG 'user ligand, propane' ligand 11 3 .

data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   LIG C5   C  -0.5694   1.4093   0.0000
   LIG C17  C   0.0000   0.0000   0.0000
   LIG C20  C   1.5200   0.0000   0.0000
   LIG H171 H  -0.3131  -0.4641   0.7921
   LIG H172 H  -0.3131  -0.4641  -0.7921
   LIG H5A  H   0.1571   2.0521   0.0000
   LIG H5B  H  -1.1146   1.5383  -0.7919
   LIG H5C  H  -1.1146   1.5383   0.7919
   LIG H20A H   1.8438  -0.9144   0.0000
   LIG H20B H   1.8438   0.4572  -0.7919
   LIG H20C H   1.8438   0.4572   0.7919
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   LIG C5   C17  single 1.520 0.020
   LIG C17  C20  single 1.520 0.020
   LIG C17  H171 single 0.970 0.020
   LIG C17  H172 single 0.970 0.020
   LIG C5   H5A  single 0.970 0.020
   LIG C5   H5B  single 0.970 0.020
   LIG C5   H5C  single 0.970 0.020
   LIG C20  H20A single 0.970 0.020
   LIG C20  H20B single 0.970 0.020
   LIG C20  H20C single 0.970 0.020
"""

pdb_str_002 = """
HETATM    1  C5  LIG A   1      -0.569   1.409   0.000  1.00 20.00           C
HETATM    2  C17 LIG A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C20 LIG A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    4 H171 LIG A   1      -0.313  -0.464   0.792  1.00 20.00           H
HETATM    5 H172 LIG A   1      -0.313  -0.464  -0.792  1.00 20.00           H
END
"""

pdb_str_004 = """
CRYST1  201.615  201.615  123.572  90.00  90.00 120.00 H 3
HETATM    1  N   6EL A 501     -18.527  39.011   4.685  1.00 42.85           N
HETATM    2  CA  6EL A 501     -19.238  37.942   3.975  1.00 41.65           C
HETATM    3  C   6EL A 501     -19.091  38.062   2.480  1.00 39.05           C
HETATM    4  O   6EL A 501     -18.732  39.122   1.981  1.00 35.88           O
HETATM    5  CB  6EL A 501     -18.971  36.528   4.429  1.00 41.97           C
HETATM    6  CG  6EL A 501     -20.045  35.653   3.781  1.00 34.50           C
HETATM    7  CD  6EL A 501     -20.198  35.965   2.286  1.00 34.24           C
HETATM    8  OE1 6EL A 501     -20.787  35.166   1.567  1.00 32.94           O
HETATM    9  NE2 6EL A 501     -19.708  37.135   1.735  1.00 37.42           N
HETATM   10  CAE 6EL A 501     -16.110  42.317   6.704  1.00 46.24           C
HETATM   11  CAF 6EL A 501     -17.290  43.026   6.465  1.00 43.30           C
HETATM   12  CAG 6EL A 501     -15.983  40.985   6.304  1.00 44.16           C
HETATM   13  CAH 6EL A 501     -18.362  42.399   5.839  1.00 46.24           C
HETATM   14  CAN 6EL A 501     -17.263  39.108   5.113  1.00 43.31           C
HETATM   15  CAO 6EL A 501     -19.045  40.227   4.740  1.00 46.27           C
HETATM   16  CAP 6EL A 501     -17.059  40.369   5.678  1.00 43.68           C
HETATM   17  CAQ 6EL A 501     -18.247  41.059   5.505  1.00 41.14           C
HETATM   18  OAC 6EL A 501     -16.455  38.154   5.201  1.00 37.79           O
HETATM   19  OAD 6EL A 501     -20.260  40.446   4.509  1.00 36.94           O
END
"""

pdb_str_005 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   0TD A 501       9.131  12.228   9.576  1.00 20.00           N
HETATM    2  CA  0TD A 501       9.235  10.781   9.451  1.00 20.00           C
HETATM    3  C   0TD A 501       7.932  10.147  10.006  1.00 20.00           C
HETATM    4  O   0TD A 501       7.556  10.575  11.117  1.00 20.00           O
HETATM    5  CSB 0TD A 501      13.190  10.241  10.653  1.00 20.00           C
HETATM    6  SB  0TD A 501      11.968  10.858   9.476  1.00 20.00           S
HETATM    7  CB  0TD A 501      10.441  10.157  10.187  1.00 20.00           C
HETATM    8  CG  0TD A 501      10.431   8.611  10.102  1.00 20.00           C
HETATM    9  OD2 0TD A 501      10.667   8.104   8.990  1.00 20.00           O
HETATM   10  OD1 0TD A 501      10.165   8.013  11.163  1.00 20.00           O
HETATM   11  OXT 0TD A 501       7.383   9.269   9.316  1.00 20.00           O
END
"""

pdb_str_006 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   GLZ A 501       9.839   9.850  11.172  1.00 20.00           N
HETATM    2  CA  GLZ A 501       8.601   9.331  10.629  1.00 20.00           C
HETATM    3  C   GLZ A 501       8.696   7.841  10.440  1.00 20.00           C
HETATM    4  O   GLZ A 501       8.048   7.204   9.648  1.00 20.00           O
END
"""

pdb_str_007 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  O1  PEO A 501      14.917  15.127  15.671  1.00 20.00           O
HETATM    2  O2  PEO A 501      15.083  15.207  14.264  1.00 20.00           O
END
"""

pdb_str_008 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C   MOH A 501      12.578  10.567  10.234  1.00 20.00           C
HETATM    2  O   MOH A 501      13.940  10.360   9.945  1.00 20.00           O
END
"""

pdb_str_008_control = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EOH A 501      10.157  14.225  12.482  1.00 20.00           C
HETATM    2  C2  EOH A 501      11.073  15.376  12.125  1.00 20.00           C
HETATM    3  O   EOH A 501      10.877  13.288  13.257  1.00 20.00           O
END
"""

pdb_str_009 = """
CRYST1   88.679  100.287  219.749  90.00  90.00  90.00 P 2 21 21
ATOM      1  N   ILE D   3      -4.315  36.803  -4.113  1.00 99.90           N
ATOM      2  CA  ILE D   3      -2.949  37.203  -4.414  1.00 90.56           C
ATOM      3  C   ILE D   3      -2.949  38.463  -5.269  1.00 91.77           C
ATOM      4  O   ILE D   3      -3.680  38.556  -6.255  1.00 91.19           O
ATOM      5  CB  ILE D   3      -2.182  36.067  -5.114  1.00 86.55           C
ATOM      6  CG1 ILE D   3      -2.220  34.795  -4.265  1.00 87.21           C
ATOM      7  CG2 ILE D   3      -0.747  36.482  -5.399  1.00 77.11           C
ATOM      8  CD1 ILE D   3      -1.760  33.558  -5.003  1.00 82.40           C
HETATM    9  N   9FZ D   4      -2.130  39.434  -4.879  1.00 91.16           N
HETATM   10  CA  9FZ D   4      -1.988  40.695  -5.655  1.00 89.10           C
HETATM   11  C   9FZ D   4      -1.176  40.466  -6.907  1.00 90.90           C
HETATM   12  O   9FZ D   4       0.010  40.188  -6.831  1.00 91.37           O
HETATM   13  CB  9FZ D   4      -1.337  41.773  -4.798  1.00 89.28           C
HETATM   14  C1  9FZ D   4      -1.144  43.047  -5.607  1.00 90.46           C
HETATM   15  C2  9FZ D   4      -0.679  44.202  -4.730  1.00 97.10           C
HETATM   16  C7  9FZ D   4      -1.020  47.117  -6.816  1.00106.60           C
HETATM   17  C8  9FZ D   4      -1.634  45.994  -6.287  1.00107.15           C
HETATM   18  C9  9FZ D   4      -1.906  47.950  -7.705  1.00102.31           C
HETATM   19  N4  9FZ D   4      -0.689  45.413  -5.545  1.00104.45           N
HETATM   20  N5  9FZ D   4       0.515  46.128  -5.583  1.00105.46           N
HETATM   21  N6  9FZ D   4       0.267  47.232  -6.418  1.00106.68           N
ATOM     22  N   ILE D   5      -1.820  40.576  -8.065  1.00 90.55           N
ATOM     23  CA  ILE D   5      -1.159  40.296  -9.335  1.00 92.65           C
ATOM     24  C   ILE D   5      -1.295  41.456 -10.315  1.00 90.25           C
ATOM     25  O   ILE D   5      -0.639  41.471 -11.356  1.00 88.80           O
ATOM     26  CB  ILE D   5      -1.710  39.006  -9.973  1.00 84.60           C
ATOM     27  CG1 ILE D   5      -3.209  39.150 -10.250  1.00 82.15           C
ATOM     28  CG2 ILE D   5      -1.432  37.800  -9.086  1.00 79.75           C
ATOM     29  CD1 ILE D   5      -3.813  37.968 -10.971  1.00 78.88           C
ATOM     30  N   LEU D   9       0.789  42.265 -13.727  1.00 85.51           N
ATOM     31  CA  LEU D   9       0.169  42.262 -15.045  1.00 87.53           C
ATOM     32  C   LEU D   9       0.491  43.540 -15.814  1.00 93.45           C
ATOM     33  O   LEU D   9       0.677  43.508 -17.030  1.00 95.25           O
ATOM     34  CB  LEU D   9      -1.347  42.092 -14.925  1.00 86.95           C
ATOM     35  CG  LEU D   9      -1.825  40.801 -14.256  1.00 85.95           C
ATOM     36  CD1 LEU D   9      -3.346  40.738 -14.220  1.00 82.56           C
ATOM     37  CD2 LEU D   9      -1.241  39.579 -14.952  1.00 67.49           C
HETATM   38  N   9G2 D  10       0.559  44.662 -15.103  1.00 96.10           N
HETATM   39  CA  9G2 D  10       0.902  45.962 -15.743  1.00 97.17           C
HETATM   40  C   9G2 D  10       2.361  46.013 -16.127  1.00 95.05           C
HETATM   41  O   9G2 D  10       2.735  46.711 -17.055  1.00 95.78           O
HETATM   42  CB  9G2 D  10       0.575  47.130 -14.821  1.00 95.86           C
HETATM   43  C1  9G2 D  10      -0.906  47.477 -14.891  1.00 98.06           C
HETATM   44  C10 9G2 D  10      -1.116  49.129  -8.268  1.00106.65           C
HETATM   45  C2  9G2 D  10      -1.181  48.865 -14.323  1.00107.10           C
HETATM   46  C7  9G2 D  10      -0.869  48.747 -10.699  1.00104.41           C
HETATM   47  C8  9G2 D  10      -0.231  48.723 -11.936  1.00101.61           C
HETATM   48  C9  9G2 D  10      -0.181  48.654  -9.359  1.00104.27           C
HETATM   49  N4  9G2 D  10      -1.186  48.828 -12.861  1.00107.44           N
HETATM   50  N5  9G2 D  10      -2.449  48.916 -12.275  1.00105.64           N
HETATM   51  N6  9G2 D  10      -2.206  48.860 -10.894  1.00106.43           N
ATOM     52  N   GLU D  11       3.190  45.273 -15.398  1.00 91.80           N
ATOM     53  CA  GLU D  11       4.607  45.175 -15.723  1.00 97.96           C
ATOM     54  C   GLU D  11       4.775  44.382 -17.015  1.00 99.71           C
ATOM     55  O   GLU D  11       5.723  44.593 -17.772  1.00103.15           O
ATOM     56  CB  GLU D  11       5.382  44.518 -14.579  1.00104.15           C
ATOM     57  CG  GLU D  11       6.888  44.481 -14.781  1.00112.57           C
ATOM     58  CD  GLU D  11       7.602  43.692 -13.700  1.00128.77           C
ATOM     59  OE1 GLU D  11       6.919  43.144 -12.809  1.00117.73           O
ATOM     60  OE2 GLU D  11       8.848  43.618 -13.744  1.00136.20           O
END
"""

pdb_str_013 = pdb_str_009

pdb_str_010 = """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  41      14.715 -15.370  14.923  1.00 39.91           N
ATOM      2  CA  GLU A  41      14.788 -16.780  14.439  1.00 35.65           C
ATOM      3  C   GLU A  41      13.961 -16.899  13.156  1.00 36.25           C
ATOM      4  O   GLU A  41      14.465 -17.527  12.205  1.00 38.69           O
ATOM      5  CB  GLU A  41      14.321 -17.749  15.518  1.00 33.50           C
ATOM      6  CG  GLU A  41      15.333 -17.993  16.600  1.00 35.41           C
ATOM      7  CD  GLU A  41      16.636 -18.604  16.105  1.00 43.11           C
ATOM      8  OE1 GLU A  41      16.630 -19.181  15.035  1.00 45.69           O
ATOM      9  OE2 GLU A  41      17.647 -18.494  16.799  1.00 51.56           O
ATOM     37  N   LEU A 107       4.234 -18.401  11.202  1.00 29.52           N
ATOM     38  CA  LEU A 107       5.506 -19.021  10.766  1.00 31.21           C
ATOM     39  C   LEU A 107       5.985 -20.013  11.824  1.00 29.33           C
ATOM     40  O   LEU A 107       5.845 -19.764  13.033  1.00 32.67           O
ATOM     41  CB  LEU A 107       6.554 -17.940  10.507  1.00 31.28           C
ATOM     42  CG  LEU A 107       6.196 -16.889   9.461  1.00 33.41           C
ATOM     43  CD1 LEU A 107       7.334 -15.910   9.340  1.00 35.68           C
ATOM     44  CD2 LEU A 107       5.884 -17.500   8.088  1.00 32.97           C
END
"""

pdb_str_011 = """
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM     10  N   ALA A  42      12.752 -16.313  13.130  1.00 37.05           N
ATOM     11  CA  ALA A  42      11.861 -16.326  11.950  1.00 41.68           C
ATOM     12  C   ALA A  42      12.606 -15.740  10.747  1.00 42.80           C
ATOM     13  O   ALA A  42      12.597 -16.391   9.684  1.00 36.14           O
ATOM     14  CB  ALA A  42      10.567 -15.603  12.226  1.00 43.61           C
ATOM     61  N   LYS A  87      14.521  55.475  65.158  1.00 16.69           N
ATOM     62  CA  LYS A  87      13.600  54.884  64.201  1.00 17.21           C
ATOM     63  C   LYS A  87      12.210  55.516  64.232  1.00 15.78           C
ATOM     64  O   LYS A  87      11.790  56.116  65.226  1.00 14.71           O
ATOM     65  CB  LYS A  87      13.500  53.368  64.403  1.00 18.88           C
ATOM     66  CG  LYS A  87      13.338  52.939  65.839  1.00 24.38           C
ATOM     67  CD  LYS A  87      11.961  52.424  66.104  1.00 30.22           C
ATOM     68  CE  LYS A  87      11.700  51.121  65.394  1.00 26.16           C
ATOM     69  NZ  LYS A  87      10.264  50.804  65.531  1.00 31.31           N
END
"""

pdb_str_012 = """
CRYST1   52.650   54.499   66.365  93.37  97.68 102.44 P 1
HETATM    1  C1' AMP A1201       1.751   9.985  45.277  1.00 24.11           C
HETATM    2  C2  AMP A1201       3.845   8.243  48.752  1.00 26.36           C
HETATM    3  C2' AMP A1201       2.736  11.117  44.952  1.00 24.56           C
HETATM    4  C3' AMP A1201       2.584  11.238  43.440  1.00 24.53           C
HETATM    5  C4  AMP A1201       3.059   8.342  46.558  1.00 23.44           C
HETATM    6  C4' AMP A1201       1.117  10.895  43.217  1.00 26.37           C
HETATM    7  C5  AMP A1201       3.693   7.052  46.280  1.00 24.64           C
HETATM    8  C5' AMP A1201       0.762  10.378  41.845  1.00 27.37           C
HETATM    9  C6  AMP A1201       4.423   6.413  47.392  1.00 24.61           C
HETATM   10  C8  AMP A1201       2.694   7.786  44.489  1.00 23.99           C
HETATM   11  N1  AMP A1201       4.457   7.065  48.572  1.00 22.89           N
HETATM   12  N3  AMP A1201       3.166   8.890  47.789  1.00 25.59           N
HETATM   13  N6  AMP A1201       5.046   5.216  47.246  1.00 26.07           N
HETATM   14  N7  AMP A1201       3.442   6.772  44.990  1.00 23.32           N
HETATM   15  N9  AMP A1201       2.467   8.713  45.433  1.00 24.84           N
HETATM   16  O1P AMP A1201      -0.476  10.260  39.151  1.00 31.62           O
HETATM   17  O2' AMP A1201       2.299  12.295  45.610  1.00 27.03           O
HETATM   18  O2P AMP A1201      -1.283   8.318  40.533  1.00 33.06           O
HETATM   19  O3' AMP A1201       2.930  12.509  42.922  1.00 26.18           O
HETATM   20  O3P AMP A1201      -2.715  10.330  40.359  1.00 28.38           O
HETATM   21  O4' AMP A1201       0.851   9.868  44.194  1.00 25.17           O
HETATM   22  O5' AMP A1201      -0.643  10.350  41.667  1.00 29.44           O
HETATM   23  P   AMP A1201      -1.307   9.798  40.318  1.00 33.89           P
END
"""

pdb_str_012_control = """
CRYST1   17.880   31.420   43.900  90.00  90.00  90.00 P 21 21 21    8
ATOM      1  O5'  DC A   1      19.545  18.136  17.917  1.00  3.07           O
ATOM      2  C5'  DC A   1      19.769  17.119  18.884  1.00  2.46           C
ATOM      3  C4'  DC A   1      18.610  16.148  19.001  1.00  2.19           C
ATOM      4  O4'  DC A   1      17.462  16.852  19.514  1.00  2.62           O
ATOM      5  C3'  DC A   1      18.161  15.506  17.674  1.00  2.26           C
ATOM      6  O3'  DC A   1      17.782  14.139  17.875  1.00  2.47           O
ATOM      7  C2'  DC A   1      16.906  16.282  17.315  1.00  2.41           C
ATOM      8  C1'  DC A   1      16.340  16.624  18.692  1.00  2.34           C
ATOM      9  N1   DC A   1      15.516  17.837  18.704  1.00  2.30           N
ATOM     10  C2   DC A   1      14.145  17.720  18.492  1.00  2.34           C
ATOM     11  O2   DC A   1      13.658  16.581  18.329  1.00  3.07           O
ATOM     12  N3   DC A   1      13.385  18.831  18.454  1.00  2.35           N
ATOM     13  C4   DC A   1      13.943  20.039  18.611  1.00  2.37           C
ATOM     14  N4   DC A   1      13.161  21.108  18.580  1.00  2.79           N
ATOM     15  C5   DC A   1      15.357  20.189  18.812  1.00  2.71           C
ATOM     16  C6   DC A   1      16.103  19.061  18.841  1.00  2.68           C
ATOM     17  P    DG A   2      18.825  12.942  17.684  1.00  2.51           P
ATOM     18  OP1  DG A   2      19.788  13.206  16.573  1.00  3.24           O
ATOM     19  OP2  DG A   2      17.976  11.710  17.621  1.00  3.25           O
ATOM     20  O5'  DG A   2      19.719  12.937  19.002  1.00  2.66           O
ATOM     21  C5'  DG A   2      19.103  12.733  20.284  1.00  2.85           C
ATOM     22  C4'  DG A   2      20.140  13.045  21.335  1.00  2.73           C
ATOM     23  O4'  DG A   2      20.546  14.399  21.207  1.00  2.76           O
ATOM     24  C3'  DG A   2      19.598  12.910  22.753  1.00  2.85           C
ATOM     25  O3'  DG A   2      19.812  11.563  23.230  1.00  3.51           O
ATOM     26  C2'  DG A   2      20.430  13.919  23.526  1.00  3.30           C
ATOM     27  C1'  DG A   2      20.834  14.964  22.481  1.00  2.85           C
ATOM     28  N9   DG A   2      20.140  16.222  22.572  1.00  2.76           N
ATOM     29  C8   DG A   2      20.744  17.451  22.654  1.00  3.56           C
ATOM     30  N7   DG A   2      19.903  18.441  22.626  1.00  3.75           N
ATOM     31  C5   DG A   2      18.658  17.840  22.510  1.00  2.68           C
ATOM     32  C6   DG A   2      17.359  18.399  22.397  1.00  2.62           C
ATOM     33  O6   DG A   2      17.067  19.612  22.382  1.00  3.52           O
ATOM     34  N1   DG A   2      16.371  17.435  22.313  1.00  2.19           N
ATOM     35  C2   DG A   2      16.610  16.084  22.263  1.00  1.99           C
ATOM     36  N2   DG A   2      15.540  15.294  22.108  1.00  2.48           N
ATOM     37  N3   DG A   2      17.814  15.534  22.356  1.00  2.10           N
ATOM     38  C4   DG A   2      18.782  16.466  22.477  1.00  2.26           C
ATOM     39  P    DC A   3      18.598  10.645  23.664  1.00  3.54           P
ATOM     40  OP1  DC A   3      19.187   9.314  23.999  1.00  5.42           O
ATOM     41  OP2  DC A   3      17.526  10.682  22.656  1.00  5.19           O
ATOM     42  O5'  DC A   3      18.133  11.412  24.998  1.00  3.68           O
ATOM     43  C5'  DC A   3      17.562  10.688  26.104  1.00  3.12           C
ATOM     44  C4'  DC A   3      16.175  11.192  26.422  1.00  2.67           C
ATOM     45  O4'  DC A   3      16.263  12.557  26.885  1.00  2.93           O
ATOM     46  C3'  DC A   3      15.179  11.206  25.269  1.00  2.83           C
ATOM     47  O3'  DC A   3      13.854  10.971  25.771  1.00  3.13           O
ATOM     48  C2'  DC A   3      15.208  12.650  24.803  1.00  3.07           C
ATOM     49  C1'  DC A   3      15.379  13.367  26.137  1.00  2.64           C
ATOM     50  N1   DC A   3      15.979  14.696  26.016  1.00  2.59           N
ATOM     51  C2   DC A   3      15.126  15.784  25.823  1.00  2.60           C
ATOM     52  O2   DC A   3      13.894  15.575  25.734  1.00  3.08           O
ATOM     53  N3   DC A   3      15.658  17.016  25.725  1.00  2.61           N
ATOM     54  C4   DC A   3      16.987  17.184  25.770  1.00  2.65           C
ATOM     55  N4   DC A   3      17.464  18.422  25.674  1.00  3.18           N
ATOM     56  C5   DC A   3      17.881  16.072  25.908  1.00  3.00           C
ATOM     57  C6   DC A   3      17.330  14.855  26.026  1.00  2.98           C
END
"""

pdb_str_012_linked = """
CRYST1   56.300   33.400   63.000  90.00  90.25  90.00 P 1 21 1      2
ATOM      1  OP3 5GP A   1      23.215   5.145  51.161  1.00  0.00           O
ATOM      2  P   5GP A   1      24.650   4.594  51.620  1.00  0.00           P
ATOM      3  OP1 5GP A   1      24.810   3.219  51.082  1.00  0.00           O
ATOM      4  OP2 5GP A   1      25.707   5.598  51.344  1.00  0.00           O
ATOM      5  O5' 5GP A   1      24.476   4.546  53.204  1.00  0.00           O
ATOM      6  C5' 5GP A   1      25.156   5.525  53.978  1.00  0.00           C
ATOM      7  C4' 5GP A   1      25.808   4.801  55.152  1.00  0.00           C
ATOM      8  O4' 5GP A   1      24.983   3.707  55.565  1.00  0.00           O
ATOM      9  C3' 5GP A   1      27.185   4.241  54.785  1.00  0.00           C
ATOM     10  O3' 5GP A   1      28.130   5.240  55.163  1.00  0.00           O
ATOM     11  C2' 5GP A   1      27.282   3.124  55.808  1.00  0.00           C
ATOM     12  O2' 5GP A   1      27.691   3.684  57.045  1.00  0.00           O
ATOM     13  C1' 5GP A   1      25.833   2.612  55.903  1.00  0.00           C
ATOM     14  N9  5GP A   1      25.681   1.623  54.832  1.00  0.00           N
ATOM     15  C8  5GP A   1      24.865   1.712  53.755  1.00  0.00           C
ATOM     16  N7  5GP A   1      25.012   0.580  52.989  1.00  0.00           N
ATOM     17  C5  5GP A   1      25.874  -0.248  53.697  1.00  0.00           C
ATOM     18  C6  5GP A   1      26.336  -1.541  53.483  1.00  0.00           C
ATOM     19  O6  5GP A   1      25.983  -2.159  52.466  1.00  0.00           O
ATOM     20  N1  5GP A   1      27.127  -2.122  54.405  1.00  0.00           N
ATOM     21  C2  5GP A   1      27.481  -1.468  55.532  1.00  0.00           C
ATOM     22  N2  5GP A   1      28.055  -2.114  56.534  1.00  0.00           N
ATOM     23  N3  5GP A   1      27.078  -0.243  55.747  1.00  0.00           N
ATOM     24  C4  5GP A   1      26.260   0.388  54.839  1.00  0.00           C
ATOM     25  P     C A   2      29.640   5.301  54.591  1.00  0.00           P
ATOM     26  OP1   C A   2      30.307   6.511  55.152  1.00  0.00           O
ATOM     27  OP2   C A   2      29.605   5.142  53.117  1.00  0.00           O
ATOM     28  O5'   C A   2      30.332   4.036  55.292  1.00  0.00           O
ATOM     29  C5'   C A   2      30.763   4.137  56.639  1.00  0.00           C
ATOM     30  C4'   C A   2      31.424   2.812  56.982  1.00  0.00           C
ATOM     31  O4'   C A   2      30.483   1.753  56.828  1.00  0.00           O
ATOM     32  C3'   C A   2      32.562   2.541  56.002  1.00  0.00           C
ATOM     33  O3'   C A   2      33.772   3.124  56.479  1.00  0.00           O
ATOM     34  C2'   C A   2      32.637   1.038  56.077  1.00  0.00           C
ATOM     35  O2'   C A   2      33.390   0.643  57.211  1.00  0.00           O
ATOM     36  C1'   C A   2      31.162   0.636  56.248  1.00  0.00           C
ATOM     37  N1    C A   2      30.570   0.343  54.911  1.00  0.00           N
ATOM     38  C2    C A   2      30.799  -0.845  54.338  1.00  0.00           C
ATOM     39  O2    C A   2      31.777  -1.511  54.728  1.00  0.00           O
ATOM     40  N3    C A   2      30.243  -1.143  53.141  1.00  0.00           N
ATOM     41  C4    C A   2      29.496  -0.265  52.481  1.00  0.00           C
ATOM     42  N4    C A   2      28.782  -0.652  51.425  1.00  0.00           N
ATOM     43  C5    C A   2      29.284   0.988  53.032  1.00  0.00           C
ATOM     44  C6    C A   2      29.844   1.264  54.273  1.00  0.00           C
END
"""

pdb_str_014 = """
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N   MEA A   1      73.283  60.707   2.566  1.00 47.48           N
HETATM    2  CA  MEA A   1      74.453  61.240   1.817  1.00 47.82           C
HETATM    3  C   MEA A   1      74.572  60.601   0.426  1.00 47.58           C
HETATM    4  O   MEA A   1      74.374  61.266  -0.601  1.00 48.32           O
HETATM    5  CB  MEA A   1      75.759  60.973   2.578  1.00 50.44           C
HETATM    6  CG  MEA A   1      75.799  61.532   3.983  1.00 53.06           C
HETATM    7  CD1 MEA A   1      76.807  61.121   4.869  1.00 53.98           C
HETATM    8  CD2 MEA A   1      74.845  62.444   4.432  1.00 52.76           C
HETATM    9  CE1 MEA A   1      76.863  61.605   6.175  1.00 53.65           C
HETATM   10  CE2 MEA A   1      74.898  62.932   5.742  1.00 53.06           C
HETATM   11  CZ  MEA A   1      75.906  62.511   6.612  1.00 54.01           C
ATOM     19  N   LEU A   3      73.863  57.566  -2.661  1.00 43.69           N
ATOM     20  CA  LEU A   3      72.718  56.937  -3.309  1.00 43.82           C
ATOM     21  C   LEU A   3      72.246  55.652  -2.643  1.00 43.49           C
ATOM     22  O   LEU A   3      71.088  55.537  -2.256  1.00 44.24           O
ATOM     23  CB  LEU A   3      73.043  56.613  -4.762  1.00 45.49           C
ATOM     24  CG  LEU A   3      72.022  55.671  -5.400  1.00 45.57           C
ATOM     25  CD1 LEU A   3      70.684  56.381  -5.444  1.00 45.67           C
ATOM     26  CD2 LEU A   3      72.460  55.263  -6.795  1.00 46.15           C
END
"""

pdb_str_014_ht7 = """
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N   HT7 H   3     -41.786   7.174  23.475  1.00 20.96           N
HETATM    2  CA  HT7 H   3     -41.367   8.455  25.530  1.00 24.74           C
HETATM    3  C   HT7 H   3     -39.954   8.864  25.125  1.00 28.36           C
HETATM    4  O   HT7 H   3     -39.671   9.203  23.968  1.00 25.39           O
HETATM    5  CB  HT7 H   3     -42.289   8.244  24.314  1.00 23.13           C
HETATM    6  CG  HT7 H   3     -43.736   7.905  24.673  1.00 24.83           C
HETATM    7  CD  HT7 H   3     -44.619   7.834  23.465  1.00 26.29           C
HETATM    8  CE1 HT7 H   3     -45.324   8.860  22.907  1.00 26.79           C
HETATM    9  CE2 HT7 H   3     -44.830   6.695  22.610  1.00 28.69           C
HETATM   10  NZ1 HT7 H   3     -45.987   8.430  21.785  1.00 28.51           N
HETATM   11  CZ2 HT7 H   3     -45.699   7.114  21.574  1.00 29.24           C
HETATM   12  CZ3 HT7 H   3     -44.383   5.366  22.628  1.00 29.19           C
HETATM   13  CH2 HT7 H   3     -46.149   6.246  20.583  1.00 32.02           C
HETATM   14  CH3 HT7 H   3     -44.805   4.510  21.637  1.00 29.79           C
HETATM   15  CT2 HT7 H   3     -45.670   4.949  20.628  1.00 32.39           C
END
"""

pdb_str_015 = """
HETATM  421  P   BGM A  14     -11.641  -6.888   0.692  1.00  0.00           P
HETATM  422  OP1 BGM A  14     -12.726  -7.004  -0.316  1.00  0.00           O
HETATM  423  OP2 BGM A  14     -11.949  -7.102   2.128  1.00  0.00           O
HETATM  424  O5' BGM A  14     -10.470  -7.884   0.276  1.00  0.00           O
HETATM  425  C5' BGM A  14      -9.316  -8.044   1.100  1.00  0.00           C
HETATM  426  C4' BGM A  14      -8.293  -8.912   0.405  1.00  0.00           C
HETATM  427  O4' BGM A  14      -7.610  -8.119  -0.598  1.00  0.00           O
HETATM  428  C1' BGM A  14      -6.214  -8.121  -0.343  1.00  0.00           C
HETATM  429  N9  BGM A  14      -5.686  -6.811  -0.712  1.00  0.00           N
HETATM  430  C8  BGM A  14      -5.048  -6.468  -1.881  1.00  0.00           C
HETATM  431  N7  BGM A  14      -4.689  -5.214  -1.921  1.00  0.00           N
HETATM  432  C5  BGM A  14      -5.114  -4.697  -0.705  1.00  0.00           C
HETATM  433  C4  BGM A  14      -5.731  -5.670   0.053  1.00  0.00           C
HETATM  434  N3  BGM A  14      -6.274  -5.562   1.283  1.00  0.00           N
HETATM  435  C2  BGM A  14      -6.165  -4.336   1.766  1.00  0.00           C
HETATM  436  N2  BGM A  14      -6.653  -4.052   2.982  1.00  0.00           N
HETATM  437  N1  BGM A  14      -5.570  -3.296   1.093  1.00  0.00           N
HETATM  438  C6  BGM A  14      -5.003  -3.385  -0.175  1.00  0.00           C
HETATM  439  O6  BGM A  14      -4.490  -2.385  -0.693  1.00  0.00           O
HETATM  440  C2' BGM A  14      -6.064  -8.472   1.123  1.00  0.00           C
HETATM  441  C3' BGM A  14      -7.194  -9.468   1.312  1.00  0.00           C
HETATM  442  O3' BGM A  14      -6.817 -10.779   0.874  1.00  0.00           O
HETATM  443 BR   BGM A  14      -4.719  -7.696  -3.314  1.00  0.00          BR
HETATM  454  P   BGM A  15      -5.374 -11.371   1.266  1.00  0.00           P
HETATM  455  OP1 BGM A  15      -5.290 -12.745   0.711  1.00  0.00           O
HETATM  456  OP2 BGM A  15      -5.159 -11.152   2.719  1.00  0.00           O
HETATM  457  O5' BGM A  15      -4.354 -10.448   0.463  1.00  0.00           O
HETATM  458  C5' BGM A  15      -2.945 -10.667   0.541  1.00  0.00           C
HETATM  459  C4' BGM A  15      -2.335 -10.632  -0.841  1.00  0.00           C
HETATM  460  O4' BGM A  15      -2.416  -9.278  -1.353  1.00  0.00           O
HETATM  461  C1' BGM A  15      -1.159  -8.905  -1.881  1.00  0.00           C
HETATM  462  N9  BGM A  15      -1.032  -7.455  -1.771  1.00  0.00           N
HETATM  463  C8  BGM A  15      -0.465  -6.598  -2.684  1.00  0.00           C
HETATM  464  N7  BGM A  15      -0.488  -5.351  -2.301  1.00  0.00           N
HETATM  465  C5  BGM A  15      -1.112  -5.385  -1.061  1.00  0.00           C
HETATM  466  C4  BGM A  15      -1.455  -6.676  -0.721  1.00  0.00           C
HETATM  467  N3  BGM A  15      -2.073  -7.118   0.394  1.00  0.00           N
HETATM  468  C2  BGM A  15      -2.357  -6.130   1.226  1.00  0.00           C
HETATM  469  N2  BGM A  15      -2.977  -6.394   2.385  1.00  0.00           N
HETATM  470  N1  BGM A  15      -2.057  -4.811   0.983  1.00  0.00           N
HETATM  471  C6  BGM A  15      -1.420  -4.334  -0.159  1.00  0.00           C
HETATM  472  O6  BGM A  15      -1.194  -3.123  -0.279  1.00  0.00           O
HETATM  473  C2' BGM A  15      -0.142  -9.693  -1.085  1.00  0.00           C
HETATM  474  C3' BGM A  15      -0.854 -11.020  -0.902  1.00  0.00           C
HETATM  475  O3' BGM A  15      -0.612 -11.878  -2.021  1.00  0.00           O
HETATM  476 BR   BGM A  15       0.273  -7.165  -4.357  1.00  0.00          BR
END
"""

pdb_str_016 = """
HETATM    1  P   35G A   1      40.589  33.310  86.588  1.00 20.00           P
HETATM    2  O1P 35G A   1      41.439  33.145  85.410  1.00 20.00           O
HETATM    3  O2P 35G A   1      39.279  34.185  86.398  1.00 20.00           O
HETATM    4  O5' 35G A   1      41.295  34.060  87.788  1.00 20.00           O
HETATM    5  C5' 35G A   1      40.748  34.048  89.142  1.00 20.00           C
HETATM    6  C4' 35G A   1      40.485  32.596  89.477  1.00 20.00           C
HETATM    7  O4' 35G A   1      39.809  32.402  90.710  1.00 20.00           O
HETATM    8  C3' 35G A   1      39.519  31.984  88.475  1.00 20.00           C
HETATM    9  O3' 35G A   1      40.103  31.901  87.171  1.00 20.00           O
HETATM   10  C2' 35G A   1      39.158  30.672  89.148  1.00 20.00           C
HETATM   11  O2' 35G A   1      40.202  29.745  88.940  1.00 20.00           O
HETATM   12  C1' 35G A   1      39.105  31.145  90.637  1.00 20.00           C
HETATM   13  N9  35G A   1      37.779  31.316  91.172  1.00 20.00           N
HETATM   14  C8  35G A   1      36.714  32.012  90.621  1.00 20.00           C
HETATM   15  N7  35G A   1      35.655  31.993  91.380  1.00 20.00           N
HETATM   16  C5  35G A   1      36.022  31.256  92.493  1.00 20.00           C
HETATM   17  C6  35G A   1      35.282  30.898  93.657  1.00 20.00           C
HETATM   18  O6  35G A   1      34.119  31.163  93.950  1.00 20.00           O
HETATM   19  N1  35G A   1      36.074  30.137  94.539  1.00 20.00           N
HETATM   20  C2  35G A   1      37.382  29.777  94.342  1.00 20.00           C
HETATM   21  N2  35G A   1      37.987  29.092  95.330  1.00 20.00           N
HETATM   22  N3  35G A   1      38.055  30.104  93.256  1.00 20.00           N
HETATM   23  C4  35G A   1      37.340  30.832  92.383  1.00 20.00           C
END
"""

pdb_str_018 = """
MODEL        1
ATOM      1  N   GLY A 202      17.011  19.566   6.604  1.00  0.00           N
ATOM      2  CA  GLY A 202      15.941  18.568   6.783  1.00  0.00           C
ATOM      3  C   GLY A 202      15.797  18.196   8.252  1.00  0.00           C
ATOM      4  O   GLY A 202      16.055  19.048   9.103  1.00  0.00           O
ATOM     10  N   PRO A 203      15.393  16.952   8.573  1.00  0.00           N
ATOM     11  CA  PRO A 203      15.401  16.399   9.933  1.00  0.00           C
ATOM     12  C   PRO A 203      16.792  15.872  10.343  1.00  0.00           C
ATOM     13  O   PRO A 203      17.746  15.930   9.561  1.00  0.00           O
ATOM     14  CB  PRO A 203      14.365  15.269   9.880  1.00  0.00           C
ATOM     15  CG  PRO A 203      14.551  14.718   8.469  1.00  0.00           C
ATOM     16  CD  PRO A 203      14.860  15.967   7.639  1.00  0.00           C
ATOM     24  N   HIS A 204      16.902  15.311  11.556  1.00  0.00           N
ATOM     25  CA  HIS A 204      18.108  14.658  12.103  1.00  0.00           C
ATOM     26  C   HIS A 204      19.382  15.513  12.044  1.00  0.00           C
ATOM     27  O   HIS A 204      20.489  14.992  11.921  1.00  0.00           O
ATOM     28  CB  HIS A 204      18.319  13.285  11.452  1.00  0.00           C
ATOM     29  CG  HIS A 204      17.114  12.383  11.482  1.00  0.00           C
ATOM     30  ND1 HIS A 204      16.433  11.921  10.375  1.00  0.00           N1+
ATOM     31  CD2 HIS A 204      16.521  11.852  12.592  1.00  0.00           C
ATOM     32  CE1 HIS A 204      15.446  11.111  10.804  1.00  0.00           C
ATOM     33  NE2 HIS A 204      15.476  11.063  12.151  1.00  0.00           N
ENDMDL
MODEL        2
ATOM      1  N   GLY A 202      20.851  15.633   2.605  1.00  0.00           N
ATOM      2  CA  GLY A 202      21.227  14.437   3.379  1.00  0.00           C
ATOM      3  C   GLY A 202      20.033  13.890   4.148  1.00  0.00           C
ATOM      4  O   GLY A 202      18.905  14.325   3.908  1.00  0.00           O
ATOM     10  N   PRO A 203      20.237  12.925   5.066  1.00  0.00           N
ATOM     11  CA  PRO A 203      19.142  12.235   5.753  1.00  0.00           C
ATOM     12  C   PRO A 203      18.563  13.008   6.950  1.00  0.00           C
ATOM     13  O   PRO A 203      17.344  13.033   7.125  1.00  0.00           O
ATOM     14  CB  PRO A 203      19.736  10.885   6.174  1.00  0.00           C
ATOM     15  CG  PRO A 203      21.221  11.190   6.377  1.00  0.00           C
ATOM     16  CD  PRO A 203      21.504  12.252   5.313  1.00  0.00           C
ATOM     24  N   HIS A 204      19.419  13.622   7.780  1.00  0.00           N
ATOM     25  CA  HIS A 204      19.057  14.325   9.023  1.00  0.00           C
ATOM     26  C   HIS A 204      20.059  15.445   9.334  1.00  0.00           C
ATOM     27  O   HIS A 204      21.244  15.297   9.043  1.00  0.00           O
ATOM     28  CB  HIS A 204      19.006  13.330  10.204  1.00  0.00           C
ATOM     29  CG  HIS A 204      18.099  12.139  10.029  1.00  0.00           C
ATOM     30  ND1 HIS A 204      18.479  10.807   9.996  1.00  0.00           N1+
ATOM     31  CD2 HIS A 204      16.739  12.183   9.926  1.00  0.00           C
ATOM     32  CE1 HIS A 204      17.366  10.061   9.859  1.00  0.00           C
ATOM     33  NE2 HIS A 204      16.297  10.881   9.817  1.00  0.00           N
ENDMDL
END
"""

pdb_str_019 = """
ATOM   5392  N   ASN U 172      28.919 -35.712  34.830  1.00 43.06           N
ATOM   5393  CA  ASN U 172      29.803 -36.865  34.659  1.00 42.35           C
ATOM   5394  C   ASN U 172      30.623 -36.917  33.379  1.00 40.69           C
ATOM   5395  O   ASN U 172      30.581 -37.901  32.648  1.00 39.42           O
ATOM   5396  CB  ASN U 172      28.999 -38.145  34.780  1.00 43.40           C
ATOM   5397  CG  ASN U 172      27.809 -38.142  33.876  1.00 50.40           C
ATOM   5398  OD1 ASN U 172      27.350 -37.080  33.448  1.00 50.72           O
ATOM   5399  ND2 ASN U 172      27.289 -39.330  33.572  1.00 54.60           N
HETATM 6170  C1  NAG C   1      26.142 -38.639  32.099  1.00 77.20           C
HETATM 6171  C2  NAG C   1      25.260 -39.875  32.134  1.00 77.28           C
HETATM 6172  C3  NAG C   1      24.134 -39.690  31.125  1.00 76.77           C
HETATM 6173  C4  NAG C   1      23.347 -38.412  31.448  1.00 75.42           C
HETATM 6174  C5  NAG C   1      24.298 -37.207  31.547  1.00 76.20           C
HETATM 6175  C6  NAG C   1      23.588 -35.950  32.029  1.00 76.67           C
HETATM 6176  C7  NAG C   1      27.298 -41.161  32.257  1.00 77.02           C
HETATM 6177  C8  NAG C   1      28.395 -41.056  31.211  1.00 77.44           C
HETATM 6178  N2  NAG C   1      26.041 -41.060  31.830  1.00 77.79           N
HETATM 6179  O1  NAG C   1      27.173 -38.796  33.005  1.00 78.97           O
HETATM 6180  O3  NAG C   1      23.271 -40.814  31.167  1.00 80.30           O
HETATM 6181  O4  NAG C   1      22.383 -38.172  30.432  1.00 71.73           O
HETATM 6182  O5  NAG C   1      25.376 -37.484  32.475  1.00 76.45           O
HETATM 6183  O6  NAG C   1      24.190 -35.421  33.203  1.00 79.30           O
HETATM 6184  O7  NAG C   1      27.592 -41.308  33.445  1.00 75.78           O
END
"""

pdb_str_019_control = """
ATOM    859  N   ASN A 114      13.568   0.065  35.530  1.00  8.61           N
ATOM    860  CA  ASN A 114      13.918  -0.275  36.879  1.00  9.29           C
ATOM    861  C   ASN A 114      15.288   0.210  37.331  1.00 10.96           C
ATOM    862  O   ASN A 114      15.734  -0.079  38.422  1.00 15.40           O
ATOM    863  CB  ASN A 114      13.746  -1.754  37.229  1.00 11.35           C
ATOM    864  CG  ASN A 114      12.350  -2.251  37.154  1.00 11.55           C
ATOM    865  OD1 ASN A 114      11.372  -1.469  37.146  1.00 11.66           O
ATOM    866  ND2 ASN A 114      12.173  -3.567  37.053  1.00 12.97           N
HETATM 1710  C1  NAG A 401      10.902  -4.146  37.038  1.00 13.79           C
HETATM 1711  C2  NAG A 401      10.956  -5.537  36.326  1.00 15.71           C
HETATM 1712  C3  NAG A 401       9.504  -6.092  36.328  1.00 17.70           C
HETATM 1713  C4  NAG A 401       9.043  -6.173  37.788  1.00 20.69           C
HETATM 1714  C5  NAG A 401       9.072  -4.723  38.394  1.00 21.37           C
HETATM 1715  C6  NAG A 401       8.589  -4.620  39.782  1.00 28.32           C
HETATM 1716  C7  NAG A 401      12.396  -5.712  34.327  1.00 17.48           C
HETATM 1717  C8  NAG A 401      12.509  -5.462  32.899  1.00 22.32           C
HETATM 1718  N2  NAG A 401      11.252  -5.279  34.917  1.00 15.00           N
HETATM 1719  O3  NAG A 401       9.545  -7.411  35.765  1.00 21.47           O
HETATM 1720  O4  NAG A 401       7.678  -6.625  37.777  1.00 25.59           O
HETATM 1721  O5  NAG A 401      10.488  -4.329  38.368  1.00 17.19           O
HETATM 1722  O6  NAG A 401       9.364  -5.399  40.659  1.00 38.39           O
HETATM 1723  O7  NAG A 401      13.161  -6.362  34.883  1.00 19.28           O
END
"""

pdb_str_020 = """
CRYST1   35.150   61.360   88.940  90.00  90.00  90.00 C 2 2 21     16
ATOM   2094  N  ASER A 102      -0.505 -11.072   2.316  0.62 14.20           N
ATOM   2095  N  BSER A 102      -0.440 -11.462   2.296  0.38  8.25           N
ATOM   2096  CA ASER A 102       0.562 -10.397   1.608  0.62 14.24           C
ATOM   2097  CA BSER A 102       0.488 -10.636   1.545  0.38  9.67           C
ATOM   2098  C  ASER A 102       0.420  -8.947   2.004  0.62 13.70           C
ATOM   2099  C  BSER A 102       0.444  -9.171   1.999  0.38  8.99           C
ATOM   2100  O  ASER A 102       1.376  -8.180   2.169  0.62 14.83           O
ATOM   2101  O  BSER A 102       1.496  -8.531   2.119  0.38  9.78           O
ATOM   2102  CB ASER A 102       0.481 -10.668   0.054  0.62 12.00           C
ATOM   2103  CB BSER A 102       0.155 -10.652   0.056  0.38 14.81           C
ATOM   2104  OG ASER A 102       0.909 -11.911  -0.426  0.62 12.63           O
ATOM   2105  OG BSER A 102       1.222 -10.020  -0.610  0.38 17.17           O
ATOM   2116  N   ALA A 103      -0.777  -8.531   2.155  1.00 12.01           N
ATOM   2117  CA AALA A 103      -1.028  -7.141   2.663  0.62 12.55           C
ATOM   2118  CA BALA A 103      -0.799  -7.167   2.561  0.38 11.32           C
ATOM   2119  C   ALA A 103      -0.434  -6.917   4.050  1.00 11.54           C
ATOM   2120  O   ALA A 103       0.175  -5.884   4.351  1.00 12.98           O
ATOM   2121  CB AALA A 103      -2.475  -6.621   2.598  0.62 13.63           C
ATOM   2122  CB BALA A 103      -2.223  -6.688   2.343  0.38 12.27           C
END
"""

pdb_str_021 = """
CRYST1   57.482   57.482   54.469  90.00  90.00  90.00 P 42 21 2     8
ATOM     84  N   GLY A  10       3.511  23.530   4.999  1.00  6.37           N
ATOM     85  CA  GLY A  10       3.794  22.322   4.237  1.00  6.80           C
ATOM     86  C   GLY A  10       3.735  22.565   2.744  1.00  6.99           C
ATOM     87  O   GLY A  10       3.789  23.705   2.292  1.00  9.58           O
ATOM     88  N   ASP A  11       3.619  21.484   1.977  1.00  6.21           N
ATOM     89  CA  ASP A  11       3.612  21.594   0.527  1.00  6.44           C
ATOM     90  C   ASP A  11       2.281  22.091  -0.033  1.00  6.88           C
ATOM     91  O   ASP A  11       2.200  22.414  -1.208  1.00  7.79           O
ATOM     92  CB  ASP A  11       4.111  20.307  -0.163  1.00  6.77           C
ATOM     93  CGX ASP A  11       3.109  19.154  -0.141  1.00  7.29           C
ATOM     94  OD1 ASP A  11       2.010  19.254   0.435  1.00  6.86           O
ATOM     95  OD2 ASP A  11       3.456  18.124  -0.766  1.00 10.63           O
ATOM     96  N   GLY A  12       1.252  22.153   0.807  1.00  7.18           N
ATOM     97  CA  GLY A  12      -0.056  22.616   0.386  1.00  8.11           C
ATOM     98  C   GLY A  12      -0.787  21.679  -0.546  1.00  8.25           C
ATOM     99  O   GLY A  12      -1.761  22.086  -1.177  1.00  9.47           O
END
"""

pdb_str_022 = """
HETATM    1  C   ACB A   1       2.446   1.907  -2.265  1.00  0.00           C
HETATM    2  O   ACB A   1       1.420   2.124  -1.611  1.00  0.00           O
HETATM    3  OXT ACB A   1       2.511   0.914  -3.003  1.00  0.00           O
HETATM    4  CA  ACB A   1       3.666   2.852  -2.160  1.00  0.00           C
HETATM    5  N   ACB A   1       4.493   2.334  -1.074  1.00  0.00           N
HETATM    6  CB  ACB A   1       3.309   4.300  -1.999  1.00  0.00           C
HETATM    7  CG  ACB A   1       2.709   4.844  -0.712  1.00  0.00           C
HETATM    8  C4  ACB A   1       4.466   5.181  -2.274  1.00  0.00           C
HETATM    9  OD2 ACB A   1       3.383   4.704   0.319  1.00  0.00           O
ATOM     16  N   VAL A   2       1.508   5.383  -0.621  1.00  0.00           N
ATOM     17  CA  VAL A   2       0.831   5.820   0.601  1.00  0.00           C
ATOM     18  C   VAL A   2      -0.168   4.765   1.046  1.00  0.00           C
ATOM     19  O   VAL A   2      -1.267   5.122   1.495  1.00  0.00           O
ATOM     20  CB  VAL A   2       0.462   7.363   0.614  1.00  0.00           C
ATOM     21  CG1 VAL A   2       0.534   8.119  -0.686  1.00  0.00           C
ATOM     22  CG2 VAL A   2      -0.874   7.823   1.149  1.00  0.00           C
HETATM   32  C1  1ZN A   3      -7.398  -0.933  -3.062  1.00  0.00           C
HETATM   33  O1  1ZN A   3      -6.743  -0.971  -1.826  1.00  0.00           O
HETATM   34  C2  1ZN A   3      -5.407  -1.333  -1.916  1.00  0.00           C
HETATM   35  C3  1ZN A   3      -5.116  -2.581  -2.743  1.00  0.00           C
HETATM   36  C4  1ZN A   3      -6.159  -3.648  -3.029  1.00  0.00           C
HETATM   37  C5  1ZN A   3      -7.053  -4.123  -2.061  1.00  0.00           C
HETATM   38  C6  1ZN A   3      -8.144  -4.905  -2.440  1.00  0.00           C
HETATM   39  C7  1ZN A   3      -8.327  -5.247  -3.779  1.00  0.00           C
HETATM   40  C8  1ZN A   3      -7.409  -4.834  -4.738  1.00  0.00           C
HETATM   41  C9  1ZN A   3      -6.323  -4.053  -4.364  1.00  0.00           C
HETATM   42  C10 1ZN A   3      -4.576  -0.178  -2.405  1.00  0.00           C
HETATM   43  C11 1ZN A   3      -5.053   0.691  -3.523  1.00  0.00           C
HETATM   44  C12 1ZN A   3      -4.103   0.686  -1.290  1.00  0.00           C
HETATM   45  C13 1ZN A   3      -2.835   1.056  -1.225  1.00  0.00           C
HETATM   46  C14 1ZN A   3      -1.920   0.679  -2.352  1.00  0.00           C
HETATM   47  C15 1ZN A   3      -2.429   1.676  -0.130  1.00  0.00           C
HETATM   48  C16 1ZN A   3      -1.153   1.597   0.239  1.00  0.00           C
HETATM   49  CA  1ZN A   3      -0.669   2.408   1.336  1.00  0.00           C
HETATM   50  N   1ZN A   3       0.187   3.520   0.932  1.00  0.00           N
HETATM   51  C18 1ZN A   3       0.018   1.612   2.437  1.00  0.00           C
HETATM   52  C19 1ZN A   3      -0.875   0.772   3.283  1.00  0.00           C
HETATM   53  C   1ZN A   3       0.830   2.579   3.261  1.00  0.00           C
HETATM   54  O   1ZN A   3       0.489   3.730   3.569  1.00  0.00           O
HETATM   82  N   FGA A   4       1.905   2.010   3.721  1.00  0.00           N
HETATM   83  CA  FGA A   4       2.861   2.665   4.609  1.00  0.00           C
HETATM   84  C   FGA A   4       2.978   1.832   5.914  1.00  0.00           C
HETATM   85  O   FGA A   4       4.055   1.309   6.252  1.00  0.00           O
HETATM   86  CB  FGA A   4       4.163   2.893   3.900  1.00  0.00           C
HETATM   87  CG  FGA A   4       4.920   1.738   3.242  1.00  0.00           C
HETATM   88  CD  FGA A   4       5.387   1.427   1.819  1.00  0.00           C
HETATM   89  OE1 FGA A   4       4.268   1.471   1.286  1.00  0.00           O
HETATM   96  N   MDH A   5       6.206   1.058   0.852  1.00  0.00           N
HETATM   97  CM  MDH A   5       7.489   1.578   1.487  1.00  0.00           C
HETATM   98  CA  MDH A   5       6.631   0.777  -0.516  1.00  0.00           C
HETATM   99  C   MDH A   5       5.615   1.660  -1.243  1.00  0.00           C
HETATM  100  O   MDH A   5       6.060   1.412  -2.374  1.00  0.00           O
HETATM  101  CB  MDH A   5       7.635  -0.135  -0.796  1.00  0.00           C
HETATM  102  CG  MDH A   5       8.551  -0.894   0.158  1.00  0.00           C
END
"""

pdb_str_022_control = """
MODEL        1
HETATM    1  N   CSU A   1     -17.038  -8.727  -8.532  1.00  1.00           N
HETATM    2  CA  CSU A   1     -17.363  -7.323  -8.738  1.00  1.00           C
HETATM    3  CB  CSU A   1     -18.453  -7.183  -9.802  1.00  1.00           C
HETATM    4  SG  CSU A   1     -18.918  -5.438  -9.894  1.00  1.00           S
HETATM    5  S   CSU A   1     -20.321  -5.607 -11.354  1.00  1.00           S
HETATM    6  C   CSU A   1     -17.784  -6.690  -7.414  1.00  1.00           C
HETATM    7  O   CSU A   1     -17.096  -5.821  -6.882  1.00  1.00           O
HETATM    8  OD1 CSU A   1     -20.774  -4.175 -11.430  1.00  1.00           O
HETATM    9  OD2 CSU A   1     -21.438  -6.382 -10.849  1.00  1.00           O
HETATM   10  OD3 CSU A   1     -19.675  -5.949 -12.608  1.00  1.00           O
ATOM     18  N   SER A   2     -18.871  -7.186  -6.818  1.00  1.00           N
ATOM     19  CA  SER A   2     -19.308  -6.763  -5.491  1.00  1.00           C
ATOM     20  C   SER A   2     -18.499  -7.499  -4.413  1.00  1.00           C
ATOM     21  O   SER A   2     -19.057  -7.929  -3.403  1.00  1.00           O
ATOM     22  CB  SER A   2     -20.819  -7.011  -5.355  1.00  1.00           C
ATOM     23  OG  SER A   2     -21.285  -6.562  -4.099  1.00  1.00           O
ENDMDL
MODEL        2
HETATM    1  N   CSU A   1     -16.367  -8.724  -8.446  1.00  1.00           N
HETATM    2  CA  CSU A   1     -17.541  -8.173  -9.115  1.00  1.00           C
HETATM    3  CB  CSU A   1     -18.544  -9.280  -9.455  1.00  1.00           C
HETATM    4  SG  CSU A   1     -17.916 -10.298 -10.816  1.00  1.00           S
HETATM    5  S   CSU A   1     -19.479 -11.591 -10.932  1.00  1.00           S
HETATM    6  C   CSU A   1     -18.190  -7.113  -8.230  1.00  1.00           C
HETATM    7  O   CSU A   1     -18.527  -6.023  -8.686  1.00  1.00           O
HETATM    8  OD1 CSU A   1     -19.506 -12.426  -9.745  1.00  1.00           O
HETATM    9  OD2 CSU A   1     -20.659 -10.869 -11.368  1.00  1.00           O
HETATM   10  OD3 CSU A   1     -19.020 -12.433 -12.091  1.00  1.00           O
ATOM     18  N   SER A   2     -18.353  -7.420  -6.943  1.00  1.00           N
ATOM     19  CA  SER A   2     -18.701  -6.432  -5.922  1.00  1.00           C
ATOM     20  C   SER A   2     -18.054  -6.828  -4.593  1.00  1.00           C
ATOM     21  O   SER A   2     -18.567  -6.548  -3.511  1.00  1.00           O
ATOM     22  CB  SER A   2     -20.227  -6.316  -5.805  1.00  1.00           C
ATOM     23  OG  SER A   2     -20.801  -6.046  -7.072  1.00  1.00           O
ENDMDL
END
"""

pdb_str_023 = """
ATOM   1944  N   VAL C 504      25.617   7.043  38.901  1.00 21.91           N
ATOM   1945  CA  VAL C 504      25.298   6.520  37.593  1.00 25.45           C
ATOM   1946  C   VAL C 504      26.023   5.185  37.459  1.00 27.80           C
ATOM   1947  O   VAL C 504      26.017   4.447  38.451  1.00 32.85           O
ATOM   1948  CB  VAL C 504      25.759   7.467  36.446  1.00 22.73           C
ATOM   1949  CG1 VAL C 504      25.062   8.814  36.558  1.00 21.23           C
ATOM   1950  CG2 VAL C 504      27.268   7.641  36.479  1.00 21.12           C
HETATM 1951  C   ASJ C 505      25.757   3.265  34.352  1.00 21.67           C
HETATM 1952  N   ASJ C 505      25.878   4.624  36.266  1.00 30.84           N
HETATM 1953  O   ASJ C 505      25.368   2.035  33.737  1.00 21.73           O
HETATM 1954  CA  ASJ C 505      25.449   3.280  35.850  1.00 24.76           C
HETATM 1955  CB  ASJ C 505      23.949   3.073  36.085  1.00 25.50           C
HETATM 1956  CG  ASJ C 505      23.651   1.958  37.095  1.00 22.63           C
HETATM 1957  OD1 ASJ C 505      24.572   1.633  37.874  1.00 18.73           O
HETATM 1958  OD2 ASJ C 505      22.458   1.603  37.218  1.00 20.52           O
END
"""

if __name__ == '__main__':
  run()
  print("OK")
