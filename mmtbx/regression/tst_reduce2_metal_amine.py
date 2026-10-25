"""reduce2 metal_bound_amines_as_nh2 (A5): placement and geometry. A primary
amine that the dictionary makes NH3+ (Lys NZ, an N-terminal NH3) and that
coordinates a metal of AMINE_METAL_ELEMENTS (Mg included: 1MNZ Tris; not Na,
K, Ca) is placed as NH2, the two H oriented with the metal in place of the H
not placed; the classifier, metal-bound chain starts, geometry checks,
symmetry images, hd_and_h, the qmready charge and the option off. Supplied H
and reruns: tst_reduce2_metal_amine_supplied.py.
"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression.reduce2_metal_amine_utils import *
from mmtbx.regression.reduce2_metal_amine_utils import shared, shared_text

def test_lys_cd():
  h, description, log = shared('lys_on')
  check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])
  for text in (description, log):
    assert section(text, CONVERTED) == ['NZ  LYS A 359: HZ1 not placed  '
      '(NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'], text
    assert REPORTED not in text
  # the coordinates place_metal_nh2() accepts are the final ones
  from mmtbx.hydrogens.reduce_hydrogen import place_metal_nh2
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  placed, reason = place_metal_nh2(n, ce, cd, lys_targets())
  assert reason is None, reason
  for name, xyz in zip(('HZ2', 'HZ3'), placed):
    d = (matrix.col(atom(h, 359, name).xyz) - xyz).length()
    assert d < 2.e-3, (name, d)

def test_option_off():
  h, description, log = shared('lys_off')
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  assert CONVERTED not in description and CONVERTED not in log
  assert section(description, REPORTED)[-1].startswith(
    'NZ  LYS A 359 (HZ1, HZ2, HZ3)'), description

def test_tris_mg():
  '''Mg counts (NH3+ has no lone pair, and an N-Mg contact is too short for
  an N-H...O bond): the 1MNZ Tris amine on Mg is placed as NH2, HN1 not
  placed'''
  h, description, log = run('tst_reduce2_metal_amine_trs', pdb_trs_mg, [ON])
  check_nh2(h, 392, 'N', 'C', 389, 'MG', ['HN2', 'HN3'])
  assert section(description, CONVERTED) == ['N   TRS A 392: HN1 not placed  '
    '(N   TRS A 392 - MG    MG A 389 2.26 A (geometric))'], description

def test_lys_na_k_ca():
  '''Na, K and Ca do not count: with the option, Lys NZ 2.21 A from one of
  them (the Cd of pdb_lys_cd replaced) keeps HZ1, HZ2 and HZ3'''
  for e in ('NA', 'K', 'CA'):
    pdb_str = pdb_lys_cd.replace('CD    CD A 396', '%-2s    %-2s A 396' % (
      e, e)).replace('10.07          CD', '10.07          %2s' % e)
    assert ' CD A' not in pdb_str, pdb_str
    h, description, log = run('tst_reduce2_metal_amine_lys_' + e, pdb_str,
      [ON])
    assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3'], e
    assert CONVERTED not in description, (e, description)

def test_terminus_cu():
  # resseq 1: placement builds NH3 (n_terminal_charge=residue_one)
  h, description, log = run('tst_reduce2_metal_amine_his', pdb_his_cu, [ON])
  check_nh2(h, 1, 'N', 'CA', 301, 'CU', ['H2', 'H3'])
  assert section(description, CONVERTED) == ['N   HIS A   1: H1 not placed  '
    '(N   HIS A   1 - CU    CU A 301 2.13 A (geometric))'], description

def test_neutral_nh2_kept():
  # negative control through reduce2: Cu(en) keeps two H per N
  h, description, log = run('tst_reduce2_metal_amine_edn', pdb_edn_cu, [ON])
  assert [n for n, alt in h_on(h, 1, 'NAA')] == ['HAA1', 'HAA2']
  assert [n for n, alt in h_on(h, 1, 'NAD')] == ['HAD1', 'HAD2']
  assert CONVERTED not in description

def test_classifier():
  # which N the conversion takes (metal_amine_sites.amine_class(), as
  # convert_metal_bound_amines() uses it; no reduce2 run): a primary amine
  # (Lys NZ) that its dictionary makes NH3+; Cu(en) is a primary amine but
  # neutral NH2 (both H kept); Asn ND2 is planar; YI1 NAV is an O-NH3+
  from mmtbx.hydrogens import reduce_hydrogen
  m = classifier_model()
  found = {}
  class probe(reduce_hydrogen.metal_amine_sites):
    def _candidates(self, atoms, elements, heavy, h_on):
      found['heavy'] = heavy
      return []
  sites = probe(m)
  atoms = m.get_hierarchy().atoms()
  srv = m.get_mon_lib_srv()
  def judge(chain, name):
    i = [a.i_seq for a in atoms if a.name.strip() == name and
         a.parent().parent().parent().id == chain][0]
    branch = sites.amine_class(atoms, found['heavy'], i, i).branch
    return branch, reduce_hydrogen.dictionary_nh3_plus(m, atoms[i], i, srv)
  assert judge('D', 'NZ') == (reduce_hydrogen.PRIMARY_AMINE, True)
  for name in ('NAA', 'NAD'):
    # Cu(en): neutral NH2, kept with two H
    assert judge('A', name) == (reduce_hydrogen.PRIMARY_AMINE, False), name
  branch, plus = judge('B', 'ND2')
  assert branch != reduce_hydrogen.PRIMARY_AMINE, branch
  branch, plus = judge('C', 'NAV')
  assert branch == 'heavy neighbour OAU is not C', branch
  # chain starts (place_n_terminal_propeller()): a primary amine by the
  # dictionary; not DHA (CA=CB), SAR (N-methyl) or PRO
  p = reduce_hydrogen.place_hydrogens(m)
  p.query_mon_lib_srv = srv
  starts = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1       1.000   1.000   1.000  1.00 20.00           N
ATOM      2  N   DHA A   2       5.000   1.000   1.000  1.00 20.00           N
ATOM      3  N   SAR A   3       9.000   1.000   1.000  1.00 20.00           N
ATOM      4  N   PRO A   4      13.000   1.000   1.000  1.00 20.00           N
""").construct_hierarchy()
  result = dict((a.parent().resname, p._dictionary_primary_amine(a))
    for a in starts.atoms())
  assert result['GLY'] == (True, ''), result
  assert result['DHA'] == (False, 'not a primary amine (CA in a multiple or '
    'aromatic bond)'), result
  assert result['SAR'][0] is False and result['PRO'][0] is False, result

def test_chain_start_settings():
  # Gly 288 starts the chain in the model: NH2 whatever n_terminal_charge
  # says (no sequence records are read)
  for setting in ('residue_one', 'no_charge', 'all_chain_starts'):
    if setting == 'residue_one':
      # the default
      h, description, log = shared('x1_on')
    else:
      h, description, log = run('tst_reduce2_metal_amine_gly_' + setting,
        pdb_1xmk, [ON, 'n_terminal_charge=' + setting])
    assert h_on(h, 288, 'N') == [('H2', ''), ('H3', '')], (setting,
      h_on(h, 288, 'N'))
    check_nh2(h, 288, 'N', 'CA', 398, 'NI', ['H2', 'H3'])
    expected = (GLY288 % setting if setting != 'all_chain_starts' else
      'N GLY A 288: metal-bound chain start, placed as NH2')
    for text in (description, log):
      assert section(text, CHAIN_START) == [expected], (setting, text)

def test_chain_start_preceding_residue():
  # a residue before Gly 288 with its C bonded to N: no chain start, an amide
  # on the Ni (step B: no H). Passes without C1c.
  pdb = pdb_1xmk.replace('ATOM      1  N   GLY A 288', """\
ATOM     41  N   GLY A 287      11.272  25.005  34.136  1.00 20.00           N
ATOM     42  CA  GLY A 287      11.618  24.171  32.989  1.00 20.00           C
ATOM     43  C   GLY A 287      12.237  24.201  31.601  1.00 20.00           C
ATOM     44  O   GLY A 287      12.648  25.262  31.134  1.00 20.00           O
ATOM      1  N   GLY A 288""")
  h, description, log = run('tst_reduce2_metal_amine_gly_prev', pdb, [ON])
  assert h_on(h, 288, 'N') == [], h_on(h, 288, 'N')
  assert 'metal-bound chain start' not in description
  assert 'H   GLY A 288  (N   GLY A 288 - NI' in '\n'.join(
    section(description, 'N-H not placed')), description

def test_chain_start_refused():
  # two Ni positions (A, B): no unique frame, NH3 kept with the reason
  pdb = pdb_1xmk.replace(
    'HETATM   21 NI    NI A 398      13.849  21.950  30.990  1.00',
    'HETATM   21 NI  A NI A 398      13.849  21.950  30.990  0.50').replace(
    'END\n', 'HETATM   22 NI  B NI A 398      13.849  22.450  30.990  0.50 20.00'
    '          NI\nEND\n')
  h, description, log = run('tst_reduce2_metal_amine_gly_alt', pdb, [ON])
  assert [x for x, alt in h_on(h, 288, 'N')] == ['H1', 'H2', 'H3'], \
    h_on(h, 288, 'N')
  lines = section(description, CHAIN_START)
  assert len(lines) == 1 and lines[0].startswith('N GLY A 288: metal-bound '
    'chain start, NH3 kept: no unique metal position') and lines[0].endswith(
    '(n_terminal_charge=residue_one would place 1 H)'), lines

def test_chain_start_without_metal():
  # no metal: as n_terminal_charge places it (passes without C1c)
  pdb = '\n'.join(l for l in pdb_1xmk.splitlines() if ' NI ' not in l) + '\n'
  for setting in ('residue_one', 'no_charge'):
    h, description, log = run('tst_reduce2_metal_amine_gly_nometal', pdb,
      [ON, 'n_terminal_charge=' + setting])
    assert h_on(h, 288, 'N') == [('H', '')], (setting, h_on(h, 288, 'N'))
    assert CHAIN_START not in description

def test_chain_start_option_off():
  # without the option: one H on Gly 288 (residue_one), as on master
  h, description, log = shared('x1_off')
  assert h_on(h, 288, 'N') == [('H', '')]
  assert CHAIN_START not in description and CONVERTED not in description

def test_chain_start_conformers():
  # Gly 288 N in A and B, the Ni in one of them: that conformer NH2, the
  # other one H (residue_one), the line names the conformer
  for metal_alt, other in (('B', 'A'), ('A', 'B')):
    h, d, l = run('tst_reduce2_metal_amine_gc' + metal_alt,
      gly288_conformers(metal_alt), [ON])
    n = dict((alt, atom(h, 288, 'N', alt)) for alt in ('A', 'B'))
    on = lambda alt: sorted(x.name.strip() for x in h.atoms()
      if x.element.strip() == 'H' and x.parent().altloc == alt
      and x.distance(n[alt]) < 1.2)
    assert on(metal_alt) == ['H2', 'H3'], (metal_alt, on(metal_alt))
    assert on(other) == ['H'], (metal_alt, on(other))
    assert section(d, CHAIN_START) == [GLY288.replace('N GLY A 288:',
      'N GLY A 288 conformer %s:' % metal_alt) % 'residue_one'], d
  # a blank N with the Ni only in B: no common frame, unchanged
  h, d, l = run('tst_reduce2_metal_amine_gcn', gly288_conformers('B', False),
    [ON])
  assert [x for x, alt in h_on(h, 288, 'N')] == ['H', 'H'], h_on(h, 288, 'N')
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: metal only in conformer B (no common frame)'], d

def test_unresolved():
  # two alternate Cd positions: full coverage, but no single frame
  pdb = pdb_lys_cd.replace(
    'HETATM   10 CD    CD A 396      10.293  10.651  44.885  1.00',
    'HETATM   10 CD  A CD A 396      10.293  10.651  44.885  0.50').replace(
    'END\n', 'HETATM   11 CD  B CD A 396      11.011  10.159  44.542  0.50 10.07'
    '          CD\nEND\n')
  h, description, log = run('tst_reduce2_metal_amine_alt', pdb, [ON])
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  lines = section(description, UNRESOLVED)
  assert len(lines) == 1 and lines[0].startswith(
    'NZ  LYS A 359: no unique metal position'), lines

def test_geometry_not_placeable():
  # C-N-Cd 140 and 145 deg: no valid NH2 (place_metal_nh2() refuses both).
  # In the pipeline (140): no H is removed, new H as without the option,
  # supplied H as they are
  from mmtbx.hydrogens.reduce_hydrogen import place_metal_nh2
  h_lys = processed(pdb_lys_cd).get_hierarchy()
  for angle in (140, 145):
    h = iotbx.pdb.input(source_info=None,
      lines=lys_cd_at(pdb_lys_cd, angle)).construct_hierarchy()
    n, ce, cd = [matrix.col(atom(h, r, x).xyz)
      for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
    placed, reason = place_metal_nh2(n, ce, cd, lys_targets())
    assert placed is None and reason.startswith('geometry not placeable, '
      'C-N-M %d deg' % angle), (angle, reason)
  nh2 = shared_text('lys_on')
  angle = 140
  pdb = lys_cd_at(pdb_lys_cd, angle)
  off, d_off, l_off = run('tst_reduce2_metal_amine_g%d_off' % angle, pdb,
    [OFF])
  h, description, log = run('tst_reduce2_metal_amine_g%d' % angle, pdb, [ON])
  assert CONVERTED not in description
  lines = section(description, UNRESOLVED)
  assert len(lines) == 1 and lines[0].startswith('NZ  LYS A 359: geometry '
    'not placeable, C-N-M %d deg' % angle), lines
  assert h_on(h, 359, 'NZ') == h_on(off, 359, 'NZ') == [('HZ1', ''),
    ('HZ2', ''), ('HZ3', '')]
  for name in ('HZ1', 'HZ2', 'HZ3'):
    assert atom(h, 359, name).distance(atom(off, 359, name)) < 1.e-3, name
  # supplied NH2
  supplied = lys_cd_at(nh2, angle)
  h0 = iotbx.pdb.input(source_info=None, lines=supplied).construct_hierarchy()
  h, description, log = run('tst_reduce2_metal_amine_g%d_keep' % angle,
    supplied, [ON, 'keep_existing_H=True'])
  assert h_on(h, 359, 'NZ') == [('HZ2', ''), ('HZ3', '')], h_on(h, 359, 'NZ')
  for name in ('HZ2', 'HZ3'):
    assert atom(h, 359, name).distance(atom(h0, 359, name)) < 1.e-3, name
  # kept as supplied: not placed, not checked
  assert UNRESOLVED not in description and CONVERTED not in description
  assert section(description, SUPPLIED_NH2)[0].startswith(
    'NZ  LYS A 359 (HZ2, HZ3)'), description

def test_only_reference_angles_removed():
  # an angle already on C-N-M (here 107 deg, weight 123, origin edits)
  # survives placement unchanged
  import mmtbx.model
  from libtbx.utils import null_out
  from cctbx import geometry_restraints
  from mmtbx.hydrogens import reduce_hydrogen
  class placement(reduce_hydrogen.place_hydrogens):
    def add_link_h_restraints(self):
      reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)
      atoms = self.model.get_hierarchy().atoms()
      i = dict((a.name.strip(), a.i_seq) for a in atoms)
      self.model.get_restraints_manager().geometry.add_angles_in_place([
        geometry_restraints.angle_proxy(i_seqs = [i['CE'], i['NZ'], i['CD']],
          angle_ideal = 107, weight = 123, origin_id = 4)])
  m = mmtbx.model.manager(model_input = iotbx.pdb.input(source_info=None,
    lines=pdb_lys_cd), log = null_out())
  p = placement(m, metal_bound_amines_as_nh2 = True)
  p.run()
  assert p.metal_amine_converted_text, p.metal_amine_unresolved
  m = p.get_model()
  i = dict((a.name.strip(), a.i_seq) for a in m.get_hierarchy().atoms())
  found = [(round(x.angle_ideal, 3), x.weight, x.origin_id)
    for x in m.get_restraints_manager().geometry.angle_proxies
    if set(x.i_seqs) == set([i['CE'], i['NZ'], i['CD']])]
  assert found == [(107, 123, 4)], found
  # and the injected H-N-M angles are gone
  assert len([x for x in m.get_restraints_manager().geometry.angle_proxies
    if i['CD'] in x.i_seqs]) == 1

def test_lists_after_gc():
  # the lists are built at the end of placement from what step B saw; the
  # atoms of that earlier hierarchy may be collected by then (shared x1_off
  # runs with a gc.collect() just before the lists are built)
  h, description, log = shared('x1_off')
  assert section(description, REPORTED) == [
    'N   GLY A 288 (H)  (N   GLY A 288 - NI    NI A 398 1.89 A (geometric))'
    ], description

def test_other_h_unchanged():
  # with a converted amine, every H away from it as without the option
  # (placement's riding H manager is kept; a new one moved the glycan H)
  off, d, l = run('tst_reduce2_metal_amine_gly_off', pdb_his_cu_glycan,
    [OFF])
  on, d, l = run('tst_reduce2_metal_amine_gly_on', pdb_his_cu_glycan, [ON])
  assert section(d, CONVERTED)[0].startswith('N   HIS A   1: H1 not placed')
  n = atom(on, 1, 'N')
  xyz = lambda h: dict((a.id_str(), a.xyz) for a in h.atoms()
    if a.element.strip() == 'H' and a.distance(n) > 1.2)
  a, b = xyz(off), xyz(on)
  assert sorted(a) == sorted(b)
  moved = [k for k in a if (matrix.col(a[k]) - matrix.col(b[k])).length()
           > 1.e-3]
  assert not moved, moved

def test_unequal_targets():
  # C-N-HZ2 105, C-N-HZ3 115 deg and unequal N-H: place_metal_nh2() accepts
  # them with each H's own targets (test_lys_cd: the accepted coordinates are
  # the final ones)
  from mmtbx.hydrogens.reduce_hydrogen import place_metal_nh2
  h = shared('lys_on')[0]
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  placed, reason = place_metal_nh2(n, ce, cd, lys_targets(d = [0.89, 0.95],
    cnh = [105., 115.]))
  assert reason is None, reason
  assert abs((placed[0] - n).length() - 0.89) < 1.e-6
  assert abs((placed[1] - n).length() - 0.95) < 1.e-6
  assert (ce - n).angle(placed[0] - n, deg=True) < \
    (ce - n).angle(placed[1] - n, deg=True)

def test_unequal_targets_placement():
  # the same through reduce2: a LYS dictionary with C-N-HZ2 105, C-N-HZ3 115
  # deg and N-HZ3 0.95 A; each H keeps its own targets and the accepted
  # coordinates are the final ones (riding H would give both H one set)
  lys_cif('tst_reduce2_metal_amine_lys.cif', (
      (' LYS  NZ  HZ3  coval  0.890  0.020  1.040',
       ' LYS  NZ  HZ3  coval  0.950  0.020  1.100'),
      (' LYS  CE   NZ  HZ2  110.000  3.000', ' LYS  CE   NZ  HZ2  105.000  3.000'),
      (' LYS  CE   NZ  HZ3  110.000  3.000', ' LYS  CE   NZ  HZ3  115.000  3.000')))
  h, description, log = run('tst_reduce2_metal_amine_unequal', pdb_lys_cd,
    [ON, 'tst_reduce2_metal_amine_lys.cif'])
  assert section(description, CONVERTED) == ['NZ  LYS A 359: HZ1 not placed  '
    '(NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'], description
  from mmtbx.hydrogens.reduce_hydrogen import place_metal_nh2
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  placed, reason = place_metal_nh2(n, ce, cd, lys_targets(d = [0.89, 0.95],
    cnh = [105., 115.]))
  assert reason is None, reason
  for name, xyz in zip(('HZ2', 'HZ3'), placed):
    d = (matrix.col(atom(h, 359, name).xyz) - xyz).length()
    assert d < 2.e-3, (name, d)
  assert abs(atom(h, 359, 'HZ3').distance(atom(h, 359, 'NZ')) - 0.95) < 2.e-3
  assert angle(atom(h, 359, 'CE'), atom(h, 359, 'NZ'), atom(h, 359, 'HZ2')) < \
    angle(atom(h, 359, 'CE'), atom(h, 359, 'NZ'), atom(h, 359, 'HZ3'))

def test_symmetry_contact():
  # the same Cd through a translation and through a rotation
  # (-x+1/2,-y,z+1/2): the contact's image is the Cd itself, and the H placed
  # from it are those of the Cd itself; the rotation also through reduce2
  from cctbx import crystal, sgtbx
  from mmtbx.hydrogens import reduce_hydrogen
  h0 = shared('lys_on')[0]
  cs = crystal.symmetry(unit_cell=(35.557, 43.526, 45.471, 90, 90, 90),
    space_group_symbol='P 21 21 21')
  uc = cs.unit_cell()
  cd = matrix.col(atom(h0, 396, 'CD').xyz)
  n, ce = [matrix.col(atom(h0, 359, x).xyz) for x in ('NZ', 'CE')]
  for tag, op in (('t', 'x+1,y,z'), ('r', '-x+1/2,-y,z+1/2')):
    rt = sgtbx.rt_mx(op)
    xyz = uc.orthogonalize(rt.inverse() * uc.fractionalize(cd.elems))
    pdb = pdb_lys_cd.replace('CD A 396      10.293  10.651  44.885',
      'CD A 396    %8.3f%8.3f%8.3f' % tuple(xyz))
    hp = iotbx.pdb.input(source_info=None, lines=pdb).construct_hierarchy()
    atoms = hp.atoms()
    el = [a.element.strip().upper() for a in atoms]
    i_n = [i for i, a in enumerate(atoms) if a.name.strip() == 'NZ'][0]
    i_m = [i for i, a in enumerate(atoms) if el[i] == 'CD'][0]
    pairs = reduce_hydrogen.metal_contact_pairs([a.xyz for a in atoms], el,
      [i_n], [i_m], cs)
    images = []
    for i_x, j, r in pairs:
      image = matrix.col(uc.orthogonalize(r * uc.fractionalize(atoms[j].xyz)))
      if not [x for x in images if (x - image).length() < 0.01]:
        images.append(image)
    assert len(images) == 1 and (images[0] - cd).length() < 1.e-3, (tag,
      images)
    placed, reason = reduce_hydrogen.place_metal_nh2(n, ce, images[0],
      lys_targets())
    assert reason is None, (tag, reason)
    for name, x in zip(('HZ2', 'HZ3'), placed):
      d = (matrix.col(atom(h0, 359, name).xyz) - x).length()
      assert d < 2.e-3, (tag, name, d)
    if tag == 'r':
      h, description, log = run('tst_reduce2_metal_amine_sym' + tag, pdb,
        [ON])
      lines = section(description, CONVERTED)
      assert len(lines) == 1 and ' sym ' in lines[0], (tag, lines)
      for name in ('HZ2', 'HZ3'):
        d = atom(h, 359, name).distance(atom(h0, 359, name))
        assert d < 2.e-3, (tag, name, d)

def test_shared_frame_restraints():
  # blank N, C and Cd, H groups A and B (O split): both converted, and the
  # angle restraints of the rest, with one already on C-N-M, unchanged in
  # count and metadata
  import mmtbx.model
  from libtbx.utils import null_out
  from cctbx import geometry_restraints
  from mmtbx.hydrogens import reduce_hydrogen
  class placement(reduce_hydrogen.place_hydrogens):
    def add_link_h_restraints(self):
      reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)
      i = dict((a.name.strip(), a.i_seq)
        for a in self.model.get_hierarchy().atoms())
      self.model.get_restraints_manager().geometry.add_angles_in_place([
        geometry_restraints.angle_proxy(i_seqs = [i['CE'], i['NZ'], i['CD']],
          angle_ideal = 107, weight = 123, origin_id = 4)])
  pdb = pdb_lys_cd.replace(
    'ATOM      4  O   LYS A 359       7.955   4.152  45.067  1.00',
    'ATOM      4  O  ALYS A 359       7.955   4.152  45.067  0.50').replace(
    'ATOM      5  CB', 'ATOM      4  O  BLYS A 359       7.855   4.252  45.067'
    '  0.50  7.65           O\nATOM      5  CB')
  result = []
  for on in (False, True):
    m = mmtbx.model.manager(model_input = iotbx.pdb.input(source_info=None,
      lines=pdb), log = null_out())
    p = placement(m, metal_bound_amines_as_nh2 = on)
    p.run()
    m = p.get_model()
    atoms = m.get_hierarchy().atoms()
    nz_h = set(a.i_seq for a in atoms if a.name.strip().startswith('HZ'))
    label = lambda i: atoms[i].id_str()
    result.append((sorted((tuple(label(i) for i in x.i_seqs), x.angle_ideal,
      x.weight, x.slack, x.origin_id, x.sym_ops is None)
      for x in m.get_restraints_manager().geometry.angle_proxies
      if not nz_h & set(x.i_seqs)), p, nz_h, atoms))
  assert result[0][0] == result[1][0]
  p, nz_h, atoms = result[1][1:]
  assert sorted((atoms[i].name.strip(), atoms[i].parent().altloc)
    for i in nz_h) == [('HZ2', 'A'), ('HZ2', 'B'), ('HZ3', 'A'), ('HZ3', 'B')]
  assert len(p.metal_amine_converted_text) == 2, p.metal_amine_converted_text

def test_alternate_metals_per_conformer():
  # NZ (or only CE) and the Cd in A and B: each conformer NH2, oriented to
  # its own Cd (a guard: passes without C1d)
  for tag, names in (('nz', ('CE', 'NZ')), ('ce', ('CE',))):
    h, d, l = run('tst_reduce2_metal_amine_alt_' + tag, lys_cd_alternates(
      names), [ON])
    assert len(section(d, CONVERTED)) == 2, (tag, d)
    a = dict(((x.name.strip(), x.parent().altloc), x) for x in h.atoms())
    for alt in 'AB':
      n = a.get(('NZ', alt)) or a[('NZ', '')]
      ce, cd = a[('CE', alt)], a[('CD', alt)]
      hs = [a[(x, alt)] for x in ('HZ2', 'HZ3')]
      assert not ('HZ1', alt) in a, (tag, alt)
      for x in hs:
        assert 104 < angle(cd, n, x) < 115, (tag, alt, angle(cd, n, x))
        assert 104 < angle(ce, n, x) < 115, (tag, alt, angle(ce, n, x))
        assert x.distance(cd) > 2.3, (tag, alt, x.distance(cd))
      assert 104 < angle(hs[0], n, hs[1]) < 115, (tag, alt)

def test_special_position():
  # a Cd on a 2-fold (and 0.04 A off it): one metal position, one contact,
  # and a valid NH2 from it (a guard). Model processing puts an atom whose
  # symmetry-equivalent images are closer than min_distance_sym_equiv (0.5 A;
  # for a 2-fold 0.25 A from the axis) onto the special position, so the
  # conversion sees one image
  from mmtbx.hydrogens import reduce_hydrogen
  for tag, pdb in (('on', pdb_lys_cd_2fold), ('off', pdb_lys_cd_2fold.replace(
      'CD A 396       0.000', 'CD A 396       0.040'))):
    m = processed(pdb)
    cs = m.crystal_symmetry()
    uc = cs.unit_cell()
    atoms = m.get_hierarchy().atoms()
    el = [a.element.strip().upper() for a in atoms]
    i_n = [i for i, a in enumerate(atoms) if a.name.strip() == 'NZ'][0]
    i_m = [i for i, a in enumerate(atoms) if el[i] == 'CD'][0]
    images = []
    for i_x, j, r in reduce_hydrogen.metal_contact_pairs(
        [a.xyz for a in atoms], el, [i_n], [i_m], cs):
      image = matrix.col(uc.orthogonalize(r * uc.fractionalize(atoms[j].xyz)))
      if not [x for x in images if (x - image).length() < 0.01]:
        images.append(image)
    assert len(images) == 1, (tag, images)
    assert abs(images[0][0]) < 1.e-3, (tag, images[0])
    n, ce = [matrix.col(atoms[i].xyz) for i in (i_n, [i for i, a in
      enumerate(atoms) if a.name.strip() == 'CE'][0])]
    placed, reason = reduce_hydrogen.place_metal_nh2(n, ce, images[0],
      lys_targets())
    assert reason is None, (tag, reason)
    for x in placed:
      assert (x - images[0]).length() > 2.3, (tag, (x - images[0]).length())
  # through reduce2 (its own contacts): the Cd on the axis, one contact, NH2
  h, d, l = run('tst_reduce2_metal_amine_2fold_on', pdb_lys_cd_2fold, [ON])
  lines = section(d, CONVERTED)
  assert len(lines) == 1 and lines[0].startswith('NZ  LYS A 359: HZ1 not '
    'placed  (NZ  LYS A 359 - CD    CD A 396 2.') and ';' not in lines[0], lines
  check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])

def test_hd_and_h():
  h, description, log = run('tst_reduce2_metal_amine_hd', pdb_lys_cd,
    [ON, 'neutron_option=hd_and_h'])
  assert h_on(h, 359, 'NZ') == [('DZ2', 'B'), ('DZ3', 'B'), ('HZ2', 'A'),
    ('HZ3', 'A')], h_on(h, 359, 'NZ')

def test_qmready_charge():
  # one converted site: total charge one lower
  pdb = pdb_lys_cd.replace('          CD\n', '          CD2+\n')
  charges = []
  for args in ([OFF], [ON]):
    h, description, log, r = run('tst_reduce2_metal_amine_qm', pdb,
      ['qmready=True'] + args, result=True)
    charges.append(r.total_charge)
  assert charges[1] == charges[0] - 1, charges

if __name__ == '__main__':
  test_lys_cd()
  test_option_off()
  test_tris_mg()
  test_lys_na_k_ca()
  test_terminus_cu()
  test_neutral_nh2_kept()
  test_classifier()
  test_chain_start_settings()
  test_chain_start_preceding_residue()
  test_chain_start_refused()
  test_chain_start_without_metal()
  test_chain_start_option_off()
  test_chain_start_conformers()
  test_unresolved()
  test_geometry_not_placeable()
  test_only_reference_angles_removed()
  test_lists_after_gc()
  test_other_h_unchanged()
  test_unequal_targets()
  test_unequal_targets_placement()
  test_symmetry_contact()
  test_shared_frame_restraints()
  test_alternate_metals_per_conformer()
  test_special_position()
  test_hd_and_h()
  test_qmready_charge()
  print('OK')
