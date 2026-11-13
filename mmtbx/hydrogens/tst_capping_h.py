"""
Capping H as first-class atoms in cctbx: recognition by pdb_interpretation,
the `capping` origin id, restraints and consumers of the cap set.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
import mmtbx.refinement.refinement_flags
from libtbx.utils import null_out
from libtbx import group_args
from six.moves import cStringIO as StringIO
from scitbx import matrix
from scitbx.array_family import flex
from mmtbx.hydrogens import capping
from mmtbx.hydrogens import tst_capping as fixtures

def test_origin_id():
  from cctbx.geometry_restraints.linking_class import linking_class
  from cctbx.geometry_restraints.auto_linking_types import iterate_covalent
  ids = linking_class()
  i = ids.get_origin_id('capping')
  assert i > 0
  assert ids.data['capping'][3] == ['Capping H', 'Capping H']
  assert 'capping' not in list(iterate_covalent())

def h_names_on(model, resseq, name, chain_id='A'):
  # fixtures.h_names_on builds the bond table of the model's restraints,
  # which a model with cap H appended no longer matches: use the geometry
  parent = fixtures.atom_named(model, resseq, name, chain_id)
  return [a.name.strip() for a in parent.parent().atoms()
          if a.element_is_hydrogen() and
          (fixtures.site(a) - fixtures.site(parent)).length() < 1.3]

def test_shared_server():
  s1 = capping.shared_mon_lib_srv()
  assert s1 is capping.shared_mon_lib_srv()

def test_cap_occupancy_zero():
  # occupancy 0 only for caps whose name is no dictionary atom of the
  # residue (HCG); the dictionary-named HB2, HB3 take CB's occupancy, as any
  # reduce2 H. CB at 0.7 so that the two cases differ.
  cb_line = fixtures.pdb_arg.splitlines()[5]
  assert ' CB  ARG A   3' in cb_line and cb_line[54:60] == '  1.00', cb_line
  pdb_str = fixtures.pdb_arg.replace(cb_line,
    cb_line[:54] + '  0.70' + cb_line[60:])
  model, result = fixtures.cap(pdb_str)
  cb = fixtures.atom_named(model, 3, 'CB')
  assert cb.occ == 0.7, cb.occ
  names = h_names_on(model, 3, 'CB')
  assert names == ['HB2', 'HB3', 'HCG'], names
  for name in names:
    h = fixtures.atom_named(model, 3, name)
    expected = 0.0 if name == 'HCG' else cb.occ
    assert h.occ == expected, (name, h.occ)
    assert (h.b, h.segid, h.hetero) == (cb.b, cb.segid, cb.hetero)

def test_cap_name_avoids_dictionary():
  # DC 4 of pdb_dna without O3': the cap on C3' must not be HO3' (the
  # dictionary H of the absent O3'); the fallback is HC1
  model, result = fixtures.cap(fixtures.without(fixtures.pdb_dna, 4, ["O3'"]))
  names = h_names_on(model, 4, "C3'", 'C')
  assert "HO3'" not in names, names
  assert 'HC1' in names, names

def test_cap_name_keeps_free_dictionary_h():
  model, result = fixtures.cap(fixtures.pdb_arg)
  assert h_names_on(model, 3, 'CB') == ['HB2', 'HB3', 'HCG']

def test_budget_refactor_unchanged():
  # find_dangling_atoms before and after the refactor: same records (the
  # values are those of the unrefactored code; CD has no H in the fixture)
  r = fixtures.find(fixtures.pdb_lys)
  by = fixtures.records_by_name(r)
  assert list(by) == [('LYS', '2', '', 'CD')], list(by)
  assert (by[('LYS', '2', '', 'CD')].n_fill,
          by[('LYS', '2', '', 'CD')].missing_names) == (3, ['CE'])

def capped_hierarchy(pdb_str):
  model, result = fixtures.cap(pdb_str)
  return model, result, model.get_hierarchy()

def added_i_seqs(hierarchy, result):
  """
  i_seqs of the caps cap_dangling_atoms placed: result.added is
  [(atom_key of the parent, [cap names])]; the caps are the atoms of the
  parent's atom group with those names.
  """
  want = set()
  for key, names in result.added:
    model_id, chain_id, resseq, icode, resname, altloc, parent_name = key
    for a in hierarchy.atoms():
      ag = a.parent(); rg = ag.parent()
      if ((rg.parent().id.strip(), rg.resseq.strip(), ag.altloc.strip())
          == (chain_id, resseq, altloc) and a.name.strip() in names):
        want.add(a.i_seq)
  return want

def test_find_matches_added():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  added = added_i_seqs(h, result)
  found = set(r.selection.iselection())
  assert found and found <= added, (found, added)
  assert not r.rejected, r.rejected
  # HB2 and HB3 on the truncated CB are dictionary atoms: not caps
  names = set(h.atoms()[i].name.strip() for i in found)
  assert names == {'HCG'}, names
  assert [(x.i_seq, h.atoms()[x.parent_i_seq].name.strip(), x.coordination,
           x.steric) for x in r.records] == [(list(found)[0], 'CB', 4, 4)], \
    r.records

def test_find_lys_cap():
  # CD without CE: HD2 and HD3 are dictionary H, the cap is HCE
  model, result, h = capped_hierarchy(fixtures.pdb_lys)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  rec = [x for x in r.records if h.atoms()[x.parent_i_seq].name.strip() == 'CD']
  assert [h.atoms()[x.i_seq].name.strip() for x in rec] == ['HCE'], rec
  assert all(x.coordination == 4 and x.steric == 4 for x in rec), rec

def test_find_old_ho3_name():
  # a file written before Task 2: the cap on C3' named HO3' (dictionary H of
  # the absent O3') is recognized by condition 2
  model, result, h = capped_hierarchy(
    fixtures.without(fixtures.pdb_dna, 4, ["O3'"]))
  cap = [a for a in h.atoms() if a.name.strip() == 'HC1'][0]
  cap.set_name(new_name="HO3'")
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection[cap.i_seq]
  assert not r.rejected, r.rejected

def test_find_rejects_wrong_name_complete_residue():
  model = fixtures.model_with_h(fixtures.pdb_ala_ala)
  h = model.get_hierarchy()
  a = fixtures.atom_named(model, 2, 'HB1'); a.set_name(new_name=' HX1')
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert len(r.rejected) == 1, r.rejected
  assert 'misses no heavy atom' in r.rejected[0][1], r.rejected

def test_find_rejects_excess():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cb = fixtures.atom_named(model, 3, 'CB')
  extra = fixtures.atom_named(model, 3, 'HCG').detached_copy()
  extra.set_name(new_name=' HX9')
  extra.set_xyz(new_xyz=(cb.xyz[0]+1.0, cb.xyz[1], cb.xyz[2]))
  cb.parent().append_atom(extra); h.atoms().reset_i_seq()
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert len(r.rejected) == 2, r.rejected
  assert all('valence allows' in reason for label, reason in r.rejected), \
    r.rejected

def test_find_rejects_no_parent():
  # condition 3: an H 3 A from every heavy atom of its residue
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cap = fixtures.atom_named(model, 3, 'HCG')
  cap.set_xyz(new_xyz=(cap.xyz[0]+3.0, cap.xyz[1], cap.xyz[2]))
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert [reason for label, reason in r.rejected] == [
    'no heavy atom of its residue within bonding distance'], r.rejected

pdb_lone_cb_hx1 = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  CB  ALA A   1      25.000  25.000  25.000  1.00  7.41           C
ATOM      2  HX1 ALA A   1      26.000  25.000  25.000  1.00  7.41           H
END
"""

def test_find_rejects_lone_parent():
  # an unknown-named H next to a CB with no other atom: the parent has no
  # bonded neighbour to orient a cap against (capping rule 5), so no cap,
  # with or without a bond table; processing does not assert
  pdb_inp = iotbx.pdb.input(lines=pdb_lone_cb_hx1.split('\n'),
                            source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  r = capping.find_capping_h(model.get_hierarchy(), model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert r.rejected == [('ALA A 1 HX1', 'parent CB has no bonded neighbour')], \
    r.rejected
  log = StringIO()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=pdb_lone_cb_hx1.split('\n'), source_info=None),
    log=log, stop_for_unknowns=False)
  m.process(make_restraints=True)
  bonds, angles = cap_proxies(m)
  assert not bonds and not angles
  assert 'Capping H: ALA A 1 HX1 is not a capping H: parent CB has no ' \
    'bonded neighbour' in log.getvalue(), log.getvalue()

# ALA with CB and O only: O is 1.4 A from CB, within bonding distance, but
# not a dictionary partner of CB (the NAD_594_HD case: an invented ligand
# whose name collides with an unrelated CCD entry)
pdb_cb_near_o = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  CB  ALA A   1      25.000  25.000  25.000  0.75  7.41           C
ATOM      2  O   ALA A   1      26.400  25.000  25.000  1.00  7.41           O
ATOM      3  HX1 ALA A   1      24.000  25.000  25.000  0.75  7.41           H
END
"""

def test_distance_neighbours_are_dictionary_partners():
  # without a bond table a same-residue atom counts as a neighbour only if
  # it is a dictionary partner of the parent, as pdb_interpretation bonds
  # nothing else inside a residue: same answer with and without the table
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=pdb_cb_near_o.split('\n'), source_info=None),
    log=null_out(), stop_for_unknowns=False)
  m.process(make_restraints=True)
  h, srv = m.get_hierarchy(), m.get_mon_lib_srv()
  table = capping.bond_table(m)
  assert not table.get(0), table
  with_table = capping.find_capping_h(h, srv, bond_table=table)
  without = capping.find_capping_h(h, srv)
  assert list(without.selection) == list(with_table.selection)
  assert without.selection.count(True) == 0
  assert without.rejected == with_table.rejected == [
    ('ALA A 1 HX1', 'parent CB has no bonded neighbour')], (
    without.rejected, with_table.rejected)
  assert m.selection('capping_h').count(True) == 0

def test_two_models():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  atom_lines = [l for l in model.model_as_pdb().splitlines()
                if l.startswith(('ATOM', 'HETATM'))]
  two = iotbx.pdb.input(lines=['MODEL        1'] + atom_lines +
    ['ENDMDL', 'MODEL        2'] + atom_lines + ['ENDMDL', 'END'],
    source_info=None).construct_hierarchy()
  assert len(two.models()) == 2
  r = capping.find_capping_h(two, model.get_mon_lib_srv())
  assert r.selection.count(True) == 2 and not r.rejected, r.rejected

def test_altloc_parent():
  # CB B without OG: HB2 and HB3 are dictionary H, the cap is HOG (altloc B)
  model, result, h = capped_hierarchy(fixtures.pdb_ser_altloc)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  found = set(r.selection.iselection())
  assert found and found <= added_i_seqs(h, result), found
  assert [(h.atoms()[i].name.strip(), h.atoms()[i].parent().altloc)
          for i in found] == [('HOG', 'B')], found
  assert not r.rejected, r.rejected

def test_deuterium_and_old_names_not_candidates():
  # D for H (DB2) and old digit-first names (1HB) are dictionary names to
  # pdb_interpretation, as are the chain-end H (H1, H3 of NH3)
  model = fixtures.model_with_h(fixtures.pdb_ala_ala)
  h = model.get_hierarchy()
  fixtures.atom_named(model, 2, 'HB1').set_name(new_name='1HB ')
  d = fixtures.atom_named(model, 2, 'HB2')
  d.set_name(new_name=' DB2'); d.element = ' D'
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0 and not r.rejected, r.rejected

def test_stale_i_seq():
  # atom.i_seq cannot be set from Python: the copied atoms of a selection
  # have i_seq 0
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  keep = flex.bool(h.atoms_size(), True)
  keep[fixtures.atom_named(model, 4, 'HA').i_seq] = False
  stale = h.select(keep, copy_atoms=True)
  i_seqs = list(stale.atoms().extract_i_seq())
  assert i_seqs != list(range(len(i_seqs))), i_seqs
  fresh = stale.deep_copy(); fresh.atoms().reset_i_seq()
  srv = model.get_mon_lib_srv()
  r = capping.find_capping_h(stale, srv)
  assert list(r.selection) == list(capping.find_capping_h(fresh, srv).selection)
  assert r.selection.count(True) == 1
  # the caller's hierarchy keeps its i_seqs
  assert list(stale.atoms().extract_i_seq()) == i_seqs

def test_selection_keyword_helper():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  assert (list(capping.capping_h_selection(h)) ==
          list(capping.find_capping_h(h, model.get_mon_lib_srv()).selection))

def big_hierarchy_without_caps(n_atoms=20000):
  """
  pdb_lys_full with H, copied into separate chains until n_atoms.
  """
  base = fixtures.model_with_h(fixtures.pdb_lys_full).get_hierarchy()
  chain = base.only_model().only_chain()
  n_copies = -(-n_atoms // chain.atoms_size())
  root = iotbx.pdb.hierarchy.root()
  m = iotbx.pdb.hierarchy.model()
  root.append_model(m)
  for k in range(n_copies):
    c = chain.detached_copy()
    c.id = '%d' % k
    m.append_chain(c)
  root.atoms().reset_i_seq()
  assert root.atoms_size() >= n_atoms
  return root

def test_speed():
  import time
  h = big_hierarchy_without_caps()
  srv = capping.shared_mon_lib_srv()
  t0 = time.time(); r = capping.find_capping_h(h, srv); dt = time.time() - t0
  assert r.selection.count(True) == 0 and not r.rejected, r.rejected
  assert dt < 2.0, dt

def test_speed_no_h():
  # a model without H or D: one flex element test, no per-atom Python work
  import time
  h = big_hierarchy_without_caps(n_atoms=50000)
  h = h.select(~h.atom_selection_cache().selection('element H or element D'))
  assert h.atoms_size() >= 20000, h.atoms_size()
  h.atoms().reset_i_seq()
  srv = capping.shared_mon_lib_srv()
  t0 = time.time(); r = capping.capping_h_selection(h, srv)
  dt = time.time() - t0
  assert r.count(True) == 0
  assert dt < 0.1, dt
  t0 = time.time()
  sel = h.atom_selection_cache().selection('capping_h')
  dt = time.time() - t0
  assert sel.count(True) == 0
  assert dt < 0.1, dt

def test_dictionary_h_aliases():
  # the name aliases pdb_interpretation applies (_get_mappings)
  srv = capping.shared_mon_lib_srv()
  for resname, name, expected in [
      ('DG', 'H1*',  ("H1'", "C1'")),  # prime/star
      ('DG', "h1'",  ("H1'", "C1'")),  # case folding
      ('DG', 'D1*',  ("H1'", "C1'")),  # D, star
      ('DG', '1H5*', ("H5'", "C5'")),  # backbone alias
      ('DG', "H2'1", ("H2'", "C2'")),  # backbone alias
      ('DG', 'H5T',  ("HO5'", "O5'")), # backbone alias of a chain-end H
      ('G',  '2HO*', ("HO2'", "O2'")), # backbone alias
      ('ALA', 'hb1', ('HB1', 'CB')),   # case folding
      ('ALA', '1HB', ('HB1', 'CB')),   # digit rotation
      ('ALA', 'H1*', None),            # no backbone aliases for a peptide
      ]:
    names = capping._dictionary_names(resname, srv, {})
    synonyms = srv.comp_synonym_atom_list_dict.get(resname, {})
    found = capping._dictionary_h(name, names, synonyms)
    assert found == expected, (resname, name, found)

def star_model(log):
  """
  (unprocessed model, [real cap i_seq], H1* i_seq): DG 3 of pdb_dna without
  its base, capped (HN9 on C1'), with H1' renamed H1*; pdb_interpretation
  attaches H1* through its prime/star aliases. log: the model's log.
  """
  base = ["N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"]
  model, result = fixtures.cap(fixtures.without(fixtures.pdb_dna, 3, base))
  assert [names for key, names in result.added] == [['HN9']], result.added
  fixtures.atom_named(model, 3, "H1'", 'C').set_name(new_name=" H1*")
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=log)
  atoms = m.get_hierarchy().atoms()
  caps = [a.i_seq for a in atoms if a.name.strip() == 'HN9']
  star = [a.i_seq for a in atoms if a.name.strip() == 'H1*']
  assert len(caps) == 1 and len(star) == 1
  return m, caps, star[0]

def check_star_processed(m, log, caps, star):
  # exactly the real cap: one capping bond (to HN9), capping angles only of
  # HN9 (H1*, an attached neighbour of C1', is their outer atom)
  assert 'Capping H: 1 recognized' in log.getvalue(), log.getvalue()
  bonds, angles = cap_proxies(m)
  assert [p.i_seqs[1] for p in bonds] == caps, [p.i_seqs for p in bonds]
  assert angles and all(p.i_seqs[2] in caps for p in angles), \
    [p.i_seqs for p in angles]

def test_star_alias_not_a_cap():
  # I1: an H attached by pdb_interpretation through an alias is no cap,
  # neither for the keyword (unprocessed or processed) nor for the hook
  log = StringIO()
  m, caps, star = star_model(log)
  assert list(m.selection('capping_h').iselection()) == caps
  asc = m.get_hierarchy().atom_selection_cache()
  assert list(asc.selection('capping_h').iselection()) == caps
  m.process(make_restraints=True)
  check_star_processed(m, log, caps, star)
  assert list(m.selection('capping_h').iselection()) == caps

def test_hook_acts_on_unattached_only():
  # I1 (a): even when _dictionary_h misses the alias (simulated), the hook
  # does not make the attached H1* a cap: pdb_interpretation attached it
  log = StringIO()
  m, caps, star = star_model(log)
  saved = capping._dictionary_h
  def no_star(name, names, synonyms):
    if name == 'H1*': return None
    return saved(name, names, synonyms)
  capping._dictionary_h = no_star
  try:
    # the finder alone is fooled
    assert m.get_hierarchy().atom_selection_cache().selection(
      'capping_h')[star]
    m.process(make_restraints=True)
    check_star_processed(m, log, caps, star)
  finally:
    capping._dictionary_h = saved

# ------------------------------------------------------------------------------
# pdb_interpretation restrains the caps

def processed(pdb_str, **kw):
  """
  (model, cap_dangling_atoms result): pdb_str with H placed and capped,
  written and read back, processed with make_restraints=True.
  """
  model, result, h = capped_hierarchy(pdb_str)
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out(), **kw)
  m.process(make_restraints=True)
  return m, result

def cap_proxies(m):
  from cctbx.geometry_restraints.linking_class import linking_class
  oid = linking_class().get_origin_id('capping')
  grm = m.get_restraints_manager().geometry
  bonds = [p for p in grm.get_all_bond_proxies(sites_cart=m.get_sites_cart())[0]
           if p.origin_id == oid]
  angles = [p for p in grm.get_all_angle_proxies() if p.origin_id == oid]
  return bonds, angles

def angle_names(m, angles):
  atoms = m.get_hierarchy().atoms()
  return sorted(tuple(atoms[i].name.strip() for i in p.i_seqs) for p in angles)

def test_find_with_bond_table_equals_without():
  m, result = processed(fixtures.pdb_arg)
  h = m.get_hierarchy()
  a = capping.find_capping_h(h, m.get_mon_lib_srv())
  b = capping.find_capping_h(h, m.get_mon_lib_srv(),
    bond_table=capping.bond_table(m))
  assert a.selection.count(True) == 1
  assert list(a.selection) == list(b.selection)
  assert ([(x.i_seq, x.parent_i_seq, x.coordination) for x in a.records] ==
          [(x.i_seq, x.parent_i_seq, x.coordination) for x in b.records])

def test_process_capped_model():
  m, result = processed(fixtures.pdb_arg)      # default stop_for_unknowns
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1
  assert angle_names(m, angles) == [('CA', 'CB', 'HCG'), ('HB2', 'CB', 'HCG'),
    ('HB3', 'CB', 'HCG')], angle_names(m, angles)
  assert all(abs(a.angle_ideal - 109.5) < 1e-6 for a in angles)
  assert m.selection('capping_h').count(True) == 1
  # default pdb_interpretation: X-ray distances
  assert abs(bonds[0].distance_ideal - capping.xh_length('C', 4, False)) < 1e-6
  assert abs(bonds[0].weight - 1/0.02**2) < 1e-6

def test_process_lys_angles():
  # one cap, HCE on CD; CD also has CG and the dictionary HD2, HD3: 3 angles
  m, result = processed(fixtures.pdb_lys)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1
  assert angle_names(m, angles) == [('CG', 'CD', 'HCE'), ('HD2', 'CD', 'HCE'),
    ('HD3', 'CD', 'HCE')], angle_names(m, angles)
  assert all(abs(a.angle_ideal - 109.5) < 1e-6 for a in angles)

def test_process_two_caps_on_one_parent():
  # ALA A 1 C without O: its caps get one angle with each other, once
  m, result = processed(fixtures.pdb_ala_ala)
  bonds, angles = cap_proxies(m)
  atoms = m.get_hierarchy().atoms()
  caps = sorted(atoms[p.i_seqs[1]].name.strip() for p in bonds)
  assert len(caps) == 2, caps
  names = angle_names(m, angles)
  between = [x for x in names if x[0] in caps and x[2] in caps]
  assert len(between) == 1, names
  # C of ALA 1: CA, N of ALA 2 and the two caps, steric 4
  assert len(names) == 5, names
  assert all(abs(a.angle_ideal - 109.5) < 1e-6 for a in angles)

def test_statistics_exclude_caps():
  # the capped model and the same model with the caps removed, reprocessed:
  # same covalent bond and angle statistics
  m, result = processed(fixtures.pdb_arg)
  sel = capping.capping_h_selection(m.get_hierarchy())
  assert sel.count(True) == 1
  text = m.select(~sel).model_as_pdb()
  uncapped = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=text.split('\n'), source_info=None), log=null_out())
  uncapped.process(make_restraints=True)
  s_c = m.geometry_statistics()
  s_u = uncapped.geometry_statistics()
  for stat in ('bond', 'angle'):
    c, u = getattr(s_c, stat)(), getattr(s_u, stat)()
    assert c.n == u.n, (stat, c.n, u.n)
    assert abs(c.mean - u.mean) < 1e-6, (stat, c.mean, u.mean)
  geo = m.restraints_as_geo()
  assert 'Capping H' in geo

def test_mmcif_no_link():
  # a cap bond is no struct_conn and no cif link
  m, result = processed(fixtures.pdb_arg)
  text = m.model_as_mmcif()
  assert 'HCG' in text
  assert '_struct_conn' not in text and 'link_capping' not in text

def test_idempotent():
  m, result = processed(fixtures.pdb_arg)
  r = capping.find_dangling_atoms(m)
  assert not r.records, r.records

def test_still_stops_wrong_name():
  # an H pdb_interpretation cannot attach that is no cap still stops the run,
  # and the log says why it is no cap
  model = fixtures.model_with_h(fixtures.pdb_ala_ala)
  fixtures.atom_named(model, 2, 'HB1').set_name(new_name=' HX1')
  text = model.model_as_pdb()
  log = StringIO()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=log)
  fixtures.assert_sorry(lambda: m.process(make_restraints=True), 'nknown')
  assert 'Capping H: ALA A 2 HX1 is not a capping H: parent CB misses no ' \
    'heavy atom' in log.getvalue(), log.getvalue()

def test_still_stops_excess():
  # spec 11.7: one H too many on a truncated residue (test_find_rejects_excess)
  # keeps the unknown-atom Sorry of process()
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cb = fixtures.atom_named(model, 3, 'CB')
  extra = fixtures.atom_named(model, 3, 'HCG').detached_copy()
  extra.set_name(new_name=' HX9')
  extra.set_xyz(new_xyz=(cb.xyz[0]+1.0, cb.xyz[1], cb.xyz[2]))
  cb.parent().append_atom(extra); h.atoms().reset_i_seq()
  text = h.as_pdb_string()
  log = StringIO()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=log)
  fixtures.assert_sorry(lambda: m.process(make_restraints=True), 'nknown')
  assert 'valence allows' in log.getvalue(), log.getvalue()

def test_process_no_stop():
  m, result = processed(fixtures.pdb_arg, stop_for_unknowns=False)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1

def test_no_caps_no_work():
  # all H have bond proxies and all atoms known types: find_capping_h not run
  model = fixtures.model_with_h(fixtures.pdb_arg)
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out())
  def fail(*args, **kw): raise AssertionError('find_capping_h called')
  saved, capping.find_capping_h = capping.find_capping_h, fail
  try: m.process(make_restraints=True)
  finally: capping.find_capping_h = saved
  bonds, angles = cap_proxies(m)
  assert not bonds and not angles

def test_neutron_and_deuterium():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cap = fixtures.atom_named(model, 3, 'HCG')
  cap.set_element('D'); cap.set_name(new_name=' DCG')
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out())
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = True
  m.process(make_restraints=True, pdb_interpretation_params=p)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1 and len(angles) == 3
  assert abs(bonds[0].distance_ideal - capping.xh_length('C', 4, True)) < 1e-6

def test_riding_minimization():
  """
  The riding manager does not parameterize the H of LYS CD without CE: the
  dictionary HD2, HD3 and the cap HCE alike. phenix.geometry_minimization
  (mmtbx.refinement.geometry_minimization.run2) refines an H without a
  parameterization as a free atom whether riding_h is on or off, so the cap
  behaves as the dictionary H and all three end at the X-H length (default:
  X-ray).
  """
  from mmtbx.command_line import geometry_minimization
  m, result = processed(fixtures.pdb_lys)
  m.setup_riding_h_manager()
  m.idealize_h_riding()
  with open('capped_lys.pdb', 'w') as f: f.write(m.model_as_pdb())
  for riding in (True, False):
    prefix = 'capped_lys_min_riding_%s' % riding
    geometry_minimization.run(args=['capped_lys.pdb',
      'output_file_name_prefix=%s' % prefix,
      'minimization.riding_h=%s' % riding], log=null_out())
    out = iotbx.pdb.input(file_name=prefix+'.pdb').construct_hierarchy()
    by_name = dict((a.name.strip(), a) for a in out.atoms())
    cd = matrix.col(by_name['CD'].xyz)
    d = dict((n, (matrix.col(by_name[n].xyz) - cd).length())
             for n in ('HD2', 'HCE'))
    assert abs(d['HCE'] - d['HD2']) < 0.01, (riding, d)
    assert abs(d['HCE'] - capping.xh_length('C', 4, False)) < 0.01, (riding, d)

def test_keyword_selection():
  m, result = processed(fixtures.pdb_arg)
  expected = capping.capping_h_selection(m.get_hierarchy())
  assert list(m.selection('capping_h')) == list(expected)
  asc = m.get_hierarchy().atom_selection_cache()
  assert list(asc.selection('capping_h')) == list(expected)
  assert asc.selection('not capping_h').count(False) == 1
  assert asc.selection('capping_h and resseq 3').count(True) == 1

def test_keyword_no_h():
  h = iotbx.pdb.input(lines=fixtures.pdb_arg.split('\n'),
    source_info=None).construct_hierarchy()
  assert h.atom_selection_cache().selection('capping_h').count(True) == 0

def test_pdbtools_remove():
  from libtbx import easy_run
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  uncapped = fixtures.model_with_h(fixtures.pdb_arg)
  with open('capped_arg.pdb', 'w') as fh:
    fh.write(model.model_as_pdb())
  assert easy_run.call('phenix.pdbtools capped_arg.pdb remove=capping_h '
                       'output.file_name=stripped.pdb > pdbtools.log') == 0
  out = iotbx.pdb.input(file_name='stripped.pdb').construct_hierarchy()
  names = [(a.parent().parent().resseq.strip(), a.name.strip())
    for a in out.atoms()]
  want = [(a.parent().parent().resseq.strip(), a.name.strip())
    for a in uncapped.get_hierarchy().atoms()]
  # capping also completes the dictionary H of the cut atom (HB2, HB3 of CB);
  # those are ordinary H, only HCG is a cap
  i = want.index(('3', 'HA')) + 1
  want[i:i] = [('3', 'HB2'), ('3', 'HB3')]
  assert sorted(names) == sorted(want), (names, want)
  assert easy_run.call('phenix.pdbtools capped_arg.pdb remove="not capping_h" '
                       'output.file_name=only_caps.pdb > pdbtools2.log') == 0
  only = iotbx.pdb.input(file_name='only_caps.pdb').construct_hierarchy()
  assert [a.name.strip() for a in only.atoms()] == ['HCG']

def test_occupancy_survives():
  # pdb_ser_altloc caps CB B: the cap HOG is in the altloc B atom group, so
  # occupancy_groups_simple would put it into conformer B's constrained group
  from mmtbx.refinement import occupancies
  m, result = processed(fixtures.pdb_ser_altloc)
  caps = m.selection('capping_h')
  assert caps.count(True) > 0
  assert [m.get_hierarchy().atoms()[i].parent().altloc
          for i in caps.iselection()] == ['B']
  m.reset_occupancy_for_hydrogens_simple()
  assert all(o == 0.0 for o in m.get_occ().select(caps))
  # reset_occupancies_for_hydrogens reads refinement_flags.s_occupancies;
  # give it the hierarchy's own groups, in which the cap is in conformer B's
  sog = m.get_hierarchy().occupancy_groups_simple(
    common_residue_name_class_only=None, always_group_adjacent=False,
    ignore_hydrogens=False)
  assert set(caps.iselection()) <= set(i for g in sog for gg in g for i in gg)
  m.refinement_flags = mmtbx.refinement.refinement_flags.manager(
    occupancies=True,
    s_occupancies=[[flex.size_t(gg) for gg in g] for g in sog])
  m.reset_occupancies_for_hydrogens()
  assert all(o == 0.0 for o in m.get_occ().select(caps))
  cap_set = set(caps.iselection())
  def assert_no_caps(groups):
    flat = set(i for g in groups for gg in g for i in gg)
    assert flat and not (flat & cap_set), flat & cap_set
  assert_no_caps(occupancies.occupancy_selections(model=m))
  # user selections covering the capped residue do not put the cap back
  assert cap_set <= set(m.selection('resseq 4').iselection())
  assert_no_caps(occupancies.occupancy_selections(model=m,
    other_individual_selection_strings=['resseq 4']))
  assert_no_caps(occupancies.occupancy_selections(model=m,
    other_constrained_groups=[group_args(
      selection=['resseq 4 and altloc A', 'resseq 4 and altloc B'])]))
  # structure factors: caps add nothing. The reference is this model without
  # the caps, not fixtures.model_with_h(pdb_ser_altloc): that one has
  # unrounded coordinates (processed() goes through PDB text) and lacks HB2,
  # HB3 of CB B, the dictionary H capping completes, which the resets above
  # rightly gave CB B's occupancy
  uncapped = m.select(~caps)
  for mm in (m, uncapped):
    mm.setup_scattering_dictionaries(scattering_table='n_gaussian')
  fc_c = m.get_xray_structure().structure_factors(d_min=2.0).f_calc()
  fc_u = uncapped.get_xray_structure().structure_factors(d_min=2.0).f_calc()
  assert flex.max(flex.abs(fc_c.data() - fc_u.common_set(fc_c).data())) < 1e-6

# SER B 64 of tst_add_hydrogen_12 SER_1S6_B_401_4jxg_H without N, ligand
# and H: a residue missing its N
pdb_ser_no_n = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  CA  SER B  64       8.812   8.398  16.900  1.00 15.06           C
ATOM      2  C   SER B  64       7.911   7.427  16.145  1.00 14.26           C
ATOM      3  O   SER B  64       6.732   7.735  16.004  1.00 14.43           O
ATOM      4  CB  SER B  64       8.998   9.607  16.069  1.00 14.19           C
ATOM      5  OG  SER B  64       9.585  10.629  16.889  1.00 18.69           O
END
"""

def test_placement_does_not_cap_placeholders():
  """
  reduce2 placement adds every dictionary H of SER B 64 without N and removes
  the backbone H, whose parent is absent, as isolated; with cap recognition on
  during placement the hook bonded it to CA (CA misses N) and it survived
  with occupancy 0 (tst_add_hydrogen_12 tst_04). The finished model, capped
  and processed with the default parameters, has its real cap on CA
  restrained.
  """
  model = fixtures.model_with_h(pdb_ser_no_n)
  names = sorted(a.name.strip() for a in model.get_hierarchy().atoms())
  assert names == ['C', 'CA', 'CB', 'HA', 'HB2', 'HB3', 'HG', 'O', 'OG'], \
    names
  assert h_names_on(model, 64, 'CA', chain_id='B') == ['HA']
  m, result = processed(pdb_ser_no_n)         # default recognize_capping_h
  bonds, angles = cap_proxies(m)
  atoms = m.get_hierarchy().atoms()
  pairs = [tuple(atoms[i].name.strip() for i in p.i_seqs) for p in bonds]
  # HN for the absent N, HOXT on C (no OXT)
  assert sorted(pairs) == [('C', 'HOXT'), ('CA', 'HN')], pairs
  assert sorted(m.selection('capping_h').iselection()) == \
    sorted(p.i_seqs[1] for p in bonds)

def test_recognition_off_leaves_caps_unknown():
  # recognize_capping_h=False: a capped model's cap is an unknown atom again
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out())
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  assert p.pdb_interpretation.recognize_capping_h is True
  p.pdb_interpretation.recognize_capping_h = False
  fixtures.assert_sorry(lambda: m.process(make_restraints=True,
    pdb_interpretation_params=p), 'nknown')

def run():
  test_origin_id()
  test_shared_server()
  test_cap_occupancy_zero()
  test_cap_name_avoids_dictionary()
  test_cap_name_keeps_free_dictionary_h()
  test_budget_refactor_unchanged()
  test_find_matches_added()
  test_find_lys_cap()
  test_find_old_ho3_name()
  test_find_rejects_wrong_name_complete_residue()
  test_find_rejects_excess()
  test_find_rejects_no_parent()
  test_find_rejects_lone_parent()
  test_distance_neighbours_are_dictionary_partners()
  test_two_models()
  test_altloc_parent()
  test_deuterium_and_old_names_not_candidates()
  test_stale_i_seq()
  test_selection_keyword_helper()
  test_speed()
  test_speed_no_h()
  test_dictionary_h_aliases()
  test_star_alias_not_a_cap()
  test_hook_acts_on_unattached_only()
  test_find_with_bond_table_equals_without()
  test_process_capped_model()
  test_process_lys_angles()
  test_process_two_caps_on_one_parent()
  test_statistics_exclude_caps()
  test_mmcif_no_link()
  test_idempotent()
  test_still_stops_wrong_name()
  test_still_stops_excess()
  test_process_no_stop()
  test_no_caps_no_work()
  test_neutron_and_deuterium()
  test_riding_minimization()
  test_keyword_selection()
  test_keyword_no_h()
  test_pdbtools_remove()
  test_occupancy_survives()
  test_placement_does_not_cap_placeholders()
  test_recognition_off_leaves_caps_unknown()

if __name__ == '__main__':
  run()
  print('OK')
