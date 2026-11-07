"""reduce2: IUPAC names for the planar NH2 H of Asn, Gln and Arg (A49).
IUPAC 1998 (PAC 70:117, Fig. 1): H "1" is cis, "2" trans to the heavy atom
three bonds closer to the main chain (CB for Asn, CG for Gln, NE for Arg).
Placement named them by the bogus start site the riding torsion (periodicity
2) snaps to. Supplied H/D keep their names, sites and values; so do placed H
whose restraints fix the side (a user CIF).
Fragment: 4wsc Arg A 430 - Asn A 433, no H.
"""
from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
import iotbx.pdb
from iotbx.data_manager import DataManager
from scitbx import matrix
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

pdb_4wsc = """\
CRYST1  135.640  259.390  143.310  90.00 105.12  90.00 P 1 21 1
ATOM   3142  N   ARG A 430      78.087  18.464  16.235  1.00 34.35           N
ATOM   3143  CA  ARG A 430      77.510  19.646  16.862  1.00 22.53           C
ATOM   3144  C   ARG A 430      77.128  20.698  15.847  1.00 32.56           C
ATOM   3145  O   ARG A 430      77.594  20.673  14.708  1.00 22.11           O
ATOM   3146  CB  ARG A 430      78.451  20.249  17.889  1.00 23.78           C
ATOM   3147  CG  ARG A 430      78.708  19.354  19.056  1.00 24.25           C
ATOM   3148  CD  ARG A 430      77.538  19.321  20.018  1.00 23.98           C
ATOM   3149  NE  ARG A 430      77.784  18.353  21.079  1.00 24.50           N
ATOM   3150  CZ  ARG A 430      77.618  17.035  20.952  1.00 24.02           C
ATOM   3151  NH1 ARG A 430      77.163  16.505  19.818  1.00 22.98           N
ATOM   3152  NH2 ARG A 430      77.898  16.242  21.974  1.00 24.68           N
ATOM   3153  N   GLY A 431      76.253  21.605  16.287  1.00 35.90           N
ATOM   3154  CA  GLY A 431      75.797  22.758  15.519  1.00 43.33           C
ATOM   3155  C   GLY A 431      76.312  24.102  16.037  1.00 44.95           C
ATOM   3156  O   GLY A 431      77.322  24.159  16.760  1.00 24.06           O
ATOM   3157  N   GLN A 432      75.625  25.184  15.647  1.00 50.88           N
ATOM   3158  CA  GLN A 432      76.005  26.558  16.018  1.00 36.67           C
ATOM   3159  C   GLN A 432      75.494  27.003  17.391  1.00 24.53           C
ATOM   3160  O   GLN A 432      75.995  27.956  17.981  1.00 25.57           O
ATOM   3161  CB  GLN A 432      75.540  27.524  14.924  1.00 24.00           C
ATOM   3162  CG  GLN A 432      76.467  27.522  13.708  1.00 48.40           C
ATOM   3163  CD  GLN A 432      75.864  28.157  12.460  1.00 57.17           C
ATOM   3164  OE1 GLN A 432      74.642  28.214  12.298  1.00 54.17           O
ATOM   3165  NE2 GLN A 432      76.729  28.650  11.574  1.00 62.36           N
ATOM   3166  N   ASN A 433      74.542  26.258  17.924  1.00 33.78           N
ATOM   3167  CA  ASN A 433      73.985  26.546  19.234  1.00 33.90           C
ATOM   3168  C   ASN A 433      73.245  25.324  19.741  1.00 23.49           C
ATOM   3169  O   ASN A 433      73.183  24.312  19.046  1.00 22.71           O
ATOM   3170  CB  ASN A 433      73.083  27.781  19.159  1.00 35.14           C
ATOM   3171  CG  ASN A 433      71.961  27.623  18.157  1.00 23.25           C
ATOM   3172  OD1 ASN A 433      71.002  26.904  18.395  1.00 22.56           O
ATOM   3173  ND2 ASN A 433      72.094  28.279  17.017  1.00 23.24           N
END
"""

# (resname, parent, heavy neighbour, reference, H "1", H "2"), H names
# without the element letter
GROUPS = [('ARG', 'NH1', 'CZ', 'NE', 'H11', 'H12'),
          ('ARG', 'NH2', 'CZ', 'NE', 'H21', 'H22'),
          ('GLN', 'NE2', 'CD', 'CG', 'E21', 'E22'),
          ('ASN', 'ND2', 'CG', 'CB', 'D21', 'D22')]
FLIPS = ['add_flip_movers=True']
NEUTRON = ['neutron_option=hd_and_h', 'use_neutron_distances=True']

def run(prefix, pdb_str, args=(), log=None):
  '''reduce2 on pdb_str. :return: output hierarchy.'''
  fn, out = prefix + '_in.pdb', prefix + '_out.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  if os.path.exists(out): os.remove(out)
  if log is None: log = StringIO()
  run_program(program_class=reduce2.Program, logger=log, args=[fn,
    'output.filename=' + out, 'output.description_file_name=' + prefix +
    '.txt', 'output.overwrite=True'] + list(args))
  return iotbx.pdb.input(out).construct_hierarchy()

def residues(hierarchy):
  '''{(resname, altloc): {name: atom}}, blank atoms in every altloc.'''
  result = {}
  for rg in hierarchy.residue_groups():
    blank = [a for ag in rg.atom_groups() if ag.altloc == '' for a in ag.atoms()]
    for ag in rg.atom_groups():
      d = result.setdefault((ag.resname, ag.altloc), {})
      for a in blank + list(ag.atoms()):
        d[a.name.strip()] = a
  return result

def sides(hierarchy):
  '''{(resname, altloc, H name): "cis" or "trans" to the reference}.'''
  result = {}
  for (resname, altloc), at in residues(hierarchy).items():
    for g in GROUPS:
      if g[0] != resname: continue
      for e in ('H', 'D'):
        for s in g[4:]:
          if e + s not in at: continue
          t = matrix.dihedral_angle(sites=[at[g[3]].xyz, at[g[2]].xyz,
            at[g[1]].xyz, at[e + s].xyz], deg=True)
          result[(resname, altloc, e + s)] = 'cis' if abs(t) < 90 else 'trans'
  return result

def assert_iupac(hierarchy, n_expected):
  s = sides(hierarchy)
  assert len(s) == n_expected, (len(s), sorted(s))
  for (resname, altloc, name), side in s.items():
    expected = 'cis' if name[1:] in [g[4] for g in GROUPS] else 'trans'
    assert side == expected, (resname, altloc, name, side)

def swapped_names(hierarchy):
  '''pdb text of hierarchy with every NH2 pair named opposite IUPAC ("1"
  trans), whatever the placement named: the two sites exchanged where "1" is
  cis. Names and atom order stay (as deposited in 7T5E: HD21 first, trans).'''
  h = hierarchy.deep_copy()
  s = sides(h)
  for (resname, altloc), at in residues(h).items():
    for g in GROUPS:
      if g[0] != resname: continue
      for e in ('H', 'D'):
        a, b = at.get(e + g[4]), at.get(e + g[5])
        if a is None or b is None: continue
        if a.parent().altloc != altloc: continue
        if s[(resname, altloc, e + g[4])] == 'cis':
          a.xyz, b.xyz = b.xyz, a.xyz
  return h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(source_info=None,
    lines=pdb_4wsc).crystal_symmetry())

def nh2_atoms(hierarchy):
  '''{(resname, altloc, name): atom} of the NH2 H/D.'''
  result = {}
  for (resname, altloc), at in residues(hierarchy).items():
    for g in GROUPS:
      if g[0] != resname: continue
      for e in ('H', 'D'):
        for s in g[4:]:
          a = at.get(e + s)
          if a is not None and a.parent().altloc == altloc:
            result[(resname, altloc, e + s)] = a
  return result

def assert_kept(h_in, h_out, values=False):
  '''Every supplied NH2 H/D written at its input site, with its values.'''
  a_in, a_out = nh2_atoms(h_in), nh2_atoms(h_out)
  assert len(a_in) > 0
  for k, a in a_in.items():
    assert k in a_out, k
    d = (matrix.col(a.xyz) - matrix.col(a_out[k].xyz)).length()
    assert d < 0.05, (k, d)
    if values:
      assert (a.occ, a.b) == (a_out[k].occ, a_out[k].b), (k, a.occ, a.b,
        a_out[k].occ, a_out[k].b)

def test_placed_names():
  # add with and without flip movers, and place: every placed NH2 H named by
  # IUPAC (Arg 4 H, Gln 2, Asn 2)
  for prefix, args in (('add', []), ('add_flips', FLIPS),
                       ('place', ['approach=place'])):
    assert_iupac(run('tst_r2_nh2_' + prefix, pdb_4wsc, args), 8)

def test_placed_names_neutron():
  # hd_and_h: the H and the D of each site take the same IUPAC name
  h = run('tst_r2_nh2_neutron', pdb_4wsc, NEUTRON + FLIPS)
  assert_iupac(h, 16)

def test_supplied_h_kept():
  # supplied H named opposite IUPAC: names and sites kept, with and without
  # flip movers
  placed = run('tst_r2_nh2_src', pdb_4wsc)
  src = swapped_names(placed)
  h_in = iotbx.pdb.input(source_info=None, lines=src).construct_hierarchy()
  for prefix, args in (('keep', []), ('keep_flips', FLIPS)):
    out = run('tst_r2_nh2_' + prefix, src, ['keep_existing_H=True'] + args)
    assert_kept(h_in, out)

def test_partial_pair():
  # a supplied HD21 on the trans site (non-IUPAC) without HD22: HD21 stays
  # at its site, the placed H takes HD22 at the other site
  placed = run('tst_r2_nh2_psrc', pdb_4wsc)
  src = [l for l in swapped_names(placed).splitlines()
         if not (l[12:16] == 'HD22' and l[17:20] == 'ASN')]
  src = '\n'.join(src) + '\n'
  h_in = iotbx.pdb.input(source_info=None, lines=src).construct_hierarchy()
  hd21 = nh2_atoms(h_in)[('ASN', '', 'HD21')]
  for prefix, args in (('partial', []), ('partial_flips', FLIPS)):
    out = run('tst_r2_nh2_' + prefix, src, ['keep_existing_H=True'] + args)
    a = nh2_atoms(out)
    d = (matrix.col(a[('ASN', '', 'HD21')].xyz) - matrix.col(hd21.xyz)).length()
    assert d < 0.05, (prefix, d)
    s = sides(out)
    assert s[('ASN', '', 'HD21')] == 'trans', prefix
    assert s[('ASN', '', 'HD22')] == 'cis', prefix

def test_supplied_d_kept():
  # supplied H/D pairs named opposite IUPAC, each site with its own
  # occupancy and B: names, sites and values kept (hd_and_h)
  placed = run('tst_r2_nh2_dsrc', pdb_4wsc, NEUTRON)
  h = iotbx.pdb.input(source_info=None,
    lines=swapped_names(placed)).construct_hierarchy()
  at = nh2_atoms(h)
  hs = sorted(k for k in at if k[2][0] == 'H')
  assert len(hs) == 8
  for i, k in enumerate(hs):
    d = [x for x in at if x[0] == k[0] and x[2] == 'D' + k[2][1:]]
    assert len(d) == 1, k
    a_h, a_d = at[k], at[d[0]]
    a_h.occ, a_d.occ = 0.30 + 0.02 * i, 0.70 - 0.02 * i
    a_h.b, a_d.b = 20 + i, 40 + i
  src = h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(source_info=None,
    lines=pdb_4wsc).crystal_symmetry())
  h_in = iotbx.pdb.input(source_info=None, lines=src).construct_hierarchy()
  for prefix, args in (('dkeep', []), ('dkeep_flips', FLIPS)):
    out = run('tst_r2_nh2_' + prefix, src,
      ['keep_existing_H=True'] + NEUTRON + args)
    assert_kept(h_in, out, values=True)

def test_restraints_fix_side():
  # a user ASN CIF (GeoStd with hh1 at 0 degrees, periodicity 1) fixes HD22
  # cis to CB: kept against IUPAC and reported; a fresh riding setup moves
  # nothing. Arg and Gln (GeoStd, periodicity 2) are still renamed.
  geostd = libtbx.env.find_in_repositories(
    relative_path=os.path.join('chem_data', 'geostd', 'a', 'data_ASN.cif'),
    test=os.path.isfile)
  if geostd is None:
    print('Skipping test_restraints_fix_side: no chem_data/geostd')
    return
  with open(geostd) as f:
    lines = f.read().splitlines()
  hh1 = [i for i, l in enumerate(lines) if l.split()[:2] == ['ASN', 'hh1']]
  assert len(hh1) == 1
  assert lines[hh1[0]].split()[2:] == ['CB', 'CG', 'ND2', 'HD22', '180.000',
    '.', '30.000', '2'], lines[hh1[0]]
  lines[hh1[0]] = ' ASN  hh1     CB  CG  ND2  HD22  0.000  .  30.000  1'
  cif = 'tst_r2_nh2_asn_p1.cif'
  with open(cif, 'w') as f:
    f.write('\n'.join(lines) + '\n')
  log = StringIO()
  h = run('tst_r2_nh2_fixed', pdb_4wsc, [cif, 'approach=place'], log=log)
  s = sides(h)
  assert s[('ASN', '', 'HD21')] == 'trans', s
  assert s[('ASN', '', 'HD22')] == 'cis', s
  for k, side in s.items():
    if k[0] == 'ASN': continue
    assert side == ('cis' if k[2][1:] in [g[4] for g in GROUPS] else 'trans'
                    ), (k, side)
  line = ('ASN A 433 HD21 trans to CB (IUPAC: cis), kept: torsion '
          'CB-CG-ND2-HD22 0.0, periodicity 1')
  with open('tst_r2_nh2_fixed.txt') as f:
    description = f.read()
  for text in (log.getvalue(), description):
    assert text.count(line) == 1, text
  # a fresh riding setup with the same restraints moves no NH2 H
  dm = DataManager()
  dm.process_model_file('tst_r2_nh2_fixed_out.pdb')
  dm.process_restraint_file(cif)
  m = dm.get_model()
  m.process(make_restraints=True)
  nh2 = flex.bool([a.name.strip()[1:] in [x for g in GROUPS for x in g[4:]]
                   for a in m.get_hierarchy().atoms()])
  before = m.get_sites_cart().select(nh2)
  m.setup_riding_h_manager(use_ideal_dihedral=True)
  m.idealize_h_riding()
  d = flex.sqrt((m.get_sites_cart().select(nh2) - before).dot())
  assert nh2.count(True) == 8 and flex.max(d) < 0.01, (nh2.count(True),
    flex.max(d))

if __name__ == '__main__':
  test_placed_names()
  test_placed_names_neutron()
  test_supplied_h_kept()
  test_partial_pair()
  test_supplied_d_kept()
  test_restraints_fix_side()
  print('OK')
