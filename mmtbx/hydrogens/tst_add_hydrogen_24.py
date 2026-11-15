from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import mmtbx.model
from cctbx import geometry_restraints
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from scitbx import matrix

# ------------------------------------------------------------------------------

def place(pdb_str):
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=pdb_str), log=null_out())
  o = reduce_hydrogen.place_hydrogens(model=m)
  o.run()
  return o.get_model()

def residuals(model):
  """
  {H name: (largest |residual| of its torsions, the same over those with a
  heavy third atom or None)}, by name, on a fresh model of the placed atoms
  with the restraints placement uses (CCD fallback included)
  """
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=model.get_hierarchy().as_pdb_string()), log=null_out())
  reduce_hydrogen.add_ccd_restraint_objects(m)
  m.process(make_restraints=True)
  g = m.get_restraints_manager().geometry
  atoms = m.get_hierarchy().atoms(); sc = m.get_sites_cart()
  is_h = lambda i: atoms[i].element.strip().upper() in ('H', 'D')
  result = {}
  for kind, proxies in (('var', g.dihedral_proxies),
                        ('const', g.const_dihedral_proxies or [])):
    for p in proxies:
      for ih, i_third in ((p.i_seqs[0], p.i_seqs[3]),
                          (p.i_seqs[3], p.i_seqs[0])):
        if not is_h(ih): continue
        d = geometry_restraints.dihedral(sites=[sc[i] for i in p.i_seqs],
          angle_ideal=0, weight=1).angle_model
        per = p.periodicity if kind == 'var' else 1
        r = abs(geometry_restraints.angle_delta_deg(angle_1=d,
          angle_2=p.angle_ideal, periodicity=max(1, per)))
        a, h = result.get(atoms[ih].name.strip(), (0., None))
        if not is_h(i_third): h = r if h is None else max(h, r)
        result[atoms[ih].name.strip()] = (max(a, r), h)
  return result

# ------------------------------------------------------------------------------

def test_000():
  '''
    1NA (CCD fallback): every torsion is period 1 and names its H. The methyl
    HM1-3 (on CM) and H81-83 (on C8) and every other H meet their own
    torsions by name. Before, reduce2's propeller naming swapped the second
    and third names by handedness only, while riding gives the target to
    whichever H's torsion comes last: with the riding sign fixed (A62) the
    methyl names came out 120 deg off; on master its sign error cancelled.
  '''
  res = residuals(place(pdb_str_000))
  for n in ('HM1', 'HM2', 'HM3', 'H81', 'H82', 'H83'):
    assert n in res, (n, sorted(res))
  bad = dict((n, r) for n, r in res.items() if r[0] > 1.)
  assert not bad, bad

# ------------------------------------------------------------------------------

def test_001():
  '''
    TRE: HO6 has a torsion with a heavy third atom (C5-C6-O6-HO6) and one with
    an H third (H61/H62-C6-O6-HO6). In placement the H third is still at its
    seed: the heavy one sets the riding target whatever the order (A62c).
    Before, the H-third torsion came last and won; with the sign fixed (A62)
    HO6 was 120 deg off C5-C6-O6-HO6.
  '''
  res = residuals(place(pdb_str_001))
  for n in ('HO6', 'HO6P'):
    assert res[n][1] is not None and res[n][1] < 1., (n, res[n])

# ------------------------------------------------------------------------------

def test_002():
  '''
    Standard methyls (ALA, ILE, LEU, MET, THR, VAL from 1us0, DT from 5gsk):
    their torsions are period 3 (VAL 2), which fit each H at every site, so
    reduce2 names them by handedness as before. Every H line is as master
    placed it, to the byte.
  '''
  m = place(pdb_str_002)
  lines = [l for l in m.get_hierarchy().as_pdb_string().splitlines()
    if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() == 'H']
  answer = answer_002.strip().splitlines()
  assert len(lines) == len(answer), (len(lines), len(answer))
  diff = [(a, b) for a, b in zip(answer, lines) if a != b]
  assert not diff, '\n'.join(a + '\n' + b for a, b in diff[:10])

# ------------------------------------------------------------------------------

def nerf(a, b, c, d, ang, tor):
  a, b, c = [matrix.col(v) for v in (a, b, c)]
  bc = (c - b).normalize(); n = (b - a).cross(bc).normalize(); m = n.cross(bc)
  A, T = math.radians(ang), math.radians(tor)
  return tuple(c + (-d * math.cos(A)) * bc + (d * math.sin(A) * math.cos(T)) * m
    + (d * math.sin(A) * math.sin(T)) * n)

def ethanol_methyl(torsions):
  """ethanol's methyl C1 with H11, H12, H13 at O1-C2-C1-H 180, 60, -60
  (residue atoms by name) and reference sites with the same names"""
  s = {'C1': (10., 10., 10.), 'C2': (11.52, 10., 10.)}
  s['O1'] = nerf((10., 11., 10.), s['C1'], s['C2'], 1.43, 109.5, 180.)
  for n, t in torsions.items():
    s[n] = nerf(s['O1'], s['C2'], s['C1'], 0.97, 109.5, t)
  lines = ['HETATM%5d %-4s LG9 A   1    %8.3f%8.3f%8.3f  1.00 20.00'
    '          %2s' % (k + 1, ' ' + n, s[n][0], s[n][1], s[n][2],
    'H' if n.startswith('H') else n[0]) for k, n in enumerate(sorted(s))]
  h = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  return dict((a.name.strip(), a) for a in h.atoms()), s

class _comp(object):
  def __init__(self, tors):
    from libtbx import group_args
    self.atom_list = [group_args(atom_id=n, type_symbol=n[0])
      for n in ('C1', 'C2', 'O1', 'H11', 'H12', 'H13', 'H21')]
    self.tor_list = [group_args(atom_id_1=a[0], atom_id_2=a[1], atom_id_3=a[2],
      atom_id_4=a[3], value_angle=v, value_angle_esd=esd, period=p)
      for a, v, esd, p in tors]

class _srv(object):
  def __init__(self, tors): self.comp = _comp(tors)
  def get_comp_comp_id_direct(self, resname): return self.comp

def test_003():
  '''
    The torsions that name a propeller's H: period 1 or CONST (esd 0), H at an
    end, heavy third atom. Not period 3 (GeoStd methyls, all standard
    residues), not an H third atom (still at its seed in placement).
  '''
  from mmtbx.hydrogens.reduce_hydrogen import _deciding_torsions
  t = _deciding_torsions('LG9', _srv([
    (('H11', 'C1', 'C2', 'O1'), 180., 10., 1),
    (('O1', 'C2', 'C1', 'H12'), 60., 0., 0),
    (('H13', 'C1', 'C2', 'O1'), -60., 10., 3),
    (('H13', 'C1', 'C2', 'H21'), 0., 10., 1)]), {})
  assert sorted(t) == ['H11', 'H12'], t
  assert t['H11'] == [(('H11', 'C1', 'C2', 'O1'), 180., 1)], t
  assert t['H12'] == [(('O1', 'C2', 'C1', 'H12'), 60., 0)], t

def test_004():
  '''
    Names from the deciding torsions over the three sites. With one H decided
    (H11 at 180), the other two follow the reference's handedness; with all
    three decided, each H goes to its own site whatever the reference.
  '''
  from mmtbx.hydrogens.reduce_hydrogen import _torsion_names
  hs = ('H11', 'H12', 'H13')
  # the model has the names one step round from the reference
  atoms, model = ethanol_methyl({'H11': 60., 'H12': -60., 'H13': 180.})
  ref = ethanol_methyl({'H11': 180., 'H12': 60., 'H13': -60.})[1]
  one = {'H11': [(('H11', 'C1', 'C2', 'O1'), 180., 1)]}
  new = dict((a.name.strip(), n) for a, n in
    _torsion_names(atoms, 'C1', ('C2',), hs, ref, one))
  # H13 (at 180) becomes H11; H11 (60) H12; H12 (-60) H13: the reference's
  # turning sense, not H11 at 180 with the other two swapped
  assert new == {'H13': 'H11', 'H11': 'H12', 'H12': 'H13'}, new
  # mirror reference handedness, the same decided H: the other two swap
  ref_m = ethanol_methyl({'H11': 180., 'H12': -60., 'H13': 60.})[1]
  new = dict((a.name.strip(), n) for a, n in
    _torsion_names(atoms, 'C1', ('C2',), hs, ref_m, one))
  assert new == {'H13': 'H11', 'H12': 'H12', 'H11': 'H13'}, new
  # all three decided: the torsions, whatever the reference's hand
  three = dict(one, H12=[(('H12', 'C1', 'C2', 'O1'), -60., 1)],
    H13=[(('H13', 'C1', 'C2', 'O1'), 60., 1)])
  new = dict((a.name.strip(), n) for a, n in
    _torsion_names(atoms, 'C1', ('C2',), hs, ref, three))
  assert new == {'H13': 'H11', 'H12': 'H12', 'H11': 'H13'}, new
  # no deciding torsion: None (the handedness rule applies)
  assert _torsion_names(atoms, 'C1', ('C2',), hs, ref, {}) is None

# ------------------------------------------------------------------------------

pdb_str_000 = """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  1NA A   1      15.887  14.694  15.119  1.00 20.00           C
HETATM    2  C2  1NA A   1      14.569  15.126  15.764  1.00 20.00           C
HETATM    3  C3  1NA A   1      13.407  14.697  14.862  1.00 20.00           C
HETATM    4  C4  1NA A   1      13.654  15.245  13.452  1.00 20.00           C
HETATM    5  C5  1NA A   1      15.042  14.805  12.982  1.00 20.00           C
HETATM    6  C6  1NA A   1      15.284  15.320  11.562  1.00 20.00           C
HETATM    7  C7  1NA A   1      14.902  15.107  18.178  1.00 20.00           C
HETATM    8  C8  1NA A   1      14.761  14.452  19.527  1.00 20.00           C
HETATM    9  CM  1NA A   1      18.172  14.629  15.320  1.00 20.00           C
HETATM   10  N2  1NA A   1      14.432  14.490  17.076  1.00 20.00           N
HETATM   11  O1  1NA A   1      16.975  15.060  15.970  1.00 20.00           O
HETATM   12  O3  1NA A   1      12.180  15.222  15.371  1.00 20.00           O
HETATM   13  O4  1NA A   1      12.663  14.737  12.558  1.00 20.00           O
HETATM   14  O5  1NA A   1      16.038  15.333  13.855  1.00 20.00           O
HETATM   15  O6  1NA A   1      16.576  14.904  11.118  1.00 20.00           O
HETATM   16  O7  1NA A   1      15.437  16.191  18.083  1.00 20.00           O
END
"""

pdb_str_001 = """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  TRE A   1      14.188  14.118  15.763  1.00 20.00           C
HETATM    2  C2  TRE A   1      13.558  13.064  14.849  1.00 20.00           C
HETATM    3  C3  TRE A   1      12.646  13.759  13.833  1.00 20.00           C
HETATM    4  C4  TRE A   1      11.595  14.581  14.586  1.00 20.00           C
HETATM    5  C5  TRE A   1      12.301  15.576  15.510  1.00 20.00           C
HETATM    6  C6  TRE A   1      11.258  16.359  16.310  1.00 20.00           C
HETATM    7  O2  TRE A   1      14.587  12.353  14.158  1.00 20.00           O
HETATM    8  O3  TRE A   1      11.996  12.778  13.024  1.00 20.00           O
HETATM    9  O4  TRE A   1      10.783  15.292  13.649  1.00 20.00           O
HETATM   10  O5  TRE A   1      13.156  14.866  16.408  1.00 20.00           O
HETATM   11  O6  TRE A   1      11.910  17.367  17.085  1.00 20.00           O
HETATM   12  C1P TRE A   1      15.812  15.882  15.762  1.00 20.00           C
HETATM   13  C2P TRE A   1      16.442  16.936  14.847  1.00 20.00           C
HETATM   14  C3P TRE A   1      17.353  16.240  13.832  1.00 20.00           C
HETATM   15  C4P TRE A   1      18.405  15.419  14.585  1.00 20.00           C
HETATM   16  C5P TRE A   1      17.699  14.424  15.510  1.00 20.00           C
HETATM   17  C6P TRE A   1      18.743  13.643  16.310  1.00 20.00           C
HETATM   18  O1  TRE A   1      15.000  15.000  14.984  1.00 20.00           O
HETATM   19  O2P TRE A   1      15.412  17.647  14.158  1.00 20.00           O
HETATM   20  O3P TRE A   1      18.003  17.221  13.022  1.00 20.00           O
HETATM   21  O4P TRE A   1      19.217  14.708  13.649  1.00 20.00           O
HETATM   22  O5P TRE A   1      16.844  15.135  16.408  1.00 20.00           O
HETATM   23  O6P TRE A   1      18.091  12.635  17.086  1.00 20.00           O
END
"""

pdb_str_002 = """
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM     18  N   ALA A   1      15.555  12.152  39.893  1.00 20.00           N
ATOM     19  CA  ALA A   1      14.153  11.844  39.611  1.00 20.00           C
ATOM     20  C   ALA A   1      14.014  10.631  38.728  1.00 20.00           C
ATOM     21  O   ALA A   1      14.644  10.406  37.752  1.00 20.00           O
ATOM     22  CB  ALA A   1      13.343  12.950  38.894  1.00 20.00           C
ATOM     63  N   ILE A   4      10.519  10.546  33.584  1.00 20.00           N
ATOM     64  CA  ILE A   4      10.571  11.633  32.619  1.00 20.00           C
ATOM     65  C   ILE A   4       9.172  11.830  32.066  1.00 20.00           C
ATOM     66  O   ILE A   4       8.408  10.887  31.787  1.00 20.00           O
ATOM     67  CB  ILE A   4      11.629  11.351  31.519  1.00 20.00           C
ATOM     68  CG1 ILE A   4      11.902  12.558  30.625  1.00 20.00           C
ATOM     69  CG2 ILE A   4      11.309  10.131  30.694  1.00 20.00           C
ATOM     70  CD1 ILE A   4      13.175  12.395  29.776  1.00 20.00           C
ATOM     82  N   LEU A   5       8.793  13.088  31.859  1.00 20.00           N
ATOM     83  CA  LEU A   5       7.495  13.415  31.309  1.00 20.00           C
ATOM     84  C   LEU A   5       7.475  13.258  29.812  1.00 20.00           C
ATOM     85  O   LEU A   5       8.283  13.868  29.108  1.00 20.00           O
ATOM     86  CB  LEU A   5       7.178  14.879  31.691  1.00 20.00           C
ATOM     87  CG  LEU A   5       5.778  15.362  31.308  1.00 20.00           C
ATOM     88  CD1 LEU A   5       4.608  14.558  31.862  1.00 20.00           C
ATOM     89  CD2 LEU A   5       5.657  16.811  31.757  1.00 20.00           C
ATOM      1  N   MET A   0      18.709  11.104  41.491  1.00 20.00           N
ATOM      2  CA  MET A   0      17.937  11.685  40.368  1.00 20.00           C
ATOM      3  C   MET A   0      16.479  11.264  40.309  1.00 20.00           C
ATOM      4  O   MET A   0      16.161  10.141  40.705  1.00 20.00           O
ATOM      5  CB  MET A   0      18.672  11.229  39.111  1.00 20.00           C
ATOM      6  CG  MET A   0      18.004  11.534  37.807  1.00 20.00           C
ATOM      7  SD  MET A   0      18.852  11.085  36.316  1.00 20.00           S
ATOM      8  CE  MET A   0      20.086  12.267  35.949  1.00 20.00           C
ATOM    321  N   THR A  19      14.553  -7.200  27.782  1.00 20.00           N
ATOM    322  CA  THR A  19      14.176  -8.603  27.740  1.00 20.00           C
ATOM    323  C   THR A  19      14.578  -9.324  26.462  1.00 20.00           C
ATOM    324  O   THR A  19      14.284 -10.524  26.332  1.00 20.00           O
ATOM    325  CB  THR A  19      12.675  -8.802  27.996  1.00 20.00           C
ATOM    326  OG1 THR A  19      11.927  -8.210  26.942  1.00 20.00           O
ATOM    327  CG2 THR A  19      12.241  -8.267  29.356  1.00 20.00           C
ATOM    452  N   VAL A  27      13.214 -15.629  32.266  1.00 20.00           N
ATOM    453  CA  VAL A  27      12.767 -14.249  32.254  1.00 20.00           C
ATOM    454  C   VAL A  27      11.641 -14.003  33.257  1.00 20.00           C
ATOM    455  O   VAL A  27      11.487 -12.880  33.732  1.00 20.00           O
ATOM    456  CB  VAL A  27      12.357 -13.813  30.831  1.00 20.00           C
ATOM    457  CG1 VAL A  27      11.059 -14.438  30.367  1.00 20.00           C
ATOM    458  CG2 VAL A  27      12.301 -12.293  30.740  1.00 20.00           C
ATOM     65  P    DT B   3      10.011  19.045   4.108  1.00 20.00           P
ATOM     66  OP1  DT B   3      10.235  19.681   5.425  1.00 20.00           O
ATOM     67  OP2  DT B   3      10.906  17.959   3.643  1.00 20.00           O
ATOM     68  O5'  DT B   3       8.532  18.470   4.070  1.00 20.00           O
ATOM     69  C5'  DT B   3       7.452  19.242   4.612  1.00 20.00           C
ATOM     70  C4'  DT B   3       6.153  18.502   4.416  1.00 20.00           C
ATOM     71  O4'  DT B   3       5.899  18.337   3.006  1.00 20.00           O
ATOM     72  C3'  DT B   3       6.130  17.086   4.976  1.00 20.00           C
ATOM     73  O3'  DT B   3       5.798  17.096   6.359  1.00 20.00           O
ATOM     74  C2'  DT B   3       5.039  16.439   4.154  1.00 20.00           C
ATOM     75  C1'  DT B   3       5.228  17.096   2.794  1.00 20.00           C
ATOM     76  N1   DT B   3       6.037  16.289   1.874  1.00 20.00           N
ATOM     77  C2   DT B   3       5.382  15.302   1.178  1.00 20.00           C
ATOM     78  O2   DT B   3       4.185  15.095   1.285  1.00 20.00           O
ATOM     79  N3   DT B   3       6.183  14.567   0.341  1.00 20.00           N
ATOM     80  C4   DT B   3       7.540  14.721   0.135  1.00 20.00           C
ATOM     81  O4   DT B   3       8.121  13.998  -0.668  1.00 20.00           O
ATOM     82  C5   DT B   3       8.169  15.755   0.922  1.00 20.00           C
ATOM     83  C7   DT B   3       9.642  15.978   0.784  1.00 20.00           C
ATOM     84  C6   DT B   3       7.398  16.472   1.749  1.00 20.00           C
END
"""

answer_002 = """
ATOM      6  H   ALA A   1      15.810  12.966  39.784  1.00 20.00           H
ATOM      7  HA  ALA A   1      13.782  11.722  40.499  1.00 20.00           H
ATOM      8  HB1 ALA A   1      13.341  13.748  39.446  1.00 20.00           H
ATOM      9  HB3 ALA A   1      12.434  12.639  38.759  1.00 20.00           H
ATOM     10  HB2 ALA A   1      13.758  13.141  38.038  1.00 20.00           H
ATOM     19  H   ILE A   4      11.139  10.576  34.179  1.00 20.00           H
ATOM     20  HA  ILE A   4      10.851  12.462  33.038  1.00 20.00           H
ATOM     21  HB  ILE A   4      12.446  11.166  32.008  1.00 20.00           H
ATOM     22 HG13 ILE A   4      12.011  13.345  31.182  1.00 20.00           H
ATOM     23 HG12 ILE A   4      11.152  12.680  30.022  1.00 20.00           H
ATOM     24 HG21 ILE A   4      11.270   9.357  31.277  1.00 20.00           H
ATOM     25 HG23 ILE A   4      12.004  10.009  30.028  1.00 20.00           H
ATOM     26 HG22 ILE A   4      10.452  10.260  30.258  1.00 20.00           H
ATOM     27 HD11 ILE A   4      13.295  13.189  29.232  1.00 20.00           H
ATOM     28 HD13 ILE A   4      13.077  11.616  29.207  1.00 20.00           H
ATOM     29 HD12 ILE A   4      13.936  12.280  30.367  1.00 20.00           H
ATOM     38  H   LEU A   5       9.280  13.775  32.033  1.00 20.00           H
ATOM     39  HA  LEU A   5       6.816  12.821  31.664  1.00 20.00           H
ATOM     40  HB3 LEU A   5       7.817  15.457  31.245  1.00 20.00           H
ATOM     41  HB2 LEU A   5       7.264  14.969  32.653  1.00 20.00           H
ATOM     42  HG  LEU A   5       5.697  15.251  30.348  1.00 20.00           H
ATOM     43 HD11 LEU A   5       3.778  14.955  31.555  1.00 20.00           H
ATOM     44 HD13 LEU A   5       4.644  14.575  32.831  1.00 20.00           H
ATOM     45 HD12 LEU A   5       4.675  13.644  31.545  1.00 20.00           H
ATOM     46 HD21 LEU A   5       5.784  16.856  32.718  1.00 20.00           H
ATOM     47 HD23 LEU A   5       4.775  17.141  31.525  1.00 20.00           H
ATOM     48 HD22 LEU A   5       6.336  17.338  31.308  1.00 20.00           H
ATOM     57  H   MET A   0      19.305  11.632  41.816  1.00 20.00           H
ATOM     58  HA  MET A   0      17.900  12.650  40.455  1.00 20.00           H
ATOM     59  HB2 MET A   0      19.540  11.661  39.095  1.00 20.00           H
ATOM     60  HB3 MET A   0      18.781  10.266  39.157  1.00 20.00           H
ATOM     61  HG3 MET A   0      17.858  12.492  37.766  1.00 20.00           H
ATOM     62  HG2 MET A   0      17.153  11.068  37.794  1.00 20.00           H
ATOM     63  HE1 MET A   0      19.668  13.135  35.833  1.00 20.00           H
ATOM     64  HE3 MET A   0      20.539  12.007  35.132  1.00 20.00           H
ATOM     65  HE2 MET A   0      20.720  12.300  36.683  1.00 20.00           H
ATOM     73  H   THR A  19      14.578  -6.868  28.575  1.00 20.00           H
ATOM     74  HA  THR A  19      14.672  -9.044  28.448  1.00 20.00           H
ATOM     75  HB  THR A  19      12.483  -9.753  28.017  1.00 20.00           H
ATOM     76  HG1 THR A  19      12.092  -7.387  26.899  1.00 20.00           H
ATOM     77 HG21 THR A  19      12.724  -8.725  30.062  1.00 20.00           H
ATOM     78 HG23 THR A  19      11.290  -8.410  29.481  1.00 20.00           H
ATOM     79 HG22 THR A  19      12.426  -7.317  29.413  1.00 20.00           H
ATOM     87  H   VAL A  27      12.708 -16.164  31.822  1.00 20.00           H
ATOM     88  HA  VAL A  27      13.511 -13.694  32.536  1.00 20.00           H
ATOM     89  HB  VAL A  27      13.037 -14.139  30.221  1.00 20.00           H
ATOM     90 HG11 VAL A  27      11.156 -15.403  30.366  1.00 20.00           H
ATOM     91 HG13 VAL A  27      10.860 -14.124  29.471  1.00 20.00           H
ATOM     92 HG12 VAL A  27      10.348 -14.178  30.974  1.00 20.00           H
ATOM     93 HG21 VAL A  27      11.649 -11.962  31.377  1.00 20.00           H
ATOM     94 HG23 VAL A  27      12.042 -12.040  29.840  1.00 20.00           H
ATOM     95 HG22 VAL A  27      13.177 -11.932  30.946  1.00 20.00           H
ATOM    116 H5''  DT B   3       7.408  20.100   4.162  1.00 20.00           H
ATOM    117  H5'  DT B   3       7.602  19.390   5.559  1.00 20.00           H
ATOM    118  H4'  DT B   3       5.447  19.025   4.826  1.00 20.00           H
ATOM    119  H3'  DT B   3       6.978  16.632   4.852  1.00 20.00           H
ATOM    120  H2'  DT B   3       5.158  15.478   4.102  1.00 20.00           H
ATOM    121 H2''  DT B   3       4.161  16.629   4.520  1.00 20.00           H
ATOM    122  H1'  DT B   3       4.360  17.291   2.407  1.00 20.00           H
ATOM    123  H3   DT B   3       5.797  13.941  -0.105  1.00 20.00           H
ATOM    124 H5M1  DT B   3       9.797  16.887   0.483  1.00 20.00           H
ATOM    125 H5M2  DT B   3       9.997  15.350   0.136  1.00 20.00           H
ATOM    126 H5M3  DT B   3      10.065  15.839   1.646  1.00 20.00           H
ATOM    127  H6   DT B   3       7.803  17.130   2.268  1.00 20.00           H
"""

if (__name__ == "__main__"):
  test_000()
  test_001()
  test_002()
  test_003()
  test_004()
  print('OK')
