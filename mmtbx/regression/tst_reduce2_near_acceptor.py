"""reduce2 near-acceptor sampling for single-H rotators. Next to a close
acceptor, Reduce (RotDonor.cpp) adds coarse orientations 15 and 30 degrees to
either side of the acceptor direction and searches about 10 degrees around the
chosen one; reduce2 had only the acceptor direction and +-5 degrees, so an OH
pointing straight at a close acceptor bumped it and the Optimizer settled
elsewhere. The added angles come after the primary candidates; when one of
them wins the coarse step, the Optimizer also refines the best primary
candidate and keeps the better pose.

Tests:
- candidate generation on constructed rotators (thresholds, wrap-aware drop,
  order, primary candidates kept, fine ranges);
- 7PRG FUC A 203 O3, bridging two Ca, hd_and_h: H-bond to Asp99 OD2, away from
  both Ca; the same excerpt in X-ray mode as a guard (unchanged there);
- 1R27 Tyr C 1201 HH (no metal): out of the too-short contact with Asp C 1047 OD1;
- 1R27 Tyr D 415 HH (hd_and_h, singleton): an added angle wins the coarse step,
  the primary candidate's better pose is kept; the stored score and pose match a
  fresh rescoring of the final coordinates.
"""
from __future__ import absolute_import, division, print_function
import math, os, re
import iotbx.pdb
from iotbx import pdb
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
import mmtbx_probe_ext as probe
from mmtbx.reduce import Movers, Optimizers
from mmtbx_reduce_ext import RotateAtomDegreesAroundAxisDir

#-----------------------------------------------------------------------------
# Constructed rotators: OG-HG on CB with two friends (as Movers.Test()); the H
# starts anti to CA. Each acceptor is a water O (radius 1.4 A) at distance d
# from OG, theta degrees around the CB->OG axis from the H. Reduce's overlap is
# d - (1.05 + 1.4 + 1.0): +-15 degrees below d = 2.75 A, +-30 below d = 2.55 A.

def rotator(acceptors):
  """acceptors: list of (resseq, theta, d); returns the MoverSingleHydrogenRotator"""
  ch = pdb.hierarchy.chain(); ch.id = 'A'
  def res(n, name, atoms):
    ag = pdb.hierarchy.atom_group(); ag.resname = name
    for a in atoms: ag.append_atom(a)
    rg = pdb.hierarchy.residue_group(); rg.resseq = '%4d' % n
    rg.append_atom_group(ag); ch.append_residue_group(rg)
  def at(name, el, xyz):
    a = pdb.hierarchy.atom(); a.name = name; a.element = el; a.xyz = xyz
    return a
  h = at(' HG ', 'H', (1.0, 1.0, 1.0)); o = at(' OG ', 'O', (0.0, 0.0, 0.0))
  c = at(' CB ', 'C', (0.0, 0.0, -1.0)); f1 = at(' CA ', 'C', (1.0, 0.0, -2.0))
  f2 = at(' C  ', 'C', (-1.0, 0.0, -2.0))
  res(1, 'SER', [h, o, c, f1, f2])
  axis = flex.vec3_double([c.xyz, (0.0, 0.0, 1.0)])
  acc = []
  for n, theta, d in acceptors:
    a = at(' O  ', 'O', (-d, 0.0, 0.0))
    a.xyz = RotateAtomDegreesAroundAxisDir(axis[0], axis[1], a, theta)
    res(n, 'HOH', [a]); acc.append(a)
  m = pdb.hierarchy.model(); m.append_chain(ch); m.atoms().reset_i_seq()
  infos = []
  for a in m.atoms():
    e = probe.ExtraAtomInfo(); e.vdwRadius = 1.4 if a.parent().resname == 'HOH' else 1.0
    infos.append(e)
  extras = probe.ExtraAtomInfoMap(m.atoms(), infos)
  bonded = {h: [o], o: [h, c], c: [o, f1, f2], f1: [c], f2: [c]}
  class Item:
    def __init__(self, n, a2): self.n = n; self.a2 = a2
  mover = Movers.MoverSingleHydrogenRotator(h, bonded, extras,
    {h.i_seq: Item(0, f1.i_seq)}, acc)
  mover._model = m  # keep the hierarchy that the atoms belong to
  return mover

def split(mover):
  """(primary, added) coarse angles"""
  n = mover.NumPrimaryCoarsePositions()
  return [round(a, 6) for a in mover._coarseAngles[:n]], [round(a, 6) for a in mover._coarseAngles[n:]]

def test_candidates():
  far = split(rotator([(2, 40, 3.0)]))
  assert far[1] == [], far
  # thresholds; the primary candidates do not depend on the acceptor distance
  for d, added in ((2.751, []), (2.749, [25, 55]), (2.551, [25, 55]), (2.549, [25, 55, 10, 70])):
    p, a = split(rotator([(2, 40, d)]))
    assert p == far[0], (d, p, far)
    assert a == added, (d, a)
  # wrap-aware drop: 178 + 30 = 208 degrees lies 2 degrees from the primary -150
  p, a = split(rotator([(2, 178, 2.5)]))
  assert -150 in p and a == [163, 193, 148], (p, a)
  # HEAD's near-duplicate candidates (0 and 3 degrees) are kept
  p, a = split(rotator([(2, 3, 2.5)]))
  assert p[:2] == [0, 3] and p == split(rotator([(2, 3, 3.0)]))[0], p
  assert a == [-12, 18, -27, 33], a
  # two acceptors: offsets in acceptor order; the second one's (28, 58, 13, 73) all lie
  # within 5 degrees of the first one's and are dropped; the input order does not matter
  one = split(rotator([(2, 40, 2.5), (3, 43, 2.5)]))
  two = split(rotator([(3, 43, 2.5), (2, 40, 2.5)]))
  assert one == two, (one, two)
  assert one[1] == [25, 55, 10, 70], one
  one = split(rotator([(2, 40, 2.5), (3, -100, 2.6)]))
  assert one[1] == [25, 55, 10, 70, -115, -85], one
  # fine ranges: -10..+9 degrees for the single-H rotator, half the coarse step otherwise
  mover = rotator([(2, 40, 2.5)])
  assert mover._fineAngles == [-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7, -8, 8, -9, 9, -10], (
    mover._fineAngles)
  assert all(len(mover.FinePositions(i).positions) == 19 for i in range(len(mover._coarseAngles)))
  atoms = [pdb.hierarchy.atom() for i in range(2)]
  atoms[1].xyz = (1.0, 0.0, 0.0)
  base = Movers._MoverRotator(atoms, flex.vec3_double([(0, 0, 0), (0, 0, 1)]), 0, 0, 180,
    coarseStepDegrees=15.0)
  assert base._fineAngles == [-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7], base._fineAngles
  print("test_candidates OK")

#-----------------------------------------------------------------------------
# Excerpts, run through reduce2 (approach=add)

def run(prefix, pdb_str, args):
  fn = prefix + '.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  for x in (prefix + '_out.pdb', prefix + '_out.txt'):
    if os.path.exists(x): os.remove(x)
  run_program(program_class=reduce2.Program, logger=null_out(), args=[fn,
    'output.filename=' + prefix + '_out.pdb', 'output.description_file_name=' + prefix + '_out.txt',
    'overwrite=True', 'approach=add'] + args)
  with open(prefix + '_out.txt') as f:
    desc = f.read()
  h = iotbx.pdb.input(prefix + '_out.pdb').construct_hierarchy()
  return h, desc

def atoms(h, chain, resseq, name):
  return [a for a in h.atoms() if a.parent().parent().parent().id == chain and
    a.parent().parent().resseq.strip() == resseq and a.name.strip() == name]

def angle(a, b, c):
  return math.degrees((matrix.col(a.xyz) - matrix.col(b.xyz)).angle(matrix.col(c.xyz) - matrix.col(b.xyz)))

def pose(desc, label):
  """final score and angle of a single-H rotator from the report"""
  l = [l for l in desc.splitlines() if 'SingleHydrogenRotator at ' + label + ' ' in l]
  assert len(l) == 1, (label, desc)
  m = re.search(r'final score: (\S+) pose Angle (\S+) deg', l[0])
  return float(m.group(1)), float(m.group(2))

def check_fuc_o3(h, xh):
  o3 = atoms(h, 'A', '203', 'O3')[0]
  c3 = atoms(h, 'A', '203', 'C3')[0]
  od2 = atoms(h, 'A', '99', 'OD2')[0]
  od1 = atoms(h, 'A', '104', 'OD1')[0]
  hs = [a for a in h.atoms() if a.parent().resname == 'FUC' and a.element.strip() in ('H', 'D')
        and a.distance(o3) < 1.2]
  assert len(hs) >= 1
  for ho3 in hs:
    # the H-bond partner is Asp99 OD2 (Asp104 OD1 before)
    assert ho3.distance(od2) < 2.0 and ho3.distance(od2) < ho3.distance(od1), (
      ho3.distance(od2), ho3.distance(od1))
    for ca in (atoms(h, 'A', '201', 'CA')[0], atoms(h, 'A', '202', 'CA')[0]):
      assert angle(ca, o3, ho3) > 90, angle(ca, o3, ho3)
    # riding geometry
    assert abs(ho3.distance(o3) - xh) < 0.005, ho3.distance(o3)
    assert abs(angle(c3, o3, ho3) - 108.7) < 1.0, angle(c3, o3, ho3)

def test_7prg_fuc_o3_neutron():
  h, desc = run('tst_reduce2_near_acceptor_7prg_hd', pdb_7prg_fuc, ['neutron_option=hd_and_h'])
  check_fuc_o3(h, 0.98)
  print("test_7prg_fuc_o3_neutron OK")

def test_7prg_fuc_o3_xray():
  # guard: with the X-ray O-H the H already found Asp99 OD2
  h, desc = run('tst_reduce2_near_acceptor_7prg_x', pdb_7prg_fuc, [])
  check_fuc_o3(h, 0.85)
  print("test_7prg_fuc_o3_xray OK")

def test_tyr_short_contact():
  # no metal: HH pointed straight at Asp C 1047 OD1 (1.53 A, score 2.11)
  h, desc = run('tst_reduce2_near_acceptor_tyr1201', pdb_1r27_tyr1201, ['neutron_option=hd_and_h'])
  score, ang = pose(desc, 'chain C TYR 1201 HH')
  od1 = atoms(h, 'C', '1047', 'OD1')[0]
  hs = atoms(h, 'C', '1201', 'HH') + atoms(h, 'C', '1201', 'DH')
  assert len(hs) >= 1
  for hh in hs:
    assert 1.65 < hh.distance(od1) < 2.0, hh.distance(od1)
  assert score > 12.0, score
  print("test_tyr_short_contact OK")

def test_safeguard_singleton():
  # an added angle (38 degrees, coarse 12.50) wins the coarse step but ends at
  # 28.2 degrees (13.16) after its fine search; the best primary candidate
  # (-0.1 degrees) ends at -5.1 degrees (14.49), which is kept
  h, desc = run('tst_reduce2_near_acceptor_tyr415', pdb_1r27_tyr415, ['neutron_option=hd_and_h'])
  score, ang = pose(desc, 'chain D TYR 415 HH')
  assert abs(ang - (-5.1)) < 3.0, ang
  assert score >= 14.49 - 0.005, score
  print("test_safeguard_singleton OK")

def test_fallback_state_cache():
  # after the switch to the primary candidate's pose, the coordinates, the pose and
  # the stored score match a fresh rescoring of the final coordinates
  real = Optimizers.OptimizerC
  seen = []
  class Capture(object):
    def __init__(self, *args):
      self._r = real(*args)
      (self.verb, self.mag, self.maxvdw, self.minocc, self.probe, self.density, self.atoms,
       self.excl, self.scorer, self.cache, self.aml, self.sq, self.extra, self.dele) = args
    def __getattr__(self, n):
      return getattr(self._r, n)
    def OptimizeSingleMoverFine(self, m):
      before = self._r.GetCoarseLocation(m)
      ret = self._r.OptimizeSingleMoverFine(m)
      coarse = m.CoarsePositions()
      if coarse.atoms[0].parent().resname == 'TYR' and coarse.atoms[0].name.strip() == 'HH':
        c, f = self._r.GetCoarseLocation(m), self._r.GetFineLocation(m)
        poses = m.FinePositions(c) if f >= 0 else coarse
        i = f if f >= 0 else c
        fresh = self.mag * poses.preferenceEnergies[i]
        for a in coarse.atoms:
          fresh += self.scorer.score_dots(a, self.minocc, self.sq,
            self.extra.getMappingFor(a).vdwRadius + self.maxvdw, self.probe, self.excl[a.i_seq],
            self._r.GetDots(a.i_seq), self.density, False, False).totalScore()
        moved = max((matrix.col(a.xyz) - matrix.col(p)).length()
                    for a, p in zip(coarse.atoms, poses.positions[i]))
        seen.append(dict(before=before, coarse=c, fine=f, primary=m.NumPrimaryCoarsePositions(),
          stored=self._r.GetHighScore(m), fresh=fresh, moved=moved,
          pose=m.PoseDescription(c, f, True)))
      return ret
  Optimizers.OptimizerC = Capture
  try:
    h, desc = run('tst_reduce2_near_acceptor_tyr415_cache', pdb_1r27_tyr415, ['neutron_option=hd_and_h'])
  finally:
    Optimizers.OptimizerC = real
  assert len(seen) == 1, seen
  s = seen[0]
  # the coarse winner was an added angle; the kept pose is around a primary one
  assert s['before'] >= s['primary'] and s['coarse'] < s['primary'], s
  assert s['moved'] < 1e-6, s
  assert abs(s['stored'] - s['fresh']) < 1e-6, s
  score, ang = pose(desc, 'chain D TYR 415 HH')
  assert abs(score - s['stored']) < 0.005 and s['pose'] == 'Angle %.1f deg .' % ang, (score, ang, s)
  print("test_fallback_state_cache OK")

pdb_7prg_fuc = """\
CRYST1   52.900   73.870   55.003  90.00  94.58  90.00 P 1 21 1
ATOM    165  N   ASN A  21      32.119  20.515  53.277  1.00 11.27           N
ATOM    166  CA  ASN A  21      31.916  21.253  54.508  1.00 11.28           C
ATOM    167  C   ASN A  21      31.317  20.327  55.563  1.00 12.41           C
ATOM    168  O   ASN A  21      30.166  20.510  55.948  1.00 13.56           O
ATOM    169  CB  ASN A  21      31.001  22.468  54.351  1.00 12.02           C
ATOM    170  CG  ASN A  21      31.628  23.599  53.524  1.00 12.72           C
ATOM    171  OD1 ASN A  21      32.400  24.418  54.007  1.00 14.22           O
ATOM    172  ND2 ASN A  21      31.253  23.640  52.268  1.00 17.28           N
ATOM    173  N   SER A  22      32.120  19.376  56.065  1.00 14.99           N
ATOM    174  CA  SER A  22      31.602  18.375  56.980  1.00 14.35           C
ATOM    175  C   SER A  22      32.766  17.609  57.614  1.00 12.79           C
ATOM    176  O   SER A  22      33.853  17.519  57.047  1.00 12.86           O
ATOM    177  CB  SER A  22      30.683  17.368  56.269  1.00 15.38           C
ATOM    178  OG  SER A  22      30.323  16.286  57.127  1.00 16.76           O
ATOM    179  N  ASER A  23      32.489  17.005  58.771  0.57 12.35           N
ATOM    180  CA ASER A  23      33.458  16.101  59.386  0.57 15.29           C
ATOM    181  C  ASER A  23      33.533  14.778  58.665  0.57 16.89           C
ATOM    182  O  ASER A  23      34.509  14.040  58.857  0.57 15.21           O
ATOM    183  CB ASER A  23      33.113  15.839  60.852  0.57 16.98           C
ATOM    184  OG ASER A  23      33.508  16.943  61.655  0.57 14.04           O
ATOM    185  N  BSER A  23      32.514  17.033  58.790  0.43 12.43           N
ATOM    186  CA BSER A  23      33.511  16.123  59.351  0.43 15.30           C
ATOM    187  C  BSER A  23      33.643  14.877  58.497  0.43 17.16           C
ATOM    188  O  BSER A  23      34.744  14.309  58.403  0.43 11.34           O
ATOM    189  CB BSER A  23      33.155  15.724  60.783  0.43 16.88           C
ATOM    190  OG BSER A  23      32.017  14.902  60.768  0.43 18.25           O
ATOM    733  N   GLU A  95      22.924  17.094  52.603  1.00 15.51           N
ATOM    734  CA  GLU A  95      22.949  16.979  54.060  1.00 11.34           C
ATOM    735  C   GLU A  95      23.246  15.543  54.492  1.00 13.20           C
ATOM    736  O   GLU A  95      22.497  14.625  54.170  1.00 13.18           O
ATOM    737  CB  GLU A  95      21.630  17.457  54.646  1.00  9.87           C
ATOM    738  CG  GLU A  95      21.556  17.380  56.166  1.00 14.98           C
ATOM    739  CD  GLU A  95      22.683  18.064  56.853  1.00 17.67           C
ATOM    740  OE1 GLU A  95      22.940  19.233  56.541  1.00 16.81           O
ATOM    741  OE2 GLU A  95      23.348  17.425  57.702  1.00 15.59           O
ATOM    742  N   ASP A  96      24.328  15.365  55.232  1.00 15.95           N
ATOM    743  CA  ASP A  96      24.765  14.045  55.669  1.00 17.88           C
ATOM    744  C   ASP A  96      24.415  13.764  57.124  1.00 16.16           C
ATOM    745  O   ASP A  96      24.712  12.664  57.617  1.00 16.49           O
ATOM    746  CB  ASP A  96      26.282  13.868  55.456  1.00 15.00           C
ATOM    747  CG  ASP A  96      27.160  14.781  56.330  1.00 14.72           C
ATOM    748  OD1 ASP A  96      26.776  15.225  57.436  1.00 16.98           O
ATOM    749  OD2 ASP A  96      28.277  15.062  55.881  1.00 14.71           O
ATOM    750  N   GLY A  97      23.805  14.728  57.810  1.00 15.66           N
ATOM    751  CA  GLY A  97      23.614  14.658  59.236  1.00 19.77           C
ATOM    752  C   GLY A  97      22.308  15.221  59.683  1.00 20.00           C
ATOM    753  O   GLY A  97      21.274  14.989  59.048  1.00 17.76           O
ATOM    761  N   ASP A  99      22.049  18.680  61.261  1.00 16.76           N
ATOM    762  CA  ASP A  99      21.978  20.115  61.480  1.00 14.41           C
ATOM    763  C   ASP A  99      21.317  20.878  60.319  1.00 17.31           C
ATOM    764  O   ASP A  99      21.015  22.066  60.468  1.00 17.75           O
ATOM    765  CB  ASP A  99      23.371  20.687  61.785  1.00 21.61           C
ATOM    766  CG  ASP A  99      24.320  20.651  60.599  1.00 21.79           C
ATOM    767  OD1 ASP A  99      24.052  19.932  59.642  1.00 16.80           O
ATOM    768  OD2 ASP A  99      25.364  21.297  60.670  1.00 16.89           O
ATOM    777  N   ASP A 101      22.362  22.111  57.488  1.00 19.68           N
ATOM    778  CA  ASP A 101      23.010  23.234  56.791  1.00 15.52           C
ATOM    779  C   ASP A 101      23.118  22.977  55.281  1.00 19.73           C
ATOM    780  O   ASP A 101      23.294  23.927  54.520  1.00 15.79           O
ATOM    781  CB  ASP A 101      24.377  23.612  57.413  1.00 13.41           C
ATOM    782  CG  ASP A 101      25.417  22.500  57.391  1.00 10.61           C
ATOM    783  OD1 ASP A 101      25.104  21.342  57.003  1.00 13.37           O
ATOM    784  OD2 ASP A 101      26.552  22.826  57.774  1.00 14.83           O
ATOM    797  N   ASN A 103      25.248  21.969  53.265  1.00 14.96           N
ATOM    798  CA  ASN A 103      26.469  22.363  52.592  1.00 11.94           C
ATOM    799  C   ASN A 103      27.540  21.269  52.645  1.00 14.17           C
ATOM    800  O   ASN A 103      28.699  21.518  52.343  1.00 12.33           O
ATOM    801  CB  ASN A 103      27.003  23.643  53.240  1.00 13.00           C
ATOM    802  CG  ASN A 103      27.607  23.386  54.617  1.00 16.25           C
ATOM    803  OD1 ASN A 103      27.597  22.268  55.102  1.00 15.19           O
ATOM    804  ND2 ASN A 103      28.145  24.394  55.224  1.00 16.14           N
ATOM    805  N   ASP A 104      27.179  20.071  53.066  1.00 14.09           N
ATOM    806  CA  ASP A 104      28.190  19.107  53.439  1.00 12.05           C
ATOM    807  C   ASP A 104      28.983  18.669  52.223  1.00 11.65           C
ATOM    808  O   ASP A 104      30.216  18.554  52.282  1.00 12.49           O
ATOM    809  CB  ASP A 104      27.518  17.942  54.163  1.00 12.29           C
ATOM    810  CG  ASP A 104      26.834  18.398  55.472  1.00 12.48           C
ATOM    811  OD1 ASP A 104      27.363  19.322  56.161  1.00 15.88           O
ATOM    812  OD2 ASP A 104      25.759  17.851  55.748  1.00 12.48           O
TER
ATOM   1768  N   GLY B 114      31.365  25.232  56.738  1.00 16.78           N
ATOM   1769  CA  GLY B 114      31.585  24.498  57.962  1.00 17.57           C
ATOM   1770  C   GLY B 114      30.686  23.303  58.135  1.00 18.46           C
ATOM   1771  O   GLY B 114      29.634  23.118  57.479  1.00 15.73           O
ATOM   1772  OXT GLY B 114      31.046  22.486  58.995  1.00 19.44           O
TER
HETATM 3559 CA    CA A 201      25.177  19.096  57.681  1.00 15.05          CA
HETATM 3560 CA    CA A 202      28.204  21.197  57.073  1.00 13.66          CA
HETATM 3561  C1  FUC A 203      28.806  17.074  59.972  1.00 13.89           C
HETATM 3562  C2  FUC A 203      27.962  18.021  59.116  1.00 11.73           C
HETATM 3563  C3  FUC A 203      27.746  19.309  59.807  1.00 14.66           C
HETATM 3564  C4  FUC A 203      29.025  19.990  60.179  1.00 15.20           C
HETATM 3565  C5  FUC A 203      29.912  19.006  60.970  1.00 14.32           C
HETATM 3566  C6  FUC A 203      31.307  19.533  61.159  1.00 14.25           C
HETATM 3567  O1  FUC A 203      28.177  16.809  61.181  1.00 18.99           O
HETATM 3568  O2  FUC A 203      26.677  17.445  58.898  1.00 12.84           O
HETATM 3569  O3  FUC A 203      26.909  20.162  58.970  1.00 15.03           O
HETATM 3570  O4  FUC A 203      29.631  20.343  58.941  1.00 13.21           O
HETATM 3571  O5  FUC A 203      30.071  17.728  60.293  1.00 14.70           O
HETATM 3677  O   HOH A 362      26.756  20.780  63.163  1.00 44.59           O
HETATM 3678  O   HOH A 363      27.293  13.526  59.733  1.00 27.38           O
HETATM 3690  O   HOH A 375      25.022  16.997  61.522  1.00 26.03           O
END
"""

pdb_1r27_tyr1201 = """\
CRYST1  127.401  298.300  296.853  90.00  90.00  90.00 I 2 2 2
ATOM  16017  N   ARG C 356       7.575  71.661  23.558  1.00 39.83           N
ATOM  16018  CA  ARG C 356       8.927  72.211  23.601  1.00 40.00           C
ATOM  16019  C   ARG C 356       8.839  73.723  23.795  1.00 40.77           C
ATOM  16020  O   ARG C 356       7.801  74.331  23.528  1.00 39.79           O
ATOM  16021  CB  ARG C 356       9.680  71.902  22.301  1.00 39.66           C
ATOM  16022  CG  ARG C 356       9.591  70.445  21.843  1.00 40.99           C
ATOM  16023  CD  ARG C 356      10.028  69.472  22.933  1.00 39.51           C
ATOM  16024  NE  ARG C 356       9.788  68.083  22.549  1.00 39.83           N
ATOM  16025  CZ  ARG C 356      10.465  67.436  21.604  1.00 38.82           C
ATOM  16026  NH1 ARG C 356      11.437  68.047  20.941  1.00 38.06           N
ATOM  16027  NH2 ARG C 356      10.156  66.180  21.311  1.00 37.75           N
ATOM  21379  N   LEU C1038      16.054  63.328  26.780  1.00 33.34           N
ATOM  21380  CA  LEU C1038      15.414  64.531  26.255  1.00 33.21           C
ATOM  21381  C   LEU C1038      15.049  65.503  27.377  1.00 33.24           C
ATOM  21382  O   LEU C1038      15.339  66.695  27.293  1.00 34.52           O
ATOM  21383  CB  LEU C1038      14.144  64.155  25.486  1.00 32.34           C
ATOM  21384  CG  LEU C1038      13.225  65.300  25.045  1.00 31.54           C
ATOM  21385  CD1 LEU C1038      13.949  66.200  24.057  1.00 30.59           C
ATOM  21386  CD2 LEU C1038      11.965  64.726  24.414  1.00 31.06           C
ATOM  21455  N   ASP C1047       8.324  59.914  23.481  1.00 37.88           N
ATOM  21456  CA  ASP C1047       7.566  60.655  22.481  1.00 37.10           C
ATOM  21457  C   ASP C1047       8.089  60.607  21.054  1.00 36.68           C
ATOM  21458  O   ASP C1047       7.603  61.340  20.197  1.00 38.13           O
ATOM  21459  CB  ASP C1047       7.430  62.116  22.922  1.00 38.09           C
ATOM  21460  CG  ASP C1047       8.663  62.942  22.604  1.00 39.82           C
ATOM  21461  OD1 ASP C1047       9.788  62.412  22.731  1.00 39.75           O
ATOM  21462  OD2 ASP C1047       8.507  64.130  22.238  1.00 39.30           O
ATOM  22656  N   TYR C1201      17.944  61.058  18.108  1.00 31.93           N
ATOM  22657  CA  TYR C1201      16.950  61.625  19.015  1.00 30.51           C
ATOM  22658  C   TYR C1201      17.359  61.605  20.491  1.00 28.70           C
ATOM  22659  O   TYR C1201      16.768  60.886  21.299  1.00 28.44           O
ATOM  22660  CB  TYR C1201      15.619  60.886  18.853  1.00 31.95           C
ATOM  22661  CG  TYR C1201      14.418  61.616  19.425  1.00 32.23           C
ATOM  22662  CD1 TYR C1201      13.967  62.810  18.857  1.00 33.61           C
ATOM  22663  CD2 TYR C1201      13.718  61.102  20.518  1.00 33.06           C
ATOM  22664  CE1 TYR C1201      12.846  63.473  19.361  1.00 33.76           C
ATOM  22665  CE2 TYR C1201      12.596  61.756  21.032  1.00 33.12           C
ATOM  22666  CZ  TYR C1201      12.166  62.939  20.446  1.00 34.19           C
ATOM  22667  OH  TYR C1201      11.052  63.581  20.937  1.00 34.40           O
TER
HETATM28714  O   HOH C8320       9.280  63.671  18.743  1.00 59.06           O
END
"""

pdb_1r27_tyr415 = """\
CRYST1  127.401  298.300  296.853  90.00  90.00  90.00 I 2 2 2
ATOM  26327  N   TYR D 415      14.965  40.320  64.778  1.00 27.53           N
ATOM  26328  CA  TYR D 415      15.673  41.262  65.631  1.00 27.06           C
ATOM  26329  C   TYR D 415      16.469  42.274  64.797  1.00 28.08           C
ATOM  26330  O   TYR D 415      16.404  43.484  65.040  1.00 28.08           O
ATOM  26331  CB  TYR D 415      16.596  40.474  66.563  1.00 26.98           C
ATOM  26332  CG  TYR D 415      17.589  41.298  67.342  1.00 29.66           C
ATOM  26333  CD1 TYR D 415      17.170  42.218  68.306  1.00 29.06           C
ATOM  26334  CD2 TYR D 415      18.958  41.142  67.132  1.00 29.85           C
ATOM  26335  CE1 TYR D 415      18.095  42.960  69.043  1.00 30.10           C
ATOM  26336  CE2 TYR D 415      19.891  41.878  67.861  1.00 32.39           C
ATOM  26337  CZ  TYR D 415      19.453  42.785  68.815  1.00 31.26           C
ATOM  26338  OH  TYR D 415      20.379  43.509  69.530  1.00 33.01           O
ATOM  26498  N   GLU D 437      20.460  39.379  73.511  1.00 32.61           N
ATOM  26499  CA  GLU D 437      20.860  40.274  72.434  1.00 32.89           C
ATOM  26500  C   GLU D 437      22.306  40.052  72.036  1.00 32.64           C
ATOM  26501  O   GLU D 437      22.642  40.083  70.852  1.00 31.37           O
ATOM  26502  CB  GLU D 437      20.702  41.735  72.840  1.00 34.39           C
ATOM  26503  CG  GLU D 437      19.310  42.293  72.697  1.00 36.51           C
ATOM  26504  CD  GLU D 437      19.301  43.808  72.811  1.00 38.76           C
ATOM  26505  OE1 GLU D 437      19.730  44.491  71.847  1.00 38.38           O
ATOM  26506  OE2 GLU D 437      18.879  44.312  73.874  1.00 38.44           O
TER
HETATM29064  O   HOH D7887      18.656  46.536  68.336  1.00 40.54           O
HETATM29153  O   HOH D7976      22.116  44.147  67.400  1.00 61.89           O
END
"""

if __name__ == "__main__":
  test_candidates()
  test_7prg_fuc_o3_neutron()
  test_7prg_fuc_o3_xray()
  test_tyr_short_contact()
  test_safeguard_singleton()
  test_fallback_state_cache()
  print("OK")
