"""
Tests of mmtbx.refinement.rotamer_bayes:
  priors: the inverted Top8000 masses sum to 1 and give the known
    frequencies;
  side chains: chi angles are set as asked;
  evidence: the exact dLL of a move matches the target evaluated on the
    F_model of the moved structure; the prediction agrees with it;
  decisions: on data from a model with some side chains in other wells,
    those are changed to the right wells and the others are not changed.
"""
from __future__ import absolute_import, division, print_function
import math, os, random
import numpy as np
import iotbx.pdb
import libtbx.load_env
import mmtbx.f_model
from cctbx.array_family import flex
from scitbx.math import dihedral_angle
from libtbx import group_args
from libtbx.utils import null_out
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.refinement.rotamer_bayes as rb
from mmtbx.refinement.rotamer_bayes import priors, side_chain

def exercise_priors():
  p = priors.rotamer_priors()
  for aa in ("leu", "val", "ser", "phe", "lys", "arg", "glu", "asn", "trp"):
    w = p.wells(aa)
    total = sum([x.mass for x in w.values()])
    assert abs(total - 1) < 1.e-3, (aa, total)
    assert abs(w[priors.OUTLIER].mass - priors.OUTLIER_THRESHOLD) < 1.e-4
  leu = p.wells("leu")
  assert 0.60 < leu["mt"].mass < 0.67, leu["mt"]
  assert 0.29 < leu["tp"].mass < 0.35, leu["tp"]
  val = p.wells("val")
  assert 0.70 < val["t"].mass < 0.76, val["t"]
  for aa in ("leu", "phe", "lys", "asn"):
    for w in p.wells(aa).values():
      if(w.name in (priors.OUTLIER, priors.UNCLASSIFIED)): continue
      if(w.mass < 0.01): continue
      assert p.classify(aa, w.chis) == w.name, (aa, w.name, w.chis)

def residue_files():
  d = libtbx.env.find_in_repositories(os.path.join("mmtbx",
    "idealized_aa_residues", "data"))
  return d

def exercise_side_chain():
  sa = SidechainAngles(False)
  for aa in ("leu", "arg", "phe", "lys", "met", "asn", "ile", "thr"):
    h = iotbx.pdb.input(os.path.join(residue_files(), aa + ".pdb")
      ).construct_hierarchy()
    ag = h.only_model().chains()[0].residue_groups()[0].atom_groups()[0]
    sc = side_chain.side_chain(ag, sa)
    assert sc.ok, (aa, sc.reason)
    start = sc.chis()
    target = [c + 47. * (k + 1) for k, c in enumerate(start)]
    xyz = sc.set_chis(target)
    got = sc.chis(xyz)
    for a, b in zip(target, got):
      assert abs((a - b + 180) % 360 - 180) < 1.e-3, (aa, target, got)
    # bond lengths unchanged
    d0 = np.sqrt(((sc.xyz[:, None] - sc.xyz[None]) ** 2).sum(axis=2))
    d1 = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(axis=2))
    assert np.abs((d0 - d1)[d0 < 1.7]).max() < 1.e-6
    # same as scitbx's dihedral
    q = sc.chi_atoms[0]
    assert abs(dihedral_angle(sites=[tuple(x) for x in xyz[q]], deg=True)
      - got[0]) < 1.e-6

# residue, (truth well, model well)
LAYOUT = [("leu", "tp", "mt"), ("phe", "m-80", "m-80"), ("lys", "mttt",
  "mttt"), ("arg", "mtt180", "ttt180"), ("val", "t", "t"), ("met", "mmm",
  "mtm"), ("leu", "mt", "mt"), ("glu", "mt-10", "mt-10"), ("ile", "mt",
  "mt"), ("tyr", "t80", "m-80"), ("asn", "m-40", "m-40"), ("lys", "tttt",
  "tttt")]

def build_pair(d_min=2.0, seed=1):
  """ truth and model hierarchies in a P1 box, residues in the wells of
  LAYOUT """
  sa = SidechainAngles(False)
  p = priors.rotamer_priors()
  lines = ["CRYST1   45.000   45.000   45.000  90.00  90.00  90.00 P 1"]
  pdb_t, pdb_m = [], []
  n = 0
  serial = 1
  for k, (aa, wt, wm) in enumerate(LAYOUT):
    h = iotbx.pdb.input(os.path.join(residue_files(), aa + ".pdb")
      ).construct_hierarchy()
    ag = h.only_model().chains()[0].residue_groups()[0].atom_groups()[0]
    sc = side_chain.side_chain(ag, sa)
    shift = np.array([(k % 3) * 14. + 2, ((k // 3) % 3) * 14. + 2,
      (k // 9) * 14. + 2])
    for well_name, out in ((wt, pdb_t), (wm, pdb_m)):
      xyz = sc.set_chis(p.wells(aa)[well_name].chis) + shift
      for a, x in zip(ag.atoms(), xyz):
        out.append("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 20.00"
          "          %2s" % (serial + len(out), a.name, aa.upper(),
          chr(ord("A") + k), 1, x[0], x[1], x[2], a.element.strip()))
  t = iotbx.pdb.input(source_info=None, lines=lines + pdb_t)
  m = iotbx.pdb.input(source_info=None, lines=lines + pdb_m)
  return t, m

def make_fmodel(truth_inp, model_inp, d_min=2.0, seed=1):
  random.seed(seed); flex.set_random_seed(seed)
  xt = truth_inp.xray_structure_simple()
  xm = model_inp.xray_structure_simple()
  f = abs(xt.structure_factors(d_min=d_min, algorithm="direct").f_calc())
  noise = flex.double([random.gauss(0, 1) for i in range(f.size())])
  f = f.customized_copy(data=f.data() * (1 + 0.05 * noise),
    sigmas=f.data() * 0.05)
  f = f.select(f.data() > 0)
  flags = f.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(f_obs=f, r_free_flags=flags,
    xray_structure=xm, target_name="ml")
  fmodel.update_all_scales(remove_outliers=False)
  return fmodel

def exercise_evidence_and_decisions():
  t, m = build_pair()
  fmodel = make_fmodel(t, m)
  hm = m.construct_hierarchy()
  params = rb.master_params().extract()
  params.exact = True
  r = rb.run(fmodel=fmodel, pdb_hierarchy=hm, params=params, log=null_out())
  r.show(out=null_out())
  ev = r.ev
  # exact dLL against the target on the f_model of the moved structure
  rec = r.records[0]
  c = rec.best
  xrs = fmodel.xray_structure.deep_copy_scatterers()
  sites = xrs.sites_cart()
  for i, x in zip(rec.sc.i_seqs, c.xyz):
    sites[int(i)] = tuple(x)
  xrs.set_sites_cart(sites)
  fm2 = fmodel.deep_copy()
  fm2.update_xray_structure(xray_structure=xrs, update_f_calc=True,
    update_f_mask=True, force_update_f_mask=True)
  t2 = ev.functor.core(f_calc=fm2.f_model(), compute_gradients=False)
  brute = -ev.n_work * (t2.target_work() - ev.t0)
  # f_model's F_calc is by FFT, dF by direct summation (~0.2%), and the
  # mask change is recomputed locally
  assert abs(brute - c.exact) < 1.e-2 * max(1., abs(brute)), (brute, c.exact)
  # prediction against exact
  pe = [(cd.predicted, cd.exact) for rr in r.records for cd in rr.candidates
    if cd.exact is not None]
  pr, ex = np.array(pe).T
  cc = np.corrcoef(pr, ex)[0, 1]
  assert cc > 0.9, cc
  # decisions
  by_chain = dict((rr.label.split()[0], rr) for rr in r.records)
  for k, (aa, wt, wm) in enumerate(LAYOUT):
    rr = by_chain[chr(ord("A") + k)]
    if(wt == wm):
      assert not rr.accept, (rr.label, rr.best.well, rr.best.log_odds)
    else:
      assert rr.accept and rr.best.well == wt, (rr.label, rr.best.well,
        rr.best.log_odds)
  n = r.apply(hm, fmodel.xray_structure)
  assert n == len([1 for x in LAYOUT if x[1] != x[2]])
  return r

def exercise_stubs():
  """ the wrong side chains made into stubs (zero occupancy beyond CB) are
  rebuilt in the right wells """
  t, m = build_pair()
  hm = m.construct_hierarchy()
  wrong = [chr(ord("A") + k) for k, x in enumerate(LAYOUT) if x[1] != x[2]]
  for ch in hm.only_model().chains():
    if(ch.id not in wrong): continue
    for a in ch.atoms():
      if(a.name.strip() not in ("N", "CA", "C", "O", "CB")): a.set_occ(0.)
  m2 = iotbx.pdb.input(source_info=None, lines=hm.as_pdb_string(
    crystal_symmetry=m.crystal_symmetry()))
  fmodel = make_fmodel(t, m2)
  h2 = m2.construct_hierarchy()
  r = rb.run(fmodel=fmodel, pdb_hierarchy=h2, log=null_out())
  stubs = [rr for rr in r.records if rr.stub]
  assert len(stubs) == len(wrong)
  for rr in stubs:
    k = ord(rr.label.split()[0]) - ord("A")
    assert rr.accept and rr.best.well == LAYOUT[k][1], (rr.label,
      rr.best.well, rr.posterior_best)
  n = r.apply(h2, fmodel.xray_structure)
  assert n == len(wrong)
  occ = fmodel.xray_structure.scatterers().extract_occupancies()
  assert occ.all_eq(1.)

def exercise_local_refinement():
  """ gradients of the local refinement target (data: sites and the B
  shift; chi well restraint) against finite differences """
  from mmtbx.refinement.rotamer_bayes import refine
  t, m = build_pair()
  fmodel = make_fmodel(t, m)
  r = rb.run(fmodel=fmodel, pdb_hierarchy=m.construct_hierarchy(),
    log=null_out())
  rec = r.records[0]; sc = rec.sc; mv = sc.moving
  lr = refine.local_refinement(r.ev, None,
    fmodel.xray_structure.sites_cart(), sc.i_seqs[mv])
  x0 = rec.candidates[0].xyz[mv].copy(); db = 3.
  t0, gx, gb = lr.data(x0, db)
  h = 1.e-4
  for i, j in ((0, 0), (1, 2), (2, 1)):
    xp = x0.copy(); xp[i, j] += h; xm = x0.copy(); xm[i, j] -= h
    fd = (lr.data(xp, db, False)[0] - lr.data(xm, db, False)[0]) / (2 * h)
    assert abs(fd - gx[i, j]) < 1.e-3 * max(1., abs(fd)), (fd, gx[i, j])
  fd = (lr.data(x0, db + 1.e-3, False)[0] - lr.data(x0, db - 1.e-3,
    False)[0]) / 2.e-3
  assert abs(fd - gb) < 1.e-3 * max(1., abs(fd)), (fd, gb)
  # well restraint: centres 40 degrees away, so it is active
  lr.set_well(sc.xyz, mv, sc.chi_atoms, [c + 40. for c in sc.chis(
    sc.xyz)])
  e, g = lr.well_penalty(sc.xyz[mv])
  assert e > 0
  xs = sc.xyz[mv].copy()
  for i, j in ((0, 0), (1, 1)):
    xp = xs.copy(); xp[i, j] += h; xm = xs.copy(); xm[i, j] -= h
    fd = (lr.well_penalty(xp)[0] - lr.well_penalty(xm)[0]) / (2 * h)
    assert abs(fd - g[i, j]) < 1.e-3 * max(1., abs(fd)), (fd, g[i, j])

def exercise_unsupported_targets():
  """ least-squares targets are refused (no likelihood scale) """
  from libtbx.test_utils import Exception_expected
  t, m = build_pair()
  fmodel = make_fmodel(t, m)
  fmodel.set_target_name("ls_wunit_k1")
  try:
    rb.run(fmodel=fmodel, pdb_hierarchy=m.construct_hierarchy(),
      log=null_out())
  except Exception as e:
    assert "least squares" in str(e), str(e)
  else: raise Exception_expected

class _candidate(object):
  def __init__(self, well, dll):
    self.well = well
    self.clash = False
    self.predicted = dll
  def dll(self, exact):
    return self.predicted

class _record(object):
  """ a Leu in well mt, with the given evidence for the other wells """
  def __init__(self, dlls, wells):
    self.sc = group_args(resname="LEU")
    self.current_well = "mt"
    self.current_prior = wells["mt"].mass
    self.candidates = [_candidate(w, dlls.get(w, -50.)) for w in wells
      if w not in ("mt", priors.OUTLIER, priors.UNCLASSIFIED)]

def exercise_ncs_pooling():
  """ ncs.pooled with given per-copy evidence. Copies that agree weakly
  (dLL 0.5 each for tp, too little alone against the prior): rho is high
  and the pooled odds favour tp. Copies that disagree strongly (two for
  tp, two for staying in mt): rho is low and each keeps its own choice. """
  from mmtbx.refinement.rotamer_bayes import ncs
  p = priors.rotamer_priors()
  wells = p.wells("leu")
  own_odds = 0.5 + math.log(wells["tp"].mass / wells["mt"].mass)
  assert own_odds < 0
  agree = [[_record({"tp": 0.5}, wells) for i in range(4)] for g in range(2)]
  pool = ncs.pooled(agree, p, exact=False)
  assert pool.rho > 0.9, pool.rho
  for rec, lp, r_shared, others, e in pool.posteriors():
    assert max(lp, key=lp.get) == "tp", lp
    assert lp["tp"] - lp["mt"] > own_odds + 1, (lp["tp"] - lp["mt"], own_odds)
  disagree = [[_record({"tp": 20.}, wells), _record({"tp": 20.}, wells),
    _record({"tp": -20.}, wells), _record({"tp": -20.}, wells)]
    for g in range(2)]
  pool = ncs.pooled(disagree, p, exact=False)
  assert pool.rho < 0.1, pool.rho
  for k, (rec, lp, r_shared, others, e) in enumerate(pool.posteriors()):
    best = max(lp, key=lp.get)
    assert best == ("tp" if k % 4 < 2 else "mt"), (k, best)

TRIPEPTIDE = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM    440  N   LEU A  59      -7.510 122.749  38.503  1.00 89.42           N
ATOM    441  CA  LEU A  59      -8.480 121.839  37.894  1.00 90.13           C
ATOM    442  C   LEU A  59      -9.914 122.306  38.123  1.00 90.59           C
ATOM    443  O   LEU A  59     -10.778 122.134  37.254  1.00 88.28           O
ATOM    444  CB  LEU A  59      -8.304 120.428  38.434  1.00 82.68           C
ATOM    445  CG  LEU A  59      -9.142 119.474  37.595  1.00 86.20           C
ATOM    446  CD1 LEU A  59      -8.562 119.366  36.180  1.00 93.60           C
ATOM    447  CD2 LEU A  59      -9.250 118.118  38.258  1.00 93.32           C
ATOM    448  N   PHE A  60     -10.187 122.893  39.292  1.00 89.57           N
ATOM    449  CA  PHE A  60     -11.493 123.501  39.534  1.00 89.43           C
ATOM    450  C   PHE A  60     -11.793 124.577  38.497  1.00 88.98           C
ATOM    451  O   PHE A  60     -12.907 124.648  37.970  1.00 87.16           O
ATOM    452  CB  PHE A  60     -11.538 124.079  40.950  1.00 86.28           C
ATOM    453  CG  PHE A  60     -12.899 124.556  41.382  1.00 85.30           C
ATOM    454  CD1 PHE A  60     -13.947 123.661  41.533  1.00 88.93           C
ATOM    455  CD2 PHE A  60     -13.120 125.892  41.665  1.00 81.10           C
ATOM    456  CE1 PHE A  60     -15.193 124.092  41.937  1.00 81.71           C
ATOM    457  CE2 PHE A  60     -14.358 126.329  42.074  1.00 80.56           C
ATOM    458  CZ  PHE A  60     -15.397 125.427  42.209  1.00 83.74           C
ATOM    459  N   GLU A  61     -10.801 125.413  38.183  1.00 87.78           N
ATOM    460  CA  GLU A  61     -10.971 126.431  37.152  1.00 89.55           C
ATOM    461  C   GLU A  61     -11.309 125.807  35.805  1.00 93.53           C
ATOM    462  O   GLU A  61     -12.214 126.275  35.104  1.00 93.94           O
ATOM    463  CB  GLU A  61      -9.703 127.270  37.029  1.00 89.15           C
ATOM    464  CG  GLU A  61      -9.646 128.501  37.893  1.00 81.08           C
ATOM    465  CD  GLU A  61      -8.344 129.250  37.693  1.00 92.28           C
ATOM    466  OE1 GLU A  61      -7.298 128.739  38.156  1.00 87.99           O
ATOM    467  OE2 GLU A  61      -8.362 130.324  37.045  1.00 95.81           O
"""

def exercise_backrub():
  """ side_chain.backrub keeps the flanking C-alphas and the bond lengths,
  turns CA-CB by about the angle and moves the carbonyl oxygens less than
  CA; a run with backrub on finds the residue's neighbours """
  inp = iotbx.pdb.input(source_info=None, lines=TRIPEPTIDE)
  h = inp.construct_hierarchy()
  rgs = list(h.only_model().chains()[0].residue_groups())
  prev, res, nxt = [dict((a.name.strip(), np.array(a.xyz)) for a in
    rg.atoms()) for rg in rgs]
  for theta in (10., -20.):
    p1, r1, n1 = side_chain.backrub(prev, res, nxt, theta)
    for a, b in (("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")):
      assert abs(np.linalg.norm(r1[a] - r1[b]) -
        np.linalg.norm(res[a] - res[b])) < 1.e-6
    assert abs(np.linalg.norm(p1["C"] - r1["N"]) -
      np.linalg.norm(prev["C"] - res["N"])) < 1.e-6
    assert abs(np.linalg.norm(r1["C"] - n1["N"]) -
      np.linalg.norm(res["C"] - nxt["N"])) < 1.e-6
    v0 = res["CB"] - res["CA"]; v1 = r1["CB"] - r1["CA"]
    turn = math.degrees(math.acos(v0.dot(v1) / np.linalg.norm(v0) /
      np.linalg.norm(v1)))
    assert abs(turn - abs(theta)) < 2., (theta, turn)
    ca = np.linalg.norm(r1["CA"] - res["CA"])
    for o in (np.linalg.norm(p1["O"] - prev["O"]),
              np.linalg.norm(r1["O"] - res["O"])):
      assert o < ca, (o, ca)
  fmodel = make_fmodel(inp, inp)
  p = rb.master_params().extract()
  p.backrub = True
  r = rb.run(fmodel=fmodel, pdb_hierarchy=h, params=p, log=null_out())
  assert list(r._neighbours) == [("A", "60", "")], list(r._neighbours)

if(__name__ == "__main__"):
  exercise_priors()
  exercise_side_chain()
  exercise_evidence_and_decisions()
  exercise_stubs()
  exercise_local_refinement()
  exercise_unsupported_targets()
  exercise_ncs_pooling()
  exercise_backrub()
  print("OK")
