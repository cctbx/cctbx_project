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
import os, math, random
import numpy as np
import iotbx.pdb
import libtbx.load_env
import mmtbx.f_model
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.math import dihedral_angle
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

if(__name__ == "__main__"):
  exercise_priors()
  exercise_side_chain()
  exercise_evidence_and_decisions()
  exercise_stubs()
  exercise_local_refinement()
  print("OK")
