from __future__ import absolute_import, division, print_function
""" mmtbx.refinement.llgi_tncs: tNCS groups from a model, the pair
modulation A(h), renormalised data, and the rho/tau fits, on a synthetic
structure with a known tNCS relationship. """
import math, random
import iotbx.pdb
from cctbx import crystal, miller, sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx import matrix
import mmtbx.refinement.llgi_tncs as llgi_tncs
import mmtbx.refinement.llgi_hybrid as llgi_hybrid

CS = crystal.symmetry(unit_cell=(62, 71, 83, 90, 90, 90),
  space_group_symbol="P 21 21 21")

def build_model(n_res=60, seed=0):
  """ chain A (poly-ALA, random sites in a 12 A sphere); B = A rotated by
  3 degrees and translated, then moved by a symmetry operator so that the
  tNCS relationship needs it; C a different molecule; D = A rotated by
  60 degrees (not tNCS). """
  rnd = random.Random(seed)
  def random_sites(n):
    out = []
    while len(out) < n:
      v = [rnd.uniform(-12, 12) for k in range(3)]
      if sum(x * x for x in v) < 144: out.append(v)
    return out
  names = ["N", "CA", "C", "O", "CB"]
  a_sites = random_sites(5 * n_res)
  r3 = matrix.col((1, 1, 0.3)).normalize().axis_and_angle_as_r3_rotation_matrix(3, deg=True)
  r60 = matrix.col((0, 1, 1)).normalize().axis_and_angle_as_r3_rotation_matrix(60, deg=True)
  uc = CS.unit_cell()
  centre_a = matrix.col(uc.orthogonalize((0.25, 0.25, 0.25)))
  t_b = matrix.col(uc.orthogonalize((0.02, 0.49, 0.03)))
  op = sgtbx.rt_mx("-x+1/2,-y,z+1/2")
  inv = op.inverse()
  b_cart = []
  for s in a_sites:
    x = r3 * matrix.col(s) + centre_a + t_b
    f = inv * uc.fractionalize(tuple(x))
    b_cart.append(uc.orthogonalize(f))
  lines = []
  serial = [0]
  def add(chain, resname, i_res, name, xyz, b):
    serial[0] += 1
    el = name[0]
    lines.append("ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f           %s" % (
      serial[0], name, resname, chain, i_res, xyz[0], xyz[1], xyz[2], b, el))
  for k, s in enumerate(a_sites):
    add("A", "ALA", k // 5 + 1, names[k % 5], tuple(matrix.col(s) + centre_a), 20 + (k % 7))
  for k, s in enumerate(b_cart):
    add("B", "ALA", k // 5 + 1, names[k % 5], s, 22 + (k % 5))
  centre_c = matrix.col(uc.orthogonalize((0.7, 0.2, 0.6)))
  for k, s in enumerate(random_sites(4 * n_res)):
    add("C", "GLY", k // 4 + 1, ["N", "CA", "C", "O"][k % 4], tuple(matrix.col(s) + centre_c), 25.)
  centre_d = matrix.col(uc.orthogonalize((0.6, 0.7, 0.15)))
  for k, s in enumerate(a_sites):
    add("D", "ALA", k // 5 + 1, names[k % 5], tuple(r60 * matrix.col(s) + centre_d), 20.)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=flex.split_lines("\n".join(lines)))
  h = pdb_inp.construct_hierarchy()
  return h, h.extract_xray_structure(crystal_symmetry=CS)

def exercise_groups(h):
  groups = llgi_tncs.tncs_groups(h, CS)
  assert len(groups) == 1
  g = groups.groups[0]
  assert g.chain_ids == ["A", "B"]
  assert g.i_seqs[0].size() == 300 and g.i_seqs[1].size() == 300
  assert approx_equal(g.angles[1], 3.0, eps=0.01)
  return groups

def direct_pair_term(ms, sites_m, sites_n, z_sq, b_sum):
  """ (1/Nsym) sum_k sum_j w_j exp(2 pi i (T_k^T h).d_j) / sum_j w_j over
  all the space group's operators (inversion and centring included) """
  d = sites_m - sites_n
  out = flex.double()
  rots = [op.r().as_double() for op in ms.space_group().all_ops()]
  for h, ss in zip(ms.indices(), ms.d_star_sq().data()):
    w = z_sq * flex.exp(-b_sum * (ss / 4.))
    tot = 0
    for r in rots:
      hr = matrix.row(h) * matrix.sqr(r)
      phase = d.dot(tuple(hr)) * (2 * math.pi)
      tot += flex.sum(w * flex.cos(phase))
    out.append(tot / (len(rots) * flex.sum(w)))
  return out

def direct_modulation(ms, groups, xrs):
  """ A(h) as the normalised model-model block: 2 Re sum_{m<n} (1/Nsym)
  sum_k sum_j Z_j^2 exp(-(B_jm + B_jn) s^2/4) exp(2 pi i (T_k^T h).d_j)
  over sum_all Z^2 exp(-B s^2/2) (here every chain atom of a copy is a
  common atom, so no scaling by the copy fractions is needed) """
  z_all = llgi_tncs._atomic_weight_sq(xrs)
  b_all = xrs.extract_u_iso_or_u_equiv() * (8 * math.pi ** 2)
  sf = xrs.sites_frac()
  out = flex.double(ms.size(), 0)
  for g in groups.groups:
    placed = [llgi_tncs.apply_op(op, sf.select(sel)) for op, sel in zip(g.ops, g.i_seqs)]
    z = z_all.select(g.i_seqs[0])
    for m in range(len(placed)):
      for n in range(m + 1, len(placed)):
        b_sum = b_all.select(g.i_seqs[m]) + b_all.select(g.i_seqs[n])
        p = direct_pair_term(ms, placed[m], placed[n], z, b_sum)
        for i, ss in enumerate(ms.d_star_sq().data()):
          out[i] += 2 * p[i] * flex.sum(z * flex.exp(-b_sum * (ss / 4.)))
  for i, ss in enumerate(ms.d_star_sq().data()):
    out[i] /= flex.sum(z_all * flex.exp(-b_all * (ss / 2.)))
  return out

def exercise_pair_term(groups, xrs):
  ms = miller.build_set(CS, anomalous_flag=False, d_min=6.0)
  g = groups.groups[0]
  sf = xrs.sites_frac()
  placed = [llgi_tncs.apply_op(op, sf.select(sel)) for op, sel in zip(g.ops, g.i_seqs)]
  z_sq = llgi_tncs._atomic_weight_sq(xrs).select(g.i_seqs[0])
  b = xrs.extract_u_iso_or_u_equiv() * (8 * math.pi ** 2)
  b_sum = b.select(g.i_seqs[0]) + b.select(g.i_seqs[1])
  p_fft = llgi_tncs.pair_term(ms, placed[0], placed[1], z_sq, b_sum)
  p_dir = direct_pair_term(ms, placed[0], placed[1], z_sq, b_sum)
  err = flex.max(flex.abs(p_fft - p_dir))
  assert err < 0.01, err
  assert flex.max(flex.abs(p_dir)) > 0.5  # the copies do interfere
  mod = llgi_tncs.pair_modulation(ms, groups, xrs)
  f = mod.fractions[0]
  assert approx_equal(f[0], f[1], eps=1e-6)  # both 300 C/N/O atoms
  assert mod.f_total >= sum(f) - 1e-12
  assert flex.max(flex.abs(mod.a - direct_modulation(ms, groups, xrs))) < 0.01
  # copy B with B-factors 60 A^2 higher: it weighs less at high resolution
  # and the pair modulation falls off as sech((B_B - B_A) s^2/4)
  x2 = xrs.deep_copy_scatterers()
  b = x2.extract_u_iso_or_u_equiv()
  b.set_selected(g.i_seqs[1], b.select(g.i_seqs[1]) + 60 / (8 * math.pi ** 2))
  x2.set_u_iso(values=b)
  mod2 = llgi_tncs.pair_modulation(ms, groups, x2)
  a2 = direct_modulation(ms, groups, x2)
  assert flex.max(flex.abs(mod2.a - a2)) < 0.01
  assert flex.max(flex.abs(mod2.a)) < flex.max(flex.abs(mod.a))
  # the consensus of two copies related exactly by the operator is the
  # copies themselves; with the second copy perturbed, the consensus pair
  # term is the stronger one (the perturbation is not in the consensus)
  cons = llgi_tncs.pair_modulation(ms, groups, xrs, basis="consensus")
  assert flex.max(flex.abs(cons.a - mod.a)) < 0.01
  x3 = xrs.deep_copy_scatterers()
  rnd = random.Random(2)
  sites = x3.sites_cart()
  for i in g.i_seqs[1]:
    sites[i] = tuple(v + rnd.gauss(0, 0.6) for v in sites[i])
  x3.set_sites_cart(sites)
  noisy = llgi_tncs.pair_modulation(ms, groups, x3)
  noisy_cons = llgi_tncs.pair_modulation(ms, groups, x3, basis="consensus")
  hi = ms.d_spacings().data() < 7.
  rms = lambda x: math.sqrt(flex.mean(flex.pow2(x.select(hi))))
  assert rms(noisy_cons.a) > rms(noisy.a) * 1.05, (rms(noisy_cons.a), rms(noisy.a))

def exercise_pair_term_centrosymmetric():
  """ P-1: the pair term counts the inversion (Re is the same for R and -R) """
  cs = crystal.symmetry(unit_cell=(40, 45, 50, 80, 95, 100),
    space_group_symbol="P -1")
  rnd = random.Random(4)
  n = 80
  m = flex.vec3_double([(rnd.uniform(0.3, 0.5), rnd.uniform(0.3, 0.5),
    rnd.uniform(0.3, 0.5)) for i in range(n)])
  copy = m + flex.vec3_double(n, (0.03, 0.5, 0.02)) + flex.vec3_double([
    tuple(rnd.gauss(0, 0.004) for k in range(3)) for i in range(n)])
  z_sq = flex.double(n, 64.); b_sum = flex.double(n, 40.)
  ms = miller.build_set(cs, anomalous_flag=False, d_min=5.0)
  p_fft = llgi_tncs.pair_term(ms, copy, m, z_sq, b_sum)
  p_dir = direct_pair_term(ms, copy, m, z_sq, b_sum)
  slope = flex.sum(p_fft * p_dir) / flex.sum(p_dir * p_dir)
  rms = math.sqrt(flex.mean(flex.pow2(p_fft - p_dir)))
  assert abs(slope - 1) < 0.02 and rms < 0.01, (slope, rms)  # not 0.5
  assert flex.max(flex.abs(p_dir)) > 0.5

def exercise_renormalise():
  n = 50
  flex.set_random_seed(3)
  i_obs = flex.random_double(n) * 4 - 0.3
  sig = flex.random_double(n) * 0.5 + 0.05
  resn = flex.random_double(n) + 0.5
  # E_obs^2 = 20, sigma 5: no Rice solution (the posterior is more
  # dispersed than any Rice distribution with the same <E^2>)
  i_obs[0], sig[0] = 20 * resn[0] ** 2, 5 * resn[0] ** 2
  cen = flex.bool(n, False)
  t = flex.double(n, 1.)
  d = llgi_tncs.renormalise(i_obs, sig, resn, t, cen)
  e2, s2 = llgi_hybrid.e_obs_sq_and_sigma(i_obs=i_obs, sig_i_obs=sig, resn=resn, teps=t)
  assert approx_equal(d.e_obs_sq, e2) and approx_equal(d.sig_e_obs_sq, s2)
  assert d.force_exact[0] and d.dobs[0] == llgi_hybrid.RICE_DOBS_NONE
  assert approx_equal(d.feff, d.eeff * resn)
  t2 = flex.random_double(n) + 0.3
  d2 = llgi_tncs.renormalise(i_obs, sig, resn, t2, cen)
  assert approx_equal(d2.e_obs_sq, i_obs / (resn * resn * t2))

def synthetic_e_data(a, rho_or_tau, sigmaa=None, seed=5, sig_frac=0.05):
  """ acentric E with Wilson variance t = 1 + rho A (sigmaa None), or
  E = sigmaa Ec + error with <|Ec|^2> = 1 + A and error variance
  (1 - sigmaa^2)(1 + tau A); intensities with Gaussian error. """
  rnd = random.Random(seed)
  e_obs_sq, sig, e_model = flex.double(), flex.double(), flex.double()
  for ai in a:
    def cn(var):
      s = math.sqrt(var / 2)
      return complex(rnd.gauss(0, s), rnd.gauss(0, s))
    if sigmaa is None:
      expected = 1 + rho_or_tau * ai
      e = cn(expected); ec = 0
    else:
      ec = cn(1 + ai)
      v = (1 - sigmaa ** 2) * (1 + rho_or_tau * ai)
      e = sigmaa * ec + cn(v)
      expected = abs(sigmaa * ec) ** 2 + v
    i_true = abs(e) ** 2
    # sigma from the expected, not the realised, intensity: otherwise a
    # known sigma carries information about the realised intensity
    s = sig_frac * (1 + expected)
    e_obs_sq.append(i_true + rnd.gauss(0, s)); sig.append(s); e_model.append(abs(ec))
  return e_obs_sq, sig, e_model

def exercise_fits(groups, xrs):
  ms = miller.build_set(CS, anomalous_flag=False, d_min=2.8).select_acentric()
  mod = llgi_tncs.pair_modulation(ms, groups, xrs)
  a = mod.a
  isel = flex.size_t(sorted(random.Random(1).sample(range(a.size()), 6000)))
  a = a.select(isel)
  dss = ms.d_star_sq().data().select(isel)
  n = a.size()
  cen = flex.bool(n, False)
  one = flex.double(n, 1.)
  # data-only: rho
  e2, sig, _ = synthetic_e_data(a, 0.8 / mod.f_total)
  rho = llgi_tncs.fit_wilson_modulation(e2, sig, cen, a, dss, mod.f_total)
  assert approx_equal(rho.p0, 0.8 / mod.f_total, eps=0.1), rho.p0
  # tau, with sigmaA fixed: constant (default) and linear in d*^2. With
  # 6000 reflections, sigmaA 0.5 and |A| <= 0.5, tau has a standard
  # deviation of about 0.1 over random data sets.
  e2, sig, em = synthetic_e_data(a, 0.7, sigmaa=0.5, seed=6)
  args = (e2, sig, one, cen, em, flex.double(n, 0.5), a, dss, mod.f_total,
          flex.bool(n, True))
  tau = llgi_tncs.fit_tau(*args)
  lo, hi = tau.end_values()
  assert lo == hi and approx_equal(lo, 0.7, eps=0.25), (lo, hi)
  tau = llgi_tncs.fit_tau(*args, linear=True)
  assert approx_equal(sum(tau.end_values()) / 2, 0.7, eps=0.3), tau.end_values()
  # an error uncorrelated between the copies gives tau near 0
  e2, sig, em = synthetic_e_data(a, 0.0, sigmaa=0.5, seed=7)
  tau = llgi_tncs.fit_tau(e2, sig, one, cen, em, flex.double(n, 0.5), a, dss,
    mod.f_total, flex.bool(n, True))
  assert abs(tau.p0) < 0.25, tau.p0
  # coherence: with the consensus equal to the copies (A_cons = A), errors
  # with u = 1 + tau A correspond to Gamma = tau (1 - sigmaA^2) + sigmaA^2
  gamma_true = 0.9
  tau_eq = (gamma_true - 0.25) / 0.75
  e2, sig, em = synthetic_e_data(a, tau_eq, sigmaa=0.5, seed=8)
  coherence = llgi_tncs.fit_coherence(e2, sig, cen, em, flex.double(n, 0.5),
    a, a, dss, flex.bool(n, True))
  assert abs(coherence.gamma - gamma_true) < 0.15, coherence.gamma
  assert coherence.v < 0.5, coherence.v
  u, n_low = llgi_tncs.error_variance_modulation(flex.double(n, 0.5), a, a,
    coherence.at(dss))
  assert n_low == 0 and flex.min(u) > 0

def synthetic_fmodel(h, xrs, sigma_common=0.6, sigma_indep=0.15, d_min=3.5,
      seed=11):
  """ fmodel (target llgi) for the model xrs, with intensities from a
  "true" structure whose coordinate errors are shared by chains A and B
  (sigma_common; B's error is A's moved by the inverse of the symmetry
  operator that relates the copies, as in build_model) plus independent
  errors (sigma_indep); nacelle-like LLGI data with TEPS = 1 (RESN from
  binned mean intensities, Dobs/Feff by moment matching). """
  import mmtbx.f_model
  rnd = random.Random(seed)
  uc = CS.unit_cell()
  o = matrix.sqr(uc.orthogonalization_matrix())
  f = matrix.sqr(uc.fractionalization_matrix())
  r_inv = o * matrix.sqr(sgtbx.rt_mx("-x+1/2,-y,z+1/2").inverse().r().as_double()) * f
  chain_of = [a.chain().id for a in h.atoms()]
  res_atom = [(a.parent().parent().resseq_as_int(), a.name) for a in h.atoms()]
  common = {}
  shifts = flex.vec3_double()
  for cid, key in zip(chain_of, res_atom):
    g = lambda s: tuple(rnd.gauss(0, s / math.sqrt(3)) for k in range(3))
    if cid in "AB":
      c = common.setdefault(key, g(sigma_common))
      if cid == "B": c = tuple(r_inv * matrix.col(c))
    else:
      c = g(sigma_common)
    e = g(sigma_indep)
    shifts.append(tuple(c[k] + e[k] for k in range(3)))
  truth = xrs.deep_copy_scatterers()
  truth.set_sites_cart(xrs.sites_cart() + shifts)
  ms = miller.build_set(CS, anomalous_flag=False, d_min=d_min).resolution_filter(
    d_max=20).remove_systematic_absences()
  f_true = abs(ms.structure_factors_from_scatterers(truth,
    algorithm="direct").f_calc())
  fmodel = mmtbx.f_model.manager(xray_structure=xrs, f_obs=f_true,
    r_free_flags=f_true.generate_r_free_flags(fraction=0.1))
  fmodel.update_all_scales()
  f_obs = fmodel.f_obs()   # fmodel's order
  i_true = flex.pow2(f_obs.data())
  eps = f_obs.epsilons().data().as_double()
  f_obs.setup_binner(n_bins=8)
  sigma_n = flex.double(f_obs.size(), 0)
  for i_bin in f_obs.binner().range_used():
    sel = f_obs.binner().selection(i_bin)
    sigma_n.set_selected(sel, flex.mean((i_true / eps).select(sel)))
  sig = 0.05 * sigma_n * eps + 0.02 * i_true
  i_obs = i_true + sig * flex.double([rnd.gauss(0, 1) for i in range(f_obs.size())])
  resn = flex.sqrt(eps * sigma_n)
  centric = f_obs.centric_flags().data()
  n = f_obs.size()
  d = llgi_tncs.renormalise(i_obs, sig, resn, flex.double(n, 1.), centric)
  from libtbx import group_args
  params = llgi_hybrid.llgi_hybrid_params.extract()
  llgi_data = group_args(dobs=f_obs.array(data=d.dobs),
    feff=f_obs.array(data=d.feff), teps=f_obs.array(data=flex.double(n, 1.)),
    resn=f_obs.array(data=resn), e_obs_sq=f_obs.array(data=d.e_obs_sq),
    sig_e_obs_sq=f_obs.array(data=d.sig_e_obs_sq),
    force_exact=f_obs.array(data=d.force_exact), hybrid_params=params)
  fmodel.set_llgi_data(llgi_data)
  fmodel.update(f_obs=llgi_data.feff)
  fmodel.set_target_name("llgi")
  return fmodel

def exercise_fmodel_integration(h, xrs):
  from six.moves import cStringIO as StringIO
  import mmtbx.refinement.llgi_maps as llgi_maps
  xrs = xrs.deep_copy_scatterers()
  xrs.scatterers().flags_set_grads(state=False)
  xrs.scatterers().flags_set_grad_site(
    iselection=flex.size_t(range(xrs.scatterers().size())))
  fmodel = synthetic_fmodel(h, xrs)
  # no tNCS treatment when switched off or without tNCS groups
  p = llgi_tncs.llgi_tncs_params.extract()
  p.enabled = False
  assert llgi_tncs.setup(fmodel.llgi_data(), fmodel.f_obs(), h, xrs,
    params=p) is None
  no_b = h.atom_selection_cache().selection("not chain B")
  h_no_b = h.select(no_b)
  h_no_b.atoms().reset_i_seq()
  assert llgi_tncs.setup(fmodel.llgi_data(), fmodel.f_obs(), h_no_b,
    xrs.select(no_b)) is None
  log = StringIO()
  new = llgi_tncs.setup(fmodel.llgi_data(), fmodel.f_obs(), h, xrs, log=log)
  assert new is not None
  assert "tNCS group: chains A B" in log.getvalue()
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=new.feff)
  fmodel.update_all_scales()
  result = fmodel.update_llgi_sigmaa_scatfrac()
  ld = fmodel.llgi_data()
  # errors shared by the copies: the true copies are coherent
  assert result.coherence.gamma > 0.3, result.coherence.gamma
  assert result.tncs_summary.startswith("coherence")
  assert approx_equal(fmodel.f_obs().data(), ld.feff.data())
  assert flex.min(ld.tncs_t.data()) > 0 and flex.min(ld.tncs_u.data()) > 0
  u, t, sa = ld.tncs_u.data(), ld.tncs_t.data(), ld.sigmaa.data()
  # t exceeds (1 - sigmaA^2) u: the frame is valid everywhere
  assert flex.min(t - (1 - sa * sa) * u) > 0
  assert fmodel.llgi_sigmaa_is_current()
  # the coordinate target
  tf = fmodel.target_functor()(compute_gradients=True)
  assert math.isfinite(tf.target_work())
  g = tf.gradients_wrt_atomic_parameters(site=True).packed()
  assert g.size() == 3 * xrs.scatterers().size()
  # (d target/d Fcalc is checked by finite differences in
  # cctbx/xray/targets/tst_llgi_tncs.py)
  assert flex.max(flex.abs(g)) > 0 and all(math.isfinite(x) for x in g)
  # the Python frame (llgi_tncs.frame, as in llgi_hybrid_counts) and the
  # C++ functor make the same exact/Rice choice (the functor counts the
  # working set)
  exact = llgi_hybrid.get_hybrid(ld).exact_selection(
    llgi_tncs.frame(sa, u, t)[0])
  work = ~fmodel.r_free_flags().data()
  assert (exact & work).count(True) == tf.core_result.n_exact()
  assert fmodel.llgi_hybrid_counts().n_exact == exact.count(True)
  # statistics/map helper and LLGI maps
  mch = fmodel.map_calculation_helper_llgi()
  assert mch.n_exact > 0
  assert flex.min(mch.fom) >= 0 and flex.max(mch.fom) < 1
  maps = llgi_maps.compute(fmodel)
  for label, a in maps.arrays:
    assert flex.max(flex.abs(a.data())) < 1e10, label
  # sigmaA and ScatFrac refitted with atoms removed (as for the ADP weight,
  # mmtbx.refinement.weights.adp_gradients): the tNCS frame and the data
  # are kept; refitting the tNCS treatment itself needs the full model
  fd = fmodel.deep_copy()
  keep = h.atom_selection_cache().selection("not chain D")
  fd.update_xray_structure(xray_structure=fd.xray_structure.select(keep),
    update_f_calc=True)
  u0, f0 = fd.llgi_data().tncs_u.data().deep_copy(), fd.f_obs().data().deep_copy()
  r = fd.update_llgi_sigmaa_scatfrac(refit_tncs=False)
  assert r.sigmaa.size() == u0.size() and flex.max(r.sigmaa) > 0
  assert fd.llgi_data().tncs_u.data().all_eq(u0)
  assert fd.f_obs().data().all_eq(f0)
  try:
    fd.update_llgi_sigmaa_scatfrac()
  except AssertionError as e:
    assert "tNCS groups were found for a model with" in str(e)
  else:
    raise RuntimeError("expected the atom-count check")
  # a selection keeps the tNCS arrays aligned
  sel = fmodel.f_obs().d_spacings().data() > 4.0
  fs = fmodel.select(sel)
  assert fs.llgi_data().tncs_u.size() == sel.count(True)
  tf = fs.target_functor()(compute_gradients=False)
  assert math.isfinite(tf.target_work())
  # coherence = tau, u = 1 + tau A, fitted on all reflections (the free
  # set of this small data set is too small for tau)
  ld.tncs_state.params.coherence = "tau"
  ld.tncs_state.params.fit_reflections = "all"
  result = fmodel.update_llgi_sigmaa_scatfrac()
  assert result.tau.p0 > 0.3, result.tau.p0
  assert result.tncs_summary.startswith("tau")

def run():
  h, xrs = build_model()
  groups = exercise_groups(h)
  exercise_pair_term(groups, xrs)
  exercise_pair_term_centrosymmetric()
  exercise_renormalise()
  exercise_fits(groups, xrs)
  exercise_fmodel_integration(h, xrs)
  print("OK")

if __name__ == "__main__":
  run()
