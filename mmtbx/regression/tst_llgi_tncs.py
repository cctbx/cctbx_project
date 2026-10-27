from __future__ import absolute_import, division, print_function
""" mmtbx.refinement.llgi_tncs: tNCS groups from a model, the pair
modulation A(h), renormalised data, and the rho/tau fits, on a synthetic
structure with a known tNCS relationship. """
import math, random
import iotbx.pdb
from cctbx import crystal, miller, sgtbx
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.math import superpose
import iotbx.phil
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
  # the upper bound: data whose optimum is above 1 (gamma_true = 1.15)
  tau_eq = (1.15 - 0.25) / 0.75
  e2, sig, em = synthetic_e_data(a, tau_eq, sigmaa=0.5, seed=8)
  args = (e2, sig, cen, em, flex.double(n, 0.5), a, a, dss,
    flex.bool(n, True))
  default = llgi_tncs.fit_coherence(*args)
  assert 0.95 < default.gamma <= 1., default.gamma
  wide = llgi_tncs.fit_coherence(*args, gamma_max=1.05)
  assert 1. < wide.gamma <= 1.05, wide.gamma
  # the start is below the bound
  start = llgi_tncs.coherence_curve(1.05, 0.2)
  wide = llgi_tncs.fit_coherence(*args, start=start, gamma_max=1.05)
  assert 1. < wide.gamma <= 1.05, wide.gamma

def synthetic_fmodel(h, xrs, sigma_common=0.6, sigma_indep=0.15, d_min=3.5,
      seed=11, source="nacelle"):
  """ fmodel (target llgi) for the model xrs, with intensities from a
  "true" structure whose coordinate errors are shared by chains A and B
  (sigma_common; B's error is A's moved by the inverse of the symmetry
  operator that relates the copies, as in build_model) plus independent
  errors (sigma_indep); nacelle-like LLGI data with TEPS = 1 (RESN from
  binned mean intensities, Dobs/Feff by moment matching). source="cctbx":
  as phenix.refinement.llgi_data.get_llgi_data_cctbx without tNCS groups
  (RESN from llgi_normalisation.fit, TEPS = 1), with .i_obs_base and
  .normalisation_source = "cctbx". """
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
  extra = {}
  if source == "cctbx":
    import mmtbx.refinement.llgi_normalisation as llgi_normalisation
    i_obs_base = f_obs.array(data=i_obs, sigmas=sig
      ).set_observation_type_xray_intensity()
    resn = llgi_normalisation.fit(i_obs_base).resn.data()
    extra = dict(i_obs_base=i_obs_base, normalisation_source="cctbx")
  else:
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
    force_exact=f_obs.array(data=d.force_exact), hybrid_params=params,
    **extra)
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
  # are kept; refitting the tNCS treatment itself needs the model with the
  # current atoms
  from libtbx.utils import Sorry
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
  except Sorry as e:
    assert "tNCS groups were found" in str(e)
  else:
    raise RuntimeError("expected the atom-count check")
  exercise_atoms_added(fmodel, h)
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

def exercise_atoms_added(fmodel, h):
  """ waters added after setup (as by ordered solvent): the tNCS groups are
  found again in the new model, and the treatment continues """
  from libtbx.utils import Sorry
  from six.moves import cStringIO as StringIO
  lines = [l for l in h.as_pdb_string().splitlines()
           if l.startswith("ATOM")]
  uc = CS.unit_cell()
  for k in range(20):
    xyz = uc.orthogonalize((0.05 * k, 0.9 - 0.04 * k, 0.5))
    lines.append("HETATM%5d  O   HOH W%4d    %8.3f%8.3f%8.3f  1.00 30.00"
      "           O" % (len(lines) + 1, k + 1, xyz[0], xyz[1], xyz[2]))
  h_w = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines("\n".join(lines))).construct_hierarchy()
  xrs_w = h_w.extract_xray_structure(crystal_symmetry=CS)
  fw = fmodel.deep_copy()
  fw.update_xray_structure(xray_structure=xrs_w, update_f_calc=True)
  # without the new model, a Sorry, not a failed assertion
  try:
    fw.update_llgi_sigmaa_scatfrac()
  except Sorry as e:
    assert "no model with the current atoms" in str(e), str(e)
  else:
    raise RuntimeError("expected a Sorry")
  # or with a model that does not match the structure
  try:
    fw.update_llgi_sigmaa_scatfrac(pdb_hierarchy=h)
  except Sorry as e:
    assert "no model with the current atoms" in str(e), str(e)
  else:
    raise RuntimeError("expected a Sorry")
  result = fw.update_llgi_sigmaa_scatfrac(pdb_hierarchy=h_w)
  groups = fw.llgi_data().tncs_state.groups
  assert groups.n_atoms == xrs_w.scatterers().size()
  assert [g.chain_ids for g in groups.groups] == [["A", "B"]]
  assert "groups found again" in result.tncs_summary
  assert result.coherence.gamma > 0.3, result.coherence.gamma
  # the next refit uses the new groups without finding them again
  result = fw.update_llgi_sigmaa_scatfrac()
  assert "groups found again" not in result.tncs_summary
  assert fw.llgi_data().tncs_state.groups is groups

# _nacelle_run with the code before the cctbx normalisation source was
# added (exercise_nacelle_source_unchanged; the new code reproduces them
# exactly on the machine that made them); the values from gamma on, which
# follow the sigmaA fit, since the model's E there is normalised by the
# moments ScatFrac (llgi_e_sigmaa.build_e_model)
NACELLE_REFERENCE = (
  ("rho", 0.787717865840551),
  ("t0_sum", 4893.009102968329),
  ("feff0_sum", 1320034.8770827255),
  ("gamma", 0.923527482773334),
  ("v", 0.15150228930869325),
  ("u_sum", 4891.683615049821),
  ("t_sum", 4891.996110755512),
  ("sigmaa_mean", 0.8972941943236018),
  ("feff_sum", 1320008.6653574002))

def _nacelle_run(h, xrs, frame):
  # the R-free flags use the global random state
  flex.set_random_seed(0)
  random.seed(0)
  fmodel = synthetic_fmodel(h, xrs)
  p = llgi_tncs.llgi_tncs_params.extract()
  p.frame = frame
  new = llgi_tncs.setup(fmodel.llgi_data(), fmodel.f_obs(), h, xrs, params=p)
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=new.feff)
  fmodel.update_all_scales()
  values = [("rho", new.tncs_state.tau.p0),
    ("t0_sum", flex.sum(new.tncs_t.data())),
    ("feff0_sum", flex.sum(new.feff.data()))]
  r = fmodel.update_llgi_sigmaa_scatfrac()
  ld = fmodel.llgi_data()
  values += [("gamma", r.coherence.gamma), ("v", r.coherence.v),
    ("u_sum", flex.sum(ld.tncs_u.data())),
    ("t_sum", flex.sum(ld.tncs_t.data())),
    ("sigmaa_mean", flex.mean(ld.sigmaa.data())),
    ("feff_sum", flex.sum(fmodel.f_obs().data()))]
  return fmodel, r, values

def exercise_nacelle_source_unchanged(h, xrs):
  """ the nacelle source (llgi_data without normalisation_source) gives
  the results of the code before the cctbx source was added, and the
  frame option does not apply to it """
  fm, rm, values = _nacelle_run(h, xrs, "model")
  for (name, v), (ref_name, ref) in zip(values, NACELLE_REFERENCE):
    assert name == ref_name
    # tolerance for other platforms and library versions (the fits stop
    # at xatol = 1e-4)
    assert abs(v - ref) <= 1e-5 * abs(ref), (name, v, ref)
  s = fm.llgi_data().tncs_state
  assert s.teps_data is None and s.normalisation is None
  assert s.t_model.all_eq(fm.llgi_data().tncs_t.data())
  assert not [w for w in rm.warnings if "u clamped" in w]
  assert "would change" not in rm.tncs_summary
  fd, rd, values_d = _nacelle_run(h, xrs, "data")
  assert values_d == values
  for name in ["tncs_u", "tncs_t", "sigmaa", "feff", "resn"]:
    assert getattr(fd.llgi_data(), name).data().all_eq(
      getattr(fm.llgi_data(), name).data()), name
  assert fd.f_obs().data().all_eq(fm.f_obs().data())

def exercise_gamma_max_passed(h, xrs):
  """ fit_and_renormalise bounds the coherence by gamma_max of the
  normalisation parameters with the cctbx source and by 1 with nacelle's """
  seen = []
  original = llgi_tncs.fit_coherence
  def recorder(*args, **kwds):
    seen.append(kwds.get("gamma_max"))
    return original(*args, **kwds)
  llgi_tncs.fit_coherence = recorder
  try:
    _nacelle_run(h, xrs, "model")
    assert seen and set(seen) == set([1.0]), seen
    del seen[:]
    fmodel = cctbx_fmodel(h, xrs)
    fmodel.update_llgi_sigmaa_scatfrac()
    assert seen and set(seen) == set([1.05]), seen
    state = fmodel.llgi_data().tncs_state
    state.normalisation_params.gamma_max = 1.2
    del seen[:]
    fmodel.update_llgi_sigmaa_scatfrac()
    assert seen and set(seen) == set([1.2]), seen
  finally:
    llgi_tncs.fit_coherence = original

def _max_rel(a, b):
  """ largest difference relative to the largest value """
  return flex.max(flex.abs(a - b)) / flex.max(flex.abs(b))

def _check_cctbx_base(llgi_data):
  """ RESN, E_obs^2 base and sigma base from the normalisation fit kept
  in tncs_state: E^2 = I/RESN^2 """
  s = llgi_data.tncs_state
  i = llgi_data.i_obs_base
  resn = s.normalisation.resn.data()
  assert _max_rel(llgi_data.e_obs_sq_base.data(),
    i.data() / flex.pow2(resn)) < 1e-12
  assert _max_rel(llgi_data.sig_e_obs_sq_base.data(),
    i.sigmas() / flex.pow2(resn)) < 1e-12
  assert s.teps_data.all_eq(s.normalisation.teps.data())

def exercise_cctbx_setup_frame(h, xrs):
  """ the cctbx source: setup refits the normalisation with the consensus
  pair terms of the tNCS groups and starts the frame at u = t = TEPS_data,
  with E_obs^2 base = I/RESN^2 from the new RESN """
  from six.moves import cStringIO as StringIO
  import mmtbx.refinement.llgi_normalisation as llgi_normalisation
  fmodel = synthetic_fmodel(h, xrs, source="cctbx")
  ld = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  # no tNCS treatment without tNCS groups, as for the nacelle source
  no_b = h.atom_selection_cache().selection("not chain B")
  h_no_b = h.select(no_b)
  h_no_b.atoms().reset_i_seq()
  assert llgi_tncs.setup(ld, f_obs, h_no_b, xrs.select(no_b)) is None
  log = StringIO()
  n_params = llgi_normalisation.normalisation_params.extract()
  new = llgi_tncs.setup(ld, f_obs, h, xrs, log=log,
    normalisation_params=n_params)
  assert "tNCS group: chains A B" in log.getvalue()
  assert "cctbx normalisation" in log.getvalue(), log.getvalue()
  s = new.tncs_state
  assert s.normalisation_params is n_params
  assert s.t_model is None
  # the per-group consensus pair terms (one group: their sum, exactly)
  pm = llgi_tncs.pair_modulations(f_obs, s.groups, xrs, per_group=True)
  a_groups = pm.a_groups["consensus"]
  assert len(a_groups) == 1 and len(pm.a_groups["model"]) == 1
  assert a_groups[0].all_eq(pm.a["consensus"])
  assert pm.a_groups["model"][0].all_eq(pm.a["model"])
  plain = llgi_tncs.pair_modulations(f_obs, s.groups, xrs)
  assert plain.a_groups is None and plain.a["model"].all_eq(pm.a["model"])
  fit = llgi_normalisation.fit(ld.i_obs_base, tncs_a_groups=a_groups,
    params=n_params)
  t = new.tncs_t.data()
  assert flex.max(flex.abs(t - fit.teps.data())) < 1e-12
  assert new.tncs_u.data().all_eq(t)
  assert flex.max(flex.abs(t - 1)) > 0.1   # the data are modulated
  assert fit.coherence.gamma > 0.3, fit.coherence.gamma
  resn = fit.resn.data()
  i = ld.i_obs_base
  assert _max_rel(new.e_obs_sq_base.data(), i.data() / flex.pow2(resn)) < 1e-12
  assert _max_rel(new.resn.data(), resn) < 1e-12
  assert _max_rel(new.e_obs_sq.data(), i.data() / (flex.pow2(resn) * t)) < 1e-12
  _check_cctbx_base(new)
  assert new.normalisation_selection.data().all_eq(fit.selection)
  return new

def cctbx_fmodel(h, xrs, **kwds):
  """ fmodel with the cctbx source after setup and scaling (kwds:
  synthetic_fmodel's) """
  flex.set_random_seed(0)
  random.seed(0)
  fmodel = synthetic_fmodel(h, xrs, source="cctbx", **kwds)
  new = llgi_tncs.setup(fmodel.llgi_data(), fmodel.f_obs(), h, xrs)
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=new.feff)
  fmodel.update_all_scales()
  return fmodel

def _refit(fmodel, frame):
  """ fit_and_renormalise on a copy of fmodel (the copies share the
  tNCS parameters) """
  f = fmodel.deep_copy()
  f.llgi_data().tncs_state.params.frame = frame
  return f, f.update_llgi_sigmaa_scatfrac()

def _frame_valid(llgi_data):
  sa = llgi_data.sigmaa.data()
  sel = sa > 0
  assert flex.max(sa) < 1
  a_p = llgi_tncs.frame(sa, llgi_data.tncs_u.data(),
    llgi_data.tncs_t.data())[0]
  return flex.min(a_p.select(sel)) > 0

def exercise_frame_data_vs_model(fmodel):
  """ after one fit_and_renormalise on the cctbx source, frame = data
  gives t = TEPS_data and frame = model t = t_model; the fits of the
  macrocycle do not depend on the frame """
  runs = dict((frame, _refit(fmodel, frame)) for frame in ("data", "model"))
  for frame, (f, r) in runs.items():
    ld = f.llgi_data()
    s = ld.tncs_state
    assert ("frame %s" % frame) in r.tncs_summary, r.tncs_summary
    _check_cctbx_base(ld)
    assert ld.resn.data().all_eq(s.normalisation.resn.data())
    assert approx_equal(f.f_obs().data(), ld.feff.data())
    assert "would change" in r.tncs_summary, r.tncs_summary
    assert _frame_valid(ld)
    assert math.isfinite(f.target_functor()(
      compute_gradients=False).target_work())
  ld_d, ld_m = runs["data"][0].llgi_data(), runs["model"][0].llgi_data()
  s_d, s_m = ld_d.tncs_state, ld_m.tncs_state
  assert ld_d.tncs_t.data().all_eq(s_d.teps_data)
  assert ld_m.tncs_t.data().all_eq(s_m.t_model)
  assert flex.max(flex.abs(s_d.t_model - s_m.t_model)) < 1e-10
  assert flex.max(flex.abs(ld_d.sigmaa.data() - ld_m.sigmaa.data())) < 1e-10
  assert flex.max(flex.abs(s_m.teps_data - s_m.t_model)) > 0.01
  # frame = model: t_model exceeds (1 - sigmaA^2) u, nothing clamped
  assert not [w for w in runs["model"][1].warnings if "u clamped" in w]
  assert runs["model"][1].n_u_clamped == 0
  # a selection keeps the normalisation's intensities and fit selection
  # aligned
  f = runs["data"][0]
  f.llgi_data().tncs_state.params.frame = "data"
  sel = f.f_obs().d_spacings().data() > 4.0
  fs = f.select(sel)
  r = fs.update_llgi_sigmaa_scatfrac()
  ls = fs.llgi_data()
  assert ls.i_obs_base.size() == sel.count(True)
  assert ls.normalisation_selection.size() == sel.count(True)
  assert ls.tncs_t.data().all_eq(ls.tncs_state.teps_data)
  _check_cctbx_base(ls)
  return runs

def exercise_u_clamp_reported(fmodel):
  """ the frame-validity clamp u = (1 - U_MIN) t/(1 - sigmaA^2) where
  t <= (1 - sigmaA^2) u: on arrays, and in fit_and_renormalise with
  TEPS_data forced to 0.5 t_model (fmodel: with sigmaA low enough, below
  about 0.7, for that to need the clamp) """
  import mmtbx.refinement.llgi_normalisation as llgi_normalisation
  k = 1 - llgi_tncs.U_MIN
  sa = flex.double([0., 0.5, 0.5, 0.5, 0.9, 1.0])
  u = flex.double([1., 1., 1., 1., 2., 1.])
  t = flex.double([1., 1., 0.75, 0.7, 0.1, 0.5])
  u2, n = llgi_tncs.clamp_u_to_frame(sa, u, t)
  assert n == 3
  assert approx_equal(u2, [1., 1., k * 0.75 / 0.75, k * 0.7 / 0.75,
    k * 0.1 / 0.19, 1.])
  assert u.all_eq(flex.double([1., 1., 1., 1., 2., 1.]))   # not in place
  a_p = llgi_tncs.frame(sa, u2, t)[0]
  assert flex.min(a_p[1:5]) > 0
  # fit_and_renormalise
  fm, rm = _refit(fmodel, "model")
  t_model = fm.llgi_data().tncs_state.t_model
  u_model = fm.llgi_data().tncs_u.data()
  real_fit = llgi_normalisation.fit
  def forced_fit(*args, **kwds):
    r = real_fit(*args, **kwds)
    r.teps = r.teps.customized_copy(data=0.5 * t_model)
    return r
  llgi_normalisation.fit = forced_fit
  try:
    fc, rc = _refit(fmodel, "data")
  finally:
    llgi_normalisation.fit = real_fit
  ld = fc.llgi_data()
  s = ld.tncs_state
  assert flex.max(flex.abs(s.t_model - t_model)) < 1e-10
  sa, u, t = ld.sigmaa.data(), ld.tncs_u.data(), ld.tncs_t.data()
  assert flex.max(flex.abs(t - 0.5 * t_model)) < 1e-12
  expected = (sa > 0) & (sa < 1) & (t <= (1 - sa * sa) * u_model)
  n = expected.count(True)
  assert 0 < n < sa.size(), n
  assert (llgi_tncs.U_CLAMP_WARNING % n) in rc.warnings, rc.warnings
  assert "tNCS: u clamped for %d reflections where t <= " \
    "(1 - sigmaA^2) u" % n in rc.warnings
  assert _frame_valid(ld)
  sel = expected.iselection()
  assert approx_equal(u.select(sel),
    k * t.select(sel) / (1 - flex.pow2(sa.select(sel))))
  assert flex.max(flex.abs(u.select(~expected) - u_model.select(~expected))) \
    < 1e-10
  # the scale of RESN (put on the model's scale by mmtbx.f_model.manager.
  # apply_scale_k1_to_f_obs) is kept when RESN is refitted
  fs = fmodel.deep_copy()
  ld0 = fs.llgi_data()
  ld0.tncs_state.params.frame = "data"
  fs.update(f_obs=fs.f_obs().customized_copy(data=fs.f_obs().data() / 2))
  fs.set_llgi_data(llgi_hybrid.replace_llgi_data(ld0,
    feff=ld0.feff.customized_copy(data=ld0.feff.data() / 2),
    resn=ld0.resn.customized_copy(data=ld0.resn.data() / 2)))
  fs.update_llgi_sigmaa_scatfrac()
  ls = fs.llgi_data()
  ratio = ls.resn.data() / ls.tncs_state.normalisation.resn.data()
  assert flex.max(flex.abs(ratio - 0.5)) < 1e-9
  _check_cctbx_base(ls)

def exercise_precomputed_pair_terms(h, xrs):
  """ setup with the consensus pair terms given (pair_terms=) gives the
  same llgi_data as without, and evaluates none; select_pair_terms maps
  pair terms onto a reordered subset of indices, None if one is missing """
  fmodel = synthetic_fmodel(h, xrs, source="cctbx")
  ld = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  base = llgi_tncs.setup(ld, f_obs, h, xrs)
  groups = base.tncs_state.groups
  # only the consensus basis is computed without pair_terms
  real = llgi_tncs.pair_modulations
  bases = []
  def spy(*args, **kwds):
    bases.append(args[3] if len(args) > 3 else kwds.get("bases"))
    return real(*args, **kwds)
  llgi_tncs.pair_modulations = spy
  try:
    llgi_tncs.setup(ld, f_obs, h, xrs)
  finally:
    llgi_tncs.pair_modulations = real
  assert bases == [("consensus",)], bases
  # given: nothing is evaluated, and the numbers are the same
  pm = real(f_obs, groups, xrs, ("consensus",), per_group=True)
  calls = []
  def fail(*args, **kwds):
    calls.append(1)
    return real(*args, **kwds)
  llgi_tncs.pair_modulations = fail
  real_pair_term = llgi_tncs.pair_term
  llgi_tncs.pair_term = fail
  try:
    given = llgi_tncs.setup(ld, f_obs, h, xrs, pair_terms=pm)
  finally:
    llgi_tncs.pair_modulations = real
    llgi_tncs.pair_term = real_pair_term
  assert not calls
  for name in ("tncs_u", "tncs_t", "resn", "feff", "e_obs_sq_base",
      "sig_e_obs_sq_base", "e_obs_sq", "dobs"):
    assert flex.max(flex.abs(getattr(given, name).data() -
      getattr(base, name).data())) < 1e-12, name
  assert abs(given.tncs_state.coherence.gamma -
    base.tncs_state.coherence.gamma) < 1e-12
  assert abs(given.tncs_state.f_total - base.tncs_state.f_total) < 1e-12
  # wrong size
  try:
    llgi_tncs.setup(ld, f_obs, h, xrs, pair_terms=llgi_tncs.select_pair_terms(
      pm, f_obs, f_obs.select(f_obs.d_spacings().data() > 4.0)))
  except AssertionError:
    pass
  else:
    raise RuntimeError("expected the size check")
  # the mapping
  sub = f_obs.select(f_obs.d_spacings().data() > 4.0)
  perm = flex.size_t(list(reversed(range(sub.size()))))
  sub = sub.select(perm)
  mapped = llgi_tncs.select_pair_terms(pm, f_obs, sub)
  assert mapped.a["consensus"].size() == sub.size()
  assert mapped.f_total == pm.f_total and mapped.fractions is pm.fractions
  # the values are those of pm at the reordered indices ...
  lookup = dict((h_, i) for i, h_ in enumerate(f_obs.indices()))
  isel = flex.size_t([lookup[h_] for h_ in sub.indices()])
  assert mapped.a["consensus"].all_eq(pm.a["consensus"].select(isel))
  assert len(mapped.a_groups["consensus"]) == len(pm.a_groups["consensus"])
  for m, g in zip(mapped.a_groups["consensus"], pm.a_groups["consensus"]):
    assert m.all_eq(g.select(isel))
  # ... which agree with the terms computed on sub itself to the accuracy of
  # the FFT (quality_factor 1e5: ~0.001; the grid depends on the set)
  direct = real(sub, groups, xrs, ("consensus",), per_group=True)
  assert flex.max(flex.abs(mapped.a["consensus"] - direct.a["consensus"])) \
    < 0.01
  # without a_groups
  plain = real(f_obs, groups, xrs, ("consensus",))
  assert llgi_tncs.select_pair_terms(plain, f_obs, sub).a_groups is None
  # an index of the target missing from the source
  src = f_obs.select(f_obs.d_spacings().data() > 4.0)
  src_pm = llgi_tncs.select_pair_terms(pm, f_obs, src)
  assert llgi_tncs.select_pair_terms(src_pm, src, f_obs) is None

def exercise_refit_clamps_u(fmodel):
  """ refit_sigmaa_scatfrac (refit_tncs=False, as adp_gradients) with
  sigmaA lowered so that the stored (u, t) would be frame-invalid: u is
  clamped, stored, counted and reported; t and the data are kept """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  f, _ = _refit(fmodel, "data")
  ld0 = f.llgi_data()
  t0, u0 = ld0.tncs_t.data().deep_copy(), ld0.tncs_u.data().deep_copy()
  e0 = ld0.e_obs_sq.data().deep_copy()
  real = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel
  def lowered(*args, **kwds):
    r = real(*args, **kwds)
    r.sigmaa = 0.6 * r.sigmaa
    return r
  llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel = lowered
  try:
    r = llgi_tncs.refit_sigmaa_scatfrac(f)
  finally:
    llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel = real
  ld = f.llgi_data()
  sa, u, t = ld.sigmaa.data(), ld.tncs_u.data(), ld.tncs_t.data()
  expected = (sa > 0) & (sa < 1) & (t <= (1 - sa * sa) * u0)
  n = expected.count(True)
  assert n > 0
  assert r.n_u_clamped == n
  assert (llgi_tncs.U_CLAMP_WARNING % n) in r.warnings, r.warnings
  assert t.all_eq(t0) and ld.e_obs_sq.data().all_eq(e0)
  assert _frame_valid(ld)
  sel = (sa > 0) & (sa < 1)
  a_p = llgi_tncs.frame(sa, u, t)[0]
  assert flex.min(a_p.select(sel)) > 0
  assert u.select(~expected).all_eq(u0.select(~expected))
  # no clamp needed: nothing changes
  f2, _ = _refit(fmodel, "model")
  u2 = f2.llgi_data().tncs_u.data().deep_copy()
  r2 = llgi_tncs.refit_sigmaa_scatfrac(f2)
  assert r2.n_u_clamped == 0
  assert not [w for w in r2.warnings if "u clamped" in w]
  assert f2.llgi_data().tncs_u.data().all_eq(u2)

def monitor_inputs(seed=11, d_min=3.0, sig_frac=0.1):
  """ (f_obs, free flags, i_obs, resn, A_1, A_2) on random P 21 21 21 data """
  ms = miller.build_set(CS, anomalous_flag=False, d_min=d_min)
  n = ms.size()
  flex.set_random_seed(seed)
  a1 = flex.random_double(n) * 0.8
  a2 = flex.random_double(n) * 0.8
  resn = flex.random_double(n) + 0.5
  i_data = -flex.log(flex.random_double(n) + 1e-3) * resn * resn
  sigmas = i_data * sig_frac + 0.05
  i_obs = ms.array(data=i_data, sigmas=sigmas).set_observation_type_xray_intensity()
  free = flex.random_double(n) < 0.1
  return ms.array(data=flex.double(n, 1.)), free, i_obs, resn, a1, a2

def exercise_monitor_identical():
  f_obs, free, i_obs, resn, a1, a2 = monitor_inputs()
  t = 1 + a1 + a2
  c = llgi_tncs.coherence_curve(0.9, 0.3)
  r = llgi_tncs.monitor(t, t.deep_copy(), f_obs, free, [a1, a2], i_obs, resn,
    c, c)
  assert len(r.table) > 0
  assert sum(row[3] for row in r.table) == f_obs.size()
  n_big = 0
  for d_max, d_min, cls, n, corr, mean_ratio, rms in r.table:
    assert d_max >= d_min
    assert approx_equal(mean_ratio, 1.) and approx_equal(rms, 0.)
    if(n >= 10):
      n_big += 1
      assert approx_equal(corr, 1.), (n, corr)
  assert n_big > 0
  assert r.ll_free_data == r.ll_free_model
  assert approx_equal(r.ll_free_data - r.ll_free_model, 0.)
  assert r.warnings == [], r.warnings
  assert max(abs(x) for x in r.group_coefficients) < 1e-12
  assert set(row[2] for row in r.table) <= set(["weak", "neutral", "strong"])
  from six.moves import cStringIO as StringIO
  out = StringIO()
  llgi_tncs.show_monitor(r, out)
  assert "gamma = 0.900" in out.getvalue() and "corr_ln" in out.getvalue()

def exercise_monitor_classes():
  f_obs, free, i_obs, resn, a1, a2 = monitor_inputs()
  t_model = flex.double(f_obs.size(), 1.)
  teps = flex.double(f_obs.size(), 1.)
  teps.set_selected(flex.size_t_range(0, 100), 0.4)
  teps.set_selected(flex.size_t_range(100, 200), 2.0)
  c = llgi_tncs.coherence_curve(0.9, 0.3)
  r = llgi_tncs.monitor(teps, t_model, f_obs, free, [a1], i_obs, resn, c, c)
  count = dict((k, sum(row[3] for row in r.table if row[2] == k))
    for k in ("weak", "neutral", "strong"))
  assert count == dict(weak=100, strong=100,
    neutral=f_obs.size() - 200), count
  for row in r.table:
    if(row[2] == "weak"): assert approx_equal(row[5], 0.4)
    if(row[2] == "strong"): assert approx_equal(row[5], 2.0)

def exercise_monitor_group_stands_out():
  f_obs, free, i_obs, resn, a1, a2 = monitor_inputs()
  t_model = 1 + a1 + a2
  teps = 1 + a1 + 1.3 * a2
  c = llgi_tncs.coherence_curve(0.9, 0.3)
  r = llgi_tncs.monitor(teps, t_model, f_obs, free, [a1, a2], i_obs, resn,
    c, c)
  assert abs(r.group_coefficients[1]) > 0.05, r.group_coefficients
  assert [w for w in r.warnings if w.startswith("tNCS monitor: group 2 "
    "stands out (coefficient ")], (r.warnings, r.group_coefficients)
  assert not [w for w in r.warnings if "group 1 " in w], r.warnings
  assert ("tNCS monitor: group %d stands out (coefficient %.3f)" % (
    2, r.group_coefficients[1])) in r.warnings

def exercise_monitor_reports_once():
  f_obs, free, i_obs, resn, a1, a2 = monitor_inputs()
  t_model = 1 + a1 + a2
  teps = 1 + a1 + 1.3 * a2
  c = llgi_tncs.coherence_curve(0.9, 0.3)
  def run(teps, previous=None):
    return llgi_tncs.monitor(teps, t_model, f_obs, free, [a1, a2], i_obs,
      resn, c, c, previous=previous)
  def stands_out(r):
    return [w for w in r.warnings if "stands out" in w]
  first = run(teps)
  assert len(stands_out(first)) == 1, first.warnings
  assert list(first.reported.keys()) == [1], first.reported
  assert first.reported[1] == first.group_coefficients[1]
  # the same inputs again: not reported again
  second = run(teps, first)
  assert not stands_out(second), second.warnings
  assert second.reported == first.reported
  assert len(second.group_coefficients) == 2
  # a coefficient that moved by more than MONITOR_REPORT_CHANGE is reported
  # again; the reported value is updated
  def with_reported(reported):
    return group_args(ll_free_data=first.ll_free_data,
      ll_free_model=first.ll_free_model, reported=reported)
  moved = with_reported({
    1: first.group_coefficients[1] - 2 * llgi_tncs.MONITOR_REPORT_CHANGE})
  third = run(teps, moved)
  assert len(stands_out(third)) == 1, third.warnings
  assert third.reported == first.reported, third.reported
  # a change below the threshold is not reported and the old value is kept
  near = with_reported({
    1: first.group_coefficients[1] - 0.5 * llgi_tncs.MONITOR_REPORT_CHANGE})
  fourth = run(teps, near)
  assert not stands_out(fourth), fourth.warnings
  assert fourth.reported == near.reported, fourth.reported
  # a group that no longer stands out is dropped, and reported when it
  # stands out again
  quiet = run(1 + a1 + a2, second)
  assert not stands_out(quiet) and quiet.reported == {}, quiet.reported
  again = run(teps, quiet)
  assert len(stands_out(again)) == 1, again.warnings
  assert again.reported == first.reported

def exercise_monitor_gap_growth():
  f_obs, free, i_obs, resn, a1, a2 = monitor_inputs()
  t_model = 1 + a1 + a2
  teps = 1 + a1 + 1.3 * a2
  c = llgi_tncs.coherence_curve(0.9, 0.3)
  first = llgi_tncs.monitor(teps, t_model, f_obs, free, [a1, a2], i_obs,
    resn, c, c)
  gap = first.ll_free_data - first.ll_free_model
  assert abs(gap) > 1e-6
  prev = group_args(ll_free_data=gap - 0.01, ll_free_model=0.)
  r = llgi_tncs.monitor(teps, t_model, f_obs, free, [a1, a2], i_obs, resn,
    c, c, previous=prev)
  assert ("tNCS monitor: free-set gap grew from %.4f to %.4f nats/reflection"
    % (gap - 0.01, gap)) in r.warnings, r.warnings
  prev = group_args(ll_free_data=gap - 0.001, ll_free_model=0.)
  r = llgi_tncs.monitor(teps, t_model, f_obs, free, [a1, a2], i_obs, resn,
    c, c, previous=prev)
  assert not [w for w in r.warnings if "free-set gap" in w], r.warnings

def exercise_monitor_in_fit_and_renormalise(fmodel):
  """ the cctbx source stores the monitor, passes it on as previous, and
  reports its warnings """
  f = fmodel.deep_copy()
  r1 = f.update_llgi_sigmaa_scatfrac()
  m1 = f.llgi_data().tncs_state.monitor
  assert m1 is not None and len(m1.table) > 0
  assert len(m1.group_coefficients) == len(f.llgi_data().tncs_state.groups)
  assert [w for w in r1.warnings if w.startswith("tNCS monitor")] == \
    m1.warnings
  r2 = f.update_llgi_sigmaa_scatfrac()
  m2 = f.llgi_data().tncs_state.monitor
  assert m2 is not m1
  gap1 = m1.ll_free_data - m1.ll_free_model
  gap2 = m2.ll_free_data - m2.ll_free_model
  assert (gap2 - gap1 > 0.002) == bool(
    [w for w in r2.warnings if "free-set gap" in w]), (gap1, gap2)

def build_two_domain(rot2_deg, seed=3):
  """ chain A (domain 1: residues 1-40, domain 2: residues 41-80, poly-ALA)
  and chain B, a copy translated, with domain 1 rotated by 2 degrees and
  domain 2 by rot2_deg about its own centre. """
  rnd = random.Random(seed)
  def random_sites(n):
    out = []
    while len(out) < n:
      v = [rnd.uniform(-10, 10) for k in range(3)]
      if sum(x * x for x in v) < 100: out.append(matrix.col(v))
    return out
  names = ["N", "CA", "C", "O", "CB"]
  uc = CS.unit_cell()
  centres = [matrix.col(uc.orthogonalize(f))
             for f in ((0.2, 0.2, 0.25), (0.2, 0.2, 0.6))]
  t_b = matrix.col(uc.orthogonalize((0.5, 0.45, 0.02)))
  axis = matrix.col((1, 2, 0.5)).normalize()
  lines = []
  def add(chain, i_res, name, xyz):
    lines.append("ATOM  %5d  %-3s ALA %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f"
      "           %s" % (len(lines) + 1, name, chain, i_res, xyz[0], xyz[1],
      xyz[2], 20., name[0]))
  domains = [random_sites(200), random_sites(200)]
  for chain in "AB":
    for i_dom, (sites, c) in enumerate(zip(domains, centres)):
      if chain == "B":
        r = axis.axis_and_angle_as_r3_rotation_matrix(
          (2., rot2_deg)[i_dom], deg=True)
        sites = [r * s for s in sites]
      for k, s in enumerate(sites):
        add(chain, 40 * i_dom + k // 5 + 1, names[k % 5],
            s + c + (t_b if chain == "B" else matrix.col((0, 0, 0))))
  pdb_inp = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines("\n".join(lines)))
  return pdb_inp.construct_hierarchy()

DOMAIN_SELECTIONS = [
  ["chain A and resseq 1:40", "chain B and resseq 1:40"],
  ["chain A and resseq 41:80", "chain B and resseq 41:80"]]

def exercise_domain_groups():
  h = build_two_domain(20.)
  groups = llgi_tncs.tncs_groups.from_selections(h, CS, DOMAIN_SELECTIONS,
    angle_max=30.)
  assert len(groups) == 2 and groups.n_atoms == h.atoms().size()
  uc = CS.unit_cell()
  sites_frac = uc.fractionalize(h.atoms().extract_xyz())
  for g, applied, sels in zip(groups.groups, (2., 20.), DOMAIN_SELECTIONS):
    assert g.chain_ids == sels
    assert len(g.i_seqs) == 2 and len(g.chain_i_seqs) == 2
    assert g.i_seqs[0].size() == g.i_seqs[1].size() == 200
    assert g.chain_i_seqs[0].size() == 200
    assert g.ops[0].is_unit_mx() and g.ops[1].is_unit_mx()
    assert g.angles[0] == 0.
    assert approx_equal(g.angles[1], applied, eps=0.01)
    # the superposition of copy 2 on copy 1 leaves less than 0.01 A
    lsq = superpose.least_squares_fit(
      reference_sites=uc.orthogonalize(
        llgi_tncs.apply_op(g.ops[1], sites_frac.select(g.i_seqs[1]))),
      other_sites=uc.orthogonalize(sites_frac.select(g.i_seqs[0])))
    assert approx_equal(lsq.r.rotation_angle(), applied, eps=0.01)
    moved = lsq.other_sites_best_fit()
    ref = uc.orthogonalize(
      llgi_tncs.apply_op(g.ops[1], sites_frac.select(g.i_seqs[1])))
    assert approx_equal(flex.mean((moved - ref).dot()) ** 0.5, 0., eps=0.01)
  # the same through the phil parameters
  params = llgi_tncs.llgi_tncs_params.fetch(iotbx.phil.parse(
    "group { copy = 'chain A and resseq 1:40'\n copy = 'chain B and resseq 1:40' }\n"
    "group { copy = 'chain A and resseq 41:80'\n copy = 'chain B and resseq 41:80' }\n"
    "angle_max = 30\n")).extract()
  g2 = llgi_tncs.groups_from_params(h, CS, params)
  assert len(g2) == 2
  assert [list(g.i_seqs[1]) for g in g2.groups] == \
         [list(g.i_seqs[1]) for g in groups.groups]
  from six.moves import cStringIO as StringIO
  log = StringIO()
  g2.show(log)
  assert "'chain A and resseq 1:40'" in log.getvalue()
  # no group given: automatic, from whole chains
  auto = llgi_tncs.groups_from_params(h, CS,
    llgi_tncs.llgi_tncs_params.extract())
  assert len(auto) == 1 and auto.groups[0].chain_ids == ["A", "B"]

def exercise_domain_group_angle_error():
  from libtbx.utils import Sorry
  h = build_two_domain(20.)
  try:
    llgi_tncs.tncs_groups.from_selections(h, CS, DOMAIN_SELECTIONS,
      angle_max=10.)
  except Sorry as e:
    assert str(e) == ("tNCS group 2: copy 'chain B and resseq 41:80' is "
      "related to 'chain A and resseq 41:80' by 20.0 degrees, more than "
      "angle_max"), str(e)
  else:
    raise AssertionError("Sorry expected")

def exercise_selected_hierarchy():
  h = build_two_domain(20.)
  full = llgi_tncs.tncs_groups.from_selections(h, CS, DOMAIN_SELECTIONS,
    angle_max=30.)
  # drop the first 10 atoms of chain A: i_seq values of the rest are stale
  sel = flex.bool(h.atoms().size(), True)
  for i in range(10): sel[i] = False
  sub = h.select(sel)
  stale = list(sub.atoms().extract_i_seq())
  assert stale[0] == 10
  groups = llgi_tncs.tncs_groups.from_selections(sub, CS, [
    ["chain A and resseq 3:40", "chain B and resseq 3:40"],
    DOMAIN_SELECTIONS[1]], angle_max=30.)
  assert list(sub.atoms().extract_i_seq()) == stale   # input unchanged
  assert groups.n_atoms == sub.atoms().size() == h.atoms().size() - 10
  xyz_full, xyz_sub = h.atoms().extract_xyz(), sub.atoms().extract_xyz()
  g_full, g = full.groups[1], groups.groups[1]
  for a, b in zip(g_full.i_seqs, g.i_seqs):
    assert approx_equal(list(xyz_full.select(a).as_double()),
                        list(xyz_sub.select(b).as_double()))
  assert approx_equal(g.angles[1], g_full.angles[1], eps=1e-6)
  assert [list(v) for v in g.chain_i_seqs] == \
         [[i - 10 for i in v] for v in g_full.chain_i_seqs]
  # group 1: residues 3:40 keep their atoms (residues 1-2 hold the 10 dropped
  # atoms of chain A and are not in the selections)
  assert groups.groups[0].i_seqs[0].size() == 190
  assert approx_equal(groups.groups[0].angles[1], 2., eps=0.01)

def run():
  h, xrs = build_model()
  groups = exercise_groups(h)
  exercise_pair_term(groups, xrs)
  exercise_pair_term_centrosymmetric()
  exercise_renormalise()
  exercise_fits(groups, xrs)
  exercise_fmodel_integration(h, xrs)
  exercise_nacelle_source_unchanged(h, xrs)
  exercise_gamma_max_passed(h, xrs)
  exercise_cctbx_setup_frame(h, xrs)
  exercise_frame_data_vs_model(cctbx_fmodel(h, xrs))
  exercise_u_clamp_reported(cctbx_fmodel(h, xrs, sigma_common=1.0,
    sigma_indep=0.5))
  exercise_precomputed_pair_terms(h, xrs)
  exercise_refit_clamps_u(cctbx_fmodel(h, xrs, sigma_common=1.0,
    sigma_indep=0.5))
  exercise_monitor_identical()
  exercise_monitor_classes()
  exercise_monitor_group_stands_out()
  exercise_monitor_reports_once()
  exercise_monitor_gap_growth()
  exercise_monitor_in_fit_and_renormalise(cctbx_fmodel(h, xrs))
  exercise_domain_groups()
  exercise_domain_group_angle_error()
  exercise_selected_hierarchy()
  print("OK")

if __name__ == "__main__":
  run()
