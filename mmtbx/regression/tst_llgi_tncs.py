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
  """ (1/Nsym) sum_k sum_j w_j exp(2 pi i (T_k^T h).d_j) / sum_j w_j """
  d = sites_m - sites_n
  out = flex.double()
  rots = [op.r().as_double() for op in ms.space_group().smx()]
  for h, ss in zip(ms.indices(), ms.d_star_sq().data()):
    w = z_sq * flex.exp(-b_sum * (ss / 4.))
    tot = 0
    for r in rots:
      hr = matrix.row(h) * matrix.sqr(r)
      phase = d.dot(tuple(hr)) * (2 * math.pi)
      tot += flex.sum(w * flex.cos(phase))
    out.append(tot / (len(rots) * flex.sum(w)))
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
  assert approx_equal(mod.f_total, sum(f))
  assert approx_equal(mod.a, 2 * math.sqrt(f[0] * f[1]) * p_fft, eps=1e-12)

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

def run():
  h, xrs = build_model()
  groups = exercise_groups(h)
  exercise_pair_term(groups, xrs)
  exercise_renormalise()
  exercise_fits(groups, xrs)
  print("OK")

if __name__ == "__main__":
  run()
