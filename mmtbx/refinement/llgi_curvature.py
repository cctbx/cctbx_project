"""
Diagonal curvatures for preconditioning the LBFGS minimisation of
coordinates with the llgi target (mmtbx.refinement.minimization.lbfgs).

When the model is poor the data constrain it weakly (sigmaA is small), so
the data term's curvature is small next to that of the restraints, and
LBFGS, with only a scalar estimate of the inverse Hessian, makes little
progress along the directions that only the data determine (collective and
rigid-body shifts). A diagonal estimate of the Hessian, given to LBFGS as
its initial inverse Hessian, puts the shifts on the right scale for each
atom without changing the target.

Data term: the expected (Fisher) information of each reflection about the
model's structure factor, in the Gaussian approximation of the target:
with mean D Ecalc and variance V in the E scale (D = Dobs sigmaA,
V = TEPS - D^2; in the tNCS frame of mmtbx.refinement.llgi_tncs D = Dobs
a', V = 1 - D^2 and Ecalc scaled by s), I_E = D^2/V, halved for acentric
reflections (only the component along the model phase is informative),
and put on the scale of Fcalc by (k/(RESN sqrt(ScatFrac)))^2 times the
squared overall/anisotropic scale. For atom j and each Cartesian
component, averaging over directions and symmetry copies,

  H_j = (wx/N_work) (4 pi^2/3) n_sym occ_j^2
          sum_h I_h d*^2 f_j(d*)^2 exp(-B_j d*^2/2),

evaluated with sums over resolution bins.

Restraints: one curvature for all atoms, the mean of the diagonal of the
Hessian of the (normalised) restraint target, estimated with random +-1
probes and finite differences of the gradients (E[v (Hv)] = diag H).
"""
from __future__ import absolute_import, division, print_function
import math
from cctbx.array_family import flex

def fisher_information(fmodel):
  """ Per-reflection Fisher information of the llgi target about the
  complex Fcalc (on f_calc's scale), zero for test reflections and where the
  target has no contribution. """
  import mmtbx.refinement.llgi_tncs as llgi_tncs
  ld = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  sa = ld.sigmaa.data()
  sf = ld.scatfrac.data()
  dobs = ld.dobs.data()
  resn = ld.resn.data()
  teps = ld.teps.data()
  k = fmodel.scale_ml_wrapper()
  if(k <= 0): k = 1.
  u, t = llgi_tncs.arrays(ld)
  ok = ~fmodel.r_free_flags().data() & (sf > 0) & (resn > 0) & (sa > 0) & \
    (sa < 1)
  if(u is not None):
    # the t-normalised frame (llgi.h tncs_frame): a' and s
    a_p, s = llgi_tncs.frame(sa, u, t)
    ok &= a_p > 0
    d = dobs * a_p
    i_e = flex.pow2(d * s)
    v = 1. - d * d
  else:
    d = dobs * sa
    i_e = d * d
    v = teps - d * d
  v.set_selected(v < 1.e-3, 1.e-3)
  half = flex.double(n, 0.5)                 # acentric: half per component
  half.set_selected(f_obs.centric_flags().data(), 1.)
  i_e = i_e / v * half
  sf_safe = sf.deep_copy()
  sf_safe.set_selected(~ok, 1.)
  resn_safe = resn.deep_copy()
  resn_safe.set_selected(~ok, 1.)
  ec_per_fc = k / (resn_safe * flex.sqrt(sf_safe))
  k_total = fmodel.k_isotropic() * fmodel.k_anisotropic()
  result = i_e * flex.pow2(ec_per_fc * k_total)
  result.set_selected(~ok, 0.)
  return result

def data_curvatures(fmodel, xray_structure, n_bins=40):
  """ Per-scatterer data curvature (per Cartesian component) of the mean
  -log-likelihood over the working set, without the weight wx. """
  info = fisher_information(fmodel)
  dss = fmodel.f_obs().d_star_sq().data()
  n_work = max(1, (~fmodel.r_free_flags().data()).count(True))
  order = flex.sort_permutation(dss)
  bins = []
  n = dss.size()
  for b in range(n_bins):
    sel = order[b * n // n_bins:(b + 1) * n // n_bins]
    if(sel.size() == 0): continue
    s2 = dss.select(sel)
    bins.append((flex.mean(s2), flex.sum(info.select(sel) * s2)))
  reg = xray_structure.scattering_type_registry()
  f_sq = {}
  for t in reg.type_index_pairs_as_dict().keys():
    g = reg.gaussian_not_optional(scattering_type=t)
    f_sq[t] = [g.at_d_star_sq(s2) ** 2 for s2, _ in bins]
  n_sym = xray_structure.space_group().order_z()
  c = 4 * math.pi ** 2 / 3 * n_sym / n_work
  b_iso = xray_structure.extract_u_iso_or_u_equiv() * (8 * math.pi ** 2)
  result = flex.double()
  for i, sc in enumerate(xray_structure.scatterers()):
    fs = f_sq[sc.scattering_type]
    h = 0.
    for (s2, w), f2 in zip(bins, fs):
      h += w * f2 * math.exp(-b_iso[i] * s2 / 2)
    result.append(c * sc.occupancy ** 2 * h)
  return result

def restraint_curvature(model, sites_cart, selection, n_probes=8, eps=0.01,
      seed=0):
  """ Mean diagonal element of the Hessian of the restraint target (as the
  minimiser uses it, normalised) over the selected atoms' coordinates. """
  rm = model.restraints_manager
  hd = model.get_hd_selection()
  def gradients(sites):
    return rm.energies_sites(sites_cart=sites, compute_gradients=True,
      hd_selection=hd).gradients.as_double()
  mt = flex.mersenne_twister(seed=seed)
  n = sites_cart.size()
  mask = flex.double(n * 3, 0.)
  for i in selection.iselection():
    mask[3 * i] = mask[3 * i + 1] = mask[3 * i + 2] = 1.
  total, count = 0., 0
  for p in range(n_probes):
    v = (mt.random_double(size=n * 3) < 0.5).as_double() * 2 - 1
    v *= mask
    dv = flex.vec3_double(v * eps)
    hv = (gradients(sites_cart + dv) - gradients(sites_cart - dv)) / (2 * eps)
    total += flex.sum(v * hv)
    count += mask.count(1.)
  return max(total / max(count, 1), 0.)

def xyz_inverse_diagonal(fmodel, model, weights, xray_structure,
      restraints_weight):
  """ The inverse diagonal Hessian for LBFGS (flex.double, the minimiser's
  parameter order: three per scatterer refining its site) of
  wx*target + w*restraints. """
  sc = xray_structure.scatterers()
  refine = flex.bool([s.flags.grad_site() for s in sc])
  wx = weights.wx * weights.wx_scale
  h_data = data_curvatures(fmodel, xray_structure) * wx
  h_restr = 0.
  if(restraints_weight > 0 and model is not None and
     model.restraints_manager is not None):
    sites, sel = xray_structure.sites_cart(), refine
    ias = model.ias_manager.get_ias_selection() \
      if model.ias_manager is not None else None
    if(ias is not None and ias.count(True) > 0):
      # ias are not in the geometry restraints (as in the minimiser)
      sites, sel = sites.select(~ias), sel.select(~ias)
    h_restr = restraints_weight * restraint_curvature(model, sites, sel)
  result = flex.double()
  floor = max(1.e-6, 1.e-3 * flex.max(h_data + h_restr))
  for i in refine.iselection():
    h = max(h_data[i] + h_restr, floor)
    result.extend(flex.double(3, 1. / h))
  return result, h_data.select(refine), h_restr
