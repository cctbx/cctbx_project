""" Anisotropic normalisation of the model's intensities for the LLGI target.

The LLGI likelihood compares the data and the model on the E scale, and
assumes both have unit mean square everywhere in reciprocal space. The data
are normalised by nacelle's RESN, which includes its anisotropy correction,
an ellipsoid that may change with resolution (phasertng.macaniso.
resolution_dependent). The model's normalisation (SigmaP for the sigmaA fit,
ScatFrac for the target) was isotropic, so wherever the model's own
intensities fall off differently with direction, for instance because TLS
or anisotropic ADPs carry part of the anisotropy, or because the anisotropic
scale in F_model and nacelle's correction differ, the model's E had the wrong
mean square by direction, which sigmaA cannot repair.

Here the model normalisation is

  Sigma(h) = Sigma_iso(s) A(h),   A(h) = exp(-sum_k b_k u_k(h)),

with u_k(h) = q^T T_k q / 4 (q the Cartesian reciprocal vector) for a basis
T_k of the symmetric, trace-free tensors invariant under the point group
(mode "ellipsoid"), and, for mode "resolution_dependent", also u_k(h) d*^2/
d*^2_max (the same family as nacelle's). The isotropic part stays in
Sigma_iso, which the caller computes as before from I/A (iso_fn). b is
fitted to the model's intensities by the Wilson likelihood
  sum_h w_h [ln Sigma(h) + I_h/Sigma(h)]   (w = 1 acentric, 1/2 centric),
by Newton steps alternating with Sigma_iso; at the solution the normalised
intensities I/Sigma have unit mean in each resolution shell and no residual
trend of the family in direction. Nothing finer than the family is fitted:
directional structure beyond it (the molecular transform, tNCS modulation)
is shared by the model and the true structure and stays in both.

Mode "auto" (resolve) takes the family of nacelle's correction:
resolution_dependent if its ANISOBETA changes with resolution, otherwise
ellipsoid. The model's family matching the data's works best in tests;
a model normalised with the resolution-dependent family against data
normalised with an ellipsoid gains freedom without correcting anything.
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from cctbx.array_family import flex
from libtbx import group_args

MODES = ("off", "ellipsoid", "resolution_dependent")

def tensor_basis(unit_cell, space_group):
  """ Orthonormal basis (Frobenius) of the symmetric, trace-free Cartesian
  tensors invariant under the point group of space_group: a list of 3x3
  numpy arrays (empty for cubic groups). """
  fm = np.array(unit_cell.fractionalization_matrix()).reshape(3, 3)
  finv = np.linalg.inv(fm.T)
  mats = [fm.T.dot(np.array(op.r().as_double()).reshape(3, 3).T).dot(finv)
          for op in space_group.all_ops()]
  basis = []
  for i in range(3):
    for j in range(i, 3):
      e = np.zeros((3, 3)); e[i, j] = e[j, i] = 1
      t = np.mean([m.T.dot(e).dot(m) for m in mats], axis=0)
      t = 0.5 * (t + t.T)
      t -= np.trace(t) / 3 * np.eye(3)
      for b in basis:
        t -= np.sum(t * b) * b
      norm = math.sqrt(np.sum(t * t))
      if(norm > 1.e-8): basis.append(t / norm)
  return basis

class anisotropic_factor(object):
  """ A(h) = exp(-sum_k b_k u_k(h)) for any reflections of the crystal. """

  def __init__(self, unit_cell, space_group, mode, d_star_sq_max, b=None):
    assert mode in MODES[1:]
    self.unit_cell = unit_cell
    self.basis = tensor_basis(unit_cell, space_group)
    self.mode = mode
    self.d_star_sq_max = d_star_sq_max
    self.b = np.zeros(self.n_parameters()) if b is None else np.asarray(b, float)

  def n_parameters(self):
    return len(self.basis) * (2 if self.mode == "resolution_dependent" else 1)

  def design(self, indices):
    """ u_k(h) for each reflection (rows) and parameter (columns) """
    h = np.array(indices, dtype=float).reshape(-1, 3)
    q = h.dot(np.array(self.unit_cell.fractionalization_matrix()).reshape(3, 3))
    u = np.column_stack([0.25 * np.einsum("ij,jk,ik->i", q, t, q)
      for t in self.basis]) if self.basis else np.zeros((h.shape[0], 0))
    if(self.mode == "resolution_dependent"):
      x = (q * q).sum(1) / self.d_star_sq_max
      u = np.column_stack([u, u * x[:, None]])
    return u

  def evaluate(self, indices):
    """ A(h) as flex.double """
    if(self.n_parameters() == 0): return flex.double(len(indices), 1)
    return flex.double(np.exp(-self.design(indices).dot(self.b)))

def fit(intensity, indices, centric_flags, unit_cell, space_group, iso_fn,
      mode, selection=None, n_outer=6, max_newton=30):
  """ Fit A(h) to the intensities (flex.double, epsilon already divided
  out) by the Wilson likelihood, alternating with the isotropic part
  iso_fn(I/A) -> Sigma_iso per reflection (flex.double, all reflections).
  selection: the reflections used in the fit (flex.bool, or None for all).

  Returns group_args(sigma (Sigma_iso*A, flex.double), sigma_iso, a (A),
  factor (anisotropic_factor, or None for mode "off")). """
  n = intensity.size()
  if(mode == "off"):
    sigma_iso = iso_fn(intensity)
    return group_args(sigma=sigma_iso, sigma_iso=sigma_iso,
      a=flex.double(n, 1), factor=None)
  sel = np.ones(n, bool) if selection is None else np.array(selection)
  i_np = np.array(intensity)
  d_star_sq = unit_cell.d_star_sq(indices).as_numpy_array()
  factor = anisotropic_factor(unit_cell, space_group, mode,
    d_star_sq_max=float(d_star_sq[sel].max()))
  if(factor.n_parameters() == 0):
    sigma_iso = iso_fn(intensity)
    return group_args(sigma=sigma_iso, sigma_iso=sigma_iso,
      a=flex.double(n, 1), factor=factor)
  u = factor.design(indices)
  w = np.where(np.array(centric_flags), 0.5, 1.0)
  us, ws = u[sel], w[sel]
  b = np.zeros(u.shape[1])
  a = np.ones(n)
  for outer in range(n_outer):
    sigma_iso = np.array(iso_fn(flex.double(i_np / a)))
    y = (i_np / sigma_iso)[sel]
    for it in range(max_newton):
      e = np.exp(us.dot(b))            # 1/A
      g = us.T.dot(ws * (y * e - 1))
      hess = (us * (ws * y * e)[:, None]).T.dot(us)
      step = np.linalg.solve(hess + 1.e-9 * np.eye(len(b)), g)
      b -= step
      if(np.max(np.abs(step)) < 1.e-7): break
    a = np.exp(-u.dot(b))
  factor.b = b
  sigma_iso = iso_fn(flex.double(i_np / a))
  a = flex.double(a)
  return group_args(sigma=sigma_iso * a, sigma_iso=sigma_iso, a=a,
    factor=factor)

# below this (in ln ANISOBETA, over the reflections), the d*^2 term of
# nacelle's anisotropy is rounding: about 1e-8 without
# resolution-dependent anisotropy, 0.4 to 8 with it in tests
RESOLUTION_DEPENDENT_THRESHOLD = 1.e-4

def data_anisotropy_is_resolution_dependent(anisobeta):
  """ True if nacelle's ANISOBETA (a miller array) changes with resolution
  (phasertng.macaniso.resolution_dependent kept its beta1): the d*^2 part
  of its exact form exp(h^T B0 h + d*^2 h^T B1 h) (llgi_maps.fit_anisobeta)
  is not negligible """
  import mmtbx.refinement.llgi_maps as llgi_maps
  fit = llgi_maps.fit_anisobeta(anisobeta.indices(), anisobeta.data())
  full = llgi_maps.fit_anisobeta(anisobeta.indices(), anisobeta.data(),
    unit_cell=anisobeta.unit_cell())
  b0 = np.log(np.array(fit.evaluate(anisobeta.indices())))
  b01 = np.log(np.array(full.evaluate(anisobeta.indices())))
  return bool(np.max(np.abs(b01 - b0)) > RESOLUTION_DEPENDENT_THRESHOLD)

def resolve(mode, llgi_data):
  """ mode, with "auto" replaced by the family of the data's anisotropy
  correction (llgi_data.anisobeta): resolution_dependent or ellipsoid
  (ellipsoid also without llgi_data or ANISOBETA) """
  if(mode != "auto"): return mode
  anisobeta = getattr(llgi_data, "anisobeta", None)
  if(anisobeta is not None and anisobeta.size() > 0 and
     data_anisotropy_is_resolution_dependent(anisobeta)):
    return "resolution_dependent"
  return "ellipsoid"

def mode_for(llgi_data):
  """ The model_anisotropy (resolved) of the sigmaA fit that made
  llgi_data's sigmaA (llgi_data.fit_params.e_scale, the default if None),
  "off" if there is none """
  fp = getattr(llgi_data, "fit_params", None)
  if(fp is None): return "off"
  e = getattr(fp, "e_scale", None)
  if(e is None):
    import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
    e = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
  return resolve(getattr(e, "model_anisotropy", "off"), llgi_data)
