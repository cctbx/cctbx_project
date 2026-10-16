from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx.array_family import flex
from libtbx import group_args
import iotbx.phil

""" LLGI map coefficients from the exact posterior of E. On the E scale,
with Ec = |Emodel| along the model phase, sigmaA (clamped to SIGMAA_MAX)
and the information fraction phi(sigmaA, sigma(E_obs^2))
(mmtbx.refinement.llgi_phi_table), three coefficients along the model phase:

  bias-reduced   M = (<E> - c*Ec)/(1 - c*sigmaA),  c = (1 - phi)*sigmaA
                 (tends to 0 as the information in a reflection goes to 0)
  fill           2<E> - sigmaA*Ec (acentric), <E> (centric)
                 (tends to sigmaA*Ec, the LLGI D*Fc, as the information goes
                 to 0: a generalisation of the filled 2mFo-DFc map)
  difference     <E> - sigmaA*Ec

<E> is the exact posterior mean along the model phase (cctbx.xray.
llgi_exact_evaluate) for every reflection, whatever form the target uses.
Its model bias is sigmaA*(1 - phi); M0 = <E> - c*Ec removes exactly that
much of Ec. Then E[M0*conj(E)] = phi and E[|M0|^2] = phi*(1 - c*sigmaA)
(from the definition of phi, with E[<E>*conj(Ec)] = sigmaA), so the scale
1/(1 - c*sigmaA) minimises the mean-square error E[|M - E|^2] among
multiples of M0 (checked by Monte Carlo).

F scale: sqrt(TEPS)*RESN/sqrt(ANISOBETA), the data's own normalisation with
nacelle's anisotropic factor divided out (the overall isotropic falloff of
the data kept), or, for the _aniso maps (only with anisotropic_maps=True),
sqrt(TEPS)*RESN with the anisotropy left in. Maps (labels in LABELS):

  IWT            bias-reduced            IWT_aniso       (anisotropy left in)
  IWT_fill       fill                    IWT_fill_aniso
  dIWT           difference              dIWT_aniso

The anisotropy-removed maps correlate better with the true map and the true
difference map, most clearly for anisotropic data. Without an ANISOBETA
column the anisotropy cannot be removed, and the maps without suffix keep it.

IWT and dIWT have the reflections used by the LLGI target. IWT_fill adds
every other reflection to the resolution limit (unmeasured, or measured but
left out of the target as uninformative) with sigmaA*Ec, i.e. the value the
fill coefficient tends to; there the F scale is RESN/sqrt(ANISOBETA)
(RESN for IWT_fill_aniso), with nacelle's RESN for every reflection in the
nacelle file and nacelle's model of it for the others (resn_model), and
ANISOBETA(h) = exp(h^T B0 h + d*^2 h^T B1 h) fitted to nacelle's ANISOBETA
(exact for nacelle's anisotropy model, with or without its resolution
dependence, phasertng.macaniso.resolution_dependent); TEPS is taken as 1.
"""

from mmtbx.refinement.llgi_e_dmodel import SIGMAA_MAX

llgi_map_params = iotbx.phil.parse("""\
  enabled = True
    .type = bool
    .short_caption = Write LLGI map coefficients
    .help = "Write the LLGI map coefficients (IWT, IWT_fill, dIWT; see " \
            "mmtbx.refinement.llgi_maps) to the output MTZ file when " \
            "target=llgi, in place of the usual ones."
  anisotropic_maps = False
    .type = bool
    .short_caption = Also write the LLGI maps with the anisotropy left in
    .help = "Also write IWT_aniso, IWT_fill_aniso and dIWT_aniso: the LLGI " \
            "maps on the data's anisotropic scale (sqrt(TEPS)*RESN), without " \
            "nacelle's anisotropy (ANISOBETA) divided out."
""")

LABELS = group_args(
  bias_reduced="IWT", fill="IWT_fill", difference="dIWT",
  bias_reduced_aniso="IWT_aniso", fill_aniso="IWT_fill_aniso",
  difference_aniso="dIWT_aniso")

def posterior_mean_e(e_obs_sq, sig_e_obs_sq, e_calc, sigmaa, centric_flags):
  """ Exact <E> along the model phase. sig_e_obs_sq <= 0 (no measurement
  error): |E| = sqrt(max(E_obs^2, 0)) and only the phase is uncertain. """
  from cctbx.xray import ext as xray_ext
  n = e_obs_sq.size()
  result = flex.double(n, 0)
  measured = sig_e_obs_sq > 0
  isel = measured.iselection()
  if(isel.size() > 0):
    r = xray_ext.llgi_exact_evaluate(
      e_obs_sq=e_obs_sq.select(isel),
      sig_e_obs_sq=sig_e_obs_sq.select(isel),
      e_calc=e_calc.select(isel),
      sigmaa=sigmaa.select(isel),
      centric_flags=centric_flags.select(isel),
      null_log_z=flex.double(isel.size(), 0))  # LLGI value not needed
    result.set_selected(isel, r.e_expected)
  isel = (~measured).iselection()
  if(isel.size() > 0):
    import scitbx.math
    e2 = e_obs_sq.select(isel)
    e = flex.sqrt(e2.set_selected(e2 < 0, 0))
    sa = sigmaa.select(isel)
    ec = e_calc.select(isel)
    cen = centric_flags.select(isel)
    x = 2 * sa * ec * e / (1 - sa * sa)
    ratio = scitbx.math.bessel_i1_over_i0(x)
    ratio.set_selected(cen, flex.tanh((x / 2).select(cen)))
    result.set_selected(isel, e * ratio)
  return result

def e_scale_coefficients(e_expected, e_calc, sigmaa, phi, centric_flags):
  """ The E-scale coefficient magnitudes (signed, along the model phase),
  as flex.double arrays: bias_reduced, fill and difference. """
  c = (1 - phi) * sigmaa
  difference = e_expected - sigmaa * e_calc
  fill = e_expected + difference
  fill.set_selected(centric_flags, e_expected.select(centric_flags))
  return group_args(
    bias_reduced=(e_expected - c * e_calc) / (1 - c * sigmaa),
    fill=fill,
    difference=difference)

def _unit_phase(f):
  # exp(i*phase); 1 where f = 0
  return flex.polar(flex.double(f.size(), 1), flex.arg(f))

def _indices_as_numpy(indices):
  import numpy as np
  if(hasattr(indices, "as_vec3_double")):
    return indices.as_vec3_double().as_numpy_array()
  return np.array(indices, dtype=float).reshape(-1, 3)

def fit_anisobeta(indices, anisobeta, unit_cell=None):
  """ An anisotropic factor exp(h^T B h) as fitted (least squares) to
  log(anisobeta): the six coefficients of the quadratic form in h, and a
  function evaluating it for other indices. With unit_cell, the form is
  exp(h^T B0 h + d*^2 h^T B1 h) (twelve coefficients, B0's first): nacelle's
  ANISOBETA with resolution-dependent anisotropy, in which the anisotropy
  changes linearly with d*^2 (phasertng.macaniso.resolution_dependent); it
  includes the quadratic form (B1 = 0). """
  import numpy as np
  def design(idx):
    h = _indices_as_numpy(idx)
    quad = np.column_stack([h[:, 0]**2, h[:, 1]**2, h[:, 2]**2,
      h[:, 0]*h[:, 1], h[:, 0]*h[:, 2], h[:, 1]*h[:, 2]])
    if(unit_cell is None): return quad
    g = unit_cell.reciprocal_metrical_matrix()   # (G11, G22, G33, G12, G13, G23)
    dss = quad.dot(np.array([g[0], g[1], g[2], 2*g[3], 2*g[4], 2*g[5]]))
    return np.column_stack([quad, quad * dss[:, None]])
  y = np.log(np.asarray(anisobeta, dtype=float))
  coef = np.linalg.lstsq(design(indices), y, rcond=None)[0]
  rms = float(np.sqrt(np.mean((y - design(indices).dot(coef))**2)))
  def evaluate(idx):
    return flex.double(np.exp(design(idx).dot(coef)))
  return group_args(coefficients=coef, rms_log_residual=rms, evaluate=evaluate)

def fit_k_anisotropic(miller_array):
  """ fmodel's anisotropic scale k_anisotropic (a miller array of positive
  values) as one of the forms mmtbx.bulk_solvent.scaler chooses between:
  exp(-h^T U* h)-type (fit_anisobeta) or the 12-parameter polynomial
  1 + sum a_i p_i(h) of mmtbx.f_model.ext.k_anisotropic(indices, a,
  unit_cell), whichever fits better (in ln). Returns group_args(form,
  rms_log_residual, evaluate(miller_set)). """
  import numpy as np
  uc = miller_array.unit_cell()
  ras, rbs, rcs = uc.reciprocal_parameters()[:3]
  def poly_design(ms):
    h = _indices_as_numpy(ms.indices())
    stol = np.array(ms.sin_theta_over_lambda_sq().data())
    s = np.where(stol != 0, 1. / np.where(stol != 0, stol, 1), 0)
    hh, kk, ll = h[:, 0]**2*ras*ras, h[:, 1]**2*rbs*rbs, h[:, 2]**2*rcs*rcs
    kl, hl = 2*h[:, 1]*h[:, 2]*rbs*rcs, 2*h[:, 0]*h[:, 2]*ras*rcs
    hk = 2*h[:, 0]*h[:, 1]*ras*rbs
    return np.column_stack([hh*s, hh, kk*s, kk, ll*s, ll, kl*s, kl, hl*s, hl,
      hk*s, hk])
  k = np.array(miller_array.data())
  y = np.log(k)
  quad = fit_anisobeta(miller_array.indices(), miller_array.data())
  a = np.linalg.lstsq(poly_design(miller_array), k - 1, rcond=None)[0]
  k_poly = 1 + poly_design(miller_array).dot(a)
  rms_poly = float(np.sqrt(np.mean((np.log(np.clip(k_poly, 1e-12, None)) - y)**2)))
  if(quad.rms_log_residual <= rms_poly):
    return group_args(form="exp", rms_log_residual=quad.rms_log_residual,
      evaluate=lambda ms: quad.evaluate(ms.indices()))
  from mmtbx.f_model import ext as f_model_ext
  def evaluate(ms):
    k = f_model_ext.k_anisotropic(ms.indices(), flex.double(a), uc)
    k.set_selected(k < 1e-12, 1e-12)
    return k
  return group_args(form="polynomial", rms_log_residual=rms_poly,
    evaluate=evaluate)

class resn_model(object):
  """ nacelle's RESN for any reflections. nacelle computes
    RESN^2 = epsilon * SigmaN(s^2) * ANISOBETA(h),
  SigmaN its isotropic fit to the data as a function of resolution
  (smooth, with steps between its resolution bins) and ANISOBETA(h) =
  exp(h^T B0 h + d*^2 h^T B1 h) its anisotropy model (B1 = 0 without
  resolution-dependent anisotropy). So G = RESN^2/(epsilon*ANISOBETA)
  depends on resolution only, and log(ANISOBETA) is exactly of that form
  (fit_anisobeta with the unit cell). RESN is nacelle's value where nacelle
  has one; otherwise from G at the nearest reflection in d*^2 and the
  fitted anisotropy.

  Built from llgi_data.nacelle_all (every reflection in the nacelle file,
  including those excluded from refinement) if present, otherwise from
  llgi_data.resn and .anisobeta (the reflections refined). nacelle_all is
  as read, so it is put on the scale of llgi_data.resn (which
  mmtbx.f_model.manager.apply_scale_k1_to_f_obs rescales) by their ratio
  on the reflections they share. Without ANISOBETA, G = RESN^2/epsilon and
  the model is isotropic. """

  def __init__(self, llgi_data):
    import numpy as np
    source = getattr(llgi_data, "nacelle_all", None)
    if(source is not None):
      resn, anisobeta = source.resn.map_to_asu(), source.anisobeta
      refined = llgi_data.resn.map_to_asu()
      pairs = miller.match_indices(refined.indices(), resn.indices()).pairs()
      if(pairs.size() > 0):
        k = float(np.median((refined.data().select(pairs.column(0)) /
          resn.data().select(pairs.column(1))).as_numpy_array()))
        resn = resn.customized_copy(data=resn.data() * k)
    else:
      resn, anisobeta = llgi_data.resn, getattr(llgi_data, "anisobeta", None)
    resn = resn.map_to_asu()
    keep = resn.data() > 0
    if(anisobeta is not None):
      resn, anisobeta = resn.common_sets(anisobeta.map_to_asu())
      keep = (resn.data() > 0) & (anisobeta.data() > 0)
      anisobeta = anisobeta.select(keep)
    resn = resn.select(keep)
    self.resn = resn
    self.anisotropic = anisobeta is not None
    ab = anisobeta.data() if self.anisotropic else flex.double(resn.size(), 1)
    self.fit = fit_anisobeta(resn.indices(), ab,
      unit_cell=resn.unit_cell()) if self.anisotropic else None
    self._g = flex.pow2(resn.data()) / (
      resn.epsilons().data().as_double() * ab)
    self._dss = resn.d_star_sq().data()

  def g(self, d_star_sq):
    """ G = RESN^2/(epsilon*ANISOBETA) at the nearest reflection in d*^2 """
    from mmtbx.map_tools import nearest_in_d_star_sq
    return nearest_in_d_star_sq(self._dss, self._g, d_star_sq)

  def anisobeta(self, miller_set):
    if(not self.anisotropic): return flex.double(miller_set.size(), 1)
    return self.fit.evaluate(miller_set.map_to_asu().indices())

  def resn_for(self, miller_set):
    """ RESN for the reflections of miller_set (flex.double, in its
    order), and a flex.bool marking nacelle's own values """
    ms = miller_set.map_to_asu()
    result = flex.sqrt(ms.epsilons().data().as_double() *
      self.g(ms.d_star_sq().data()) * self.anisobeta(ms))
    pairs = miller.match_indices(ms.indices(), self.resn.indices()).pairs()
    result.set_selected(pairs.column(0), self.resn.data().select(pairs.column(1)))
    exact = flex.bool(ms.size(), False)
    exact.set_selected(pairs.column(0), True)
    return result, exact

def get_resn_model(llgi_data):
  """ resn_model(llgi_data), kept on llgi_data (._resn_model) while
  llgi_data.resn is the same array (it is replaced when f_obs is
  rescaled, mmtbx.f_model.manager.apply_scale_k1_to_f_obs) """
  cached = getattr(llgi_data, "_resn_model", None)
  if(cached is not None and cached[0] is llgi_data.resn):
    return cached[1]
  result = resn_model(llgi_data)
  llgi_data._resn_model = (llgi_data.resn, result)
  return result

def compute(fmodel, params=None, log=None):
  """ LLGI map coefficients for fmodel (target=llgi, llgi_data with sigmaa
  and the intensities attached). Returns group_args(arrays=[(label,
  complex miller.array)], ...) or None (with a message on log) if the
  data needed are not available. The maps have nacelle's anisotropy
  removed (if there is an ANISOBETA column); with
  params.anisotropic_maps, the _aniso versions are added. """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  import mmtbx.refinement.llgi_phi_table as llgi_phi_table
  if(params is None):
    params = llgi_map_params.extract()
  def say(msg):
    if(log is not None): print("LLGI maps: " + msg, file=log)
  llgi_data = fmodel.llgi_data()
  if(llgi_data is None or getattr(llgi_data, "sigmaa", None) is None):
    say("no LLGI data or sigmaA; not written.")
    return None
  if(getattr(llgi_data, "e_obs_sq", None) is None):
    say("the LLGI data file has no intensities; not written.")
    return None
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  fmnas = fmodel.f_model_no_aniso_scale()
  em = llgi_e_sigmaa.build_e_model(fmnas.data(), epsilons, d_star_sq)
  e_calc = flex.abs(em.e_model)
  phase = _unit_phase(fmnas.data())
  sigmaa = llgi_data.sigmaa.data().deep_copy()
  sigmaa.set_selected(sigmaa < 0, 0)
  sigmaa.set_selected(sigmaa > SIGMAA_MAX, SIGMAA_MAX)
  import mmtbx.refinement.llgi_tncs as llgi_tncs
  tncs_u, tncs_t = llgi_tncs.arrays(llgi_data)
  if(tncs_u is not None):
    # tNCS: the t-normalised frame (llgi.h tncs_frame), in which the data
    # are normalised: sigmaA -> a', Ec -> s*Ec; sqrt(t) in the F scale
    sigmaa, s_frame = llgi_tncs.frame(sigmaa, tncs_u, tncs_t)
    e_calc = e_calc * s_frame
  e_obs_sq = llgi_data.e_obs_sq.data()
  sig = llgi_data.sig_e_obs_sq.data()
  e_expected = posterior_mean_e(e_obs_sq, sig, e_calc, sigmaa, centric)
  phi = llgi_phi_table.phi(sigmaa, sig, centric)
  ce = e_scale_coefficients(e_expected, e_calc, sigmaa, phi, centric)
  # F scales: f_scale_aniso keeps the anisotropy, f_scale removes it
  resn = llgi_data.resn.data()
  teps = llgi_data.teps.data() if tncs_u is None else tncs_t
  f_scale_aniso = flex.sqrt(teps) * resn
  anisobeta = getattr(llgi_data, "anisobeta", None)
  ab = flex.double(resn.size(), 1)
  if(anisobeta is None):
    say("no ANISOBETA column; the maps keep the anisotropy of the data.")
  else:
    ab = anisobeta.data().deep_copy()
    ab_bad = ab <= 0
    if(ab_bad.count(True) > 0):
      say("%d reflection(s) with ANISOBETA <= 0 keep their anisotropy."
        % ab_bad.count(True))
      ab.set_selected(ab_bad, 1.)
  f_scale = f_scale_aniso / flex.sqrt(ab)
  def as_complex(values, scale):
    v = values * scale
    return f_obs.customized_copy(
      data=flex.complex_double(v, flex.double(v.size(), 0)) * phase,
      sigmas=None)
  # the reflections left out: sigmaA*Ec on the F scale (model_missing_
  # reflections_llgi supplies Ec and sigmaA on the complete set; complete_with
  # takes those not in f_obs), with nacelle's RESN (resn_model)
  from mmtbx.map_tools import model_missing_reflections_llgi
  mro = model_missing_reflections_llgi(fmodel=fmodel)
  mro.get_missing(isotropize=False)   # for its e_scale_missing
  mis = mro.e_scale_missing
  sa_m = mis.sigmaa.deep_copy()
  sa_m.set_selected(sa_m < 0, 0)
  sa_m.set_selected(sa_m > SIGMAA_MAX, SIGMAA_MAX)
  ms_m = mis.miller_set
  rm = get_resn_model(llgi_data)
  def fill_array(scale):
    return ms_m.array(
      data=flex.complex_double(sa_m * mis.e_model_abs * scale,
        flex.double(sa_m.size(), 0)) * _unit_phase(mis.f_model_no_aniso))
  # added reflections: the F scale RESN/sqrt(ANISOBETA) = sqrt(epsilon*G)
  # (RESN itself in the _aniso map), as on the observed reflections
  arrays = [
    (LABELS.bias_reduced, as_complex(ce.bias_reduced, f_scale)),
    (LABELS.fill, as_complex(ce.fill, f_scale).complete_with(
      other=fill_array(mis.resn / flex.sqrt(rm.anisobeta(ms_m))), scale=False)),
    (LABELS.difference, as_complex(ce.difference, f_scale))]
  if(params.anisotropic_maps and anisobeta is not None and rm.anisotropic):
    if(rm.fit.rms_log_residual > 1.e-3):
      say("ANISOBETA is not a quadratic form in h (rms residual of the fit "
        "%.3g in ln); the anisotropy of the filled reflections is approximate."
        % rm.fit.rms_log_residual)
    arrays += [
      (LABELS.bias_reduced_aniso, as_complex(ce.bias_reduced, f_scale_aniso)),
      (LABELS.fill_aniso, as_complex(ce.fill, f_scale_aniso).complete_with(
        other=fill_array(mis.resn), scale=False)),
      (LABELS.difference_aniso, as_complex(ce.difference, f_scale_aniso))]
  if(getattr(llgi_data, "from_amplitudes", False)):
    say("WARNING: computed from intensities reconstructed from amplitudes "
      "(approximate). Use the original intensities if at all possible.")
  n_filled = arrays[1][1].size() - f_obs.size()
  say("wrote %s; the fill maps add %d reflections (unmeasured or left out "
    "of the target) as sigmaA*Ec." % (", ".join(a[0] for a in arrays), n_filled))
  return group_args(arrays=arrays, e_expected=e_expected, phi=phi,
    sigmaa=sigmaa, e_calc=e_calc, coefficients_e=ce, f_scale=f_scale,
    f_scale_aniso=f_scale_aniso, n_filled=n_filled)

def _is_two_fofc(map_type):
  import mmtbx
  mn = mmtbx.map_names(map_name_string=map_type)
  return mn.k == 2 and mn.n == -1 and mn.ml_map is True

class filled_posterior_map_server(object):
  """ Stands in for the electron_density_map object emap of the LLGI helper
  (fmodel.electron_density_map_for_model_building): for the real-space
  steps (fft_map with use_all_data=False), the 2mFo-DFc-type maps
  (mmtbx.map_names k = 2, n = -1, ML weights) are the LLGI posterior map
  2<E> - sigmaA*Ec of the working reflections with every other reflection
  to the resolution limit, free or without data, filled with D*Fc
  (mmtbx.map_tools.model_missing_reflections_llgi: sigmaA*Ec, the limit of
  the posterior coefficient as the information goes to 0), i.e. emap's
  map_coefficients with exclude_free_r_reflections and fill_missing.
  Everything else comes from emap, map_coefficients included: filled where
  fill_missing is asked for (ion placement), not otherwise (water picking,
  which looks for atoms the model lacks). """

  def __init__(self, fmodel, emap):
    self.fmodel = fmodel
    self.emap = emap

  def __getattr__(self, name):
    if(name == "emap"): raise AttributeError(name)
    return getattr(self.emap, name)

  def _filled(self, map_type, isotropize=True):
    return self.emap.map_coefficients(map_type=map_type,
      exclude_free_r_reflections=True, fill_missing=True,
      isotropize=isotropize)

  def fft_map(self, map_type, resolution_factor=1/3., symmetry_flags=None,
        other_fft_map=None, use_all_data=True, **kwds):
    if(use_all_data or not _is_two_fofc(map_type)):
      return self.emap.fft_map(map_type=map_type,
        resolution_factor=resolution_factor, symmetry_flags=symmetry_flags,
        other_fft_map=other_fft_map, use_all_data=use_all_data, **kwds)
    assert kwds.get("map_coefficients") is None
    assert kwds.get("acentrics_scale", 2.0) == 2.0
    assert kwds.get("centrics_pre_scale", 1.0) == 1.0
    assert not kwds.get("force_anomalous_flag_false")
    coeffs = self._filled(map_type)
    if(other_fft_map is not None):
      return miller.fft_map(crystal_gridding=other_fft_map,
        fourier_coefficients=coeffs)
    return coeffs.fft_map(resolution_factor=resolution_factor,
      symmetry_flags=symmetry_flags)

class bias_reduced_map_server(object):
  """ Stands in for the electron_density_map object emap of the LLGI helper
  (fmodel.electron_density_map_for_model_building): the 2mFo-DFc-type maps
  (mmtbx.map_names k = 2, n = -1, ML weights) are the bias-reduced map
  IWT (missing reflections 0, whatever fill_missing; the free reflections
  too with exclude_free_r_reflections=True, or use_all_data=False in
  fft_map, as in mmtbx.map_tools); everything else comes from emap, whose
  mFo-DFc form <E> - sigmaA*Ec is the LLGI difference map. """

  def __init__(self, fmodel, emap):
    self.fmodel = fmodel
    self.emap = emap
    self._coeffs = None

  def __getattr__(self, name):
    if(name == "emap"): raise AttributeError(name)
    return getattr(self.emap, name)

  def _bias_reduced(self):
    if(self._coeffs is None):
      from libtbx.utils import Sorry
      result = compute(fmodel=self.fmodel,
        params=getattr(self.fmodel.llgi_data(), "map_params", None))
      if(result is None):
        raise Sorry("LLGI bias-reduced map not available "
          "(llgi_data.real_space_map=bias_reduced); "
          "llgi_data.real_space_map=posterior_2mfo_dfc uses the LLGI "
          "2mFo-DFc map instead.")
      self._coeffs = dict(result.arrays)[LABELS.bias_reduced]
    return self._coeffs

  def _work(self, coeffs):
    return coeffs.select(self.fmodel.arrays.work_sel)

  def map_coefficients(self, map_type, exclude_free_r_reflections=False,
        **kwds):
    if(_is_two_fofc(map_type)):
      coeffs = self._bias_reduced()
      return self._work(coeffs) if exclude_free_r_reflections else coeffs
    return self.emap.map_coefficients(map_type=map_type,
      exclude_free_r_reflections=exclude_free_r_reflections, **kwds)

  def fft_map(self, map_type, resolution_factor=1/3., symmetry_flags=None,
        other_fft_map=None, use_all_data=True, **kwds):
    if(_is_two_fofc(map_type)):
      coeffs = self._bias_reduced()
      if(not use_all_data):
        coeffs = self._work(coeffs)
      if(other_fft_map is not None):
        return miller.fft_map(crystal_gridding=other_fft_map,
          fourier_coefficients=coeffs)
      return coeffs.fft_map(resolution_factor=resolution_factor,
        symmetry_flags=symmetry_flags)
    return self.emap.fft_map(map_type=map_type,
      resolution_factor=resolution_factor, symmetry_flags=symmetry_flags,
      other_fft_map=other_fft_map, use_all_data=use_all_data, **kwds)
