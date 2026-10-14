from __future__ import absolute_import, division, print_function
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
fill coefficient tends to; there the F scale is the data's isotropic falloff
RESN/sqrt(epsilon*ANISOBETA) at the nearest observed resolution, times
sqrt(epsilon), and for IWT_fill_aniso also sqrt(ANISOBETA(h)), with
ANISOBETA(h) = exp(h^T B h) fitted to nacelle's ANISOBETA (an exact
quadratic form in h); TEPS is taken as 1.
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

def fit_anisobeta(indices, anisobeta):
  """ nacelle's ANISOBETA as exp(h^T B h): the six coefficients of the
  quadratic form in h fitted (least squares) to log(ANISOBETA), and a
  function evaluating it for other indices. """
  import numpy as np
  def design(idx):
    h = np.array(idx, dtype=float).reshape(-1, 3)
    return np.column_stack([h[:, 0]**2, h[:, 1]**2, h[:, 2]**2,
      h[:, 0]*h[:, 1], h[:, 0]*h[:, 2], h[:, 1]*h[:, 2]])
  y = np.log(np.asarray(anisobeta, dtype=float))
  coef = np.linalg.lstsq(design(indices), y, rcond=None)[0]
  rms = float(np.sqrt(np.mean((y - design(indices).dot(coef))**2)))
  def evaluate(idx):
    return flex.double(np.exp(design(idx).dot(coef)).tolist())
  return group_args(coefficients=coef, rms_log_residual=rms, evaluate=evaluate)

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
  if(anisobeta is None):
    say("no ANISOBETA column; the maps keep the anisotropy of the data.")
    f_scale = f_scale_aniso
  else:
    f_scale = f_scale_aniso / flex.sqrt(anisobeta.data())
  def as_complex(values, scale):
    v = values * scale
    return f_obs.customized_copy(
      data=flex.complex_double(v, flex.double(v.size(), 0)) * phase,
      sigmas=None)
  # the reflections left out: sigmaA*Ec on the F scale (model_missing_
  # reflections_llgi supplies Ec and sigmaA on the complete set; complete_with
  # takes those not in f_obs)
  from mmtbx.map_tools import model_missing_reflections_llgi
  mro = model_missing_reflections_llgi(fmodel=fmodel,
    coeffs=as_complex(ce.fill, f_scale_aniso))
  mro.get_missing()
  mis = mro.e_scale_missing
  sa_m = mis.sigmaa.deep_copy()
  sa_m.set_selected(sa_m < 0, 0)
  sa_m.set_selected(sa_m > SIGMAA_MAX, SIGMAA_MAX)
  ms_m = mis.miller_set
  eps_m = ms_m.epsilons().data().as_double()
  # the data's isotropic falloff RESN/sqrt(epsilon*ANISOBETA) (smooth in
  # resolution) at the nearest observed resolution
  ab = anisobeta.data() if anisobeta is not None else flex.double(resn.size(), 1)
  iso_m = mro.nearest_observed(ms_m.d_star_sq().data(),
    resn / flex.sqrt(epsilons * ab))
  if(anisobeta is not None):
    fit = fit_anisobeta(f_obs.indices(), anisobeta.data())
    if(fit.rms_log_residual > 1.e-3):
      say("ANISOBETA is not a quadratic form in h (rms residual of the fit "
        "%.3g in ln); the anisotropy of the filled reflections is approximate."
        % fit.rms_log_residual)
    ab_m = fit.evaluate(ms_m.indices())
  else:
    ab_m = flex.double(ms_m.size(), 1)
  def fill_array(scale):
    return ms_m.customized_copy(
      data=flex.complex_double(sa_m * mis.e_model_abs * scale,
        flex.double(sa_m.size(), 0)) * _unit_phase(mis.f_model_no_aniso))
  # added reflections: the isotropic falloff (times sqrt(ANISOBETA(h)) in the
  # _aniso map, as on its observed reflections)
  arrays = [
    (LABELS.bias_reduced, as_complex(ce.bias_reduced, f_scale)),
    (LABELS.fill, as_complex(ce.fill, f_scale).complete_with(
      other=fill_array(iso_m * flex.sqrt(eps_m)), scale=False)),
    (LABELS.difference, as_complex(ce.difference, f_scale))]
  if(params.anisotropic_maps and anisobeta is not None):
    arrays += [
      (LABELS.bias_reduced_aniso, as_complex(ce.bias_reduced, f_scale_aniso)),
      (LABELS.fill_aniso, as_complex(ce.fill, f_scale_aniso).complete_with(
        other=fill_array(iso_m * flex.sqrt(eps_m * ab_m)), scale=False)),
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
