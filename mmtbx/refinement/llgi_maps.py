from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import group_args
import iotbx.phil

""" LLGI map coefficients from the exact posterior of E (map-coefficient
handoff, Oct 2026). On the E scale, with Ec = |Emodel| along the model
phase, sigmaA (clamped to SIGMAA_MAX) and the information fraction
phi(sigmaA, sigma(E_obs^2)) (mmtbx.refinement.llgi_phi_table):

  filled         <E>                                    (missing: sigmaA*Ec)
  bias-reduced   M = (<E> - c*Ec)/(1 - c*sigmaA),  c = (1 - phi)*sigmaA
                                                        (missing: 0)
  difference     <E> - sigmaA*Ec  (proportional to dLLGI/dEc; missing: 0)

<E> is the exact posterior mean along the model phase (cctbx.xray.
llgi_exact_evaluate) for every reflection, whatever form the target uses.
Its model bias is sigmaA*(1 - phi); M subtracts exactly that much of Ec,
and the 1/(1 - c*sigmaA) scale is the handoff's "proposed" normalisation.

F scale: multiply by sqrt(TEPS)*RESN, the data's own normalisation, with
nacelle's anisotropic factor divided out (RESN/sqrt(ANISOBETA): the
isotropic overall falloff is kept, the anisotropy removed) unless
remove_anisotropy=False.
"""

SIGMAA_MAX = 0.995

llgi_map_params = iotbx.phil.parse("""\
  enabled = True
    .type = bool
    .short_caption = Write LLGI map coefficients
    .help = "Write the LLGI map coefficients (filled, bias-reduced and " \
            "difference; see mmtbx.refinement.llgi_maps) to the output " \
            "MTZ file when target=llgi, in addition to the usual ones."
  remove_anisotropy = True
    .type = bool
    .short_caption = Remove anisotropy from F-scale LLGI maps
    .help = "Put F-scale LLGI map coefficients on the isotropic overall " \
            "falloff (RESN/sqrt(ANISOBETA)) rather than the anisotropic " \
            "one (RESN). Needs the ANISOBETA column."
  e_scale = True
    .type = bool
    .short_caption = Also write E-scale bias-reduced and difference maps
""")

LABELS = group_args(
  filled="LLGI_FILLED", bias_reduced="LLGI_BIASRED", difference="LLGI_DIFF",
  bias_reduced_e="LLGI_BIASRED_E", difference_e="LLGI_DIFF_E")

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

def e_scale_coefficients(e_expected, e_calc, sigmaa, phi):
  """ The three E-scale coefficient magnitudes (signed, along the model
  phase), as flex.double arrays. """
  c = (1 - phi) * sigmaa
  return group_args(
    filled=e_expected,
    bias_reduced=(e_expected - c * e_calc) / (1 - c * sigmaa),
    difference=e_expected - sigmaa * e_calc)

def _unit_phase(f):
  # exp(i*phase); 1 where f = 0
  return flex.polar(flex.double(f.size(), 1), flex.arg(f))

def compute(fmodel, params=None, log=None):
  """ LLGI map coefficients for fmodel (target=llgi, llgi_data with sigmaa
  and the intensities attached). Returns group_args(arrays=[(label,
  complex miller.array)], ...) or None (with a message on log) if the
  data needed are not available. """
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
  e_obs_sq = llgi_data.e_obs_sq.data()
  sig = llgi_data.sig_e_obs_sq.data()
  e_expected = posterior_mean_e(e_obs_sq, sig, e_calc, sigmaa, centric)
  phi = llgi_phi_table.phi(sigmaa, sig, centric)
  ce = e_scale_coefficients(e_expected, e_calc, sigmaa, phi)
  # F scale
  resn = llgi_data.resn.data()
  teps = llgi_data.teps.data()
  anisobeta = getattr(llgi_data, "anisobeta", None)
  remove_aniso = params.remove_anisotropy and anisobeta is not None
  if(params.remove_anisotropy and anisobeta is None):
    say("no ANISOBETA column; F-scale maps keep the anisotropy.")
  aniso = anisobeta.data() if remove_aniso else flex.double(resn.size(), 1)
  f_scale = flex.sqrt(teps) * resn / flex.sqrt(aniso)
  def as_complex(values, scale=None):
    v = values if scale is None else values * scale
    return f_obs.customized_copy(
      data=flex.complex_double(v, flex.double(v.size(), 0)) * phase,
      sigmas=None)
  arrays = []
  filled_obs = as_complex(ce.filled, f_scale)
  # Filled: unmeasured reflections get sigmaA*Ec, the posterior mean with
  # no data (model_missing_reflections_llgi supplies Ec, sigmaA and RESN
  # there, from the nearest observed resolution).
  from mmtbx.map_tools import model_missing_reflections_llgi
  mro = model_missing_reflections_llgi(fmodel=fmodel, coeffs=filled_obs)
  mro.get_missing()
  mis = mro.e_scale_missing
  sa_m = mis.sigmaa.deep_copy()
  sa_m.set_selected(sa_m < 0, 0)
  sa_m.set_selected(sa_m > SIGMAA_MAX, SIGMAA_MAX)
  scale_m = mis.resn
  if(remove_aniso):
    scale_m = scale_m / flex.sqrt(mro.nearest_observed(
      mis.miller_set.d_star_sq().data(), aniso))
  fill = mis.miller_set.customized_copy(
    data=flex.complex_double(sa_m * mis.e_model_abs * scale_m,
      flex.double(sa_m.size(), 0)) * _unit_phase(mis.f_model_no_aniso))
  filled = filled_obs.complete_with(other=fill, scale=False)
  arrays.append((LABELS.filled, filled))
  arrays.append((LABELS.bias_reduced, as_complex(ce.bias_reduced, f_scale)))
  arrays.append((LABELS.difference, as_complex(ce.difference, f_scale)))
  if(params.e_scale):
    arrays.append((LABELS.bias_reduced_e, as_complex(ce.bias_reduced)))
    arrays.append((LABELS.difference_e, as_complex(ce.difference)))
  if(getattr(llgi_data, "from_amplitudes", False)):
    say("WARNING: computed from intensities reconstructed from amplitudes "
      "(approximate). Use a nacelle file made from intensities if at all "
      "possible.")
  say("wrote %s (F scale%s)%s." % (", ".join(a[0] for a in arrays),
    ", anisotropy removed" if remove_aniso else "",
    "; filled map completed with %d unmeasured reflections" % (
      filled.size() - f_obs.size())))
  return group_args(arrays=arrays, e_expected=e_expected, phi=phi,
    sigmaa=sigmaa, e_calc=e_calc, coefficients_e=ce, f_scale=f_scale,
    n_filled=filled.size() - f_obs.size())
