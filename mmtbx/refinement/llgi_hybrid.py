from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from libtbx import group_args
import iotbx.phil

""" Hybrid LLGI: the Rice approximation with moment-matched (Dobs, Eeff)
for most reflections, and the exact likelihood (cctbx/xray/targets/
llgi_exact.h) where the approximation is known to be poor:

  - no 2nd/4th-moment Rice solution exists for the reflection: its
    French-Wilson posterior is more dispersed than any Rice distribution
    with the same <E^2> (nacelle gives these Dobs = 0.01 and Eeff equal to
    the posterior rms; validity is recomputed here from the intensities,
    not inferred from those values),
  - measurement variance not small compared with the square of the model
    variance: 1 - D^2 > rice_kappa (1 - sigmaA^2)^2 (evaluated per call,
    since the effective sigmaA changes during refinement).

The exact likelihood needs the observed intensities on the E^2 scale,
E_obs^2 = I/(RESN^2 TEPS) and sigma(E_obs^2) = SIGI/(RESN^2 TEPS),
carried in llgi_data as e_obs_sq/sig_e_obs_sq (miller arrays, so they are
selected/common_set with everything else), plus the per-reflection
force_exact flags. The C++ hybrid object itself (which caches the
null-hypothesis term per reflection) is built on first use and cached on
the llgi_data instance.
"""

llgi_hybrid_params = iotbx.phil.parse("""\
  enabled = True
    .type = bool
    .short_caption = Hybrid exact/Rice LLGI
    .help = "Evaluate the exact LLGI (numerical integration over the " \
            "observed intensity's error) instead of the Rice " \
            "approximation for reflections where the approximation is " \
            "not accurate enough (see rice_kappa, and reflections with no " \
            "2nd/4th-moment Rice solution). Needs the intensities and sigmas in the LLGI data " \
            "file."
  rice_kappa = 0.1
    .type = float
    .short_caption = Exact LLGI unless the measurement error is small
    .help = "The Rice approximation has variance 1 - D^2 sigmaA^2, of which " \
            "1 - D^2 comes from the measurement error and D^2 (1 - sigmaA^2) " \
            "from the model error. It is exact as the measurement error " \
            "vanishes, and its worst case (a calculated E far from what the " \
            "observation implies) grows as the model error shrinks, so the " \
            "exact likelihood is used where 1 - D^2 > rice_kappa*(1 - " \
            "sigmaA^2)^2. With 0.1 the worst-case Rice error is about 0.15 " \
            "in the LLGI of a reflection, at any sigmaA; 0 evaluates every " \
            "reflection exactly."
  sigma_e_obs_sq_cutoff = 8.5
    .type = float
    .short_caption = Exclude reflections with sigma(E_obs^2) above
    .help = "Reflections with sigma(E_obs^2) above this carry less " \
            "than ~0.01 bits of information and are excluded from " \
            "refinement altogether (replaces info_cutoff when the " \
            "intensities and sigmas are available)."
""")

def e_obs_sq_and_sigma(i_obs, sig_i_obs, resn, teps):
  """ E_obs^2 and sigma(E_obs^2) from intensities and nacelle's RESN/TEPS
  (flex.double arrays, same order). """
  esn = resn * resn * teps
  return i_obs / esn, sig_i_obs / esn

def force_exact_flags(e_obs_sq, sig_e_obs_sq, centric_flags):
  """ Reflections that use the exact LLGI whatever sigmaA is: those with
  no 2nd/4th-moment Rice solution. (Reflections approaching that limit
  have D -> 0 and large Eeff; the rice_kappa rule makes them exact at any
  sigmaA, since 1 - D^2 > rice_kappa*(1 - sigmaA^2)^2 whenever D^2 < 1 -
  rice_kappa.) Returns group_args(force_exact, rice): a flex.bool array
  and the llgi_rice_moments object (whose dsqr/eeff can be compared with
  nacelle's DOBS/FEFF). """
  rice = xray_ext.llgi_rice_moments(
    e_obs_sq=e_obs_sq, sig_e_obs_sq=sig_e_obs_sq, centric_flags=centric_flags)
  invalid = ~rice.valid
  # Reflections without an intensity error estimate always use Rice
  invalid.set_selected(~(sig_e_obs_sq > 0), False)
  return group_args(force_exact=invalid, rice=rice)

RICE_DOBS_NONE = 0.01

def rice_dobs_eeff(rice):
  """ (Dobs, Eeff) from an llgi_rice_moments object, as phasertng.nacelle
  writes them: the moment-matched values, or, with no Rice solution,
  Dobs = RICE_DOBS_NONE (negligible weight in the Rice approximation; the
  exact likelihood is needed) and Eeff = sqrt(<E^2>), the posterior rms,
  so that Feff stays a sensible amplitude. """
  dobs = flex.sqrt(rice.dsqr)
  eeff = rice.eeff.deep_copy()
  invalid = ~rice.valid
  dobs.set_selected(invalid, RICE_DOBS_NONE)
  eeff.set_selected(invalid,
    flex.sqrt(flex.max(rice.mu2.select(invalid), 0)) if invalid.count(True)
    else flex.double())
  return dobs, eeff

def get_exact_data(llgi_data):
  """ The cctbx.xray.llgi_hybrid object for this llgi_data (E_obs^2,
  sigma(E_obs^2) and the cached null-hypothesis terms), whether or not
  the hybrid target is enabled, or None if the intensities are not
  available. The LLGI map coefficients use it for the exact posterior <E>
  of every measured reflection. """
  if(llgi_data is None): return None
  e_obs_sq = getattr(llgi_data, "e_obs_sq", None)
  if(e_obs_sq is None): return None
  result = getattr(llgi_data, "_hybrid_object", None)
  if(result is None):
    params = getattr(llgi_data, "hybrid_params", None)
    if(params is None): params = llgi_hybrid_params.extract()
    result = xray_ext.llgi_hybrid(
      e_obs_sq=e_obs_sq.data(),
      sig_e_obs_sq=llgi_data.sig_e_obs_sq.data(),
      force_exact=llgi_data.force_exact.data(),
      centric_flags=e_obs_sq.centric_flags().data(),
      rice_kappa=params.rice_kappa)
    llgi_data._hybrid_object = result
  return result

def get_hybrid(llgi_data):
  """ The cctbx.xray.llgi_hybrid object for this llgi_data, or None if the
  hybrid is disabled or the intensities are not available. """
  if(llgi_data is None): return None
  params = getattr(llgi_data, "hybrid_params", None)
  if(params is None or not params.enabled): return None
  return get_exact_data(llgi_data)

def _e_scale_ok(e_params):
  # The E-scale quantities normalise Feff by RESN only; with
  # renormalise_e_eff the Rice Eeff is rescaled and E_obs^2 would have to
  # be too, so the exact likelihood is not used there.
  return not (e_params is not None
              and getattr(e_params, "renormalise_e_eff", False))

def get_e_scale_hybrid(llgi_data, e_params):
  """ get_hybrid() for the E-scale fits (see _e_scale_ok). """
  return get_hybrid(llgi_data) if _e_scale_ok(e_params) else None

def get_e_scale_exact_data(llgi_data, e_params):
  """ get_exact_data() for E-scale map coefficients (see _e_scale_ok). """
  return get_exact_data(llgi_data) if _e_scale_ok(e_params) else None

def map_llgi_data(llgi_data, op):
  """ Copy of llgi_data with op applied to every miller-array component
  (e.g. select, common_set); other components are carried over as is,
  except cached objects (names starting with '_'), which depend on the
  index set. """
  if(llgi_data is None): return None
  d = {}
  for k, v in llgi_data.__dict__.items():
    if(k.startswith("_")): continue
    if(v is not None and hasattr(v, "indices") and hasattr(v, "data")):
      d[k] = op(v)
    else:
      d[k] = v
  return group_args(**d)

def replace_llgi_data(llgi_data, **changes):
  """ Copy of llgi_data with some components replaced, keeping everything
  else (including cached objects: the index set is unchanged). """
  d = dict(llgi_data.__dict__)
  d.update(changes)
  result = group_args(**dict((k, v) for k, v in d.items()
    if not k.startswith("_")))
  for k, v in d.items():
    if(k.startswith("_")): setattr(result, k, v)
  return result

AMPLITUDE_WARNING = """\
******************************************************************************
WARNING: the LLGI data were derived from amplitudes, not intensities.

The LLGI likelihood, the hybrid exact/Rice evaluation and the bias-reduced
map coefficients all need the measured intensities and their errors. These
have been reconstructed approximately from the amplitudes%s. Weak
reflections, which French-Wilson treatment changes most, are reconstructed
least accurately (tested on 2G38: within 0.2 sigma(I) for 95%% of
reflections, but up to a few sigma(I) for the weakest).

For best results, run phasertng.nacelle on the original intensities
(I, SIGI) instead, and use that file here.
******************************************************************************"""

def is_french_wilson(f, sigf, centric_flags, max_violation_fraction=0.005):
  """ Same test as phasertng's french_wilson::is_FrenchWilsonF: French-Wilson
  amplitudes are all positive with SIGF/F below its large-sigma limits,
  sqrt(4/pi - 1) = 0.523 (acentric) and sqrt(pi/2 - 1) = 0.756 (centric). """
  if(f.size() == 0): return False
  if((f <= 0).count(True) > 0 or (sigf <= 0).count(True) > 0): return False
  ratio = sigf / f
  if((ratio > 1).count(True) > 0): return False
  limit = flex.double(f.size(), 0.523)
  limit.set_selected(centric_flags, 0.756)
  return (ratio > limit).count(True) <= max_violation_fraction * f.size()

def prior_mean_intensity(f_obs, n_per_bin=200, max_bins=60):
  """ Estimate of the prior <I> a French-Wilson calculation used, as the
  binned mean of F^2 + SIGF^2 (the posterior <J>, whose average is the
  prior mean), interpolated linearly in d*^2. Iterating with the
  reconstructed intensities instead was tested and is unstable for the
  weakest shells. """
  import numpy as np
  d = f_obs.d_star_sq().data().as_numpy_array()
  y = (flex.pow2(f_obs.data()) + flex.pow2(f_obs.sigmas())).as_numpy_array()
  order = np.argsort(d)
  nb = max(1, min(max_bins, len(d)//n_per_bin))
  edges = np.linspace(0, len(d), nb + 1).astype(int)
  xc = [d[order[a:b]].mean() for a, b in zip(edges[:-1], edges[1:])]
  yc = [y[order[a:b]].mean() for a, b in zip(edges[:-1], edges[1:])]
  return flex.double(np.interp(d, xc, yc))

def intensities_from_amplitudes(f_obs, log=None):
  """ Approximate (I, SIGI) from amplitudes with sigmas. French-Wilson
  amplitudes are inverted exactly given the prior <I> (cctbx.xray.
  llgi_french_wilson_inverse), with <I> estimated by prior_mean_intensity;
  other amplitudes are converted as I = F^2, SIGI = 2*F*SIGF. Returns
  group_args(i_obs (miller array, xray intensity), french_wilson (bool),
  n_prior_dominated). """
  f = f_obs.data()
  sigf = f_obs.sigmas()
  assert sigf is not None
  centric = f_obs.centric_flags().data()
  fw = is_french_wilson(f, sigf, centric)
  n_prior = 0
  if(fw):
    inv = xray_ext.llgi_french_wilson_inverse(
      f=f, sigf=sigf, mean_intensity=prior_mean_intensity(f_obs),
      centric_flags=centric)
    i_data, sig_data = inv.i_obs, inv.sig_i_obs
    n_prior = inv.prior_dominated.count(True)
  else:
    i_data = flex.pow2(f)
    sig_data = 2 * f * sigf
  i_obs = f_obs.customized_copy(data=i_data, sigmas=sig_data)
  i_obs.set_observation_type_xray_intensity()
  if(log is not None):
    print(AMPLITUDE_WARNING % (
      " by inverting the French-Wilson calculation" if fw
      else " as I = F^2, SIGI = 2*F*SIGF (the amplitudes do not look like "
           "French-Wilson output)"), file=log)
  return group_args(i_obs=i_obs, french_wilson=fw,
    n_prior_dominated=n_prior)
