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
  form = rice *hybrid exact
    .type = choice
    .short_caption = Form of the LLGI likelihood
    .help = "rice: the Rice approximation with moment-matched (Dobs, Feff) " \
            "for every reflection, in the target and in the sigmaA and " \
            "ScatFrac fits. hybrid: the exact likelihood (numerical " \
            "integration over the observed intensity's error) where the " \
            "Rice approximation is not accurate enough (see rice_kappa, " \
            "and reflections with no 2nd/4th-moment Rice solution), Rice " \
            "elsewhere; the sigmaA and ScatFrac fits use the exact form " \
            "throughout. exact: the exact likelihood for every reflection " \
            "(slower). hybrid and exact need intensities with sigmas; " \
            "without them the Rice form is used. The map coefficients use " \
            "the exact posterior whatever the form."
  rice_kappa = 0.1
    .type = float
    .short_caption = Hybrid: exact LLGI unless the measurement error is small
    .help = "For form=hybrid. The Rice approximation has variance 1 - D^2 " \
            "sigmaA^2, of which 1 - D^2 comes from the measurement error " \
            "and D^2 (1 - sigmaA^2) from the model error. It is exact as " \
            "the measurement error vanishes, and its worst case (a " \
            "calculated E far from what the observation implies) grows as " \
            "the model error shrinks, so the exact likelihood is used where " \
            "1 - D^2 > rice_kappa*(1 - sigmaA^2)^2. With 0.1 the worst-case " \
            "Rice error is about 0.15 in the LLGI of a reflection, at any " \
            "sigmaA."
""")

def effective_rice_kappa(params):
  """ The rice_kappa of the cctbx.xray.llgi_hybrid object for the likelihood
  form of params (extracted llgi_hybrid_params): 0 (every reflection
  exact) for form=exact, params.rice_kappa otherwise. """
  if(params.form == "exact"): return 0.0
  return params.rice_kappa

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
  mu2 = rice.mu2.select(invalid)
  mu2.set_selected(mu2 < 0, 0)
  eeff.set_selected(invalid, flex.sqrt(mu2))
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
      rice_kappa=effective_rice_kappa(params))
    llgi_data._hybrid_object = result
  return result

def get_hybrid(llgi_data):
  """ The cctbx.xray.llgi_hybrid object for this llgi_data, or None if the
  likelihood form is rice or the intensities are not available. """
  if(llgi_data is None): return None
  params = getattr(llgi_data, "hybrid_params", None)
  if(params is None or params.form == "rice"): return None
  return get_exact_data(llgi_data)

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
