from cctbx_xray_observations_ext import *

def custom_copy(obs, twin_fractions=None, twin_components=None):
  """ Creates a copy of the observation object and with new twin fractions
  and twin components
  """
  twf = ()
  if twin_fractions is not None:
    twf = twin_fractions
  twc = ()
  if twin_components is not None:
    twc = twin_components
  result = observations(obs, twf, twc)
  result.fo_sq = obs.fo_sq
  result.ref_twin_fractions = twin_fractions
  result.ref_twin_components = twin_components
  return result
