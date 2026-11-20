from __future__ import absolute_import, division, print_function
from cctbx_xray_observations_ext import *

def customized_copy(obs, space_group, twin_fractions=None, twin_components=None):
  """ Creates a copy of the observation object and with new twin fractions
  and twin components
  """
  twf = ()
  if twin_fractions is not None:
    twf = twin_fractions
  twc = ()
  if twin_components is not None:
    twc = twin_components
  result = observations(obs, space_group, twf, twc)
  result.fo_sq = obs.fo_sq
  result.ref_twin_fractions = twin_fractions
  result.ref_twin_components = twin_components
  return result

#convenience method for HKLF5 format to get a full miller index set unique under symmetry
def full_set(obs):
  fc_set = []
  for i,h in enumerate(obs.indices):
    fc_set.append(h)
    itr = obs.iterator(i)
    while itr.has_next():
      fc_set.append(itr.next().h)
  fc_set = miller.set(crystal_symmetry=obs.fo_sq.crystal_symmetry(),
    indices=flex.miller_index(fc_set),
    anomalous_flag=obs.fo_sq.anomalous_flag())
  return fc_set.select(fc_set.unique_under_symmetry_selection())
