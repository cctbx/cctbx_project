from __future__ import division
from mmtbx.secondary_structure import build as ssb
import time

def exercise_process_params():
  ppars = ssb.process_params(params=None)
  assert ppars.enabled
  assert not ppars.restrain_torsion_angles
  #
  params = ssb.master_phil.fetch().extract()
  ppars = ssb.process_params(params=params)
  assert not ppars.enabled
  assert not ppars.restrain_torsion_angles
  params.model_idealization.enabled = True
  ppars = ssb.process_params(params=params)
  assert ppars.enabled
  ppars = ssb.process_params(params=params.model_idealization)
  assert ppars.enabled
  params.model_idealization.sigma_on_reference_non_ss = -1
  try:
    ppars = ssb.process_params(params=params.model_idealization)
  except AssertionError as e:
    assert e.args[0] == 'Bad sigma_on_reference_non_ss parameter'

if (__name__ == "__main__"):
  t0=time.time()
  exercise_process_params()
  print "Time: %6.4f"%(time.time()-t0)
  print "OK"
