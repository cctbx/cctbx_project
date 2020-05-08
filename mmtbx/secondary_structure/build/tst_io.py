from __future__ import absolute_import, division, print_function
from mmtbx.secondary_structure.build import ss_idealization as ssb
from mmtbx.regression import model_1yjp
import iotbx.pdb
import mmtbx.model
import time

def exercise_process_params():
  pdb_inp = iotbx.pdb.input(lines=model_1yjp, source_info='model_1yjp')
  model = mmtbx.model.manager(model_input=pdb_inp)
  ssb_obj = ssb.substitute_ss(model)
  ppars = ssb_obj._process_params(params=None)
  assert ppars.enabled
  assert not ppars.restrain_torsion_angles
  #
  params = ssb.master_phil.fetch().extract()
  ppars = ssb_obj._process_params(params=params)
  assert not ppars.enabled
  assert not ppars.restrain_torsion_angles
  params.ss_idealization.enabled = True
  ppars = ssb_obj._process_params(params=params)
  assert ppars.enabled
  ppars = ssb_obj._process_params(params=params.ss_idealization)
  assert ppars.enabled
  params.ss_idealization.sigma_on_reference_non_ss = -1
  try:
    ppars = ssb_obj._process_params(params=params.ss_idealization)
  except AssertionError as e:
    assert e.args[0] == 'Bad sigma_on_reference_non_ss parameter'

if (__name__ == "__main__"):
  t0=time.time()
  exercise_process_params()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
