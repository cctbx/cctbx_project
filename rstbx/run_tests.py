from libtbx import test_utils
import libtbx.load_env

def run():
  from rstbx.indexing.tst_auto_monoscan import test_automatic_monoscan
  from rstbx.indexing.tst_dataset1 import test_case_obs_data
  from rstbx.indexing.tst_dataset1 import test_case_synthetic_data

  from rstbx.dps_core.sampling import test_simple_sampler

  for functional_test in [test_automatic_monoscan,
                          test_case_obs_data,test_case_synthetic_data]:
    dps,groups = functional_test(verbose=False)
    assert groups[0].reference_lookup_symbol() == "F m -3 m"
  assert test_simple_sampler()
  print "OK"

if (__name__ == "__main__"):
  run()
