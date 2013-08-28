
from __future__ import division
import libtbx.load_env
from libtbx.test_utils import approx_equal
from cStringIO import StringIO
import os.path

def exercise () :
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_examples/p9-sad/p9.sca",
    test=os.path.isfile)
  if (hkl_file is None) :
    print "phenix_examples/p9-sad/p9.sca not available, skipping test"
    return
  # XXX very hard to avoid cross-imports here if we want to test on actual
  # real data...
  from iotbx import file_reader
  from scitbx.array_family import flex
  hkl_in = file_reader.any_file(hkl_file)
  i_obs = hkl_in.file_server.miller_arrays[0]
  i_obs.setup_binner(n_bins=50)
  meas = i_obs.measurability()
  assert approx_equal(meas, 0.091044)
  meas_binned = i_obs.measurability(use_binning=True)
  bijvoet_ratios = i_obs.bijvoet_ratios(measurable_only=True)
  assert (bijvoet_ratios.size() == 1713)
  assert approx_equal(flex.mean(bijvoet_ratios), 0.476718)
  bijvoet_ratios = i_obs.bijvoet_ratios(measurable_only=False)
  assert (bijvoet_ratios.size() == 18815)
  assert approx_equal(flex.mean(bijvoet_ratios), 0.2727699)
  wilson_mean = i_obs.wilson_plot()
  assert approx_equal(wilson_mean, 18481.1252)
  i_obs.setup_binner(n_bins=10)
  wilson_binned = i_obs.wilson_plot(use_binning=True)
  out = StringIO()
  wilson_binned.show(f=out)
  assert (out.getvalue() == """\
unused:         - 28.4910 [   0/7   ]
bin  1: 28.4910 -  3.7555 [4161/4167] 79994.6693
bin  2:  3.7555 -  2.9819 [4160/4164] 47274.9530
bin  3:  2.9819 -  2.6052 [4175/4177] 17682.1962
bin  4:  2.6052 -  2.3672 [4133/4139] 10444.3793
bin  5:  2.3672 -  2.1976 [4115/4133] 9079.7295
bin  6:  2.1976 -  2.0680 [4180/4210] 6529.7430
bin  7:  2.0680 -  1.9645 [4076/4123] 5050.6608
bin  8:  1.9645 -  1.8790 [4109/4156] 3617.7585
bin  9:  1.8790 -  1.8067 [4115/4157] 2096.0724
bin 10:  1.8067 -  1.7443 [3957/4151] 1472.1231
unused:  1.7443 -         [   0/0   ]
""")
  # TODO lots more stuff - eventually most of the functions used in Xtriage
  # should be tested here
  print "OK"

if (__name__ == "__main__") :
  exercise()
