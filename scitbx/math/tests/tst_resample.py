from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx.array_family import flex
import math
from six.moves import range

def test_resample(seed=0):
  obs_ori=flex.double(range(20))

  npb_draw = scitbx.math.non_parametric_bootstrap( obs_ori, -seed-1 )
  obs = npb_draw.draw( 100 )

  npb = scitbx.math.non_parametric_bootstrap( obs, -seed-2)
  sbs = scitbx.math.smooth_bootstrap( obs, -seed-3 )

  mean_t = flex.mean( obs_ori )
  var_t =  flex.mean( obs_ori*obs_ori ) - mean_t*mean_t

  mean_of_mean = 0
  var_of_mean = 0

  mean_sbs = 0
  std_sbs = 0

  n_sample=1e3
  size = 100.0
  for iteration in range(int(n_sample)):
    sample = npb.draw(int(size))
    sample_2 = sbs.draw(int(size))

    single_mean = flex.mean( sample )
    mean_of_mean += single_mean
    var_of_mean += single_mean*single_mean

    tmp = flex.mean( sample_2 )
    mean_sbs += tmp
    std_sbs += tmp*tmp

  mean_of_mean = mean_of_mean/n_sample
  var_of_mean = var_of_mean/(n_sample) - mean_of_mean*mean_of_mean
  var_of_mean  = math.sqrt( var_of_mean )

  mean_sbs /=n_sample
  std_sbs = (std_sbs/(n_sample))-mean_sbs*mean_sbs
  std_sbs = math.sqrt( std_sbs )

  assert math.fabs(9.5-mean_of_mean)/var_of_mean < 4
  assert math.fabs(9.5-mean_sbs)/std_sbs < 4


def run():
  for ii in range(10):
    test_resample(ii)
  print("OK")

if (__name__ == "__main__"):
  run()
