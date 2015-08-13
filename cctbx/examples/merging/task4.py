from __future__ import division
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell
import pickle

from cctbx.examples.merging import intensity_data
def prepare_simulation_with_noise(sim, transmittance,
                                       apply_noise,
                                       ordered_intensities=None,
                                       half_data_flag = 0):
  result = intensity_data()
  result.frame = sim["frame_lookup"]
  result.miller= sim['miller_lookup']
  raw_obs_no_noise = transmittance * sim['observed_intensity']
  if apply_noise:
    import scitbx.random
    from scitbx.random import variate, normal_distribution
         # bernoulli_distribution, gamma_distribution, poisson_distribution
    scitbx.random.set_random_seed(321)
    g = variate(normal_distribution())
    noise = flex.sqrt(raw_obs_no_noise) * g(len(raw_obs_no_noise))
    # adds in Gauss noise to signal
  else:
    noise = flex.double(len(raw_obs_no_noise),0.)

  raw_obs = raw_obs_no_noise + noise

  if half_data_flag in [1,2]:  # apply selection after random numbers have been applied
    half_data_selection = (sim["frame_lookup"]%2)==(half_data_flag%2)
    result.frame  = sim["frame_lookup"].select(half_data_selection)
    result.miller = sim['miller_lookup'].select(half_data_selection)
    raw_obs       = raw_obs.select(half_data_selection)

  mean_signal = flex.mean(raw_obs)

  sigma_obs = flex.sqrt(flex.abs(raw_obs))
  mean_sigma = flex.mean(sigma_obs)
  print "<I> / <sigma>", (mean_signal/ mean_sigma)

  scale_factor = mean_signal/10.

  #most important line; puts input data on a numerically reasonable scale
  result.raw_obs = raw_obs / scale_factor
  scaled_sigma = sigma_obs / scale_factor

  result.exp_var = scaled_sigma * scaled_sigma

  #ordered intensities gets us the unit cell & miller indices to
  # gain a static array of (sin theta over lambda)**2
  if ordered_intensities is not None:
    uc = ordered_intensities.unit_cell()
    stol_sq = flex.double()
    for i in xrange(len(result.miller)):
      this_hkl = ordered_intensities.indices()[result.miller[i]]
      stol_sq_item = uc.stol_sq(this_hkl)
      stol_sq.append(stol_sq_item)
    result.stol_sq = stol_sq
  return result

if __name__=="__main__":

  import cPickle as pickle
  ordered_intensities = pickle.load(open("intensities.pickle","rb"))
  frames = pickle.load(open("frames.pickle","rb"))
  case = 1000
  sim = pickle.load(open("simulated%05d_0.pickle"%case,"rb"))
  print
  #print "accepted %d obs"%(len(sim["miller_lookup"]))
  #print "accepted frames %d"%(len(sim["frame_lookup"]))
  print "accepted obs %d"%(len(sim["observed_intensity"]))
  #print list(sim["frame_lookup"])
  #print list(sim["miller_lookup"])
  #print list(sim["observed_intensity"])

  transmittance = 0.00001
  apply_noise = True
  FSIM = prepare_simulation_with_noise(sim, transmittance=transmittance,
                                       apply_noise=apply_noise)


  print "OK"
