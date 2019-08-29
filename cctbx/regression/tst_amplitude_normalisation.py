from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import group_args
import math
from six.moves import range

class Empty(object):
  pass

def as_normalised_array(miller_array,
                        asu_contents,
                        wilson_plot=None):
  """Old python code replaced by faster C++ code."""
  from cctbx import statistics
  from cctbx import eltbx
  if not wilson_plot:
    wilson_plot = statistics.wilson_plot(miller_array, asu_contents)

  # cache scattering factor info
  gaussians = {}
  for chemical_type in asu_contents:
    gaussians.setdefault(chemical_type, eltbx.xray_scattering.wk1995(
      chemical_type).fetch())

  stol_sq = miller_array.sin_theta_over_lambda_sq()
  epsilons = miller_array.epsilons()
  e_sq_minus_1 = 0
  n_e_greater_than_2 = 0

  normalised_f_obs = flex.double()
  space_group = miller_array.space_group()
  tau = space_group.n_ltr()
  for i in range(0,miller_array.size()):
    s_sq = stol_sq.data()[i]
    f_sq = math.pow(miller_array.data()[i], 2)
    epsilon = epsilons.data()[i]

    sum_fj_sq = 0
    for chemical_type, n_atoms in asu_contents.items():
      n_atoms *= space_group.order_z()
      f0 = gaussians[chemical_type].at_stol_sq(s_sq)
      sum_fj_sq += f0 * f0 * n_atoms

    e_sq = f_sq\
         /(wilson_plot.wilson_intensity_scale_factor*math.exp(-2*wilson_plot.wilson_b*s_sq)
           *epsilon
           *tau
           *sum_fj_sq)
    normalised_f_obs.append(math.sqrt(e_sq))
    e_sq_minus_1 += abs(e_sq - 1)
    if (e_sq > 4.0): n_e_greater_than_2 += 1

  r = Empty()
  r.array = miller.array(
    miller_set=miller.set(
      crystal_symmetry=miller_array.crystal_symmetry(),
      indices=miller_array.indices()).auto_anomalous(),
    data=normalised_f_obs,
    sigmas=miller_array.sigmas())
  r.mean_e_sq_minus_1 = e_sq_minus_1/r.array.size()
  r.percent_e_sq_gt_2 = (100.0*n_e_greater_than_2)/r.array.size()

  return r

def exercise_normalised_amplitudes():
  f = {'C':10,'H':10,'N':1}
  k = 1.322
  b = 5.672
  u = uctbx.unit_cell((3,4,5,90,81,90))
  sgi = sgtbx.space_group_info("P 21")
  i = flex.miller_index(((1,-2,3), (-3,5,0), (1,0,0)))
  d = flex.double((1,2,0.1))
  s = flex.double((2,3,1.1))
  expected = flex.double(
    (0.24360789372276667, 10.918685276828379, 0.0052268772340737469))

  ms = miller.set(
      crystal_symmetry=crystal.symmetry(
        space_group_info=sgi,
        unit_cell=u),
      indices=i)
  ma = miller.array(miller_set=ms, data=d, sigmas=s)
  ma.set_observation_type_xray_amplitude()
  normalised = ma.normalised_amplitudes(
    asu_contents=f,
    wilson_plot=group_args(wilson_intensity_scale_factor=k, wilson_b=b))

  python_normalised = as_normalised_array(
    miller_array=ma,
    asu_contents=f,
    wilson_plot=group_args(wilson_intensity_scale_factor=k, wilson_b=b))

  assert normalised.array().size() == ma.size()
  assert approx_equal(normalised.array().data(),
                      python_normalised.array.data())
  assert approx_equal(normalised.array().data(), expected)

  normalisations = ms.amplitude_normalisations(
    asu_contents=f,
    wilson_plot=group_args(wilson_intensity_scale_factor=k, wilson_b=b))
  assert approx_equal((normalised.array()*normalisations).data(), d)

def run():
  exercise_normalised_amplitudes()
  print("OK")

if __name__ == '__main__':
  run()
