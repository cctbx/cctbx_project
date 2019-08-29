# LIBTBX_SET_DISPATCHER_NAME cluster.simulate_partial_data
"""
This module contains tools for simulating partial integration data. In
particular, it is intended to help test XFEL data merging tools.
"""
from __future__ import absolute_import, division, print_function

from iotbx import mtz
import cctbx.miller
from cctbx.crystal_orientation import crystal_orientation, basis_type
from scitbx.matrix import sqr, col
#from xfel.Toy_Network.generate_toy_data import ImageNode
from xfel.clustering.singleframe import ImageNode
from dials.array_family import flex
from dxtbx.model import DetectorFactory

import random
from six.moves import cPickle as pickle
import logging
from six.moves import range

eps = 0.001  # Tolerance for assertions
p_threshold = 0.1  # Partiality threshold for inclusion

def process_mtz(filename):
  mtzfile = mtz.object(filename)
  miller_dict = mtzfile.as_miller_arrays_dict()
  if ('crystal', 'dataset', 'F(ake)obs') in miller_dict:
    return miller_dict[('crystal', 'dataset', 'F(ake)obs')]

def get_pix_coords(wavelength, A, mill_arr, detector, delta_i=0.02):
    """ Code copied from sim.py courtesy of Aaron and Tara """
    s0=col((0,0,-1/wavelength))
    q=flex.vec3_double([A*col(idx) for idx in  mill_arr.indices().as_vec3_double()])
    s0_hat=flex.vec3_double([s0.normalize()]*len(q))
    q_hat=q.each_normalize()
    #q_hat.cross(flex.vec3_double([s0_hat]*len(q_hat)))
    e1_hat = q_hat.cross(s0_hat)
    c0_hat = s0_hat.cross(e1_hat)
    q_len_sq = flex.double([col(v).length_sq() for v in q])
    a_side=q_len_sq*wavelength/2
    b_side=flex.sqrt(q_len_sq)-a_side**2
    #flex.vec3_double([sqrt(q.length_sq()-a_side**2 for idx in mill_arr)])
    r_vec=flex.vec3_double(-a_side*s0_hat+b_side*c0_hat)
    s1=r_vec+s0

    EQ=q+s0
    len_EQ=flex.double([col(v).length() for v in EQ])
    ratio=len_EQ*wavelength

    indices = flex.miller_index()
    coords =flex.vec2_double()
    for i in range(len(s1)):
        if ratio[i] > 1 - delta_i and ratio[i] < 1 + delta_i:
            indices.append(mill_arr.indices()[i])
            pix = detector[0].get_ray_intersection_px(s1[i])
            if detector[0].is_coord_valid(pix):
                coords.append(pix)

    return coords, indices

def run(args):

  distance = 125
  centre = (97.075, 97.075)
  pix_size = (0.11, 0.11)
  image_size = (1765, 1765)
  wavelength = args.w or 1.0
  # 1. Make a dummy detector
  detector = DetectorFactory.simple('SENSOR_UNKNOWN', # Sensor
                                          distance,
                                          centre,
                                          '+x','-y', # fast/slow direction
                                          pix_size,
                                          image_size)

  # 2. Get the miller array!
  mill_array = process_mtz(args.mtzfile[0]).as_intensity_array()
  ortho = sqr(mill_array.crystal_symmetry().unit_cell().reciprocal() \
            .orthogonalization_matrix())

  # 3.Create some image_pickle dictionairies that contain 'full' intensities,
  # but are otherwise complete.
  im = 0
  while im < args.n:
    im += 1
    A = sqr(flex.random_double_r3_rotation_matrix()) * ortho
    orientation = crystal_orientation(A, basis_type.reciprocal)
    pix_coords, miller_set = get_pix_coords(wavelength, A, mill_array, detector)
    if len(miller_set) > 10:  # at least 10 reflections
        miller_set = cctbx.miller.set(mill_array.crystal_symmetry(), miller_set,
                anomalous_flag=False)
        obs = mill_array.common_set(miller_set)
        temp_dict = {'observations': [obs],
                     'mapped_predictions': [pix_coords],
                     'pointgroup': None,
                     'current_orientation': [orientation],
                     'xbeam': centre[0],
                     'ybeam': centre[1],
                     'wavelength': wavelength}
        old_node = ImageNode(dicti=temp_dict, scale=False)
        # Remove all reflection that are not at least p partial
        partial_sel = (old_node.partialities > p_threshold)

        temp_dict['full_observations'] = [obs.select(partial_sel)]
        temp_dict['observations'] = [obs.select(partial_sel)
                            * old_node.partialities.select(partial_sel)]
        temp_dict['mapped_predictions'] = \
                    [temp_dict['mapped_predictions'][0].select(partial_sel)]

        if logging.Logger.root.level <= logging.DEBUG:  # debug!
          before = temp_dict['full_observations'][0]
          after = temp_dict['observations'][0] / old_node.partialities
          assert sum(abs(before.data() - after.data())) < eps

        if args.r:
          partials = list(temp_dict['observations'][0].data())
          jiggled_partials = flex.double([random.gauss(obs, args.r * obs)
                                          for obs in partials])
          temp_dict['observations'][0] = temp_dict['observations'][0] \
                                      .customized_copy(data=jiggled_partials)

        pkl_name = "simulated_data_{0:04d}.pickle".format(im)
        with(open(pkl_name, 'wb')) as pkl:
          pickle.dump(temp_dict, pkl)

        ''' Only works with no noise:
        if logging.Logger.root.level <= logging.DEBUG:  # debug!
          new_node = ImageNode(pkl_name, scale=False)
          assert sum(old_node.partialities != new_node.partialities) == 0
          new_node.G = 1
          assert sum(new_node.miller_array.indices() !=
                     old_node.miller_array.indices()) == 0

          old_arr = old_node.miller_array
          new_arr = new_node.miller_array / (new_node.partialities *
                                                    new_node.scales)
          assert sum(old_arr.indices() != new_arr.indices()) == 0
          assert sum(abs(old_arr.data() - new_arr.data())) < eps, \
            "delta is {}".format(sum(abs(old_arr.data() - new_arr.data())) )
          assert new_node.G == 1
          assert new_node.minus_2B == 0
          assert len(set(list(new_node.scales))) == 1, "Scales: {}"\
                      .format(new_node.scales)
          assert old_node.G == 1
          assert old_node.minus_2B == 0
          assert len(set(list(old_node.scales))) == 1
          # Testing more:
          d_spacings = list(new_node.miller_array.d_spacings().data())
          miller_indeces = list(new_node.miller_array.indices())
          miller_intensities = list(new_node.miller_array.data()
                                    / (new_node.partialities * new_node.scales))

          for observation in zip(miller_indeces, miller_intensities):
            try:
              test_dict[observation[0]].append(observation[1])
            except KeyError:
              test_dict[observation[0]] = [observation[1]]

  if logging.Logger.root.level <= logging.DEBUG:  # debug!
    for miller in test_dict:
      assert max(test_dict[miller]) - min(test_dict[miller]) < eps,  \
        "Miller index {} were not all equal: {}".format(miller,
                                                        test_dict[miller])
        '''


  # 5. Optional: add some random scale to the images.

  # Optional. Perturb the orientation matrices a bit. (Gaussian, s.d. 0.05deg?)

  # 6. Make a Graph. :D


if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description=('Generate still image data from'
                                                ' an mtz file made by '
                                                'phenix.fake_f_obs.'))
  parser.add_argument('mtzfile', type=str, nargs=1,
                      help='MTZ filename.')
  parser.add_argument('-n', type=int, default=500,
                      help='Number of pickles to generate')
  parser.add_argument('-w', type=float, default=1,
                      help='wavelength of simulated data.')
  parser.add_argument('-r', type=float, default=None,
                      help='Random noise to be applied to partials. Parameter '
                      'is the standard deviation of normally distributed'
                      'noise as a fraction of the partial intensity.')
  args = parser.parse_args()
  run(args)
