"""
This module contains tools for simulating partial integration data. In
particular, it is intended to help test XFEL data merging tools.
"""

from iotbx import mtz
import cctbx.miller
from cctbx.crystal_orientation import crystal_orientation, basis_type
from scitbx.matrix import sqr, col
from xfel.Toy_Network.generate_toy_data import ImageNode
from dials.array_family import flex
from dxtbx.model.detector import detector_factory
import cPickle

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
  detector = detector_factory.simple('SENSOR_UNKNOWN', # Sensor
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
  dict_list = []
  for im in range(args.n):
    A = sqr(flex.random_double_r3_rotation_matrix()) * ortho
    orientation = crystal_orientation(A, basis_type.reciprocal)
    pix_coords, miller_set = get_pix_coords(wavelength, A, mill_array, detector)
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
    node = ImageNode(dicti=temp_dict)
    temp_dict['full_observations'] = [obs]
    temp_dict['observations'] = [obs*node.partialities]

    with(open("simulated_data_{}.pickle".format(im), 'wb')) as pkl:
        cPickle.dump(temp_dict, pkl)
  # 4. Calculate 'corrected' intensities for each orientation
  # i.e find the partialities for each, and 'partialise'

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
  args = parser.parse_args()
  run(args)


