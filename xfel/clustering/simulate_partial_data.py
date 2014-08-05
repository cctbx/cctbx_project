"""
This module contains tools for simulating partial integration data. In
particular, it is intended to help test XFEL data merging tools.
"""

from iotbx import mtz
from cctbx.crystal_orientation import crystal_orientation, basis_type
from scitbx.matrix import sqr
import random
import math
from xfel.Toy_Network.generate_toy_data import ImageNode, Graph


def make_minimal_dict(miller_array, orientation, wavelength=1.0,
                             use_b=True):

  d = {'observations': [miller_array],
       'mapped_predictions': [None], # <<--- Need to have something here! Yikes!
       'pointgroup': None,
       'current_orientation': [orientation],
       'xbeam': None,
       'ybeam': None,
       'wavelength': wavelength}

  return d


def process_mtz(filename):
  mtzfile = mtz.object(filename)
  miller_dict = mtzfile.as_miller_arrays_dict()
  if ('crystal', 'dataset', 'F(ake)obs') in miller_dict:
    return miller_dict[('crystal', 'dataset', 'F(ake)obs')]


def generate_orientations(n, random_seed=None):
  """

  :param n: number of random orientations to create.
  :param random_seed: random seed (optional)
  :return: a list of crystal orientation objects, with random orientations.
  """

  random.seed(random_seed)
  random_numbers = []  # list of uniform random 3-tuples in the interval (0,1)
  for r in range(n * 3):
    random_numbers.append((random.random(), random.random(), random.random()))

  rand_rotations = []
  for ran in random_numbers:
     rand_rotations.append(calc_rotations(ran))

  return rand_rotations


def calc_rotations(x):
  """
  DEPRECATED: Found that this duplicates flex.random_double_r3_rotation_matric()
  ==============================
  Generates random rotation matrices using the method from Arvo, James (1992),
  "Fast random rotation matrices", in David Kirk, Graphics Gems III, San Diego:
  Academic Press, pages 117-120.
  :param x: a 3 element tuple of random variables. For uniform random variables
  about the unit sphere, use range (1,1,1). For small pertubations,
  use (d,1,d) where d < 1.
  :return: a 9-tuple rotation matrix.
  """
  theta = x[0] * 2 * math.pi
  phi = x[1] * 2 * math.pi
  z = x[2] * 2

  r = math.sqrt(z)
  Vx = math.sin(phi) * r
  Vy = math.cos(phi) * r
  Vz = math.sqrt(2.0 - z)

  st = math.sin(theta)
  ct = math.cos(theta)
  Sx = Vx * ct - Vy * st
  Sy = Vx * st + Vy * ct

  return (Vx * Sx - ct, Vx * Sy - st, Vx * Vz,
          Vy * Sx + st, Vy * Sy - ct, Vy * Vz,
          Vz * Sx, Vz * Sy, 1.0 - z)


def make_orientation(uc, rot=None):
  if rot:
    crystal_matrix = sqr(uc.orthogonalization_matrix()).transpose() \
                     * rot.transpose()
  else:
    crystal_matrix = sqr(uc.orthogonalization_matrix()).transpose()

  return crystal_orientation(crystal_matrix, basis_type.direct)


def run(args):

  # 0. Get the miller array!
  mill_array = process_mtz(args.mtzfile[0]).as_intensity_array()

  # 1. Make a list of (random) crystal orientation objects.
  rotations = generate_orientations(args.n, args.r)
  initial_members = []
  for rot in rotations:
    orr = make_orientation(mill_array.unit_cell(), sqr(rot))
    d = make_minimal_dict(mill_array, orr)
    min_node = ImageNode(dicti=d, use_b=False)
    print 'arg'
  # 2. Make (limited) ImageNode objects!


  # 2. Calculate effective partialities for each crystal orientation



  # 3. Generate random scales (Gaussian distribution).

  # 4. Calculate 'corrected' intensities for each orientation

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
  parser.add_argument('-r', type=float,
                      help='Random seed')
  args = parser.parse_args()
  run(args)


