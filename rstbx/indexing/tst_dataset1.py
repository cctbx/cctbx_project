from __future__ import absolute_import, division, print_function
import math
import libtbx.load_env
from libtbx.test_utils import approx_equal
from scitbx import matrix
from rstbx.array_family import flex
from rstbx.indexing.do_index import do_index

# reciprocal space vectors expressed internally in units of inverse Angstroms:
# eV per KeV *
# Joules per eV /
# Joule-seconds (Planck's const) /
# meters/second (speed of light) /
# Angstroms per meter
keV_to_inv_Angstrom = 1E3             * \
                      1.602176487E-19 / \
                      6.62606896E-34  / \
                      299792458.E0    / \
                      1E10

def parse_input(filename):
  with open(filename,"r") as G:
    lines = G.readlines()
  reciprocal_vectors = flex.vec3_double()
  for line in lines[0:len(lines)-0]:
    tokens = line.strip().split('\t')
    assert len(tokens)==7
    q_vector = matrix.col([float(i) for i in tokens[4:7]])
    checklength = q_vector.dot(q_vector)
    assert approx_equal(checklength, 1.0, eps=0.001)
    energy_keV = float(tokens[2])
    theta = (math.pi/180.) * float(tokens[3])
    inv_lambda_Angstrom = keV_to_inv_Angstrom * energy_keV
    inv_d_Angstrom = 2. * math.sin(theta) * inv_lambda_Angstrom
    extended_q_Angstrom = inv_d_Angstrom * q_vector
    E=extended_q_Angstrom
    reciprocal_vectors.append(extended_q_Angstrom.elems)
  return reciprocal_vectors

def parse_synthetic(filename):
  with open(filename,"r") as G:
    lines = G.readlines()
  reciprocal_vectors = flex.vec3_double()
  #Example of Silicon F d 3-bar m unit cell oriented along lab axes
  A_mat = matrix.sqr((5.43,0.,0.,0.,5.43,0.,0.,0.,5.43))
  A_star = A_mat.inverse()
  for line in lines[1:]:
    tokens = line.strip().split('\t')
    assert len(tokens)==5
    miller = matrix.col([float(i) for i in tokens[0:3]])
    reciprocal_vector = A_star*miller
    reciprocal_vectors.append(reciprocal_vector.elems)
  return reciprocal_vectors

def test_case_obs_data(verbose=True):
  R = parse_input(libtbx.env.under_dist("rstbx", "indexing/si_brief.dat"))
  return do_index(R,verbose)

def test_case_synthetic_data(verbose=True):
  R = parse_synthetic(libtbx.env.under_dist("rstbx", "indexing/si_synthetic.dat"))
  return do_index(R,verbose)

if __name__=='__main__':
  print("Test autoindexing on synthetic reciprocal space positions")
  test_case_synthetic_data()
  print("Test autoindexing on monoscan data")
  test_case_obs_data()
