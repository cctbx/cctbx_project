from __future__ import absolute_import, division, print_function
from six.moves import range
import math
import libtbx.load_env
from libtbx.test_utils import approx_equal
from scitbx import matrix
from rstbx.array_family import flex
from rstbx.indexing.do_index import do_index

def get_token(name,buf):
  token_positions = {'dspacing':[16,26],'qx':[102,112],'qy':[112,122],'qz':[122,132],}
  return float(buf[token_positions[name][0]:token_positions[name][1]])

def parse_input(filename):
  with open(filename,"r") as G:
    lines = G.readlines()
  reciprocal_vectors = flex.vec3_double()
  qvec = ('qx','qy','qz')
  for line in lines[1:len(lines)]:
    l_buffer = line.rstrip()
    assert len(l_buffer)==132
    extended_q_nm = matrix.col(
     [get_token(qvec[i],l_buffer) for i in range(3)]
    )
    checklength = math.sqrt(extended_q_nm.dot(extended_q_nm))
    assert approx_equal(1./checklength, get_token('dspacing',l_buffer), eps=0.0001)
    #convert from nanometers to Angstroms
    extended_q_Angstrom = 0.1 * extended_q_nm
    reciprocal_vectors.append(extended_q_Angstrom.elems)
  return reciprocal_vectors

def test_automatic_monoscan(verbose=True):
  R = parse_input(libtbx.env.under_dist('rstbx', 'indexing/results2.dat'))
  return do_index(R,verbose)

if __name__=='__main__':
  test_automatic_monoscan()
