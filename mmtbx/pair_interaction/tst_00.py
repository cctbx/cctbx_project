from __future__ import absolute_import, division, print_function
import iotbx.pdb
from mmtbx.pair_interaction import pair_interaction
from libtbx.test_utils import approx_equal

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_pair_interaction_ext")

pdb_str = """
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
END
"""

def test_density_props_cpp():
  # All numbers verified against Java and Python implementation: all numbers
  # match.
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  ph = pdb_inp.construct_hierarchy()
  atoms=ph.atoms()
  a=atoms[0]
  e=a.element.strip(" ").lower()
  if(len(e)==1):e=e+"_"
  wfc_obj = pair_interaction.load_wfc(e)
  a_xyz = [10.448638218462468  , 3.71765835974282866E-002 , 1.0816225869801266]
  p=[10.149391801428232, 1.4386837495289542, 2.4831297529116525]
  r = ext.atom_density_props(p = p, a_xyz = a_xyz, wfc_obj = wfc_obj)
  #print ("density", r.density)
  #print ("gradient_vector", r.gradient_vector)
  #print ("hessian", r.hessian)
  return r

if __name__=="__main__":
  r=test_density_props_cpp()
  h = [-0.02470403, -0.01459848, -0.01459848, -0.01459848,  0.04055025,
        0.06837132, -0.01459848,  0.06837132,  0.04055025]
  assert approx_equal(r.density, 0.027461463669643152)
  assert approx_equal(r.gradient_vector, [0.00832535, -0.03899142, -0.03899142])
  assert approx_equal(r.hessian, h)
