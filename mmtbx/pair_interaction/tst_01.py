from __future__ import absolute_import, division, print_function
from mmtbx.pair_interaction import pair_interaction
from libtbx.test_utils import approx_equal
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_pair_interaction_ext")

def test_density_props_cpp():
  # This matches numbers from NCIPLOT (FORTRAN code).
  e="C_"
  wfc_obj = pair_interaction.load_wfc(e)
  a_xyz = [0,0,0]
  p=[-0.11027379670470783, -0.11027379670470783, -1.0551368983523539]
  r = ext.atom_density_props(p = p, a_xyz = a_xyz, wfc_obj = wfc_obj)
  #print ("density", r.density)
  #print ("gradient_vector", r.gradient_vector)
  #print ("hessian", r.hessian)
  #print ("has_silva_interaction (dori)", r.has_silva_interaction("dori"))
  #print ("has_silva_interaction (sedd)", r.has_silva_interaction("sedd"))
  #silva = r.cal_silva()
  #print ("cal_silva", silva)
  #print ("get_dori_value", r.get_dori_value())
  #print ("get_sedd_value", r.get_sedd_value())
  return r

if __name__=="__main__":
  r=test_density_props_cpp()
  h = [-0.21713634909425739, 3.4647776644023028E-003, 3.3152161869310914E-002,
       3.4647776644023028E-003, -0.21713634909425739, 3.3152161869310914E-002,
       3.3152161869310914E-002, 3.3152161869310914E-002, 9.6609945080372378E-002]
  assert approx_equal(r.density, 0.15909578867162319)
  assert approx_equal(r.gradient_vector,
    [ 2.4326523805013919E-002, 2.4326523805013919E-002,0.23276438866116667])
  assert approx_equal(r.hessian, h)
  silva = r.cal_silva()
  assert approx_equal(silva, 8.37330299844e-05)
  assert approx_equal(r.get_dori_value(), 0.66372941504308625)
  assert approx_equal(r.get_sedd_value(), 6.7056325871680365)
  assert not r.has_silva_interaction("dori")
  assert not r.has_silva_interaction("sedd")
