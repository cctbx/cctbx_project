from __future__ import division
from scitbx.array_family import flex
import numpy as np

import boost.python
ext = boost.python.import_ext("mmtbx_pair_interaction_ext")
from mmtbx_pair_interaction_ext import *



def run():
#  h = ext.hessian(
#    distanceVector    = [1,2,3],
#    distanceReciprocal = 1,
#    distanceUnitVector = [4,5,6],
#    fac1 = 2,
#    fac2 = 5)
#  print list(h)
#  print h
#  print dir(ext)

  #o = ext.wfc()
  #print o

  print list(ext.points_and_pairs(a=1,b=1,c=1))
  STOP()

#  o = ext.wfc(
#    node_offsets = flex.vec3_int(10),
#    coefficients_of_first_derivative = flex.vec3_int(10),
#    coefficients_of_second_derivative = flex.vec3_int(10),
#    prefactor_of_first_derivative = 2,
#    prefactor_of_second_derivative = 4,
#    core_cutdens = 9)
#  o.a = 1.9
#  print o.a
#  print
#  print list(o.node_offsets)
#  #
#  print

  h = [0,0,0, 0,0,0, 0,0,0]
  np.array(h).reshape(3,3)
  print h

  #print flex.double(flex.grid(2))
  #density_props_obj = ext.density_props()
  #print dir(density_props_obj)

  #density_props_obj = ext.density_props(
  #  density = 0,
  #  gradient_vector = [0,0,0],
  #  hessian = h)

  density_props_obj = ext.density_props()

  h = [1,2,3, 4,5,6, 7,8,9]
  h = np.array(h).reshape(3,3)

  h = h.reshape(-1)


  density_props_obj_2 = ext.density_props(
    density = 1,
    gradient_vector = [2,3,4],
    hessian = h)

  density_props_obj.add(density_props = density_props_obj_2)
  print "density_props_obj.density        ", density_props_obj.density
  print "density_props_obj.gradient_vector", density_props_obj.gradient_vector
  print "density_props_obj.hessian        ", density_props_obj.hessian
  ###
  density_props_obj.gradient=100
  print "has_silva_interaction:",density_props_obj.has_silva_interaction()



if (__name__ == "__main__"):
  run()
