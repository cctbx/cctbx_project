import math
"""assume data collection at two-theta angle of zero.
   assume a flat square detector normal to beam.
   calculate the fraction of a resolution shell covered by the detector.
"""
from scitbx import matrix
from spotfinder.core_toolbox import geometry_2d_base

def point(a,b):
  return matrix.col((a,b,0))

def magnitude(a):
  return math.sqrt( a.dot(a) )

def unit_vector(a):
  return a.normalize()

def dot(a,b):
  return a.dot(b)

def polar(a):
  return magnitude(a), math.atan2(a[1],a[0]) # r,theta.  z direction ignored

def cross_product(a,b):
  return a.cross(b)

def radius_to_resol(radius,parameter_dictionary):
  # assumes radius & distance in same units (e.g., mm)
  # wavelength and resolution in same units (e.g., Angstroms)
  distance = float(parameter_dictionary['distance'])
  wavelength = float(parameter_dictionary['wavelength'])
  theta = math.atan2(radius,distance)/2.0
  return wavelength/(2.0*math.sin(theta))

class Geom2d(geometry_2d_base):
  def __init__(self,pd):
    self.pd = pd
    geometry_2d_base.__init__(self,pixel_size=float(pd["pixel_size"]),
    size1=float(pd["size1"]),
    size2=float(pd["size2"]),
    xbeam=float(pd["xbeam"]),
    ybeam=float(pd["ybeam"]),
    distance=float(pd["distance"]),
    wavelength=float(pd["wavelength"]),
  )

if __name__=='__main__':
  from libtbx.test_utils import approx_equal
  pd = {'pixel_size':'1.0','xbeam':'3','ybeam':'6','size1':'10','size2':'10',
        'distance':'200','wavelength':'1'}
  f = Geom2d(pd)
  test_values = [(200,1.0), (100,1.0), (50,0.7643789), (30,0.2827813), (20,0.0)]
  for item in test_values:
    assert approx_equal(f(item[0]),item[1])
  print "OK"
