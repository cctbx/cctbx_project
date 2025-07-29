from __future__ import division
from math import pi, sin, cos, atan2, sqrt, atan

from mmtbx.conformation_dependent_library.cdl_utils import round_to_ten
from mmtbx.validation.mean_devs_PRO_phi_psi import mean_devs as PRO_phi_psi
from mmtbx.validation.mean_devs_others_phi_psi import mean_devs as others_phi_psi
from mmtbx.validation.mean_devs_VAL_THR_ILE_phi_psi import mean_devs as VAL_THR_ILE_phi_psi

class radial_deviation:
  def __init__(self, r, t, input_in_radians=False, input_in_cartesian=False):
    if input_in_cartesian:
      self.r = sqrt(r*r + t*t)
      self.t = atan(t/r)
      input_in_radians=True
    self.r = r
    self.t = t
    if not input_in_radians: self.t *= pi/180

  def __repr__(self):
    #return u'%0.2f \u2220 %6.1f' % (self.r, self.t/pi*180)
    return "'%4.1f/_%6.1f'" % (self.r*100, self.t/pi*180)

  def __add__(self, other):
    r = sqrt(self.r**2 + other.r**2 + 2*self.r*other.r*cos(other.t-self.t))
    numer = other.r * sin(other.t-self.t)
    denom = self.r + other.r * cos(other.t-self.t)
    t = self.t + atan2(numer, denom)
    new = radial_deviation(0,0)
    new.r = r
    new.t = t
    return new

  def __sub__(self, other):
    new = radial_deviation(0,0)
    new.r = other.r
    new.t = other.t + pi
    return self.__add__(new)

  def __truediv__(self, other):
    if type(other)==type(1):
      self.r/=other
      return self
    else:
      assert 0

  def __cmp__(self, other):
    if self.r<other.r: return -1
    return 1

  def __eq__(self, other):
    return self.r == other.r

  def __ne__(self, other):
    return self.r != other.r

  def __lt__(self, other):
    return self.r < other.r

  def __le__(self, other):
    return self.r <= other.r

  def __gt__ (self, other):
    return self.r > other.r

  def __ge__(self, other):
    return self.r >= other.r

def get_phi_psi_correction(result,
                           residue,
                           phi_psi,
                           display_phi_psi_correction=False,
                           verbose=False):
  rc = None
  key = (round_to_ten(phi_psi[0]), round_to_ten(phi_psi[1]))
  if residue.resname=='PRO':
    correction = PRO_phi_psi.get(key, None)
  elif residue.resname in ['VAL', 'THR', 'ILE']:
    correction = VAL_THR_ILE_phi_psi.get(key, None)
  else:
    correction = others_phi_psi.get(key, None)
  if correction:
    correction = radial_deviation(*tuple(correction), input_in_radians=True)
    current = radial_deviation(result.deviation, result.dihedral)
    rc = current-correction
    start = (current.r>=0.25)
    finish = (rc.r>=0.25)
    if verbose:
      print('current %s is corrected with %s' % (current, correction)),
      print(' to give %s\n' % rc)
    if display_phi_psi_correction and (start or finish):
      show_phi_psi_correction(residue, phi_psi, current, correction, rc)
    return rc.r, rc.t/pi*180, start, finish
  else:
    return None

def show_phi_psi_correction(residue, phi_psi, current, correction, rc, units=100):
  import matplotlib.pyplot as plt

  xs = [current.t, correction.t, rc.t]
  ys = [current.r, correction.r, rc.r]
  labels = ['default', 'correction', 'result']

  for i, (x, y) in enumerate(zip(xs, ys)):
      y*=units
      plt.polar(x, y, 'ro')
      plt.text(x, y,
               '%s' % (labels[i]),
               )
  key = (round_to_ten(phi_psi[0]), round_to_ten(phi_psi[1]))
  plt.title('%s (%5.1f, %5.1f) (%4.0f, %4.0f)\n%s minus %s equals %s' %(
    residue.id_str(),
    phi_psi[0],
    phi_psi[1],
    key[0],
    key[1],
    current,
    correction,
    rc,
   )
  )
  plt.show()

if __name__ == '__main__':
  # tests
  v1 = radial_deviation(1,0)
  v2 = radial_deviation(-1,0)
  print(v1,v2)
  print(v1-v2)
  v1 = radial_deviation(1,1)
  print(v1)
  v1 = radial_deviation(1,1, input_in_radians=True)
  print(v1)
  v1 = radial_deviation(1,1, input_in_cartesian=True)
  print(v1)
