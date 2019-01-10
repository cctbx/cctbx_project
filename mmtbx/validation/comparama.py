from __future__ import division

import iotbx.phil
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from mmtbx.validation.ramalyze import ramalyze
import math

master_phil_str = '''
comparama {
  nproc = 1
    .type = int
}
'''

def get_to_0_360(angle):
  if angle < 0:
    return angle + 360
  return angle

def get_distance(a1, a2):
  a = a1-a2
  if a < -180:
    a = a + 360
  if a > 180:
    a = a-360
  return a

class two_rama_points(object):
  def __init__(self, a, b):
    self.a = a
    self.b = b
    self.min_len = None
    self.best_xy_multipliers = None

  def length(self, abeg, aend):
    return math.sqrt((abeg[0]-aend[0])**2 + (abeg[1]-aend[1])**2)
  def get_xy_multipliers(self):
    if self.best_xy_multipliers is not None:
      self.min_length()
    return self.best_xy_multipliers

  def min_length(self, plot_ranges=((-180, 180),(-180, 180))):
    if self.min_len is not None:
      return self.min_len
    base_len = self.length(self.a, self.b)
    self.min_len = base_len
    # print("base_len", base_len)
    xspan = plot_ranges[0][1] - plot_ranges[0][0]
    yspan = plot_ranges[1][1] - plot_ranges[1][0]
    self.best_xy_multipliers = [0,0]
    for x in [-1,0,1]:
      for y in [-1,0,1]:
        new_x = self.b[0] + x*xspan
        new_y = self.b[1] + y*yspan
        tlen = self.length(self.a, (new_x, new_y))
        if tlen < self.min_len:
          self.min_len = tlen
          self.best_xy_multipliers = [x,y]
          # print("  min_len", min_len, best_xy_multipliers)
    return self.min_len

def determine_validation_change_text(r1, r2):
  rt1 = r1.ramalyze_type()
  rt2 = r2.ramalyze_type()
  # if r1.is_outlier() and not r2.is_outlier():
  #   return "fixed"
  # if not r1.is_outlier() and r2.is_outlier():
  #   return "broke"
  if rt1 == rt2:
    return rt1
  return "%s -> %s" % (rt1, rt2)


class rcompare(object):
  def __init__(self, model1, model2, params=None, log=null_out()):
    self.params = params
    if self.params is None:
      self.params = rcompare.get_default_params().comparama
    self.rama1 = ramalyze(model1.get_hierarchy(), out=null_out())
    self.rama2 = ramalyze(model2.get_hierarchy(), out=null_out())
    self.results = []
    print dir(self.rama1.results[0])
    for r1, r2 in zip(self.rama1.results, self.rama2.results):
      assert r1.id_str() == r2.id_str()
      diff = math.sqrt((r1.phi-r2.phi)**2 + (r1.psi-r2.psi)**2)
      v = determine_validation_change_text(r1, r2)
      diff2 = math.sqrt(get_distance(r1.phi, r2.phi)**2 +
          get_distance(r1.psi, r2.psi)**2)
      diff3 = two_rama_points((r1.phi, r1.psi), (r2.phi, r2.psi)).min_length()
      assert approx_equal(diff2, diff3), "%s, %s" % ((r1.phi, r1.psi), (r2.phi, r2.psi))
      self.results.append((r1.id_str(), diff2, r1.phi, r1.psi, r2.phi, r2.psi, v, r2.res_type, r1.score/100, r2.score/100))
      # print "Score:", r1.score, r2.score

  def get_results(self):
    return self.results

  def get_ramalyze_objects(self):
    return self.rama1, self.rama2

  @staticmethod
  def get_default_params():
    """
    Get extracted params
    """
    return iotbx.phil.parse(
          input_string=master_phil_str,
          process_includes=True).extract()
