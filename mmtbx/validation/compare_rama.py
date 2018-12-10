from __future__ import division

import iotbx.phil
from libtbx.utils import null_out
from mmtbx.validation.ramalyze import ramalyze, RAMALYZE_OUTLIER, \
    RAMALYZE_ALLOWED, RAMALYZE_FAVORED, RAMALYZE_ANY, RAMALYZE_NOT_FAVORED
import math

master_phil_str = '''
compare_rama {
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
      self.params = rcompare.get_default_params().compare_rama
    rama1 = ramalyze(model1.get_hierarchy(), out=null_out())
    rama2 = ramalyze(model2.get_hierarchy(), out=null_out())
    self.results = []
    for r1, r2 in zip(rama1.results, rama2.results):
      assert r1.id_str() == r2.id_str()
      diff = math.sqrt((r1.phi-r2.phi)**2 + (r1.psi-r2.psi)**2)
      v = determine_validation_change_text(r1, r2)
      diff2 = math.sqrt(get_distance(r1.phi, r2.phi)**2 +
          get_distance(r1.psi, r2.psi)**2)
      self.results.append((r1.id_str(), diff2, r1.phi, r1.psi, r2.phi, r2.psi, v))

  def get_results(self):
    return self.results

  @staticmethod
  def get_default_params():
    """
    Get extracted params
    """
    return iotbx.phil.parse(
          input_string=master_phil_str,
          process_includes=True).extract()
