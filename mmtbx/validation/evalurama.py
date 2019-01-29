from __future__ import division

from mmtbx.validation.ramalyze import ramalyze, find_region_max_value, get_favored_peaks, \
    RAMALYZE_FAVORED, res_types
from libtbx.utils import null_out
from mmtbx.rotamer import ramachandran_eval
import numpy as np


master_phil_str = """\
evalurama {
  nproc = 1
    .type = int
}
"""

class square(object):
  def __init__(self, low_left, high_right, probability=0):
    self.low_left = low_left
    self.high_right = high_right
    self.probability = probability
    self.n_points = 0
    self.prob_points = 0

  def add_point_if_inside(self, point):
    if (self.low_left[0] <= point[0] <= self.high_right[0] and
        self.low_left[1] <= point[1] <= self.high_right[1]):
      self.n_points += 1
      return True
    return False

  def _show(self):
    print self.low_left, self.high_right, self.probability, self.n_points, self.prob_points


class grid(list):
  def __init__(self, *args):
    super(grid, self).__init__(*args)

  def scale_to_1(self):
    prob_sum = 0
    for x in self:
      prob_sum += x.probability
    for x in self:
      x.probability /= prob_sum
    print "prob_sum", prob_sum
    prob_sum = 0
    for x in self:
      prob_sum += x.probability
    print "after scaling:", prob_sum

  def add_points(self, points, normalize=True):
    """ points: [(xy),(xy)...] """
    self.total_points = len(points)
    for p in points:
      placed = False
      for x in self:
        placed = x.add_point_if_inside(p)
        if placed:
          break
      # print "what's up", p
      if not placed:
        assert 0, "point was not placed in any square (%f, %f)" % (p[0], p[1])
    for x in self:
      x.prob_points = x.n_points/self.total_points

  def get_corr(self):
    a = []
    b = []
    for x in self:
      a.append(x.probability)
      b.append(x.prob_points)
    print "a:", a
    print "b:", b
    corr = np.corrcoef(a,b)
    print "Correlation", corr
    return corr

  def _show(self):
    for x in self:
      x._show()

class grid_over_favored(object):
  def __init__(self, rama_type, rama_region, start_point, grid_size, corners_inside):
    self.rama_type = rama_type
    self.rama_region = rama_region
    self.start_point = start_point
    self.grid_size = grid_size
    self.corners_inside = corners_inside
    self.r = ramachandran_eval.RamachandranEval()
    self.grid = grid() # list of square objects
    points_x = grid_over_favored._get_grid_points(start_point[0], grid_size)
    points_y = grid_over_favored._get_grid_points(start_point[1], grid_size)
    # print points_x
    # print points_y
    for x in points_x:
      for y in points_y:
        n_inside = 0
        for dx in [0,1]:
          for dy in [0,1]:
            reg = find_region_max_value(0, x+grid_size*dx, y+grid_size*dy)
            if reg is not None and reg[0] == rama_region:
              n_inside += 1
        if n_inside >= corners_inside:
          v = self.r.evaluate(rama_type, [x+0.5*grid_size, y+0.5*grid_size])
          self.grid.append(square((x,y), (x+grid_size, y+grid_size), v))
    print "len grid", len(self.grid)

    self.grid.scale_to_1()
    # print "=" *80
    # self.grid.add_points([(-123,0), (-100,0), (-85,0)])
    # self.grid._show()
    # self.grid.get_corr()

  def add_points(self, points):
    self.grid.add_points(points)
    self.grid._show()

  def get_corr(self):
    c = self.grid.get_corr()
    return c[0,1]

  @staticmethod
  def _get_grid_points(start_point, size):
    """ returns values in [-180, 180] one of them is start point and """
    x = start_point
    result = [x]
    x -= size
    while x >= -180:
      result.append(x)
      x -= size
    x = start_point
    x += size
    while x <= 180:
      result.append(x)
      x += size
    return sorted(result)

class eval(object):
  def __init__(self, model, params, log):
    self.model = model
    self.log = log
    self.rama = ramalyze(self.model.get_hierarchy(), out=null_out())

    self.results = {}
    for rama_type in range(5):
      print "Working with", res_types[rama_type]
      for peak in get_favored_peaks(rama_type):
        print "  Working with peak:", peak
        gof = grid_over_favored(
            rama_type=rama_type,
            rama_region=peak[0],
            start_point=(0,0),
            grid_size=10,
            corners_inside=1)
        points = []
        for r in self.rama.results:
          if (r.res_type == rama_type and r.rama_type == RAMALYZE_FAVORED and
              find_region_max_value(rama_type, r.phi, r.psi) == peak):
            points.append((r.phi, r.psi))
        print points
        # STOP()
        c = None
        if len(points) > 30:
          gof.add_points(points)
          c = gof.get_corr()
        self.results[(rama_type, peak)] = c

  def get_results(self):
    return self.results
