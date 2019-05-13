from __future__ import division
from __future__ import print_function

from mmtbx.validation.ramalyze import ramalyze, find_region_max_value, get_favored_peaks, \
    RAMALYZE_FAVORED, res_types
from libtbx.utils import null_out
from mmtbx.rotamer import ramachandran_eval
import numpy as np
from libtbx import easy_mp


master_phil_str = """\
evalurama {
  use_allowed = False
    .type = bool
    .help = include residues that are allowed but still falls into grid squares
  corners_inside = 1
    .type = int
    .help = for grid, how many corners of the square need to be inside favored
  grid_size = 10
    .type = float
    .help = for grid, what size of the squares
  grid_center_on_peak = False
    .type = bool
    .help = center grid on peak. Otherwise, center will be (0,0)
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
    inside = False
    for i in [0,1]:
      for j in [0,1]:
        if (self.low_left[0] <= point[0] + 360*i <= self.high_right[0] and
            self.low_left[1] <= point[1] + 360*j <= self.high_right[1]):
          inside = True
    if inside:
      self.n_points += 1
    return inside

  def _show(self):
    print(self.low_left, self.high_right, self.probability, self.n_points, self.prob_points)


class grid(list):
  def __init__(self, *args):
    super(grid, self).__init__(*args)

  def scale_to_1(self):
    prob_sum = 0
    for x in self:
      prob_sum += x.probability
    for x in self:
      x.probability /= prob_sum
    prob_sum = 0
    for x in self:
      prob_sum += x.probability

  def add_points(self, points, stop_on_outsiders=True):
    """ points: [(xy),(xy)...] """
    self.total_points = len(points)
    print("    Using", self.total_points, "points")

    self.not_placed_points = 0
    for p in points:
      placed = False
      for x in self:
        placed = x.add_point_if_inside(p)
        if placed:
          break
      if not placed:
        # basically all points in the vicinity of this peak, including outliers
        # print "      Rejected point:", p
        self.not_placed_points += 1
      if not placed and not stop_on_outsiders:
        assert 0, "point was not placed in any square (%f, %f)" % (p[0], p[1])
    self.total_points -= self.not_placed_points
    print("    Total residues placed:", self.total_points)
    if self.total_points > 0:
      for x in self:
        x.prob_points = x.n_points/self.total_points

  def get_corr(self):
    a = []
    b = []
    for x in self:
      a.append(x.probability)
      b.append(x.prob_points)
    if len(a) == 0 or len(b) == 0:
      return None
    if self.total_points < 30:
      return None
    corr = np.corrcoef(a,b)
    print("    Correlation", corr[0,1])
    return corr

  def get_total_points(self):
    t = 0
    for x in self:
      t += x.n_points
    return t

  def _show(self):
    for x in self:
      x._show()

class grid_over_favored(object):
  def __init__(self, rama_type, rama_region, start_point, grid_size, corners_inside, use_allowed=False):
    self.rama_type = rama_type
    self.rama_region = rama_region
    self.start_point = start_point
    self.grid_size = grid_size
    self.corners_inside = corners_inside
    self.r = ramachandran_eval.RamachandranEval()
    self.grid = grid() # list of square objects
    points_x = grid_over_favored._get_grid_points(start_point[0], grid_size)
    points_y = grid_over_favored._get_grid_points(start_point[1], grid_size)
    for x in points_x:
      for y in points_y:
        n_inside = 0
        for dx in [0,1]:
          for dy in [0,1]:
            reg = find_region_max_value(rama_type, x+grid_size*dx, y+grid_size*dy)
            if reg is not None and reg[0] == rama_region:
              n_inside += 1
        if n_inside >= corners_inside:
          v = self.r.evaluate(rama_type, [x+0.5*grid_size, y+0.5*grid_size])
          self.grid.append(square((x,y), (x+grid_size, y+grid_size), v))
    print("    Number of grid cells", len(self.grid))

    self.grid.scale_to_1()

  def add_points(self, points, stop_on_outsiders=True):
    self.grid.add_points(points, stop_on_outsiders=stop_on_outsiders)

  def get_corr(self):
    c = self.grid.get_corr()
    if c is not None:
      return c[0,1]
    return None

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

def ramalyze_parallel(hierarchy):
  return ramalyze(hierarchy, out=null_out())

class eval(object):
  def __init__(self, models, params, log):
    self.models = models
    self.log = log
    self.params = params
    self.total_rama = []
    rama_res = easy_mp.pool_map(
        processes=self.params.nproc,
        fixed_func=ramalyze_parallel,
        args=[m.get_hierarchy() for m in self.models])
    for r in rama_res:
      for res in r.results:
        self.total_rama.append(res)

    self.results = {}
    for rama_type in range(5):
      print("Working with", res_types[rama_type], "plot")
      for peak in get_favored_peaks(rama_type):
        print("  Working with peak:", peak)
        start_point = peak[0] if self.params.grid_center_on_peak else (0,0)
        gof = grid_over_favored(
            rama_type=rama_type,
            rama_region=peak[0],
            start_point=start_point,
            grid_size=self.params.grid_size,
            corners_inside=self.params.corners_inside)
        points = []
        for r in self.total_rama:
          if (r.res_type == rama_type
              and (r.rama_type == RAMALYZE_FAVORED or self.params.use_allowed)
              and find_region_max_value(rama_type, r.phi, r.psi,
                  allow_outside=self.params.use_allowed) == peak):
            points.append((r.phi, r.psi))
        gof.add_points(points, stop_on_outsiders=self.params.use_allowed)
        c = gof.get_corr(), gof.grid.get_total_points()
        if c.count(None) > 0:
          c = None
        self.results[(rama_type, peak)] = c

  def get_results(self):
    return self.results
