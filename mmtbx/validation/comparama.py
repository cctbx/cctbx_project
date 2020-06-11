from __future__ import absolute_import, division, print_function

import iotbx.phil
from libtbx import group_args
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from mmtbx.validation.ramalyze import ramalyze, find_region_max_value
import math
import numpy as np
from collections import Counter
from scitbx.array_family import flex
import six
from six.moves import zip

master_phil_str = '''
comparama {
  nproc = 1
    .type = int
  show_labels = True
    .type = bool
    .help = show labels for outlier residues
  outlier_favored = lime
    .type = str
    .help = Color of outlier->favored arrows. None if not needed
  outlier_allowed = lime
    .type = str
    .help = Color of outlier->allowed arrows. None if not needed
  allowed_outlier = red
    .type = str
    .help = Color of allowed->outlier arrows. None if not needed
  allowed_favored = green
    .type = str
    .help = Color of allowed->favored arrows. None if not needed
  favored_outlier = red
    .type = str
    .help = Color of favored->outlier arrows. None if not needed
  favored_allowed = orange
    .type = str
    .help = Color of favored->allowed arrows. None if not needed
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
    self.plots = None
    self.params = params
    if self.params is None:
      self.params = rcompare.get_default_params().comparama
    self.rama1 = ramalyze(model1.get_hierarchy(), out=null_out())
    self.rama2 = ramalyze(model2.get_hierarchy(), out=null_out())
    self.results = []
    # looping technique trying to recover when 1 or several residues are
    # missing
    i1 = i2 = 0
    self.skipped_1 = []
    self.skipped_2 = []
    while i1 < len(self.rama1.results) and i2 < len(self.rama2.results):
      r1 = self.rama1.results[i1]
      r2 = self.rama2.results[i2]
      if r1.id_str() == r2.id_str():
        # regular calculations
        diff = math.sqrt((r1.phi-r2.phi)**2 + (r1.psi-r2.psi)**2)
        v = determine_validation_change_text(r1, r2)
        diff2 = math.sqrt(get_distance(r1.phi, r2.phi)**2 +
            get_distance(r1.psi, r2.psi)**2)
        diff3 = two_rama_points((r1.phi, r1.psi), (r2.phi, r2.psi)).min_length()
        assert approx_equal(diff2, diff3), "%s, %s" % ((r1.phi, r1.psi), (r2.phi, r2.psi))
        self.results.append((r1.id_str(), diff2, r1.phi, r1.psi, r2.phi, r2.psi, v, r2.res_type, r1.score/100, r2.score/100))
        i1 += 1
        i2 += 1
      else:
        skip_1 = False
        # figure out what to skip
        if r1.chain_id == r2.chain_id:
          if r1.resseq_as_int() < r2.resseq_as_int():
            skip_1 = True
        else:
          if r1.resseq_as_int() > r2.resseq_as_int():
            skip_1 = True
        if skip_1:
          i1 += 1
          self.skipped_1.append(r1)
        else:
          i2 += 1
          self.skipped_2.append(r2)
    self.res_columns = None
    if len(self.results) > 0:
      self.res_columns = list(zip(*self.get_results()))

  def get_results(self):
    return self.results

  def get_results_as_vec3(self):
    r1 = flex.vec3_double()
    r2 = flex.vec3_double()
    assert len(self.results) > 0
    for r in self.results:
      r1.append((r[2], r[3], 0))
      r2.append((r[4], r[5], 0))
    return r1, r2

  def get_ramalyze_objects(self):
    return self.rama1, self.rama2

  def get_skipped(self):
    return self.skipped_1, self.skipped_2

  def get_number_results(self):
    if len(self.results) > 0:
      v1, v2 = rama_rescale(self.results)
      rv1 = [1-x for x in v1]
      rv2 = [1-x for x in v2]
      return group_args(
          mean_diff=np.mean(self.res_columns[1]),
          std_diff=np.std(self.res_columns[1]),
          sum_1 = np.sum(self.res_columns[-2]),
          sum_2 = np.sum(self.res_columns[-1]),
          n_res = len(self.res_columns[-1]),
          scaled_sum_1 = np.sum(v1),
          scaled_sum_2 = np.sum(v2),
          rev_scaled_sum_1 = np.sum(rv1),
          rev_scaled_sum_2 = np.sum(rv2),
          counts = Counter(self.res_columns[-4]),
          )
    return None

  def get_plots(self, wrap_arrows=True):
    if self.plots is not None:
      return self.plots
    self.plots = self.rama2.get_plots(
        show_labels=self.params.show_labels,
        point_style='bo',
        markersize=3,
        markeredgecolor="black",
        dpi=300,
        markerfacecolor="white")
    for pos, plot in six.iteritems(self.plots):
      # prepare data
      arrows_info = []
      if self.params.allowed_outlier is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "Allowed -> OUTLIER")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.allowed_outlier))
      if self.params.allowed_favored is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "Allowed -> Favored")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.allowed_favored))
      if self.params.favored_outlier is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "Favored -> OUTLIER")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.favored_outlier))
      if self.params.favored_allowed is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "Favored -> Allowed")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.favored_allowed))
      if self.params.outlier_favored is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "OUTLIER -> Favored")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.outlier_favored))
      if self.params.outlier_allowed is not None:
        arrs = [x for x in self.results if (x[-3]==pos and x[-4] == "OUTLIER -> Allowed")]
        arrs.sort(key=lambda x:x[1], reverse=True)
        arrows_info.append((arrs, self.params.outlier_allowed))

      for data, color in arrows_info:
        # print (len(data))
        if data and len(data) < 0: continue
        ad = [((x[2], x[3]),(x[4], x[5])) for x in data]
        add_arrows_on_plot(
            p=plot,
            arrows_data=ad,
            wrap_arrows=wrap_arrows,
            color=color)
    return self.plots

  @staticmethod
  def get_default_params():
    """
    Get extracted params
    """
    return iotbx.phil.parse(
          input_string=master_phil_str,
          process_includes=True).extract()

def rama_rescale(results):
  res1 = []
  res2 = []
  for (r1_id_str, diff2, r1_phi, r1_psi, r2_phi, r2_psi, v, r2_res_type,
      r1_score, r2_score) in results:
    for phi, psi, score, res in [(r1_phi, r1_psi, r1_score, res1),
                                 (r2_phi, r2_psi, r2_score, res2)]:
      max_value = find_region_max_value(r2_res_type, phi, psi)
      if max_value is None:
        res.append(score)
      else:
        res.append(score/max_value[1])
  return res1, res2

def breake_arrow_if_needed(abeg, aend, plot_ranges):
  eps = 1e-3
  tp = two_rama_points(abeg, aend)
  actual_len = tp.length(abeg, aend)
  min_len = tp.min_length()
  best_xy_multipliers = tp.get_xy_multipliers()
  result = []
  if best_xy_multipliers == [0,0]:
    return [(abeg,aend)]
  # Now we figure out how to brake it.
  result = [ [abeg, (0,0)], [(0,0), aend] ]
  ix = 0 if best_xy_multipliers[0] == -1 else 1
  iy = 0 if best_xy_multipliers[1] == -1 else 1
  if approx_equal(abeg[0], aend[0], eps, out=None):
    # case where x1 == x2
    result[0][1] = (abeg[0], plot_ranges[0][iy])
    result[1][0] = (abeg[0], plot_ranges[0][1-iy])
  elif best_xy_multipliers.count(0) == 1:
    # general case, 1 border crossing
    # y = ax + b
    n_aend = (aend[0]+360*best_xy_multipliers[0], aend[1]+360*best_xy_multipliers[1])
    a = (n_aend[1]-abeg[1]) / (n_aend[0] - abeg[0])
    b = n_aend[1] - a*n_aend[0]
    if best_xy_multipliers[0] != 0:
      # x wrapping, calculating y
      y = a*(plot_ranges[0][ix]) + b
      y = get_distance(y, 0)
      result[0][1] = (plot_ranges[0][ix],   y)
      result[1][0] = (plot_ranges[0][1-ix], y)
    else:
      # y wrapping, calculating x
      x = (plot_ranges[1][iy] - b) / a
      x = get_distance(x, 0)
      result[0][1] = (x, plot_ranges[1][iy])
      result[1][0] = (x, plot_ranges[1][1-iy])
  else:
    # both sides cutting. just go to the corner to make things simple
    result[0][1] = (plot_ranges[0][ix], plot_ranges[1][iy])
    result[1][0] = (plot_ranges[0][1-ix], plot_ranges[1][1-iy])
  return result

def add_arrows_on_plot(
    p,
    arrows_data,
    color='green',
    wrap_arrows=True,
    plot_ranges=[(-180, 180), (-180, 180)]):
  """
  p - pyplot
  arrows_data - [((x,y beginning), (x,y end)), ... ((xy),(xy))]
  wrap_arrows - draw shortest possible arrow - wrap around plot edges
  ranges - ranges of the plot
  """
  import matplotlib.patches as patches
  import matplotlib.lines as lines

  style="Simple,head_length=10,head_width=5,tail_width=1"
  for arrow in arrows_data:
    r = [(arrow[0], arrow[1])]
    if wrap_arrows:
      r = breake_arrow_if_needed(arrow[0], arrow[1], plot_ranges)
      for l_coors in r[:-1]:
        l = lines.Line2D(
            xdata = [l_coors[0][0], l_coors[1][0]],
            ydata = [l_coors[0][1], l_coors[1][1]],
            linewidth=1.7, color=color, zorder=10)
        p.plot.add_line(l)
    p.plot.add_patch(patches.FancyArrowPatch(
        r[-1][0],
        r[-1][1],
        arrowstyle=style,
        color = color,
        linewidth=0.5,
        zorder=10,
        ))
