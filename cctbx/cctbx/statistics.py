import math
from cctbx_boost.arraytbx import shared
from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995

mean_number_of_atoms_per_amino_acid = {'C': 5, 'N': 3, 'O': 1}

class empty: pass

class wilson_plot:

  def __init__(self, fobs_set, asu_contents):
    self.info = fobs_set.info
    # compute <fobs^2> in resolution shells
    self.mean_fobs_sq = fobs_set.mean_sq(use_binning=1, use_multiplicities=1)
    # compute <s^2> = <(sin(theta)/lambda)^2> in resolution shells
    stol_sq_set = fobs_set.sin_theta_over_lambda_sq()
    stol_sq_set.binner = fobs_set.binner
    self.mean_stol_sq = stol_sq_set.mean(use_binning=1, use_multiplicities=1)
    # cache scattering factor info
    caasf = {}
    for chemical_type in asu_contents.keys():
      caasf[chemical_type] = CAASF_WK1995(chemical_type)
    # compute expected fcalc^2 in resolution shells
    self.expected_f_sq = shared.double()
    for stol_sq in self.mean_stol_sq:
      sum_fj_sq = 0
      for chemical_type, n_atoms in asu_contents.items():
        f0 = caasf[chemical_type].stol2(stol_sq)
        sum_fj_sq += f0 * f0 * n_atoms
      self.expected_f_sq.append(sum_fj_sq)
    self.expected_f_sq *= fobs_set.SgOps.OrderZ() * fobs_set.SgOps.nLTr()
    # fit to straight line
    self.x = self.mean_stol_sq
    self.y = shared.log(self.mean_fobs_sq / self.expected_f_sq)
    self.fit = shared.linear_regression(self.x, self.y)
    self.wilson_k = math.exp(self.fit.b())
    self.wilson_b = -self.fit.m() / 2

  def get_xy_plot_info(self):
    r = empty()
    r.title = "Wilson Plot"
    if (self.info != 0):
      r.title += ": " + str(self.info)
    r.x = self.x
    r.y = self.y
    r.xLegend = "(sin(theta)/lambda)^2"
    r.yLegend = "ln(<Fobs^2>/<Fcalc^2>)"
    r.b = self.fit.b()
    r.m = self.fit.m()
    r.cc = self.fit.cc()
    r.overlayLegend = ("k=%f, b=%f, corr=%f" % (
      self.wilson_k, self.wilson_b, self.fit.cc()))
    return r

  def dump_raw_xy(self):
    assert self.x.size() == self.y.size()
    for i,x in self.x.items():
      print x, self.y[i]
