import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from cctbx.array_family import flex
from libtbx.utils import plural_s
import math

mean_number_of_atoms_per_amino_acid = {'C': 5, 'N': 3, 'O': 1}

class empty: pass

class wilson_plot(object):

  def __init__(self, f_obs, asu_contents):
    assert f_obs.is_real_array()
    self.info = f_obs.info()
    f_obs_selected = f_obs.select(f_obs.data() > 0)
    f_obs_selected.use_binning_of(f_obs)
    # compute <fobs^2> in resolution shells
    self.mean_fobs_sq = f_obs_selected.mean_sq(
      use_binning=True,
      use_multiplicities=True).data[1:-1]
    n_none = self.mean_fobs_sq.count(None)
    if (n_none > 0):
      error_message = "wilson_plot error: %d empty bin%s:" % plural_s(n_none)
      if (self.info is not None):
        error_message += "\n  Info: " + str(self.info)
      error_message += "\n  Number of bins: %d" % len(self.mean_fobs_sq)
      error_message += "\n  Number of f_obs > 0: %d" % (
        f_obs_selected.indices().size())
      error_message += "\n  Number of f_obs <= 0: %d" % (
        f_obs.indices().size() - f_obs_selected.indices().size())
      raise RuntimeError(error_message)
    self.mean_fobs_sq = flex.double(self.mean_fobs_sq)
    # compute <s^2> = <(sin(theta)/lambda)^2> in resolution shells
    stol_sq = f_obs_selected.sin_theta_over_lambda_sq()
    stol_sq.use_binner_of(f_obs_selected)
    self.mean_stol_sq = flex.double(stol_sq.mean(
      use_binning=True,
      use_multiplicities=True).data[1:-1])
    # cache scattering factor info
    gaussians = {}
    for chemical_type in asu_contents.keys():
      gaussians[chemical_type] = eltbx.xray_scattering.wk1995(
        chemical_type).fetch()
    # compute expected f_calc^2 in resolution shells
    self.expected_f_sq = flex.double()
    for stol_sq in self.mean_stol_sq:
      sum_fj_sq = 0
      for chemical_type, n_atoms in asu_contents.items():
        f0 = gaussians[chemical_type].at_stol_sq(stol_sq)
        sum_fj_sq += f0 * f0 * n_atoms
      self.expected_f_sq.append(sum_fj_sq)
    self.expected_f_sq *= f_obs_selected.space_group().order_z() \
                        * f_obs_selected.space_group().n_ltr()
    # fit to straight line
    self.x = self.mean_stol_sq
    self.y = flex.log(self.mean_fobs_sq / self.expected_f_sq)
    fit = flex.linear_regression(self.x, self.y)
    assert fit.is_well_defined()
    self.fit_y_intercept = fit.y_intercept()
    self.fit_slope = fit.slope()
    self.wilson_k = math.sqrt(math.exp(self.fit_y_intercept))
    self.wilson_b = -self.fit_slope / 2
    self.fit_correlation = flex.linear_correlation(self.x,self.y).coefficient()

  def xy_plot_info(self):
    r = empty()
    r.title = "Wilson Plot"
    if (self.info != 0):
      r.title += ": " + str(self.info)
    r.x = self.x
    r.y = self.y
    r.xLegend = "(sin(theta)/lambda)^2"
    r.yLegend = "ln(<Fobs^2>/<Fcalc^2>)"
    r.fit_y_intercept = self.fit_y_intercept
    r.fit_slope = self.fit_slope
    r.fit_correlation = self.fit_correlation
    r.overlayLegend = ("k=%f, b=%f, corr=%f" % (
      self.wilson_k, self.wilson_b, self.fit_correlation))
    return r
