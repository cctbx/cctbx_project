from cctbx.eltbx.caasf import wk1995
from cctbx.array_family import flex
import math

mean_number_of_atoms_per_amino_acid = {'C': 5, 'N': 3, 'O': 1}

class empty: pass

class wilson_plot:

  def __init__(self, f_obs_array, asu_contents):
    self.info = f_obs_array.info()
    # compute <fobs^2> in resolution shells
    self.mean_fobs_sq = f_obs_array.mean_sq(
      use_binning=0001,
      use_multiplicities=0001)
    # compute <s^2> = <(sin(theta)/lambda)^2> in resolution shells
    stol_sq = f_obs_array.sin_theta_over_lambda_sq()
    stol_sq.use_binner_of(f_obs_array)
    self.mean_stol_sq = stol_sq.mean(
      use_binning=0001,
      use_multiplicities=0001)
    # cache scattering factor info
    caasf = {}
    for chemical_type in asu_contents.keys():
      caasf[chemical_type] = wk1995(chemical_type)
    # compute expected f_calc^2 in resolution shells
    self.expected_f_sq = flex.double()
    for stol_sq in self.mean_stol_sq:
      sum_fj_sq = 0
      for chemical_type, n_atoms in asu_contents.items():
        f0 = caasf[chemical_type].at_stol_sq(stol_sq)
        sum_fj_sq += f0 * f0 * n_atoms
      self.expected_f_sq.append(sum_fj_sq)
    self.expected_f_sq *= f_obs_array.space_group().order_z() \
                        * f_obs_array.space_group().n_ltr()
    # fit to straight line
    self.x = self.mean_stol_sq
    self.y = flex.log(self.mean_fobs_sq / self.expected_f_sq)
    fit = flex.linear_regression(self.x, self.y)
    assert fit.is_well_defined()
    self.fit_y_intercept = fit.y_intercept()
    self.fit_slope = fit.slope()
    self.fit_correlation = fit.correlation()
    self.wilson_k = math.exp(self.fit_y_intercept)
    self.wilson_b = -self.fit_slope / 2

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
