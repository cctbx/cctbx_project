import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")

from cctbx.array_family import flex

import math


class _mainstream_shelx_weighting(boost.python.injector,
                                  ext.mainstream_shelx_weighting):

  def __str__(self):
    """ A string representation of the weighting scheme in a format that is
        appropriate for the CIF item _refine_ls_weighting_details.
    """
    if round(self.a, 4) in (0.1, 0.2):
      a = "%.1f" %self.a
    else:
      a = "%.4f" %self.a
    if round(self.b, 4) == 0: b_part=""
    else: b_part = "+%.4fP" %self.b
    return ("w=1/[\s^2^(Fo^2^)+(%sP)^2^%s]"
            " where P=(Fo^2^+2Fc^2^)/3" %(a, b_part))

  def type(self):
    return "calc"

  def optimise_parameters(self, fo_sq, f_calc,
                          scale_factor, n_independent_params):
    """ Find optimal values of a and b that give a flat analysis of the variance
        when binned by Fc/max(Fc), and a goodness of fit close to 1.

        This is done in a grid search fashion similar to Shelxl.

        self is not modified in place; instead a new instance of the weighting
        scheme is returned.

        It is intended that f_calc should already contain the contribution from
        f_mask (if a solvent mask is used).
    """
    assert f_calc.is_complex_array()
    weighting = ext.mainstream_shelx_weighting(a=self.a, b=self.b)

    def compute_chi_sq(fo_sq, fc_sq, a,b):
      weighting.a = a
      weighting.b = b
      weights = weighting(
        fo_sq.data(), fo_sq.sigmas(), fc_sq.data(), scale_factor)
      return (flex.sum(
        weights * flex.pow2(fo_sq.data() - scale_factor * fc_sq.data())))

    fo_sq = fo_sq.deep_copy()
    fo_sq.data().set_selected(fo_sq.data() < 0, 0)
    fc_sq = f_calc.as_intensity_array()

    fo2 = fo_sq.data().deep_copy()
    fo2 /= scale_factor
    sigmas = fo_sq.sigmas() / scale_factor
    sigmas_sq = flex.pow2(sigmas)
    fc2 = fc_sq.data()

    # determine starting values for a and b, formulae taken from shelxl code
    p = (fo2 + 2 * fc2)/3
    p_sq = flex.pow2(p)
    x = flex.sum((flex.pow2(fo2-fc2)-sigmas) * (p_sq/sigmas_sq))
    y = flex.sum( flex.pow2(p_sq)/sigmas_sq)
    z = flex.sum(p)
    start_a = math.sqrt(max(0.0001, 0.64*x/max(1e-8, y)))
    start_b = 0.5 * z * start_a**2 /fo_sq.size()
    a_step = 0.2 * start_a
    b_step = 0.4 * start_b

    # sort data and setup binning by fc/fc_max
    fc_sq_over_fc_sq_max = fc_sq.data()/flex.max(fc_sq.data())
    permutation = flex.sort_permutation(fc_sq_over_fc_sq_max)
    fc_sq_over_fc_sq_max = fc_sq.customized_copy(
      data=fc_sq_over_fc_sq_max).select(permutation)
    fc_sq = fc_sq.select(permutation)
    fo_sq = fo_sq.select(permutation)
    n_bins = 10
    bins = []
    bin_max = 0
    bin_limits = flex.size_t(1, 0)
    bin_count = flex.size_t()
    for i in range(n_bins):
      bin_limits.append(int(math.ceil((i+1) * fc_sq.size()/n_bins)))
      bin_count.append(bin_limits[i+1] - bin_limits[i])

    n = fo_sq.size()/(fo_sq.size()-n_independent_params)

    # search on a 9x9 grid to determine best values of a and b
    gridding = flex.grid(9,9)
    while (a_step > 1e-4 and b_step > 5e-3):
      tmp = flex.double(gridding, 0)
      binned_chi_sq = [tmp.deep_copy() for i in range(n_bins)]
      start_a = max(start_a, 4*a_step) - 4*a_step
      start_b = max(start_b, 4*b_step) - 4*b_step
      for i_bin in range(n_bins):
        sel = flex.size_t_range(bin_limits[i_bin], bin_limits[i_bin+1])
        fc2 = fc_sq.select(sel)
        fo2 = fo_sq.select(sel)
        b = start_b
        for j in range(9):
          a = start_a
          b += b_step
          for k in range(9):
            a += a_step
            binned_chi_sq[i_bin][j,k] += compute_chi_sq(fo2, fc2, a, b)
      min_variance = 9e9
      j_min, k_min = (0, 0)
      for j in range(9):
        for k in range(9):
          variance = 0
          for i_bin in range(n_bins):
            if bin_count[i_bin] == 0: continue
            goof = math.sqrt(binned_chi_sq[i_bin][j,k]*n/bin_count[i_bin])
            variance += (goof-1)**2
          min_variance = min(variance, min_variance)
          if variance == min_variance:
            j_min = j
            k_min = k
      start_a += k_min*a_step
      start_b += j_min*b_step
      if k_min == 8:
        a_step *= 2
        continue
      elif k_min != 0:
        a_step /= 4
      if j_min == 8:
        b_step *= 2
        continue
      elif j_min != 0:
        b_step /=4
      if start_a <= 1e-4: a_step /= 4
      if start_b <= 1e-3: b_step /= 4
    if start_a > 0.2:
      start_a = 0.2
      start_b = 0
    weighting.a = start_a
    weighting.b = start_b
    return weighting

class _unit_weighting(boost.python.injector,
                      ext.unit_weighting):

  def __str__(self):
    return "w=1"

  def type(self):
    return "unit"

  def optimise_parameters(self, fo_sq, fc_sq,
                          scale_factor, n_independent_params):
    # no parameters to optimise!
    return self
