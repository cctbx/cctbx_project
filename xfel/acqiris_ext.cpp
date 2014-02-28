/* -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*- */

/*-
 * $Id$
 */

#include <boost/python.hpp>

#include <scitbx/array_family/shared.h>

//#define DEBUG 1


/*
 * Integrate Acqiris trace from @p i_sig to @p i_bkg.  This
 * implementation trades accuracy for speed, by only adding corrected
 * numbers which leads to cache misses.  Note that because the signal
 * is negative, signal is subtracted from background.  Should be
 * faster than apd_hitfind.  For reading out the Silicon Avalanche
 * Photodiode SAR 3000 E1 (http://XXX).
 *
 * @param w     Waveform
 * @param i_sig Index of first signal point
 * @param i_bkg Index of first background point
 * @param nmemb Number of points in integral
 * @return      Integral
 */
double
acqiris_integrate(
  const scitbx::af::const_ref<double>& w,
  std::size_t i_sig,
  std::size_t i_bkg,
  std::size_t nmemb)
{
  double integral = w[i_bkg] - w[i_sig];
  for (std::size_t i = 1; i < nmemb; i++)
    integral += w[i_bkg + i] - w[i_sig + i];
  return (integral);
}


/*
 * Optimal window width depends on the decay time of the diode.  At
 * LCLS the pulse length (tenths of fs) is much shorter than the
 * interval between digitised samples (hundreds of ps).  Should be
 * longer than w.size() to make any sense.  This is really more suited
 * for the Optodiode SXUV100 (http://www.optodiode.com/pro_09.html)

 * @param w      Waveform
 * @param wwidth Window width
 * @return       Hit
 */
#ifdef DEBUG
scitbx::af::shared<double>
#else
bool
#endif
apd_hitfind(const scitbx::af::const_ref<double>& w, std::size_t wwidth)
{
#ifdef DEBUG
  scitbx::af::shared<double> window3(
    std::max(std::size_t(0), w.size() - wwidth));
#endif
  double rsum[3] = { 0, 0, 0 }; // running sums for expectations (E)
  double skew_amax, skew_mean, skew_ssq, skew_stddev, skew_sum;

  /*
   * Compute non-central moments up to order three for all windows.
   * This implementation will accumulate error, but requires only one
   * flat O(n) pass over the data.
   *
   * Compute skewness, the third central moment, in each window and
   * store in window3.  Establish baseline (mean and standard
   * deviation) over the first samples.  This is really a threshold.
   * Maybe it makes more sense to use the median and make use of
   * trailing samples as well?
   *
   * Remember absolute maximum skewness above the mean outside the
   * baseline region.
   */
  skew_amax = skew_ssq = skew_sum = 0;
  for (std::size_t i = 0; i < std::min(wwidth, w.size()); i++) {
    rsum[0] += w[i];
    rsum[1] += w[i] * w[i];
    rsum[2] += w[i] * w[i] * w[i];
  }
  for (std::size_t i = wwidth; i < w.size(); i++) {
    const std::size_t j  = i - wwidth;
    const double      c1 = rsum[0] / wwidth;
    const double      c2 = std::max(0.0, rsum[1] / wwidth - c1 * c1);
    const double      c3 = c2 > 0
      ? (rsum[2] / wwidth - 3 * c1 * c2 - c1 * c1 * c1) / std::pow(c2, 1.5)
      : 0;

    if (j < wwidth) {
      skew_sum += c3;
      skew_ssq += c3 * c3;
    } else {
      skew_amax = std::max(skew_amax, std::abs(c3 - skew_sum / wwidth));
    }

#ifdef DEBUG
    window3[j] = c3;
#endif

    rsum[0] += w[i]               - w[j];
    rsum[1] += w[i] * w[i]        - w[j] * w[j];
    rsum[2] += w[i] * w[i] * w[i] - w[j] * w[j] * w[j];
  }

  skew_mean   = skew_sum / wwidth;
  skew_stddev = std::max(0.0, skew_ssq - skew_sum * skew_mean);
  skew_stddev = std::sqrt(
    wwidth > 1 ? skew_stddev / (wwidth - 1) : skew_stddev);

  /*
   * For actual APD, it may make a lot more sense to look at one-step
   * transitions.  XXX Read specs for the two diodes.
   */
  if (skew_amax > 15 * skew_stddev) { // XXX magic
#ifdef DEBUG
    printf("It's a HIT (amax %f, 15stddev %f)\n", skew_amax, 15 * skew_stddev);
#else
    return (true);
#endif
  }

#ifdef DEBUG
  return window3;
#else
  return (false);
#endif
}


BOOST_PYTHON_MODULE(acqiris_ext)
{
  // XXX named arguments!
  using namespace boost::python;
  def("acqiris_integrate", acqiris_integrate);
  def("apd_hitfind", apd_hitfind);
}
