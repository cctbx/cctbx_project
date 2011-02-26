#ifndef CCTBX_XRAY_TARGETS_CORRELATION_H
#define CCTBX_XRAY_TARGETS_CORRELATION_H

#include <cctbx/xray/targets/common_results.h>
#include <boost/scoped_array.hpp>
#include <cmath>

namespace cctbx { namespace xray { namespace targets {

  struct correlation : common_results
  {
    protected:
      char obs_type_;
      boost::optional<double> cc_;
      public:

    correlation() {}

    correlation(
      char obs_type,
      af::const_ref<double> const& obs,
      af::const_ref<double> const& weights,
      af::const_ref<bool> const& r_free_flags,
      af::const_ref<std::complex<double> > const& f_calc,
      int derivatives_depth)
    :
      common_results(obs.size()),
      obs_type_(obs_type)
    {
      TBXX_ASSERT(obs_type == 'F' || obs_type == 'I');
      TBXX_ASSERT(weights.size() == 0 || weights.size() == obs.size());
      TBXX_ASSERT(r_free_flags.size() == 0
               || r_free_flags.size() == obs.size());
      TBXX_ASSERT(f_calc.size() == obs.size());
      TBXX_ASSERT(derivatives_depth >= 0 && derivatives_depth <= 2);
      const double* wghts = weights.begin();
      const bool* rff = r_free_flags.begin();
      std::size_t n_work = 0;
      double sum_w = 0;
      double sum_wx = 0;
      double sum_wy = 0;
      boost::scoped_array<double> wb(new double[obs.size()]);
      boost::scoped_array<double> xb(new double[obs.size()]);
      boost::scoped_array<double> yb(new double[obs.size()]);
      double w = 1;
      for(std::size_t i=0;i<obs.size();i++) {
        if (rff == 0 || !rff[i]) {
          if (wghts != 0) {
            w = wghts[i];
            sum_w += wghts[i];
          }
          double x = obs[i];
          double y = (obs_type == 'F'
            ? std::abs(f_calc[i])
            : std::norm(f_calc[i]));
          sum_wx += w * x;
          sum_wy += w * y;
          wb[n_work] = w;
          xb[n_work] = x;
          yb[n_work] = y;
          n_work++;
        }
      }
      if (wghts == 0) sum_w = static_cast<double>(n_work);
      TBXX_ASSERT(sum_w != 0);
      double wxm = sum_wx / sum_w;
      double wym = sum_wy / sum_w;
      double sum_wxy = 0;
      double sum_wxx = 0;
      double sum_wyy = 0;
      for(std::size_t i_work=0;i_work<n_work;i_work++) {
        double w = wb[i_work];
        double& xc = xb[i_work];
        double& yc = yb[i_work];
        xc -= wxm;
        yc -= wym;
        sum_wxy += w * xc * yc;
        sum_wxx += w * xc * xc;
        sum_wyy += w * yc * yc;
      }
      double cc_den_sq = sum_wxx * sum_wyy;
      if (cc_den_sq != 0) {
        double cc_den = std::sqrt(cc_den_sq);
        TBXX_ASSERT(cc_den != 0);
        double cc = sum_wxy / cc_den;
        if (derivatives_depth != 0) {
          gradients_work_.reserve(n_work);
          TBXX_ASSERT(derivatives_depth != 2); // not implemented
          if (derivatives_depth == 2) {
            hessians_work_.reserve(n_work);
          }
          std::size_t i_work = 0;
          for(std::size_t i=0;i<obs.size();i++) {
            if (rff == 0 || !rff[i]) {
              double w = wb[i_work];
              double xc = xb[i_work];
              double yc = yb[i_work];
              i_work++;
              double d_sum_wxy = w * xc;
              double d_sum_wyy = 2 * w * yc;
              double d_cc_den_sq = sum_wxx * d_sum_wyy;
              double d_cc_den = 1 / (2 * cc_den) * d_cc_den_sq;
              double d_cc = (d_sum_wxy * cc_den - sum_wxy * d_cc_den)
                          / cc_den_sq;
              double gf;
              if (obs_type == 'F') {
                gf = std::abs(f_calc[i]);
                if (gf != 0) gf = d_cc / gf;
              }
              else {
                gf = 2 * d_cc;
              }
              gradients_work_.push_back(-gf * f_calc[i]);
            }
          }
        }
        target_work_ = 1 - cc;
        cc_ = cc;
      }
    }

    char
    obs_type() const { return obs_type_; }

    boost::optional<double>
    cc() const { return cc_; }
  };

}}} // namespace cctbx::xray::targets

#endif // GUARD
