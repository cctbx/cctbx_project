#ifndef CCTBX_XRAY_TARGETS_LEAST_SQUARES_H
#define CCTBX_XRAY_TARGETS_LEAST_SQUARES_H

#include <cctbx/xray/targets/common_results.h>
#include <cmath>

namespace cctbx { namespace xray { namespace targets {

  struct least_squares : common_results
  {
    protected:
      bool compute_scale_using_all_data_;
      char obs_type_;
      double scale_factor_;
      public:

    //! derivatives_depth == -1: compute only scale factor
    least_squares(
      bool compute_scale_using_all_data,
      char obs_type,
      af::const_ref<double> const& obs,
      af::const_ref<double> const& weights,
      af::const_ref<bool> const& r_free_flags,
      af::const_ref<std::complex<double> > const& f_calc,
      int derivatives_depth,
      double scale_factor)
    :
      common_results(obs.size()),
      compute_scale_using_all_data_(compute_scale_using_all_data),
      obs_type_(obs_type)
    {
      TBXX_ASSERT(obs_type == 'F' || obs_type == 'I');
      TBXX_ASSERT(weights.size() == 0 || weights.size() == obs.size());
      TBXX_ASSERT(r_free_flags.size() == 0
               || r_free_flags.size() == obs.size());
      TBXX_ASSERT(f_calc.size() == obs.size());
      TBXX_ASSERT(derivatives_depth >= -1 && derivatives_depth <= 2);
      TBXX_ASSERT(scale_factor >= 0);
      TBXX_ASSERT(!(derivatives_depth == -1 && scale_factor != 0));
      const double* wghts = weights.begin();
      const bool* rff = r_free_flags.begin();
      if (!rff) compute_scale_using_all_data = true;
      std::size_t n_work = 0;
      double sum_w_o_sq_work = 0;
      double sum_w_o_sq_test = 0;
      {
        double num = 0;
        double denom = 0;
        for(std::size_t i=0;i<obs.size();i++) {
          double o = obs[i];
          double w = (wghts ? wghts[i] : 1);
          double a = f_calc[i].real();
          double b = f_calc[i].imag();
          double aabb = a*a + b*b;
          double c = (obs_type == 'F' ? std::sqrt(aabb) : aabb);
          double w_o = w * o;
          double w_o_sq = w_o * o;
          if (compute_scale_using_all_data || !rff[i]) {
            num += w_o * c;
            denom += w * c * c;
          }
          if (!rff || !rff[i]) {
            n_work++;
            sum_w_o_sq_work += w_o_sq;
          }
          else {
            sum_w_o_sq_test += w_o_sq;
          }
        }
        if (scale_factor == 0) {
          TBXX_ASSERT(denom > 0);
          scale_factor_ = num / denom;
          if (derivatives_depth == -1) return;
        }
        else {
          scale_factor_ = scale_factor;
        }
      }
      TBXX_ASSERT(sum_w_o_sq_work > 0);
      if (derivatives_depth != 0) {
        if (compute_scale_using_all_data && n_work != obs.size()) {
          throw std::runtime_error(
            "Sorry: cctbx::xray::targets::least_squares:"
            " derivatives for compute_scale_using_all_data"
            " not implemented.");
        }
        gradients_work_.reserve(n_work);
        if (derivatives_depth == 2) {
          hessian_work_.reserve(n_work);
        }
      }
      double target_test = 0;
      double grad_factor = (obs_type == 'F' ? -2 : -4)
                         * scale_factor_ / sum_w_o_sq_work;
      for(std::size_t i=0;i<obs.size();i++) {
        double o = obs[i];
        double w = (wghts ? wghts[i] : 1);
        double a = f_calc[i].real();
        double b = f_calc[i].imag();
        double aabb = a*a + b*b;
        double c = (obs_type == 'F' ? std::sqrt(aabb) : aabb);
        double delta = o - scale_factor_ * c;
        double wd =  w * delta;
        double t = wd * delta;
        target_per_reflection_[i] = t;
        if (rff && rff[i]) {
          target_test += t;
        }
        else {
          target_work_ += t;
          if (derivatives_depth != 0) {
            double c_cub = c * c * c;
            if (c == 0 || c_cub == 0) {
              gradients_work_.push_back(std::complex<double>(0,0));
              if (derivatives_depth == 2) {
                hessian_work_.push_back(scitbx::vec3<double>(1,1,1));
              }
            }
            else {
              double gf = grad_factor * wd;
              if (obs_type == 'F') gf /= c;
              gradients_work_.push_back(gf * f_calc[i]);
              if (derivatives_depth == 2) {
                scitbx::vec3<double> cw;
                if (obs_type == 'F') {
                  double term = -grad_factor * w;
                  double oocc = o / c_cub;
                  /*daa*/ cw[0] = term * (scale_factor_ - b * b * oocc);
                  /*dbb*/ cw[1] = term * (scale_factor_ - a * a * oocc);
                  /*dab*/ cw[2] = term * a * b * oocc;
                }
                else {
                  double term = -2 * grad_factor * scale_factor_ * w;
                  /*daa*/ cw[0] = term * a * a + gf;
                  /*dbb*/ cw[1] = term * b * b + gf;
                  /*dab*/ cw[2] = term * a * b;
                }
                hessian_work_.push_back(cw);
              }
            }
          }
        }
      }
      target_work_ /= sum_w_o_sq_work;
      if (rff && sum_w_o_sq_test > 0) {
        target_test_ = boost::optional<double>(
          target_test / sum_w_o_sq_test);
      }
    }

    bool
    compute_scale_using_all_data() const
    {
      return compute_scale_using_all_data_;
    }

    char
    obs_type() const { return obs_type_; }

    double
    scale_factor() const { return scale_factor_; }

    /* Mathematica input for gradients, hessian, analytical scale_factor (k) :

       obs_type == 'F':
         t=w(fo-k Sqrt[a^2+b^2])^2
         D[t,a]
         D[t,b]
         D[D[t,a],a]
         D[D[t,b],b]
         D[D[t,a],b]

         t=w1(fo1-k Sqrt[a1^2+b1^2])^2 + w2(fo2-k Sqrt[a2^2+b2^2])^2
         Solve[D[t,k]==0,k]

       obs_type == 'I':
         t=w(io-k (a^2+b^2))^2
         D[t,a]
         D[t,b]
         D[D[t,a],a]
         D[D[t,b],b]
         D[D[t,a],b]

         t=w1(io1-k(a1^2+b1^2))^2 + w2(io2-k(a2^2+b2^2))^2
         Solve[D[t,k]==0,k]
     */
  };

}}} // namespace cctbx::xray::targets

#endif // GUARD
