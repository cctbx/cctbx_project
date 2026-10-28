#ifndef SCITBX_MATH_VON_MISES_H
#define SCITBX_MATH_VON_MISES_H

#include <scitbx/math/quadrature.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/shared.h>
#include <cmath>

namespace scitbx { namespace math {

  //! <|phi|> in radians for the von Mises distribution
  //! p(phi) ~ exp(x cos(phi)) on [-pi, pi], x >= 0: an n-point
  //! Gauss-Legendre rule over [0, min(pi, 15/sqrt(x))]; beyond that the
  //! weight is below exp(-45), since 1 - cos(phi) >= 2 phi^2/pi^2. Within
  //! 1e-10 of the integral for x from 0 to 1e12 with 64 points.
  class von_mises_mean_abs_angle
  {
    public:
      explicit
      von_mises_mean_abs_angle(int n_nodes = 64)
      {
        quadrature::gauss_legendre_engine<double> gle(n_nodes);
        t_ = gle.x();
        w_ = gle.w();
      }

      //! a shared 64-point instance
      static von_mises_mean_abs_angle const& get()
      {
        static const von_mises_mean_abs_angle instance;
        return instance;
      }

      double operator()(double x) const
      {
        double pi = constants::pi;
        double upper = pi;
        if (x*pi*pi > 225.) upper = 15./std::sqrt(x);
        if (!(upper > 0.)) return 0.;
        double num = 0., den = 0.;
        for(std::size_t k=0; k < t_.size(); k++) {
          double phi = 0.5 * upper * (t_[k] + 1.);
          double wk = w_[k] * std::exp(x * (std::cos(phi) - 1.));
          num += wk * phi;
          den += wk;
        }
        return num / den;
      }

    protected:
      af::shared<double> t_, w_;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_VON_MISES_H
