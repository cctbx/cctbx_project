#if !defined(CCTBX_ELTBX_XRAY_SCATTERING_H)
#error "Do not include this file directly."
#endif
      //! Analytical approximation to the scattering factor.
      /*! stol_sq = sin-theta-over-lambda-squared: (sin(theta)/lambda)^2
          <p>
          See also: uctbx::unit_cell::at_stol_sq()
       */
      double
      at_stol_sq(double stol_sq) const
      {
        double sf = c();
        for(std::size_t i=0;i<n_ab();i++)
          sf += a(i) * std::exp(-b(i) * stol_sq);
        return sf;
      }

      //! Analytical approximation to the scattering factor.
      /*! See also: at_stol_sq(), uctbx::unit_cell::stol()
       */
      double
      at_stol(double stol) const
      {
        return at_stol_sq(stol * stol);
      }

      //! Analytical approximation to the scattering factor.
      /*! See also: at_stol_sq(), uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star_sq(double d_star_sq) const
      {
        return at_stol_sq(d_star_sq / 4.);
      }

      //! Analytical approximation to the scattering factor.
      /*! See also: at_stol_sq(), uctbx::unit_cell::d_star_sq()
       */
      af::shared<double>
      at_d_star_sq(af::const_ref<double> const& d_star_sq) const
      {
        af::shared<double> result(
          d_star_sq.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<d_star_sq.size();i++) {
          result[i] = at_d_star_sq(d_star_sq[i]);
        }
        return result;
      }

      //! Analytical gradient of the scattering factor at the point d_star.
      double
      gradient_at_d_star(double d_star) const
      {
        double d_star_sq = d_star * d_star;
        double result = 0;
        for(std::size_t i=0;i<n_ab();i++) {
          result -= (a(i)*b(i)*d_star)/(2*std::exp((b(i)*d_star_sq/4)));
        }
        return result;
      }

      //! Analytical integral of the scattering factor from zero to d_star.
      double
      integral_at_d_star(double d_star) const
      {
        double result = c() * d_star;
        for(std::size_t i=0;i<n_ab();i++) {
          result += one_gaussian_term_integral_at_d_star(
            static_cast<double>(a(i)),
            static_cast<double>(b(i)),
            d_star);
        }
        return result;
      }

      //! Analytical gradients w.r.t. a,b,c at the point d_star.
      gaussian
      gradients_at_d_star_sq(double d_star_sq) const
      {
        double stol_sq = d_star_sq / 4.;
        af::small<float, gaussian::max_n_ab> gr_a;
        af::small<float, gaussian::max_n_ab> gr_b;
        for(std::size_t i=0;i<n_ab();i++) {
          gr_a.push_back(std::exp(-b(i) * stol_sq));
        }
        for(std::size_t i=0;i<n_ab();i++) {
          gr_b.push_back(-stol_sq * a(i) * gr_a[i]);
        }
        return gaussian(gr_a, gr_b, 1);
      }
