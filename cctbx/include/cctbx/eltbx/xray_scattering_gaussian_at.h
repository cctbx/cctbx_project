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
        for (std::size_t i = 0; i < n_ab(); i++)
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
