#ifndef CCTBX_ELTBX_XRAY_SCATTERING_FORM_FACTOR_H
#define CCTBX_ELTBX_XRAY_SCATTERING_FORM_FACTOR_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace eltbx {

//! X-ray scattering tables.
namespace xray_scattering {

  //! Helper.
  template<class Derived>
  class isotropic_form_factor_mixin
  {
    Derived const & heir() const {
      return static_cast<Derived const &>(*this);
    }

    public:
      /*! \brief Value at the point stol
          (sin-theta-over-lambda), given stol^2.
       */
      /*! See also: at_stol(), at_d_star(), at_d_star_sq(),
                    uctbx::unit_cell::stol()
       */
      double
      at_stol_sq(double stol_sq) const
      {
        return heir().at_x_sq(stol_sq);
      }

      /*! \brief Value at the point stol
          (sin-theta-over-lambda).
       */
      /*! See also: at_stol_sq(), at_d_star(), at_d_star_sq(),
       */
      double
      at_stol(double stol) const
      {
        return heir().at_x_sq(stol * stol);
      }

      /*! \brief Value at the point d_star
          (1/d), given d_star^2.
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star_sq(double d_star_sq) const
      {
        return heir().at_x_sq(d_star_sq / 4);
      }

      /*! \brief Value at the points d_star
          (1/d), given d_star^2.
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star(),
                    uctbx::unit_cell::d_star_sq()
       */
      af::shared<double>
      at_d_star_sq(af::const_ref<double> const& d_star_sq) const
      {
        af::shared<double>
          result(d_star_sq.size(), af::init_functor_null<double>());
        for(std::size_t i=0;i<d_star_sq.size();i++) {
          result[i] = at_d_star_sq(d_star_sq[i]);
        }
        return result;
      }

      /*! \brief Value at the point d_star
          (1/d).
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star_sq(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star(double d_star) const
      {
        return heir().at_x_sq(d_star * d_star / 4);
      }
  };

}}} // cctbx::eltbx::xray_scattering

#endif // GUARD
