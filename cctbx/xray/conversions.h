#ifndef CCTBX_XRAY_CONVERSIONS_H
#define CCTBX_XRAY_CONVERSIONS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cmath>

namespace cctbx { namespace xray {

  /*! From the Xtal ADDREF documentation (http://xtal.crystal.uwa.edu.au/):

      The calculation of (F) is straightforward except where F is close
      to zero or negative. To avoid the asymptotic form of this
      conversion, the following expression is used in ADDREF.

      sigma_F = (RLP sigma_I) / (F + {F2 + RLP sigma_I}1/2)
   */
  template <typename FloatType=double>
  struct f_sq_as_f_xtal_3_7
  {
    f_sq_as_f_xtal_3_7(
      FloatType const& f_sq,
      FloatType const& sigma_f_sq,
      FloatType const& tolerance=1.e-6)
    {
      init_f(f_sq);
      if ((f < tolerance && sigma_f_sq < tolerance) || sigma_f_sq < 0) sigma_f=0;
      else sigma_f = sigma_f_sq / (f + std::sqrt(sigma_f_sq + f * f));
    }

    f_sq_as_f_xtal_3_7(
      FloatType const& f_sq)
    {
      init_f(f_sq);
    }

    void
    init_f(FloatType f_sq)
    {
      if (f_sq > 0) f = std::sqrt(f_sq);
      else f = 0;
    }

    FloatType f;
    FloatType sigma_f;
  };

  template<typename FloatType=double>
  struct f_sq_as_f_crystals
  {
    f_sq_as_f_crystals(
      FloatType f_sq,
      FloatType sigma_f_sq,
      FloatType tolerance=1e-6)
    {
      init_f(f_sq);
      // I assume sigma_f_sq > 0 here
      if (f_sq < sigma_f_sq) {
        sigma_f = sigma_f_sq;
      }
      else {
        sigma_f = sigma_f_sq/(2*f);
      }
    }

    f_sq_as_f_crystals(
      FloatType f_sq)
    {
      init_f(f_sq);
    }

    void
    init_f(FloatType f_sq)
    {
      if (f_sq > 0) {
        f = std::sqrt(f_sq);
      }
      else {
        f = -std::sqrt(-f_sq);
      }
    }

    FloatType f;
    FloatType sigma_f;
  };

  template <typename FloatType=double>
  struct f_as_f_sq
  {
    f_as_f_sq(FloatType const& f,
              FloatType const& sigma_f)
    {
      f_sq = f * f;
      if (f == 0) {
        sigma_f_sq = sigma_f * sigma_f;
      }
      else {
        sigma_f_sq = 2 * f * sigma_f;
      }
    }

    f_as_f_sq(FloatType const& f)
    {
      f_sq = f * f;
    }

    FloatType f_sq;
    FloatType sigma_f_sq;
  };

  template <template<typename> class FsqAsF, typename FloatType=double>
  struct array_f_sq_as_f
  {
    array_f_sq_as_f(
      af::const_ref<FloatType> const& f_sq,
      af::const_ref<FloatType> const& sigma_f_sq,
      FloatType const& tolerance=1.e-6)
    {
      CCTBX_ASSERT(sigma_f_sq.size() == f_sq.size());
      f.reserve(f_sq.size());
      sigma_f.reserve(f_sq.size());
      for(std::size_t i=0;i<f_sq.size();i++) {
        FsqAsF<FloatType> r(f_sq[i], sigma_f_sq[i], tolerance);
        f.push_back(r.f);
        sigma_f.push_back(r.sigma_f);
      }
    }

    array_f_sq_as_f(
      af::const_ref<FloatType> const& f_sq)
    {
      f.reserve(f_sq.size());
      for(std::size_t i=0;i<f_sq.size();i++) {
        FsqAsF<FloatType> r(f_sq[i]);
        f.push_back(r.f);
      }
    }

    af::shared<FloatType> f;
    af::shared<FloatType> sigma_f;
  };

  template <typename FloatType=double>
  struct array_f_as_f_sq
  {
    array_f_as_f_sq(af::const_ref<FloatType> const& f,
                    af::const_ref<FloatType> const& sigma_f)
    {
      CCTBX_ASSERT(sigma_f.size() == f.size());
      f_sq.reserve(f.size());
      sigma_f_sq.reserve(f.size());
      for(std::size_t i=0;i<f.size();i++) {
        f_as_f_sq<FloatType> r(f[i], sigma_f[i]);
        f_sq.push_back(r.f_sq);
        sigma_f_sq.push_back(r.sigma_f_sq);
      }
    }

    array_f_as_f_sq(af::const_ref<FloatType> const& f)
    {
      f_sq.reserve(f.size());
      for(std::size_t i=0;i<f.size();i++) {
        f_as_f_sq<FloatType> r(f[i]);
        f_sq.push_back(r.f_sq);
      }
    }

    af::shared<FloatType> f_sq;
    af::shared<FloatType> sigma_f_sq;
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_CONVERSIONS_H
