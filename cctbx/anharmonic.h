#ifndef CCTBX_GRAM_CHARLIER_H
#define CCTBX_GRAM_CHARLIER_H

#include <vector>
#include <scitbx/matrix/tensors.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <scitbx/constants.h>

namespace cctbx { namespace adptbx { namespace anharmonic {
  using namespace scitbx::matrix::tensors;

  template <typename FloatType=double>
  struct GramCharlier4 {
    tensor_rank_3<FloatType> C;
    tensor_rank_4<FloatType> D;

    GramCharlier4() {}

    GramCharlier4(const af::shared<FloatType> &Cijk,
      const af::shared<FloatType> &Dijkl)
      : C(Cijk), D(Dijkl)
    {}

    std::complex<FloatType> calculate(const miller::index<> &h) const {
      using namespace scitbx::constants;
      FloatType c = C.sum_up(h), d = D.sum_up(h);
      return std::complex<FloatType>(
        1 + d * pi_sq*pi_sq * 2 / 3,
        -c * pi_sq*pi * 4 / 3);
    }

    af::shared<std::complex<FloatType> >
      gradient_coefficients(const miller::index<> &h) const
    {
      using namespace scitbx::constants;
      af::shared<std::complex<FloatType> > r(25);
      af::shared<FloatType> c = C.gradient_coefficients(h);
      const FloatType c_factor = -pi_sq * pi * 4 / 3,
        d_factor = pi_sq * pi_sq * 2 / 3;
      for (size_t i = 0; i < 10; i++) {
        r[i] = std::complex<FloatType>(0, c[i]*c_factor);
      }
      c = D.gradient_coefficients(h);
      for (size_t i = 0; i < 15; i++) {
        r[i + 10] = c[i] * d_factor;
      }
      return r;
    }

    af::shared<FloatType> data() const {
      af::shared<FloatType> rv(25);
      for (size_t i = 0; i < 10; i++) {
        rv[i] = C[i];
      }
      for (size_t i = 0; i < 15; i++) {
        rv[i + 10] = D[i];
      }
      return rv;
    }

  };

}}} // end of cctbx::adptbx::anharmonic
#endif
