#ifndef SMTBX_AB_INITIO_DENSITY_MODIFICATION_H
#define SMTBX_AB_INITIO_DENSITY_MODIFICATION_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/loops.h>

namespace smtbx { namespace ab_initio { namespace density_modification {

  namespace af = scitbx::af;

  template<class FloatType, class AccessorType>
  void flip_charges_in_place(af::ref<FloatType, AccessorType> const &rho,
                             FloatType delta)
  {
    for (std::size_t i=0; i<rho.size(); ++i) {
      if (rho[i] < delta) rho[i] = -rho[i];
    }
  }

  template<class FloatType, class AccessorType>
  void low_density_elimination_in_place_tanaka_et_al_2001(
    af::ref<FloatType, AccessorType> const &rho,
    FloatType rho_s)
  {
    for (std::size_t i=0; i<rho.size(); ++i) {
      if (rho[i] <= 0) {
        rho[i] = 0;
      }
      else {
        FloatType x = rho[i]/rho_s;
        rho[i] *= 1 - std::exp(-0.5*x*x);
      }
    }
  }

}}} // namespace smtbx::ab_initio:: density_modification

#endif
