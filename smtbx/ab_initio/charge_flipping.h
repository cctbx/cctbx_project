#ifndef SMTBX_AB_INITIO_CHARGE_FLIPPING_H
#define SMTBX_AB_INITIO_CHARGE_FLIPPING_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/loops.h>

namespace smtbx { namespace ab_initio {

  namespace af = scitbx::af;

  template<class FloatType, class AccessorType>
  void flip_charges_in_place(af::ref<FloatType, AccessorType> rho,
                             FloatType delta)
  {
    typedef typename AccessorType::index_type index_type;
    for (af::nested_loop<index_type> i(rho.accessor().focus());
         !i.over(); i.incr())
    {
      if (rho(i()) < delta) rho(i()) = -rho(i());
    }
  }

}}

#endif
