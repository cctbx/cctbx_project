#ifndef CCTBX_XRAY_SIGMAA_H
#define CCTBX_XRAY_SIGMAA_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/packing_order.h>
#include <scitbx/array_family/block_iterator.h>
#include <scitbx/sym_mat3.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace xray { namespace sigmaa {

  template <typename FloatType, typename ComplexType>
  af::shared<FloatType>
  compute(
    af::const_ref<FloatType> const& f_obs,
    af::const_ref<ComplexType> const& f_calc)
  {
    CCTBX_ASSERT(f_obs.size() == f_calc.size());
    af::shared<FloatType> result;
    //for(std::size_t i=0; i<f_obs.size(); i++) {
    //  FloatType f_obs_i = f_obs[i];
    //}
    return result;
  }

}}} // namespace cctbx::xray::sigmaa

#endif // CCTBX_XRAY_SIGMAA_H
