#ifndef CCTBX_XRAY_ASU_MAPPINGS_H
#define CCTBX_XRAY_ASU_MAPPINGS_H

#include <cctbx/xray/scatterer.h>

namespace cctbx { namespace xray {

  template <typename AsuMappingsType, typename ScattererType>
  void
  asu_mappings_process(
    AsuMappingsType& asu_mappings,
    af::const_ref<ScattererType> const& scatterers)
  {
    asu_mappings.reserve(asu_mappings.mappings().size() + scatterers.size());
    for(std::size_t i=0;i<scatterers.size();i++) {
      asu_mappings.process(scatterers[i].site);
    }
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_ASU_MAPPINGS_H
