#ifndef CCTBX_XRAY_SCATTERER_UTILS_H
#define CCTBX_XRAY_SCATTERER_UTILS_H

#include <cctbx/xray/scatterer.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>

namespace cctbx { namespace xray {

  template <typename FloatType,
            typename CaasfType,
            typename LabelType>
  af::shared<std::size_t>
  apply_symmetry(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group,
    af::ref<scatterer<FloatType,CaasfType,LabelType> > const& scatterers,
    double min_distance_sym_equiv=0.5,
    double u_star_tolerance=0,
    bool assert_is_positive_definite=false,
    bool assert_min_distance_sym_equiv=true)
  {
    af::shared<std::size_t> special_position_indices;
    for(std::size_t i=0;i<scatterers.size();i++) {
      sgtbx::site_symmetry site_symmetry = scatterers[i].apply_symmetry(
        unit_cell,
        space_group,
        min_distance_sym_equiv,
        u_star_tolerance,
        assert_is_positive_definite,
        assert_min_distance_sym_equiv);
      if (!site_symmetry.is_point_group_1()) {
        special_position_indices.push_back(i);
      }
    }
    return special_position_indices;
  }

  template <typename FloatType,
            typename CaasfType,
            typename LabelType>
  af::shared<scatterer<FloatType,CaasfType,LabelType> >
  rotate(
    uctbx::unit_cell const& unit_cell,
    scitbx::mat3<double> const& rotation_matrix,
    af::const_ref<scatterer<FloatType,CaasfType,LabelType> > const& scatterers)
  {
    af::shared<scatterer<FloatType,CaasfType,LabelType> >
      rot_scatterers(af::reserve(scatterers.size()));
    for(std::size_t i=0;i<scatterers.size();i++) {
      CCTBX_ASSERT(!scatterers[i].anisotropic_flag);
      cartesian<> c = unit_cell.orthogonalize(scatterers[i].site);
      cartesian<> rc = rotation_matrix * c;
      rot_scatterers.push_back(scatterers[i]);
      rot_scatterers[i].site = unit_cell.fractionalize(rc);
    }
    return rot_scatterers;
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_UTILS_H
