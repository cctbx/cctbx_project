/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (R.W. Grosse-Kunstleve)
 */

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

  template <typename FloatType,
            typename CaasfType,
            typename LabelType,
            typename FloatTypePacked,
            typename FlexGridIndexType>
  std::size_t
  pack_parameters(
    const uctbx::unit_cell* unit_cell,
    af::const_ref<scatterer<FloatType,CaasfType,LabelType> > const& scatterers,
    af::versa<FloatTypePacked, af::flex_grid<FlexGridIndexType> >& x,
    bool site_flag,
    bool u_iso_flag,
    bool occupancy_flag)
  {
    CCTBX_ASSERT(x.accessor().is_trivial_1d());
    af::shared_plain<FloatTypePacked> xb = x.as_base_array();
    CCTBX_ASSERT(xb.size() == x.size());
    if (site_flag) {
      if (!unit_cell) {
        for(std::size_t i=0;i<scatterers.size();i++) {
          fractional<FloatType> const&
            site = scatterers[i].site;
          xb.insert(xb.end(), site.begin(), site.end());
        }
      }
      else {
        for(std::size_t i=0;i<scatterers.size();i++) {
          cartesian<FloatType>
            site = unit_cell->orthogonalize(scatterers[i].site);
          xb.insert(xb.end(), site.begin(), site.end());
        }
      }
    }
    if (u_iso_flag) {
      for(std::size_t i=0;i<scatterers.size();i++) {
        if (!scatterers[i].anisotropic_flag) {
          xb.push_back(scatterers[i].u_iso);
        }
      }
    }
    if (occupancy_flag) {
      for(std::size_t i=0;i<scatterers.size();i++) {
        xb.push_back(scatterers[i].occupancy);
      }
    }
    x.resize(af::flex_grid<>(xb.size()));
    return x.size();
  }

  template <typename FloatTypePacked,
            typename FloatType,
            typename CaasfType,
            typename LabelType>
  std::size_t
  unpack_parameters(
    const uctbx::unit_cell* unit_cell,
    std::size_t space_group_order_z,
    af::const_ref<FloatTypePacked> const& x,
    std::size_t start,
    af::ref<scatterer<FloatType,CaasfType,LabelType> > const& scatterers,
    bool site_flag,
    bool u_iso_flag,
    bool occupancy_flag)
  {
    if (site_flag) {
      CCTBX_ASSERT(x.size() >= start + scatterers.size() * 3);
      if (!unit_cell) {
        for(std::size_t i=0;i<scatterers.size();i++) {
          scatterers[i].site
            = fractional<double>(&x[start]);
          start += 3;
        }
      }
      else {
        for(std::size_t i=0;i<scatterers.size();i++) {
          scatterers[i].site
            = unit_cell->fractionalize(cartesian<double>(&x[start]));
          start += 3;
        }
      }
    }
    if (u_iso_flag) {
      CCTBX_ASSERT(x.size() >= start + scatterers.size());
      for(std::size_t i=0;i<scatterers.size();i++) {
        if (!scatterers[i].anisotropic_flag) {
          scatterers[i].u_iso = x[start];
          start++;
        }
      }
    }
    if (occupancy_flag) {
      CCTBX_ASSERT(space_group_order_z > 0);
      CCTBX_ASSERT(x.size() >= start + scatterers.size());
      for(std::size_t i=0;i<scatterers.size();i++) {
        scatterers[i].occupancy = x[start];
        scatterers[i].update_weight(space_group_order_z);
        start++;
      }
    }
    return start;
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_UTILS_H
