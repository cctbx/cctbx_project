#ifndef CCTBX_XRAY_SCATTERER_UTILS_H
#define CCTBX_XRAY_SCATTERER_UTILS_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>

namespace cctbx { namespace xray {

  template <typename ScattererType>
  void
  apply_symmetry_site(
    sgtbx::site_symmetry_table const& site_symmetry_table,
    af::ref<ScattererType> const& scatterers)
  {
    CCTBX_ASSERT(
      scatterers.size() == site_symmetry_table.indices_const_ref().size());
    af::const_ref<std::size_t>
      sp_indices = site_symmetry_table.special_position_indices().const_ref();
    for(std::size_t i_sp=0;i_sp<sp_indices.size();i_sp++) {
      std::size_t i_seq = sp_indices[i_sp];
      scatterers[i_seq].apply_symmetry_site(site_symmetry_table.get(i_seq));
    }
  }

  template <typename ScattererType>
  void
  apply_symmetry_u_star(
    uctbx::unit_cell const& unit_cell,
    sgtbx::site_symmetry_table const& site_symmetry_table,
    af::ref<ScattererType> const& scatterers,
    double u_star_tolerance=0,
    bool assert_is_positive_definite=false,
    bool assert_min_distance_sym_equiv=true)
  {
    CCTBX_ASSERT(
      scatterers.size() == site_symmetry_table.indices_const_ref().size());
    af::const_ref<std::size_t>
      sp_indices = site_symmetry_table.special_position_indices().const_ref();
    for(std::size_t i_sp=0;i_sp<sp_indices.size();i_sp++) {
      std::size_t i_seq = sp_indices[i_sp];
      scatterers[i_seq].apply_symmetry_u_star(
        unit_cell,
        site_symmetry_table.get(i_seq),
        u_star_tolerance,
        assert_is_positive_definite,
        assert_min_distance_sym_equiv);
    }
  }

  template <typename ScattererType>
  void
  add_scatterers_ext(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group,
    af::ref<ScattererType> const& scatterers,
    sgtbx::site_symmetry_table& site_symmetry_table,
    sgtbx::site_symmetry_table const& site_symmetry_table_for_new,
    double min_distance_sym_equiv,
    double u_star_tolerance,
    bool assert_is_positive_definite,
    bool assert_min_distance_sym_equiv)
  {
    if (site_symmetry_table_for_new.indices_const_ref().size() == 0) {
      CCTBX_ASSERT(scatterers.size()
                >= site_symmetry_table.indices_const_ref().size());
      for(std::size_t i=site_symmetry_table.indices_const_ref().size();
                      i<scatterers.size();
                      i++) {
        sgtbx::site_symmetry site_symmetry = scatterers[i].apply_symmetry(
          unit_cell,
          space_group,
          min_distance_sym_equiv,
          u_star_tolerance,
          assert_is_positive_definite,
          assert_min_distance_sym_equiv);
        site_symmetry_table.process(site_symmetry);
      }
    }
    else {
      CCTBX_ASSERT(scatterers.size()
                == site_symmetry_table.indices_const_ref().size()
                 + site_symmetry_table_for_new.indices_const_ref().size());
      std::size_t j = 0;
      for(std::size_t i=site_symmetry_table.indices_const_ref().size();
                      i<scatterers.size();
                      i++,j++) {
        sgtbx::site_symmetry_ops const&
          site_symmetry_ops = site_symmetry_table_for_new.get(j);
        if (!site_symmetry_ops.is_point_group_1()) {
          scatterers[i].apply_symmetry_site(
            site_symmetry_ops);
          scatterers[i].apply_symmetry_u_star(
            unit_cell,
            site_symmetry_ops,
            u_star_tolerance,
            assert_is_positive_definite,
            assert_min_distance_sym_equiv);
        }
        site_symmetry_table.process(site_symmetry_ops);
      }
    }
  }

  template <typename AsuMappingsType, typename ScattererType>
  void
  asu_mappings_process(
    AsuMappingsType& asu_mappings,
    af::const_ref<ScattererType> const& scatterers,
    sgtbx::site_symmetry_table const& site_symmetry_table)
  {
    CCTBX_ASSERT(site_symmetry_table.indices_const_ref().size()
              == scatterers.size());
    asu_mappings.reserve(asu_mappings.mappings().size() + scatterers.size());
    for(std::size_t i=0;i<scatterers.size();i++) {
      asu_mappings.process(scatterers[i].site, site_symmetry_table.get(i));
    }
  }

  template <typename ScattererType>
  af::shared<ScattererType>
  rotate(
    uctbx::unit_cell const& unit_cell,
    scitbx::mat3<double> const& rotation_matrix,
    af::const_ref<ScattererType> const& scatterers)
  {
    af::shared<ScattererType>
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
