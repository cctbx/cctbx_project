#ifndef CCTBX_XRAY_SCATTERER_UTILS_H
#define CCTBX_XRAY_SCATTERER_UTILS_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cstdio>

namespace cctbx { namespace xray {

  template <typename ScattererType>
  void
  apply_symmetry_sites(
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
  apply_symmetry_u_stars(
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
        scatterers[i].apply_symmetry(
          unit_cell,
          site_symmetry_ops,
          u_star_tolerance,
          assert_is_positive_definite,
          assert_min_distance_sym_equiv);
        site_symmetry_table.process(site_symmetry_ops);
      }
    }
  }

  template <typename ScattererType>
  af::shared<ScattererType>
  change_basis(
    af::const_ref<ScattererType> const& scatterers,
    sgtbx::change_of_basis_op const& cb_op)
  {
    af::shared<ScattererType> new_scatterers(
      scatterers.begin(), scatterers.end());
    af::ref<ScattererType> new_scatterers_ref = new_scatterers.ref();
    scitbx::mat3<double> c = cb_op.c().r().as_double();
    for(std::size_t i_seq=0;i_seq<new_scatterers_ref.size();i_seq++) {
      ScattererType& new_scatterer = new_scatterers_ref[i_seq];
      new_scatterer.site = cb_op(new_scatterer.site);
      if (new_scatterer.anisotropic_flag) {
        new_scatterer.u_star = new_scatterer.u_star.tensor_transform(c);
      }
    }
    return new_scatterers;
  }

  template <typename ScattererType>
  af::shared<ScattererType>
  expand_to_p1(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group const& space_group,
    af::const_ref<ScattererType> const& scatterers,
    sgtbx::site_symmetry_table const& site_symmetry_table,
    bool append_number_to_labels)
  {
    af::shared<ScattererType> new_scatterers((af::reserve(scatterers.size())));
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      ScattererType const& scatterer = scatterers[i_seq];
      const char* fmt = 0;
      if (append_number_to_labels) {
        if      (scatterer.multiplicity() >= 1000) fmt = "_%04u";
        else if (scatterer.multiplicity() >= 100)  fmt = "_%03u";
        else if (scatterer.multiplicity() >= 10)   fmt = "_%02u";
        else                                       fmt = "_%u";
      }
      sgtbx::sym_equiv_sites<> equiv_sites(
        unit_cell,
        space_group,
        scatterer.site,
        site_symmetry_table.get(i_seq));
      af::const_ref<scitbx::vec3<double> >
        coordinates = equiv_sites.coordinates().ref();
      ScattererType new_scatterer = scatterer;
      for(unsigned i_coor=0;i_coor<coordinates.size();i_coor++) {
        if (fmt) {
          char buf[40];
          std::sprintf(buf, fmt, i_coor);
          new_scatterer.label = scatterer.label + buf;
        }
        new_scatterer.site = coordinates[i_coor];
        if (new_scatterer.anisotropic_flag) {
          scitbx::mat3<double> c = equiv_sites.sym_op(i_coor).r().as_double();
          new_scatterer.u_star = scatterer.u_star.tensor_transform(c);
        }
        new_scatterers.push_back(new_scatterer);
      }
    }
    return new_scatterers;
  }

  template <typename ScattererType>
  std::size_t
  n_undefined_multiplicities(
    af::const_ref<ScattererType> const& scatterers)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<scatterers.size();i++) {
      if (scatterers[i].multiplicity() <= 0) result += 1;
    }
    return result;
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
