#ifndef CCTBX_XRAY_SCATTERER_UTILS_H
#define CCTBX_XRAY_SCATTERER_UTILS_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/fp_fdp.h>
#include <cctbx/eltbx/wavelengths.h>
#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/adptbx.h>
#include <cstdio>

namespace cctbx { namespace xray {

  template <typename ScattererType>
  af::shared<bool>
  is_positive_definite_u(
    af::const_ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell)
  {
    af::shared<bool> result((af::reserve(scatterers.size())));
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      result.push_back(scatterers[i_seq].is_positive_definite_u(
        unit_cell));
    }
    return result;
  }

  template <typename ScattererType>
  af::shared<bool>
  is_positive_definite_u(
    af::const_ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell,
    double u_cart_tolerance)
  {
    af::shared<bool> result((af::reserve(scatterers.size())));
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      result.push_back(scatterers[i_seq].is_positive_definite_u(
        unit_cell, u_cart_tolerance));
    }
    return result;
  }

  template <typename ScattererType>
  void
  tidy_us(
    af::ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell,
    sgtbx::site_symmetry_table const& site_symmetry_table,
    double u_min,
    double u_max,
    double anisotropy_min)
  {
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      scatterers[i_seq].tidy_u(
        unit_cell,
        site_symmetry_table.get(i_seq),
        u_min,
        u_max,
        anisotropy_min);
    }
  }

  template <typename ScattererType>
  void
  u_star_plus_u_iso(
    af::ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell)
  {
    typedef typename ScattererType::float_type float_type;
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      ScattererType& sc = scatterers[i_seq];
      if (sc.flags.use_u_iso() && sc.flags.use_u_aniso()) {
        sc.u_star += adptbx::u_iso_as_u_star(unit_cell, sc.u_iso);
      }
    }
  }

  template <typename ScattererType>
  void
  shift_us(
    af::ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell,
    double u_shift)
  {
    typedef typename ScattererType::float_type float_type;
    scitbx::sym_mat3<float_type>
      u_star_shift = adptbx::u_iso_as_u_star(unit_cell, u_shift);
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      ScattererType& sc = scatterers[i_seq];
      if (sc.flags.use_u_iso()) {
        sc.u_iso += u_shift;
      }
      else if (sc.flags.use_u_aniso()) {
        sc.u_star += u_star_shift;
      }
    }
  }

  template <typename ScattererType>
  void
  shift_us(
    af::ref<ScattererType> const& scatterers,
    uctbx::unit_cell const& unit_cell,
    double u_shift,
    af::const_ref<std::size_t> const& selection)
  {
    typedef typename ScattererType::float_type float_type;
    for(std::size_t j=0;j<selection.size();j++) {
      std::size_t i_seq=selection[j];
      ScattererType& sc = scatterers[i_seq];
      if (sc.flags.use_u_iso() & sc.flags.use_u_aniso()) {
        scitbx::sym_mat3<float_type>
          u_star_shift = adptbx::u_iso_as_u_star(unit_cell, u_shift);
        scitbx::sym_mat3<float_type>
          u_cart_total = adptbx::u_star_as_u_cart(unit_cell,
                                                  sc.u_star+u_star_shift);
        u_cart_total[0] += sc.u_iso;
        u_cart_total[1] += sc.u_iso;
        u_cart_total[2] += sc.u_iso;
        if(adptbx::is_positive_definite(u_cart_total)) {
          sc.u_iso += u_shift;
        }
      }
      else if (sc.flags.use_u_iso()) {
        double new_u_iso = sc.u_iso + u_shift;
        if(new_u_iso >= 0.0) sc.u_iso = new_u_iso;
      }
      else if (sc.flags.use_u_aniso()) {
        scitbx::sym_mat3<float_type>
          u_star_shift = adptbx::u_iso_as_u_star(unit_cell, u_shift);
        scitbx::sym_mat3<float_type> new_u_star = sc.u_star + u_star_shift;
        if(adptbx::is_positive_definite(new_u_star)) {
          sc.u_star = new_u_star;
        }
      }
    }
  }

  template <typename ScattererType>
  void
  shift_occupancies(
    af::ref<ScattererType> const& scatterers,
    double q_shift,
    af::const_ref<std::size_t> const& selection)
  {
    typedef typename ScattererType::float_type float_type;
    for(std::size_t j=0;j<selection.size();j++) {
      std::size_t i_seq=selection[j];
      ScattererType& sc = scatterers[i_seq];
      sc.occupancy = sc.occupancy + q_shift;
    }
  }

  template <typename ScattererType>
  void
  shift_occupancies(
    af::ref<ScattererType> const& scatterers,
    double q_shift)
  {
    typedef typename ScattererType::float_type float_type;
    for(std::size_t i=0;i<scatterers.size();i++) {
      ScattererType& sc = scatterers[i];
      sc.occupancy = sc.occupancy + q_shift;
    }
  }

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
    sgtbx::site_symmetry_table const& site_symmetry_table,
    af::ref<ScattererType> const& scatterers,
    double u_star_tolerance=0)
  {
    CCTBX_ASSERT(
      scatterers.size() == site_symmetry_table.indices_const_ref().size());
    af::const_ref<std::size_t>
      sp_indices = site_symmetry_table.special_position_indices().const_ref();
    for(std::size_t i_sp=0;i_sp<sp_indices.size();i_sp++) {
      std::size_t i_seq = sp_indices[i_sp];
      scatterers[i_seq].apply_symmetry_u_star(
        site_symmetry_table.get(i_seq),
        u_star_tolerance);
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
    bool assert_min_distance_sym_equiv,
    bool non_unit_occupancy_implies_min_distance_sym_equiv_zero)
  {
    if (site_symmetry_table_for_new.indices_const_ref().size() == 0) {
      CCTBX_ASSERT(scatterers.size()
                >= site_symmetry_table.indices_const_ref().size());
      for(std::size_t i=site_symmetry_table.indices_const_ref().size();
                      i<scatterers.size();
                      i++) {
        ScattererType& sc = scatterers[i];
        sgtbx::site_symmetry site_symmetry = sc.apply_symmetry(
          unit_cell,
          space_group,
          (   sc.occupancy != 1
           && non_unit_occupancy_implies_min_distance_sym_equiv_zero
             ? 0 :  min_distance_sym_equiv),
          u_star_tolerance,
          assert_min_distance_sym_equiv);
        site_symmetry_table.process(site_symmetry);
      }
    }
    else {
      CCTBX_ASSERT(!non_unit_occupancy_implies_min_distance_sym_equiv_zero);
      CCTBX_ASSERT(scatterers.size()
                == site_symmetry_table.indices_const_ref().size()
                 + site_symmetry_table_for_new.indices_const_ref().size());
      std::size_t j = 0;
      for(std::size_t i=site_symmetry_table.indices_const_ref().size();
                      i<scatterers.size();
                      i++,j++) {
        sgtbx::site_symmetry_ops const&
          site_symmetry_ops = site_symmetry_table_for_new.get(j);
        scatterers[i].apply_symmetry(site_symmetry_ops, u_star_tolerance);
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
      if (new_scatterer.flags.use_u_aniso()) {
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
          std::snprintf(buf, sizeof(buf), fmt, i_coor);
          new_scatterer.label = scatterer.label + buf;
        }
        new_scatterer.site = coordinates[i_coor];
        if (new_scatterer.flags.use_u_aniso()) {
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
      CCTBX_ASSERT(!scatterers[i].flags.use_u_aniso());
      cartesian<> c = unit_cell.orthogonalize(scatterers[i].site);
      cartesian<> rc = rotation_matrix * c;
      rot_scatterers.push_back(scatterers[i]);
      rot_scatterers[i].site = unit_cell.fractionalize(rc);
    }
    return rot_scatterers;
  }

template <typename FloatType=double>
class apply_rigid_body_shift
{
  public:
    scitbx::vec3<FloatType> center_of_mass;
    af::shared<scitbx::vec3<FloatType> > sites_cart;
    af::shared<scitbx::vec3<FloatType> > sites_frac;

    apply_rigid_body_shift() {}

    apply_rigid_body_shift(
                 af::shared<scitbx::vec3<FloatType> > sites_cart_,
                 af::shared<scitbx::vec3<FloatType> > sites_frac_,
                 scitbx::mat3<FloatType> const& rot,
                 scitbx::vec3<FloatType> const& trans,
                 af::const_ref<FloatType> const& atomic_weights,
                 uctbx::unit_cell const& unit_cell,
                 af::const_ref<std::size_t> const& selection)
    :
    center_of_mass(0,0,0), sites_cart(sites_cart_), sites_frac(sites_frac_)
    {
      CCTBX_ASSERT(sites_cart.size() == sites_frac.size());
      CCTBX_ASSERT(sites_cart.size() == atomic_weights.size());
      FloatType xcm = 0, ycm = 0, zcm = 0, weight = 0;
      for(std::size_t j=0;j<selection.size();j++) {
          std::size_t i=selection[j];
          CCTBX_ASSERT(i < sites_cart.size());
          scitbx::vec3<FloatType> const& site_cart = sites_cart[i];
          xcm += site_cart[0]*atomic_weights[i];
          ycm += site_cart[1]*atomic_weights[i];
          zcm += site_cart[2]*atomic_weights[i];
          weight += atomic_weights[i];
      }
      if (weight != 0) {
        center_of_mass = scitbx::vec3<FloatType>(
          xcm/weight,ycm/weight,zcm/weight);
      }
      scitbx::vec3<FloatType> tcm = trans + center_of_mass;
      for(std::size_t j=0;j<selection.size();j++) {
          std::size_t i=selection[j];
          scitbx::vec3<FloatType> new_site_cart =
                               rot * (sites_cart[i] - center_of_mass) + tcm;
          sites_cart[i] = new_site_cart;
          sites_frac[i] = unit_cell.fractionalize(new_site_cart);
      }
    }
};

  template <class TableType>
  struct inelastic_form_factors
  {
    template <class ScattererType>
    static void set(af::ref<ScattererType> const &scatterers,
                    eltbx::wavelengths::characteristic photon,
                    bool set_use_fp_fdp)
    {
      set(scatterers, photon.as_angstrom(), set_use_fp_fdp);
    }

    template <class ScattererType>
    static void set(af::ref<ScattererType> const &scatterers,
                    float wavelength, // in angstrom
                    bool set_use_fp_fdp)
    {
      for (int i=0; i < scatterers.size(); ++i) {
        ScattererType &sc = scatterers[i];
        if (sc.scattering_type == "H" || sc.scattering_type == "D") continue;
        TableType tb(sc.scattering_type);
        CCTBX_ASSERT(tb.is_valid());
        eltbx::fp_fdp ff_inel = tb.at_angstrom(wavelength);
        sc.fp = ff_inel.fp();
        sc.fdp = ff_inel.fdp();
        if (set_use_fp_fdp) sc.flags.set_use_fp_fdp(true);
      }
    }
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_UTILS_H
