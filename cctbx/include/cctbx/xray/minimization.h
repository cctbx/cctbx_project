#ifndef CCTBX_XRAY_MINIMIZATION_H
#define CCTBX_XRAY_MINIMIZATION_H

#include <cctbx/xray/scattering_dictionary.h>
#include <cctbx/xray/gradient_flags.h>
#include <cctbx/xray/packing_order.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/block_iterator.h>

namespace cctbx { namespace xray { namespace minimization {

  template <typename XrayScattererType,
            typename FloatType>
  af::shared<XrayScattererType>
  apply_shifts(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group_type const& space_group_type,
    af::const_ref<XrayScattererType> const& scatterers,
    scattering_dictionary const& scattering_dict,
    xray::gradient_flags const& gradient_flags,
    af::const_ref<FloatType> const& shifts,
    FloatType const& d_min)
  {
    CCTBX_ASSERT(scattering_dict.n_scatterers() == scatterers.size());
    BOOST_STATIC_ASSERT(packing_order_convention == 1);
    typedef typename XrayScattererType::float_type sc_f_t;
    af::shared<XrayScattererType> result((af::reserve(scatterers.size())));
    sc_f_t d_star_sq_max = 0;
    if (d_min > 0) d_star_sq_max = 1 / (d_min * d_min);
    scitbx::af::const_block_iterator<FloatType> next_shifts(
      shifts, "Array of shifts is too small.");
    for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
      XrayScattererType sc = scatterers[i_sc];
      if (gradient_flags.site) {
        sc.site += unit_cell.fractionalize(cartesian<sc_f_t>(next_shifts(3)));
      }
      if (!sc.anisotropic_flag) {
        if (gradient_flags.u_iso) {
          sc.u_iso += next_shifts();
          if (sc.u_iso < 0) sc.u_iso = 0;
        }
      }
      else {
        if (gradient_flags.u_aniso) {
          scitbx::sym_mat3<sc_f_t> u_cart = adptbx::u_star_as_u_cart(
            unit_cell, sc.u_star);
          u_cart += scitbx::sym_mat3<sc_f_t>(next_shifts(6));
          sc.u_star = adptbx::u_cart_as_u_star(
            unit_cell, adptbx::eigenvalue_filtering(u_cart));
        }
      }
      if (gradient_flags.occupancy) {
        sc.occupancy += next_shifts();
        if (sc.occupancy < 0) sc.occupancy = 0;
      }
      if (gradient_flags.fp) {
        sc.fp += next_shifts();
        if (d_star_sq_max != 0) {
          sc_f_t f0 = scattering_dict
            .lookup(sc.scattering_type)
            .gaussian.at_d_star_sq(d_star_sq_max);
          if (f0 + sc.fp < 0) {
            sc.fp = -f0;
          }
        }
      }
      if (gradient_flags.fdp) {
        sc.fdp += next_shifts();
      }
      result.push_back(sc);
    }
    if (!next_shifts.is_at_end()) {
      throw error("Array of shifts is too large.");
    }
    return result;
  }

  template <typename XrayScattererType,
            typename FloatType>
  void
  add_geometry_restraints_site_gradients(
    af::const_ref<XrayScattererType> const& scatterers,
    xray::gradient_flags const& gradient_flags,
    af::ref<FloatType> const& xray_gradients,
    af::const_ref<scitbx::vec3<FloatType> > const&
      geometry_restraints_site_gradients)
  {
    CCTBX_ASSERT(
      geometry_restraints_site_gradients.size() == scatterers.size());
    CCTBX_ASSERT(gradient_flags.site == true);
    BOOST_STATIC_ASSERT(packing_order_convention == 1);
    scitbx::af::block_iterator<FloatType> next_xray_gradients(
      xray_gradients, "Array of xray gradients is too small.");
    for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
      XrayScattererType const& sc = scatterers[i_sc];
      FloatType* xg = next_xray_gradients(3);
      scitbx::vec3<FloatType> const&
        grsg = geometry_restraints_site_gradients[i_sc];
      for(std::size_t i=0;i<3;i++) xg[i] += grsg[i];
      if (!sc.anisotropic_flag) {
        if (gradient_flags.u_iso) {
          next_xray_gradients();
        }
      }
      else {
        if (gradient_flags.u_aniso) {
          next_xray_gradients(6);
        }
      }
      if (gradient_flags.occupancy) {
        next_xray_gradients();
      }
      if (gradient_flags.fp) {
        next_xray_gradients();
      }
      if (gradient_flags.fdp) {
        next_xray_gradients();
      }
    }
    if (!next_xray_gradients.is_at_end()) {
      throw error("Array of xray gradients is too large.");
    }
  }

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_MINIMIZATION_H
