#ifndef CCTBX_XRAY_MINIMIZATION_H
#define CCTBX_XRAY_MINIMIZATION_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/packing_order.h>
#include <scitbx/array_family/block_iterator.h>
#include <scitbx/sym_mat3.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace xray { namespace minimization {

  template <typename FloatType>
  void
  damp_shifts(
    af::const_ref<FloatType> const& previous,
    af::ref<FloatType> current,
    FloatType const& max_value)
  {
    CCTBX_ASSERT(previous.size() == current.size());
    for(std::size_t i=0; i<previous.size(); i++) {
      FloatType p = previous[i];
      FloatType c = current[i];
      FloatType delta = c-p;
      if(std::abs(delta)>max_value) {
        if(delta>=0) current[i] = p+max_value;
        if(delta< 0) current[i] = p-max_value;
      }
    }
  }

  template <typename FloatType>
  void
  truncate_shifts(
    af::ref<FloatType> shifts,
    FloatType const& min_value,
    FloatType const& max_value)
  {
    CCTBX_ASSERT(min_value < max_value);
    for(std::size_t i=0; i<shifts.size(); i++) {
      FloatType sh = shifts[i];
      if(sh>max_value) shifts[i] = max_value;
      if(sh<min_value) shifts[i] = min_value;
    }
  }

  template <typename XrayScattererType,
            typename FloatType>
  af::shared<FloatType>
  shift_scales(
    af::const_ref<XrayScattererType> const& scatterers,
    std::size_t n_parameters,
    FloatType const& site_cart,
    FloatType const& u_iso,
    FloatType const& u_cart,
    FloatType const& occupancy,
    FloatType const& fp,
    FloatType const& fdp)
  {
    BOOST_STATIC_ASSERT(packing_order_convention == 2);
    af::shared<FloatType> result(n_parameters);
    scitbx::af::block_iterator<FloatType> next_shifts(
      result.ref(), "n_parameters is too small.");
    for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
        XrayScattererType const& sc = scatterers[i_sc];
        if (sc.flags.grad_site()) {
          FloatType* sh = next_shifts(3);
          for(std::size_t i=0;i<3;i++) sh[i] = site_cart;
        }
        if (sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
          next_shifts() = u_iso;
        }
        if (sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
          FloatType* sh = next_shifts(6);
          for(std::size_t i=0;i<6;i++) sh[i] = u_cart;
        }
        if(sc.flags.grad_occupancy()) {
          next_shifts() = occupancy;
        }
        if(sc.flags.grad_fp()) {
          next_shifts() = fp;
        }
        if(sc.flags.grad_fdp()) {
          next_shifts() = fdp;
        }
    }
    CCTBX_ASSERT(next_shifts.is_at_end());
    return result;
  }

  template <typename XrayScattererType,
            typename FloatType>
  struct apply_shifts
  {
    af::shared<XrayScattererType> shifted_scatterers;
    af::shared<FloatType> u_iso_refinable_params;

    apply_shifts(
      uctbx::unit_cell const& unit_cell,
      af::const_ref<XrayScattererType> const& scatterers,
      af::const_ref<FloatType> const& shifts)
    {
      BOOST_STATIC_ASSERT(packing_order_convention == 2);
      typedef typename XrayScattererType::float_type sc_f_t;
      shifted_scatterers.reserve(scatterers.size());
      cctbx::xray::scatterer_grad_flags_counts grad_flags_counts(scatterers);
      if (grad_flags_counts.tan_u_iso != 0) {
        CCTBX_ASSERT(grad_flags_counts.u_iso != 0);
        u_iso_refinable_params.resize(scatterers.size(), 0);
      }

      FloatType* u_iso_refinable_params_ptr = u_iso_refinable_params.begin();
      scitbx::af::const_block_iterator<FloatType> next_shifts(
        shifts, "Array of shifts is too small.");
      for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
         XrayScattererType sc = scatterers[i_sc];
         if(sc.flags.grad_site()) {
            sc.site += unit_cell.fractionalize(cartesian<sc_f_t>(next_shifts(3)));
         }
         if(sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
           if(sc.flags.tan_u_iso() && sc.flags.param > 0) {
             if (sc.u_iso < 0) {
               throw error(sc.report_negative_u_iso(__FILE__, __LINE__));
             }
             FloatType pi = scitbx::constants::pi;
             FloatType u_iso_max=adptbx::b_as_u(sc.flags.param);
             FloatType u_iso_refinable_param = std::tan(pi*(sc.u_iso/u_iso_max-
                                           1./2.))+next_shifts();
             sc.u_iso = u_iso_max*(std::atan(u_iso_refinable_param)+pi/2.)/pi;
             u_iso_refinable_params_ptr[i_sc] = u_iso_refinable_param;
           }
           else {
             sc.u_iso += next_shifts();
           }
         }
         if(sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
           scitbx::sym_mat3<sc_f_t> u_cart = adptbx::u_star_as_u_cart(
             unit_cell, sc.u_star);
           u_cart += scitbx::sym_mat3<sc_f_t>(next_shifts(6));
           sc.u_star = adptbx::u_cart_as_u_star(unit_cell, u_cart);
         }
        if(sc.flags.grad_occupancy()) {
           sc.occupancy += next_shifts();
        }
        if(sc.flags.grad_fp()) {
           sc.fp += next_shifts();
        }
        if(sc.flags.grad_fdp()) {
           sc.fdp += next_shifts();
        }
        shifted_scatterers.push_back(sc);
      }
      if (!next_shifts.is_at_end()) {
        throw error("Array of shifts is too large.");
      }
    }
  };

  template <typename XrayScattererType,
            typename FloatType>
  void
  add_gradients(
    af::const_ref<XrayScattererType> const& scatterers,
    af::ref<FloatType> const& xray_gradients,
    af::const_ref<scitbx::vec3<FloatType> > const& site_gradients,
    af::const_ref<FloatType> const& u_iso_gradients,
    af::const_ref<scitbx::sym_mat3<FloatType> > const& u_aniso_gradients,
    af::const_ref<FloatType> const& occupancy_gradients)
  {
    BOOST_STATIC_ASSERT(packing_order_convention == 2);
    CCTBX_ASSERT(site_gradients.size() == 0
              || site_gradients.size() == scatterers.size());
    CCTBX_ASSERT(u_iso_gradients.size() == 0
              || u_iso_gradients.size() == scatterers.size());
    CCTBX_ASSERT(u_aniso_gradients.size() == 0
              || u_aniso_gradients.size() == scatterers.size());
    CCTBX_ASSERT(occupancy_gradients.size() == 0
              || occupancy_gradients.size() == scatterers.size());
    scitbx::af::block_iterator<FloatType> next_xray_gradients(
      xray_gradients, "Array of xray gradients is too small.");
    for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
        XrayScattererType const& sc = scatterers[i_sc];
        if(sc.flags.grad_site()) {
          FloatType* xg = next_xray_gradients(3);
          if (site_gradients.size() != 0) {
            scitbx::vec3<FloatType> const& grsg = site_gradients[i_sc];
            for(std::size_t i=0;i<3;i++) xg[i] += grsg[i];
          }
        }
        if(sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
          FloatType& xg = next_xray_gradients();
          if (u_iso_gradients.size() != 0) {
            xg += u_iso_gradients[i_sc];
          }
        }
        if(sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
          FloatType* xg = next_xray_gradients(6);
          if (u_aniso_gradients.size() != 0) {
            scitbx::sym_mat3<FloatType> const& gu = u_aniso_gradients[i_sc];
            for(std::size_t i=0;i<6;i++) xg[i] += gu[i];
          }
        }
        if(sc.flags.grad_occupancy()) {
          FloatType& xg = next_xray_gradients();
          if (occupancy_gradients.size() != 0) {
            xg += occupancy_gradients[i_sc];
          }
        }
        if (sc.flags.grad_fp()) {
          next_xray_gradients();
        }
        if (sc.flags.grad_fdp()) {
          next_xray_gradients();
        }
    }
    if (!next_xray_gradients.is_at_end()) {
      throw error("Array of xray gradients is too large.");
    }
  }

  template <typename XrayScattererType,
            typename FloatType>
  af::shared<scitbx::vec3<FloatType> >
  extract_site_gradients(
    af::const_ref<XrayScattererType> const& scatterers,
    af::const_ref<FloatType> const& xray_gradients)
  {
    cctbx::xray::scatterer_grad_flags_counts grad_flags_counts(scatterers);
    CCTBX_ASSERT(grad_flags_counts.site != 0);
    BOOST_STATIC_ASSERT(packing_order_convention == 2);
    af::shared<scitbx::vec3<FloatType> > result(
      (af::reserve(scatterers.size())));
    scitbx::af::const_block_iterator<FloatType> next_xray_gradients(
      xray_gradients, "Array of xray gradients is too small.");
    for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
      XrayScattererType const& sc = scatterers[i_sc];
      const FloatType* xg = next_xray_gradients(3);
      scitbx::vec3<FloatType> grsg;
      for(std::size_t i=0;i<3;i++) grsg[i] = xg[i];
      result.push_back(grsg);
      if (sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
        next_xray_gradients();
      }
      if (sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
        next_xray_gradients(6);
      }
      if (sc.flags.grad_occupancy()) {
        next_xray_gradients();
      }
      if (sc.flags.grad_fp()) {
        next_xray_gradients();
      }
      if (sc.flags.grad_fdp()) {
        next_xray_gradients();
      }
    }
    if (!next_xray_gradients.is_at_end()) {
      throw error("Array of xray gradients is too large.");
    }
    return result;
  }

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_MINIMIZATION_H
