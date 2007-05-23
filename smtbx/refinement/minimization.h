#ifndef SMTBX_XRAY_MINIMIZATION_H
#define SMTBX_XRAY_MINIMIZATION_H

#include <scitbx/array_family/block_iterator.h>
#include <cctbx/xray/packing_order.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/import_cctbx.h>
#include <smtbx/refinement/import_cctbx_xray.h>

  

namespace smtbx { namespace refinement { namespace minimization {

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
  struct reduce_gradient_as_per_special_position_constraints {
  
    af::shared<FloatType> reduced_gradients;

    reduce_gradient_as_per_special_position_constraints (
      uctbx::unit_cell const& unit_cell,
      af::const_ref<XrayScattererType> const& scatterers,
      af::ref<FloatType> const& xray_gradients
      )
    {
      BOOST_STATIC_ASSERT(packing_order_convention == 2);
      scitbx::af::block_iterator<FloatType> next_xray_gradients(
        xray_gradients, "Array of xray gradients is too small.");
      for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
          XrayScattererType const& sc = scatterers[i_sc];
          if(sc.flags.grad_site()) {
            FloatType* xg = next_xray_gradients(3);
          }
          if(sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
            FloatType& xg = next_xray_gradients();
          }
          if(sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
            FloatType* xg = next_xray_gradients(6);
          }
          if(sc.flags.grad_occupancy()) {
            FloatType& xg = next_xray_gradients();
          }
          if (sc.flags.grad_fp()) {
            FloatType& xg = next_xray_gradients();
          }
          if (sc.flags.grad_fdp()) {
            FloatType& xg = next_xray_gradients();
          }
      }
      if (!next_xray_gradients.is_at_end()) {
        throw error("Array of xray gradients is too large.");
      }
    }
  };

}}} // namespace smtbx::refinement::targets

#endif // SMTBX_XRAY_MINIMIZATION_H
