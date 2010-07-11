#ifndef SMTBX_XRAY_MINIMIZATION_H
#define SMTBX_XRAY_MINIMIZATION_H

#include <scitbx/array_family/block_iterator.h>
#include <scitbx/sym_mat3.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/xray/packing_order.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <cctbx/coordinates.h>
#include <cctbx/adptbx.h>

#include <smtbx/error.h>
#include <smtbx/import_cctbx.h>
#include <smtbx/import_scitbx_af.h>

#include <iostream>


namespace smtbx { namespace refinement {

  using cctbx::xray::scatterer;
  using cctbx::xray::packing_order_convention;

namespace minimization {

  template <typename XrayScattererType,
            typename FloatType>
  struct apply_special_position_constrained_shifts
  {
    af::shared<XrayScattererType> shifted_scatterers;
    af::shared<FloatType> u_iso_refinable_params;

    apply_special_position_constrained_shifts(
      uctbx::unit_cell const& unit_cell,
      const sgtbx::site_symmetry_table& site_symmetry_table,
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
        const sgtbx::site_symmetry_ops& op = site_symmetry_table.get(i_sc);
        if(sc.flags.grad_site()) {
          const sgtbx::site_constraints<FloatType>& site_c = op.site_constraints();
          int n = site_c.n_independent_params();
          const FloatType *xg = next_shifts(n);
          if (n < 3) {
            af::small<FloatType,3> shift(xg, xg+n);
            sc.site += site_c.all_params(shift);
          }
          else {
            sc.site += fractional<sc_f_t>(xg);
          }
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
           const sgtbx::tensor_rank_2::cartesian_constraints<FloatType>&
           adp_c = op.cartesian_adp_constraints(unit_cell);
           int n = adp_c.n_independent_params();
           const FloatType *xg = next_shifts(n);
           scitbx::sym_mat3<FloatType> shift(xg);
           if (n < 6) {
             af::small<FloatType,6> shift(xg, xg+n);
             sc.u_star += adptbx::u_cart_as_u_star(
                unit_cell,
                adp_c.all_params(shift));
           }
           else {
             sc.u_star += adptbx::u_cart_as_u_star(
               unit_cell,
               scitbx::sym_mat3<FloatType>(xg));
           }
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
  struct special_position_constrained_gradients {

    af::shared<FloatType> reduced_gradients;

    special_position_constrained_gradients (
      uctbx::unit_cell const& unit_cell,
      const sgtbx::site_symmetry_table& site_symmetry_table,
      af::const_ref<XrayScattererType> const& scatterers,
      af::ref<FloatType> const& xray_gradients
      )
    {
      BOOST_STATIC_ASSERT(packing_order_convention == 2);
      scitbx::af::block_iterator<FloatType> next_xray_gradients(
        xray_gradients, "Array of xray gradients is too small.");
      reduced_gradients.reserve(2*xray_gradients.size()/3);
      for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
        XrayScattererType const& sc = scatterers[i_sc];
        const sgtbx::site_symmetry_ops& op = site_symmetry_table.get(i_sc);
        if(sc.flags.grad_site()) {
          FloatType *xg = next_xray_gradients(3);
          scitbx::vec3<FloatType> cart_grad(xg);
          scitbx::vec3<FloatType> frac_grad
            = unit_cell.fractionalize_gradient(cart_grad);
          if (op.is_point_group_1()) {
            reduced_gradients.extend(frac_grad.begin(), frac_grad.end());
          }
          else {
            const sgtbx::site_constraints<FloatType>& site_c
              = op.site_constraints();
            af::small<FloatType,3> frac_grad1
              = site_c.independent_gradients(frac_grad.const_ref());
            reduced_gradients.extend(frac_grad1.begin(), frac_grad1.end());
          }
        }
        if(sc.flags.grad_u_iso() && sc.flags.use_u_iso()) {
          FloatType& xg = next_xray_gradients();
          reduced_gradients.push_back(xg);
        }
        if(sc.flags.grad_u_aniso() && sc.flags.use_u_aniso()) {
          FloatType* xg = next_xray_gradients(6);
          scitbx::sym_mat3<FloatType> cart_grad(xg);
          if (!op.is_point_group_1()) {
            const sgtbx::tensor_rank_2::cartesian_constraints<FloatType>&
            adp_c = op.cartesian_adp_constraints(unit_cell);
            af::small<FloatType,6> ind_cart_grad =
              adp_c.independent_gradients(cart_grad);
            reduced_gradients.extend(ind_cart_grad.begin(), ind_cart_grad.end());
          }
          else {
            reduced_gradients.extend(cart_grad.begin(), cart_grad.end());
          }
        }
        if(sc.flags.grad_occupancy()) {
          FloatType& xg = next_xray_gradients();
          reduced_gradients.push_back(xg);
        }
        if (sc.flags.grad_fp()) {
          FloatType& xg = next_xray_gradients();
          reduced_gradients.push_back(xg);
        }
        if (sc.flags.grad_fdp()) {
          FloatType& xg = next_xray_gradients();
          reduced_gradients.push_back(xg);
        }
      }
      if (!next_xray_gradients.is_at_end()) {
        throw error("Array of xray gradients is too large.");
      }
    }
  };

}}} // namespace smtbx::refinement::targets

#endif // SMTBX_XRAY_MINIMIZATION_H
