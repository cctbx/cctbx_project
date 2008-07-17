#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SPECIAL_POSITION_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SPECIAL_POSITION_H

#include <utility>
#include <map>

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <cctbx/xray/scatterer_flags.h>

#include <smtbx/error.h>
#include <smtbx/refinement/parameter_map.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/import_cctbx.h>

namespace smtbx { namespace refinement { namespace constraints {

template<typename FloatType, class XrayScattererType,
         template<class> class SharedArray1D=af::shared>
class special_positions
{
  public:
    typedef FloatType float_type;
    typedef XrayScattererType xray_scatterer_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;

    special_positions(uctbx::unit_cell const &unit_cell_,
                      sgtbx::site_symmetry_table const &site_symmetry_table_,
                      af::shared<xray_scatterer_type> scatterers_,
                      parameter_map_type const &crystallographic_parameter_map_,
                      af::ref<xray::scatterer_flags> constraint_flags)
      : unit_cell(unit_cell_),
        site_symmetry_table(site_symmetry_table_),
        scatterers(scatterers_),
        crystallographic_parameter_map(crystallographic_parameter_map_)
    {
      int j = 0;
      for(int i_sc=0; i_sc < crystallographic_parameter_map.size(); ++i_sc) {
        sgtbx::site_symmetry_ops const &ops = site_symmetry_table.get(i_sc);
        if (ops.is_point_group_1()) continue;
        xray_scatterer_type const &sc = scatterers[i_sc];
        xray::scatterer_flags &flags = constraint_flags[i_sc];
        parameter_indices const &param_ids
          = crystallographic_parameter_map[i_sc];
        if (sc.flags.grad_site()) {
          if (flags.grad_site()) {
            flags.set_grad_site(false);
            site_shift_map.push_back(scatterer_indices(i_sc, j));
            int p = ops.site_constraints().n_independent_params();
            j += p;
          }
          else already_constrained_[i_sc] = flags;
        }
        if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
          if (flags.grad_u_aniso()) {
            flags.set_grad_u_aniso(false);
            adp_shift_map.push_back(scatterer_indices(i_sc, j));
            int p = ops.adp_constraints().n_independent_params();
            j += p;
          }
          else already_constrained_[i_sc] = flags;
        }
      }
    }

    std::map<int, xray::scatterer_flags>
    already_constrained() { return already_constrained_; }

    void compute_gradients(
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparameterization_gradients)
    {
      SMTBX_ASSERT(crystallographic_gradients.size()
                   == crystallographic_parameter_map.n_parameters())
                   (crystallographic_gradients.size())
                   (crystallographic_parameter_map.n_parameters());
      begin_grad_idx = reparameterization_gradients.size();
      for (int i=0; i < site_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = site_shift_map[i];
        parameter_indices const &param_ids
          = crystallographic_parameter_map[shift_ids.i_sc];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        scitbx::vec3<float_type>
          site_grad(&crystallographic_gradients[param_ids.site]);
        af::small<float_type, 3> site_ind_grad
          = ops.site_constraints()
               .independent_gradients(site_grad.const_ref());
        reparameterization_gradients.extend(site_ind_grad.begin(),
                                            site_ind_grad.end());
      }
      for (int i=0; i < adp_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = adp_shift_map[i];
        parameter_indices const &param_ids
          = crystallographic_parameter_map[shift_ids.i_sc];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        scitbx::sym_mat3<float_type>
          adp_grad(&crystallographic_gradients[param_ids.u_aniso]);
        af::small<float_type, 6> adp_ind_grad
          = ops.adp_constraints().independent_gradients(adp_grad);
        reparameterization_gradients.extend(adp_ind_grad.begin(),
                                            adp_ind_grad.end());
      }
      end_grad_idx = reparameterization_gradients.size();
    }

    void apply_shifts(
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts)
    {
      SMTBX_ASSERT(reparametrization_shifts.size() >= end_grad_idx)
                  (reparametrization_shifts.size())(end_grad_idx);
      for (int i=0; i < site_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = site_shift_map[i];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        sgtbx::site_constraints<float_type> const &ct = ops.site_constraints();
        int p = begin_grad_idx + shift_ids.i_shift;
        int q = p + ct.n_independent_params();
        if (p == q) continue;
        af::small<float_type, 3> ind_delta(&reparametrization_shifts[p],
                                           &reparametrization_shifts[q]);
        af::small<float_type, 3> ind_site
          = ct.independent_params(scatterers[shift_ids.i_sc].site);
        ind_site += ind_delta;
        fractional<float_type> new_site = ct.all_params(ind_site);
        scatterers[shift_ids.i_sc].site = new_site;
      }
      for (int i=0; i < adp_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = adp_shift_map[i];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        sgtbx::tensor_rank_2::constraints<float_type> const &ct
          = ops.adp_constraints();
        int p = begin_grad_idx + shift_ids.i_shift;
        int q = p + ct.n_independent_params();
        if (p == q) continue;
        af::small<float_type, 6> ind_delta(&reparametrization_shifts[p],
                                           &reparametrization_shifts[q]);
        af::small<float_type, 6> ind_adp
          = ct.independent_params(scatterers[shift_ids.i_sc].u_star);
        ind_adp += ind_delta;
        scitbx::sym_mat3<float_type> new_adp = ct.all_params(ind_adp);
        scatterers[shift_ids.i_sc].u_star = new_adp;
      }
    }

    void place_constrained_scatterers()
    {
      for (int i=0; i < site_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = site_shift_map[i];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        sgtbx::rt_mx const &g = ops.special_op();
        fractional<float_type> new_site = g*scatterers[shift_ids.i_sc].site;
        scatterers[shift_ids.i_sc].site = new_site;
      }
      for (int i=0; i < adp_shift_map.size(); ++i) {
        scatterer_indices const &shift_ids = adp_shift_map[i];
        sgtbx::site_symmetry_ops const &ops
          = site_symmetry_table.get(shift_ids.i_sc);
        scitbx::sym_mat3<float_type> new_adp = ops.average_u_star(
          scatterers[shift_ids.i_sc].u_star);
        scatterers[shift_ids.i_sc].u_star = new_adp;
      }
    }

  private:
    struct scatterer_indices {
      int i_sc, i_shift;
      scatterer_indices(int scatt_idx, int shift_idx)
        : i_sc(scatt_idx), i_shift(shift_idx)
      {}
    };

    uctbx::unit_cell const &unit_cell;
    sgtbx::site_symmetry_table const &site_symmetry_table;
    af::shared<xray_scatterer_type> scatterers;
    parameter_map_type const &crystallographic_parameter_map;

    af::shared<scatterer_indices> site_shift_map, adp_shift_map;
    std::map<int, xray::scatterer_flags> already_constrained_;
    unsigned begin_grad_idx, end_grad_idx;
};

}}}

#endif // GUARD
