#ifndef CCTBX_ADP_RESTRAINTS_ANISO_H
#define CCTBX_ADP_RESTRAINTS_ANISO_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/xray/scatterer_flags.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/restraints.h>
#include <scitbx/matrix/matrix_vector_operations.h>

namespace cctbx { namespace adp_restraints {
  using namespace cctbx::geometry_restraints;
  using namespace cctbx::xray;
  namespace adptbx = cctbx::adptbx;

  inline void check_flags (
    scatterer_flags fl)
  {
    if (fl.grad_u_iso()) {
      CCTBX_ASSERT(!fl.grad_u_aniso());
      CCTBX_ASSERT(fl.use_u_iso());
      CCTBX_ASSERT(fl.use());
    } else if (fl.grad_u_aniso()) {
      CCTBX_ASSERT(!fl.grad_u_iso());
      CCTBX_ASSERT(fl.use_u_aniso());
      CCTBX_ASSERT(fl.use());
    }
  }

  class eval_adp_aniso_restraints {
    public :
      double target;
      unsigned number_of_restraints;
      af::shared<double> gradients_iso_;
      af::shared<scitbx::sym_mat3<double> > gradients_aniso_cart_;

      eval_adp_aniso_restraints (
        af::const_ref<scatterer<> > const& scatterers,
        af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
        af::const_ref<double> const& u_iso,
        af::const_ref<bond_simple_proxy> const& bond_proxies,
        af::const_ref<bool> const& selection,
        af::const_ref<bool> const& hd_selection,
        unsigned n_grad_u_iso,
        bool use_hd)
      :
        target(0.0),
        number_of_restraints(0),
        gradients_aniso_cart_(scatterers.size()),
        gradients_iso_(scatterers.size())
      {
        unsigned n_proxies = bond_proxies.size();
        for (unsigned k = 0; k < n_proxies; k++) {
          bond_simple_proxy proxy = bond_proxies[k];
          unsigned i = proxy.i_seqs[0];
          unsigned j = proxy.i_seqs[1];
          bool tmp_flag = true;
          if (!use_hd) {
            if ((hd_selection[i]) || (hd_selection[j])) {
              tmp_flag = false;
            }
          }
          if ((selection[i]) && (selection[j]) && (tmp_flag)) {
            scatterer_flags const fl_i = scatterers[i].flags;
            scatterer_flags const fl_j = scatterers[j].flags;
            check_flags(fl_i);
            check_flags(fl_j);
            if ((fl_i.use_u_aniso()) && (fl_j.use_u_aniso())) {
              scitbx::sym_mat3<double> u_i = u_cart[i];
              scitbx::sym_mat3<double> u_j = u_cart[j];
              scitbx::sym_mat3<double> g_i = gradients_aniso_cart_[i];
              scitbx::sym_mat3<double> g_j = gradients_aniso_cart_[j];
              for (unsigned i_seq = 0; i_seq < 6; i_seq++) {
                double diff = u_i[i_seq] - u_j[i_seq];
                target += diff * diff;
                if (fl_i.grad_u_aniso()) g_i[i_seq] +=  2.0 * diff;
                if (fl_j.grad_u_aniso()) g_j[i_seq] += -2.0 * diff;
                number_of_restraints++;
              }
              gradients_aniso_cart_[i] = g_i;
              gradients_aniso_cart_[j] = g_j;
            } else if ((fl_i.use_u_iso()) && (fl_j.use_u_iso())) {
              double u_i = u_iso[i];
              double u_j = u_iso[j];
              double diff = u_i - u_j;
              target += diff * diff;
              if (fl_i.grad_u_iso()) gradients_iso_[i] +=  2.0 * diff;
              if (fl_j.grad_u_iso()) gradients_iso_[j] += -2.0 * diff;
              number_of_restraints++;
            } else if ((fl_i.use_u_aniso()) && (fl_j.use_u_iso())) {
              scitbx::sym_mat3<double> u_i = u_cart[i];
              double u_j = u_iso[j];
              scitbx::sym_mat3<double> u_j_cart = adptbx::u_iso_as_u_cart(u_j);
              scitbx::sym_mat3<double> g_i = gradients_aniso_cart_[i];
              scitbx::sym_mat3<double> g_j(0,0,0,0,0,0);
              for (unsigned i_seq = 0; i_seq < 3; i_seq++) {
                double diff = u_i[i_seq] - u_j_cart[i_seq];
                target += diff * diff;
                if (fl_i.grad_u_aniso()) g_i[i_seq] +=  2.0 * diff;
                if (fl_j.grad_u_iso())   g_j[i_seq] += -2.0 * diff;
                number_of_restraints++;
              }
              if (fl_j.grad_u_iso()) {
                gradients_iso_[j] += (g_j[0] + g_j[1] + g_j[2]);
              }
              if (fl_i.grad_u_aniso()) gradients_aniso_cart_[i] = g_i;
            } else if ((fl_i.use_u_iso()) && (fl_j.use_u_aniso())) {
              double u_i = u_iso[i];
              scitbx::sym_mat3<double> u_j = u_cart[j];
              scitbx::sym_mat3<double> u_i_cart = adptbx::u_iso_as_u_cart(u_i);
              scitbx::sym_mat3<double> g_i(0,0,0,0,0,0);
              scitbx::sym_mat3<double> g_j = gradients_aniso_cart_[j];
              for (unsigned i_seq = 0; i_seq < 3; i_seq++) {
                double diff = u_i_cart[i_seq] - u_j[i_seq];
                target += diff * diff;
                if (fl_i.grad_u_iso())   g_i[i_seq] +=  2.0 * diff;
                if (fl_j.grad_u_aniso()) g_j[i_seq] += -2.0 * diff;
                number_of_restraints++;
              }
              if (fl_i.grad_u_iso()) {
                gradients_iso_[i] += (g_i[0] + g_i[1] + g_i[2]);
              }
              if (fl_j.grad_u_aniso()) gradients_aniso_cart_[j] = g_j;
            }
          }
        }
      }

      af::shared<double> gradients_iso ()
      {
        return gradients_iso_;
      }

      af::shared<scitbx::sym_mat3<double> > gradients_aniso_cart ()
      {
        return gradients_aniso_cart_;
      }
  };
}}

#endif
