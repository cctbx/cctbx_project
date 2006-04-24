#ifndef MMTBX_BULK_SOLVENT_BULK_SOLVENT_H
#define MMTBX_BULK_SOLVENT_BULK_SOLVENT_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cmath>
#include <cctbx/miller.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>


using scitbx::vec3;
namespace af=scitbx::af;
using scitbx::mat3;
using scitbx::sym_mat3;

namespace mmtbx { namespace bulk_solvent {

template <typename DataType, typename TagType>
void
  symmetrize_mask(
      af::ref<DataType, af::c_grid<3> > const& data,
      af::const_ref<TagType, af::c_grid<3> > const& tags)
  {
    CCTBX_ASSERT(tags.accessor().all_eq(data.accessor()));
    for(std::size_t i=0;i<data.size();i++) {
      if (tags[i] < 0) continue;
      if (data[i] == 0) data[tags[i]] = 0;
    }
    for(std::size_t i=0;i<data.size();i++) {
      if (tags[i] < 0) continue;
      data[i] = data[tags[i]];
    }
  }

vec3<double> ksol_bsol_grid_search(
                             af::const_ref<double> const& fo,
                             af::const_ref< std::complex<double> > const& fc,
                             af::const_ref< std::complex<double> > const& fm,
                             sym_mat3<double> const& b_cart,
                             af::const_ref<double> const& ksol_range,
                             af::const_ref<double> const& bsol_range,
                             double const& r_ref,
                             af::const_ref<cctbx::miller::index<> > const& hkl,
                             cctbx::uctbx::unit_cell const& uc);

class target_gradients_aniso {
public:
   target_gradients_aniso(af::const_ref<double> const& fo,
                          af::const_ref< std::complex<double> > const& fc,
                          af::const_ref< std::complex<double> > const& fm,
                          sym_mat3<double> const& b_cart,
                          double const& ksol,
                          double const& bsol,
                          af::const_ref<cctbx::miller::index<> > const& hkl,
                          cctbx::uctbx::unit_cell const& uc,
                          bool const& calc_grad_u,
                          bool const& calc_grad_ksol,
                          bool const& calc_grad_bsol);
   double target() const { return tgx; }
   af::shared<double> grad_b_cart() { return gtgx_u; }
   double grad_ksol() const { return gtgx_ksol; }
   double grad_bsol() const { return gtgx_bsol; }
   double scale_target() const { return scale_tgx; }
private:
   double tgx, scale_tgx, gtgx_ksol, gtgx_bsol;
   af::shared<double> gtgx_u;
};

af::shared<double> fb_cart(sym_mat3<double> const& b_cart,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc);

double scale(af::const_ref<double> const& fo,
             af::const_ref< std::complex<double> > const& fc);

double r_factor(af::const_ref<double> const& fo,
                af::const_ref< std::complex<double> > const& fc);

double r_factor_aniso_fast(af::const_ref<double> const& fo,
                           af::const_ref< std::complex<double> > const& fc,
                           af::const_ref< std::complex<double> > const& fm,
                           sym_mat3<double> const& b_cart,
                           double const& ksol,
                           double const& bsol,
                           af::const_ref<cctbx::miller::index<> > const& hkl,
                           cctbx::uctbx::unit_cell const& uc);


class target_gradients_aniso_ml {
public:
   target_gradients_aniso_ml(af::const_ref<double> const& fo,
                            af::const_ref< std::complex<double> > const& fc,
                            af::const_ref< std::complex<double> > const& fm,
                            sym_mat3<double> const& b_cart,
                            double const& ksol,
                            double const& bsol,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc,
                            cctbx::sgtbx::space_group const& sg,
                            af::const_ref<bool> const& gradient_flags,
                            af::const_ref<double> const& alpha,
                            af::const_ref<double> const& beta,
                            double k);
   double target() const { return tgx; }
   af::shared<double> grad_b_cart() { return gtgx_u; }
   double grad_ksol() const { return gtgx_ksol; }
   double grad_bsol() const { return gtgx_bsol; }
   double grad_k() const { return gtgx_k; }
private:
   double tgx, gtgx_ksol, gtgx_bsol, gtgx_k;
   af::shared<double> gtgx_u;
};

}} // namespace mmtbx::bulk_solvent

#endif // MMTBX_BULK_SOLVENT_BULK_SOLVENT_H
