#ifndef CCTBX_XRAY_CURVATURES_SIMPLE_H
#define CCTBX_XRAY_CURVATURES_SIMPLE_H

#include <cctbx/xray/scattering_type_registry.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/matrix/packed.h>

namespace cctbx { namespace xray { namespace structure_factors {
namespace curvatures_simple {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::sym_mat3;

  static const double
    mtps = -scitbx::constants::two_pi * scitbx::constants::pi;

  template <typename FloatType=double>
  struct d2f_d_params_diag
  {
    typedef FloatType f_t;
    typedef std::complex<f_t> c_t;
    vec3<f_t> tphkl;
    af::tiny<f_t, 6> tphkl_outer;
    af::tiny<f_t, 21> d2_exp_huh_d_u_star_u_star;
    af::tiny<c_t, 12> diag_values;
    unsigned diag_size;

    d2f_d_params_diag(
      miller::index<> const& hkl)
    :
      tphkl(hkl)
    {
      tphkl *= scitbx::constants::two_pi;
      af::versa<f_t, af::c_grid<2> >
        outer_full = scitbx::matrix::outer_product(
          tphkl.const_ref(),
          tphkl.const_ref());
      af::shared<f_t>
        outer_u = scitbx::matrix::symmetric_as_packed_u(
          outer_full.const_ref());
      CCTBX_ASSERT(outer_u.size() == tphkl_outer.size());
      std::copy(
        outer_u.begin(), outer_u.end(), tphkl_outer.begin());
      int h = hkl[0];
      int k = hkl[1];
      int l = hkl[2];
      af::tiny<f_t, 6> d_exp_huh_d_u_star(h*h, k*k, l*l, 2*h*k, 2*h*l, 2*k*l);
      outer_full = scitbx::matrix::outer_product(
        d_exp_huh_d_u_star.const_ref(),
        d_exp_huh_d_u_star.const_ref());
      outer_u = scitbx::matrix::symmetric_as_packed_u(
        outer_full.const_ref());
      CCTBX_ASSERT(outer_u.size() == d2_exp_huh_d_u_star_u_star.size());
      std::copy(
        outer_u.begin(), outer_u.end(), d2_exp_huh_d_u_star_u_star.begin());
    }

    template <typename ScattererType>
    void
    compute(
      sgtbx::space_group const& space_group,
      miller::index<> const& hkl,
      f_t const& d_star_sq,
      af::const_ref<ScattererType> const& scatterers,
      xray::scattering_type_registry const& scattering_type_registry,
      sgtbx::site_symmetry_table const& site_symmetry_table,
      unsigned i_scatterer)
    {
      ScattererType const& scatterer = scatterers[i_scatterer];
      CCTBX_ASSERT(scatterer.flags.use_u_iso()!=scatterer.flags.use_u_aniso());
      sgtbx::site_constraints<> const* site_constraints = 0;
      sgtbx::tensor_rank_2::constraints<> const* adp_constraints = 0;
      if (site_symmetry_table.is_special_position(i_scatterer)) {
        sgtbx::site_symmetry_ops const&
          site_symmetry_ops = site_symmetry_table.get(i_scatterer);
        site_constraints = &site_symmetry_ops.site_constraints();
        if (scatterer.flags.use_u_aniso()) {
          adp_constraints = &site_symmetry_ops.adp_constraints();
        }
      }
      f_t w = scatterer.weight();
      f_t dw = -1e20; // uninitialized
      if (!scatterer.flags.use_u_aniso()) {
        double huh = scatterer.u_iso * d_star_sq;
        dw = std::exp(mtps * huh);
      }
      typename xray::scattering_type_registry::form_factor_t const&
        gaussian = scattering_type_registry.gaussian_not_optional(
          scatterer.scattering_type);
      f_t f0 = gaussian.at_d_star_sq(d_star_sq);
      f_t ffp = f0 + scatterer.fp;
      f_t fdp = scatterer.fdp;
      c_t ff(ffp, fdp);
      unsigned d2_site_site_n = 3;
      af::tiny<c_t, 6> d2_site_site;
      d2_site_site.fill(c_t(0,0));
      c_t d2_u_iso_u_iso;
      unsigned d2_u_star_u_star_n = 6;
      af::tiny<c_t, 21> d2_u_star_u_star;
      if (!scatterer.flags.use_u_aniso()) {
        d2_u_iso_u_iso = c_t(0,0);
      }
      else {
        d2_u_star_u_star.fill(c_t(0,0));
      }
      for(std::size_t i_smx=0;i_smx<space_group.order_z();i_smx++) {
        sgtbx::rt_mx s = space_group(i_smx);
        mat3<f_t> r = s.r().as_floating_point(scitbx::type_holder<f_t>());
        vec3<f_t> s_site = s * scatterer.site;
        f_t alpha = tphkl * s_site;
        if (scatterer.flags.use_u_aniso()) {
          sym_mat3<f_t> s_u_star_s = scatterer.u_star.tensor_transform(r);
          f_t huh = (hkl * s_u_star_s) * hkl;
          dw = std::exp(mtps * huh);
        }
        c_t e(std::cos(alpha), std::sin(alpha));
        {
          af::tiny<c_t, 9> ab;
          af::tiny<c_t, 6> abat;
          scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
            /* site_gtmx */ r.transpose().begin(),
            tphkl_outer.begin(),
            3, 3,
            ab.begin(),
            abat.begin());
          abat *= (w * dw * ff * e * c_t(-1,0));
          d2_site_site += abat;
        }
        if (!scatterer.flags.use_u_aniso()) {
          d2_u_iso_u_iso += w * dw * ff * e * scitbx::fn::pow2(mtps*d_star_sq);
        }
        else {
          af::versa<f_t, af::c_grid<2> >
            u_star_gtmx = scitbx::matrix::tensor_rank_2
              ::gradient_transform_matrix(r);
          af::tiny<c_t, 36> ab;
          af::tiny<c_t, 21> abat;
          scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
            u_star_gtmx.begin(),
            d2_exp_huh_d_u_star_u_star.begin(),
            6, 6,
            ab.begin(),
            abat.begin());
          abat *= (w * dw * ff * e * mtps*mtps);
          d2_u_star_u_star += abat;
        }
      }
      unsigned i_u;
      if (site_constraints == 0) {
        i_u = 3;
      }
      else {
        i_u = site_constraints->n_independent_params();
      }
      unsigned i_occ;
      if (!scatterer.flags.use_u_aniso()) {
        i_occ = i_u + 1;
      }
      else if (adp_constraints == 0) {
        i_occ = i_u + 6;
      }
      else {
        i_occ = i_u + adp_constraints->n_independent_params();
      }
      diag_size = i_occ + 3;
      if (site_constraints != 0) {
        {
          af::const_ref<f_t, af::mat_grid>
            gsm = site_constraints->gradient_sum_matrix();
          d2_site_site_n = gsm.n_rows();
          CCTBX_ASSERT(d2_site_site_n <= 2);
          af::tiny<c_t, 6> ab;
          af::tiny<c_t, 3> abat;
          scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
            gsm.begin(),
            d2_site_site.begin(),
            d2_site_site_n, gsm.n_columns(),
            ab.begin(),
            abat.begin());
          std::copy(
            abat.begin(),
            abat.begin()+(d2_site_site_n*(d2_site_site_n+1)/2),
            d2_site_site.begin());
        }
        if (scatterer.flags.use_u_aniso()) {
          af::const_ref<f_t, af::mat_grid>
            gsm = adp_constraints->gradient_sum_matrix();
          d2_u_star_u_star_n = gsm.n_rows();
          CCTBX_ASSERT(d2_u_star_u_star_n <= 4);
          af::tiny<c_t, 24> ab;
          af::tiny<c_t, 10> abat;
          scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
            gsm.begin(),
            d2_u_star_u_star.begin(),
            d2_u_star_u_star_n, gsm.n_columns(),
            ab.begin(),
            abat.begin());
          std::copy(
            abat.begin(),
            abat.begin()+(d2_u_star_u_star_n*(d2_u_star_u_star_n+1)/2),
            d2_u_star_u_star.begin());
        }
      }
      std::fill_n(diag_values.begin(), diag_size, c_t(0,0));
      scitbx::matrix::packed_u_diagonal(
        &diag_values[0], d2_site_site.begin(), d2_site_site_n);
      if (!scatterer.flags.use_u_aniso()) {
        diag_values[i_u] = d2_u_iso_u_iso;
      }
      else {
        scitbx::matrix::packed_u_diagonal(
          &diag_values[i_u], d2_u_star_u_star.begin(), d2_u_star_u_star_n);
      }
    }

    af::shared<c_t>
    copy_diag() const
    {
      c_t const* d = diag_values.begin();
      return af::shared<c_t>(d, d+diag_size);
    }
  };

}}}} // namespace cctbx::xray::structure_factors::curvatures_simple

#endif // GUARD
