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
  struct grads_and_curvs_hkl_scatterer
  {
    typedef FloatType f_t;
    typedef std::complex<f_t> c_t;

    vec3<f_t> tphkl;
    af::tiny<f_t, 6> tphkl_outer;
    af::tiny<f_t, 6> d_exp_huh_d_u_star;
    af::tiny<f_t, 21> d2_exp_huh_d_u_star_u_star;
    unsigned n_params;
    af::tiny<c_t, 12> grad_values;
    af::tiny<c_t, 12> curv_values;

    grads_and_curvs_hkl_scatterer(
      miller::index<> const& hkl)
    :
      tphkl(hkl)
    {
      tphkl *= scitbx::constants::two_pi;
      af::tiny<f_t, 36> outer_buffer;
      scitbx::matrix::outer_product(
        outer_buffer.begin(),
        tphkl.const_ref(),
        tphkl.const_ref());
      scitbx::matrix::symmetric_as_packed_u(
        tphkl_outer.begin(), outer_buffer.begin(), 3);
      int h = hkl[0];
      int k = hkl[1];
      int l = hkl[2];
      d_exp_huh_d_u_star = af::tiny<f_t, 6>(h*h, k*k, l*l, 2*h*k, 2*h*l, 2*k*l);
      scitbx::matrix::outer_product(
        outer_buffer.begin(),
        d_exp_huh_d_u_star.const_ref(),
        d_exp_huh_d_u_star.const_ref());
      scitbx::matrix::symmetric_as_packed_u(
        d2_exp_huh_d_u_star_u_star.begin(), outer_buffer.begin(), 6);
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
      CCTBX_ASSERT(
        scatterer.flags.use_u_iso() != scatterer.flags.use_u_aniso());
      sgtbx::site_constraints<> const* site_constraints = 0;
      sgtbx::site_symmetry_ops::tensor_rank_2_constraints_t const* adp_constraints = 0;
      if (site_symmetry_table.is_special_position(i_scatterer)) {
        sgtbx::site_symmetry_ops const&
          site_symmetry_ops = site_symmetry_table.get(i_scatterer);
        site_constraints = &site_symmetry_ops.site_constraints();
        if (scatterer.flags.use_u_aniso()) {
          adp_constraints = &site_symmetry_ops.adp_constraints();
        }
      }
      f_t w = scatterer.weight();
      f_t wwo = scatterer.weight_without_occupancy();
      f_t dw = -1e20; // uninitialized
      if (!scatterer.flags.use_u_aniso()) {
        double huh = scatterer.u_iso * d_star_sq;
        dw = std::exp(mtps * huh);
      }
      typename xray::scattering_type_registry::gaussian_t const&
        gaussian = scattering_type_registry.gaussian_not_optional(
          scatterer.scattering_type);
      f_t f0 = gaussian.at_d_star_sq(d_star_sq);
      f_t ffp = f0 + scatterer.fp;
      f_t fdp = scatterer.fdp;
      c_t ff(ffp, fdp);
      unsigned n_params_site = 3;
      af::tiny<c_t, 3> d_site;
      af::tiny<c_t, 6> d2_site_site;
      d_site.fill(c_t(0,0));
      d2_site_site.fill(c_t(0,0));
      c_t d_u_iso;
      c_t d2_u_iso_u_iso;
      unsigned n_params_u_star = 6;
      af::tiny<c_t, 6> d_u_star;
      af::tiny<c_t, 21> d2_u_star_u_star;
      if (!scatterer.flags.use_u_aniso()) {
        d_u_iso = c_t(0,0);
        d2_u_iso_u_iso = c_t(0,0);
      }
      else {
        d_u_star.fill(c_t(0,0));
        d2_u_star_u_star.fill(c_t(0,0));
      }
      c_t d_occ(0,0);
      c_t d_fp(0,0);
      c_t d_fdp(0,0);
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
        c_t w_dw_e = w * dw * e;
        c_t w_dw_ff_e = w_dw_e * ff;
        {
          mat3<f_t> site_gtmx = r.transpose();
          d_site += af::tiny<c_t, 3>(site_gtmx * tphkl)
                  * (w_dw_ff_e * c_t(0,1));
          af::tiny<f_t, 9> ab;
          af::tiny<f_t, 6> abat;
          scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
            site_gtmx.begin(),
            tphkl_outer.begin(),
            3, 3,
            ab.begin(),
            abat.begin());
          d2_site_site += af::tiny<c_t, 6>(abat) * (w_dw_ff_e * c_t(-1,0));
        }
        if (!scatterer.flags.use_u_aniso()) {
          f_t mtps_d_star_sq = mtps * d_star_sq;
          c_t term = w_dw_ff_e * mtps_d_star_sq;
          d_u_iso += term;
          d2_u_iso_u_iso += term * mtps_d_star_sq;
        }
        else {
          af::tiny<f_t, 36> u_star_gtmx;
          scitbx::matrix::tensor_rank_2::gradient_transform_matrix(
            u_star_gtmx.begin(), r.begin());
          c_t w_dw_ff_e_mtps = w_dw_ff_e * mtps;
          {
            af::tiny<f_t, 6> ab;
            scitbx::matrix::multiply(
              u_star_gtmx.begin(),
              d_exp_huh_d_u_star.begin(),
              6, 6, 1,
              ab.begin());
            d_u_star += af::tiny<c_t, 6>(ab) * w_dw_ff_e_mtps;
          }
          {
            af::tiny<f_t, 36> ab;
            af::tiny<f_t, 21> abat;
            scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
              u_star_gtmx.begin(),
              d2_exp_huh_d_u_star_u_star.begin(),
              6, 6,
              ab.begin(),
              abat.begin());
            d2_u_star_u_star += af::tiny<c_t, 21>(abat)
                              * (w_dw_ff_e_mtps*mtps);
          }
        }
        d_occ += wwo * dw * ff * e;
        d_fp += w_dw_e;
        d_fdp += w_dw_e * c_t(0,1);
      }
      //
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
      unsigned i_fp = i_occ + 1;
      unsigned i_fdp = i_fp + 1;
      n_params = i_fdp + 1;
      //
      if (site_constraints != 0) {
        {
          af::const_ref<f_t, af::mat_grid>
            gsm = site_constraints->gradient_sum_matrix();
          n_params_site = gsm.n_rows();
          CCTBX_ASSERT(n_params_site <= 2);
          {
            af::tiny<c_t, 2> ab;
            scitbx::matrix::multiply(
              gsm.begin(),
              d_site.begin(),
              n_params_site, 3, 1,
              ab.begin());
            std::copy(
              ab.begin(),
              ab.begin()+n_params_site,
              d_site.begin());
          }
          {
            af::tiny<c_t, 6> ab;
            af::tiny<c_t, 3> abat;
            scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
              gsm.begin(),
              d2_site_site.begin(),
              n_params_site, gsm.n_columns(),
              ab.begin(),
              abat.begin());
            std::copy(
              abat.begin(),
              abat.begin()+(n_params_site*(n_params_site+1)/2),
              d2_site_site.begin());
          }
        }
        if (scatterer.flags.use_u_aniso()) {
          af::const_ref<f_t, af::mat_grid>
            gsm = adp_constraints->gradient_sum_matrix();
          n_params_u_star = gsm.n_rows();
          CCTBX_ASSERT(n_params_u_star <= 4);
          {
            af::tiny<c_t, 4> ab;
            scitbx::matrix::multiply(
              gsm.begin(),
              d_u_star.begin(),
              n_params_u_star, 6, 1,
              ab.begin());
            std::copy(
              ab.begin(),
              ab.begin()+n_params_u_star,
              d_u_star.begin());
          }
          {
            af::tiny<c_t, 24> ab;
            af::tiny<c_t, 10> abat;
            scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
              gsm.begin(),
              d2_u_star_u_star.begin(),
              n_params_u_star, gsm.n_columns(),
              ab.begin(),
              abat.begin());
            std::copy(
              abat.begin(),
              abat.begin()+(n_params_u_star*(n_params_u_star+1)/2),
              d2_u_star_u_star.begin());
          }
        }
      }
      //
      std::copy(&d_site[0], &d_site[n_params_site], &grad_values[0]);
      if (!scatterer.flags.use_u_aniso()) {
        grad_values[i_u] = d_u_iso;
      }
      else {
        std::copy(&d_u_star[0], &d_u_star[n_params_u_star], &grad_values[i_u]);
      }
      grad_values[i_occ] = d_occ;
      grad_values[i_fp] = d_fp;
      grad_values[i_fdp] = d_fdp;
      //
      std::fill_n(curv_values.begin(), n_params, c_t(0,0));
      scitbx::matrix::packed_u_diagonal(
        &curv_values[0], d2_site_site.begin(), n_params_site);
      if (!scatterer.flags.use_u_aniso()) {
        curv_values[i_u] = d2_u_iso_u_iso;
      }
      else {
        scitbx::matrix::packed_u_diagonal(
          &curv_values[i_u], d2_u_star_u_star.begin(), n_params_u_star);
      }
    }

    af::shared<c_t>
    copy_gradients() const
    {
      c_t const* v = grad_values.begin();
      return af::shared<c_t>(v, v+n_params);
    }

    af::shared<c_t>
    copy_curvatures() const
    {
      c_t const* v = curv_values.begin();
      return af::shared<c_t>(v, v+n_params);
    }
  };

  template <typename FloatType=double>
  struct grads_and_curvs_target
  {
    typedef FloatType f_t;
    typedef std::complex<f_t> c_t;

    af::shared<f_t> grads;
    af::shared<f_t> curvs;

    template <typename ScattererType>
    grads_and_curvs_target(
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      af::const_ref<ScattererType> const& scatterers,
      xray::scattering_type_registry const& scattering_type_registry,
      sgtbx::site_symmetry_table const& site_symmetry_table,
      af::const_ref<miller::index<> > const& miller_indices,
      af::const_ref<c_t> const& da_db,
      af::const_ref<scitbx::vec3<f_t> > const& daa_dbb_dab)
    {
      CCTBX_ASSERT(da_db.size() == miller_indices.size());
      CCTBX_ASSERT(daa_dbb_dab.size() == miller_indices.size());
      unsigned n_scatterers = boost::numeric_cast<unsigned>(scatterers.size());
      //
      unsigned n_params_total = 0;
      for(unsigned i_scatterer=0;i_scatterer<n_scatterers;i_scatterer++) {
        ScattererType const& scatterer = scatterers[i_scatterer];
        CCTBX_ASSERT(
          scatterer.flags.use_u_iso() != scatterer.flags.use_u_aniso());
        if (!site_symmetry_table.is_special_position(i_scatterer)) {
          n_params_total += 3;
          if (!scatterer.flags.use_u_aniso()) {
            n_params_total += 1;
          }
          else {
            n_params_total += 6;
          }
        }
        else {
          sgtbx::site_symmetry_ops const&
            site_symmetry_ops = site_symmetry_table.get(i_scatterer);
          n_params_total += site_symmetry_ops.site_constraints()
            .n_independent_params();
          if (!scatterer.flags.use_u_aniso()) {
            n_params_total += 1;
          }
          else {
            n_params_total += site_symmetry_ops.adp_constraints()
              .n_independent_params();
          }
        }
        n_params_total += 3; // occ, fp, fdp
      }
      grads.resize(n_params_total, 0);
      curvs.resize(n_params_total, 0);
      //
      for(std::size_t ih=0;ih<miller_indices.size();ih++) {
        miller::index<> const& hkl = miller_indices[ih];
        f_t d_star_sq = unit_cell.d_star_sq(hkl);
        f_t dah = da_db[ih].real();
        f_t dbh = da_db[ih].imag();
        f_t daah = daa_dbb_dab[ih][0];
        f_t dbbh = daa_dbb_dab[ih][1];
        f_t dabh = daa_dbb_dab[ih][2];
        grads_and_curvs_hkl_scatterer<FloatType> gachs(hkl);
        unsigned i_param_total = 0;
        for(unsigned i_scatterer=0;i_scatterer<n_scatterers;i_scatterer++) {
          gachs.compute(
            space_group, hkl, d_star_sq,
            scatterers, scattering_type_registry, site_symmetry_table,
            i_scatterer);
          for(unsigned ip=0;ip<gachs.n_params;ip++) {
            c_t const& g = gachs.grad_values[ip];
            c_t const& c = gachs.curv_values[ip];
            grads[i_param_total]
              += dah * g.real() + dbh * g.imag();
            curvs[i_param_total]
              += daah * g.real() * g.real()
              +  dbbh * g.imag() * g.imag()
              +  dabh * 2 * g.real() * g.imag()
              +  dah * c.real()
              +  dbh * c.imag();
            i_param_total++;
          }
        }
        CCTBX_ASSERT(i_param_total == n_params_total);
      }
    }
  };

}}}} // namespace cctbx::xray::structure_factors::curvatures_simple

#endif // GUARD
