#ifndef CCTBX_MAPTBX_GRID_INDICES_AROUND_SITES_H
#define CCTBX_MAPTBX_GRID_INDICES_AROUND_SITES_H

#include <cctbx/uctbx.h>
#include <scitbx/math/modulo.h>
#include <boost/scoped_array.hpp>

namespace cctbx { namespace maptbx {

  struct grid_indices_around_sites_enumerator
  {
    virtual
    ~grid_indices_around_sites_enumerator() {}

    virtual
    void
    next_point(
      std::size_t i_grid) = 0;

    void
    enumerate(
      uctbx::unit_cell const& unit_cell,
      af::tiny<int, 3> const& fft_n_real,
      af::tiny<int, 3> const& fft_m_real,
      af::const_ref<scitbx::vec3<double> > const& sites_cart,
      af::const_ref<double> const& site_radii)
    {
      CCTBX_ASSERT(unit_cell.volume() > 0);
      CCTBX_ASSERT(fft_n_real.const_ref().all_gt(0));
      CCTBX_ASSERT(fft_m_real.const_ref().all_ge(fft_n_real.const_ref()));
      CCTBX_ASSERT(site_radii.size() == sites_cart.size());
      typedef unsigned svt;
      {
        svt m[3];
        for(unsigned i=0;i<3;i++) {
          m[i] = boost::numeric_cast<svt>(fft_m_real[i]);
        }
        if (scitbx::math::unsigned_product_leads_to_overflow(m, 3U)) {
          std::ostringstream ostr;
          ostr << " Grid: " << m[0] << ", " << m[1] << ", " << m[2];
          throw std::runtime_error(
            "product of fft_m_real grid dimensions leads to unsigned"
            " integer overflow (i.e. the grid is too large)."+ostr.str());
        }
      }
      scitbx::vec3<double> fft_n_real_f;
      scitbx::vec3<double> one_over_grid_spacing;
      for(unsigned i=0;i<3;i++) {
        fft_n_real_f[i] = boost::numeric_cast<double>(fft_n_real[i]);
        double unit_cell_parameter = unit_cell.parameters()[i];
        CCTBX_ASSERT(unit_cell_parameter > 0);
        one_over_grid_spacing[i] = fft_n_real_f[i] / unit_cell_parameter;
      }
      svt n0 = boost::numeric_cast<svt>(fft_n_real[0]);
      svt n1 = boost::numeric_cast<svt>(fft_n_real[1]);
      svt n2 = boost::numeric_cast<svt>(fft_n_real[2]);
      svt m1 = boost::numeric_cast<svt>(fft_m_real[1]);
      svt m2 = boost::numeric_cast<svt>(fft_m_real[2]);
      for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
        scitbx::vec3<double> uvw = unit_cell.fractionalize(sites_cart[i_site]);
        for(unsigned i=0;i<3;i++) uvw[i] *= fft_n_real_f[i];
        double site_radius = site_radii[i_site];
        af::tiny<svt, 3> box_min;
        af::tiny<svt, 3> box_max;
        for(unsigned i=0;i<3;i++) {
          double grid_radius = std::min(
            site_radius * one_over_grid_spacing[i],
            fft_n_real_f[i] * (0.5+1.e-5));
          typedef scitbx::math::float_int_conversions<double, int> fic;
          int box_min_i = fic::ifloor(uvw[i] - grid_radius);
          int box_max_i = fic::iceil( uvw[i] + grid_radius);
          int box_min_mod = scitbx::math::mod_positive(
            box_min_i, fft_n_real[i]);
          int shift = box_min_mod - box_min_i;
          uvw[i] += boost::numeric_cast<double>(shift);
          box_min[i] = boost::numeric_cast<svt>(box_min_i + shift);
          box_max[i] = boost::numeric_cast<svt>(box_max_i + shift);
        }
        double site_radius_sq = site_radius * site_radius;
        scitbx::mat3<double> orth = unit_cell.orthogonalization_matrix();
        CCTBX_ASSERT(orth[3] == 0);
        CCTBX_ASSERT(orth[6] == 0);
        CCTBX_ASSERT(orth[7] == 0);
        orth[0] /= fft_n_real_f[0];
        orth[1] /= fft_n_real_f[1];
        orth[4] /= fft_n_real_f[1];
        orth[2] /= fft_n_real_f[2];
        orth[5] /= fft_n_real_f[2];
        orth[8] /= fft_n_real_f[2];
        svt box_n1_2 = (box_max[1] + 1 - box_min[1]) * 2;
        svt box_n2_3 = (box_max[2] + 1 - box_min[2]) * 3;
        boost::scoped_array<double> buffer(new double[box_n1_2 + box_n2_3]);
        double* o14d1_beg = buffer.get();
        double* o258d2_beg = o14d1_beg + box_n1_2;
        const double* o258d2_end = o258d2_beg + box_n2_3;
        {
          double* o14d1 = o14d1_beg;
          for(svt b1=box_min[1];b1<=box_max[1];b1++) {
            double diff1 = boost::numeric_cast<double>(b1) - uvw[1];
            *o14d1++ = orth[1] * diff1;
            *o14d1++ = orth[4] * diff1;
          }
        }
        {
          double* o258d2 = o258d2_beg;
          for(svt b2=box_min[2];b2<=box_max[2];b2++) {
            double diff2 = boost::numeric_cast<double>(b2) - uvw[2];
            *o258d2++ =                  orth[2] * diff2;
            *o258d2++ =                  orth[5] * diff2;
            *o258d2++ = scitbx::fn::pow2(orth[8] * diff2);
          }
        }
        svt g0 = box_min[0];
        svt g1b = box_min[1];
        svt g2b = box_min[2];
        for(svt b0=box_min[0];b0<=box_max[0];b0++,g0++) {
          double diff0 = boost::numeric_cast<double>(b0) - uvw[0];
          double orth0_diff0 = orth[0] * diff0;
          if (g0 == n0) g0 = 0;
          svt gm0 = g0 * m1;
          svt g1 = g1b;
          for(const double* o14d1=o14d1_beg;o14d1!=o258d2_beg;g1++) {
            double orth0_diff0_orth1_diff_1 = orth0_diff0 + *o14d1++;
            double orth4_diff1 = *o14d1++;
            if (g1 == n1) g1 = 0;
            svt gm01 = (gm0 + g1) * m2;
            svt g2 = g2b;
            for(const double* o258d2=o258d2_beg;o258d2!=o258d2_end;g2++) {
              if (g2 == n2) g2 = 0;
              double dc_sq = orth0_diff0_orth1_diff_1 + *o258d2++;
              dc_sq *= dc_sq;
              double dc1 = orth4_diff1 + *o258d2++;
              dc_sq += dc1 * dc1;
              dc_sq += *o258d2++;
              if (dc_sq > site_radius_sq) continue;
              next_point(gm01 + g2);
            }
          }
        }
      }
    }
  };

}} // namespace cctbx::maptbx

#endif // GUARD
