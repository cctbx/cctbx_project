#ifndef CCTBX_MAPTBX_GRID_INDICES_AROUND_SITES_H
#define CCTBX_MAPTBX_GRID_INDICES_AROUND_SITES_H

#include <cctbx/uctbx.h>
#include <scitbx/math/modulo.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/noncopyable.hpp>
#include <set>

namespace cctbx { namespace maptbx {

  std::auto_ptr<std::set<unsigned> >
  grid_indices_around_sites(
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
    std::auto_ptr<std::set<unsigned> > result(new std::set<unsigned>);
    scitbx::vec3<double> fft_n_real_f;
    scitbx::vec3<double> one_over_grid_spacing;
    for(unsigned i=0;i<3;i++) {
      fft_n_real_f[i] = boost::numeric_cast<double>(fft_n_real[i]);
      double unit_cell_parameter = unit_cell.parameters()[i];
      CCTBX_ASSERT(unit_cell_parameter > 0);
      one_over_grid_spacing[i] = fft_n_real_f[i] / unit_cell_parameter;
    }
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      scitbx::vec3<double> uvw = unit_cell.fractionalize(sites_cart[i_site]);
      for(unsigned i=0;i<3;i++) uvw[i] *= fft_n_real_f[i];
      double site_radius = site_radii[i_site];
      af::tiny<int, 3> box_min;
      af::tiny<int, 3> box_max;
      for(unsigned i=0;i<3;i++) {
        double grid_radius = std::min(
          site_radius * one_over_grid_spacing[i],
          fft_n_real_f[i] * (0.5+1.e-5));
        typedef scitbx::math::float_int_conversions<double, int> fic;
        box_min[i] = fic::ifloor(uvw[i] - grid_radius);
        box_max[i] = fic::iceil( uvw[i] + grid_radius);
        int box_min_mod = scitbx::math::mod_positive(box_min[i], fft_n_real[i]);
        int shift = box_min_mod - box_min[i];
        uvw[i] += boost::numeric_cast<double>(shift);
        box_min[i] += shift;
        box_max[i] += shift;
      }
      double site_radius_sq = site_radius * site_radius;
      af::tiny<int, 3> b;
      scitbx::vec3<double> bf;
      af::tiny<int, 3> g;
      g[0] = box_min[0];
      for(b[0]=box_min[0];b[0]<=box_max[0];b[0]++,g[0]++) {
        bf[0] = boost::numeric_cast<double>(b[0]);
        if (g[0] == fft_n_real[0]) g[0] = 0;
        g[1] = box_min[1];
        for(b[1]=box_min[1];b[1]<=box_max[1];b[1]++,g[1]++) {
          bf[1] = boost::numeric_cast<double>(b[1]);
          if (g[1] == fft_n_real[1]) g[1] = 0;
          g[2] = box_min[2];
          for(b[2]=box_min[2];b[2]<=box_max[2];b[2]++,g[2]++) {
            bf[2] = boost::numeric_cast<double>(b[2]);
            if (g[2] == fft_n_real[2]) g[2] = 0;
            scitbx::vec3<double> diff;
            for(unsigned i=0;i<3;i++) {
              diff[i] = (bf[i] - uvw[i]) / fft_n_real_f[i];
            }
            scitbx::vec3<double> diff_cart = unit_cell.orthogonalize(diff);
            if (diff_cart.length_sq() > site_radius_sq) continue;
            unsigned grid_index = boost::numeric_cast<unsigned>(
              (g[0] * fft_m_real[1] + g[1]) * fft_m_real[2] + g[2]);
            result->insert(grid_index);
          }
        }
      }
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // GUARD
