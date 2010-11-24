#ifndef CCTBX_MAPTBX_STANDARD_DEVIATIONS_AROUND_SITES_HPP
#define CCTBX_MAPTBX_STANDARD_DEVIATIONS_AROUND_SITES_HPP

#include <cctbx/maptbx/grid_indices_around_sites.h>
#include <scitbx/math/mean_and_variance.h>
#include <vector>

namespace cctbx { namespace maptbx {

  struct grid_indices_around_sites_std_dev_plugin :
    grid_indices_around_sites_enumerator
  {
    double const* density_map_begin;
    std::vector<double> density_values;

    void
    clear()
    {
      density_values.clear();
    }

    virtual
    void
    next_point(
      std::size_t i_grid)
    {
      density_values.push_back(density_map_begin[i_grid]);
    }

    double
    standard_deviation() const
    {
      if (density_values.size() == 0) return 0;
      return scitbx::math::mean_and_variance<double>(
        af::const_ref<double>(
          &*density_values.begin(), density_values.size()))
            .unweighted_sample_standard_deviation();
    }
  };

  af::shared<double>
  standard_deviations_around_sites(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<double, af::flex_grid<> > const& density_map,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<double> const& site_radii)
  {
    af::flex_grid<> const& grid = density_map.accessor();
    CCTBX_ASSERT(grid.nd() == 3);
    CCTBX_ASSERT(grid.is_0_based());
    CCTBX_ASSERT(site_radii.size() == sites_cart.size());
    af::shared<double> result(
      sites_cart.size(), af::init_functor_null<double>());
    af::tiny<int, 3> fft_n_real(af::adapt_with_static_cast(grid.focus()));
    af::tiny<int, 3> fft_m_real(af::adapt_with_static_cast(grid.all()));
    grid_indices_around_sites_std_dev_plugin gias;
    gias.density_map_begin = density_map.begin();
    gias.density_values.reserve(256);
      // should avoid re-allocation in most cases
    for(std::size_t i_site=0;i_site!=sites_cart.size();i_site++) {
      gias.clear();
      gias.enumerate(
        unit_cell, fft_n_real, fft_m_real,
        af::const_ref<scitbx::vec3<double> >(&sites_cart[i_site], 1),
        af::const_ref<double>(&site_radii[i_site], 1));
      result[i_site] = gias.standard_deviation();
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // GUARD
