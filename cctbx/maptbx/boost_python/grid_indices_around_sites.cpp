#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cctbx/maptbx/grid_indices_around_sites.h>
#include <boost/unordered_set.hpp>
#include <set>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  std::size_t
  grid_indices_around_sites_unordered(
    uctbx::unit_cell const& unit_cell,
    af::tiny<int, 3> const& fft_n_real,
    af::tiny<int, 3> const& fft_m_real,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<double> const& site_radii)
  {
    return grid_indices_around_sites<boost::unordered_set<unsigned> >(
      unit_cell, fft_n_real, fft_m_real, sites_cart, site_radii)->size();
  }

} // namespace <anonymous>

  void wrap_grid_indices_around_sites()
  {
    using namespace boost::python;
    def("grid_indices_around_sites",
      grid_indices_around_sites<std::set<unsigned> >, (
        arg_("unit_cell"),
        arg_("fft_n_real"),
        arg_("fft_m_real"),
        arg_("sites_cart"),
        arg_("site_radii")));
    def("grid_indices_around_sites_unordered",
      grid_indices_around_sites_unordered, (
        arg_("unit_cell"),
        arg_("fft_n_real"),
        arg_("fft_m_real"),
        arg_("sites_cart"),
        arg_("site_radii")));
  }

}}} // namespace cctbx::maptbx::boost_python
