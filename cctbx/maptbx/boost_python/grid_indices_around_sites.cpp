#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cctbx/maptbx/grid_indices_around_sites.h>
#include <boost/unordered_set.hpp>
#include <vector>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct grid_indices_around_sites_unordered_set :
         grid_indices_around_sites_enumerator
  {
    boost::unordered_set<unsigned> buffer;

    virtual
    void
    next_point(
      std::size_t i_grid)
    {
      buffer.insert(i_grid);
    }
  };

  boost::shared_ptr<std::vector<unsigned> >
  grid_indices_around_sites_wrapper(
    uctbx::unit_cell const& unit_cell,
    af::tiny<int, 3> const& fft_n_real,
    af::tiny<int, 3> const& fft_m_real,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<double> const& site_radii)
  {
    grid_indices_around_sites_unordered_set gias;
    gias.enumerate(unit_cell, fft_n_real, fft_m_real, sites_cart, site_radii);
    boost::shared_ptr<std::vector<unsigned> >
      result(new std::vector<unsigned>());
    std::size_t n = gias.buffer.size();
    result->reserve(n);
    boost::unordered_set<unsigned>::const_iterator p = gias.buffer.begin();
    for(std::size_t i=0;i!=n;i++) result->push_back(*p++);
    std::sort(result->begin(), result->end());
    return result;
  }

} // namespace <anonymous>

  void wrap_grid_indices_around_sites()
  {
    using namespace boost::python;
    def("grid_indices_around_sites",
      grid_indices_around_sites_wrapper, (
        arg("unit_cell"),
        arg("fft_n_real"),
        arg("fft_m_real"),
        arg("sites_cart"),
        arg("site_radii")));
  }

}}} // namespace cctbx::maptbx::boost_python
