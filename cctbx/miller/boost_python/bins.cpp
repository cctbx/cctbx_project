#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/bins.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct binning_wrappers
  {
    typedef binning w_t;

    static boost::python::object
    range_used(w_t const& o)
    {
      return scitbx::boost_python::range(1, o.n_bins_used());
    }

    static boost::python::object
    range_all(w_t const& o)
    {
      return scitbx::boost_python::range(o.n_bins_all());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("binning", no_init)
        .def(init<uctbx::unit_cell const&,
                  std::size_t,
                  double,
                  double,
                  optional<double> >())
        .def(init<uctbx::unit_cell const&,
                  std::size_t,
                  af::const_ref<index<> > const&,
                  optional<double, double, double> >())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("n_bins_used", &w_t::n_bins_used)
        .def("n_bins_all", &w_t::n_bins_all)
        .def("range_used", range_used)
        .def("range_all", range_all)
        .def("i_bin_d_too_large", &w_t::i_bin_d_too_large)
        .def("i_bin_d_too_small", &w_t::i_bin_d_too_small)
        .def("d_max", &w_t::d_max)
        .def("d_min", &w_t::d_min)
        .def("bin_d_range", &w_t::bin_d_range)
        .def("bin_d_min", &w_t::bin_d_min)
        .def("limits", &w_t::limits, ccr())
        .def("get_i_bin",
          (std::size_t(w_t::*)(double) const)
          &w_t::get_i_bin)
        .def("get_i_bin",
          (std::size_t(w_t::*)(index<> const&) const)
          &w_t::get_i_bin)
      ;
    }
  };

  struct binner_wrappers
  {
    typedef binner w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t, bases<binning> >("binner", no_init)
        .def(init<binning const&,
                  af::shared<index<> > const&>())
        .def("miller_indices", &w_t::miller_indices, ccr())
        .def("bin_indices", &w_t::bin_indices, ccr())
        .def("count", &w_t::count)
        .def("counts", &w_t::counts)
        .def("selection", &w_t::selection)
        .def("array_indices", &w_t::array_indices)
        .def("bin_centers",
          (af::shared<double>(w_t::*)(
            double const&) const)
          &w_t::bin_centers)
        .def("interpolate",
          (af::shared<double>(w_t::*)(
            af::const_ref<double> const&,
            double const&) const)
          &w_t::interpolate)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_bins()
  {
    binning_wrappers::wrap();
    binner_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
