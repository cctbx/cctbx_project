#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/merge_equivalents.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct merge_equivalents_wrappers
  {
    typedef merge_equivalents<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("merge_equivalents", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<double> const&,
                  optional<af::const_ref<double> const&> >())
        .def("indices", &w_t::indices)
        .def("data", &w_t::data)
        .def("sigmas", &w_t::sigmas)
        .def("redundancies", &w_t::redundancies)
      ;
    }
  };

  struct merge_equivalents_hl_wrappers
  {
    typedef merge_equivalents_hl<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("merge_equivalents_hl", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<hendrickson_lattman<> > const&>())
        .def("indices", &w_t::indices)
        .def("data", &w_t::data)
        .def("redundancies", &w_t::redundancies)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_merge_equivalents()
  {
    merge_equivalents_wrappers::wrap();
    merge_equivalents_hl_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
