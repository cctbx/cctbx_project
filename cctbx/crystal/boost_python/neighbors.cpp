#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/crystal/neighbors_simple.h>

namespace cctbx { namespace crystal { namespace neighbors {

namespace {

  struct simple_pair_generator_wrappers
  {
    typedef simple_pair_generator<> w_t;

    static direct_space_asu::asu_mapping_index_pair<>
    next(w_t& o)
    {
      if (o.at_end()) {
        PyErr_SetString(PyExc_StopIteration, "asu_mappings are exhausted.");
        boost::python::throw_error_already_set();
      }
      return o.next();
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("neighbors_simple_pair_generator", no_init)
        .def(init<direct_space_asu::asu_mappings<>*,
                  optional<double const&> >(
          (arg_("asu_mappings"), arg_("distance_cutoff")))
          [with_custodian_and_ward_postcall<0,1>()])
        .def("at_end", &w_t::at_end)
        .def("next", next)
        .def("__iter__", scitbx::boost_python::pass_through)
      ;
    }
  };

}} // namespace neighbors::<anoymous>

namespace boost_python {

  void wrap_neighbors()
  {
    neighbors::simple_pair_generator_wrappers::wrap();
  }

}}} // namespace cctbx::crystal::boost_python
