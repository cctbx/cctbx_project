#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/crystal/neighbors_fast.h>

namespace cctbx { namespace crystal { namespace neighbors {

namespace {

  template <typename PairGeneratorType>
  struct helper
  {
    static boost::python::object
    next(PairGeneratorType& o)
    {
      if (o.at_end()) {
        PyErr_SetString(PyExc_StopIteration, "asu_mappings are exhausted.");
        boost::python::throw_error_already_set();
      }
      return boost::python::object(o.next());
    }
  };

  struct simple_pair_generator_wrappers
  {
    typedef simple_pair_generator<> w_t;

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
        .def("distance_cutoff_sq", &w_t::distance_cutoff_sq)
        .def("at_end", &w_t::at_end)
        .def("next", helper<w_t>::next)
        .def("__iter__", scitbx::boost_python::pass_through)
        .def("restart", &w_t::restart)
      ;
    }
  };

  struct fast_pair_generator_wrappers
  {
    typedef fast_pair_generator<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t, bases<simple_pair_generator<> > >(
        "neighbors_fast_pair_generator", no_init)
        .def(init<direct_space_asu::asu_mappings<>*,
                  double const&,
                  optional<double const&> >(
          (arg_("asu_mappings"),
           arg_("distance_cutoff"),
           arg_("epsilon")))
          [with_custodian_and_ward_postcall<0,1>()])
        .def("epsilon", &w_t::epsilon)
        .def("n_boxes", &w_t::n_boxes, ccr())
        .def("next", helper<w_t>::next)
        .def("restart", &w_t::restart)
      ;
    }
  };

}} // namespace neighbors::<anoymous>

namespace boost_python {

  void wrap_neighbors()
  {
    neighbors::simple_pair_generator_wrappers::wrap();
    neighbors::fast_pair_generator_wrappers::wrap();
  }

}}} // namespace cctbx::crystal::boost_python
