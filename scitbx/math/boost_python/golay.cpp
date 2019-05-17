#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <scitbx/math/golay.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace scitbx { namespace math {

namespace {

  struct golay_24_12_generator_wrappers
  {
    typedef golay_24_12_generator<> w_t;

    static af::tiny<int, 24>
    next(w_t& o)
    {
      if (o.at_end()) {
        PyErr_SetString(PyExc_StopIteration,
          "golay_24_12_generator is exhausted.");
        boost::python::throw_error_already_set();
      }
      return o.next();
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("golay_24_12_generator")
        .def("at_end", &w_t::at_end)
        .def("next", next)
        .def("__next__", next)
        .def("__iter__", boost_python::pass_through)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_golay()
  {
    golay_24_12_generator_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
