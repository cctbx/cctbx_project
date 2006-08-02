#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <scitbx/math/halton.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace scitbx { namespace math {

namespace {

  struct halton_wrappers
  {
    typedef halton::halton<double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("halton",no_init)
      .def(init<int const& > (( arg_("dimension")) ))
      .def("nth_given_base", &w_t::nth_given_base)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_halton()
  {
    halton_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
