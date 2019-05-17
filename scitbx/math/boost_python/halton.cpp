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
      .def(init<int const& > (( arg("dimension")) ))
      .def("nth_given_base", &w_t::nth_given_base)
      .def("nth_all", &w_t::nth_all)
      ;
    }
  };

  struct square_halton_sampling_wrappers
  {
    typedef halton::square_halton_sampling<double> w_t;
    static void
    wrap()
    {
       using namespace boost::python;
       class_<w_t>("square_halton_sampling",no_init)
         .def(init<
              double const&,
              double const&,
              double const&,
              double const&>(( arg("low_x"),
                               arg("high_x"),
                               arg("low_y"),
                               arg("high_y") )))
         .def("next", &w_t::next )
         .def("__next__", &w_t::next )
         .def("start", &w_t::start )
         .def("state", &w_t::state )
         .def("set_state", &w_t::set_state )
         ;

    }

  };



} // namespace <anonymous>

namespace boost_python {

  void wrap_halton()
  {
    halton_wrappers::wrap();
    square_halton_sampling_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
