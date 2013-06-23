#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/math/least_squares_plane.h>

namespace scitbx { namespace math { namespace boost_python {

template<typename ScalarType=double>
struct least_squares_plane_wrapper
{
  typedef least_squares_plane<ScalarType> wt;
  typedef typename wt::vector_type vector_type;

  static void wrap() {
    using namespace boost::python;
    return_value_policy<return_by_value> rbv;
    class_<wt>("least_squares_plane", no_init)
      .add_property("normal", make_function(&wt::normal, rbv))
      .add_property("distance_to_origin", &wt::distance_to_origin)
      .def(init<af::const_ref<vector_type> const &>((arg("points"))))
      ;
  }
};

void wrap_least_squares_plane() {
  least_squares_plane_wrapper<>::wrap();
}

}}} // scitbx::math::boost_python
