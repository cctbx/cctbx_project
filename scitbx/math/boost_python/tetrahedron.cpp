#include <scitbx/math/tetrahedron.h>

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>

namespace scitbx { namespace math {

namespace {

template<typename T>
struct tetrahedron_wrapper
{
  typedef tetrahedron<T> wt;

  static void wrap() {
    using namespace boost::python;
    return_value_policy<copy_const_reference> ccr;
    class_<wt>("tetrahedron", no_init)
      .def(init<typename wt::vertices_t const&>((arg("vertices"))))
      .add_property("vertices", make_function(&wt::vertices, ccr))
      .def("volume", &wt::volume)
      .def("gradients", &wt::gradients)
      ;
  }
};

}

namespace boost_python {
  void wrap_tetrahedron() {
    tetrahedron_wrapper<double>::wrap();
  }
}

}}
