#include <boost/python/class.hpp>
#include <boost/python/def.hpp>


#include <scitbx/math/continued_fraction.h>

namespace scitbx { namespace math { namespace boost_python {

template<typename FloatType, typename IntType>
struct continued_fraction_wrapper
{
  typedef continued_fraction<IntType> wt;
  typedef typename wt::integral_type int_t;
  typedef typename wt::rational_type rational_t;

  static void wrap() {
    using namespace boost::python;
    class_<wt>("continued_fraction", no_init)
      .def(init<int_t>())
      .def("append", &wt::append)
      .def("as_rational", &wt::as_rational)
      .def("from_real", from_real_1, arg("value"))
      .def("from_real", from_real_2, (arg("value"), arg("eps")))
      .staticmethod("from_real")
      ;
  }

  static wt from_real_1(FloatType x) {
    return wt::from_real(x);
  }

  static wt from_real_2(FloatType x, FloatType eps) {
    return wt::from_real(x, eps);
  }
};

void wrap_continued_fraction() {
  continued_fraction_wrapper<double, int>::wrap();
}

}}} // scitbx::math::boost_python
