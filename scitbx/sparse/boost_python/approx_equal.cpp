#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>

#include <scitbx/sparse/approx_equal.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename T>
struct approx_equal_wrapper
{
  typedef approx_equal<T> wt;

  static bool vector_cmp(wt const &self, vector<T> const &a, vector<T> const &b)
  {
    return self(a, b);
  }

  static bool matrix_cmp(wt const &self, matrix<T> const &a, matrix<T> const &b)
  {
    return self(a, b);
  }

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt>(name, no_init)
      .def(init<T>(arg("tolerance")))
      .add_property("tolerance", make_getter(&wt::tolerance),
                                 make_setter(&wt::tolerance))
      .def("__call__", vector_cmp)
      .def("__call__", matrix_cmp);
  }
};

void wrap_approx_equal() {
  approx_equal_wrapper<double>::wrap("approx_equal");
}

}}} // scitbx::sparse::boost_python
