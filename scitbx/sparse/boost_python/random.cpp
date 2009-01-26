#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>

#include <scitbx/sparse/random.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename T>
struct random_matrix_generator_wrapper
{
  typedef random_matrix_generator<T> wt;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt>(name)
      .def("__call__", &wt::operator(), (
          arg("rows"), arg("cols"),
          arg("lower_bound"), arg("upper_bound"),
          arg("sparsity")))
      .add_property("non_zeroes", &wt::non_zeroes)
    ;
  }
};

void wrap_random() {
  random_matrix_generator_wrapper<double>::wrap("random_matrix_generator");
}

}}} // scitbx::sparse::boost_python
