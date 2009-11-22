#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/lu_factorization.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename FloatType>
struct gilbert_peierls_lu_factorization_wrapper
{
  typedef gilbert_peierls_lu_factorization<matrix<FloatType> > wt;

  static void wrap(char const *name) {
    using namespace boost::python;
    return_internal_reference<> rir;
    class_<wt>(name, no_init)
      .def(init<typename wt::matrix_type const&>(arg("matrix")))
      .def("factored", &wt::factored, rir)
      .def("l", &wt::l, rir)
      .def("u", &wt::u, rir)
      .def("rows_permutation", &wt::rows_permutation)
      ;
  }
};

void wrap_lu_factorization() {
  gilbert_peierls_lu_factorization_wrapper<double>::wrap(
    "gilbert_peierls_lu_factorization");
}

}}} // scitbx::sparse::boost_python
