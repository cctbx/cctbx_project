#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/matrix/tests.h>

namespace scitbx { namespace matrix { namespace boost_python {

    void wrap_matrix() {
      using namespace boost::python;
      double epsilon = std::numeric_limits<double>::epsilon();
      def("matrix_normality_ratio", normality_ratio<double>, (
        arg("a"), arg("epsilon")=epsilon));
      def("matrix_equality_ratio", equality_ratio<double>, (
        arg("a"), arg("b"), arg("epsilon")=epsilon));
      def("matrix_cholesky_test_ratio", cholesky_test_ratio<double>, (
        arg("a"), arg("x"), arg("b"), arg("epsilon")=epsilon));
    }

}}}
