#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost_adaptbx/easy_overloads.h>

#include <scitbx/matrix/tests.h>

namespace scitbx { namespace matrix { namespace boost_python {

    template <typename FloatType>
    BOOST_ADAPTBX_FUNCTION_OVERLOADS(matrix_normality_ratio_overloads,
                                     normality_ratio<FloatType>, 1, 2,
                                     boost::python::args("a", "epsilon"));

    template <typename FloatType>
    BOOST_ADAPTBX_FUNCTION_OVERLOADS(matrix_equality_ratio_overloads,
                                     equality_ratio<FloatType>, 2, 3,
                                     boost::python::args("a", "b", "epsilon"));

    void wrap_matrix() {
      using namespace matrix::boost_python;
      matrix_normality_ratio_overloads<double>::wrap("matrix_normality_ratio");
      matrix_equality_ratio_overloads<double>::wrap("matrix_equality_ratio");
    }

}}}
