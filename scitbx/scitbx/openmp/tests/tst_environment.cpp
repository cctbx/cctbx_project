#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <iostream>
#include <omp.h>

namespace scitbx { namespace openmp { namespace tests {
namespace boost_python { namespace {

  boost::python::tuple tst_environment() {
    int foo=0, bar=0;
    int m,n;
    #pragma omp parallel sections reduction(+:foo)
    {
      m = omp_get_num_threads();
      #pragma omp section
      {
        foo += 1;
      }
      #pragma omp section
      {
        foo += 1;
      }
      #pragma omp section
      {
        foo += 1;
      }
      #pragma omp section
      {
        foo += 1;
      }
    }
    #pragma omp parallel sections reduction(+:bar)\
                num_threads(4)
    {
      n = omp_get_num_threads();
      #pragma omp section
      {
        bar += 1;
      }
      #pragma omp section
      {
        bar += 1;
      }
      #pragma omp section
      {
        bar += 1;
      }
      #pragma omp section
      {
        bar += 1;
      }
    }
    return boost::python::make_tuple(foo,m, bar,n);
  }
} // namespace <anonymous>

void init_module() {
  using namespace boost::python;
  def("tst_environment", &tst_environment);
}

}}}} // namespace scitbx::openmp::tests::boost_python

BOOST_PYTHON_MODULE(scitbx_openmp_tests_ext) {
  scitbx::openmp::tests::boost_python::init_module();
}
