#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <omptbx/omp_or_stubs.h>

namespace omptbx { namespace boost_python { namespace {

  void
  set_num_threads(int num_threads) { omp_set_num_threads(num_threads); }

  int
  get_num_threads() { return omp_get_num_threads(); }

  int
  get_max_threads() { return omp_get_max_threads(); }

  int
  get_num_procs() { return omp_get_num_procs(); }

  int
  get_dynamic() { return omp_get_dynamic(); }

  void
  set_dynamic(int dynamic_threads) { omp_set_dynamic(dynamic_threads); }

  int
  get_nested() { return omp_get_nested(); }

  void
  set_nested(int nested) { omp_set_nested(nested); }

  void
  tst_parallel(){
    #pragma omp parallel for
    for (int ix = 0; ix<omp_get_num_threads(); ++ix){
        int tid = omp_get_thread_num();
        printf("Hello world from omp thread %d of %d\n", tid, omp_get_num_threads());
    }
  }
  boost::python::tuple
  tst_environment()
  {
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
#if defined(_OPENMP) && _OPENMP <= 199819
    omp_set_num_threads(4);
    #pragma omp parallel sections reduction(+:bar)
#else
    #pragma omp parallel sections reduction(+:bar) \
                num_threads(4)
#endif
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

  void
  init_module()
  {
    using namespace boost::python;

    scope s;
#if defined(OMPTBX_HAVE_OMP_H)
    s.attr("omp_version") = _OPENMP;
    s.attr("have_stubs_h") = false;
    s.attr("have_omp_h") = true;
#else
    s.attr("omp_version") = object();
    s.attr("have_stubs_h") = true;
    s.attr("have_omp_h") = false;
#endif

    def("omp_set_num_threads", set_num_threads);
    def("omp_get_num_threads", get_num_threads);
    def("omp_get_max_threads", get_max_threads);
    def("omp_get_num_procs", get_num_procs);
    def("omp_get_dynamic", get_dynamic);
    def("omp_set_dynamic", set_dynamic);
    def("omp_get_nested", get_nested);
    def("omp_set_nested", set_nested);

    def("tst_environment", &tst_environment);
    def("tst_parallel", &tst_parallel);
  }

}}} // namespace omptbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(omptbx_ext)
{
  omptbx::boost_python::init_module();
}
