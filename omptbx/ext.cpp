#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <omptbx/omp_or_stubs.h>

namespace omptbx { namespace boost_python { namespace {

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
    #pragma omp parallel sections reduction(+:bar) \
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

  void
  init_module()
  {
    using namespace boost::python;

    scope s;
#if defined(OMPTBX_HAVE_OMP_H)
    s.attr("have_stubs_h") = false;
    s.attr("have_omp_h") = true;
#else
    s.attr("have_stubs_h") = true;
    s.attr("have_omp_h") = false;
#endif

    typedef bool (*boolean_property_getter)();
    typedef void (*boolean_property_setter)(bool);

    def("omp_set_num_threads", &omp_set_num_threads);
    def("omp_get_max_threads", &omp_get_max_threads);
    def("omp_get_num_procs", &omp_get_num_procs);
    def("omp_get_dynamic", (boolean_property_getter)&omp_get_dynamic);
    def("omp_set_dynamic", (boolean_property_setter)&omp_set_dynamic);
    def("omp_get_nested", (boolean_property_getter)&omp_get_nested);
    def("omp_set_nested", (boolean_property_setter)&omp_set_nested);

    def("tst_environment", &tst_environment);
  }

}}} // namespace omptbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(omptbx_ext)
{
  omptbx::boost_python::init_module();
}
