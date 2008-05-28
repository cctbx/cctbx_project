#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <omp.h>

namespace scitbx { namespace openmp { namespace boost_python { namespace {

  void wrap_runtime_library_functions() {
    using namespace boost::python;
    typedef bool (*boolean_property_getter)();
    typedef void (*boolean_property_setter)(bool);

    def("omp_set_num_threads", &omp_set_num_threads);
    def("omp_get_max_threads", &omp_get_max_threads);
    def("omp_get_num_procs", &omp_get_num_procs);
    def("omp_get_dynamic", (boolean_property_getter)&omp_get_dynamic);
    def("omp_get_nested", (boolean_property_getter)&omp_get_nested);
    def("omp_set_dynamic", (boolean_property_setter)&omp_set_dynamic);
    def("omp_set_nested", (boolean_property_setter)&omp_set_nested);
  }

} // namespace <anonymous>

void init_module() {
  wrap_runtime_library_functions();
}

}}} // namespace scitbx::openmp::boost_python

BOOST_PYTHON_MODULE(scitbx_openmp_ext) {
  scitbx::openmp::boost_python::init_module();
}
