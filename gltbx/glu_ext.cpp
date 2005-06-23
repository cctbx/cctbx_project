#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>

namespace gltbx { namespace glu {

  namespace boost_python {

    void wrap_defines_00(boost::python::scope scope);
    void wrap_defines_01(boost::python::scope scope);
    void wrap_functions_00();
    void wrap_functions_01();
    void wrap_functions_02();
    void wrap_functions_03();
    void wrap_functions_04();

  }

  void
  init_module()
  {
    using namespace boost::python;
    boost_python::wrap_defines_00(scope());
    boost_python::wrap_defines_01(scope());
    boost_python::wrap_functions_00();
    boost_python::wrap_functions_01();
    boost_python::wrap_functions_02();
    boost_python::wrap_functions_03();
  }

}} // namespace gltbx::glu

BOOST_PYTHON_MODULE(gltbx_glu_ext)
{
  gltbx::glu::init_module();
}
