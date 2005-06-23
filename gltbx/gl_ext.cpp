#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>

namespace gltbx { namespace gl {

  namespace boost_python {

    void wrap_defines_00(boost::python::scope scope);
    void wrap_defines_01(boost::python::scope scope);
    void wrap_defines_02(boost::python::scope scope);
    void wrap_defines_03(boost::python::scope scope);
    void wrap_defines_04(boost::python::scope scope);
    void wrap_defines_05(boost::python::scope scope);
    void wrap_defines_06(boost::python::scope scope);
    void wrap_defines_07(boost::python::scope scope);
    void wrap_functions_00();
    void wrap_functions_01();
    void wrap_functions_02();
    void wrap_functions_03();
    void wrap_functions_04();
    void wrap_functions_05();
    void wrap_functions_06();
    void wrap_functions_07();
    void wrap_functions_08();
    void wrap_functions_09();
    void wrap_functions_10();
    void wrap_functions_11();
    void wrap_functions_12();
    void wrap_functions_13();
    void wrap_functions_14();
    void wrap_functions_15();

  }

  void
  init_module()
  {
    using namespace boost::python;
    boost_python::wrap_defines_00(scope());
    boost_python::wrap_defines_01(scope());
    boost_python::wrap_defines_02(scope());
    boost_python::wrap_defines_03(scope());
    boost_python::wrap_defines_04(scope());
    boost_python::wrap_defines_05(scope());
    boost_python::wrap_defines_06(scope());
    boost_python::wrap_defines_07(scope());
    boost_python::wrap_functions_00();
    boost_python::wrap_functions_01();
    boost_python::wrap_functions_02();
    boost_python::wrap_functions_03();
    boost_python::wrap_functions_04();
    boost_python::wrap_functions_05();
    boost_python::wrap_functions_06();
    boost_python::wrap_functions_07();
    boost_python::wrap_functions_08();
    boost_python::wrap_functions_09();
    boost_python::wrap_functions_10();
    boost_python::wrap_functions_11();
    boost_python::wrap_functions_12();
    boost_python::wrap_functions_13();
    boost_python::wrap_functions_14();
    boost_python::wrap_functions_15();
  }

}} // namespace gltbx::gl

BOOST_PYTHON_MODULE(gltbx_gl_ext)
{
  gltbx::gl::init_module();
}
