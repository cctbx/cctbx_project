#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/wigner/wigner3j.h>

namespace scitbx { namespace wigner {
namespace boost_python{

  struct wigner3j_fast_wrapper
  {
    typedef wigner3j_fast < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j_fast", no_init)
        .def( init< int const&
                  >
             (( arg("max")
             ))
            )
        .def("compute", &w_t::compute)
      ;
    }
  };

  struct wigner3j_wrapper
  {
    typedef wigner3j < double > w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("wigner3j", no_init)
        .def( init< int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&,
                    int const&
                  >
             (( arg("j1"),
                arg("j2"),
                arg("j3"),
                arg("m1"),
                arg("m2"),
                arg("m3")
             ))
            )
        .def("check", &w_t::check)
        .def("get_value", &w_t::get_value)
      ;
    }
  };

  void
  init_module()
  {
    wigner3j_fast_wrapper::wrap();
    wigner3j_wrapper::wrap();
  }

}}} // scitbx::wigner::boost_python

BOOST_PYTHON_MODULE(scitbx_wigner_ext)
{
  scitbx::wigner::boost_python::init_module();
}
