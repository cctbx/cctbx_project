#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <iotbx/shelx/hklf_simple.h>
//#define IOTBX_SHELX_HKLF_SPIRIT_DEBUG
#include <iotbx/shelx/hklf.h>

namespace iotbx { namespace shelx { namespace boost_python {

  template<class WrappedType>
  struct any_hklf_reader_wrapper
  {
    typedef WrappedType wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<std::string const &, bool>((arg("content"),
                                              arg("strict")=true)))
        .def("indices", &wt::indices)
        .def("data", &wt::data)
        .def("sigmas", &wt::sigmas)
        .def("alphas", &wt::alphas)
        .def("batch_numbers", &wt::batch_numbers)
        ;
    }
  };

  void init_module() {
    any_hklf_reader_wrapper<simple_hklf_reader>::wrap("simple_hklf_reader");
    any_hklf_reader_wrapper<hklf_reader>::wrap("hklf_reader");
  }

}}} //iotbx::shelx::boost_python


BOOST_PYTHON_MODULE(iotbx_shelx_ext)
{
        iotbx::shelx::boost_python::init_module();
}
