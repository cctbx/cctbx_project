#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <iotbx/shelx/hklf.h>

namespace iotbx { namespace shelx { namespace boost_python {

  struct hklf_reader_wrapper
  {
    typedef hklf_reader wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<wt>("hklf_reader", no_init)
        .def(init<af::const_ref<std::string> const&, bool>((
          arg("lines"),
          arg("strict")=true)))
        .def("indices", &wt::indices)
        .def("data", &wt::data)
        .def("sigmas", &wt::sigmas)
        .def("alphas", &wt::alphas)
        .def("batch_numbers", &wt::batch_numbers)
        ;
    }
  };

  void
  init_module()
  {
    hklf_reader_wrapper::wrap();
  }

}}} //iotbx::shelx::boost_python

BOOST_PYTHON_MODULE(iotbx_shelx_ext)
{
  iotbx::shelx::boost_python::init_module();
}
