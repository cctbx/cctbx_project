#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <iotbx/cif/parser.h>

namespace iotbx { namespace cif { namespace boost_python {

  struct cif_wrapper
  {
    typedef iotbx::cif::parser wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<std::string, boost::python::object&>((arg("input"), arg("builder"))))
        ;
    }
  };

  void init_module() {
    cif_wrapper::wrap("fast_reader");
  }

}}} //iotbx::cif::boost_python


BOOST_PYTHON_MODULE(iotbx_cif_ext)
{
        iotbx::cif::boost_python::init_module();
}
