#include <boost/python/cross_module.hpp>
#include <vector>

#include <iostream>
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

namespace {

  typedef std::vector<double> shared_real_array;
  typedef std::vector<std::complex<double> > shared_complex_array;

  void show_complex(shared_complex_array a) {
CheckPoint;
    std::cout << a.size() << std::endl;
CheckPoint;
  }

  void show_real(shared_real_array a) {
CheckPoint;
    std::cout << a.size() << std::endl;
CheckPoint;
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<std::vector<double> >
    py_shared_double(
      "cctbx.arraytbx.std_vector", "double");

    python::import_converters<std::vector<std::complex<double> > >
    py_shared_complex_double(
      "cctbx.arraytbx.std_vector", "complex_double");

    this_module.def(show_complex, "show");
    this_module.def(show_real, "show");
  }

}

BOOST_PYTHON_MODULE_INIT(debug_overload)
{
  boost::python::module_builder this_module("debug_overload");
  init_module(this_module);
}
