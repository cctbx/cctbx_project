// $Id$

#include <boost/python/cross_module.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/sasaki.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(sasaki)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("sasaki");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<fpfdp> fpfdp_converters("fpfdp", "fpfdp");
    class_builder<Sasaki> py_Sasaki(this_module, "Sasaki");

    py_Sasaki.def(constructor<const std::string&>());
    py_Sasaki.def(constructor<const std::string&, bool>());
    py_Sasaki.def(&Sasaki::Label, "Label");
    py_Sasaki.def(&Sasaki::Z, "Z");
    py_Sasaki.def(&Sasaki::operator(), "__call__");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
