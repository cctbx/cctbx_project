// $Id$

#include <boost/python/cross_module.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/henke.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(henke)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("henke");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<fpfdp> fpfdp_converters("fpfdp", "fpfdp");
    class_builder<Henke> py_Henke(this_module, "Henke");

    py_Henke.def(constructor<const std::string&>());
    py_Henke.def(constructor<const std::string&, bool>());
    py_Henke.def(&Henke::Label, "Label");
    py_Henke.def(&Henke::Z, "Z");
    py_Henke.def(&Henke::operator(), "__call__");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
