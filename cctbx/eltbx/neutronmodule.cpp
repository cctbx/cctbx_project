// $Id$

#include <boost/python/class_builder.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/neutron.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(neutron)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("neutron");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<NeutronNews1992Record>
    py_NeutronNews1992Record(this_module, "NeutronNews1992Record");

    py_NeutronNews1992Record.def(constructor<const std::string&>());
    py_NeutronNews1992Record.def(constructor<const std::string&, bool>());
    py_NeutronNews1992Record.def(&NeutronNews1992Record::Symbol, "Symbol");
    py_NeutronNews1992Record.def(&NeutronNews1992Record::BoundCohScattLength,
      "BoundCohScattLength");
    py_NeutronNews1992Record.def(&NeutronNews1992Record::AbsCrossSect,
      "AbsCrossSect");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
