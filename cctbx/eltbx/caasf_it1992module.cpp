// $Id$

#include <boost/python/class_builder.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/caasf.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(caasf_it1992)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("caasf_it1992");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<CAASF<4> > py_CAASF_4(this_module, "CAASF_4");
    class_builder<CAASF_IT1992> py_CAASF_IT1992(this_module, "CAASF_IT1992");

    py_CAASF_IT1992.declare_base(py_CAASF_4, python::without_downcast);

    py_CAASF_IT1992.def(constructor<const std::string&>());
    py_CAASF_IT1992.def(constructor<const std::string&, bool>());
    py_CAASF_IT1992.def(&CAASF_IT1992::Table, "Table");
    py_CAASF_IT1992.def(&CAASF_IT1992::Label, "Label");
    py_CAASF_IT1992.def(&CAASF_IT1992::n_ab, "n_ab");
    py_CAASF_IT1992.def(&CAASF_IT1992::a, "a");
    py_CAASF_IT1992.def(&CAASF_IT1992::b, "b");
    py_CAASF_IT1992.def(&CAASF_IT1992::c, "c");
    py_CAASF_IT1992.def(&CAASF_IT1992::operator(), "__call__");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
