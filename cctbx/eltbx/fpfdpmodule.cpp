// $Id$

#include <boost/python/cross_module.hpp>
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/efpfdp.h>

using namespace eltbx;

BOOST_PYTHON_MODULE_INIT(fpfdp)
{
  try
  {
#   include <cctbx/basic/from_bpl_import.h>

    python::module_builder this_module("fpfdp");

    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<fpfdp> py_fpfdp(this_module, "fpfdp");
    python::export_converters(py_fpfdp);

    py_fpfdp.def(constructor<float, float>());
    py_fpfdp.def(&fpfdp::fp, "fp");
    py_fpfdp.def(&fpfdp::fdp, "fdp");
    py_fpfdp.def(&fpfdp::isValid_fp, "fp");
    py_fpfdp.def(&fpfdp::isValid_fdp, "fdp");
  }
  catch(...)
  {
    boost::python::handle_exception(); // Deal with the exception for Python
  }
}
