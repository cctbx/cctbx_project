// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/basic/meta.h>
#include <cctbx/shared_storage_bpl.h>

using namespace cctbx;

namespace {

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    class_builder<shared_storage_handle_type>
    py_shared_storage_handle(this_module, "shared_storage_handle");

    py_shared_storage_handle.def(constructor<>());

#define WRAP_TYPE(python_name, type) \
    python::wrap_shared_storage(this_module, python_name, type_holder<type >())

    WRAP_TYPE("int", int);
    WRAP_TYPE("float", float);
    WRAP_TYPE("double", double);
    WRAP_TYPE("complex", std::complex<double>);
  }

}

BOOST_PYTHON_MODULE_INIT(shared_storage)
{
  boost::python::module_builder this_module("shared_storage");
  init_module(this_module);
}
