/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jul 2002: Created (rwgk)
 */

#include <boost/python/cross_module.hpp>

namespace {

  double waiting_loop(std::size_t n1, std::size_t n2)
  {
    double result = 0;
    for(std::size_t i=0;i<n1;i++) {
      double sum = 0;
      for(std::size_t j=0;j<n2;j++) {
        sum += 1.;
      }
      result += sum / n2;
    }
    return result;
  }

}

BOOST_PYTHON_MODULE_INIT(dev)
{
# include <cctbx/basic/from_bpl_import.h>

  python::module_builder this_module("dev");

  const std::string Revision = "$Revision$";
  this_module.add(ref(to_python(
      Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

  this_module.def(waiting_loop, "waiting_loop");
}
