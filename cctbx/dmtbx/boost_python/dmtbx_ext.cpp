/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (Ralf W. Grosse-Kunstleve)
 */

#include <boost/python/module.hpp>

namespace cctbx { namespace dmtbx { namespace boost_python {

  void wrap_triplet_generator();
  void wrap_triplet_phase_relation();

namespace {

  void init_module()
  {
    wrap_triplet_generator();
    wrap_triplet_phase_relation();
  }

} // namespace <anonymous>
}}} // namespace cctbx::dmtbx::boost_python

BOOST_PYTHON_MODULE(dmtbx_ext)
{
  cctbx::dmtbx::boost_python::init_module();
}
