/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/translation_search/symmetry_flags.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

namespace {

  struct symmetry_flags_wrappers
  {
    typedef symmetry_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<maptbx::symmetry_flags> >("symmetry_flags", no_init)
        .def(init<bool, bool>())
        .def("is_isotropic_search_model", &w_t::is_isotropic_search_model)
        .def("have_f_part", &w_t::have_f_part)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symmetry_flags()
  {
    symmetry_flags_wrappers::wrap();
  }

}}} // namespace cctbx::translation_search::boost_python
