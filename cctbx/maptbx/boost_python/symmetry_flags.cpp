/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/symmetry_flags.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct symmetry_flags_wrappers
  {
    typedef symmetry_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("symmetry_flags", no_init)
        .def(init<bool, optional<bool, bool> >())
        .def("use_space_group_symmetry", &w_t::use_space_group_symmetry)
        .def("use_normalizer_k2l", &w_t::use_normalizer_k2l)
        .def("use_structure_seminvariants", &w_t::use_structure_seminvariants)
        .def("select_sub_space_group", &w_t::select_sub_space_group)
        .def("grid_factors", &w_t::grid_factors)
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symmetry_flags()
  {
    symmetry_flags_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
