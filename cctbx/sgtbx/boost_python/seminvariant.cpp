/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/seminvariant.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct ss_vec_mod_wrappers
  {
    typedef ss_vec_mod w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("ss_vec_mod", no_init)
        .def_readonly("v", &w_t::v)
        .def_readonly("m", &w_t::m)
      ;
    }
  };

  struct structure_seminvariant_wrappers
  {
    typedef structure_seminvariant w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("structure_seminvariant", no_init)
        .def(init<space_group const&>())
        .def("vectors_and_moduli", &w_t::vectors_and_moduli, ccr())
        .def("size", &w_t::size)
        .def("is_ss", &w_t::is_ss)
        .def("apply_mod", &w_t::apply_mod)
        .def("gridding", &w_t::gridding)
        .def("refine_gridding",
          (sg_vec3(w_t::*)(sg_vec3 const&) const)
          &w_t::refine_gridding)
        .def("grid_adapted_moduli",
          (af::small<ss_vec_mod, 3>(w_t::*)(sg_vec3 const&) const)
          &w_t::grid_adapted_moduli)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_seminvariant()
  {
    ss_vec_mod_wrappers::wrap();
    structure_seminvariant_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
