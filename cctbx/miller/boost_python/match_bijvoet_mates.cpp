/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/miller/match_bijvoet_mates.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct match_bijvoet_mates_wrappers
  {
    typedef match_bijvoet_mates w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("match_bijvoet_mates", no_init)
        .def(init<sgtbx::space_group_type const&,
                  af::shared<index<> > const&>())
        .def(init<sgtbx::reciprocal_space::asu const&,
                  af::shared<index<> > const&>())
        .def(init<af::shared<index<> > const&>())
        .def("pairs", &w_t::pairs)
        .def("singles", &w_t::singles)
        .def("have_singles", &w_t::have_singles)
        .def("miller_indices_in_hemisphere",&w_t::miller_indices_in_hemisphere)
#define CCTBX_DEF(function_name) \
        .def(# function_name, \
          (af::shared<double>(w_t::*)(af::const_ref<double> const&) const) \
          &w_t::function_name)
        CCTBX_DEF(minus)
        CCTBX_DEF(additive_sigmas)
        CCTBX_DEF(average)
#undef CCTBX_DEF
      ;
    }
  };

} // namespace <anoymous>

  void wrap_match_bijvoet_mates()
  {
    match_bijvoet_mates_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
