/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/translation_search/map_gridding.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

namespace {

  struct map_gridding_wrappers
  {
    typedef map_gridding<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("map_gridding", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group_type const&,
                  maptbx::symmetry_flags const&,
                  double,
                  af::const_ref<miller::index<> > const&,
                  int>())
        .def("target", &w_t::target, ccr())
        .def("quarter", &w_t::quarter, ccr())
        .def("eighth", &w_t::eighth, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_map_gridding()
  {
    map_gridding_wrappers::wrap();
  }

}}} // namespace cctbx::translation_search::boost_python
