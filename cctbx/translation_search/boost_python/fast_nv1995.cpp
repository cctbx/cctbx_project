/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/translation_search/fast_nv1995.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

namespace {

  struct fast_nv1995_wrappers
  {
    typedef fast_nv1995<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("fast_nv1995", no_init)
        .def(init<map_gridding<> const&,
                  sgtbx::space_group const&,
                  bool,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<std::complex<double> > >())
        .def("target_map", &w_t::target_map)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_fast_nv1995()
  {
    fast_nv1995_wrappers::wrap();
  }

}}} // namespace cctbx::translation_search::boost_python
