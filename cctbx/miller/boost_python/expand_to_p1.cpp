/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/miller/expand_to_p1.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct expand_to_p1_wrappers
  {
    typedef expand_to_p1<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("expand_to_p1", no_init)
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<index<> > const&>())
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<index<> > const&,
                  af::const_ref<double> const&>())
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<index<> > const&,
                  af::const_ref<double> const&,
                  bool>())
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  optional<bool> >())
        .def(init<sgtbx::space_group const&,
                  bool,
                  af::const_ref<index<> > const&,
                  af::const_ref<std::complex<double> > const&>())
        .def("indices", &w_t::indices, ccr())
        .def("amplitudes", &w_t::amplitudes, ccr())
        .def("phases", &w_t::phases, ccr())
        .def("structure_factors", &w_t::structure_factors, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_expand_to_p1()
  {
    expand_to_p1_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
