/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/phase_interpolation.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace miller { namespace boost_python {

  void wrap_phase_interpolation()
  {
    using namespace boost::python;
    def("phase_interpolation",
      (af::shared<double>(*)
        (af::const_ref<bool> const&,
         af::const_ref<double> const&,
         af::const_ref<std::complex<double> > const&,
         af::const_ref<std::complex<double> > const&,
         bool,
         double epsilon))
           phase_interpolation);
  }

}}} // namespace cctbx::miller::boost_python
