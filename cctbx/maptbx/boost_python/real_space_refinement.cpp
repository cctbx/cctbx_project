/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/real_space_refinement.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct real_space_refinement_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      def("residual", real_space_refinement::residual<double>, (
        arg_("map"), arg_("gridding_matrix"), arg_("sites_cart")));
    }
  };

} // namespace <anoymous>

  void wrap_real_space_refinement()
  {
    real_space_refinement_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
