/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/math.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace cctbx { namespace miller { namespace boost_python {

  void wrap_asu();
  void wrap_bins();
  void wrap_expand_to_p1();
  void wrap_index_generator();
  void wrap_index_span();
  void wrap_match_bijvoet_mates();
  void wrap_match_indices();
  void wrap_phase_transfer();
  void wrap_sym_equiv();

namespace {

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping<af::shared<sym_equiv_index>, variable_capacity_policy>();
  }

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    register_tuple_mappings();

    wrap_sym_equiv(); // must be wrapped first to enable use of bases<>
    wrap_asu();
    wrap_bins();
    wrap_expand_to_p1();
    wrap_index_generator();
    wrap_index_span();
    wrap_match_bijvoet_mates();
    wrap_match_indices();
    wrap_phase_transfer();

    def("statistical_mean",
      (double(*)(sgtbx::space_group const&,
                 bool,
                 af::const_ref<index<> > const&,
                 af::const_ref<double> const&)) statistical_mean);
  }

} // namespace <anonymous>
}}} // namespace cctbx::miller::boost_python

BOOST_PYTHON_MODULE(miller_ext)
{
  cctbx::miller::boost_python::init_module();
}
