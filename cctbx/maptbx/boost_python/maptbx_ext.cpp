/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/copy.h>
#include <cctbx/maptbx/gridding.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

  void wrap_grid_tags();
  void wrap_peak_list();
  void wrap_pymol_interface();
  void wrap_statistics();
  void wrap_structure_factors();
  void wrap_symmetry_flags();

namespace {

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    wrap_grid_tags();
    wrap_peak_list();
    wrap_pymol_interface();
    wrap_statistics();
    wrap_structure_factors();
    wrap_symmetry_flags();

    def("copy",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, af::flex_grid<> > const& map,
         af::flex_grid<> const& result_grid))
           maptbx::copy);
    def("copy",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, af::flex_grid<> > const& map,
         af::flex_grid<> const& result_grid))
           maptbx::copy);

    def("determine_grid",
      (af::int3(*)(uctbx::unit_cell const&,
                   double,
                   double,
                   int,
                   af::int3 const&))
        determine_grid);
  }

} // namespace <anonymous>
}}} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(maptbx_ext)
{
  cctbx::maptbx::boost_python::init_module();
}
