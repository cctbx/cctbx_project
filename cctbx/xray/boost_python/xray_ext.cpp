/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/xray/scatterer.h>
#include <cctbx/xray/scatterer_utils.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_conversions();
  void wrap_sampled_model_density();
  void wrap_scatterer();
  void wrap_structure_factors();
  void wrap_targets();

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    apply_symmetry_overloads, apply_symmetry, 3, 7)

  // work around Visual C++ 7 internal compiler error
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    structure_factor_array_overloads, structure_factor_array, 4, 4)

namespace {

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    scitbx::boost_python::import_module(
      "cctbx_boost.eltbx.caasf_ext");

    wrap_conversions();
    wrap_sampled_model_density();
    wrap_scatterer();
    wrap_structure_factors();
    wrap_targets();

    def("apply_symmetry",
      (af::shared<std::size_t>(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group const&,
        af::ref<scatterer<> > const&,
        double, double, bool, bool)) 0, apply_symmetry_overloads());

    def("rotate",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        scitbx::mat3<double> const&,
        af::const_ref<scatterer<> > const&)) rotate);

    def("pack_parameters",
      (std::size_t(*)(
        const uctbx::unit_cell*,
        af::const_ref<scatterer<> > const&,
        af::versa<double, af::flex_grid<> >&,
        bool, bool, bool)) pack_parameters);

    def("unpack_parameters",
      (std::size_t(*)(
        const uctbx::unit_cell*,
        std::size_t,
        af::const_ref<double> const& x,
        std::size_t start,
        af::ref<scatterer<> > const&,
        bool, bool, bool)) unpack_parameters);
  }

} // namespace <anonymous>
}}} // namespace cctbx::xray::boost_python

BOOST_PYTHON_MODULE(xray_ext)
{
  cctbx::xray::boost_python::init_module();
}
