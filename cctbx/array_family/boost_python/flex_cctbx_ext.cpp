/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created, based on scitbx/flex_module.cpp, shared_bpl.h (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/small.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/boost_python/utils.h>
#include <cctbx/coordinates.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/maptbx/accessors/c_grid_p1.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_hendrickson_lattman();
  void wrap_flex_miller_index(boost::python::object const& flex_root_scope);
  void wrap_flex_vec3_double();
  void wrap_flex_sym_mat3_double();
  void wrap_flex_tiny_size_t_2();
  void wrap_flex_xray_scatterer();

namespace {

  void register_cctbx_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_fixed_size<tiny<bool, 3> >();
    tuple_mapping_fixed_size<tiny<int, 12> >();
    tuple_mapping_fixed_size<tiny<std::size_t, 2> >();
    tuple_mapping_fixed_size<tiny<float, 4> >();
    tuple_mapping_fixed_size<tiny<float, 5> >();
    tuple_mapping_fixed_size<tiny<double, 12> >();

    tuple_mapping_fixed_capacity<small<double, 6> >();

    tuple_mapping_fixed_size<cctbx::cartesian<> >();
    tuple_mapping_fixed_size<cctbx::fractional<> >();
    tuple_mapping_fixed_size<cctbx::hendrickson_lattman<> >();
    tuple_mapping_fixed_size<cctbx::miller::index<> >();
  }

  void register_cctbx_c_grid_conversions()
  {
    using scitbx::af::boost_python::c_grid_flex_conversions;

    c_grid_flex_conversions<long, cctbx::maptbx::c_grid_p1<3> >();
    c_grid_flex_conversions<long, cctbx::maptbx::c_grid_padded_p1<3> >();
    c_grid_flex_conversions<float, cctbx::maptbx::c_grid_p1<3> >();
    c_grid_flex_conversions<float, cctbx::maptbx::c_grid_padded_p1<3> >();
    c_grid_flex_conversions<double, cctbx::maptbx::c_grid_p1<3> >();
    c_grid_flex_conversions<double, cctbx::maptbx::c_grid_padded_p1<3> >();
    c_grid_flex_conversions<std::complex<double>,
                            cctbx::maptbx::c_grid_p1<3> >();
    c_grid_flex_conversions<std::complex<double>,
                            cctbx::maptbx::c_grid_padded_p1<3> >();
  }

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    object flex_root_scope(scitbx::boost_python::import_module(
      "scitbx_boost.array_family.flex_scitbx_ext"));

    wrap_flex_hendrickson_lattman();
    wrap_flex_miller_index(flex_root_scope);
    wrap_flex_vec3_double();
    wrap_flex_sym_mat3_double();
    wrap_flex_tiny_size_t_2();
    wrap_flex_xray_scatterer();

    // The flex module will be used from all cctbx extension modules.
    // Therefore it is convenient to register all tuple mappings from here.
    register_cctbx_tuple_mappings();
    register_cctbx_c_grid_conversions();
  }

} // namespace <anonymous>
}}} // namespace scitbx::af::boost_python

BOOST_PYTHON_MODULE(flex_cctbx_ext)
{
  scitbx::af::boost_python::init_module();
}
