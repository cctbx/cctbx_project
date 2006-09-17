#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/copy.h>
#include <cctbx/maptbx/average_densities.h>
#include <cctbx/maptbx/eight_point_interpolation.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>


#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/tls/tls.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>


namespace cctbx { namespace maptbx { namespace boost_python {

  void wrap_grid_tags();
  void wrap_gridding();
  void wrap_misc();
  void wrap_peak_list();
  void wrap_pymol_interface();
  void wrap_statistics();
  void wrap_structure_factors();
  void wrap_coordinate_transformers();
  void wrap_mappers();
  void wrap_basic_map();
  void wrap_real_space_refinement();

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    non_crystallographic_eight_point_interpolation_overloads,
    non_crystallographic_eight_point_interpolation, 3, 5)

  void init_module()
  {
    using namespace boost::python;

    wrap_grid_tags();
    wrap_gridding();
    wrap_misc();
    wrap_peak_list();
    wrap_pymol_interface();
    wrap_statistics();
    wrap_structure_factors();
    wrap_coordinate_transformers();
    wrap_mappers();
    wrap_basic_map();
    wrap_real_space_refinement();


    using namespace boost::python;
    typedef boost::python::arg arg_;
    class_<grid_points_in_sphere_around_atom_and_distances>("grid_points_in_sphere_around_atom_and_distances",
                           init<cctbx::uctbx::unit_cell const&,
                                af::const_ref<double, af::c_grid<3> > const&,
                                double const&,
                                double const&,
                                vec3<double> const& >(
                                   (arg_("unit_cell"),
                                                         arg_("data"),
                                                         arg_("radius"),
                                                         arg_("shell"),
                                                         arg_("site_frac"))))
      .def("data_at_grid_points", &grid_points_in_sphere_around_atom_and_distances::data_at_grid_points)
      .def("distances", &grid_points_in_sphere_around_atom_and_distances::distances)
    ;
    class_<one_gaussian_peak_approximation>("one_gaussian_peak_approximation",
                           init<af::const_ref<double> const&,
                                af::const_ref<double> const& >(
                                                         (arg_("data_at_grid_points"),
                                                         arg_("distances"))))
      .def("a_real_space", &one_gaussian_peak_approximation::a_real_space)
      .def("b_real_space", &one_gaussian_peak_approximation::b_real_space)
      .def("a_reciprocal_space", &one_gaussian_peak_approximation::a_reciprocal_space)
      .def("b_reciprocal_space", &one_gaussian_peak_approximation::b_reciprocal_space)
      .def("gof", &one_gaussian_peak_approximation::gof)
    ;


    def("copy",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, af::flex_grid<> > const&,
         af::flex_grid<> const&)) maptbx::copy, (
      arg_("map"),
      arg_("result_grid")));
    def("copy",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         af::flex_grid<> const&)) maptbx::copy, (
      arg_("map"),
      arg_("result_grid")));
    def("copy",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, c_grid_padded_p1<3> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy, (
      arg_("map_unit_cell"),
      arg_("first"),
      arg_("last")));
    def("copy",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, c_grid_padded_p1<3> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy, (
      arg_("map_unit_cell"),
      arg_("first"),
      arg_("last")));
    def("unpad_in_place",
      (void(*)(af::versa<float, af::flex_grid<> >&))
        maptbx::unpad_in_place, (arg_("map")));
    def("unpad_in_place",
      (void(*)(af::versa<double, af::flex_grid<> >&))
        maptbx::unpad_in_place, (arg_("map")));

    def("box_map_averaging",box_map_averaging);
    def("average_densities",
      (af::shared<double>(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         float)) average_densities, (
      arg_("unit_cell"),
      arg_("data"),
      arg_("sites_frac"),
      arg_("radius")));

    def("eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         fractional<double> const&)) eight_point_interpolation);
    def("closest_grid_point",
      (af::c_grid_padded<3>::index_type(*)
        (af::flex_grid<> const&,
         fractional<double> const&)) closest_grid_point);
    def("non_crystallographic_eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         scitbx::mat3<double> const&,
         scitbx::vec3<double> const&,
         bool,
         double const&))
           non_crystallographic_eight_point_interpolation,
         non_crystallographic_eight_point_interpolation_overloads((
           arg_("map"),
           arg_("gridding_matrix"),
           arg_("site_cart"),
           arg_("allow_out_of_bounds")=false,
           arg_("out_of_bounds_substitute_value")=0)));
    def("asu_eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         crystal::direct_space_asu::asu_mappings<double> &,
         fractional<double> const&)) asu_eight_point_interpolation);
  }

} // namespace <anonymous>
}}} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_maptbx_ext)
{
  cctbx::maptbx::boost_python::init_module();
}
