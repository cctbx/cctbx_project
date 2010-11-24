#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/fft.h>
#include <cctbx/maptbx/average_densities.h>
#include <cctbx/maptbx/standard_deviations_around_sites.hpp>
#include <cctbx/maptbx/real_space_gradients_simple.h>
#include <cctbx/maptbx/real_space_target_and_gradients.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

  void wrap_grid_indices_around_sites();
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

  void init_module()
  {
    using namespace boost::python;

    wrap_grid_indices_around_sites();
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

    {
      typedef target_and_gradients w_t;
      class_<w_t>("target_and_gradients", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<double, af::c_grid_padded<3> > const&,
                  af::const_ref<double, af::c_grid_padded<3> > const&,
                  double const&,
                  af::const_ref<scitbx::vec3<double> > const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("map_current"),
                                    arg("step"),
                                    arg("sites_frac"))))
        .def("target", &w_t::target)
        .def("gradients", &w_t::gradients)
      ;
    }
    {
      typedef grid_points_in_sphere_around_atom_and_distances w_t;

      class_<w_t>("grid_points_in_sphere_around_atom_and_distances", no_init)
        .def(init<cctbx::uctbx::unit_cell const&,
                  af::const_ref<double, af::c_grid<3> > const&,
                  double const&,
                  double const&,
                  scitbx::vec3<double> const& >(
                    (arg("unit_cell"),
                     arg("data"),
                     arg("radius"),
                     arg("shell"),
                     arg("site_frac"))))
        .def("data_at_grid_points", &w_t::data_at_grid_points)
        .def("data_at_grid_points_averaged", &w_t::data_at_grid_points_averaged)
        .def("distances", &w_t::distances)
      ;
    }
    {
      typedef one_gaussian_peak_approximation w_t;

      class_<w_t>("one_gaussian_peak_approximation", no_init)
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  bool const&,
                  bool const& >(
                    (arg("data_at_grid_points"),
                     arg("distances"),
                     arg("use_weights"),
                     arg("optimize_cutoff_radius"))))
        .def("a_real_space", &w_t::a_real_space)
        .def("b_real_space", &w_t::b_real_space)
        .def("a_reciprocal_space", &w_t::a_reciprocal_space)
        .def("b_reciprocal_space", &w_t::b_reciprocal_space)
        .def("gof", &w_t::gof)
        .def("cutoff_radius", &w_t::cutoff_radius)
        .def("weight_power", &w_t::weight_power)
        .def("first_zero_radius", &w_t::first_zero_radius)
      ;
    }

    def("copy",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, af::flex_grid<> > const&,
         af::flex_grid<> const&)) maptbx::copy, (
      arg("map"),
      arg("result_grid")));
    def("copy",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         af::flex_grid<> const&)) maptbx::copy, (
      arg("map"),
      arg("result_grid")));
    def("copy",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, c_grid_padded_p1<3> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy, (
      arg("map_unit_cell"),
      arg("first"),
      arg("last")));
    def("copy",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, c_grid_padded_p1<3> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy, (
      arg("map_unit_cell"),
      arg("first"),
      arg("last")));
    def("unpad_in_place",
      (void(*)(af::versa<float, af::flex_grid<> >&))
        maptbx::unpad_in_place, (arg("map")));
    def("unpad_in_place",
      (void(*)(af::versa<double, af::flex_grid<> >&))
        maptbx::unpad_in_place, (arg("map")));

    def("fft_to_real_map_unpadded",
      (af::versa<double, af::c_grid<3> >(*)(
        sgtbx::space_group const&,
        af::tiny<int, 3> const&,
        af::const_ref<miller::index<> > const&,
        af::const_ref<std::complex<double> > const&))
          maptbx::fft_to_real_map_unpadded, (
            arg("space_group"),
            arg("n_real"),
            arg("miller_indices"),
            arg("data")));

    def("box_map_averaging",box_map_averaging);
    def("average_densities",
      (af::shared<double>(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         float)) average_densities, (
      arg("unit_cell"),
      arg("data"),
      arg("sites_frac"),
      arg("radius")));

    def("rotate_translate_map",
      (af::versa<double, af::c_grid<3> >(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid<3> > const&,
         scitbx::mat3<double> const&,
         scitbx::vec3<double> const& )) rotate_translate_map, (
      arg("unit_cell"),
      arg("map_data"),
      arg("rotation_matrix"),
      arg("translation_vector")));

    def("superpose_maps",
      (af::versa<double, af::c_grid<3> >(*)
        (uctbx::unit_cell const&,
         uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::tiny<int, 3> const&,
         scitbx::mat3<double> const&,
         scitbx::vec3<double> const& )) superpose_maps, (
      arg("unit_cell_1"),
      arg("unit_cell_2"),
      arg("map_data_1"),
      arg("n_real_2"),
      arg("rotation_matrix"),
      arg("translation_vector")));

    def("combine_and_maximize_maps",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::tiny<int, 3> const& )) combine_and_maximize_maps, (
      arg("map_data_1"),
      arg("map_data_2"),
      arg("n_real")));

    def("eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&)) eight_point_interpolation);
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
           non_crystallographic_eight_point_interpolation, (
             arg("map"),
             arg("gridding_matrix"),
             arg("site_cart"),
             arg("allow_out_of_bounds")=false,
             arg("out_of_bounds_substitute_value")=0));
    def("asu_eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         crystal::direct_space_asu::asu_mappings<double> &,
         fractional<double> const&)) asu_eight_point_interpolation);
    def("real_space_target_simple",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           real_space_target_simple, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart")));
    def("real_space_gradients_simple",
      (af::shared<scitbx::vec3<double> >(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         double)) real_space_gradients_simple, (
           arg("unit_cell"),
           arg("density_map"),
           arg("sites_cart"),
           arg("delta")));

    def("standard_deviations_around_sites",
      standard_deviations_around_sites, (
        arg("unit_cell"),
        arg("density_map"),
        arg("sites_cart"),
        arg("site_radii")));
  }

} // namespace <anonymous>
}}} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_maptbx_ext)
{
  cctbx::maptbx::boost_python::init_module();
}
