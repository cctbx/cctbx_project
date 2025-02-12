#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/fft.h>
#include <cctbx/maptbx/average_densities.h>
#include <cctbx/maptbx/standard_deviations_around_sites.hpp>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/maptbx/histogram.h>
#include <cctbx/maptbx/resolution.h>
#include <cctbx/maptbx/fsc.h>
#include <cctbx/maptbx/loft.h>
#include <cctbx/maptbx/sphericity.h>
#include <cctbx/maptbx/mask.h>
#include <cctbx/maptbx/utils.h>
#include <cctbx/maptbx/connectivity.h>
#include <cctbx/maptbx/marked_grid_points.h>
#include <cctbx/maptbx/mask_utils.h>
#include <cctbx/maptbx/map_accumulator.h>
#include <cctbx/maptbx/ft_analytical_1d_point_scatterer_at_origin.h>
#include <cctbx/maptbx/target_and_gradients.h>

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

  template <typename FloatType, typename GridType>
  struct map_accumulator_wrapper
  {
    //typedef map_accumulator<FloatType, GridType> w_t; // works too
    typedef cctbx::maptbx::map_accumulator<FloatType, GridType> w_t;
    static void wrap() {
      using namespace boost::python;
      class_<w_t>("map_accumulator", no_init)
        .def(init<af::int3 const&,
             double const&,
             double const&,
             int const&,
             bool,
             bool >((arg("n_real"),
                     arg("smearing_b"),
                     arg("max_peak_scale"),
                     arg("smearing_span"),
                     arg("use_exp_table"),
                     arg("use_max_map"))))
        .def("as_median_map", &w_t::as_median_map)
        .def("add", &w_t::add, (arg("map_data")))
        .def("at_index", &w_t::at_index, (arg("n")))
        .def("int_to_float_at_index", &w_t::int_to_float_at_index, (arg("n")))
      ;
    }
  };

  template <typename FloatType>
  struct ft_analytical_1d_point_scatterer_at_origin_wrapper
  {
    typedef ft_analytical_1d_point_scatterer_at_origin<FloatType> w_t;
    static void wrap() {
      using namespace boost::python;
      class_<w_t>("ft_analytical_1d_point_scatterer_at_origin", no_init)
        .def(init<int const& >(
                    (arg("N"))))
        .def("distances", &w_t::distances)
        .def("rho", &w_t::rho)
        .def("compute", &w_t::compute,
          (arg("miller_indices"),arg("step"),arg("left"),
            arg("right"),arg("u_frac")))
      ;
    }
  };

namespace {

  void init_module()
  {
    using namespace boost::python;

    map_accumulator_wrapper<double, af::c_grid<3> >::wrap();
    ft_analytical_1d_point_scatterer_at_origin_wrapper<double>::wrap();
    wrap_grid_indices_around_sites();
    wrap_grid_tags();
    wrap_gridding();
    wrap_misc();
    wrap_peak_list();
    wrap_pymol_interface();
    wrap_statistics();
    wrap_structure_factors();
    wrap_coordinate_transformers();

// Real-space target and gradients ---------------------------------------------
    {
      typedef cctbx::maptbx::target_and_gradients::diffmap::compute w_t;
      class_<w_t>("target_and_gradients_diffmap", no_init)
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

    //- Magnification begin ----------------------------------------------------
    {
      typedef cctbx::maptbx::target_and_gradients::simple::magnification<double> w_t;
      class_<w_t>("target_and_gradients_simple_magnification", no_init)
        .def(init<uctbx::unit_cell const&,
             af::const_ref<double, af::c_grid_padded<3> > const&,
             af::const_ref<scitbx::vec3<double> > const&,
             scitbx::mat3<double> const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("sites_cart"),
                                    arg("K"))))
        .def(init<uctbx::unit_cell const&,
             af::const_ref<double, af::c_grid_padded<3> > const&,
             af::const_ref<scitbx::vec3<double> > const&,
             scitbx::vec3<double> const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("sites_cart"),
                                    arg("K"))))
        .def("target", &w_t::target)
        .def("gradients", &w_t::gradients)
      ;
    }

    def("magnification_isotropic",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           cctbx::maptbx::target_and_gradients::simple::magnification_isotropic, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart")));

    def("magnification_anisotropic",
      (scitbx::vec3<double>(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           cctbx::maptbx::target_and_gradients::simple::magnification_anisotropic, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart")));
    //- Magnification end ------------------------------------------------------

    {
      typedef cctbx::maptbx::target_and_gradients::simple::compute<double> w_t;
      class_<w_t>("target_and_gradients_simple", no_init)
        .def(init<uctbx::unit_cell const&,
             af::const_ref<double, af::c_grid_padded<3> > const&,
             af::const_ref<scitbx::vec3<double> > const&,
             double const&,
             af::const_ref<bool> const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("sites_cart"),
                                    arg("delta"),
                                    arg("selection"))))
        .def(init<uctbx::unit_cell const&,
             af::const_ref<double, af::c_grid_padded<3> > const&,
             af::const_ref<scitbx::vec3<double> > const&,
             double const&,
             af::const_ref<bool> const&,
             af::const_ref<double> const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("sites_cart"),
                                    arg("delta"),
                                    arg("selection"),
                                    arg("rsr_weight"))))
        .def(init<uctbx::unit_cell const&,
             af::const_ref<double, af::c_grid_padded<3> > const&,
             af::const_ref<scitbx::vec3<double> > const&,
             af::const_ref<bool> const&,
             std::string const& >((
                                    arg("unit_cell"),
                                    arg("map_target"),
                                    arg("sites_cart"),
                                    arg("selection"),
                                    arg("interpolation"))))
        .def("target", &w_t::target)
        .def("gradients", &w_t::gradients)
      ;
    }

    def("real_space_target_simple",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         af::const_ref<bool> const&))
           cctbx::maptbx::target_and_gradients::simple::target, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart"),
             arg("selection")));

    def("real_space_target_simple_with_adjacent_similarity",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         af::const_ref<std::size_t > const&,
         af::const_ref<double> const&))
           cctbx::maptbx::target_and_gradients::simple::target_with_adjacent_similarity, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart"),
             arg("selection"),
             arg("weights")));

    def("real_space_target_simple",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           cctbx::maptbx::target_and_gradients::simple::target, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart")));

    def("real_space_target_simple",
      (double(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         af::const_ref<std::size_t> const&))
           cctbx::maptbx::target_and_gradients::simple::target, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart"),
             arg("selection")));

    def("real_space_target_simple",
      (double(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           cctbx::maptbx::target_and_gradients::simple::target, (
             arg("density_map"),
             arg("sites_frac")));

    def("real_space_target_simple_per_site",
      (af::shared<double>(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&))
           cctbx::maptbx::target_and_gradients::simple::target_per_site, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart")));

    def("real_space_gradients_simple",
      (af::shared<scitbx::vec3<double> >(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&,
         double,
         af::const_ref<bool> const&))
           cctbx::maptbx::target_and_gradients::simple::gradients, (
             arg("unit_cell"),
             arg("density_map"),
             arg("sites_cart"),
             arg("delta"),
             arg("selection")));

// -----------------------------------------------------------------------------

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
      typedef non_linear_map_modification_to_match_average_cumulative_histogram w_t;

      class_<w_t>("non_linear_map_modification_to_match_average_cumulative_histogram", no_init)
        .def(init<af::const_ref<double, af::c_grid<3> > const&,
                  af::const_ref<double, af::c_grid<3> > const& >(
                    (arg("map_1"),
                     arg("map_2"))))
        .def("map_1", &w_t::map_1)
        .def("map_2", &w_t::map_2)
        .def("histogram_1", &w_t::histogram_1)
        .def("histogram_2", &w_t::histogram_2)
        .def("histogram_12", &w_t::histogram_12)
        .def("histogram_values", &w_t::histogram_values)
      ;
    }

    {
      typedef d99 w_t;

      class_<w_t>("d99", no_init)
        .def(init<
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<miller::index<> > const&,
          double const& >(
                    (arg("f"),
                     arg("d_spacings"),
                     arg("hkl"),
                     arg("cutoff"))))
        .def("d_min", &w_t::d_min)
      ;
    }

    {
      typedef fsc w_t;

      class_<w_t>("fsc", no_init)
        .def(init<
          af::const_ref<std::complex<double> > const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          int const& >(
                    (arg("f1"),
                     arg("f2"),
                     arg("d_spacings"),
                     arg("step"))))
        .def("fsc",   &w_t::cc)
        .def("d",     &w_t::d)
        .def("d_inv", &w_t::d_inv)
      ;
    }

    {
      typedef histogram w_t;

      class_<w_t>("histogram", no_init)
        .def(init<af::const_ref<double, af::c_grid<3> > const&,
                  int const& >(
                    (arg("map"),
                     arg("n_bins"))))
        .def(init<af::const_ref<double> const&,
                  int const& >(
                    (arg("map"),
                     arg("n_bins"))))
        .def("values",    &w_t::values)
        .def("c_values",  &w_t::c_values)
        .def("v_values",  &w_t::v_values)
        .def("bin_width", &w_t::bin_width)
        .def("arguments", &w_t::arguments)
      ;
    }

    {
      typedef sphericity2 w_t;

      class_<w_t>("sphericity2", no_init)
        .def(init<af::const_ref<double, af::c_grid<3> > const&,
                  cctbx::cartesian<double> const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  cctbx::uctbx::unit_cell const& >(
                    (arg("map_data"),
                     arg("center_cart"),
                     arg("points_on_sphere_cart"),
                     arg("unit_cell"))))
        .def("rho_min_max_mean", &w_t::rho_min_max_mean)
        .def("ccs_min_max_mean", &w_t::ccs_min_max_mean)
      ;
    }

    {
      typedef fit_point_3d_grid_search w_t;

      class_<w_t>("fit_point_3d_grid_search", no_init)
        .def(init<cctbx::cartesian<> const&,
                  af::const_ref<double, af::c_grid<3> > const&,
                  double const&,
                  cctbx::uctbx::unit_cell const&,
                  double const&,
                  double const& >(
                    (arg("site_cart"),
                     arg("map_data"),
                     arg("map_min"),
                     arg("unit_cell"),
                     arg("amplitude"),
                     arg("increment"))))
        .def("has_peak",        &w_t::has_peak)
        .def("map_best",        &w_t::map_best)
        .def("map_start",       &w_t::map_start)
        .def("site_cart_moved", &w_t::site_cart_moved)
      ;
    }

    {
      typedef volume_scale w_t;

      class_<w_t>("volume_scale", no_init)
        .def(init<af::const_ref<double, af::c_grid<3> > const&,
                  int const& >(
                    (arg("map"),
                     arg("n_bins"))))
        .def("map_data", &w_t::map_data)
        .def("v_values", &w_t::v_values)
      ;
    }

    {
      typedef volume_scale_1d w_t;

      class_<w_t>("volume_scale_1d", no_init)
        .def(init<af::const_ref<double> const&,
                  int const& >(
                    (arg("map"),
                     arg("n_bins"))))
        .def("map_data", &w_t::map_data)
        .def("v_values", &w_t::v_values)
      ;
    }

    {
      typedef loft w_t;

      class_<w_t>("loft", no_init)
        .def(init<af::const_ref<cctbx::miller::index<> > const&,
                  af::const_ref<int, af::c_grid<3> > const&,
                  af::shared<double> const&,
                  double const& >(
                    (arg("miller_indices"),
                     arg("map_data"),
                     arg("abc"),
                     arg("d_min"))))
        .def("structure_factors", &w_t::structure_factors)
        .def("indices", &w_t::indices)
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

    {
      typedef connectivity w_t;

      class_<w_t>("connectivity", no_init)
        .def(init<af::ref<float, af::c_grid<3> >,
                  float const&,
                  bool,
                  bool >(
                    (arg("map_data"),
                     arg("threshold"),
                     arg("wrapping")=true,
                     arg("preprocess_against_shallow")=false)))
        .def(init<af::ref<double, af::c_grid<3> > ,
                  double const&,
                  bool,
                  bool >(
                    (arg("map_data"),
                     arg("threshold"),
                     arg("wrapping")=true,
                     arg("preprocess_against_shallow")=false)))
        .def(init<af::ref<int, af::c_grid<3> >,
                  int const&,
                  bool,
                  bool >(
                    (arg("map_data"),
                     arg("threshold"),
                     arg("wrapping")=true,
                     arg("preprocess_against_shallow")=false)))

        // Symmetry-aware constructors
        // .def(init< //af::ref<float, af::c_grid<3> >
        //           float const&
        //           // const cctbx::sgtbx::space_group &,
        //           // const cctbx::sgtbx::space_group_type &,
        //           // int3,
        //           // bool,
        //           // bool
        //           >(
        //             // (arg("map_data")
        //              arg("threshold")
        //              // arg("space_group"),
        //              // arg("uc_dimensions"),
        //              // arg("wrapping")=true,
        //              // arg("preprocess_against_shallow")=false
        //             )))

        .def(init< af::ref<double, af::c_grid<3> >
                  // double const&
                  // const cctbx::sgtbx::space_group &,
                  // const cctbx::sgtbx::space_group_type &,
                  // int3,
                  // bool,
                  // bool
                  >(
                    (arg("map_data")
                      // arg("threshold")
                     // arg("space_group"),
                     // arg("uc_dimensions"),
                     // arg("wrapping")=true,
                     // arg("preprocess_against_shallow")=false
                    )))

        // amap constructor
        .def(init< asymmetric_map &,
                  double const&,
                  bool
                  >(
                    (arg("amap"),
                     arg("threshold"),
                     arg("preprocess_against_shallow")=false
                    )))



        .def("result",    &w_t::result)
        .def("regions",   &w_t::regions)
        .def("experiment_with_symmetry", &w_t::experiment_with_symmetry,
            (arg("space_group")))
        .def("merge_symmetry_related_regions", &w_t::merge_symmetry_related_regions,
            (arg("space_group")))
        .def("volume_cutoff_mask", &w_t::volume_cutoff_mask,
                    (arg("volume_cutoff")))
        .def("get_blobs_boundaries", &w_t::get_blobs_boundaries)
        .def("expand_mask", &w_t::expand_mask,
                  (arg("id_to_expand"),
                   arg("expand_size")))
        .def("noise_elimination_two_cutoffs",
             &w_t::noise_elimination_two_cutoffs,
                    (arg("connectivity_object_at_t1"),
                     arg("elimination_volume_threshold_at_t1"),
                     arg("zero_all_interblob_region")=true))
        .def("maximum_coors", &w_t::maximum_coors)
        .def("maximum_values", &w_t::maximum_values)
        .def_readonly("preprocessing_changed_voxels", &w_t::preprocessing_changed_voxels)
        .def_readonly("preprocessing_n_passes", &w_t::preprocessing_n_passes)
      ;
    }
    {
      typedef marked_grid_points w_t;

      class_<w_t>("marked_grid_points", no_init)
        .def(init<af::const_ref<bool, af::flex_grid<> > const&,
                  int const&
                  >(
                    (arg("map_data"),
                     arg("every_nth_point")
                     )))
        .def("result",    &w_t::result)
      ;
    }

    {
      typedef sample_all_mask_regions w_t;
      class_<w_t>("sample_all_mask_regions", no_init)
        .def(init<af::const_ref<int, af::flex_grid<> > const&,
            af::shared<int> const&,
            af::shared<int> const&,
            cctbx::uctbx::unit_cell const& > ((
              arg("mask"),
              arg("volumes"),
              arg("sampling_rates"),
              arg("unit_cell"))))
        .def("get_array", &w_t::get_array,
          (arg("n")))
      ;
    }

    {
      typedef zero_boundary_box_map w_t;
      class_<w_t>("zero_boundary_box_map", no_init)
        .def(init<af::const_ref<double, af::flex_grid<> > const&,
            int const& > ((
              arg("mask"),
              arg("boundary"))))
        .def("result",    &w_t::result)
      ;
    }

    {
      typedef binary_filter w_t;
      class_<w_t>("binary_filter", no_init)
        .def(init<af::const_ref<double, af::flex_grid<> > const&,
            float const& > ((
              arg("map"),
              arg("threshold"))))
        .def("result",    &w_t::result)
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
    def("copy_box",
      (af::versa<float, af::flex_grid<> >(*)
        (af::const_ref<float, af::flex_grid<> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy_box, (
      arg("map"),
      arg("first"),
      arg("last")));
    def("copy_box",
      (af::versa<double, af::flex_grid<> >(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         af::int3 const&,
         af::int3 const&)) maptbx::copy_box, (
      arg("map"),
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

    def("direct_summation_at_point",
      (std::complex<double>(*)(
        af::const_ref<miller::index<> > const&,
        af::const_ref<std::complex<double> > const&,
        scitbx::vec3<double>))
          maptbx::direct_summation_at_point, (
            arg("miller_indices"),
            arg("data"),
            arg("site_frac")));

    def("cc_weighted_maps",cc_weighted_maps);

    def("kuwahara_filter",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         int const&)) kuwahara_filter, (
      arg("map_data"),
      arg("index_span")));

    def("median_filter",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         int const&)) median_filter, (
      arg("map_data"),
      arg("index_span")));

    def("remove_single_node_peaks",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         af::ref<double, af::c_grid_padded<3> >,
         double const&,
         int const&)) remove_single_node_peaks, (
      arg("map_data"),
      arg("mask_data"),
      arg("cutoff"),
      arg("index_span")));

    def("map_box_average",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double const&,
         int const&)) map_box_average, (
      arg("map_data"),
      arg("cutoff"),
      arg("index_span")));

    def("map_box_average",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         int const&)) map_box_average, (
      arg("map_data"),
      arg("index_span")));

    def("center_of_mass",
      (cctbx::cartesian<>(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         uctbx::unit_cell const&,
         double const&)) center_of_mass, (
      arg("map_data"),
      arg("unit_cell"),
      arg("cutoff")));

    def("sharpen",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         int const&,
         int const&,
         bool)) sharpen, (
      arg("map_data"),
      arg("index_span"),
      arg("n_averages"),
      arg("allow_negatives")));

    def("gamma_compression",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double const&)) gamma_compression, (
      arg("map_data"),
      arg("gamma")));

    def("map_box_average",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         cctbx::uctbx::unit_cell const&,
         double const&)) map_box_average, (
      arg("map_data"),
      arg("unit_cell"),
      arg("radius")));

    def("hoppe_gassman_modification",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double,
         int)) hoppe_gassman_modification, (
      arg("data"),
      arg("mean_scale"),
      arg("n_iterations")));

    def("hoppe_gassman_modification2",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         int)) hoppe_gassman_modification2, (
      arg("data"),
      arg("n_iterations")));

    def("sphericity_tensor",
      (scitbx::sym_mat3<double>(*)
        (af::const_ref<double, af::c_grid<3> > const&,
          uctbx::unit_cell const&,
         double const&,
         cctbx::fractional<> const&)) sphericity_tensor, (
      arg("map_data"),
      arg("unit_cell"),
      arg("radius"),
      arg("site_frac")));

    def("sphericity",
      (af::shared<double>(*)
        (af::const_ref<double, af::c_grid<3> > const&,
          uctbx::unit_cell const&,
         double const&,
         af::const_ref<scitbx::vec3<double> > const&)) sphericity, (
      arg("map_data"),
      arg("unit_cell"),
      arg("radius"),
      arg("sites_frac")));

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

    def("mask",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<scitbx::vec3<double> > const&,
          uctbx::unit_cell const&,
          af::tiny<int, 3> const&,
          double const&,
          double const&,
          af::const_ref<double> const&,
          bool const&)) mask, (
      arg("sites_frac"),
      arg("unit_cell"),
      arg("n_real"),
      arg("mask_value_inside_molecule"),
      arg("mask_value_outside_molecule"),
      arg("radii"),
      arg("wrapping")=true
       ));

    def("convert_to_non_negative",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double)) convert_to_non_negative, (
      arg("data"),
      arg("substitute_value")));

    def("flexible_boundary_mask",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         af::ref<double, af::c_grid<3> >)) flexible_boundary_mask, (
      arg("data"),
      arg("mask")));

    def("negate_selected_in_place",
      (af::versa<double, af::c_grid_padded<3> >(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         std::vector<unsigned> const&)) negate_selected_in_place, (
      arg("map_data"),
      arg("selection")));

    def("resample",
      (void(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::ref<double, af::c_grid<3> >,
         cctbx::uctbx::unit_cell const&)) resample, (
      arg("map_data"),
      arg("map_data_new"),
      arg("unit_cell")));

    def("reset",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double,
         double,
         double,
         bool)) reset, (
      arg("data"),
      arg("substitute_value"),
      arg("less_than_threshold"),
      arg("greater_than_threshold"),
      arg("use_and")));

    def("intersection",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         af::ref<double, af::c_grid<3> >,
         af::ref<double> const&,
         bool)) intersection, (
      arg("map_data_1"),
      arg("map_data_2"),
      arg("thresholds"),
      arg("average")));

    def("rotate_translate_map",
      (af::versa<double, af::c_grid<3> >(*)
        (uctbx::unit_cell const&,
         af::const_ref<double, af::c_grid<3> > const&,
         scitbx::mat3<double> const&,
         scitbx::vec3<double> const&,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) rotate_translate_map, (
      arg("unit_cell"),
      arg("map_data"),
      arg("rotation_matrix"),
      arg("translation_vector"),
      arg("start"),
      arg("end")));

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
         scitbx::vec3<double> const&,
         bool wrap )) superpose_maps, (
      arg("unit_cell_1"),
      arg("unit_cell_2"),
      arg("map_data_1"),
      arg("n_real_2"),
      arg("rotation_matrix"),
      arg("translation_vector"),
      arg("wrapping")=true));

    def("combine_and_maximize_maps",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::tiny<int, 3> const& )) combine_and_maximize_maps, (
      arg("map_data_1"),
      arg("map_data_2"),
      arg("n_real")));

    def("denmod_simple",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::tiny<int, 3> const&,
         double,double)) denmod_simple, (
      arg("map_data"),
      arg("n_real"),
      arg("cutoffp"),
      arg("cutoffm")));

    def("binarize",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double const&,
         double const&,
         double const&)) binarize, (
      arg("map_data"),
      arg("threshold"),
      arg("substitute_value_below"),
      arg("substitute_value_above")));

    def("truncate_special",
      (void(*)
        (af::ref<int, af::c_grid<3> >,
         af::ref<double, af::c_grid<3> >)) truncate_special, (
      arg("mask"),
      arg("map_data")));

    def("combine_1",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         af::ref<double, af::c_grid<3> >)) combine_1, (
      arg("map_data"),
      arg("diff_map")));

    def("truncate_between_min_max",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double const&,
         double const&)) truncate_between_min_max, (
      arg("map_data"),
      arg("min"),
      arg("max")));

    def("truncate",
      (void(*)
        (af::ref<double, af::c_grid<3> >,
         double const&,
         double const&,
         double const&,
         double const&)) truncate, (
      arg("map_data"),
      arg("standard_deviation"),
      arg("by_sigma_less_than"),
      arg("scale_by"),
      arg("set_value")));

    def("fem_averaging_loop",
      (af::shared<std::complex<double> >(*)
        (af::const_ref<std::complex<double> > const&,
         af::const_ref<double> const&,
         af::const_ref<double> const&,
         double const&,
         int const&,
         int const&)) fem_averaging_loop, (
      arg("map_coefficients"),
      arg("r_factors"),
      arg("sigma_over_f_obs"),
      arg("random_scale"),
      arg("random_seed"),
      arg("n_cycles")));

    def("conditional_solvent_region_filter",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         double const&)) conditional_solvent_region_filter, (
      arg("bulk_solvent_mask"),
      arg("map_data"),
      arg("threshold")));

    def("update_f_part1_helper",
      (af::versa<double, af::c_grid<3> >(*)
        (af::const_ref<int, af::c_grid_padded<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         int const&)) update_f_part1_helper, (
      arg("connectivity_map"),
      arg("map_data"),
      arg("region_id")));

    def("map_sum_at_sites_frac",
      (double(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<scitbx::vec3<double> > const&)) map_sum_at_sites_frac, (
      arg("map_data"),
      arg("sites_frac")));

    def("discrepancy_function",
      (af::shared<double>(*)
        (af::const_ref<double> const&,
         af::const_ref<double> const&,
         af::const_ref<double> const&)) discrepancy_function, (
      arg("map_1"),
      arg("map_2"),
      arg("cutoffs")));

    def("discrepancy_function",
      (af::shared<double>(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<double> const&)) discrepancy_function, (
      arg("map_1"),
      arg("map_2"),
      arg("cutoffs")));

    def("cc_complex_complex",
      (af::shared<double>(*)
        (af::const_ref<std::complex<double> > const&,
         af::const_ref<std::complex<double> > const&,
         af::const_ref<double> const&,
         af::const_ref<double> const&,
         af::const_ref<double> const&,
         double const&)) cc_complex_complex, (
      arg("f_1"),
      arg("f_2"),
      arg("d_spacings"),
      arg("ss"),
      arg("d_mins"),
      arg("b_iso")));

    def("cc_complex_complex",
      (double(*)
        (af::const_ref<std::complex<double> > const&,
         af::const_ref<std::complex<double> > const&)) cc_complex_complex, (
      arg("f_1"),
      arg("f_2")));

    def("cc_peak",
      (double(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::const_ref<double, af::c_grid<3> > const&,
         double const&)) cc_peak, (
      arg("map_1"),
      arg("map_2"),
      arg("cutoff")));

    def("set_box_copy",
      (af::versa<double, af::c_grid<3> >(*)
        (double const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) set_box_copy, (
      arg("value"),
      arg("map_data_to"),
      arg("start"),
      arg("end")));

    def("set_box_copy_inside",
      (af::versa<double, af::c_grid<3> >(*)
        (double const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) set_box_copy_inside, (
      arg("value"),
      arg("map_data_to"),
      arg("start"),
      arg("end")));

    def("set_box",
      (void(*)
        (double const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) set_box, (
      arg("value"),
      arg("map_data_to"),
      arg("start"),
      arg("end")));

    def("set_box",
      (void(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) set_box, (
      arg("map_data_from"),
      arg("map_data_to"),
      arg("start"),
      arg("end")));

    def("set_box_with_symmetry",
      (void(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&,
         cctbx::uctbx::unit_cell const& unit_cel,
         af::shared<scitbx::mat3<double> > const&,
         af::shared<scitbx::vec3<double> > const&
         )) set_box_with_symmetry, (
      arg("map_data_from"),
      arg("map_data_to"),
      arg("start"),
      arg("end"),
      arg("unit_cell"),
      arg("rotation_matrix"),
      arg("translation_vector")
      ));

    def("copy_box",
      (void(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         af::ref<double, af::c_grid<3> >,
         af::tiny<int, 3> const&,
         af::tiny<int, 3> const&)) copy_box, (
      arg("map_data_from"),
      arg("map_data_to"),
      arg("start"),
      arg("end")));

    def("eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::flex_grid<> > const&,
         scitbx::vec3<double> const&)) eight_point_interpolation);
    def("eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&)) eight_point_interpolation);
    def("eight_point_interpolation",
      (double(*)
        (af::const_ref<double, af::c_grid<3> > const&,
         scitbx::vec3<double> const&)) eight_point_interpolation);
    def("eight_point_interpolation_with_gradients",
      (af::tiny<double, 4>(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&,
         scitbx::vec3<double> const&)) eight_point_interpolation_with_gradients);
    def("quadratic_interpolation_with_gradients",
      (af::tiny<double, 4>(*)
        (//af::const_ref<double, af::flex_grid<> > const&,
         af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&,
         scitbx::vec3<double> const&)) quadratic_interpolation_with_gradients);
    def("closest_grid_point",
      (af::c_grid_padded<3>::index_type(*)
        (af::flex_grid<> const&,
         fractional<double> const&)) closest_grid_point);
    def("tricubic_interpolation_with_gradients",
      (af::tiny<double, 4>(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&,
         scitbx::vec3<double> const&)) tricubic_interpolation_with_gradients);
    def("tricubic_interpolation",
      (double(*)
        (af::const_ref<double, af::c_grid_padded<3> > const&,
         scitbx::vec3<double> const&)) tricubic_interpolation);
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
