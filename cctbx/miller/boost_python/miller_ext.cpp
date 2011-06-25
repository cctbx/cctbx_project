#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/union_of_indices.h>
#include <cctbx/miller/math.h>
#include <cctbx/miller/image_simple.hpp>
#include <boost/python.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace cctbx { namespace miller { namespace boost_python {

  void wrap_asu();
  void wrap_bins();
  void wrap_change_basis();
  void wrap_expand_to_p1();
  void wrap_index_generator();
  void wrap_index_span();
  void wrap_match_bijvoet_mates();
  void wrap_match_indices();
  void wrap_merge_equivalents();
  void wrap_phase_transfer();
  void wrap_phase_integrator();
  void wrap_sym_equiv();
  void wrap_f_calc_map();
  // miller_lookup_utils
  void wrap_lookup_tensor();
  void wrap_local_neighbourhood();
  void wrap_local_area();
  void wrap_amplitude_normalisation();
  void wrap_slices();

namespace {

  struct union_of_indices_registry_wrappers
  {
    typedef union_of_indices_registry wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<wt>("union_of_indices_registry", no_init)
        .def(init<>())
        .def("update", &wt::update, (arg("indices")))
        .def("as_array", &wt::as_array)
      ;
    }
  };

  hendrickson_lattman<>
  as_hendrickson_lattman(
    bool centric_flag,
    std::complex<double> const& phase_integral,
    double const& max_figure_of_merit)
  {
    return hendrickson_lattman<>(
      centric_flag, phase_integral, max_figure_of_merit);
  }

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_variable_capacity<af::shared<sym_equiv_index> >();
  }

  void init_module()
  {
    using namespace boost::python;

    register_tuple_mappings();

    wrap_sym_equiv(); // must be wrapped first to enable use of bases<>
    wrap_asu();
    wrap_bins();
    wrap_change_basis();
    wrap_expand_to_p1();
    wrap_index_generator();
    wrap_index_span();
    wrap_match_bijvoet_mates();
    wrap_match_indices();
    wrap_merge_equivalents();
    wrap_phase_integrator();
    wrap_phase_transfer();
    wrap_f_calc_map();

    wrap_lookup_tensor();
    wrap_local_neighbourhood();
    wrap_local_area();
    wrap_amplitude_normalisation();
    wrap_slices();

    def("statistical_mean",
      (double(*)(sgtbx::space_group const&,
                 bool,
                 af::const_ref<index<> > const&,
                 af::const_ref<double> const&)) statistical_mean);

    union_of_indices_registry_wrappers::wrap();

    def("as_hendrickson_lattman", as_hendrickson_lattman, (
      arg("centric_flag"),
      arg("phase_integral"),
      arg("max_figure_of_merit")));

    {
      typedef image_simple wt;
      typedef return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;
      class_<wt>("image_simple", no_init)
        .def(init<bool, bool, bool, bool, bool, bool>((
          arg("apply_proximity_filter")=true,
          arg("apply_proximity_factor")=true,
          arg("store_miller_index_i_seqs")=false,
          arg("store_spots")=false,
          arg("store_signals")=false,
          arg("set_pixels")=false)))
        .def("compute", &wt::compute, (
          arg("unit_cell"),
          arg("miller_indices"),
          arg("spot_intensity_factors"),
          arg("crystal_rotation_matrix"),
          arg("ewald_radius"),
          arg("ewald_proximity"),
          arg("signal_max"),
          arg("detector_distance"),
          arg("detector_size"),
          arg("detector_pixels"),
          arg("point_spread"),
          arg("gaussian_falloff_scale")), rir())
        .add_property("miller_index_i_seqs",
          make_getter(&wt::miller_index_i_seqs, rbv()))
        .add_property("spots", make_getter(&wt::spots, rbv()))
        .add_property("signals", make_getter(&wt::signals, rbv()))
        .def_readonly("pixels", &wt::pixels)
      ;
    }
  }

} // namespace <anonymous>
}}} // namespace cctbx::miller::boost_python

BOOST_PYTHON_MODULE(cctbx_miller_ext)
{
  cctbx::miller::boost_python::init_module();
}
