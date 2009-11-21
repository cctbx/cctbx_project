#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/sampling_base.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct sampling_base_wrappers
  {
    typedef sampling_base<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("sampling_base", no_init)
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("u_base", &w_t::u_base)
        .def("u_extra", &w_t::u_extra)
        .def("u_min", &w_t::u_min)
        .def("wing_cutoff", &w_t::wing_cutoff)
        .def("exp_table_one_over_step_size",
          &w_t::exp_table_one_over_step_size)
        .def("tolerance_positive_definite",
          &w_t::tolerance_positive_definite)
        .def("n_scatterers_passed", &w_t::n_scatterers_passed)
        .def("n_contributing_scatterers", &w_t::n_contributing_scatterers)
        .def("n_anomalous_scatterers", &w_t::n_anomalous_scatterers)
        .def("anomalous_flag", &w_t::anomalous_flag)
        .def("exp_table_size", &w_t::exp_table_size)
        .def("max_sampling_box_n_points", &w_t::max_sampling_box_n_points)
        .def("sum_sampling_box_n_points", &w_t::sum_sampling_box_n_points)
        .def("ave_sampling_box_n_points", &w_t::ave_sampling_box_n_points)
        .def("max_sampling_box_edges", &w_t::max_sampling_box_edges, ccr())
        .def("max_sampling_box_edges_frac", &w_t::max_sampling_box_edges_frac)
        .def("excessive_sampling_radius_i_seqs",
          &w_t::excessive_sampling_radius_i_seqs, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sampling_base()
  {
    using namespace boost::python;

    def("calc_u_base", calc_u_base, (
      arg("d_min"),
      arg("grid_resolution_factor"),
      arg("quality_factor")=100,
      arg("max_u_base")=adptbx::b_as_u(1000)));

    def("apply_u_extra",
      (void(*)(
        uctbx::unit_cell const&,
        double const&,
        af::const_ref<miller::index<> > const&,
        af::ref<std::complex<double> > const&,
        double const&)) apply_u_extra, (
          arg("unit_cell"),
          arg("u_extra"),
          arg("miller_indices"),
          arg("structure_factors"),
          arg("multiplier")=1));

    sampling_base_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
