/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/sampled_model_density.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct sampled_model_density_wrappers
  {
    typedef sampled_model_density<> w_t;
    typedef w_t::grid_point_type grid_point_type;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("sampled_model_density", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  grid_point_type const&,
                  grid_point_type const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool,
                           bool> >())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("u_extra", &w_t::u_extra)
        .def("wing_cutoff", &w_t::wing_cutoff)
        .def("exp_table_one_over_step_size",
          &w_t::exp_table_one_over_step_size)
        .def("n_scatterers_passed", &w_t::n_scatterers_passed)
        .def("n_contributing_scatterers", &w_t::n_contributing_scatterers)
        .def("n_anomalous_scatterers", &w_t::n_anomalous_scatterers)
        .def("anomalous_flag", &w_t::anomalous_flag)
        .def("real_map", &w_t::real_map)
        .def("complex_map", &w_t::complex_map)
        .def("exp_table_size", &w_t::exp_table_size)
        .def("max_shell_radii", &w_t::max_shell_radii, ccr())
        .def("max_shell_radii_frac", &w_t::max_shell_radii_frac)
        .def("apply_symmetry",
          (void(w_t::*)(maptbx::grid_tags<> const&))
            &w_t::apply_symmetry)
        .def("eliminate_u_extra_and_normalize",
          &w_t::eliminate_u_extra_and_normalize)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    calc_u_extra_overloads, calc_u_extra, 2, 4)

  void
  eliminate_u_extra_double(
    uctbx::unit_cell const& unit_cell,
    double const& u_extra,
    af::const_ref<miller::index<> > const& miller_indices,
    af::ref<std::complex<double> > const& structure_factors,
    double const& norm=1)
  {
    eliminate_u_extra(unit_cell, u_extra, miller_indices, structure_factors,
                      norm);
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    eliminate_u_extra_double_overloads, eliminate_u_extra_double, 4, 5)

} // namespace <anoymous>

  void wrap_sampled_model_density()
  {
    using namespace boost::python;

    def("calc_u_extra", calc_u_extra, calc_u_extra_overloads());
    def("eliminate_u_extra", eliminate_u_extra_double,
                             eliminate_u_extra_double_overloads());

    sampled_model_density_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
