#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/asu_mappings.h>
#include <cctbx/crystal/direct_space_asu.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_wrappers
  {
    typedef scatterer<> w_t;
    typedef w_t::float_type flt_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      apply_symmetry_overloads, apply_symmetry, 2, 6)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      apply_symmetry_u_star_overloads, apply_symmetry_u_star, 2, 5)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("scatterer", no_init)
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  flt_t const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>((
          arg_("label"),
          arg_("site"),
          arg_("u_iso"),
          arg_("occupancy"),
          arg_("scattering_type"),
          arg_("fp"),
          arg_("fdp"))))
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  scitbx::sym_mat3<flt_t> const&,
                  flt_t const&,
                  std::string const&,
                  flt_t const&,
                  flt_t const&>((
          arg_("label"),
          arg_("site"),
          arg_("u_star"),
          arg_("occupancy"),
          arg_("scattering_type"),
          arg_("fp"),
          arg_("fdp"))))
        .add_property("label", make_getter(&w_t::label, rbv()),
                               make_setter(&w_t::label, dcp()))
        .add_property("scattering_type",
          make_getter(&w_t::scattering_type, rbv()),
          make_setter(&w_t::scattering_type, dcp()))
        .add_property("fp", make_getter(&w_t::fp, rbv()),
                            make_setter(&w_t::fp, dcp()))
        .add_property("fdp", make_getter(&w_t::fdp, rbv()),
                             make_setter(&w_t::fdp, dcp()))
        .add_property("site", make_getter(&w_t::site, rbv()),
                              make_setter(&w_t::site, dcp()))
        .add_property("occupancy", make_getter(&w_t::occupancy, rbv()),
                                   make_setter(&w_t::occupancy, dcp()))
        .add_property("anisotropic_flag",
          make_getter(&w_t::anisotropic_flag, rbv()),
          make_setter(&w_t::anisotropic_flag, dcp()))
        .add_property("u_iso", make_getter(&w_t::u_iso, rbv()),
                               make_setter(&w_t::u_iso, dcp()))
        .add_property("u_star", make_getter(&w_t::u_star, rbv()),
                                make_setter(&w_t::u_star, dcp()))
        .def("apply_symmetry", &w_t::apply_symmetry,
          apply_symmetry_overloads((
            arg_("unit_cell"),
            arg_("space_group"),
            arg_("min_distance_sym_equiv")=0.5,
            arg_("u_star_tolerance")=0,
            arg_("assert_is_positive_definite")=false,
            arg_("assert_min_distance_sym_equiv")=true)))
        .def("apply_symmetry_site", &w_t::apply_symmetry_site, (
          arg_("site_symmetry_ops")))
        .def("apply_symmetry_u_star", &w_t::apply_symmetry_u_star,
          apply_symmetry_u_star_overloads((
            arg_("unit_cell"),
            arg_("site_symmetry_ops"),
            arg_("u_star_tolerance")=0,
            arg_("assert_is_positive_definite")=false,
            arg_("assert_min_distance_sym_equiv")=true)))
        .def("multiplicity", &w_t::multiplicity)
        .def("weight_without_occupancy", &w_t::weight_without_occupancy)
        .def("weight", &w_t::weight)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scatterer()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround

    scatterer_wrappers::wrap();

    def("asu_mappings_process",
      (void(*)(
        crystal::direct_space_asu::asu_mappings<>&,
        af::const_ref<scatterer<> > const&,
        sgtbx::site_symmetry_table const&)) asu_mappings_process, (
      arg_("asu_mappings"), arg_("scatterers"), arg_("site_symmetry_table")));
  }

}}} // namespace cctbx::xray::boost_python
